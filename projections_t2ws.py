#!/usr/bin/env python
"""
Thomas Klijnsma
"""

#____________________________________________________________________
# Imports

import logging
import os, sys, re, copy
from OptionHandler import flag_as_option
import differentials
import differentialutils
import LatestBinning
import LatestPaths

import random
random.seed(1002)
from differentials.theory.theory_utils import FileFinder


#____________________________________________________________________
cards = differentials.core.AttrDict()
cards.hgg = 'suppliedInput/fromVittorio/pT_newBins_Feb28/projection_hgg_Jun28.txt'
cards.hzz = 'suppliedInput/fromDavid/PTH_Jan24_newBinning/smH/projection_hzz_Jun28.txt'
cards.hbb = 'suppliedInput/fromJavier/bernstein_r7428/projection_hbb_Jun28.txt'
cards.combWithHbb = 'suppliedInput/projection_combWithHbb_Jun28.txt'

@flag_as_option
def projection_pth_smH_t2ws(args):
    differentials.combine.combine_utils.set_global_outdir('projections')
    decay_channel = differentialutils.get_decay_channel_tag(args)
    t2ws = differentials.combine.t2ws.T2WS(
        cards[decay_channel],
        name = 'ws_pth_smH_' + decay_channel
        )
    if args.hzz:
        t2ws.add_maps(
            '.*/smH_PTH_0_15:r_smH_PTH_0_15[1.0,-1.0,4.0]',
            '.*/smH_PTH_15_30:r_smH_PTH_15_30[1.0,-1.0,4.0]',
            '.*/smH_PTH_30_45:r_smH_PTH_30_80[1.0,-1.0,4.0]',
            '.*/smH_PTH_45_80:r_smH_PTH_30_80[1.0,-1.0,4.0]',
            '.*/smH_PTH_80_120:r_smH_PTH_80_200[1.0,-1.0,4.0]',
            '.*/smH_PTH_120_200:r_smH_PTH_80_200[1.0,-1.0,4.0]',
            '.*/smH_PTH_200_350:r_smH_PTH_GT200[1.0,-1.0,4.0]',
            '.*/smH_PTH_350_600:r_smH_PTH_GT200[1.0,-1.0,4.0]',
            '.*/smH_PTH_GT600:r_smH_PTH_GT200[1.0,-1.0,4.0]',
            )
    else:
        t2ws.make_maps_from_processes(scale_ggH_xH_with_smH=True)

    t2ws.run()


#____________________________________________________________________
class T2WSKBKC(differentials.combine.t2ws.T2WS):
    """docstring for T2WSKBKC"""

    def __init__(self):
        super(T2WSKBKC, self).__init__()

        self.exp_binning = [ 0., 15., 30., 45., 80., 120. ]
        self.model_file = 'physicsModels/CouplingModel.py'
        self.model_name = 'couplingModel'

        self.extra_options.extend([
            '--PO linearTerms=True',
            '--PO splitggH=True',
            ])

        self.obs = copy.deepcopy(LatestBinning.obs_pth_ggH)
        self.obs.drop_bins_up_to_value(self.exp_binning[-1]+1.0)
        self.extra_options.append(
            '--PO binBoundaries={0}'
            .format(','.join([str(b) for b in self.obs.binning]))
            )
        
        self.is_reweighted = False
        self.tags.append('unreweighted')

    def add_theory(self):
        coupling_variations = FileFinder(
            muR=1.0, muF=1.0, Q=1.0, directory=LatestPaths.theory.yukawa.filedir
            ).get()

        self.sm = [ v for v in coupling_variations if v.kappab==1.0 and v.kappac==1.0 ][0]
        coupling_variations.pop(coupling_variations.index(self.sm))

        self.extra_options.append('--PO SM=[kappab=1.0,kappac=1.0,file={0}]'.format(self.sm.theory_file))
        for variation in random.sample(coupling_variations, 7):
            if variation.kappab == 0.0 and variation.kappac == 0.0: continue
            self.extra_options.append(
                '--PO theory=[kappab={0},kappac={1},file={2}]'
                .format(variation.kappab, variation.kappac, variation.theory_file)
                )

    def apply_reweighting(self):
        crosssections = self.obs.crosssection()
        crosssections_str = ','.join([str(v) for v in crosssections])
        self.extra_options.append('--PO SMXS_of_input_ws={0}'.format(crosssections_str))
        self.tags[self.tags.index('unreweighted')] = 'reweighted'
        self.is_reweighted = True

    def add_theory_uncertainties(self, uncorrelated=False):
        # Theory uncertainties
        correlation_matrix   = LatestPaths.theory.yukawa.correlation_matrix
        theory_uncertainties = LatestPaths.theory.yukawa.uncertainties
        if uncorrelated:
            correlation_matrix = LatestPaths.theory.yukawa.correlation_matrix_uncorrelated
            self.tags.append('uncorrelatedTheoryUnc')
        self.extra_options.append('--PO correlationMatrix={0}'.format(correlation_matrix))
        self.extra_options.append('--PO theoryUncertainties={0}'.format(theory_uncertainties))


cards.yukawa = differentials.core.AttrDict()
cards.yukawa.hgg         = 'suppliedInput/projection_yukawa_hgg_Jul09.txt'
cards.yukawa.hzz         = 'suppliedInput/projection_yukawa_hzz_Jul09.txt'
cards.yukawa.combination = 'suppliedInput/projection_yukawa_combination_Jul09.txt'

@flag_as_option
def projection_kbkc_t2ws_couplingdependentBRs(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    t2ws = T2WSKBKC()
    t2ws.card = cards.yukawa[decay_channel]
    if args.hgg: t2ws.extra_options.append('--PO isOnlyHgg=True' )
    if args.hzz: t2ws.extra_options.append('--PO isOnlyHZZ=True' )
    t2ws.add_theory()
    t2ws.apply_reweighting()
    t2ws.add_theory_uncertainties()
    t2ws.extra_options.append('--PO BRs_kappa_dependent=True')
    t2ws.tags.append('couplingdependentBRs')
    t2ws.run()

@flag_as_option
def projection_kbkc_t2ws_floatingBRs(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    t2ws = T2WSKBKC()
    t2ws.card = cards.yukawa[decay_channel]
    if args.hgg: t2ws.extra_options.append('--PO isOnlyHgg=True' )
    if args.hzz: t2ws.extra_options.append('--PO isOnlyHZZ=True' )
    t2ws.add_theory()
    t2ws.apply_reweighting()
    t2ws.add_theory_uncertainties()
    t2ws.extra_options.append('--PO freely_floating_BRs=True')
    t2ws.tags.append('floatingBRs')
    t2ws.run()




#____________________________________________________________________
class T2WSKTCGKB(differentials.combine.t2ws.T2WS):
    """docstring for T2WSKTCGKB"""

    def __init__(self):
        super(T2WSKTCGKB, self).__init__()

        self.model_file = 'physicsModels/CouplingModel.py'
        self.model_name = 'couplingModel'

        self.extra_options.extend([
            '--PO linearTerms=False',
            '--PO splitggH=True',
            ])

        self.obs = copy.deepcopy(LatestBinning.obs_pth_ggH)
        self.extra_options.append(
            '--PO binBoundaries={0}'
            .format(','.join([str(b) for b in self.obs.binning]))
            )
        
        self.is_reweighted = False
        self.tags.append('unreweighted')

    def add_theory_ktcg(self):
        logging.warning('Doing variations for kappa_t / kappa_g')
        coupling_variations = FileFinder(
            cb=1.0, muR=1.0, muF=1.0, Q=1.0, directory=LatestPaths.theory.top.filedir
            ).get()
        self.sm = [ v for v in coupling_variations if v.ct==1.0 and v.cg==0.0 ][0]
        coupling_variations.pop(coupling_variations.index(self.sm))

        self.extra_options.append('--PO SM=[ct=1,cg=0,file={0}]'.format(self.sm.theory_file))
        for variation in coupling_variations:
            if variation.ct == 1.0 or variation.cg == 0.0: continue
            self.extra_options.append(
                '--PO theory=[ct={0},cg={1},file={2}]'
                .format(variation.ct, variation.cg, variation.theory_file)
                )

    def add_theory_ktkb(self):
        logging.warning('Doing variations for kappa_t / kappa_b')
        coupling_variations = FileFinder(
            cg=0.0, muR=1.0, muF=1.0, Q=1.0, directory=LatestPaths.theory.top.filedir
            ).get()
        self.sm = [ v for v in coupling_variations if v.ct==1.0 and v.cb==1.0 ][0]
        coupling_variations.pop(coupling_variations.index(self.sm))

        self.extra_options.append('--PO SM=[ct=1,cb=1,file={0}]'.format(self.sm.theory_file))
        for variation in coupling_variations:
            if variation.ct == 1.0 or variation.cb == 1.0: continue
            self.extra_options.append(
                '--PO theory=[ct={0},cb={1},file={2}]'
                .format(variation.ct, variation.cb, variation.theory_file)
                )

    def apply_reweighting(self):
        crosssections = self.obs.crosssection()
        crosssections_str = ','.join([str(v) for v in crosssections])
        self.extra_options.append('--PO SMXS_of_input_ws={0}'.format(crosssections_str))
        self.tags[self.tags.index('unreweighted')] = 'reweighted'
        self.is_reweighted = True

    def add_theory_uncertainties(self, uncorrelated=False):
        correlation_matrix   = 'out/scalecorrelations_Mar06/corrMat_tophighpt.txt'
        theory_uncertainties = 'out/scalecorrelations_Mar06/errors_tophighpt.txt'
        if uncorrelated:
            correlation_matrix = LatestPaths.correlationMatrix_Yukawa_Uncorrelated # Uncorrelated
            self.tags.append('uncorrelatedTheoryUnc')
        self.extra_options.append('--PO correlationMatrix={0}'.format(correlation_matrix))
        self.extra_options.append('--PO theoryUncertainties={0}'.format(theory_uncertainties))

    def add_scaling_ttH(self):
        self.extra_options.append('--PO add_scaling_ttH=True')
        # self.tags.append('scalingttH')

    def add_scaling_bbH(self):
        self.extra_options.append('--PO add_scaling_bbH=True')
        # self.tags.append('scalingbbH')

    def constrain_bbZZ_ratio(self):
        self.extra_options.append('--PO constrain_ratio_bb_ZZ=True')
        self.tags.append('constrainedbbZZ')


cards.ktcgkb = differentials.core.AttrDict()
cards.ktcgkb.hgg         = 'suppliedInput/fromVittorio/pT_newBins_Feb28/projection_ktcgkb_hgg_Jul10.txt'
cards.ktcgkb.hzz         = 'suppliedInput/fromDavid/PTH_Jan24_newBinning/ggH/projection_ktcgkb_hzz_Jul10.txt'
cards.ktcgkb.combWithHbb = 'suppliedInput/projection_ktcgkb_combWithHbb_Jul10.txt'

@flag_as_option
def projection_ktcg_t2ws_couplingdependentBRs(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    t2ws = T2WSKTCGKB()
    t2ws.card = cards.ktcgkb[decay_channel]
    if args.hgg: t2ws.extra_options.append('--PO isOnlyHgg=True' )
    if args.hzz: t2ws.extra_options.append('--PO isOnlyHZZ=True' )
    t2ws.add_theory_ktcg()
    t2ws.apply_reweighting()
    t2ws.add_theory_uncertainties()
    t2ws.add_scaling_ttH()
    t2ws.extra_options.append('--PO BRs_kappa_dependent=True')
    t2ws.tags.append('couplingdependentBRs')
    t2ws.run()


@flag_as_option
def projection_ktcg_t2ws_floatingBRs(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    t2ws = T2WSKTCGKB()
    t2ws.card = cards.ktcgkb[decay_channel]
    if args.hgg: t2ws.extra_options.append('--PO isOnlyHgg=True' )
    if args.hzz: t2ws.extra_options.append('--PO isOnlyHZZ=True' )
    t2ws.add_theory_ktcg()
    t2ws.apply_reweighting()
    t2ws.add_theory_uncertainties()
    t2ws.add_scaling_ttH()
    if args.combWithHbb:
        t2ws.constrain_bbZZ_ratio()
    t2ws.extra_options.append('--PO freely_floating_BRs=True')
    t2ws.tags.append('floatingBRs')
    t2ws.run()


