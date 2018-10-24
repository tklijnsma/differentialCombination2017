#!/usr/bin/env python
"""
Thomas Klijnsma
"""

#____________________________________________________________________
# Imports

import logging
import os, sys, re, copy
from OptionHandler import flag_as_option, flag_as_parser_options
import differentials
import differentialutils
import LatestBinning
import LatestPaths

import random
random.seed(1002)
from differentials.theory.theory_utils import FileFinder

import projections_paths as paths

#____________________________________________________________________
@flag_as_parser_options
def add_projection_t2ws_options(parser):
    parser.add_argument('--scenario2', action='store_true', help='boolean')

#____________________________________________________________________
s2_scalingfunctions = [
    "--X-nuisance-group-function 'pBTag' '1.0'",
    "--X-nuisance-group-function 'pBTagStat' 'expr::scaleBTagStat(\"1/sqrt(@0)\",lumiscale[1])'",
    "--X-nuisance-group-function 'pEleID' 'expr::scaleEleID(\"max(0.5,1/sqrt(@0))\",lumiscale[1])'",
    "--X-nuisance-group-function 'pMuonID' 'expr::scaleMuonID(\"max(0.5,1/sqrt(@0))\",lumiscale[1])'",
    "--X-nuisance-group-function 'pPhotonID' 'expr::scalePhotonID(\"max(0.5,1/sqrt(@0))\",lumiscale[1])'",
    "--X-nuisance-group-function 'pTauID' '1.0'",
    "--X-nuisance-group-function 'pScaleJ' 'expr::scaleScaleJ(\"max(0.5,1/sqrt(@0))\",lumiscale[1])'",
    "--X-nuisance-group-function 'pScaleJAbs' 'expr::scaleScaleJAbs(\"max(0.3,1/sqrt(@0))\",lumiscale[1])'",
    "--X-nuisance-group-function 'pScaleJFlav' 'expr::scaleScaleJFlav(\"max(0.5,1/sqrt(@0))\",lumiscale[1])'",
    "--X-nuisance-group-function 'pScaleJPileup' '1.0'",
    "--X-nuisance-group-function 'pScaleJRel' 'expr::scaleScaleJRel(\"max(0.2,1/sqrt(@0))\",lumiscale[1])'",
    "--X-nuisance-group-function 'pScaleJTime' 'expr::scaleScaleJTime(\"1/sqrt(@0)\",lumiscale[1])'",
    "--X-nuisance-group-function 'pScaleJMethod' 'expr::scaleScaleJMethod(\"1/sqrt(@0)\",lumiscale[1])'",
    "--X-nuisance-group-function 'pResJ' 'expr::scaleResJ(\"max(0.5,1/sqrt(@0))\",lumiscale[1])'",
    "--X-nuisance-group-function 'pScaleMet' 'expr::scaleScaleMet(\"max(0.5,1/sqrt(@0))\",lumiscale[1])'",
    "--X-nuisance-group-function 'pLumi' '0.4'",
    "--X-nuisance-group-function 'pOther' 'expr::scaleOther(\"1/sqrt(@0)\",lumiscale[1])'",
    "--X-nuisance-group-function 'sigTheory' '0.5'",
    "--X-nuisance-group-function 'bkgTheory' '0.5'",
    ]

#____________________________________________________________________
cards_nongrouped = differentials.core.AttrDict()
cards_nongrouped.hgg = paths.pt_cards_nongrouped.hgg
cards_nongrouped.hzz = paths.pt_cards_nongrouped.hzz
cards_nongrouped.hbb = paths.pt_cards_nongrouped.hbb
cards_nongrouped.combWithHbb = paths.pt_cards_nongrouped.combWithHbb

cards = differentials.core.AttrDict()
cards.hgg         = paths.pt_cards.hgg
cards.hzz         = paths.pt_cards.hzz
cards.hbb         = paths.pt_cards.hbb
cards.combWithHbb = paths.pt_cards.combWithHbb

@flag_as_option
def projection_pth_smH_t2ws(args):
    differentials.combine.combine_utils.set_global_outdir('projections')
    decay_channel = differentialutils.get_decay_channel_tag(args)
    t2ws = differentials.combine.t2ws.T2WS(
        cards[decay_channel],
        name = 'ws_pth_smH_' + decay_channel
        )

    if args.scenario2:
        t2ws.name += '_s2'
        t2ws.extra_options.extend(s2_scalingfunctions)
    else:
        t2ws.name += '_s1'

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
@flag_as_option
def projection_pth_smH_t2ws_GT200(args):
    differentials.combine.combine_utils.set_global_outdir('projections')
    decay_channel = differentialutils.get_decay_channel_tag(args)
    t2ws = differentials.combine.t2ws.T2WS(
        cards[decay_channel],
        name = 'ws_pth_smH_' + decay_channel + '_GT200'
        )

    if args.scenario2:
        t2ws.name += '_s2'
        t2ws.extra_options.extend(s2_scalingfunctions)
    else:
        t2ws.name += '_s1'

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
        for i, map_option in enumerate(t2ws.map_options):
            poi = re.search(r':(r_.*?)\[', map_option).group(1)
            left = differentials.processinterpreter.Process(poi).left
            if left >= 200.:
                new_poi = 'r_smH_PTH_GT200'
                t2ws.map_options[i] = map_option.replace(poi, new_poi)

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

    def apply_s2_scaling(self):
        self.tags.append('scenario2')
        self.extra_options.extend(s2_scalingfunctions)

cards.yukawa = differentials.core.AttrDict()
cards.yukawa.hgg         = 'suppliedInput/projection_yukawa_hgg_Jul09.txt'
cards.yukawa.hzz         = 'suppliedInput/projection_yukawa_hzz_Jul09.txt'
cards.yukawa.combination = 'suppliedInput/projection_yukawa_combination_Jul09.txt'

cards.yukawa.s2grouping = differentials.core.AttrDict()
cards.yukawa.s2grouping.hgg = 'suppliedInput/projection_yukawa_hgg_s2groups_Jul17.txt'
cards.yukawa.s2grouping.hzz = 'suppliedInput/projection_yukawa_hzz_s2groups_Jul17.txt'
cards.yukawa.s2grouping.combination = 'suppliedInput/projection_yukawa_combination_s2groups_Jul17.txt'

@flag_as_option
def projection_kbkc_t2ws_couplingdependentBRs(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    t2ws = T2WSKBKC()
    t2ws.card = cards.yukawa.s2grouping[decay_channel]
    if args.hgg: t2ws.extra_options.append('--PO isOnlyHgg=True' )
    if args.hzz: t2ws.extra_options.append('--PO isOnlyHZZ=True' )
    t2ws.add_theory()
    t2ws.apply_reweighting()
    t2ws.add_theory_uncertainties()
    if args.scenario2: t2ws.apply_s2_scaling()
    t2ws.extra_options.append('--PO BRs_kappa_dependent=True')
    t2ws.tags.append('couplingdependentBRs')
    t2ws.run()

@flag_as_option
def projection_kbkc_t2ws_floatingBRs(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    t2ws = T2WSKBKC()
    t2ws.card = cards.yukawa.s2grouping[decay_channel]
    if args.hgg: t2ws.extra_options.append('--PO isOnlyHgg=True' )
    if args.hzz: t2ws.extra_options.append('--PO isOnlyHZZ=True' )
    t2ws.add_theory()
    t2ws.apply_reweighting()
    t2ws.add_theory_uncertainties()
    if args.scenario2: t2ws.apply_s2_scaling()
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

    def apply_s2_scaling(self):
        self.tags.append('scenario2')
        self.extra_options.extend(s2_scalingfunctions)

cards.ktcgkb = differentials.core.AttrDict()
cards.ktcgkb.hgg         = 'suppliedInput/fromVittorio/pT_newBins_Feb28/projection_ktcgkb_hgg_Jul10.txt'
cards.ktcgkb.hzz         = 'suppliedInput/fromDavid/PTH_Jan24_newBinning/ggH/projection_ktcgkb_hzz_Jul10.txt'
cards.ktcgkb.combWithHbb = 'suppliedInput/projection_ktcgkb_combWithHbb_Jul10.txt'

cards.ktcgkb.s2grouping = differentials.core.AttrDict()
cards.ktcgkb.s2grouping.hgg         = 'suppliedInput/fromVittorio/pT_newBins_Feb28/projection_ktcgkb_hgg_s2groups_Jul18.txt'
cards.ktcgkb.s2grouping.hzz         = 'suppliedInput/fromDavid/PTH_Jan24_newBinning/ggH/projection_ktcgkb_hzz_s2groups_Jul18.txt'
cards.ktcgkb.s2grouping.combWithHbb = 'suppliedInput/projection_ktcgkb_combWithHbb_s2groups_Jul18.txt'

@flag_as_option
def projection_ktcg_t2ws_couplingdependentBRs(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    t2ws = T2WSKTCGKB()
    t2ws.card = cards.ktcgkb.s2grouping[decay_channel]
    if args.hgg: t2ws.extra_options.append('--PO isOnlyHgg=True' )
    if args.hzz: t2ws.extra_options.append('--PO isOnlyHZZ=True' )
    t2ws.add_theory_ktcg()
    t2ws.apply_reweighting()
    t2ws.add_theory_uncertainties()
    t2ws.add_scaling_ttH()
    if args.scenario2: t2ws.apply_s2_scaling()
    t2ws.extra_options.append('--PO BRs_kappa_dependent=True')
    t2ws.tags.append('couplingdependentBRs')
    t2ws.run()

@flag_as_option
def projection_ktcg_t2ws_floatingBRs(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    t2ws = T2WSKTCGKB()
    t2ws.card = cards.ktcgkb.s2grouping[decay_channel]
    if args.hgg: t2ws.extra_options.append('--PO isOnlyHgg=True' )
    if args.hzz: t2ws.extra_options.append('--PO isOnlyHZZ=True' )
    t2ws.add_theory_ktcg()
    t2ws.apply_reweighting()
    t2ws.add_theory_uncertainties()
    t2ws.add_scaling_ttH()
    # if args.combWithHbb:
    #     t2ws.constrain_bbZZ_ratio()
    if args.scenario2: t2ws.apply_s2_scaling()
    t2ws.extra_options.append('--PO freely_floating_BRs=True')
    t2ws.tags.append('floatingBRs')
    t2ws.run()


