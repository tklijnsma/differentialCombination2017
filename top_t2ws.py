#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

from OptionHandler import flag_as_option

import LatestPaths
import LatestBinning

import differentials
import differentialutils

import logging
import copy
import re
import random
random.seed(1002)

from differentials.theory.theory_utils import FileFinder

from time import strftime
datestr = strftime('%b%d')

########################################
# Main
########################################

top_exp_binning = [ 0., 15., 30., 45., 80., 120., 200., 350., 600. ]


def base_t2ws(args,
        apply_theory_uncertainties=True,
        apply_reweighting=True,
        drop_last_bin=False,
        do_kappat_kappag=True,
        add_scaling_ttH=False,
        ):

    decay_channel = differentialutils.get_decay_channel_tag(args)
    t2ws = differentials.combine.t2ws.T2WS()
    t2ws.card = LatestPaths.card.pth_ggH[decay_channel]
    t2ws.model_file = 'physicsModels/CouplingModel.py'
    t2ws.model_name = 'couplingModel'

    if do_kappat_kappag:
        t2ws.name = decay_channel + '_Top'
    else:
        t2ws.name = decay_channel + '_TopCtCb'

    t2ws.extra_options.extend([
        '--PO linearTerms=False',
        '--PO splitggH=True',
        ])

    if args.hzz:
        t2ws.extra_options.append('--PO isOnlyHZZ=True' )
    if args.hgg:
        t2ws.extra_options.append('--PO isOnlyHgg=True' )

    obs = LatestBinning.obstuple_pth_ggH[decay_channel]
    # if drop_last_bin:
    #     obs.drop_last_bin()
    #     t2ws.tags.append('lastBinDropped')

    t2ws.extra_options.append(
        '--PO binBoundaries={0}'
        .format(','.join([str(b) for b in obs.binning]))
        )

    # Add the theory

    if do_kappat_kappag:
        logging.warning('Doing variations for kappa_t / kappa_g')
        coupling_variations = FileFinder(
            cb=1.0, muR=1.0, muF=1.0, Q=1.0, directory=LatestPaths.theory.top.filedir
            ).get()
        sm = [ v for v in coupling_variations if v.ct==1.0 and v.cg==0.0 ][0]
        coupling_variations.pop(coupling_variations.index(sm))

        t2ws.extra_options.append('--PO SM=[ct=1,cg=0,file={0}]'.format(sm.theory_file))
        for variation in coupling_variations:
            if variation.ct == 1.0 or variation.cg == 0.0: continue
            t2ws.extra_options.append(
                '--PO theory=[ct={0},cg={1},file={2}]'
                .format(variation.ct, variation.cg, variation.theory_file)
                )
    else:
        logging.warning('Doing variations for kappa_t / kappa_b')
        coupling_variations = FileFinder(
            cg=0.0, muR=1.0, muF=1.0, Q=1.0, directory=LatestPaths.theory.top.filedir
            ).get()
        sm = [ v for v in coupling_variations if v.ct==1.0 and v.cb==1.0 ][0]
        coupling_variations.pop(coupling_variations.index(sm))

        t2ws.extra_options.append('--PO SM=[ct=1,cb=1,file={0}]'.format(sm.theory_file))
        for variation in coupling_variations:
            if variation.ct == 1.0 or variation.cb == 1.0: continue
            t2ws.extra_options.append(
                '--PO theory=[ct={0},cb={1},file={2}]'
                .format(variation.ct, variation.cb, variation.theory_file)
                )

    if apply_reweighting:
        t2ws.extra_options.append(
            '--PO SMXS_of_input_ws={0}'.format(','.join([str(v) for v in obs.crosssection()]))
            )
        t2ws.tags.append('reweighted')
    if apply_theory_uncertainties:
        add_theory_uncertainties(t2ws)
    if add_scaling_ttH:
        t2ws.extra_options.append('--PO add_scaling_ttH=True')
        t2ws.tags.append('scalingttH')

    return t2ws

def add_theory_uncertainties(t2ws, uncorrelated=False):
    # Theory uncertainties
    correlation_matrix   = 'out/scalecorrelations_Mar06/corrMat_tophighpt.txt'
    theory_uncertainties = 'out/scalecorrelations_Mar06/errors_tophighpt.txt'
    if uncorrelated:
        correlation_matrix = LatestPaths.correlationMatrix_Yukawa_Uncorrelated # Uncorrelated
        t2ws.tags.append('uncorrelatedTheoryUnc')

    t2ws.extra_options.append('--PO correlationMatrix={0}'.format(correlation_matrix))
    t2ws.extra_options.append('--PO theoryUncertainties={0}'.format(theory_uncertainties))


#____________________________________________________________________
@flag_as_option
def t2ws_Top_scalingttH(args):
    # args = differentialutils.set_one_decay_channel(args, 'combWithHbb', asimov=True)
    t2ws = base_t2ws(args, add_scaling_ttH=True)
    # t2ws.card = LatestPaths.card.pth_ggH.combWithHbb
    t2ws.run()

@flag_as_option
def t2ws_Top_scalingttH_floatingBRs(args):
    # args = differentialutils.set_one_decay_channel(args, 'combWithHbb', asimov=True)
    t2ws = base_t2ws(args, add_scaling_ttH=True)
    # t2ws.card = LatestPaths.card.pth_ggH.combWithHbb
    t2ws.extra_options.append('--PO freely_floating_BRs=True')
    t2ws.tags.append('floatingBRs')
    t2ws.run()

@flag_as_option
def t2ws_Top_scalingttH_floatingBRs_constrainedbbZZ(args):
    # args = differentialutils.set_one_decay_channel(args, 'combWithHbb', asimov=True)
    t2ws = base_t2ws(args, add_scaling_ttH=True)
    # t2ws.card = LatestPaths.card.pth_ggH.combWithHbb
    t2ws.extra_options.append('--PO freely_floating_BRs=True')
    t2ws.extra_options.append('--PO constrain_ratio_bb_ZZ=True')
    t2ws.tags.append('floatingBRs')
    t2ws.tags.append('constrainedbbZZ')
    t2ws.run()

@flag_as_option
def t2ws_Top_scalingttH_couplingdependentBRs(args):
    # args = differentialutils.set_one_decay_channel(args, 'combWithHbb', asimov=True)
    t2ws = base_t2ws(args, add_scaling_ttH=True)
    # t2ws.card = LatestPaths.card.pth_ggH.combWithHbb
    t2ws.extra_options.append('--PO BRs_kappa_dependent=True')
    t2ws.tags.append('couplingdependentBRs')
    t2ws.run()

#____________________________________________________________________
@flag_as_option
def t2ws_ktcg_kbinterference(args):
    t2ws = base_t2ws(args, add_scaling_ttH=True)
    t2ws.extra_options.append('--PO BRs_kappa_dependent=True')
    t2ws.tags.append('kbinterference')

    # Remove previous options
    remove = lambda option: (
                option.startswith('--PO SM=')
                or option.startswith('--PO theory=')
                or option.startswith('--PO linearTerms=')
                )
    t2ws.extra_options = [ option for option in t2ws.extra_options if not(remove(option)) ]

    # Add the linear terms
    t2ws.extra_options.append('--PO linearTerms=True')

    # Add theory
    coupling_variations = FileFinder(
        cb=1.0, muR=1.0, muF=1.0, Q=1.0, directory=LatestPaths.theory.top.filedir
        ).get()
    sm = [ v for v in coupling_variations if v.ct==1.0 and v.cg==0.0 ][0]
    coupling_variations.pop(coupling_variations.index(sm))

    t2ws.extra_options.append('--PO SM=[ct=1,cg=0,file={0}]'.format(sm.theory_file))
    for variation in coupling_variations:
        if variation.ct == 1.0 or variation.cg == 0.0: continue
        t2ws.extra_options.append(
            '--PO theory=[ct={0},cg={1},file={2}]'
            .format(variation.ct, variation.cg, variation.theory_file)
            )

    t2ws.run()

#____________________________________________________________________
@flag_as_option
def t2ws_TopCtCb_scalingbbHttH(args):
    # args = differentialutils.set_one_decay_channel(args, 'combWithHbb', asimov=True)
    t2ws = base_t2ws(args, do_kappat_kappag=False)
    t2ws.extra_options.append('--PO add_scaling_ttH=True')
    t2ws.extra_options.append('--PO add_scaling_bbH=True')
    t2ws.tags.append('scalingbbHttH')
    t2ws.run()

@flag_as_option
def t2ws_TopCtCb_scalingbbHttH_couplingdependentBRs(args):
    # args = differentialutils.set_one_decay_channel(args, 'combWithHbb', asimov=True)
    t2ws = base_t2ws(args, do_kappat_kappag=False)
    t2ws.extra_options.append('--PO add_scaling_ttH=True')
    t2ws.extra_options.append('--PO add_scaling_bbH=True')
    t2ws.tags.append('scalingbbHttH')
    t2ws.extra_options.append('--PO BRs_kappa_dependent=True')
    t2ws.tags.append('couplingdependentBRs')
    t2ws.run()

@flag_as_option
def t2ws_TopCtCb_scalingbbHttH_floatingBRs(args):
    # args = differentialutils.set_one_decay_channel(args, 'combWithHbb', asimov=True)
    t2ws = base_t2ws(args, do_kappat_kappag=False)
    t2ws.extra_options.append('--PO add_scaling_ttH=True')
    t2ws.extra_options.append('--PO add_scaling_bbH=True')
    t2ws.tags.append('scalingbbHttH')
    t2ws.extra_options.append('--PO freely_floating_BRs=True')
    t2ws.tags.append('floatingBRs')
    t2ws.run()

@flag_as_option
def t2ws_TopCtCb_scalingbbHttH_floatingBRs_constrainedbbZZ(args):
    # args = differentialutils.set_one_decay_channel(args, 'combWithHbb', asimov=True)
    t2ws = base_t2ws(args, do_kappat_kappag=False)
    t2ws.extra_options.append('--PO add_scaling_ttH=True')
    t2ws.extra_options.append('--PO add_scaling_bbH=True')
    t2ws.tags.append('scalingbbHttH')
    t2ws.extra_options.append('--PO freely_floating_BRs=True')
    t2ws.extra_options.append('--PO constrain_ratio_bb_ZZ=True')
    t2ws.tags.append('floatingBRs')
    t2ws.tags.append('constrainedbbZZ')
    t2ws.run()


#____________________________________________________________________
@flag_as_option
def t2ws_Top_nominal(args):
    t2ws = base_t2ws(args)
    t2ws.card = LatestPaths.card.top.nominal[differentialutils.get_decay_channel_tag(args)]
    t2ws.tags.append('noBinsDropped')
    t2ws.run()

@flag_as_option
def t2ws_Top_lumiScale(args):
    t2ws = base_t2ws(args)
    t2ws.card = LatestPaths.card.top.nominal[differentialutils.get_decay_channel_tag(args)]
    t2ws.tags.append('noBinsDropped')
    t2ws.tags.append('lumiScale')
    t2ws.extra_options.append('--PO lumiScale=True')
    t2ws.run()

@flag_as_option
def t2ws_TopCtCb_nominal(args):
    t2ws = base_t2ws(args, do_kappat_kappag=False)
    t2ws.card = LatestPaths.card.top.nominal[differentialutils.get_decay_channel_tag(args)]
    t2ws.tags.append('noBinsDropped')
    t2ws.run()

@flag_as_option
def t2ws_TopCtCb_lumiScale(args):
    t2ws = base_t2ws(args, do_kappat_kappag=False)
    t2ws.card = LatestPaths.card.top.nominal[differentialutils.get_decay_channel_tag(args)]
    t2ws.tags.append('noBinsDropped')
    t2ws.tags.append('lumiScale')
    t2ws.extra_options.append('--PO lumiScale=True')
    t2ws.run()


@flag_as_option
def t2ws_Top_lastBinDroppedHgg(args):
    t2ws = base_t2ws(args)
    if args.combination:
        t2ws.card = 'suppliedInput/combination_pth_ggH_lastBinDroppedHgg_Mar07.txt'
    elif args.combWithHbb:
        t2ws.card = 'suppliedInput/combWithHbb_pth_ggH_lastBinDroppedHgg_Mar07.txt'
    else:
        raise NotImplementedError('Use --combination or --combWithHbb')
    t2ws.tags.append('lastBinDroppedHgg')
    t2ws.run()


@flag_as_option
def t2ws_Top_last2BinsDropped(args):
    t2ws = base_t2ws(args)
    if args.combWithHbb:
        t2ws.card = LatestPaths.card.top.combWithHbb_last2BinsDropped
    if args.combination:
        t2ws.card = LatestPaths.card.top.combination_last2BinsDropped
    else:
        raise NotImplementedError('Use --combination or --combWithHbb')
    t2ws.tags.append('last2BinsDropped')
    t2ws.run()

@flag_as_option
def t2ws_Top_BRcouplingDependency(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('BRcouplingDependency')
    t2ws.extra_options.append('--PO FitBR=True')
    # do_BR_uncertainties = False
    # if do_BR_uncertainties:
    #     t2ws.extra_options.append('--PO DoBRUncertainties=True')
    #     t2ws.tags.append('withBRUnc')
    t2ws.run()






#____________________________________________________________________

@flag_as_option
def t2ws_Top_unreweighted(args):
    t2ws = base_t2ws(args, apply_reweighting=False)
    t2ws.tags.append('nominal')
    t2ws.run()

@flag_as_option
def t2ws_Top_noTheoryUnc(args):
    t2ws = base_t2ws(args, apply_theory_uncertainties=False)
    t2ws.tags.append('noTheoryUnc')
    t2ws.run()

@flag_as_option
def t2ws_Top_uncorrelatedTheoryUnc(args):
    t2ws = base_t2ws(args, apply_theory_uncertainties=False)
    add_theory_uncertainties(t2ws, uncorrelated=True)
    t2ws.tags.append('uncorrelatedTheoryUnc')
    t2ws.run()


# @flag_as_option
# def t2ws_Top_lumiScale(args):
#     t2ws = base_t2ws(args)
#     t2ws.tags.append('lumiScale')
#     t2ws.extra_options.append('--PO lumiScale=True')
#     t2ws.run()

@flag_as_option
def t2ws_Top_profiledTotalXS(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('profiledTotalXS')
    t2ws.extra_options.append('--PO ProfileTotalXS=True')
    t2ws.run()


#____________________________________________________________________
# Not yet repeated

@flag_as_option
def t2ws_Top_fitRatioOfBRs(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('fitRatioOfBRs')
    t2ws.extra_options.append('--PO FitRatioOfBRs=True')
    t2ws.run()

@flag_as_option
def t2ws_Top_fitOnlyNormalization(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('fitOnlyNormalization')
    add_theory_uncertainties(t2ws)
    t2ws.extra_options.append('--PO FitOnlyNormalization=True')
    t2ws.run()



