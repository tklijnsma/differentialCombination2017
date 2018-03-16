#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

from OptionHandler import flag_as_option

import LatestPaths
import LatestPathsGetters
import LatestBinning

import differentials
import differentialutils

import logging
import copy
import re
import random
random.seed(1002)

import sys
sys.path.append('src')
import TheoryFileInterface

from differentials.theory.theory_utils import FileFinder

from time import strftime
datestr = strftime('%b%d')

########################################
# Main
########################################

top_exp_binning = [ 0., 15., 30., 45., 80., 120., 200., 350., 600. ]


def base_t2ws(args, apply_theory_uncertainties=True, apply_reweighting=True, drop_last_bin=False):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    card = LatestPaths.card.pth_ggH[decay_channel]
    t2ws = differentials.combine.t2ws.T2WS(card)
    t2ws.model_file = 'physicsModels/CouplingModel.py'
    t2ws.model_name = 'couplingModel'

    t2ws.name = decay_channel + '_Top'
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

    if apply_reweighting:
        t2ws.extra_options.append(
            '--PO ReweightCrossSections={0}'.format(','.join([str(v) for v in obs.crosssection()]))
            )
        t2ws.tags.append('reweighted')
    if apply_theory_uncertainties:
        add_theory_uncertainties(t2ws)
    return t2ws

def add_theory_uncertainties(t2ws, uncorrelated=False):
    coupling_variations = FileFinder(
        cb=1.0, muR=1.0, muF=1.0, Q=1.0, directory='out/theories_Mar05_tophighpt'
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

    # Theory uncertainties
    correlation_matrix   = 'out/scalecorrelations_Mar06/corrMat_tophighpt.txt'
    theory_uncertainties = 'out/scalecorrelations_Mar06/errors_tophighpt.txt'
    if uncorrelated:
        correlation_matrix = LatestPaths.correlationMatrix_Yukawa_Uncorrelated # Uncorrelated
        t2ws.tags.append('uncorrelatedTheoryUnc')

    t2ws.extra_options.append('--PO correlationMatrix={0}'.format(correlation_matrix))
    t2ws.extra_options.append('--PO theoryUncertainties={0}'.format(theory_uncertainties))


#____________________________________________________________________
def set_decay_channel(args, given_channel):
    for decay_channel in ['hgg', 'hzz', 'combination', 'hbb', 'combWithHbb']:
        setattr(args, decay_channel, False)
    setattr(args, given_channel, True)


import traceback
def try_call_function_with_args(fn, args):
    try:
        fn(args)
    except Exception as exc:
        print traceback.format_exc()
        print exc
        pass

@flag_as_option
def all_t2ws_Top(args_original):
    args = copy.deepcopy(args_original)
    for decay_channel in ['hgg', 'hzz', 'combination', 'combWithHbb']:
        logging.info('\n\n' + '*'*80)
        logging.info('Running t2ws_Top_nominal on decay channel {0}'.format(decay_channel))
        set_decay_channel(args, decay_channel)
        try_call_function_with_args(t2ws_Top_nominal, args)

    # set_decay_channel(args, 'combination')
    # fns = [
    #     t2ws_Top_noTheoryUnc,
    #     t2ws_Top_uncorrelatedTheoryUnc,
    #     t2ws_Top_lumiScale,
    #     t2ws_Top_BRcouplingDependency,
    #     t2ws_Top_profiledTotalXS,
    #     t2ws_Top_fitRatioOfBRs,
    #     t2ws_Top_fitOnlyNormalization,
    #     ]
    # for fn in fns:
    #     try_call_function_with_args(fn, args)


#____________________________________________________________________
@flag_as_option
def t2ws_Top_nominal(args):
    t2ws = base_t2ws(args)
    t2ws.card = LatestPaths.card.top.nominal[differentialutils.get_decay_channel_tag(args)]
    t2ws.tags.append('noBinsDropped')
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


@flag_as_option
def t2ws_Top_lumiScale(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('lumiScale')
    t2ws.extra_options.append('--PO lumiScale=True')
    t2ws.run()

@flag_as_option
def t2ws_Top_profiledTotalXS(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('profiledTotalXS')
    t2ws.extra_options.append('--PO ProfileTotalXS=True')
    t2ws.run()


@flag_as_option
def t2ws_Top_BRcouplingDependency(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('BRcouplingDependency')
    t2ws.extra_options.append('--PO FitBR=True')
    do_BR_uncertainties = False
    if do_BR_uncertainties:
        t2ws.extra_options.append('--PO DoBRUncertainties=True')
        t2ws.tags.append('withBRUnc')
    t2ws.run()

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



