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

import logging
import copy
import re
import random
random.seed(1002)

import sys
sys.path.append('src')
import TheoryFileInterface

from time import strftime
datestr = strftime('%b%d')

########################################
# Main
########################################

top_exp_binning = [ 0., 15., 30., 45., 80., 120., 200., 350., 600. ]

def get_nominal_card_Yukawa(args):
    datacard = LatestPaths.card_combined_ggHxH_PTH
    if args.hgg:
        datacard = LatestPaths.card_hgg_ggHxH_PTH
    if args.hzz:
        datacard = LatestPaths.card_hzz_ggHxH_PTH
    if args.hbb:
        raise NotImplementedError('This makes no sense, as Yukawa only goes up to 125!')
        datacard = LatestPaths.card_hbb_ggHxH_PTH
    if args.combWithHbb:
        raise NotImplementedError('This makes no sense, as Yukawa only goes up to 125!')
        datacard = LatestPaths.card_combinedWithHbb_ggHxH_PTH
    return datacard

def base_t2ws(args, apply_theory_uncertainties=True, apply_reweighting=True):
    decay_channel = differentials.core.get_decay_channel_tag(args)
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
    TheoryFileInterface.SetFileFinderDir(LatestPaths.theory.top.filedir)

    t2ws.extra_options.append(
        '--PO SM=[ct=1,cg=0,file={0}]'.format(
            TheoryFileInterface.FileFinder(ct=1, cg=0, muR=1, muF=1, Q=1, expectOneFile=True)
            )
        )

    theoryFiles = TheoryFileInterface.FileFinder(
        ct='*', cg='*', cb=1, muR=1, muF=1, Q=1, filter='ct_1_cg_0',
        )            
    possible_theories = []
    for theoryFile in theoryFiles:
        ct = differentials.core.str_to_float(re.search(r'ct_([\dmp]+)', theoryFile ).group(1))
        cg = differentials.core.str_to_float(re.search(r'cg_([\dmp]+)', theoryFile ).group(1))
        possible_theories.append(
            '--PO theory=[ct={0},cg={1},file={2}]'.format(
                ct, cg, theoryFile
                )                
            )
    t2ws.extra_options.extend(possible_theories)

    # Theory uncertainties
    correlationMatrix   = LatestPaths.theory.top.correlation_matrix
    theoryUncertainties = LatestPaths.theory.top.uncertainties
    if uncorrelated:
        correlationMatrix = LatestPaths.correlationMatrix_Yukawa_Uncorrelated # Uncorrelated
        t2ws.tags.append('uncorrelatedTheoryUnc')

    t2ws.extra_options.append('--PO correlationMatrix={0}'.format(correlationMatrix))
    t2ws.extra_options.append('--PO theoryUncertainties={0}'.format(theoryUncertainties))
    # t2ws.tags.append('withTheoryUncertainties') # Pretty much the default anyway

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
    t2ws.tags.append('nominal')
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
def t2ws_Top_profiledTotalXS(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('profiledTotalXS')
    t2ws.extra_options.append('--PO ProfileTotalXS=True')
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



