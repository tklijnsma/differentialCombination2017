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

yukawa_exp_binning = [ 0., 15., 30., 45., 85., 125. ]

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
    t2ws = differentials.combine.t2ws.T2WS(get_nominal_card_Yukawa(args))
    t2ws.model_file = 'physicsModels/CouplingModel.py'
    t2ws.model_name = 'couplingModel'

    decay_channel = differentials.core.get_decay_channel_tag(args, allow_default=True)
    t2ws.name = decay_channel + '_Yukawa'

    t2ws.extra_options.extend([
        '--PO linearTerms=True',
        '--PO splitggH=True',
        ])

    if args.hzz:
        t2ws.extra_options.append('--PO isOnlyHZZ=True' )
    if args.hgg:
        t2ws.extra_options.append('--PO isOnlyHgg=True' )

    obs = LatestBinning.obs_pth_ggH
    obs.drop_bins_up_to_value(yukawa_exp_binning[-1]+1.0)

    t2ws.extra_options.append(
        '--PO binBoundaries={0}'
        .format(','.join([str(b) for b in obs.binning]))
        )

    if apply_reweighting:
        crosssections = obs.crosssection()
        crosssections_str = ','.join([str(v) for v in crosssections])
        t2ws.extra_options.append('--PO ReweightCrossSections={0}'.format(crosssections_str))
        t2ws.tags.append('reweighted')
    if apply_theory_uncertainties:
        add_theory_uncertainties(t2ws)
    return t2ws


def add_theory_uncertainties(t2ws, uncorrelated=False):
    TheoryFileInterface.SetFileFinderDir(LatestPaths.derivedTheoryFiles_YukawaSummed)

    t2ws.extra_options.append(
        '--PO SM=[kappab=1,kappac=1,file={0}]'.format(
            TheoryFileInterface.FileFinder(kappab=1, kappac=1, muR=1, muF=1, Q=1, expectOneFile=True )
            )
        )
    possible_theories = []
    for kappab in [ -2, -1, 0, 1, 2 ]:
        for kappac in [ -10, -5, 0, 1, 5, 10 ]:
            if (kappab == 1 and kappac == 1 ) or (kappab == 0 and kappac == 0 ): continue
            else:
                possible_theories.append(
                    '--PO theory=[kappab={0},kappac={1},file={2}]'.format(
                        kappab, kappac,
                        TheoryFileInterface.FileFinder(kappab=kappab, kappac=kappac, muR=1, muF=1, Q=1, expectOneFile=True )
                        )
                    )
    t2ws.extra_options.extend(random.sample(possible_theories, 6))

    # Theory uncertainties
    correlationMatrix   = LatestPaths.correlationMatrix_Yukawa
    theoryUncertainties = LatestPaths.theoryUncertainties_Yukawa
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
def all_t2ws_Yukawa(args_original):
    args = copy.deepcopy(args_original)
    for decay_channel in ['hgg', 'hzz', 'combination']:
        set_decay_channel(args, decay_channel)
        # t2ws_Yukawa_nominal(args)
        try_call_function_with_args(t2ws_Yukawa_nominal, args)

    set_decay_channel(args, 'combination')

    fns = [
        t2ws_Yukawa_noTheoryUnc,
        t2ws_Yukawa_uncorrelatedTheoryUnc,
        t2ws_Yukawa_lumiScale,
        t2ws_Yukawa_BRcouplingDependency,
        t2ws_Yukawa_profiledTotalXS,
        t2ws_Yukawa_fitRatioOfBRs,
        t2ws_Yukawa_fitOnlyNormalization,
        ]
    for fn in fns:
        try_call_function_with_args(fn, args)

    # t2ws_Yukawa_noTheoryUnc(args)
    # t2ws_Yukawa_uncorrelatedTheoryUnc(args)
    # t2ws_Yukawa_lumiScale(args)
    # t2ws_Yukawa_BRcouplingDependency(args)
    # t2ws_Yukawa_profiledTotalXS(args)
    # t2ws_Yukawa_fitRatioOfBRs(args)
    # t2ws_Yukawa_fitOnlyNormalization(args)

#____________________________________________________________________
@flag_as_option
def t2ws_Yukawa_nominal(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('nominal')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_noTheoryUnc(args):
    t2ws = base_t2ws(args, apply_theory_uncertainties=False)
    t2ws.tags.append('noTheoryUnc')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_uncorrelatedTheoryUnc(args):
    t2ws = base_t2ws(args, apply_theory_uncertainties=False)
    add_theory_uncertainties(t2ws, uncorrelated=True)
    t2ws.tags.append('uncorrelatedTheoryUnc')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_lumiScale(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('lumiScale')
    t2ws.extra_options.append('--PO lumiScale=True')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_BRcouplingDependency(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('BRcouplingDependency')
    t2ws.extra_options.append('--PO FitBR=True')
    do_BR_uncertainties = False
    if do_BR_uncertainties:
        t2ws.extra_options.append('--PO DoBRUncertainties=True')
        t2ws.tags.append('withBRUnc')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_profiledTotalXS(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('profiledTotalXS')
    t2ws.extra_options.append('--PO ProfileTotalXS=True')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_fitRatioOfBRs(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('fitRatioOfBRs')
    t2ws.extra_options.append('--PO FitRatioOfBRs=True')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_fitOnlyNormalization(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('fitOnlyNormalization')
    add_theory_uncertainties(t2ws)
    t2ws.extra_options.append('--PO FitOnlyNormalization=True')
    t2ws.run()



