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
import random
random.seed(1002)

from differentials.theory.theory_utils import FileFinder

from time import strftime
datestr = strftime('%b%d')

########################################
# Main
########################################

yukawa_exp_binning = [ 0., 15., 30., 45., 80., 120. ]

@flag_as_option
def make_theory_files_Yukawa(args):
    interp = differentials.theory.kappabkappac_interpreter.KappabKappacInterpreter()
    interp.dump_gluon_induced()
    interp.dump_quark_induced()
    interp.dump_quark_induced_scaled()
    interp.dump_summed_quark_gluon_induced()


def base_t2ws(args, apply_theory_uncertainties=True, apply_reweighting=True):
    t2ws = differentials.combine.t2ws.T2WS()

    decay_channel = differentialutils.get_decay_channel_tag(args, allow_default=True)
    t2ws.card = LatestPaths.card.yukawa[decay_channel]
    t2ws.model_file = 'physicsModels/CouplingModel.py'
    t2ws.model_name = 'couplingModel'
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

    add_theory(t2ws)

    if apply_reweighting:
        crosssections = obs.crosssection()
        crosssections_str = ','.join([str(v) for v in crosssections])
        t2ws.extra_options.append('--PO SMXS_of_input_ws={0}'.format(crosssections_str))
        t2ws.tags.append('reweighted')
    else:
        t2ws.tags.append('unreweighted')

    if apply_theory_uncertainties:
        add_theory_uncertainties(t2ws)
    return t2ws

def base_t2ws_fitOnlyTotalXS(args, apply_theory_uncertainties=True, apply_reweighting=False):
    t2ws = differentials.combine.t2ws.T2WS()

    decay_channel = differentialutils.get_decay_channel_tag(args, allow_default=True)
    t2ws.card = LatestPaths.card.inclusive[decay_channel]
    t2ws.model_file = 'physicsModels/CouplingModel.py'
    t2ws.model_name = 'couplingModel'
    t2ws.name = decay_channel + '_Yukawa'

    t2ws.extra_options.extend([
        '--PO linearTerms=True',
        '--PO splitggH=True',
        '--PO FitOnlyTotalXS=True'
        ])

    if args.hzz:
        t2ws.extra_options.append('--PO isOnlyHZZ=True' )
    if args.hgg:
        t2ws.extra_options.append('--PO isOnlyHgg=True' )

    add_theory(t2ws)

    if apply_reweighting:
        logging.warning('Not sure if applying the reweighting is a good thing for the incl.')
        obs = LatestBinning.obs_pth_ggH
        obs.drop_bins_up_to_value(yukawa_exp_binning[-1]+1.0)
        inc_smxs = obs.inclusive_crosssection()
        t2ws.extra_options.append('--PO SMXS_of_input_ws={0}'.format(inc_smxs))
        t2ws.tags.append('reweighted')
    else:
        t2ws.tags.append('unreweighted')

    if apply_theory_uncertainties:
        # 0.0795 based on theory_uncertainty_on_inclusive_xs_yukawa (crosschecks.py)
        t2ws.extra_options.append('--PO inc_xs_uncertainty={0}'.format(0.079528392665))

    return t2ws


def add_theory(t2ws):
    coupling_variations = FileFinder(
        muR=1.0, muF=1.0, Q=1.0, directory=LatestPaths.theory.yukawa.filedir
        ).get()

    sm = [ v for v in coupling_variations if v.kappab==1.0 and v.kappac==1.0 ][0]
    coupling_variations.pop(coupling_variations.index(sm))

    t2ws.extra_options.append('--PO SM=[kappab=1.0,kappac=1.0,file={0}]'.format(sm.theory_file))
    for variation in random.sample(coupling_variations, 7):
        if variation.kappab == 0.0 and variation.kappac == 0.0: continue
        t2ws.extra_options.append(
            '--PO theory=[kappab={0},kappac={1},file={2}]'
            .format(variation.kappab, variation.kappac, variation.theory_file)
            )


def add_theory_uncertainties(t2ws, uncorrelated=False):
    # Theory uncertainties
    correlation_matrix   = LatestPaths.theory.yukawa.correlation_matrix
    theory_uncertainties = LatestPaths.theory.yukawa.uncertainties
    if uncorrelated:
        correlation_matrix = LatestPaths.theory.yukawa.correlation_matrix_uncorrelated
        t2ws.tags.append('uncorrelatedTheoryUnc')

    t2ws.extra_options.append('--PO correlationMatrix={0}'.format(correlation_matrix))
    t2ws.extra_options.append('--PO theoryUncertainties={0}'.format(theory_uncertainties))


#____________________________________________________________________
@flag_as_option
def t2ws_Yukawa_NOTscalingbbH(args):
    t2ws = base_t2ws(args)
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_NOTscalingbbH_floatingBRs(args):
    t2ws = base_t2ws(args)
    t2ws.extra_options.append('--PO freely_floating_BRs=True')
    t2ws.tags.append('floatingBRs')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_NOTscalingbbH_couplingdependentBRs(args):
    t2ws = base_t2ws(args)
    t2ws.extra_options.append('--PO BRs_kappa_dependent=True')
    t2ws.tags.append('couplingdependentBRs')
    t2ws.run()



#____________________________________________________________________
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
        t2ws_Yukawa_profiledTotalXS,
        # 
        # t2ws_Yukawa_BRcouplingDependency,
        # t2ws_Yukawa_fitRatioOfBRs,
        # t2ws_Yukawa_fitOnlyNormalization,
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
def t2ws_Yukawa_withBRuncertainties(args):
    t2ws = base_t2ws(args)
    t2ws.extra_options.append('--PO do_BR_uncertainties=True')
    t2ws.tags.append('withBRuncertainties')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_noReweighting(args):
    t2ws = base_t2ws(args, apply_reweighting=False)
    t2ws.tags.append('nominal')
    t2ws.run()

#____________________________________________________________________
# bbH scaling is NOT necessary for kb/kc!

# @flag_as_option
# def t2ws_Yukawa_scalingbbH(args):
#     # Basically nominal, but force asimov and hzz
#     # args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
#     t2ws = base_t2ws(args)
#     t2ws.extra_options.append('--PO add_scaling_bbH=True')
#     t2ws.tags.append('scalingbbH')
#     t2ws.run()

# @flag_as_option
# def t2ws_Yukawa_scalingbbH_floatingBRs(args):
#     # args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
#     t2ws = base_t2ws(args)
#     t2ws.extra_options.append('--PO add_scaling_bbH=True')
#     t2ws.tags.append('scalingbbH')
#     t2ws.extra_options.append('--PO freely_floating_BRs=True')
#     t2ws.tags.append('floatingBRs')
#     t2ws.run()

# @flag_as_option
# def t2ws_Yukawa_scalingbbH_couplingdependentBRs(args):
#     # Basically nominal, but force asimov and hzz
#     # args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
#     t2ws = base_t2ws(args)
#     t2ws.extra_options.append('--PO add_scaling_bbH=True')
#     t2ws.tags.append('scalingbbH')
#     t2ws.extra_options.append('--PO BRs_kappa_dependent=True')
#     t2ws.tags.append('couplingdependentBRs')
#     t2ws.run()

#____________________________________________________________________
# Scans for Giovanni
@flag_as_option
def t2ws_Yukawa_G0A(args):
    # Basically nominal, but force asimov and hzz
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    t2ws = base_t2ws(args)
    t2ws.tags.append('G0A')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_G0A_unreweighted(args):
    # This is actually not interesting
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    t2ws = base_t2ws(args, apply_reweighting=False)
    t2ws.tags.append('G0A')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_G0B(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    t2ws = base_t2ws_fitOnlyTotalXS(args)
    t2ws.tags.append('G0B')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_G0A_reweighted_scaledByMuTotalXS(args):
    # Basically nominal, but force asimov and hzz
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    t2ws = base_t2ws(args)
    t2ws.tags.append('G0A')
    t2ws.tags.append('scaledByMuTotalXS')
    t2ws.extra_options.append('--PO scale_with_mu_totalXS=True')
    t2ws.run()


@flag_as_option
def t2ws_Yukawa_G0B_reweighted(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    t2ws = base_t2ws_fitOnlyTotalXS(args, apply_reweighting=True)
    t2ws.tags.append('G0B')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_G1A(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    t2ws = base_t2ws(args, apply_reweighting=False)
    t2ws.tags.append('G1A')
    t2ws.extra_options.append('--PO BRs_kappa_dependent=True')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_G1A_unreweighted(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    t2ws = base_t2ws(args)
    t2ws.tags.append('G1A')
    t2ws.extra_options.append('--PO BRs_kappa_dependent=True')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_G1B_reweighted(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    t2ws = base_t2ws_fitOnlyTotalXS(args, apply_reweighting=True)
    t2ws.tags.append('G1B')
    t2ws.extra_options.append('--PO BRs_kappa_dependent=True')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_G1B(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    t2ws = base_t2ws_fitOnlyTotalXS(args)
    t2ws.tags.append('G1B')
    t2ws.extra_options.append('--PO BRs_kappa_dependent=True')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_G2A(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    t2ws = base_t2ws(args)
    t2ws.tags.append('G2A')
    t2ws.extra_options.append('--PO freely_floating_BRs=True')
    t2ws.run()


@flag_as_option
def t2ws_Yukawa_G1A_reweighted_noTheoryUnc(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    t2ws = base_t2ws(args, apply_reweighting=True, apply_theory_uncertainties=False)
    t2ws.tags.append('G1A')
    t2ws.tags.append('noTheoryUnc')
    t2ws.extra_options.append('--PO BRs_kappa_dependent=True')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_G1B_reweighted_noTheoryUnc(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    t2ws = base_t2ws_fitOnlyTotalXS(args, apply_reweighting=True, apply_theory_uncertainties=False)
    t2ws.tags.append('G1B')
    t2ws.tags.append('noTheoryUnc')
    t2ws.extra_options.append('--PO BRs_kappa_dependent=True')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_G1A_reweighted_noTheoryUnc_scaledByMuTotalXS(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    t2ws = base_t2ws(args, apply_reweighting=True, apply_theory_uncertainties=False)
    t2ws.tags.append('G1A')
    t2ws.tags.append('noTheoryUnc')
    t2ws.tags.append('scaledByMuTotalXS')
    t2ws.extra_options.append('--PO BRs_kappa_dependent=True')
    t2ws.extra_options.append('--PO scale_with_mu_totalXS=True')
    t2ws.run()

# 
@flag_as_option
def t2ws_Yukawa_noTheoryUnc(args):
    t2ws = base_t2ws(args, apply_theory_uncertainties=False)
    t2ws.tags.append('noTheoryUnc')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_uncorrelatedTheoryUnc(args):
    t2ws = base_t2ws(args, apply_theory_uncertainties=False)
    add_theory_uncertainties(t2ws, uncorrelated=True)
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_lumiScale(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('lumiScale')
    t2ws.extra_options.append('--PO lumiScale=True')
    t2ws.run()

@flag_as_option
def t2ws_Yukawa_profiledTotalXS(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('profiledTotalXS')
    t2ws.extra_options.append('--PO ProfileTotalXS=True')
    t2ws.run()


# @flag_as_option
# def t2ws_Yukawa_fitOnlyNormalization(args):
#     t2ws = base_t2ws(args)
#     t2ws.tags.append('fitOnlyNormalization')
#     add_theory_uncertainties(t2ws)
#     t2ws.extra_options.append('--PO FitOnlyNormalization=True')
#     t2ws.run()

@flag_as_option
def t2ws_Yukawa_BRcouplingDependency(args):
    t2ws = base_t2ws(args)
    t2ws.tags.append('BRcouplingDependency')
    t2ws.extra_options.append('--PO FitBR=True')
    do_BR_uncertainties = False
    # if do_BR_uncertainties:
    #     t2ws.extra_options.append('--PO DoBRUncertainties=True')
    #     t2ws.tags.append('withBRUnc')
    t2ws.run()

# @flag_as_option
# def t2ws_Yukawa_fitRatioOfBRs(args):
#     t2ws = base_t2ws(args)
#     t2ws.tags.append('fitRatioOfBRs')
#     t2ws.extra_options.append('--PO FitRatioOfBRs=True')
#     t2ws.run()



