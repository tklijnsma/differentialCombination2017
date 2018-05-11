from OptionHandler import flag_as_option, flag_as_parser_options

import os, logging, copy
from os.path import *

# import sys
# sys.path.append('src')
# import CombineToolWrapper
import differentials.combine.combine as combine

import LatestPaths
import differentials
import differentialutils

#____________________________________________________________________
datestr = differentials.core.datestr()

def basic_config(args, hurry=False):
    config = combine.CombineConfig(args)
    config.onBatch       = True
    config.queue         = 'all.q'
    if hurry:
        logging.warning( 'Running with quick settings' )
        config.nPointsPerJob = 5
        config.queue         = 'short.q'

    if args.hzz:
        config.nPointsPerJob = 320
        config.queue         = 'short.q'

    if args.asimov:
        config.asimov = True
    else:
        config.asimov = False

    config.decay_channel = differentialutils.get_decay_channel_tag(args)

    config.nPoints = 70*70
    kappab_ranges = [ -15., 15. ]
    kappac_ranges = [ -35., 35. ]
    if config.args.lumiScale:
        config.nPoints = 40*40
        kappab_ranges = [ -7., 7. ]
        kappac_ranges = [ -17., 17. ]

    config.POIs = [ 'kappab', 'kappac' ]
    config.PhysicsModelParameters = [ 'kappab=1.0', 'kappac=1.0' ]
    config.PhysicsModelParameterRanges = [
        'kappab={0},{1}'.format( kappab_ranges[0], kappab_ranges[1] ),
        'kappac={0},{1}'.format( kappac_ranges[0], kappac_ranges[1] )
        ]
    config.subDirectory = 'out/Scan_Yukawa_{0}_{1}'.format(datestr, config.decay_channel)
    config.deltaNLLCutOff = 70.
    return config

def set_combination_and_asimov(args):
    args = differentialutils.set_one_decay_channel(args, 'combination')
    args.asimov = True
    return args

@flag_as_option
def all_scans_Yukawa(args_original):
    args = copy.deepcopy(args_original)

    # for asimov in [ False, True ]:
    #     args.asimov = asimov
    #     for decay_channel in ['hgg', 'hzz', 'combination']:
    #         set_decay_channel(args, decay_channel)
    #         # t2ws_Yukawa_nominal(args)
    #         differentialutils.try_call_function_with_args(scan_yukawa, args)

    differentialutils.set_one_decay_channel(args, 'combination')
    args.asimov = True

    fns = [
        scan_yukawa_uncorrelatedTheoryUnc,
        # scan_yukawa_fitOnlyNormalization,
        # scan_yukawa_oneKappa,
        scan_yukawa_lumiScale,
        scan_yukawa_BRdependent,
        # scan_yukawa_BRdependent_and_profiledTotalXS,
        scan_yukawa_profiledTotalXS,
        ]
    for fn in fns:
        differentialutils.try_call_function_with_args(fn, args)



@flag_as_option
def scan_yukawa(args):
    config = basic_config(args)
    datacard_dict = LatestPaths.ws.yukawa.nominal
    config.datacard = datacard_dict[differentialutils.get_decay_channel_tag(args)]
    # if args.combination:
    #     config.datacard = 'out/workspaces_Mar09/combination_Yukawa_reweighted_nominal.root'
    # elif args.hgg:
    #     config.datacard = 'out/workspaces_Mar09/hgg_Yukawa_reweighted_nominal.root'
    differentialutils.run_postfit_fastscan_scan(config)


yukawa_G = differentials.core.AttrDict()
yukawa_G.G0A = 'out/workspaces_May11/combination_Yukawa_reweighted_G0A.root'
yukawa_G.G0B = 'out/workspaces_May11/combination_Yukawa_G0B.root'
yukawa_G.G1A = 'out/workspaces_May11/combination_Yukawa_reweighted_G1A.root'
yukawa_G.G1B = 'out/workspaces_May11/combination_Yukawa_G1B.root'
yukawa_G.G2A = 'out/workspaces_May11/combination_Yukawa_reweighted_G2A.root'

@flag_as_option
def scan_yukawa_G0B(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    config = basic_config(args)
    config.queue         = 'short.q'
    config.nPointsPerJob = 350
    config.tags.append('G0B')
    config.datacard = yukawa_G.G0B
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_yukawa_G1B(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    config = basic_config(args)
    config.queue         = 'short.q'
    config.nPointsPerJob = 350
    config.set_parameter_range('kappab', -7., 7.)
    config.set_parameter_range('kappac', -14., 14.)
    config.tags.append('G1B')
    config.datacard = yukawa_G.G1B
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_yukawa_G1BKV(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    config = basic_config(args)
    config.queue         = 'short.q'
    config.nPointsPerJob = 350
    config.set_parameter_range('kappab', -7., 7.)
    config.set_parameter_range('kappac', -14., 14.)
    config.set_parameter_range('kappa_V', -10000., 1.)
    config.hardPhysicsModelParameters.append('kappa_V=0.999')
    config.tags.append('G1BKV')
    config.datacard = yukawa_G.G1B
    differentialutils.run_postfit_fastscan_scan(config)


@flag_as_option
def scan_yukawa_G0A(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    config = basic_config(args)
    config.tags.append('G0A')
    config.datacard = yukawa_G.G0A
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_yukawa_G1A(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    config = basic_config(args)
    config.set_parameter_range('kappab', -7., 7.)
    config.set_parameter_range('kappac', -14., 14.)
    config.tags.append('G1A')
    config.datacard = yukawa_G.G1A
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_yukawa_G2A(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    config = basic_config(args)
    config.tags.append('G2A')
    config.datacard = yukawa_G.G2A
    differentialutils.run_postfit_fastscan_scan(config)



@flag_as_option
def scan_yukawa_withBRuncertainties(args):
    args = differentialutils.set_one_decay_channel(args, 'hzz', asimov=True)
    config = basic_config(args)
    config.tags.append('withBRuncertainties')
    config.datacard = 'out/workspaces_May08/hzz_Yukawa_reweighted_withBRuncertainties.root'
    differentialutils.run_postfit_fastscan_scan(config)


@flag_as_option
def scan_yukawa_unreweighted(args):
    config = basic_config(args)
    datacard_dict = LatestPaths.ws.yukawa.unreweighted
    config.datacard = datacard_dict[differentialutils.get_decay_channel_tag(args)]
    # if args.combination:
    #     config.datacard = 'out/workspaces_Mar09/combination_Yukawa_nominal.root'
    # elif args.hgg:
    #     config.datacard = 'out/workspaces_Mar09/hgg_Yukawa_nominal.root'
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_yukawa_oneKappa(args):
    args = differentialutils.set_one_decay_channel(args, 'combination')
    for kappa in [
            # 'kappac',
            'kappab'
            ]:
        config = basic_config(args)
        config.datacard = LatestPaths.ws.yukawa.nominal[config.decay_channel]
        config.nPoints       = 72
        config.nPointsPerJob = 3
        config.queue         = 'short.q'
        config.POIs          = [ kappa ]

        otherKappa = { 'kappab' : 'kappac', 'kappac' : 'kappab' }[kappa]
        config.floatNuisances.append(otherKappa)
        config.tags.append('oneKappa_' + kappa)

        for i, range_str in enumerate(config.PhysicsModelParameterRanges):
            if range_str.startswith('kappab'):
                config.PhysicsModelParameterRanges[i] = 'kappab={0},{1}'.format(-7.0, 9.0)
        if kappa == 'kappab':
            config.nPoints = 80
            config.nPointsPerJob = 2

        if args.lumiScale:
            config.datacard = 'out/workspaces_Feb12/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_lumiScale.root'
            config.freezeNuisances.append('lumiScale')
            config.tags.append('lumi300fb')
            config.hardPhysicsModelParameters.append( 'lumiScale=8.356546' )

        differentialutils.scan_directly(config)

        # DONT DO THIS:
        # differentialutils.run_postfit_fastscan_scan(config)
        # This relies on the fastscan-first-code, which is NOT implemented in combine
        # (only works for 2D now)

@flag_as_option
def scan_yukawa_uncorrelatedTheoryUnc(real_args):
    args = set_combination_and_asimov(real_args)
    config = basic_config(args)
    config.datacard = LatestPaths.ws.yukawa.uncorrelatedTheoryUnc
    config.tags.append('uncorrelatedTheoryUnc')
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_yukawa_noTheoryUnc(real_args):
    args = set_combination_and_asimov(real_args)
    config = basic_config(args)
    config.datacard = LatestPaths.ws.yukawa.noTheoryUnc
    config.tags.append('noTheoryUnc')
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_yukawa_profiledTotalXS(real_args):
    args = set_combination_and_asimov(real_args)
    config = basic_config(args)
    # config.datacard = LatestPaths.ws_combined_Yukawa_profiledTotalXS
    config.datacard = LatestPaths.ws.yukawa.profiledTotalXS
    config.tags.append('profiledTotalXS')
    config.nPoints = 60*60
    config.nPointsPerJob = 15
    config.deltaNLLCutOff = 120.
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_yukawa_lumiScale(real_args):
    args = set_combination_and_asimov(real_args)
    config = basic_config(args)
    config.datacard = LatestPaths.ws.yukawa.lumiScale
    config.freezeNuisances.append('lumiScale')
    config.hardPhysicsModelParameters.append( 'lumiScale=8.356546' )
    config.tags.append('lumi300fb')
    # config.hardPhysicsModelParameters.append( 'lumiScale=83.56546' )
    # config.tags.append('lumi3000fb')
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_yukawa_BRdependent(real_args):
    args = set_combination_and_asimov(real_args)
    config = basic_config(args)
    config.datacard = LatestPaths.ws.yukawa.BRcouplingDependency

    config.nPoints = 60*60
    config.nPointsPerJob = 10

    config.PhysicsModelParameterRanges = [
        'kappab={0},{1}'.format(-25., 25),
        'kappac={0},{1}'.format(-45., 45)
        ]

    config.tags.append('couplingDependentBR')
    config.PhysicsModelParameters.append('kappa_V=0.999')

    # config.freezeNuisances.append('kappa_V')
    config.floatNuisances.append('kappa_V')
    config.PhysicsModelParameterRanges.append('kappa_V=-1000.0,1.0')
    config.tags.append('floatKappaV')

    differentialutils.run_postfit_scan(config)


#____________________________________________________________________
# Need reimplementation in the model

@flag_as_option
def scan_yukawa_fitOnlyNormalization(real_args):
    args = set_combination_and_asimov(real_args)
    config = basic_config(args)
    config.datacard = LatestPaths.ws_combined_Yukawa_profiledTotalXS_fitOnlyNormalization
    config.tags.append('fitOnlyNormalization')
    differentialutils.run_postfit_fastscan_scan(config)


@flag_as_option
def scan_yukawa_BRdependent_and_profiledTotalXS(real_args):
    args = set_combination_and_asimov(real_args)
    config = basic_config(args)
    config.datacard = LatestPaths.ws_combined_Yukawa_couplingDependentBR_profiledTotalXS
    config.subDirectory += '_couplingDependentBR_profiledTotalXS'
    config.fix_parameter_at_value('kappa_V', 0.999)
    differentialutils.run_postfit_fastscan_scan(config)

