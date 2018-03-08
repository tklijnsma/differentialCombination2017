from OptionHandler import flag_as_option, flag_as_parser_options

import os, logging, copy
from os.path import *

import sys
sys.path.append('src')
import CombineToolWrapper

import LatestPaths
import differentials
import differentialutils

#____________________________________________________________________
datestr = differentials.core.datestr()

def basic_config(args, hurry=False):
    config = CombineToolWrapper.CombineConfig(args)
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
    config.datacard = nominal_datacard(args)
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_yukawa_uncorrelatedTheoryUnc(args):
    config = basic_config(args)
    config.datacard = 'out/workspaces_Feb23/combination_Yukawa_reweighted_uncorrelatedTheoryUnc_uncorrelatedTheoryUnc.root'
    config.subDirectory += '_uncorrelatedTheoryUnc'
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_yukawa_fitOnlyNormalization(args):
    config = basic_config(args)
    config.datacard = LatestPaths.ws_combined_Yukawa_profiledTotalXS_fitOnlyNormalization
    config.subDirectory += '_fitOnlyNormalization'
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_yukawa_oneKappa(args):
    for kappa in ['kappac', 'kappab']:
        config = basic_config(args)
        config.datacard = nominal_datacard(args)

        config.nPoints       = 72
        config.nPointsPerJob = 3
        config.queue         = 'short.q'
        config.POIs          = [ kappa ]

        otherKappa = { 'kappab' : 'kappac', 'kappac' : 'kappab' }[kappa]
        config.floatNuisances.append(otherKappa)
        config.subDirectory += '_oneKappa_' + kappa

        if args.lumiScale:
            config.datacard = 'out/workspaces_Feb12/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_lumiScale.root'
            config.subDirectory += '_lumiStudy'
            config.freezeNuisances.append('lumiScale')
            config.hardPhysicsModelParameters.append( 'lumiScale=8.356546' )
        scan_directly(config)

@flag_as_option
def scan_yukawa_lumiScale(args):
    assert_asimov(args)
    config = basic_config(args)
    # config.datacard = LatestPaths.ws_combined_Yukawa_lumiScalable
    # config.datacard = 'out/workspaces_Feb12/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_lumiScale.root'
    config.datacard = 'out/workspaces_Feb23/combination_Yukawa_reweighted_lumiScale.root'
    config.freezeNuisances.append('lumiScale')
    config.hardPhysicsModelParameters.append( 'lumiScale=8.356546' )
    config.subDirectory += '_lumiStudy'
    # config.hardPhysicsModelParameters.append( 'lumiScale=83.56546' )
    # config.subDirectory += '_lumiStudy3000fb'
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_yukawa_BRdependent(args):
    assert_asimov(args)
    config = basic_config(args)
    # config.datacard = LatestPaths.ws_combined_Yukawa_couplingDependentBR
    config.datacard = 'out/workspaces_Feb23/combination_Yukawa_reweighted_BRcouplingDependency.root'
    config.subDirectory += '_couplingDependentBR'
   
    config.PhysicsModelParameters.append('kappa_V=0.999')

    # config.freezeNuisances.append('kappa_V')

    config.floatNuisances.append('kappa_V')
    config.PhysicsModelParameterRanges.append('kappa_V=-1000.0:1.0')
    config.subDirectory += '_floatKappaV'

    differentialutils.run_postfit_fastscan_scan(config)


@flag_as_option
def scan_yukawa_BRdependent_and_profiledTotalXS(args):
    assert_asimov(args)
    config = basic_config(args)
    config.datacard = LatestPaths.ws_combined_Yukawa_couplingDependentBR_profiledTotalXS
    config.subDirectory += '_couplingDependentBR_profiledTotalXS'
    config.fix_parameter_at_value('kappa_V', 0.999)
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_yukawa_profiledTotalXS(args):
    assert_asimov(args)
    config = basic_config(args)
    # config.datacard = LatestPaths.ws_combined_Yukawa_profiledTotalXS
    config.datacard = 'out/workspaces_Feb23/combination_Yukawa_reweighted_profiledTotalXS.root'
    config.subDirectory += '_profiledTotalXS'
    config.nPoints = 60*60
    config.nPointsPerJob = 15
    config.deltaNLLCutOff = 120.
    differentialutils.run_postfit_fastscan_scan(config)

