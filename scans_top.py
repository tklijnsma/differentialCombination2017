from OptionHandler import flag_as_option, flag_as_parser_options

import logging
import copy
from copy import deepcopy
import os
from os.path import *

# import sys
# sys.path.append('src')
# import CombineToolWrapper
import differentials.combine.combine as combine

import differentials
import differentialutils
import LatestPaths


#____________________________________________________________________
datestr = differentials.core.datestr()

def set_combWithHbb_and_asimov(args):
    args = differentialutils.set_one_decay_channel(args, 'combWithHbb')
    args.asimov = True    
    return args

def basic_config(args, hurry=False):
    # assert_highpt(args)
    config = combine.CombineConfig(args)
    config.onBatch       = True
    config.queue         = 'all.q'
    if hurry:
        logging.warning('Running with quick settings')
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

    if args.combWithHbb or args.hbb:
        # raise NotImplementedError(
        #     'combWithHbb and hbb need workspaces with new binning first.'
        #     )
        config.minimizer_settings.extend([
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            ])

    config.deltaNLLCutOff = 70.
    config.nPoints = 200**2

    # default ranges; should cover general usecases
    ct_ranges = [ -8.5, 8.5 ]
    cg_ranges = [ -0.65, 0.65 ]

    # Overwrite with some optimized ranges to reduce the number of jobs
    # Can only be done with a priori knowledge
    if args.asimov:
        if args.hgg:
            ct_ranges = [ -6.0, 6.0 ]
            cg_ranges = [ -0.45, 0.45 ]
            config.nPointsPerJob = 32
        elif args.combination:
            ct_ranges = [ -6.0, 6.0 ]
            cg_ranges = [ -0.45, 0.45 ]
            config.nPointsPerJob = 20
        elif args.combWithHbb:
            ct_ranges = [ -5.0, 5.0 ]
            cg_ranges = [ -0.35, 0.35 ]
            config.nPointsPerJob = 12
    else:
        if args.hgg:
            ct_ranges = [ -7.5, 7.5 ]
            cg_ranges = [ -0.50, 0.50 ]
            config.nPointsPerJob = 32
        elif args.combination:
            ct_ranges = [ -7.0, 7.0 ]
            cg_ranges = [ -0.50, 0.50 ]
            config.nPointsPerJob = 20
        elif args.combWithHbb:
            ct_ranges = [ -6.0, 6.0 ]
            cg_ranges = [ -0.40, 0.40 ]
            config.nPointsPerJob = 12

    config.POIs = [ 'ct', 'cg' ]
    config.PhysicsModelParameters = [ 'ct=1.0', 'cg=0.0' ]
    config.PhysicsModelParameterRanges = [
        'ct={0},{1}'.format( ct_ranges[0], ct_ranges[1] ),
        'cg={0},{1}'.format( cg_ranges[0], cg_ranges[1] )
        ]
    config.subDirectory = 'out/Scan_{0}_Top_{1}'.format(datestr, config.decay_channel)
    # if config.args.highpt: config.subDirectory = 'out/Scan_TopHighPt_{0}_{1}'.format(datestr, config.decay_channel)

    return config


def basic_config_ctcb(args):
    # assert_highpt(args)
    config = combine.CombineConfig(args)
    config.onBatch       = True
    config.queue         = 'all.q'

    if args.hzz:
        config.nPointsPerJob = 320
        config.queue         = 'short.q'
    if args.asimov:
        config.asimov = True
    else:
        config.asimov = False

    config.decay_channel = differentialutils.get_decay_channel_tag(args)

    if args.combWithHbb or args.hbb:
        config.minimizer_settings.extend([
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            ])

    config.deltaNLLCutOff = 70.
    config.nPoints = 85**2

    config.nPointsPerJob = 20
    if args.hgg:
        config.nPointsPerJob = 32
    elif args.combWithHbb:
        config.nPointsPerJob = 12

    # default ranges; should cover general usecases
    # ct_ranges = [ -0.5, 2.0 ]
    # cb_ranges = [ -20., 22. ]

    # Default: Hzz ranges (the widest)
    ct_ranges = [ -1.55, 1.55 ]
    cb_ranges = [ -20., 20. ]
    if args.hzz:
        ct_ranges = [ -1.7, 1.7 ]
        cb_ranges = [ -23., 23. ]

    config.POIs = [ 'ct', 'cb' ]
    config.PhysicsModelParameters = [ 'ct=1.0', 'cb=1.0' ]
    config.PhysicsModelParameterRanges = [
        'ct={0},{1}'.format( ct_ranges[0], ct_ranges[1] ),
        'cb={0},{1}'.format( cb_ranges[0], cb_ranges[1] )
        ]
    config.subDirectory = 'out/Scan_{0}_TopCtCb_{1}'.format(datestr, config.decay_channel)

    return config


@flag_as_option
def scan_top_scalingttH(args):
    args = differentialutils.set_one_decay_channel(args, 'combWithHbb', asimov=True)
    config = basic_config(args)
    # config.datacard = 'out/workspaces_May17/combWithHbb_Top_reweighted_scalingttH.root'
    config.datacard = 'out/workspaces_May18/combWithHbb_Top_reweighted_scalingttH.root'
    config.tags.append('scalingttH')
    config.set_parameter_range('ct', -0.4, 4.4)
    config.set_parameter_range('cg', -0.28, 0.1)
    config.nPoints = 60*60
    config.nPointsPerJob = 6
    differentialutils.run_postfit_fastscan_scan(config)


#____________________________________________________________________
@flag_as_option
def scan_top_nominal_all(real_args):
    args = deepcopy(real_args)
    for asimov in [
        True,
        False
        ]:
        args.asimov = asimov
        for dc in [
            # 'combWithHbb',
            'combination',
            # 'hgg'
            ]:
            differentials.core.set_one_decay_channel(args, dc)
            scan_top(args)

@flag_as_option
def scan_top_hzz_testrun(real_args):
    args = deepcopy(real_args)
    differentials.core.set_one_decay_channel(args, 'hzz')
    args.asimov = True
    config = basic_config(args)
    config.datacard = LatestPaths.ws.top.nominal[differentialutils.get_decay_channel_tag(args)]
    config.nPoints = 10*10
    config.nPointsPerJob = 50
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_top(args):
    config = basic_config(args)
    config.datacard = LatestPaths.ws.top.nominal[differentialutils.get_decay_channel_tag(args)]
    config.nPoints = 110*110
    config.nPointsPerJob = int(0.5*config.nPointsPerJob)
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_top_lumiscale(args):
    args = differentialutils.set_one_decay_channel(args, 'combWithHbb', asimov=True)
    config = basic_config(args)
    config.datacard = LatestPaths.ws.top.lumiScale
    config.freezeNuisances.append('lumiScale')
    config.hardPhysicsModelParameters.append( 'lumiScale=8.356546' )
    config.tags.append('lumiStudy')
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_topctcb(args):
    config = basic_config_ctcb(args)
    config.datacard = LatestPaths.ws.topctcb.nominal[differentialutils.get_decay_channel_tag(args)]
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_topctcb_lumiscale(args):
    args = differentialutils.set_one_decay_channel(args, 'combWithHbb', asimov=True)
    config = basic_config_ctcb(args)
    config.datacard = LatestPaths.ws.topctcb.lumiScale
    config.freezeNuisances.append('lumiScale')
    config.hardPhysicsModelParameters.append( 'lumiScale=8.356546' )
    config.tags.append('lumiStudy')
    differentialutils.run_postfit_fastscan_scan(config)

#____________________________________________________________________

@flag_as_option
def scan_top_noBinsDropped(real_args):
    args = deepcopy(real_args)
    args.asimov = True
    config = basic_config(args)
    if args.combWithHbb:
        config.datacard = LatestPaths.ws.top.nominal.combWithHbb
    elif args.combination:
        config.datacard = LatestPaths.ws.top.nominal.combination
    else:
        raise NotImplementedError('Run with --combination or --combWithHbb')
    config.tags.append('noBinsDropped')
    config.nPoints = 120*120
    config.nPointsPerJob = int(0.5*config.nPointsPerJob)
    differentialutils.run_postfit_fastscan_scan(config)


@flag_as_option
def scan_top_last2BinsDropped(real_args):
    args = deepcopy(real_args)
    args.asimov = True
    config = basic_config(args)
    if args.combWithHbb:
        config.datacard = LatestPaths.ws.top.last2BinsDropped.combWithHbb
    elif args.combination:
        config.datacard = LatestPaths.ws.top.last2BinsDropped.combination
    else:
        raise NotImplementedError('Run with --combination or --combWithHbb')
    config.tags.append('last2BinsDropped')
    config.nPoints = 150*150
    differentialutils.run_postfit_fastscan_scan(config)


@flag_as_option
def scan_top_lastBinDropped(real_args):
    args = deepcopy(real_args)
    args.asimov = True
    config = basic_config(args)
    if args.combination:
        config.datacard = 'out/workspaces_Mar07/combination_Top_reweighted_lastBinDroppedHgg.root'
    elif args.combWithHbb:
        config.datacard = 'out/workspaces_Mar07/combWithHbb_Top_reweighted_lastBinDroppedHgg.root'
    else:
        raise NotImplementedError('Run with --combination or --combWithHbb')
    differentialutils.run_postfit_fastscan_scan(config)


@flag_as_option
def scan_top_lumiStudy(args):
    args = set_combWithHbb_and_asimov(args)
    config = basic_config(args)
    config.datacard = 'out/workspaces_Mar13/combWithHbb_Top_reweighted_lumiScale.root'
    config.freezeNuisances.append('lumiScale')
    config.hardPhysicsModelParameters.append( 'lumiScale=8.356546' )
    config.subDirectory += '_lumiStudy'
    # config.hardPhysicsModelParameters.append( 'lumiScale=83.56546' )
    # config.subDirectory += '_lumiStudy3000fb'
    differentialutils.run_postfit_fastscan_scan(config)


@flag_as_option
def scan_top_profiledTotalXS(args):
    args = set_combWithHbb_and_asimov(args)
    config = basic_config(args)
    config.datacard = 'out/workspaces_Mar13/combWithHbb_Top_reweighted_profiledTotalXS.root'
    config.tags.append('profiledTotalXS')
    config.nPoints = 150*150
    differentialutils.run_postfit_fastscan_scan(config)


#____________________________________________________________________

@flag_as_option
def scan_top_BRdependent(args):
    args = set_combWithHbb_and_asimov(args)
    config = basic_config(args)
    config.datacard = LatestPaths.ws.top.BRcouplingDependency

    config.nPoints = 80*80
    config.nPointsPerJob = 10

    config.PhysicsModelParameterRanges = [
        'ct={0},{1}'.format(-8.5, 8.5 ),
        'cg={0},{1}'.format(-0.65, 0.65 )
        ]

    config.tags.append('couplingDependentBR')
    config.PhysicsModelParameters.append('kappa_V=0.999')

    # config.freezeNuisances.append('kappa_V')
    config.floatNuisances.append('kappa_V')
    config.PhysicsModelParameterRanges.append('kappa_V=-1000.0,1.0')
    config.tags.append('floatKappaV')

    differentialutils.run_postfit_scan(config)


#____________________________________________________________________
@flag_as_option
def scan_top_fitOnlyNormalization(args):
    differentialutils.assert_asimov(args)
    config.datacard =(
        LatestPaths.ws_combined_Top_profiledTotalXS_fitOnlyNormalization
        if not args.highpt else
        LatestPaths.ws_combined_TopHighPt_profiledTotalXS_fitOnlyNormalization
        )
    config.tags.append('fitOnlyNormalization')
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_top_BRdependent_and_profiledTotalXS(args):
    differentialutils.assert_asimov(args)
    config.datacard = LatestPaths.ws_combined_Top_couplingDependentBR_profiledTotalXS
    config.subDirectory += '_couplingDependentBR_profiledTotalXS'
    config.fix_parameter_at_value('kappa_V', 0.999)
    differentialutils.run_postfit_fastscan_scan(config)

