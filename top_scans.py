from OptionHandler import flag_as_option, flag_as_parser_options

import logging
import copy
from copy import deepcopy
import os
from os.path import *

import differentials.combine.combine as combine

import differentials
import differentialutils
import LatestPaths

from math import pi

#____________________________________________________________________
datestr = differentials.core.datestr()

def basic_config(args, hurry=False):
    # assert_highpt(args)
    config = combine.CombineConfig(args)
    config.onBatch       = True
    config.queue         = 'all.q'
    config.asimov        = True if args.asimov else False
    config.decay_channel = differentialutils.get_decay_channel_tag(args)

    if args.combWithHbb or args.hbb:
        config.minimizer_settings.extend([
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            ])
    config.deltaNLLCutOff = 70.
    config.POIs = [ 'ct', 'cg' ]
    config.PhysicsModelParameters = [ 'ct=1.0', 'cg=0.0' ]
    config.subDirectory = 'out/Scan_{0}_Top_{1}'.format(datestr, config.decay_channel)
    return config

def basic_config_ctcb(args):
    # assert_highpt(args)
    config = combine.CombineConfig(args)
    config.onBatch       = True
    config.queue         = 'all.q'
    config.asimov        = True if args.asimov else False
    config.decay_channel = differentialutils.get_decay_channel_tag(args)

    if args.combWithHbb or args.hbb:
        config.minimizer_settings.extend([
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            ])

    config.deltaNLLCutOff = 70.
    config.POIs = [ 'ct', 'cb' ]
    config.PhysicsModelParameters = [ 'ct=1.0', 'cb=1.0' ]
    config.subDirectory = 'out/Scan_{0}_TopCtCb_{1}'.format(datestr, config.decay_channel)
    return config

#____________________________________________________________________
approval = differentials.core.AttrDict.create_tree(['ktcg', 'ktkb'], ['couplingdependentBRs', 'floatingBRs', 'fixedBRs'])

approval.ktcg.fixedBRs.combWithHbb = 'out/workspaces_May18/combWithHbb_Top_reweighted_scalingttH.root'
approval.ktcg.floatingBRs.combWithHbb = 'out/workspaces_May29/combWithHbb_Top_reweighted_scalingttH_floatingBRs_constrainedbbZZ.root'
approval.ktcg.floatingBRs.hgg = 'out/workspaces_Jun09/hgg_Top_reweighted_scalingttH_floatingBRs_constrainedbbZZ.root'
approval.ktcg.floatingBRs.hzz = 'out/workspaces_Jun09/hzz_Top_reweighted_scalingttH_floatingBRs_constrainedbbZZ.root'
approval.ktcg.couplingdependentBRs.combWithHbb = 'out/workspaces_May31/combWithHbb_Top_reweighted_scalingttH_couplingdependentBRs.root'
approval.ktcg.couplingdependentBRs.hgg = 'out/workspaces_Jun09/hgg_Top_reweighted_scalingttH_couplingdependentBRs.root'
approval.ktcg.couplingdependentBRs.hzz = 'out/workspaces_Jun09/hzz_Top_reweighted_scalingttH_couplingdependentBRs.root'

approval.ktkb.fixedBRs.combWithHbb = 'out/workspaces_May29/combWithHbb_TopCtCb_reweighted_scalingbbHttH.root'
# approval.ktkb.floatingBRs.combWithHbb = 'out/workspaces_May29/combWithHbb_TopCtCb_reweighted_scalingbbHttH_floatingBRs_constrainedbbZZ.root'
approval.ktkb.floatingBRs.combWithHbb = 'out/workspaces_Jun10/combWithHbb_TopCtCb_reweighted_scalingbbHttH_floatingBRs_constrainedbbZZ.root'
approval.ktkb.floatingBRs.hgg = 'out/workspaces_Jun09/hgg_TopCtCb_reweighted_scalingbbHttH_floatingBRs_constrainedbbZZ.root'
approval.ktkb.floatingBRs.hzz = 'out/workspaces_Jun09/hzz_TopCtCb_reweighted_scalingbbHttH_floatingBRs_constrainedbbZZ.root'
approval.ktkb.couplingdependentBRs.combWithHbb = 'out/workspaces_May29/combWithHbb_TopCtCb_reweighted_scalingbbHttH_couplingdependentBRs.root'
approval.ktkb.couplingdependentBRs.hgg = 'out/workspaces_Jun09/hgg_TopCtCb_reweighted_scalingbbHttH_couplingdependentBRs.root'
approval.ktkb.couplingdependentBRs.hzz = 'out/workspaces_Jun09/hzz_TopCtCb_reweighted_scalingbbHttH_couplingdependentBRs.root'


#____________________________________________________________________
@flag_as_option
def scan_top_scalingttH(args):
    config = basic_config(args)
    config.datacard = approval.ktcg.fixedBRs[config.decay_channel]
    config.tags.append('scalingttH')
    config.nPoints = 30*30
    config.set_parameter_range('ct', -0.4, 4.4)
    config.set_parameter_range('cg', -0.28, 0.1)
    config.nPointsPerJob = 6
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_top_scalingttH_couplingdependentBRs(args):
    config = basic_config(args)
    config.datacard = approval.ktcg.couplingdependentBRs[config.decay_channel]
    config.tags.append('scalingttH')
    config.tags.append('couplingdependentBRs')

    config.nPoints = 30*30
    # On combWithHbb, observed: around 4 hours for 10 points, but lots of variance
    # Good maximum per job is 10 (tiny bit of loss)
    config.nPointsPerJob = 10
    config.set_parameter_range('ct', -0.4, 4.4)
    config.set_parameter_range('cg', -0.28, 0.1)

    if args.hgg:
        config.set_parameter_range('ct', -0.3, 2.6)
        config.set_parameter_range('cg', -0.08, 0.12)
    elif args.hzz:
        config.set_parameter_range('ct', -0.8, 3.2)
        config.set_parameter_range('cg', -0.15, 0.16)
        config.nPointsPerJob = 300
        config.queue = 'short.q'

    if args.asimov and args.hgg:
        # Special scan for Vittorio
        config.nPoints = 20*20
        config.nPointsPerJob = 5
        config.set_parameter_range('ct', -0.3, 5.3)
        config.set_parameter_range('cg', -0.12, 0.2)

    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_top_scalingttH_floatingBRs_constrainedbbZZ(args):
    config = basic_config(args)
    config.datacard = approval.ktcg.floatingBRs[config.decay_channel]
    config.tags.append('scalingttH')
    config.tags.append('floatingBRs')
    config.tags.append('constrainedbbZZ')

    # config.nPoints = 30*30
    # Around 4.5 hours (good max of 6 hours) for 8 points (combWithHbb, observed)
    # Reasonably constant, can go up to 12 points
    # Keep same ranges for hgg/hzz for now
    # config.nPointsPerJob = 12
    # config.set_parameter_range('ct', -4.2, 4.2) 
    # config.set_parameter_range('cg', -0.30, 0.30)

    config.nPoints = 100*100
    config.nPointsPerJob = 12
    config.set_parameter_range('ct', -6.0, 6.0) 
    config.set_parameter_range('cg', -0.45, 0.45)

    if args.hzz:
        config.nPointsPerJob = 300
        config.queue = 'short.q'

    if args.asimov:
        config.nPointsPerJob = 10
        config.set_parameter_range('ct', -8.2, 8.2) 
        config.set_parameter_range('cg', -0.50, 0.50)

    differentialutils.run_postfit_scan(config)


#____________________________________________________________________
@flag_as_option
def scan_topctcb_scalingbbHttH(args):
    config = basic_config_ctcb(args)
    config.datacard = approval.ktkb.fixedBRs[config.decay_channel]
    config.tags.append('scalingbbHttH')
    config.set_parameter_range('ct', -1.55, 1.55)
    config.set_parameter_range('cb', -20., 20.)
    config.nPoints = 50*50
    config.nPointsPerJob = 12
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_topctcb_scalingbbHttH_couplingdependentBRs(args):
    config = basic_config_ctcb(args)
    config.datacard = approval.ktkb.couplingdependentBRs[config.decay_channel]
    config.tags.append('scalingbbHttH')
    config.tags.append('couplingdependentBRs')
    config.nPoints = 35*35
    # Up to max ~7 hours for 10 points (combWithHbb, observed), can go up to 12 points
    config.set_parameter_range('ct', -0.3, 2.5)
    config.set_parameter_range('cb', -2.8, 2.8)
    config.nPointsPerJob = 12
    if args.hgg:
        # kt-plane asym means hgg goes up to reasonably high kt
        config.nPoints = 25*25
        config.set_parameter_range('ct', 0.2, 3.3)
        config.set_parameter_range('cb', -1.8, 1.8)
        config.nPointsPerJob = 10
    if args.hzz:
        config.set_parameter_range('ct', -3.3, 3.3)
        config.set_parameter_range('cb', -3.3, 3.3)
        config.nPointsPerJob = 300
        config.queue = 'short.q'
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def scan_topctcb_scalingbbHttH_floatingBRs_constrainedbbZZ(args):
    config = basic_config_ctcb(args)
    config.datacard = approval.ktkb.floatingBRs[config.decay_channel]
    config.tags.append('scalingbbHttH')
    config.tags.append('floatingBRs')
    config.tags.append('constrainedbbZZ')

    # # around max 5, occasionally 5.5 hours for 7 points; can go up to 12
    # config.set_parameter_range('ct', -7.0, 7.0)
    # config.set_parameter_range('cb', -17., 17.)
    # config.nPoints = 25*25
    # config.nPointsPerJob = 10

    config.set_parameter_range('ct', -3.5, 3.5)
    config.set_parameter_range('cb', -15.5, 15.5)
    config.nPoints = 60*60
    config.nPointsPerJob = 10

    if args.hzz:
        config.nPointsPerJob = 300
        config.queue = 'short.q'
    differentialutils.run_postfit_scan(config)


#____________________________________________________________________
@flag_as_option
def scan_top_scalingttH_floatingBRs(args):
    raise RuntimeError('Probably not what I want!')
    args = differentialutils.set_one_decay_channel(args, 'combWithHbb')
    config = basic_config(args)
    # config.datacard = 'out/workspaces_May17/combWithHbb_Top_reweighted_scalingttH.root'
    config.datacard = 'out/workspaces_May22/combWithHbb_Top_reweighted_scalingttH_floatingBRs.root'
    config.tags.append('scalingttH')
    config.tags.append('floatingBRs')

    # Runs out of CPU
    # config.set_parameter_range('ct', -2.5, 2.5)
    # config.set_parameter_range('cg', -0.15, 0.15)
    # config.nPoints = 40*40
    # config.nPointsPerJob = 12

    # Even faster
    config.set_parameter_range('ct', -1.7, 1.7)
    config.set_parameter_range('cg', -0.10, 0.10)
    config.nPoints = 10*10
    config.nPointsPerJob = 4

    # differentialutils.run_postfit_fastscan_scan(config)
    differentialutils.run_postfit_scan(config)

@flag_as_option
def scan_top_scalingttH_floatingBRs_NONconstrainedbbZZ(args):
    args = differentialutils.set_one_decay_channel(args, 'combWithHbb')
    config = basic_config(args)
    config.datacard = 'out/workspaces_May30/combWithHbb_Top_reweighted_scalingttH_floatingBRs.root'

    config.tags.append('scalingttH')
    config.tags.append('floatingBRs')
    config.tags.append('NONconstrainedbbZZ')

    config.set_parameter_range('ct', -1.7, 1.7)
    config.set_parameter_range('cg', -0.10, 0.10)
    # config.nPoints = 10*10
    # config.nPointsPerJob = 1
    # config.queue = 'short.q'

    config.nPoints = 30*30
    config.nPointsPerJob = 7
    config.queue = 'all.q'
    differentialutils.run_postfit_scan(config)

@flag_as_option
def scan_topctcb_scalingbbHttH_floatingBRs(args):
    args = differentialutils.set_one_decay_channel(args, 'combWithHbb')
    config = basic_config_ctcb(args)
    config.datacard = 'out/workspaces_May29/combWithHbb_TopCtCb_reweighted_scalingbbHttH_floatingBRs.root'
    config.tags.append('scalingbbHttH')
    config.tags.append('floatingBRs')
    config.set_parameter_range('ct', -10.0, 10.0)
    config.set_parameter_range('cb', -20., 20.)
    config.nPoints = 40*40
    config.nPointsPerJob = 9
    differentialutils.run_postfit_scan(config)

#____________________________________________________________________

def basic_config_points(args):
    # assert_highpt(args)
    config = combine.CombineConfig(args)
    config.onBatch       = True
    config.queue         = 'short.q'
    config.asimov        = True if args.asimov else False
    config.decay_channel = differentialutils.get_decay_channel_tag(args)

    if args.combWithHbb or args.hbb:
        config.minimizer_settings.extend([
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            ])

    config.PhysicsModelParameters = [ 'ct=1.0', 'cg=0.0' ]
    config.subDirectory = 'out/Scan_{0}_TopPoints_{1}'.format(datestr, config.decay_channel)
    return config

def basic_config_radial(args):
    # assert_highpt(args)
    config = combine.CombineConfig(args)
    config.onBatch       = True
    config.queue         = 'short.q'

    config.nPointsPerJob = 5
    if args.hzz:
        config.nPointsPerJob = 320

    config.decay_channel = differentialutils.get_decay_channel_tag(args)

    if args.combWithHbb or args.hbb:
        config.minimizer_settings.extend([
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            ])

    config.deltaNLLCutOff = 60.
    config.nPoints = 50

    config.POIs = [ 'theta' ]
    config.PhysicsModelParameters = [ 'theta=0.0', 'r=1.0' ]
    config.set_parameter_range('theta', -0.25*pi, 0.75*pi)
    config.set_parameter_range('r', 0.0, 5.0)

    config.subDirectory = 'out/Scan_{0}_TopRadial_{1}'.format(datestr, config.decay_channel)

    return config


@flag_as_option
def scan_top_radialpoints(args):
    args = differentialutils.set_one_decay_channel(args, 'combWithHbb', asimov=True)    
    config = basic_config_points(args)
    config.datacard = 'out/workspaces_May22/combWithHbb_Top_reweighted_scalingttH_floatingBRs.root'
    config.make_unique_directory()

    pointsscan = differentials.combine.combine.CombineScanSinglePoints(config)
    pointsscan.add_point(ct = 1.0, cg = 0.0)
    # pointsscan.add_point(ct = 4.0, cg = -0.27)
    # pointsscan.add_point(ct = 4.0, cg = -0.25)
    # pointsscan.add_point(ct = 4.0, cg = -0.23)
    # pointsscan.add_point(ct = 4.0, cg = -0.20)
    # pointsscan.add_point(ct = 3.0, cg = -0.05)
    # pointsscan.add_point(ct = 3.0, cg = 0.05)
    # pointsscan.add_point(ct = 1.0, cg = 0.1)
    # pointsscan.add_point(ct = 0.0, cg = 0.3)
    # pointsscan.add_point(ct = -0.25, cg = 0.3)
    # pointsscan.add_point(ct = -0.5, cg = 0.3)
    # pointsscan.add_point(ct = -0.75, cg = 0.3)
    # pointsscan.add_point(ct = -1.0, cg = 0.3)
    # pointsscan.add_point(ct = -1.25, cg = 0.3)
    # pointsscan.add_point(ct = -1.5, cg = 0.3)
    # pointsscan.add_point(ct = -2.0, cg = 0.3)
    pointsscan.run()

# @flag_as_option
# def scan_top_radialktkg(args):
#     args = differentialutils.set_one_decay_channel(args, 'combWithHbb', asimov=True)    
#     config = basic_config_radial(args)
#     config.datacard = 'out/workspaces_May22/combWithHbb_Top_reweighted_scalingttH_radialctcg.root'
#     differentialutils.run_postfit_scan(config)

# @flag_as_option
# def scan_top_radialktkg_hzz(args):
#     args = differentialutils.set_one_decay_channel(args, 'hzz', asimov=True)    
#     config = basic_config_radial(args)
#     config.nPoints = 30
#     config.saveFunctions.extend(['ct', 'cg', 'r_ggH_PTH_0_15'])
#     config.datacard = 'out/workspaces_May22/hzz_Top_reweighted_scalingttH_radialctcg.root'
#     # config.nPointsPerJob = config.nPoints
#     config.nPointsPerJob = 10
#     differentialutils.scan_directly(config)

@flag_as_option
def scan_top_scalingttH_cgonly(args):
    args = differentialutils.set_one_decay_channel(args, 'combWithHbb')
    config = basic_config(args)
    # config.datacard = 'out/workspaces_May17/combWithHbb_Top_reweighted_scalingttH.root'
    config.datacard = 'out/workspaces_May22/combWithHbb_Top_reweighted_scalingttH_floatingBRs.root'
    config.tags.append('scalingttH')
    config.tags.append('cgonly')

    config.suppress_output = True
    config.queue = 'short.q'

    config.FLOAT_OTHER_POIS = False
    config.del_parameter_range('ct')
    config.del_poi('ct')
    config.set_parameter('ct', 1.0, hard=True)
    config.freeze_parameter('ct')

    if args.asimov:
        config.set_parameter_range('cg', -1.6, 0.2)
    else:
        config.set_parameter_range('cg', -0.1, 0.4)

    config.nPoints = 25
    config.nPointsPerJob = 1

    # differentialutils.run_postfit_fastscan_scan(config)
    differentialutils.run_postfit_scan(
        config,
        postfit_file = 'out/Scan_May25_Top_combWithHbb_scalingttH_cgonly_asimov/postfit_and_fastscan/higgsCombine_POSTFIT_ASIMOV_combWithHbb_Top_reweighted_scalingttH_floatingBRs.MultiDimFit.mH125.root'
        )


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

