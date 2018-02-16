from scans_yukawa import (
    run_postfit_fastscan_scan,
    scan_directly,
    assert_decay_channel,
    assert_asimov
    )

from OptionHandler import flag_as_option, flag_as_parser_options

import LatestPaths
import sys
sys.path.append('src')
import Commands
import CombineToolWrapper
import differentialTools

from time import strftime
datestr = strftime( '%b%d' )

import os
from os.path import *

#____________________________________________________________________
def basic_config(args, hurry=False):
    assert_highpt(args)
    config = CombineToolWrapper.CombineConfig(args)
    config.onBatch       = True
    config.queue         = 'all.q'
    if hurry:
        Commands.warning( 'Running with quick settings' )
        config.nPointsPerJob = 5
        config.queue         = 'short.q'

    if args.hzz:
        config.nPointsPerJob = 320
        config.queue         = 'short.q'

    if args.asimov:
        config.asimov = True
    else:
        config.asimov = False

    config.decay_channel = differentialTools.get_decay_channel_tag(args)

    if args.combWithHbb or args.hbb:
        raise NotImplementedError(
            'combWithHbb and hbb need workspaces with new binning first.'
            )
        config.default_minimizer_settings = False
        config.extraOptions.extend([
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            ])

    config.deltaNLLCutOff = 70.
    config.nPoints = 200**2

    ct_ranges = [ -8.5, 8.5 ]
    cg_ranges = [ -0.65, 0.65 ]
    config.POIs = [ 'ct', 'cg' ]
    config.PhysicsModelParameters = [ 'ct=1.0', 'cg=0.0' ]
    config.PhysicsModelParameterRanges = [
        'ct={0},{1}'.format( ct_ranges[0], ct_ranges[1] ),
        'cg={0},{1}'.format( cg_ranges[0], cg_ranges[1] )
        ]
    config.subDirectory = 'out/Scan_Top_{0}_{1}'.format(datestr, config.decay_channel)
    if config.args.highpt: config.subDirectory = 'out/Scan_TopHighPt_{0}_{1}'.format(datestr, config.decay_channel)

    return config

def nominal_datacard(args):
    if args.highpt:
        if args.hgg:
            dc = LatestPaths.ws_hgg_TopHighPt
        elif args.hzz:
            dc = LatestPaths.ws_hzz_TopHighPt
        elif args.combWithHbb:
            dc = LatestPaths.ws_combWithHbb_TopHighPt
        else:
            dc = LatestPaths.ws_combined_TopHighPt
    else:
        if args.hgg:
            dc = LatestPaths.ws_hgg_Top
        elif args.hzz:
            dc = LatestPaths.ws_hzz_Top
        elif args.combWithHbb:
            raise LookupError('No datacard defined yet')
        else:
            dc = LatestPaths.ws_combined_Top
    return dc

def assert_highpt(args):
    if not args.highpt:
        Commands.warning('Probably no reason anymore to not use --highpt; setting --highpt to True')
        args.highpt = True

#____________________________________________________________________
@flag_as_option
def scan_top(args):
    config = basic_config(args)
    config.datacard = nominal_datacard(args)
    run_postfit_fastscan_scan(config)

@flag_as_option
def scan_top_lumiStudy(args):
    assert_asimov(args)
    config.datacard = LatestPaths.ws_combined_Top_lumiScalable
    config.freezeNuisances.append('lumiScale')
    config.hardPhysicsModelParameters.append( 'lumiScale=8.356546' )
    config.subDirectory += '_lumiStudy'
    # config.hardPhysicsModelParameters.append( 'lumiScale=83.56546' )
    # config.subDirectory += '_lumiStudy3000fb'

@flag_as_option
def scan_top_profiledTotalXS(args):
    assert_asimov(args)
    config.datacard = LatestPaths.ws_combined_Top_profiledTotalXS
    config.tags.append('profiledTotalXS')
    run_postfit_fastscan_scan(config)

@flag_as_option
def scan_top_fitOnlyNormalization(args):
    assert_asimov(args)
    config.datacard =(
        LatestPaths.ws_combined_Top_profiledTotalXS_fitOnlyNormalization
        if not args.highpt else
        LatestPaths.ws_combined_TopHighPt_profiledTotalXS_fitOnlyNormalization
        )
    config.tags.append('fitOnlyNormalization')
    run_postfit_fastscan_scan(config)

@flag_as_option
def scan_top_BRdependent_and_profiledTotalXS(args):
    assert_asimov(args)
    config.datacard = LatestPaths.ws_combined_Top_couplingDependentBR_profiledTotalXS
    config.subDirectory += '_couplingDependentBR_profiledTotalXS'
    config.fix_parameter_at_value('kappa_V', 0.999)
    run_postfit_fastscan_scan(config)



#____________________________________________________________________
@flag_as_option
def couplingScan_TopCtCb(args):
    print 'Now in couplingScan_TopCtCb'

    raise NotImplementedError(
        'Re-implement this properly!'
        )

    if not args.highpt:
        Commands.warning('Probably no reason anymore to not use --highpt; setting --highpt to True')
        args.highpt = True

    TheoryCommands.set_plot_dir( 'plots_{0}_TopCtCb'.format(datestr) )

    scan = basic_scan_instance(args)
    scan.deltaNLLCutOff = 50.
    scan.nPoints = 120**2
    ct_ranges = [ -0.1, 2. ]
    cb_ranges = [ -10.0, 16.0 ]
    scan.POIs = [ 'ct', 'cb' ]
    scan.PhysicsModelParameterRanges = [
        'ct={0},{1}'.format( ct_ranges[0], ct_ranges[1] ),
        'cb={0},{1}'.format( cb_ranges[0], cb_ranges[1] )
        ]
    scan.subDirectory = 'out/Scan_TopCtCb_{0}'.format(datestr)
    if args.highpt: scan.subDirectory = 'out/Scan_TopCtCbHighPt_{0}'.format(datestr)


    # ======================================
    # Determine the physics

    def nominal_datacard(args):
        if args.highpt:
            if args.hgg:
                dc = LatestPaths.ws_hgg_TopCtCbHighPt
            elif args.hzz:
                dc = LatestPaths.ws_hzz_TopCtCbHighPt
            elif args.combWithHbb:
                raise LookupError('No datacard defined yet')
            else:
                dc = LatestPaths.ws_combined_TopCtCbHighPt
        else:
            if args.hgg:
                dc = LatestPaths.ws_hgg_TopCtCb
            elif args.hzz:
                dc = LatestPaths.ws_hzz_TopCtCb
            elif args.combWithHbb:
                raise LookupError('No datacard defined yet')
            else:
                dc = LatestPaths.ws_combined_TopCtCb
        return dc

    if args.nominal:
        datacard = nominal_datacard(args)

    else:
        print 'Pass physics option'
        return

    scan.run()