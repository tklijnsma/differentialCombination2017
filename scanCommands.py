#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

from OptionHandler import flag_as_option

import LatestPaths

import sys
sys.path.append('src')
import Commands
import PhysicsCommands
import OneOfCommands
import TheoryCommands
import CorrelationMatrices
import MergeHGGWDatacards
import TheoryFileInterface
import PlotCommands
from Container import Container
from Parametrization import Parametrization, WSParametrization
import CombineToolWrapper

from time import strftime
datestr = strftime( '%b%d' )

import ROOT
from TheoryCommands import c
from TheoryCommands import SaveC
from TheoryCommands import GetPlotBase
from TheoryCommands import SetCMargins

import os, itertools, operator, re, argparse, random
from math import isnan, isinf, sqrt
from os.path import *
from glob import glob
from copy import deepcopy
from array import array


########################################
# Main
########################################

def AppendParserOptions( parser ):

    # parser.add_argument( '--scanCommands', default=False )
    # class CustomAction(argparse.Action):
    #     def __init__(self, option_strings, dest, **kwargs):
    #         # Only 1 argument allowed (this is basically the "store_true" action)
    #         super(CustomAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
    #     def __call__(self, parser, namespace, values, option_string=None):
    #         setattr( namespace, 'scanCommands', True )
    #         setattr( namespace, self.dest, True )
    #         # super(FooAction, self).__call__( parser, namespace, values, option_string=None )
    #         # setattr( namespace, self.dest, values )

    # parser.add_argument( '--couplingScan',                             action=CustomAction )
    # parser.add_argument( '--fastscanplotmock',                         action=CustomAction )

    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument( '--yukawa',  action='store_true' )
    group.add_argument( '--top',     action='store_true' )
    group.add_argument( '--topctcb', action='store_true' )

    parser.add_argument( '--nominal',                                  action='store_true' )
    parser.add_argument( '--highpt',                                   action='store_true' )
    parser.add_argument( '--lumiStudy',                                action='store_true' )
    parser.add_argument( '--profiledTotalXS',                          action='store_true' )
    parser.add_argument( '--fitOnlyNormalization',                     action='store_true' )
    parser.add_argument( '--BRdependent',                              action='store_true' )
    parser.add_argument( '--oneKappa', type=str, default=None )


########################################
# Methods
########################################    

def basic_scan_instance(
        args,
        postfitWS=None,
        fastscanFile=None,
        hurry=False,
        ):

    scan = CombineToolWrapper.CombineScan()

    if not(postfitWS is None):
        Commands.Warning('Using existing postfitWS {0}'.format(postfitWS))
        scan.postfitWS = postfitWS
    if not(fastscanFile is None):
        Commands.Warning( 'Using existing fastscan {0}'.format(fastscanFile) )
        scan.fastscanRootFile = fastscanFile

    scan.onBatch       = True
    scan.nPoints       = 100*100
    scan.nPointsPerJob = 20
    scan.queue         = 'all.q'
    scan.APPLY_FASTSCAN_FILTER = True
    
    if hurry:
        Commands.Warning( 'Running with quick settings' )
        scan.nPointsPerJob = 5
        scan.queue         = 'short.q'

    if args.asimov:
        scan.asimov = True

    if args.hzz:
        scan.nPointsPerJob = 320
        scan.queue         = 'short.q'

    if args.hgg:
        scan.tags.append('hgg')
    elif args.hzz:
        scan.tags.append('hzz')
    elif args.combWithHbb:
        scan.tags.append('combWithHbb')
    else:
        scan.tags.append('combined')

    return scan


#____________________________________________________________________
@flag_as_option
def couplingScan_Yukawa(args):
    print 'Now in couplingScan_Yukawa'

    expBinBoundaries = [ 0., 15., 30., 45., 85., 125. ]
    print 'Hardcoded binBoundaries for Yukawa:'
    print expBinBoundaries
    print ''

    TheoryCommands.SetPlotDir( 'plots_{0}_Yukawa'.format(datestr) )

    scan = basic_scan_instance(args)
    kappab_ranges = [ -15., 15. ]
    kappac_ranges = [ -35., 35. ]
    scan.POIs = [ 'kappab', 'kappac' ]
    scan.PhysicsModelParameters = [ 'kappab=1.0', 'kappac=1.0' ]
    scan.PhysicsModelParameterRanges = [
        'kappab={0},{1}'.format( kappab_ranges[0], kappab_ranges[1] ),
        'kappac={0},{1}'.format( kappac_ranges[0], kappac_ranges[1] )
        ]
    scan.subDirectory = 'out/Scan_Yukawa_{0}'.format(datestr)
    scan.deltaNLLCutOff = 70.


    # ======================================
    # Determine the physics

    if args.hgg:
        nominal_datacard = LatestPaths.ws_hgg_Yukawa
    elif args.hzz:
        nominal_datacard = LatestPaths.ws_hzz_Yukawa
    else:
        # nominal_datacard = LatestPaths.ws_combined_Yukawa
        Commands.Warning('Picking the reweighted ws')
        nominal_datacard = LatestPaths.ws_combined_Yukawa_reweighted

    if args.nominal:
        # scan.datacard = nominal_datacard
        Commands.Warning('OVERWRITING WITH HARDCODED WORKSPACE')
        scan.datacard = 'out/workspaces_Feb02/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_reweighted.root'

    elif args.lumiStudy:
        scan.datacard = LatestPaths.ws_combined_Yukawa_lumiScalable
        scan.PhysicsModelParameters.append( 'lumiScale=8.356546' )
        scan.tags.append('lumiStudy')

    elif args.fitOnlyNormalization:
        scan.datacard = LatestPaths.ws_combined_Yukawa_profiledTotalXS_fitOnlyNormalization
        scan.tags.append('fitOnlyNormalization')

    elif args.oneKappa:
        scan.fromPostfit = False
        scan.APPLY_FASTSCAN_FILTER = False
        scan.datacard = nominal_datacard
        scan.nPoints       = 72
        scan.nPointsPerJob = 3
        scan.queue         = 'short.q'
        scan.POIs          = [ args.oneKappa ]

        otherKappa = { 'kappab' : 'kappac', 'kappac' : 'kappab' }[args.oneKappa]
        scan.floatNuisances.append(otherKappa)
        scan.tags.append('oneKappa_' + args.oneKappa)

    elif args.BRdependent:
        scan.datacard = LatestPaths.ws_combined_Yukawa_couplingDependentBR
        scan.tags.append('couplingDependentBR')

        # scan.PhysicsModelParameters.append('kappa_V=0.99')
        # scan.PhysicsModelParameterRanges.append('kappa_V=-100.0,1.0')
        
        scan.PhysicsModelParameters.append('kappa_V=0.999')
        scan.freezeNuisances.append('kappa_V')

    elif args.profiledTotalXS:
        scan.datacard = LatestPaths.ws_combined_Yukawa_profiledTotalXS
        scan.tags.append('profiledTotalXS')

    else:
        print 'Pass physics option'
        return

    scan.Run()


#____________________________________________________________________
@flag_as_option
def couplingScan_Top(args):
    print 'Now in couplingScan_Top'
    if not args.highpt:
        Commands.Warning('Probably no reason anymore to not use --highpt; setting --highpt to True')
        args.highpt = True
    
    TheoryCommands.SetPlotDir( 'plots_{0}_Top'.format(datestr) )

    scan = basic_scan_instance(args)
    scan.deltaNLLCutOff = 70.
    scan.nPoints = 200**2

    # ct_ranges = [ -1., 2. ]
    # cg_ranges = [ -0.1, 0.2 ]
    ct_ranges = [ -8.5, 8.5 ]
    cg_ranges = [ -0.65, 0.65 ]
    scan.POIs = [ 'ct', 'cg' ]
    if args.asimov: scan.PhysicsModelParameters = [ 'ct=1.0', 'cg=0.0' ]
    scan.PhysicsModelParameterRanges = [
        'ct={0},{1}'.format( ct_ranges[0], ct_ranges[1] ),
        'cg={0},{1}'.format( cg_ranges[0], cg_ranges[1] )
        ]
    scan.subDirectory = 'out/Scan_Top_{0}'.format(datestr)
    if args.highpt: scan.subDirectory = 'out/Scan_TopHighPt_{0}'.format(datestr)


    # ======================================
    # Determine the physics

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


    if args.nominal:
        scan.datacard = nominal_datacard(args)

    # In case of hbb, the minimizer settings need to be altered
    if args.combWithHbb or args.hbb:
        scan.default_minimizer_settings = False
        scan.extraOptions.extend([
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            ])

    elif args.lumiStudy:
        scan.datacard = LatestPaths.ws_combined_Top_lumiScalable
        scan.PhysicsModelParameters.append( 'lumiScale=8.356546' )
        scan.tags.append('lumiStudy')

    elif args.profiledTotalXS:
        scan.datacard = LatestPaths.ws_combined_Top_profiledTotalXS
        scan.tags.append('profiledTotalXS')

    elif args.fitOnlyNormalization:
        scan.datacard =(
            LatestPaths.ws_combined_Top_profiledTotalXS_fitOnlyNormalization
            if not args.highpt else
            LatestPaths.ws_combined_TopHighPt_profiledTotalXS_fitOnlyNormalization
            )
        scan.tags.append('fitOnlyNormalization')

    else:
        print 'Pass physics option'
        return

    scan.Run()


#____________________________________________________________________
@flag_as_option
def couplingScan_TopCtCb(args):
    print 'Now in couplingScan_TopCtCb'
    if not args.highpt:
        Commands.Warning('Probably no reason anymore to not use --highpt; setting --highpt to True')
        args.highpt = True

    TheoryCommands.SetPlotDir( 'plots_{0}_TopCtCb'.format(datestr) )

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

    scan.Run()


#____________________________________________________________________
@flag_as_option
def fastscanplotmock(args):

    TheoryCommands.SetPlotDir( 'plots_{0}_Top'.format(datestr) )

    print 'Running fastscanplotmock'

    fastscanFile = 'out/Scan_TopHighPt_Jan22_hzz_asimov_0/postfit_and_fastscan/higgsCombine_FASTSCAN_POSTFIT_ASIMOV_hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_TopHighPt_withTheoryUncertainties.MultiDimFit.mH125.root'

    mock = CombineToolWrapper.CombineScan()
    mock.POIs = [ 'ct', 'cg' ]
    mock.deltaNLLCutOff = 30.

    mock.PlotFastScan(fastscanFile)


########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'