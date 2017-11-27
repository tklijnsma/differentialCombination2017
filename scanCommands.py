#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

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

    parser.add_argument( '--scanCommands', default=False )
    class CustomAction(argparse.Action):
        def __init__(self, option_strings, dest, **kwargs):
            # Only 1 argument allowed (this is basically the "store_true" action)
            super(CustomAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            setattr( namespace, 'scanCommands', True )
            setattr( namespace, self.dest, True )
            # super(FooAction, self).__call__( parser, namespace, values, option_string=None )
            # setattr( namespace, self.dest, values )

    parser.add_argument( '--couplingScan',                             action=CustomAction )

    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument( '--yukawa', action='store_true' )
    group.add_argument( '--top',    action='store_true' )

    parser.add_argument( '--nominal',                                  action='store_true' )



########################################
# Methods
########################################    

def main( args ):

    expBinBoundaries = [ 0., 15., 30., 45., 85., 125. ]
    print 'Hardcoded binBoundaries for Yukawa:'
    print expBinBoundaries
    print ''

    TheoryCommands.SetPlotDir( 'plots_{0}_Yukawa'.format(datestr) )

    # ======================================
    # Scan


    #____________________________________________________________________
    if args.couplingScan:
        # Commands.TestMode(True)
        if not( args.yukawa or args.top ):
            print 'Use either --yukawa or --top'
            return

        RECREATE_POSTFIT = True
        # RECREATE_POSTFIT = False

        RECREATE_FASTSCAN = True
        # RECREATE_FASTSCAN = False

        # usePostfitWS = abspath( 'Scan_Top_Nov26_asimov/postfit_and_fastscan/POSTFIT_combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties.root' )
        # useFastscan  = abspath( 'Scan_Top_Nov26_asimov/postfit_and_fastscan/higgsCombine_FASTSCAN_POSTFIT_combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties.MultiDimFit.mH125.root' )

        RUN_QUICKLY = True
        # RUN_QUICKLY = False


        # ======================================
        # 

        baseContainer = Container()
        baseContainer.onBatch       = True
        baseContainer.nPointsPerJob = 20
        baseContainer.queue         = 'all.q'
        baseContainer.APPLY_FASTSCAN_FILTER = True

        if RUN_QUICKLY:
            Commands.Warning( 'Running with quick settings' )
            baseContainer.nPointsPerJob = 5
            baseContainer.queue         = 'short.q'

        if args.asimov:
            baseContainer.asimov = True

        if args.hzz:
            baseContainer.nPointsPerJob = 320
            baseContainer.queue         = 'short.q'

        if args.yukawa:
            baseContainer.nPoints       = 6400
            kappab_ranges = [ -15., 15. ]
            kappac_ranges = [ -35., 35. ]
            baseContainer.POIs = [ 'kappab', 'kappac' ]
            baseContainer.PhysicsModelParameters = [ 'kappab=1.0', 'kappac=1.0' ]
            baseContainer.PhysicsModelParameterRanges = [
                'kappab={0},{1}'.format( kappab_ranges[0], kappab_ranges[1] ),
                'kappac={0},{1}'.format( kappac_ranges[0], kappac_ranges[1] )
                ]
            baseContainer.subDirectory = 'Scan_Yukawa_{0}'.format(datestr)

        elif args.top:
            baseContainer.deltaNLLCutOff = 30.
            baseContainer.nPoints = 200**2
            # ct_ranges = [ -1., 2. ]
            # cg_ranges = [ -0.1, 0.2 ]
            ct_ranges = [ -8.5, 8.5 ]
            cg_ranges = [ -0.65, 0.65 ]
            baseContainer.POIs = [ 'ct', 'cg' ]
            baseContainer.PhysicsModelParameters = [ 'ct=1.0', 'cg=0.0' ]
            baseContainer.PhysicsModelParameterRanges = [
                'ct={0},{1}'.format( ct_ranges[0], ct_ranges[1] ),
                'cg={0},{1}'.format( cg_ranges[0], cg_ranges[1] )
                ]
            baseContainer.subDirectory = 'Scan_Top_{0}'.format(datestr)

        # Load settings into a scan baseContainer
        scan = CombineToolWrapper.CombineScan(baseContainer)

        if not RECREATE_POSTFIT:
            Commands.Warning( 'Using existing postfitWS {0}'.format(usePostfitWS) )
            scan.postfitWS = usePostfitWS

        if not RECREATE_FASTSCAN:
            Commands.Warning( 'Using existing fastscan {0}'.format(useFastscan) )
            scan.fastscanRootFile = useFastscan


        # ======================================
        # Determine the physics

        suffix = ''

        if args.hgg: suffix += '_hgg'
        if args.hzz: suffix += '_hzz'

        # ----------------------
        if args.yukawa:

            if args.nominal:
                datacard = LatestPaths.ws_combined_Yukawa
                if args.hgg:
                    datacard = LatestPaths.ws_hgg_Yukawa
                if args.hzz:
                    datacard = LatestPaths.ws_hzz_Yukawa
                scan.datacard = datacard

            else:
                print 'Pass physics option'
                return

        # ----------------------
        if args.top:

            if args.nominal:
                datacard = LatestPaths.ws_combined_Top
                if args.hgg:
                    datacard = LatestPaths.ws_hgg_Top
                if args.hzz:
                    datacard = LatestPaths.ws_hzz_Top
                scan.datacard = datacard

            else:
                print 'Pass physics option'
                return


        # ======================================
        # Run

        if args.asimov: suffix += '_asimov'

        # Make sure no other directory will be overwritten
        scan.subDirectory += suffix
        scan.subDirectory = Commands.AppendNumberToDirNameUntilItDoesNotExistAnymore(scan.subDirectory)
        scan.Run()


        

    #____________________________________________________________________
    # Literal copy of code in yukawaCommands only for easy viewing
    if False: # args.couplingBestfit_Yukawa: 

        # doFastscan = True
        # if args.notFastscan: doFastscan = False

        doFastscan = False
        if args.notFastscan: doFastscan = False
        if args.fastscan:    doFastscan = True

        ASIMOV = False
        if args.notAsimov: ASIMOV = False
        if args.asimov:    ASIMOV = True


        # ======================================
        # Flags

        SQUARE_DIST_POI_STEP              = True
        # SQUARE_DIST_POI_STEP              = False

        # LUMISTUDY = True
        LUMISTUDY = False

        # UNCORRELATED_THEORY_UNCERTAINTIES = True
        UNCORRELATED_THEORY_UNCERTAINTIES = False

        # NO_THEORY_UNCERTAINTIES           = True
        NO_THEORY_UNCERTAINTIES           = False


        PROFILE_TOTAL_XS                  = True
        # PROFILE_TOTAL_XS                  = False

        # ONEDIM_SCAN_TOTALXS               = True
        ONEDIM_SCAN_TOTALXS               = False

        FIT_ONLY_NORMALIZATION            = True
        # FIT_ONLY_NORMALIZATION            = False


        # INCLUDE_BR_COUPLING_DEPENDENCY    = True
        INCLUDE_BR_COUPLING_DEPENDENCY    = False

        # FIX_KAPPAV                        = True
        FIX_KAPPAV                        = False

        # MAX_KAPPAV_ONE                    = True
        MAX_KAPPAV_ONE                    = False


        # DO_BR_UNCERTAINTIES               = True
        DO_BR_UNCERTAINTIES               = False

        # FIT_RATIO_OF_BRS                  = True
        FIT_RATIO_OF_BRS                  = False

        # ONEDIM_SCAN_RATIO_OF_BRS          = True
        ONEDIM_SCAN_RATIO_OF_BRS          = False

        # DO_ONLY_ONE_KAPPA                 = True
        DO_ONLY_ONE_KAPPA                 = False

        # theKappa = 'kappab'
        theKappa = 'kappac'
        theOtherKappa = { 'kappab' : 'kappac', 'kappac' : 'kappab' }[theKappa]


        print
        print '(!) {0} asimov'.format( 'DOING' if ASIMOV else 'NOT DOING' )
        print '(!) {0} fastscan'.format( 'DOING' if doFastscan else 'NOT DOING' )
        print 'UNCORRELATED_THEORY_UNCERTAINTIES = ', UNCORRELATED_THEORY_UNCERTAINTIES
        print 'NO_THEORY_UNCERTAINTIES           = ', NO_THEORY_UNCERTAINTIES
        print
        print 'PROFILE_TOTAL_XS                  = ', PROFILE_TOTAL_XS
        print 'ONEDIM_SCAN_TOTALXS               = ', ONEDIM_SCAN_TOTALXS
        print 'FIT_ONLY_NORMALIZATION            = ', FIT_ONLY_NORMALIZATION
        print
        print 'INCLUDE_BR_COUPLING_DEPENDENCY    = ', INCLUDE_BR_COUPLING_DEPENDENCY
        print 'FIX_KAPPAV                        = ', FIX_KAPPAV
        print 'MAX_KAPPAV_ONE                    = ', MAX_KAPPAV_ONE
        print 'DO_BR_UNCERTAINTIES               = ', DO_BR_UNCERTAINTIES
        print
        print 'FIT_RATIO_OF_BRS                  = ', FIT_RATIO_OF_BRS
        print 'ONEDIM_SCAN_RATIO_OF_BRS          = ', ONEDIM_SCAN_RATIO_OF_BRS
        print
        print 'DO_ONLY_ONE_KAPPA                 = ', DO_ONLY_ONE_KAPPA
        if DO_ONLY_ONE_KAPPA:
            print '  Chosen kappa                      = {0}'.format(theKappa )
            print '  Other kappa                       = {0}'.format(theOtherKappa )

        if not INCLUDE_BR_COUPLING_DEPENDENCY and FIX_KAPPAV:
            Commands.ThrowError( 'INCLUDE_BR_COUPLING_DEPENDENCY == False and FIX_KAPPAV == True is not allowed' )
        if FIX_KAPPAV and MAX_KAPPAV_ONE:
            Commands.ThrowError( 'FIX_KAPPAV == True and MAX_KAPPAV_ONE == True is not allowed' )
        if not FIT_RATIO_OF_BRS and ONEDIM_SCAN_RATIO_OF_BRS:
            Commands.ThrowError( 'FIT_RATIO_OF_BRS == False and ONEDIM_SCAN_RATIO_OF_BRS == True is not allowed' )


        # ======================================
        # Set correct input

        datacard = LatestPaths.ws_combined_Yukawa
        if args.hgg:
            datacard = LatestPaths.ws_hgg_Yukawa
        if args.hzz:
            datacard = LatestPaths.ws_hzz_Yukawa

        if ( LUMISTUDY or UNCORRELATED_THEORY_UNCERTAINTIES or NO_THEORY_UNCERTAINTIES or PROFILE_TOTAL_XS ) and ( args.hzz or args.hgg ):
            print '[fixme] These studies not implemented for hgg or hzz'
            sys.exit()
        if LUMISTUDY:
            datacard = LatestPaths.ws_combined_Yukawa_lumiScalable
        if UNCORRELATED_THEORY_UNCERTAINTIES:
            datacard = LatestPaths.ws_combined_Yukawa_withUncorrelatedTheoryUncertainties
        if NO_THEORY_UNCERTAINTIES:
            datacard = LatestPaths.ws_combined_Yukawa_noTheoryUncertainties
        if PROFILE_TOTAL_XS:
            datacard = LatestPaths.ws_combined_Yukawa_profiledTotalXS
        if INCLUDE_BR_COUPLING_DEPENDENCY:
            datacard = LatestPaths.ws_combined_Yukawa_couplingDependentBR
        if FIT_RATIO_OF_BRS:
            datacard = LatestPaths.ws_combined_ratioOfBRs
        if FIT_ONLY_NORMALIZATION:
            datacard = LatestPaths.ws_combined_Yukawa_profiledTotalXS_fitOnlyNormalization

        # ======================================
        # Set some job specifics (ranges, number of points)

        jobDirectory = 'Scan_Yukawa_{0}'.format( datestr )
        if args.hgg: jobDirectory += '_hgg'
        if args.hzz: jobDirectory += '_hzz'

        kappab_ranges = [ -15., 15. ]
        kappac_ranges = [ -35., 35. ]


        if doFastscan:
            jobDirectory += '_fastscan'
            nPoints = 6400
            nPointsPerJob = 800
            queue = 'short.q'
        else:
            nPoints = 6400
            nPointsPerJob = 20
            queue = 'all.q'
            if args.hzz:
                nPointsPerJob = 320
                queue = 'short.q'


        # print '\n' + 'WARNING ' * 7
        # print 'TEMPORARY SETTINGS  - REMOVE THESE'
        # kappab_ranges = [ -8., 8. ]
        # kappac_ranges = [ -20., 20. ]
        # nPoints = 60*60
        # SQUARE_DIST_POI_STEP = False


        if DO_ONLY_ONE_KAPPA:
            nPoints       = 39
            nPointsPerJob = 3
            queue         = 'short.q'

        if ONEDIM_SCAN_RATIO_OF_BRS or ONEDIM_SCAN_TOTALXS:
            nPoints       = 50
            nPointsPerJob = 5
            queue         = 'short.q'


        # ======================================
        # Construct the fit command and process flags

        # --------------------
        # Setting POIs

        POIs = [ 'kappab', 'kappac' ]
        if DO_ONLY_ONE_KAPPA:
            POIs = [ theKappa ]
        if ONEDIM_SCAN_RATIO_OF_BRS:
            POIs = [ 'ratio_BR_hgg_hzz' ]
        if ONEDIM_SCAN_TOTALXS:
            POIs = [ 'r_totalXS' ]
        # if FIT_ONLY_NORMALIZATION:
        #     POIs = [ 'globalTotalXSmodifier' ]
            

        # --------------------
        # Setting physicsModelParameters

        physicsModelParameters = [ 'kappab=1.0', 'kappac=1.0' ]
        if LUMISTUDY:
            physicsModelParameters.append( 'lumiScale=8.356546' )
        if INCLUDE_BR_COUPLING_DEPENDENCY:
            physicsModelParameters.append( 'kappa_V=0.99' )
        if PROFILE_TOTAL_XS:
            physicsModelParameters.append( 'r_totalXS=1.0' )


        # --------------------
        # Setting physicsModelParameterRanges

        physicsModelParameterRanges = [
            'kappab={0},{1}'.format( kappab_ranges[0], kappab_ranges[1] ),
            'kappac={0},{1}'.format( kappac_ranges[0], kappac_ranges[1] )
            ]
        if INCLUDE_BR_COUPLING_DEPENDENCY and not FIX_KAPPAV:
            if MAX_KAPPAV_ONE:
                physicsModelParameterRanges.append( 'kappa_V=-100.0,1.0' )
            else:
                physicsModelParameterRanges.append( 'kappa_V=-100.0,100.0' )

        if ONEDIM_SCAN_RATIO_OF_BRS:
            # Actually overwrite list contents
            physicsModelParameterRanges = [ 'ratio_BR_hgg_hzz=0.0,0.25' ]

        if ONEDIM_SCAN_TOTALXS:
            # Actually overwrite list contents
            physicsModelParameterRanges = [ 'r_totalXS=0.0,2.0' ]

        # if FIT_ONLY_NORMALIZATION:
        #     physicsModelParameterRanges = [ 'globalTotalXSmodifier=0.0,2.0' ]


        # --------------------
        # Specify floating and frozen nuisances

        floatNuisances  = []
        freezeNuisances = []

        if INCLUDE_BR_COUPLING_DEPENDENCY:
            if FIX_KAPPAV:
                freezeNuisances.append( 'kappa_V' )
            else:
                floatNuisances.append( 'kappa_V' )

        if DO_ONLY_ONE_KAPPA:
            floatNuisances.append( theOtherKappa )


        # --------------------
        # Construct extraOptions

        extraOptions = [
            '-P ' + ' -P '.join(POIs),
            '--setPhysicsModelParameters '      + ','.join(physicsModelParameters),
            '--setPhysicsModelParameterRanges ' + ':'.join(physicsModelParameterRanges),
            ]

        if SQUARE_DIST_POI_STEP:
            extraOptions.append( '--squareDistPoiStep' )

        if len(floatNuisances) > 0:
            extraOptions.append( '--floatNuisances ' + ','.join(floatNuisances) )
        if len(freezeNuisances) > 0:
            extraOptions.append( '--freezeNuisances ' + ','.join(freezeNuisances) )

        
        # --------------------
        # Compile list of variables to save

        variablesToSave = []
        variablesToSave.extend( Commands.ListSet( datacard, 'yieldParameters' ) )
        variablesToSave.extend( [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ] )
        if INCLUDE_BR_COUPLING_DEPENDENCY:
            variablesToSave.extend( Commands.ListSet( datacard, 'hgg_yieldParameters' ) )
            variablesToSave.extend( Commands.ListSet( datacard, 'hzz_yieldParameters' ) )
            variablesToSave.extend( Commands.ListSet( datacard, 'BRvariables' ) )
        if PROFILE_TOTAL_XS:
            variablesToSave.extend([ 'r_totalXS', 'totalXSmodifier', 'totalXS_SM', 'totalXS' ])
        if FIT_RATIO_OF_BRS:
            variablesToSave.extend([ 'hgg_ratioBRmodifier', 'hzz_ratioBRmodifier', 'ratio_BR_hgg_hzz' ])
        extraOptions.append( '--saveSpecifiedFunc ' + ','.join(variablesToSave) )


        # ======================================
        # Appropriately name scan, create jobDirectoy and submit command

        if DO_ONLY_ONE_KAPPA:
            jobDirectory += '_oneKappa_{0}'.format(theKappa)
        if LUMISTUDY:
            jobDirectory += '_lumiStudy'
        if PROFILE_TOTAL_XS:
            jobDirectory += '_profiledTotalXS'
            if ONEDIM_SCAN_TOTALXS:
                jobDirectory += '_onedimTotalXSScan'
            if FIT_ONLY_NORMALIZATION:
                jobDirectory += '_fitOnlyNormalization'
        if UNCORRELATED_THEORY_UNCERTAINTIES:
            jobDirectory += '_uncorrelatedTheoryUncertainties'
        if NO_THEORY_UNCERTAINTIES:
            jobDirectory += '_noTheoryUncertainties'
        if INCLUDE_BR_COUPLING_DEPENDENCY:
            jobDirectory += '_couplingDependentBR'
            if FIX_KAPPAV:
                jobDirectory += '_fixedKappaV'
            elif MAX_KAPPAV_ONE:
                jobDirectory += '_kappaVMaxOne'
        if FIT_RATIO_OF_BRS:
            jobDirectory += '_ratioOfBRs'
            if ONEDIM_SCAN_RATIO_OF_BRS:
                jobDirectory += '_onedimRatioScan'  

        if ASIMOV:
            jobDirectory += '_asimov'

        jobDirectory = Commands.AppendNumberToDirNameUntilItDoesNotExistAnymore( jobDirectory )
        if Commands.IsTestMode(): print '\nWould now create new directory: {0}'.format( basename(jobDirectory) )

        Commands.MultiDimCombineTool(
            datacard,
            nPoints       = nPoints,
            nPointsPerJob = nPointsPerJob,
            queue         = queue,
            notOnBatch    = False,
            jobDirectory  = jobDirectory,
            fastscan      = doFastscan,
            asimov        = ASIMOV,
            jobPriority   = 0,
            extraOptions  = extraOptions
            )






########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'