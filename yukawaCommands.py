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

    parser.add_argument( '--yukawaCommands', default=False )
    class CustomAction(argparse.Action):
        def __init__(self, option_strings, dest, **kwargs):
            # Only 1 argument allowed (this is basically the "store_true" action)
            super(CustomAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            # Make sure the parser knows a combineCommand was entered
            setattr( namespace, 'yukawaCommands', True )
            setattr( namespace, self.dest, True )
            # super(FooAction, self).__call__( parser, namespace, values, option_string=None )
            # setattr( namespace, self.dest, values )

    parser.add_argument( '--createDerivedTheoryFiles_YukawaQuarkInduced',         action=CustomAction )
    parser.add_argument( '--createDerivedTheoryFiles_YukawaGluonInduced',         action=CustomAction )
    parser.add_argument( '--createDerivedTheoryFiles_YukawaQuarkInducedScaled',   action=CustomAction )

    parser.add_argument( '--mergeGluonInducedWithQuarkInduced',                   action=CustomAction )
    parser.add_argument( '--CorrelationMatrices_Yukawa',                          action=CustomAction )
    parser.add_argument( '--couplingT2WS_Yukawa',                                 action=CustomAction )
    parser.add_argument( '--couplingBestfit_Yukawa',                              action=CustomAction )

    # Plotting
    parser.add_argument( '--plot_onedimRatioOfBRs_Yukawa',                        action=CustomAction )
    parser.add_argument( '--plot_onedimTotalXS_Yukawa',                           action=CustomAction )
    parser.add_argument( '--OneKappaScanPlot_Yukawa',                             action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa',                          action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa_TheoryCrossCheck',         action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa_BRdependencyComparison',   action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa_highLumiStudy',            action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa_ProfiledTotalXS',          action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa_atfchi2',                  action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa_RatioOfBRs',               action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa_onlyNormalization',        action=CustomAction )

    parser.add_argument( '--coupling2Dplot_Yukawa',                               action=CustomAction )
    parser.add_argument( '--checkWSParametrization_Yukawa',                       action=CustomAction )
    parser.add_argument( '--combinationAndContour_Yukawa',                       action=CustomAction )


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
    # Dealing with theory IO

    if args.createDerivedTheoryFiles_YukawaQuarkInduced:
        TheoryFileInterface.CreateDerivedTheoryFiles_YukawaQuarkInduced(
            verbose = True,
            )

    if args.createDerivedTheoryFiles_YukawaGluonInduced:
        TheoryFileInterface.CreateDerivedTheoryFiles_YukawaGluonInduced(
            verbose = True,
            mainCrossSection = 'matched',
            )

    if args.createDerivedTheoryFiles_YukawaQuarkInducedScaled:
        TheoryFileInterface.ScaleQuarkInduced(
            LatestPaths.derivedTheoryFiles_YukawaQuarkInduced,
            verbose = True,
            )

    if args.mergeGluonInducedWithQuarkInduced:
        TheoryFileInterface.SumGluonAndQuarkFiles(
            gI_theoryFileDir = LatestPaths.derivedTheoryFiles_YukawaGluonInduced,
            qI_theoryFileDir = LatestPaths.derivedTheoryFiles_YukawaQuarkInducedScaled,
            verbose = True
            )


    # ======================================
    # Dealing with theory uncertainties

    if args.CorrelationMatrices_Yukawa:

        variationFiles = TheoryFileInterface.FileFinder(
            directory = LatestPaths.derivedTheoryFiles_YukawaSummed,
            kappab = 1, kappac = 1
            )

        variations = [
            TheoryFileInterface.ReadDerivedTheoryFile( variationFile, returnContainer=True )
                for variationFile in variationFiles ]

        CorrelationMatrices.GetCorrelationMatrix(
            variations,
            makeScatterPlots          = False,
            makeCorrelationMatrixPlot = True,
            outname                   = 'corrMat_theory',
            verbose                   = True,
            )

        variations_expbinning = deepcopy(variations)
        for variation in variations_expbinning:
            TheoryCommands.RebinDerivedTheoryContainer( variation, expBinBoundaries )

        CorrelationMatrices.GetCorrelationMatrix(
            variations_expbinning,
            makeScatterPlots          = True,
            makeCorrelationMatrixPlot = True,
            outname                   = 'corrMat_exp',
            verbose                   = True,
            )


    if args.couplingT2WS_Yukawa:

        # ======================================
        # Flags

        INCLUDE_THEORY_UNCERTAINTIES      = True
        # INCLUDE_THEORY_UNCERTAINTIES      = False

        # UNCORRELATED_THEORY_UNCERTAINTIES = True
        UNCORRELATED_THEORY_UNCERTAINTIES = False

        # MAKELUMISCALABLE                  = True
        MAKELUMISCALABLE                  = False

        # INCLUDE_BR_COUPLING_DEPENDENCY    = True
        INCLUDE_BR_COUPLING_DEPENDENCY    = False

        PROFILE_TOTAL_XS                  = True
        # PROFILE_TOTAL_XS                  = False

        FIT_ONLY_NORMALIZATION            = True
        # FIT_ONLY_NORMALIZATION            = False

        # DO_BR_UNCERTAINTIES               = True
        DO_BR_UNCERTAINTIES               = False

        # FIT_RATIO_OF_BRS                  = True
        FIT_RATIO_OF_BRS                  = False


        # ======================================
        # Specify needed input

        datacard = LatestPaths.card_combined_ggHxH_PTH
        if args.hgg:
            datacard = LatestPaths.card_hgg_ggHxH_PTH
        if args.hzz:
            datacard = LatestPaths.card_hzz_ggHxH_PTH

        TheoryFileInterface.SetFileFinderDir( LatestPaths.derivedTheoryFiles_YukawaSummed )

        if INCLUDE_THEORY_UNCERTAINTIES:
            correlationMatrix   = LatestPaths.correlationMatrix_Yukawa
            theoryUncertainties = LatestPaths.theoryUncertainties_Yukawa
            if UNCORRELATED_THEORY_UNCERTAINTIES:
                correlationMatrix = LatestPaths.correlationMatrix_Yukawa_Uncorrelated


        # ======================================
        # Run t2ws

        extraOptions = [
            '--PO verbose=2',
            '--PO \'higgsMassRange=123,127\'',
            '--PO linearTerms=True',
            '--PO splitggH=True',
            ]

        if args.hzz:
            extraOptions.append( '--PO isOnlyHZZ=True' )
        if args.hgg:
            extraOptions.append( '--PO isOnlyHgg=True' )

        extraOptions.append(
            '--PO binBoundaries={0}'.format( ','.join([ str(b) for b in expBinBoundaries ]) )
            )

        extraOptions.append(
            '--PO SM=[kappab=1,kappac=1,file={0}]'.format(
                TheoryFileInterface.FileFinder( kappab=1, kappac=1, muR=1, muF=1, Q=1, expectOneFile=True )
                )
            )

        possibleTheories = []
        for kappab in [ -2, -1, 0, 1, 2 ]:
            for kappac in [ -10, -5, 0, 1, 5, 10 ]:
                if ( kappab == 1 and kappac == 1 ) or ( kappab == 0 and kappac == 0 ): continue
                else:
                    possibleTheories.append(
                        '--PO theory=[kappab={0},kappac={1},file={2}]'.format(
                            kappab, kappac,
                            TheoryFileInterface.FileFinder( kappab=kappab, kappac=kappac, muR=1, muF=1, Q=1, expectOneFile=True )
                            )
                        )

        # Sample only a few needed theories
        import random
        random.seed(1002)
        extraOptions.extend( random.sample( possibleTheories, 6 ) )


        suffix = 'Yukawa'

        if INCLUDE_THEORY_UNCERTAINTIES:
            extraOptions.append( '--PO correlationMatrix={0}'.format( correlationMatrix ) )
            extraOptions.append( '--PO theoryUncertainties={0}'.format( theoryUncertainties ) )
            suffix += '_withTheoryUncertainties' if not UNCORRELATED_THEORY_UNCERTAINTIES else '_withUncorrelatedTheoryUncertainties'
        else:
            suffix += '_noTheoryUncertainties'

        if MAKELUMISCALABLE:
            extraOptions.append( '--PO lumiScale=True' )
            suffix += '_lumiScale'

        if INCLUDE_BR_COUPLING_DEPENDENCY:
            extraOptions.append( '--PO FitBR=True' )
            suffix += '_couplingDependentBR'
            if DO_BR_UNCERTAINTIES:
                extraOptions.append( '--PO DoBRUncertainties=True' )
                suffix += '_withBRUnc'

        if PROFILE_TOTAL_XS:
            extraOptions.append( '--PO ProfileTotalXS=True' )
            suffix += '_profiledTotalXS'

        if FIT_RATIO_OF_BRS:
            extraOptions.append( '--PO FitRatioOfBRs=True' )
            suffix += '_ratioOfBRs'

        if FIT_ONLY_NORMALIZATION:
            extraOptions.append( '--PO FitOnlyNormalization=True' )
            suffix += '_fitOnlyNormalization'

        Commands.BasicT2WSwithModel(
            datacard,
            'CouplingModel.py',
            suffix = suffix,
            extraOptions = extraOptions,
            )


    # ======================================
    # Scan

    if args.couplingBestfit_Yukawa:

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




    # 
    # 
    ################################################################################
    ################################################################################
    #        P L O T S
    ################################################################################
    ################################################################################
    # 
    # 


    #____________________________________________________________________
    if args.OneKappaScanPlot_Yukawa:

        scans = [
            ( 'kappac', LatestPaths.scan_combined_Yukawa_oneKappa_kappac ),
            ( 'kappac', LatestPaths.scan_combined_Yukawa_oneKappa_kappac_asimov ),
            ( 'kappab', LatestPaths.scan_combined_Yukawa_oneKappa_kappab ),
            ( 'kappab', LatestPaths.scan_combined_Yukawa_oneKappa_kappab_asimov ),
            ]

        for kappa, scanDir in scans:

            container = Container()
            container.kappa = kappa
            container.x, container.y = TheoryCommands.BasicReadScan( glob( scanDir + '/*.root' ), kappa )

            PlotCommands.PlotMultipleScans(
                [ container ],
                xTitle   = kappa,
                yTitle   = '2#Delta NLL',
                yMax     = 5.,
                plotname = 'oneKappaScan_{0}_{1}'.format( kappa, basename(scanDir).replace('/','') ),
                draw1sigmaline = True,
                translateToChi2 = True,
                printUncertainties = True,
                )



    #____________________________________________________________________
    if args.plot_onedimRatioOfBRs_Yukawa:

        scans = [
            LatestPaths.scan_combined_Yukawa_ratioOfBRs_onedimRatioScan,
            LatestPaths.scan_combined_Yukawa_ratioOfBRs_onedimRatioScan_asimov
            ]

        for scanDir in scans:

            container = Container()
            container.x, container.y = TheoryCommands.BasicReadScan( glob( scanDir + '/*.root' ), 'ratio_BR_hgg_hzz' )
            container.color = 1

            xBestfit      = container.x[ container.y.index(min(container.y)) ]
            lineAtBestfit = Container()
            lineAtBestfit.line = ROOT.TLine( xBestfit, -999, xBestfit, 999 )
            lineAtBestfit.line.SetLineWidth(2)
            lineAtBestfit.line.SetLineColor(2)

            xSM           = LatestPaths.SM_ratio_of_BRs
            lineAtSM      = Container()
            lineAtSM.line = ROOT.TLine( xSM, -999, xSM, 999 )
            lineAtSM.line.SetLineWidth(2)
            lineAtSM.line.SetLineColor(9)

            PlotCommands.PlotMultipleScans(
                [ container, lineAtSM, lineAtBestfit ],
                xTitle   = 'BR_{H#rightarrow #gamma#gamma} / BR_{H#rightarrow ZZ}',
                yTitle   = '2#DeltaNLL',
                xMin     = 0.0525,
                xMax     = 0.1425,
                yMax     = 5.,
                plotname = 'onedimRatioOfBRs_{0}'.format( basename(scanDir).replace('/','') ),
                draw1sigmaline = True,
                translateToChi2 = True
                )


    #____________________________________________________________________
    if args.plot_onedimTotalXS_Yukawa:

        scans = [
            LatestPaths.scan_combined_Yukawa_profiledTotalXS_onedimTotalXSScan,
            LatestPaths.scan_combined_Yukawa_profiledTotalXS_onedimTotalXSScan_asimov
            ]

        for scanDir in scans:

            container = Container()
            container.x, container.y = TheoryCommands.BasicReadScan( glob( scanDir + '/*.root' ), 'r_totalXS' )
            container.color = 1

            xBestfit      = container.x[ container.y.index(min(container.y)) ]
            lineAtBestfit = Container()
            lineAtBestfit.line = ROOT.TLine( xBestfit, -999, xBestfit, 999 )
            lineAtBestfit.line.SetLineWidth(2)
            lineAtBestfit.line.SetLineColor(2)

            unc = PhysicsCommands.FindMinimaAndErrors( container.x, container.y, returnContainer=True )
            uncLines = []
            for x in [ unc.leftBound, unc.rightBound ]:
                uncLine = Container()
                uncLine.line = ROOT.TLine( x, 0., x, 1.0 )
                uncLine.line.SetLineWidth(1)
                uncLine.line.SetLineColor(2)
                uncLines.append(uncLine)

            xSM           = 1.0
            lineAtSM      = Container()
            lineAtSM.line = ROOT.TLine( xSM, -999, xSM, 999 )
            lineAtSM.line.SetLineWidth(2)
            lineAtSM.line.SetLineColor(9)

            PlotCommands.PlotMultipleScans(
                [ container, lineAtSM, lineAtBestfit ] + uncLines,
                xTitle   = '#mu_{#sigma_{tot}}',
                yTitle   = '2#DeltaNLL',
                xMin     = 0.6,
                xMax     = 1.4,
                yMax     = 5.,
                plotname = 'onedimTotalXS_{0}'.format( basename(scanDir).replace('/','') ),
                draw1sigmaline = True,
                translateToChi2 = True,
                printUncertainties = True
                )


    #____________________________________________________________________
    if args.coupling2Dplot_Yukawa:

        # datacard  = combinedDatacard
        # rootfiles = combinedScanFiles
        # if args.hgg:
        #     datacard  = hggDatacard
        #     rootfiles = hggScanFiles
        # if args.hzz:
        #     datacard  = hzzDatacard
        #     rootfiles = hzzScanFiles

        rootfiles = glob( LatestPaths.scan_combined_Yukawa + '/*.root' )
        # # rootfiles = glob( 'testscan_Nov22' + '/*.root' )
        # rootfiles = glob( 'Scan_Yukawa_Nov23_0' + '/*.root' )

        res = TheoryCommands.PlotCouplingScan2D(
            rootfiles,
            xCoupling = 'kappac',
            yCoupling = 'kappab',
            xMin = -35., xMax = 35.,
            yMin = -13., yMax = 13.,
            verbose = True
            )

        print '\nBest fit:'
        print res.xCoupling, '=', res.xBestfit
        print res.yCoupling, '=', res.yBestfit



    if args.couplingContourPlot_Yukawa:

        ASIMOV = False
        if args.asimov: ASIMOV = True


        if ASIMOV: print 'Warning: plotting ASIMOV scans'

        xCoupling = 'kappac'
        yCoupling = 'kappab'
        titles = { 'kappac': '#kappa_{c}', 'kappab' : '#kappa_{b}' }

        TH2FsToPlot = []

        combined_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa ))
        if ASIMOV: combined_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_asimov))
        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name = 'combined'
        combined.title = 'Combination'
        TH2FsToPlot.append(combined)

        hgg_rootfiles      = glob( '{0}/*.root'.format( LatestPaths.scan_hgg_Yukawa ))
        if ASIMOV: hgg_rootfiles      = glob( '{0}/*.root'.format( LatestPaths.scan_hgg_Yukawa_asimov))
        hgg = TheoryCommands.GetTH2FromListOfRootFiles(
            hgg_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        hgg.color = 2
        hgg.name = 'hgg'
        hgg.title = 'H #rightarrow #gamma#gamma'
        TH2FsToPlot.append(hgg)

        hzz_rootfiles      = glob( '{0}/*.root'.format( LatestPaths.scan_hzz_Yukawa ))
        if ASIMOV: hzz_rootfiles      = glob( '{0}/*.root'.format( LatestPaths.scan_hzz_Yukawa_asimov))
        hzz = TheoryCommands.GetTH2FromListOfRootFiles(
            hzz_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        hzz.color = 4
        hzz.name = 'hzz'
        hzz.title = 'H #rightarrow ZZ'
        TH2FsToPlot.append(hzz)


        PlotCommands.BasicMixedContourPlot(
            TH2FsToPlot,
            xMin = -35.,
            xMax = 35.,
            yMin = -13.,
            yMax = 13.,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_perDecayChannel' + ( '_asimov' if ASIMOV else '' ),
            x_SM      = 1.,
            y_SM      = 1.,
            plotIndividualH2s = True,
            )


    #____________________________________________________________________
    if args.couplingContourPlot_Yukawa_RatioOfBRs:

        xCoupling = 'kappac'
        yCoupling = 'kappab'
        titles = { 'kappac': '#kappa_{c}', 'kappab' : '#kappa_{b}' }

        containers = []

        combined_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_asimov ) )
        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name  = 'regular'
        combined.title = 'SM BR'
        containers.append(combined)

        ratioOfBRs_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_ratioOfBRs_asimov ) )
        ratioOfBRs = TheoryCommands.GetTH2FromListOfRootFiles(
            ratioOfBRs_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        ratioOfBRs.color = 4
        ratioOfBRs.name  = 'ratioOfBRs'
        ratioOfBRs.title = 'floating BR_{H#gamma#gamma}/BR_{HZZ}'
        containers.append(ratioOfBRs)

        PlotCommands.BasicMixedContourPlot(
            containers,
            xMin = -35.,
            xMax = 35.,
            yMin = -13.,
            yMax = 13.,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_ratioOfBRs',
            x_SM      = 1.,
            y_SM      = 1.,
            plotIndividualH2s = True,
            )


    #____________________________________________________________________
    if args.couplingContourPlot_Yukawa_ProfiledTotalXS:

        xCoupling = 'kappac'
        yCoupling = 'kappab'
        titles = { 'kappac': '#kappa_{c}', 'kappab' : '#kappa_{b}' }

        containers = []

        combined_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_asimov ) )
        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name  = 'regular'
        combined.title = 'Nominal'
        containers.append(combined)

        profiledTotalXS_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_profiledTotalXS_asimov ) )
        profiledTotalXS = TheoryCommands.GetTH2FromListOfRootFiles(
            profiledTotalXS_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        profiledTotalXS.color = 4
        profiledTotalXS.name = 'profiledTotalXS'
        profiledTotalXS.title = '#sigma_{tot} profiled'
        containers.append(profiledTotalXS)

        PlotCommands.BasicMixedContourPlot(
            containers,
            xMin = -35.,
            xMax = 35.,
            yMin = -13.,
            yMax = 13.,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_profiledTotalXS',
            x_SM      = 1.,
            y_SM      = 1.,
            plotIndividualH2s = True,
            )


    #____________________________________________________________________
    if args.couplingContourPlot_Yukawa_TheoryCrossCheck:

        xCoupling = 'kappac'
        yCoupling = 'kappab'
        titles = { 'kappac': '#kappa_{c}', 'kappab' : '#kappa_{b}' }

        containers = []

        combined_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_asimov ) )
        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name  = 'regular'
        combined.title = 'Nominal'
        containers.append(combined)

        noTheoryUnc_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_noTheoryUncertainties_asimov ) )
        noTheoryUnc = TheoryCommands.GetTH2FromListOfRootFiles(
            noTheoryUnc_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        noTheoryUnc.color = 4
        noTheoryUnc.name = 'noTheoryUnc'
        noTheoryUnc.title = 'No unc.'
        containers.append(noTheoryUnc)

        uncTheoryUnc_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_uncorrelatedTheoryUncertainties_asimov ) )
        uncTheoryUnc = TheoryCommands.GetTH2FromListOfRootFiles(
            uncTheoryUnc_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        uncTheoryUnc.color = 2
        uncTheoryUnc.name = 'uncTheoryUnc'
        uncTheoryUnc.title = 'Uncorr. unc.'
        containers.append(uncTheoryUnc)

        PlotCommands.BasicMixedContourPlot(
            containers,
            xMin = -25.,
            xMax = 25.,
            yMin = -9.,
            yMax = 9.,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_theoryCrossCheck',
            x_SM      = 1.,
            y_SM      = 1.,
            plotIndividualH2s = True,
            )



    #____________________________________________________________________
    if args.couplingContourPlot_Yukawa_BRdependencyComparison:

        containers = []

        xCoupling = 'kappac'
        yCoupling = 'kappab'
        titles = { 'kappac': '#kappa_{c}', 'kappab' : '#kappa_{b}' }

        combined_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_asimov) )
        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name = 'regularBR'
        combined.title = 'BR constant'
        containers.append( combined )

        # scalingBR_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_couplingDependentBR_asimov ) )
        # scalingBR = TheoryCommands.GetTH2FromListOfRootFiles(
        #     scalingBR_rootfiles,
        #     xCoupling,
        #     yCoupling,
        #     verbose   = False,
        #     )
        # scalingBR.color = 2
        # scalingBR.name = 'scalingBR'
        # scalingBR.title = 'BR(#kappa_{t}, #kappa_{V})'
        # containers.append( scalingBR )

        scalingBRfixedKappaV_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_couplingDependentBR_fixedKappaV_asimov ) )
        scalingBRfixedKappaV = TheoryCommands.GetTH2FromListOfRootFiles(
            scalingBRfixedKappaV_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        scalingBRfixedKappaV.color = 4
        scalingBRfixedKappaV.name = 'scalingBRfixedKappaV'
        scalingBRfixedKappaV.title = 'BR(#kappa_{b}, #kappa_{c})' # + ' (#kappa_{V} fixed)'
        containers.append( scalingBRfixedKappaV )


        # scalingBRkappaVMaxOne_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_couplingDependentBR_kappaVMaxOne_asimov ) )
        # scalingBRkappaVMaxOne = TheoryCommands.GetTH2FromListOfRootFiles(
        #     scalingBRkappaVMaxOne_rootfiles,
        #     xCoupling,
        #     yCoupling,
        #     verbose   = False,
        #     )
        # scalingBRkappaVMaxOne.color = 8
        # scalingBRkappaVMaxOne.name = 'scalingBRkappaVMaxOne'
        # scalingBRkappaVMaxOne.title = 'BR(#kappa_{t}, #kappa_{V}#leq1)'
        # containers.append( scalingBRkappaVMaxOne )


        PlotCommands.BasicMixedContourPlot(
            containers,
            xMin = -25.,
            xMax = 25.,
            yMin = -9.,
            yMax = 13.,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_BRcouplingDependency',
            x_SM      = 1.,
            y_SM      = 1.,
            plotIndividualH2s = True,
            )

    #____________________________________________________________________
    if args.couplingContourPlot_Yukawa_highLumiStudy:

        xCoupling = 'kappac'
        yCoupling = 'kappab'
        titles = { 'kappac': '#kappa_{c}', 'kappab' : '#kappa_{b}' }

        containers = []

        combined_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_asimov ) )
        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name  = 'regular'
        combined.title = '35.9 fb^{-1}'
        containers.append(combined)

        highLumi_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_lumiStudy_asimov ) )
        highLumi = TheoryCommands.GetTH2FromListOfRootFiles(
            highLumi_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        highLumi.color = 4
        highLumi.name = 'highLumi'
        highLumi.title = '300 fb^{-1}'
        containers.append(highLumi)


        PlotCommands.BasicMixedContourPlot(
            containers,
            xMin = -25.,
            xMax = 25.,
            yMin = -9.,
            yMax = 9.,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_highLumiStudy',
            x_SM      = 1.,
            y_SM      = 1.,
            plotIndividualH2s = True,
            )



    #____________________________________________________________________
    if args.couplingContourPlot_Yukawa_onlyNormalization:

        xCoupling = 'kappac'
        yCoupling = 'kappab'
        titles = { 'kappac': '#kappa_{c}', 'kappab' : '#kappa_{b}' }

        containers = []

        combined_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_asimov ) )
        combined = TheoryCommands.GetTH2FromListOfRootFiles( combined_rootfiles, xCoupling, yCoupling, verbose   = False, )
        combined.color = 1
        combined.name  = 'combined'
        combined.title = 'Nominal'
        containers.append(combined)

        onlyNormalization_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_fitOnlyNormalization_asimov ) )
        onlyNormalization = TheoryCommands.GetTH2FromListOfRootFiles( onlyNormalization_rootfiles, xCoupling, yCoupling, verbose   = False, )
        onlyNormalization.color = 2
        onlyNormalization.name  = 'onlyNormalization'
        onlyNormalization.title = 'Only normalization'
        containers.append(onlyNormalization)

        PlotCommands.BasicMixedContourPlot(
            containers,
            xMin = -35.,
            xMax = 35.,
            yMin = -13.,
            yMax = 13.,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_onlyNormalization',
            x_SM      = 1.,
            y_SM      = 1.,
            plotIndividualH2s = True,
            )


    #____________________________________________________________________
    if args.couplingContourPlot_Yukawa_atfchi2:

        xCoupling = 'kappac'
        yCoupling = 'kappab'
        titles = { 'kappac': '#kappa_{c}', 'kappab' : '#kappa_{b}' }

        containers = []

        # combined_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_asimov ) )
        combined_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_noTheoryUncertainties_asimov ) )
        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 4
        combined.name  = 'regular'
        combined.title = 'Full likelihood (no theor. unc.)'
        containers.append(combined)


        # Get relevant info from a canvas that has the histogram and the bestfit point
        atfchi2_rootFp = ROOT.TFile.Open( 'plots_Nov14_Yukawa/AfterTheFactChi2Fit.root' )
        canvas = atfchi2_rootFp.Get('ctc')

        # Get TH2F
        atfchi2_H = canvas.GetPrimitive('H2')
        ROOT.SetOwnership( atfchi2_H, False )

        # Get values of x and y for bestfit
        bestfitpoint_Tg = canvas.GetPrimitive('bestfitpoint')
        x_Double = ROOT.Double(0)
        y_Double = ROOT.Double(0)
        bestfitpoint_Tg.GetPoint( 0, x_Double, y_Double )
        xBestfit = float(x_Double)
        yBestfit = float(y_Double)

        atfchi2_rootFp.Close()

        # Construct container
        atfchi2 = Container()
        atfchi2.H2       = atfchi2_H
        atfchi2.xBestfit = xBestfit
        atfchi2.yBestfit = yBestfit
        atfchi2.color    = 2
        atfchi2.name     = 'atfchi2'
        atfchi2.title    = 'Simple #chi^{2} fit'
        containers.append( atfchi2 )


        PlotCommands.BasicMixedContourPlot(
            containers,
            xMin = -25.,
            xMax = 25.,
            yMin = -9.,
            yMax = 9.,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_atfchi2',
            x_SM      = 1.,
            y_SM      = 1.,
            plotIndividualH2s = True,
            )

    #____________________________________________________________________
    if args.checkWSParametrization_Yukawa:

        # plotCombination = False
        plotCombination = True

        # DrawExperimentalBinLines = False
        DrawExperimentalBinLines = True


        # # wsToCheck =  LatestPaths.ws_combined_unsplit_yukawa
        # wsToCheck =  LatestPaths.ws_combined_unsplit_yukawa_onlyGluonInduced

        # # theoryDir = LatestPaths.derivedTheoryFiles_YukawaSummed
        # theoryDir = LatestPaths.derivedTheoryFiles_YukawaGluonInduced

        # wsToCheck = LatestPaths.ws_onlyhzz_split_yukawa
        wsToCheck = LatestPaths.ws_combined_Yukawa
        theoryDir = LatestPaths.derivedTheoryFiles_YukawaSummed
        

        wsParametrization = WSParametrization( wsToCheck )
        
        TheoryFileInterface.SetFileFinderDir( theoryDir )
        yukawaDerivedTheoryFiles = TheoryFileInterface.FileFinder( muF=1, muR=1, Q=1 )


        # Select subset of theory files (otherwise plot gets overcrowded)
        import random
        random.seed(1000)
        yukawaDerivedTheoryFiles = random.sample( yukawaDerivedTheoryFiles, 8 )


        containers = []
        colorCycle = itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )
        for yukawaDerivedTheoryFile in yukawaDerivedTheoryFiles:
            color = next(colorCycle)

            container = TheoryFileInterface.ReadDerivedTheoryFile( yukawaDerivedTheoryFile )
            container.name = 'kappab_{0}_kappac_{1}'.format(
                Commands.ConvertFloatToStr( container.kappab ),
                Commands.ConvertFloatToStr( container.kappac ),
                )

            container.Tg_theory = TheoryFileInterface.ReadDerivedTheoryFileToTGraph( yukawaDerivedTheoryFile, name = container.name )
            container.Tg_theory.SetLineWidth(2)
            container.Tg_theory.SetLineColor(color)
            container.Tg_theory.SetMarkerColor(color)
            container.Tg_theory.SetLineStyle(1)
            container.Tg_theory.SetMarkerStyle(8)
            container.Tg_theory.SetMarkerSize(0.8)

            container.Tg_parametrization = wsParametrization.GetOutputContainer(
                kappab = container.kappab , kappac = container.kappac,
                returnWhat='theory' ).Tg
            container.Tg_parametrization.SetLineColor(color)
            container.Tg_parametrization.SetMarkerColor(color)
            container.Tg_parametrization.SetLineStyle(1)

            container.Tg_parametrization_expBinning = wsParametrization.GetOutputContainer(
                kappab = container.kappab , kappac = container.kappac,
                returnWhat='exp' ).Tg
            container.Tg_parametrization_expBinning.SetLineColorAlpha( color, 0.5 )
            container.Tg_parametrization_expBinning.SetMarkerColor(color)
            container.Tg_parametrization_expBinning.SetLineStyle(1)

            containers.append( container )


        # ======================================
        # Additional lines to check

        extraTestCouplings = [
            # { 'kappab' : -0.565910458565, 'kappac' : 9.62666416168 },
            ]

        extraTestContainers = []
        for testCouplings in extraTestCouplings:
            color = next(colorCycle)

            container = Container(
                name = '_'.join([ '{0}_{1:.2f}'.format(key, value) for key, value in sorted(testCouplings.iteritems()) ])
                )
            container.binBoundaries = containers[0].binBoundaries

            container.Tg_parametrization = wsParametrization.GetOutputContainer(
                kappab = testCouplings['kappab'] , kappac = testCouplings['kappac'],
                returnWhat='theory' ).Tg
            container.Tg_parametrization.SetLineColor(color)
            container.Tg_parametrization.SetMarkerColor(color)
            container.Tg_parametrization.SetLineStyle(1)

            container.Tg_parametrization_expBinning = wsParametrization.GetOutputContainer(
                kappab = testCouplings['kappab'] , kappac = testCouplings['kappac'],
                returnWhat='exp' ).Tg
            container.Tg_parametrization_expBinning.SetLineColor(color)
            container.Tg_parametrization_expBinning.SetMarkerColor(color)
            container.Tg_parametrization_expBinning.SetLineStyle(1)

            extraTestContainers.append( container )



        # ======================================
        # Make plot

        doBigLegend = False

        c.cd()
        c.Clear()
        if doBigLegend:
            SetCMargins( RightMargin=0.3 )
        else:
            SetCMargins()

        xMinAbs = min([ container.Tg_theory.xMin for container in containers ])
        xMaxAbs = max([ container.Tg_theory.xMax for container in containers ])
        yMinAbs = min([ container.Tg_theory.yMin for container in containers ])
        yMaxAbs = max([ container.Tg_theory.yMax for container in containers ])

        xMin = xMinAbs - 0.1*( xMaxAbs - xMinAbs )
        xMax = xMaxAbs + 0.1*( xMaxAbs - xMinAbs )
        yMin = yMinAbs - 0.1*( yMaxAbs - yMinAbs )
        yMax = yMaxAbs + 0.1*( yMaxAbs - yMinAbs )

        base = GetPlotBase(
            xMin = xMin,
            xMax = xMax,
            yMin = yMin,
            yMax = yMax,
            xTitle = 'p_{T} [GeV]', yTitle = '#mu_{ggH}'
            )
        base.Draw('P')

        if doBigLegend:
            leg = ROOT.TLegend(
                1 - 0.3,
                c.GetBottomMargin(),
                1 - 0.02 ,
                1 - c.GetTopMargin() 
                )
        else:
            leg = ROOT.TLegend(
                1 - c.GetRightMargin() - 0.35,
                1 - c.GetTopMargin() - 0.50,
                1 - c.GetRightMargin(),
                1 - c.GetTopMargin() 
                )
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)

        if not doBigLegend:
            dummyPoint = ROOT.TGraph( 1, array( 'd', [999.] ), array( 'd', [999.] ) )
            dummyPoint.SetMarkerStyle(8)
            dummyPoint.SetLineWidth(2)
            dummyPoint.SetName('legenddummypoint')
            dummyPoint.Draw('SAMEP')
            leg.AddEntry( dummyPoint.GetName(), 'Theory calc.', 'p' )
            leg.AddEntry( dummyPoint.GetName(), 'Parametrization', 'l' )


        if plotCombination:
            # Combination scan result
            combinationPOIs = Commands.ListPOIs( 'workspaces_May15/combinedCard_May15.root' )
            combinationscans = PhysicsCommands.GetScanResults(
                combinationPOIs,
                'Scan_May15',
                pattern = 'combinedCard'
                )
            TgCombination = PhysicsCommands.GetTGraphForSpectrum( combinationPOIs, combinationscans, name='Combination' )
            TgCombination.SetLineColor(1)
            TgCombination.SetFillColorAlpha( 1, 0.2 )

            CorrelationMatrices.ConvertTGraphToLinesAndBoxes(
                TgCombination,
                drawImmediately=True, legendObject=leg, noBoxes=False, xMaxExternal=xMax )


        for container in containers:
            container.Tg_theory.Draw('XP')
            container.Tg_parametrization.Draw('XL')

            if DrawExperimentalBinLines:
                CorrelationMatrices.ConvertTGraphToLinesAndBoxes(
                    container.Tg_parametrization_expBinning,
                    drawImmediately=True,
                    legendObject= leg if doBigLegend else None,
                    noBoxes=True,
                    xMaxExternal=xMax )

            if doBigLegend:
                leg.AddEntry( container.Tg_theory.GetName(), container.name, 'p' )
            else:
                leg.AddEntry(
                    container.Tg_theory.GetName(),
                    '#kappa_{{c}} = {0:d}, #kappa_{{b}} = {1:d}'.format( int(container.kappac), int(container.kappab) ),
                    'lp'
                    )



        for container in extraTestContainers:
            container.Tg_parametrization.Draw('XL')
            if DrawExperimentalBinLines:
                container.Tg_parametrization_expBinning.SetLineStyle(2)
                CorrelationMatrices.ConvertTGraphToLinesAndBoxes(
                    container.Tg_parametrization_expBinning,
                    drawImmediately=True, legendObject=leg, noBoxes=True, xMaxExternal=xMax )

        leg.Draw()

        outname = '{0}_parametrizationCheck'.format( basename(wsToCheck) )
        SaveC( outname )


    #____________________________________________________________________
    if args.combinationAndContour_Yukawa:

        container = Container()

        container.expBinBoundaries    = expBinBoundaries
        container.ws_combination      = LatestPaths.ws_combined_smH
        container.scanDir_combination = LatestPaths.scan_combined_PTH
        container.ws_coupling         = LatestPaths.ws_combined_Yukawa
        container.scanDir_coupling    = LatestPaths.scan_combined_Yukawa

        container.xCoupling           = 'kappac'
        container.yCoupling           = 'kappab'
        container.xCouplingTitle      = '#kappa_{c}'
        container.yCouplingTitle      = '#kappa_{b}'

        container.plotTitle           = 'comparison_Yukawa'

        container.MaximaOfContour     = True
        container.BestFit             = True

        PlotCommands.PlotParametrizationsOnCombination( container )




########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'