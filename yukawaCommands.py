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
from Container import Container
from Parametrization import Parametrization, WSParametrization

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
    parser.add_argument( '--coupling2Dplot_Yukawa',                               action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa',                          action=CustomAction )
    parser.add_argument( '--couplingContourPlotAsimov_Yukawa',                    action=CustomAction )
    parser.add_argument( '--checkWSParametrization_Yukawa',                       action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa_TheoryUncertainties',      action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa_BRdependencyComparison',   action=CustomAction )


########################################
# Methods
########################################    

def main( args ):

    # expBinBoundaries = [ 0., 15., 30., 45., 85., 125. ]
    # # expBinBoundaries = [ 0., 15., 30., 45., 85., 125., 200., 350. ]
    # print 'Hardcoded binBoundaries for Yukawa:'
    # print expBinBoundaries
    # print ''

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

    if args.mergeGluonInducedWithQuarkInduced:
        # TheoryFileInterface.MergeGluonAndQuarkInduced(
        #     gI_theoryFileDir = LatestPaths.derivedTheoryFilesDirectory_YukawaGluonInduced,
        #     qI_theoryFileDir = LatestPaths.derivedTheoryFilesDirectory_YukawaQuarkInduced,
        #     verbose = True
        #     )

        TheoryFileInterface.SumGluonAndQuarkFiles(
            gI_theoryFileDir = LatestPaths.derivedTheoryFilesDirectory_YukawaGluonInduced,
            qI_theoryFileDir = LatestPaths.derivedTheoryFilesDirectory_YukawaQuarkInducedScaled,
            verbose = True
            )


    if args.createDerivedTheoryFiles_YukawaQuarkInducedScaled:
        TheoryFileInterface.ScaleQuarkInduced(
            LatestPaths.derivedTheoryFilesDirectory_YukawaQuarkInduced,
            verbose = True,
            )


    # ======================================
    # Dealing with theory uncertainties

    if args.CorrelationMatrices_Yukawa:

        variationFiles = TheoryFileInterface.FileFinder(
            directory = LatestPaths.derivedTheoryFilesDirectory_YukawaGluonInduced,
            # directory = LatestPaths.derivedTheoryFilesDirectory_YukawaSummed,
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
            TheoryCommands.RebinDerivedTheoryContainer( variation, [ 0., 15., 30., 45., 85., 125. ] )

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

        INCLUDE_THEORY_UNCERTAINTIES = True
        # INCLUDE_THEORY_UNCERTAINTIES = False

        # UNCORRELATED_THEORY_UNCERTAINTIES = True
        UNCORRELATED_THEORY_UNCERTAINTIES = False

        # MAKELUMISCALABLE = True
        MAKELUMISCALABLE = False

        # INCLUDE_BR_COUPLING_DEPENDENCY = True
        INCLUDE_BR_COUPLING_DEPENDENCY = False

        # PROFILE_TOTAL_XS                  = True
        PROFILE_TOTAL_XS                  = False

        # DO_BR_UNCERTAINTIES               = True
        DO_BR_UNCERTAINTIES               = False


        # ======================================
        # 

        datacard = LatestPaths.card_combined_split
        if args.hgg:
            datacard = LatestPaths.card_onlyhgg_split_both_renamed
        if args.hzz:
            datacard = LatestPaths.card_onlyhzz_split_renamed


        TheoryFileInterface.SetFileFinderDir( LatestPaths.derivedTheoryFilesDirectory_YukawaSummed )
        # TheoryFileInterface.SetFileFinderDir( LatestPaths.derivedTheoryFilesDirectory_YukawaGluonInduced )

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
            '--PO binBoundaries=0,15,30,45,85,125'
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


        suffix = 'Yukawa'
        if UNCORRELATED_THEORY_UNCERTAINTIES:
            extraOptions.append(
                '--PO correlationMatrix={0}'.format( LatestPaths.correlationMatrix_Yukawa_uncorrelated ) )
            extraOptions.append(
                '--PO theoryUncertainties={0}'.format( LatestPaths.theoryUncertainties_Yukawa ) )
            suffix += '_withUncorrelatedTheoryUncertainties'
        elif INCLUDE_THEORY_UNCERTAINTIES:
            extraOptions.append(
                '--PO correlationMatrix={0}'.format( LatestPaths.correlationMatrix_Yukawa ) )
            extraOptions.append(
                '--PO theoryUncertainties={0}'.format( LatestPaths.theoryUncertainties_Yukawa ) )
            suffix += '_withTheoryUncertainties'

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


        import random
        random.seed(1002)
        extraOptions.extend( random.sample( possibleTheories, 6 ) )

        Commands.BasicT2WSwithModel(
            datacard,
            'CouplingModel.py',
            suffix = suffix,
            extraOptions = extraOptions,
            )


    if args.couplingBestfit_Yukawa:

        # doFastscan = True
        # if args.notFastscan: doFastscan = False

        doFastscan = False
        if args.notFastscan: doFastscan = False
        if args.fastscan:    doFastscan = True


        # ======================================
        # Flags

        # ASIMOV = True
        ASIMOV = False
        
        # LUMISTUDY = True
        LUMISTUDY = False

        SPLIT = True
        UNCORRELATED_THEORY_UNCERTAINTIES = False
        NO_THEORY_UNCERTAINTIES           = False

        # INCLUDE_BR_COUPLING_DEPENDENCY    = True
        INCLUDE_BR_COUPLING_DEPENDENCY    = False

        # PROFILE_TOTAL_XS                  = True
        PROFILE_TOTAL_XS                  = False

        # FIX_KAPPAV                        = True
        FIX_KAPPAV                        = False

        # DO_BR_UNCERTAINTIES               = True
        DO_BR_UNCERTAINTIES               = False


        DO_ONLY_ONE_KAPPA                 = True
        # DO_ONLY_ONE_KAPPA                 = False

        # theKappa = 'kappab'
        theKappa = 'kappac'

        theOtherKappa = { 'kappab' : 'kappac', 'kappac' : 'kappab' }[theKappa]


        # ======================================
        # 

        if SPLIT:
            # combinedDatacard = LatestPaths.ws_combined_split_yukawa
            # combinedDatacard = LatestPaths.ws_combined_split_betterYukawa
            # hggDatacard      = LatestPaths.ws_onlyhgg_split_yukawa
            # hzzDatacard      = LatestPaths.ws_onlyhzz_split_yukawa

            combinedDatacard   = LatestPaths.ws_combined_yukawa
            hggDatacard        = LatestPaths.ws_hgg_yukawa
            hzzDatacard        = LatestPaths.ws_hzz_yukawa

            if LUMISTUDY:
                combinedDatacard = LatestPaths.ws_combined_split_lumiScalableWS
            if UNCORRELATED_THEORY_UNCERTAINTIES:
                combinedDatacard = LatestPaths.ws_combined_split_top_uncorrelated
            if NO_THEORY_UNCERTAINTIES:
                combinedDatacard = LatestPaths.ws_combined_split_top_notheoryunc

            if INCLUDE_BR_COUPLING_DEPENDENCY:
                # combinedDatacard = LatestPaths.ws_combined_split_betterYukawa_couplingDependentBR
                combinedDatacard = LatestPaths.ws_combined_yukawa_couplingDependentBR

                if DO_BR_UNCERTAINTIES:
                    combinedDatacard = LatestPaths.ws_combined_yukawa_couplingDependentBR_withBRunc

        else:
            combinedDatacard = LatestPaths.ws_combined_unsplit_lumiScalableWS
            hggDatacard      = LatestPaths.ws_onlyhgg_unsplit_yukawa
            hzzDatacard      = LatestPaths.ws_onlyhzz_unsplit_yukawa


        datacard = combinedDatacard
        if args.hgg:
            datacard = hggDatacard
        if args.hzz:
            datacard = hzzDatacard

        kappab_ranges = [ -20., 20. ]
        kappac_ranges = [ -50., 50. ]

        jobDirectory = 'Scan_yukawa_{0}'.format( datestr )
        if ASIMOV: jobDirectory += '_asimov'
        jobDirectory = Commands.AppendNumberToDirNameUntilItDoesNotExistAnymore( jobDirectory )


        if doFastscan:
            nPoints = 6400
            nPointsPerJob = 800
            queue = 'short.q'
        else:
            # nPoints = 6400
            # nPointsPerJob = 8
            nPoints = 4900
            nPointsPerJob = 14
            queue = 'all.q'
            if UNCORRELATED_THEORY_UNCERTAINTIES or NO_THEORY_UNCERTAINTIES or LUMISTUDY:
                nPoints = 3200
                nPointsPerJob = 8
            if args.hzz:
                nPointsPerJob = 320
                queue = 'short.q'

        if DO_ONLY_ONE_KAPPA:
            nPoints       = 39
            nPointsPerJob = 3
            queue         = 'short.q'


        extraOptions = [
            '-P kappab -P kappac' if not DO_ONLY_ONE_KAPPA else '-P {0}'.format(theKappa),
            '--squareDistPoiStep',
            # 
            ( '--setPhysicsModelParameters kappab=1.0,kappac=1.0'
                + ( ',lumiScale=8.356546' if LUMISTUDY else '' )
                + ( ',kappa_V=1.0' if INCLUDE_BR_COUPLING_DEPENDENCY else '' )
                ),
            # 
            ( '--setPhysicsModelParameterRanges kappab={0},{1}:kappac={2},{3}'.format( kappab_ranges[0], kappab_ranges[1], kappac_ranges[0], kappac_ranges[1] )
                + ( ':kappa_V=-100.0,100.0' if ( INCLUDE_BR_COUPLING_DEPENDENCY and not FIX_KAPPAV ) else '' )
                ),
            ]

        if INCLUDE_BR_COUPLING_DEPENDENCY:
            if not FIX_KAPPAV:
                extraOptions.append( '--floatNuisances kappa_V' )
            else:
                extraOptions.append( '--freezeNuisances kappa_V' )

        if DO_ONLY_ONE_KAPPA:
            extraOptions.append( '--floatNuisances {0}'.format(theOtherKappa) )

        variablesToSave = []
        variablesToSave.extend( Commands.ListSet( datacard, 'yieldParameters' ) )
        variablesToSave.extend( [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ] )
        if INCLUDE_BR_COUPLING_DEPENDENCY:
            variablesToSave.extend( Commands.ListSet( datacard, 'hgg_yieldParameters' ) )
            variablesToSave.extend( Commands.ListSet( datacard, 'hzz_yieldParameters' ) )
            variablesToSave.extend( Commands.ListSet( datacard, 'BRvariables' ) )
        if PROFILE_TOTAL_XS:
            variablesToSave.extend([ 'r_totalXS', 'totalXS', 'totalXS_SM' ])

        extraOptions.append( '--saveSpecifiedFunc ' + ','.join(variablesToSave) )


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
    #        2 D  P L O T S
    ################################################################################
    ################################################################################
    # 
    # 

    if args.coupling2Dplot_Yukawa or args.couplingContourPlot_Yukawa:

        ASIMOV = False
        SPLIT  = True

        if SPLIT:

            combinedDatacard = LatestPaths.ws_combined_split_yukawa
            hggDatacard      = LatestPaths.ws_onlyhgg_split_yukawa
            hzzDatacard      = LatestPaths.ws_onlyhzz_split_yukawa

            if not ASIMOV:
                combinedScanFiles = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined_split_profiled ) )
                hggScanFiles      = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_hgg_split_profiled ) )
                hzzScanFiles      = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_hzz_split_profiled ) )
            else:
                combinedScanFiles = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined_split_profiled_asimov ) )
                hggScanFiles      = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_hgg_split_profiled_asimov ) )
                hzzScanFiles      = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_hzz_split_profiled_asimov ) )
        else:

            combinedDatacard = LatestPaths.ws_combined_unsplit_yukawa
            hggDatacard      = LatestPaths.ws_onlyhgg_unsplit_yukawa
            hzzDatacard      = LatestPaths.ws_onlyhzz_unsplit_yukawa

            if not ASIMOV:
                combinedScanFiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_profiled ) )
                hggScanFiles      = glob( '{0}/*.root'.format( LatestPaths.scan_hgg_profiled ) )
                hzzScanFiles      = glob( '{0}/*.root'.format( LatestPaths.scan_hzz_profiled ) )
            else:
                combinedScanFiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_profiled_asimov ) )
                hggScanFiles      = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_hgg_profiled_asimov ) )
                hzzScanFiles      = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_hzz_profiled_asimov ) )


    if args.coupling2Dplot_Yukawa:

        datacard  = combinedDatacard
        rootfiles = combinedScanFiles
        if args.hgg:
            datacard  = hggDatacard
            rootfiles = hggScanFiles
        if args.hzz:
            datacard  = hzzDatacard
            rootfiles = hzzScanFiles

        res = TheoryCommands.PlotCouplingScan2D(
            # datacard,
            rootfiles,
            xCoupling = 'kappac',
            yCoupling = 'kappab',
            xMin = -35., xMax = 35.,
            yMin = -13., yMax = 13.,
            verbose = True
            # multiplyBinContents = 10.,
            )

        print '\nBest fit:'
        print res.xCoupling, '=', res.xBestfit
        print res.yCoupling, '=', res.yBestfit


    #____________________________________________________________________
    if args.couplingContourPlot_Yukawa_BRdependencyComparison:

        xCoupling = 'kappac'
        yCoupling = 'kappab'
        titles = { 'kappac': '#kappa_{c}', 'kappab' : '#kappa_{b}' }

        combined_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined_asimov ) )
        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name = 'regularBR'
        combined.title = 'BR constant'

        scalingBR_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined_asimov_couplingDependentBR ) )
        scalingBR = TheoryCommands.GetTH2FromListOfRootFiles(
            scalingBR_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        scalingBR.color = 2
        scalingBR.name = 'scalingBR'
        scalingBR.title = 'BR(#kappa_{t}, #kappa_{V})'

        scalingBRfixedKappaV_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined_asimov_couplingDependentBR_fixedKappaV ) )
        scalingBRfixedKappaV = TheoryCommands.GetTH2FromListOfRootFiles(
            scalingBRfixedKappaV_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        scalingBRfixedKappaV.color = 4
        scalingBRfixedKappaV.name = 'scalingBRfixedKappaV'
        scalingBRfixedKappaV.title = 'BR(#kappa_{t}) (#kappa_{V} fixed)'


        scalingBR_withBRunc_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined_asimov_couplingDependentBR_withBRunc ) )
        scalingBR_withBRunc = TheoryCommands.GetTH2FromListOfRootFiles(
            scalingBR_withBRunc_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        scalingBR_withBRunc.color = 8
        scalingBR_withBRunc.name = 'scalingBR_withBRunc'
        scalingBR_withBRunc.title = 'BR(#kappa_{t}) (#kappa_{V} fixed)'


        TheoryCommands.BasicMixedContourPlot(
            [
                combined,
                scalingBR,
                scalingBRfixedKappaV,
                scalingBR_withBRunc,
                ],
            xMin = -35.,
            xMax = 35.,
            yMin = -13.,
            yMax = 13.,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_BRcouplingDependency',
            x_SM      = 1.,
            y_SM      = 1.,
            plotIndividualH2s = True,
            )


    if args.couplingContourPlot_Yukawa:

        xCoupling = 'kappac'
        yCoupling = 'kappab'
        titles = { 'kappac': '#kappa_{c}', 'kappab' : '#kappa_{b}' }

        TH2FsToPlot = []

        # combined_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined_split_profiled_asimov ) )
        # hgg_rootfiles      = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_hgg_split_profiled_asimov ) )
        # hzz_rootfiles      = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_hzz_split_profiled_asimov ) )
        combined_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined_split_profiled ) )
        hgg_rootfiles      = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_hgg_split_profiled ) )
        hzz_rootfiles      = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_hzz_split_profiled ) )


        # combined_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined ) )
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

        # hgg_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_hgg ) )
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

        # hzz_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_hzz ) )
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


        TheoryCommands.BasicMixedContourPlot(
            TH2FsToPlot,
            xMin = -35.,
            xMax = 35.,
            yMin = -13.,
            yMax = 13.,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_perDecayChannel',
            x_SM      = 1.,
            y_SM      = 1.,
            plotIndividualH2s = True,
            )




    if False and args.couplingContourPlot_Yukawa:

        xCoupling = 'kappac'
        yCoupling = 'kappab'
        titles = {
            'kappac': '#kappa_{c}', 'kappab' : '#kappa_{b}',
            'hgg' : 'H #rightarrow #gamma#gamma',
            'hzz' : 'H #rightarrow 4l',
            'combined' : 'Combination',
            }

        hgg = TheoryCommands.GetTH2FromListOfRootFiles(
            hggScanFiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        hgg.color = 2
        hgg.name = 'hgg'

        hzz = TheoryCommands.GetTH2FromListOfRootFiles(
            hzzScanFiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        hzz.color = 4
        hzz.name = 'hzz'

        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combinedScanFiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name = 'combined'

        containers = [
            hzz,
            hgg,
            combined,
            ]

        for container in containers:
            container.contours_1sigma = TheoryCommands.GetContoursFromTH2( container.H2, 2.30 )
            container.contours_2sigma = TheoryCommands.GetContoursFromTH2( container.H2, 6.18 )


        c.cd()
        c.Clear()
        SetCMargins()

        xMin = -35.
        xMax = 35.
        yMin = -13.
        yMax = 13.

        base = GetPlotBase(
            xMin = xMin,
            xMax = xMax,
            yMin = yMin,
            yMax = yMax,
            xTitle = titles.get( xCoupling, xCoupling ),
            yTitle = titles.get( yCoupling, yCoupling ),
            )
        base.Draw('P')

        base.GetXaxis().SetTitleSize(0.06)
        base.GetXaxis().SetLabelSize(0.05)
        base.GetYaxis().SetTitleSize(0.06)
        base.GetYaxis().SetLabelSize(0.05)


        # leg = ROOT.TLegend(
        #     1 - c.GetRightMargin() - 0.22,
        #     c.GetBottomMargin() + 0.02,
        #     1 - c.GetRightMargin(),
        #     c.GetBottomMargin() + 0.21
        #     )

        leg = ROOT.TLegend(
            c.GetLeftMargin() + 0.01,
            c.GetBottomMargin() + 0.02,
            1 - c.GetRightMargin() - 0.01,
            c.GetBottomMargin() + 0.09
            )
        leg.SetNColumns(3)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)

        for container in containers:

            for Tg in container.contours_1sigma:
                TheoryCommands.SetExtremaOfContour(Tg)
                if abs( Tg.xMax - Tg.xMin ) < 0.1*(xMax-xMin) or abs( Tg.yMax - Tg.yMin ) < 0.1*(yMax-yMin):
                    continue
                Tg.SetLineWidth(2)
                Tg.SetLineColor( container.color )
                Tg.SetLineStyle(1)
                Tg.Draw('LSAME')
                Tg.SetName( '{0}_contour_1sigma'.format(container.name) )
                leg.AddEntry( Tg.GetName(), titles.get( container.name, container.name ), 'l' )

            for Tg in container.contours_2sigma:
                TheoryCommands.SetExtremaOfContour(Tg)
                if abs( Tg.xMax - Tg.xMin ) < 0.1*(xMax-xMin) or abs( Tg.yMax - Tg.yMin ) < 0.1*(yMax-yMin):
                    continue
                Tg.SetLineWidth(2)
                Tg.SetLineColor( container.color )
                Tg.SetLineStyle(2)
                Tg.Draw('LSAME')

            Tpoint = ROOT.TGraph( 1, array( 'd', [container.xBestfit] ), array( 'd', [container.yBestfit] ) )
            ROOT.SetOwnership( Tpoint, False )
            Tpoint.SetMarkerSize(2)
            Tpoint.SetMarkerStyle(34)
            Tpoint.SetMarkerColor( container.color )
            Tpoint.Draw('PSAME')
            Tpoint.SetName( '{0}_bestfitpoint'.format( container.name ) )

        TpointSM = ROOT.TGraph( 1, array( 'd', [1.0] ), array( 'd', [1.0] ) )
        ROOT.SetOwnership( TpointSM, False )
        TpointSM.SetMarkerSize(2)
        TpointSM.SetMarkerStyle(21)
        TpointSM.SetMarkerColor( 12 )
        TpointSM.Draw('PSAME')

        leg.Draw()

        SaveC( 'contours' )



    if args.couplingContourPlotAsimov_Yukawa:

        xCoupling = 'kappac'
        yCoupling = 'kappab'
        titles = {
            'kappac'       : '#kappa_{c}', 'kappab' : '#kappa_{b}',
            'hgg'          : 'H #rightarrow #gamma#gamma',
            'hzz'          : 'H #rightarrow 4l',
            'combined'     : 'Combination',
            'asimov_lum1'  : 'Expected at 35.9 fb^{-1}',
            'asimov_lum10' : 'Expected at 359.0 fb^{-1}',
            'asimov_lum8'  : 'Expected at 300.0 fb^{-1}',
            }

        # asimov_lum10 = TheoryCommands.GetTH2FromListOfRootFiles(
        #     glob( '{0}/*.root'.format(LatestPaths.scan_combined_profiled_asimov_lum10) ),
        #     xCoupling,
        #     yCoupling,
        #     verbose   = False,
        #     xMin = -20., xMax = 20.,
        #     yMin = -8., yMax = 8.,
        #     )
        # asimov_lum10.color = 2
        # asimov_lum10.name = 'asimov_lum10'


        asimov_lum8 = TheoryCommands.GetTH2FromListOfRootFiles(
            glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined_split_profiled_asimov_lum8 ) ),
            xCoupling,
            yCoupling,
            verbose   = False,
            xMin = -20., xMax = 20.,
            yMin = -8., yMax = 8.,
            )
        asimov_lum8.color = 2
        asimov_lum8.name = 'asimov_lum8'


        asimov_lum1 = TheoryCommands.GetTH2FromListOfRootFiles(
            glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined_split_profiled_asimov ) ),
            xCoupling,
            yCoupling,
            verbose   = False,
            xMin = -30., xMax = 30.,
            yMin = -15., yMax = 15.,
            )
        asimov_lum1.color = 1
        asimov_lum1.name = 'asimov_lum1'


        # asimov_lum10_from1 = deepcopy( asimov_lum1 )
        # asimov_lum10_from1.color = 4
        # asimov_lum10_from1.name = 'Multiplied'
        # asimov_lum10_from1.H2.SetName( TheoryCommands.GetUniqueRootName() )
        # for iX in xrange( asimov_lum10_from1.H2.GetNbinsX() ):
        #     for iY in xrange( asimov_lum10_from1.H2.GetNbinsY() ):
        #         asimov_lum10_from1.H2.SetBinContent(
        #             iX+1, iY+1,
        #             asimov_lum10_from1.H2.GetBinContent( iX+1, iY+1 ) * 10.
        #             )


        containers = [
            asimov_lum1,
            asimov_lum8,
            # asimov_lum10,
            # asimov_lum10_from1,
            ]


        for container in containers:
            container.contours_1sigma = TheoryCommands.GetContoursFromTH2( container.H2, 2.30 )
            container.contours_2sigma = TheoryCommands.GetContoursFromTH2( container.H2, 6.18 )

        c.cd()
        c.Clear()
        SetCMargins()

        xMin = -27.
        xMax = 27.
        yMin = -10.
        yMax = 10.

        base = GetPlotBase(
            xMin = xMin,
            xMax = xMax,
            yMin = yMin,
            yMax = yMax,
            xTitle = titles.get( xCoupling, xCoupling ),
            yTitle = titles.get( yCoupling, yCoupling ),
            )
        base.Draw('P')

        base.GetXaxis().SetTitleSize(0.06)
        base.GetXaxis().SetLabelSize(0.05)
        base.GetYaxis().SetTitleSize(0.06)
        base.GetYaxis().SetLabelSize(0.05)


        # leg = ROOT.TLegend(
        #     1 - c.GetRightMargin() - 0.22,
        #     c.GetBottomMargin() + 0.02,
        #     1 - c.GetRightMargin(),
        #     c.GetBottomMargin() + 0.21
        #     )

        leg = ROOT.TLegend(
            c.GetLeftMargin() + 0.01,
            c.GetBottomMargin() + 0.02,
            1 - c.GetRightMargin() - 0.01,
            c.GetBottomMargin() + 0.09
            )
        leg.SetNColumns( len(containers) )
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)

        for container in containers:

            print 'Drawing contours'
            print container

            for Tg in container.contours_1sigma:
                Tg.SetLineWidth(2)
                Tg.SetLineColor( container.color )
                Tg.SetLineStyle(1)
                Tg.Draw('LSAME')
                if Tg == container.contours_1sigma[0]:
                    Tg.SetName( '{0}_contour_1sigma'.format(container.name) )
                    leg.AddEntry( Tg.GetName(), titles.get( container.name, container.name ), 'l' )

            for Tg in container.contours_2sigma:
                Tg.SetLineWidth(2)
                Tg.SetLineColor( container.color )
                Tg.SetLineStyle(2)
                Tg.Draw('LSAME')

            # Tpoint = ROOT.TGraph( 1, array( 'd', [container.xBestfit] ), array( 'd', [container.yBestfit] ) )
            # ROOT.SetOwnership( Tpoint, False )
            # Tpoint.SetMarkerSize(2)
            # Tpoint.SetMarkerStyle(34)
            # Tpoint.SetMarkerColor( container.color )
            # Tpoint.Draw('PSAME')
            # Tpoint.SetName( '{0}_bestfitpoint'.format( container.name ) )

        TpointSM = ROOT.TGraph( 1, array( 'd', [1.0] ), array( 'd', [1.0] ) )
        ROOT.SetOwnership( TpointSM, False )
        TpointSM.SetMarkerSize(2)
        TpointSM.SetMarkerStyle(21)
        TpointSM.SetMarkerColor( 1 )
        TpointSM.Draw('PSAME')

        leg.Draw()

        SaveC( 'contoursAsimov' )



    if args.checkWSParametrization_Yukawa:

        plotCombination = False

        # DrawExperimentalBinLines = False
        DrawExperimentalBinLines = True


        # # wsToCheck =  LatestPaths.ws_combined_unsplit_yukawa
        # wsToCheck =  LatestPaths.ws_combined_unsplit_yukawa_onlyGluonInduced

        # # theoryDir = LatestPaths.derivedTheoryFilesDirectory_YukawaSummed
        # theoryDir = LatestPaths.derivedTheoryFilesDirectory_YukawaGluonInduced

        # wsToCheck = LatestPaths.ws_onlyhzz_split_yukawa
        wsToCheck = LatestPaths.ws_combined_yukawa
        theoryDir = LatestPaths.derivedTheoryFilesDirectory_YukawaSummed
        

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





    if args.couplingContourPlot_Yukawa_TheoryUncertainties:

        xCoupling = 'kappac'
        yCoupling = 'kappab'
        titles = {
            'kappac': '#kappa_{c}', 'kappab' : '#kappa_{b}',
            'hgg' : 'H #rightarrow #gamma#gamma',
            'hzz' : 'H #rightarrow 4l',
            'combined' : 'Combination',
            }

        # combinedDatacard               = LatestPaths.ws_combined_unsplit_yukawa
        # combinedDatacard_uncorrelated  = LatestPaths.ws_combined_split_top_uncorrelated
        # combinedDatacard_notheoryunc   = LatestPaths.ws_combined_split_top_notheoryunc

        ASIMOV = True
        # ASIMOV = False
       
        if ASIMOV:
            combinedScanFiles = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined_split_profiled_asimov ) )
            combinedScanFiles_uncorrelated = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined_split_profiled_uncorrelated_asimov ) )
            combinedScanFiles_notheoryunc  = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined_split_profiled_notheoryunc_asimov ) )
        else:
            combinedScanFiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_profiled ) )
            combinedScanFiles_uncorrelated = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined_split_profiled_uncorrelated ) )
            combinedScanFiles_notheoryunc  = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined_split_profiled_notheoryunc ) )


        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combinedScanFiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name = 'Regular'

        combined_uncorrelated = TheoryCommands.GetTH2FromListOfRootFiles(
            combinedScanFiles_uncorrelated,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined_uncorrelated.color = 2
        combined_uncorrelated.name = 'Uncorrelated'

        combined_notheoryunc = TheoryCommands.GetTH2FromListOfRootFiles(
            combinedScanFiles_notheoryunc,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined_notheoryunc.color = 4
        combined_notheoryunc.name = 'No_unc'


        containers = [
            combined,
            combined_uncorrelated,
            combined_notheoryunc
            ]

        for container in containers:
            container.contours_1sigma = TheoryCommands.GetContoursFromTH2( container.H2, 2.30 )
            container.contours_2sigma = TheoryCommands.GetContoursFromTH2( container.H2, 6.18 )


        c.cd()
        c.Clear()
        SetCMargins()

        xMin = -25.
        xMax = 25.
        yMin = -10.
        yMax = 10.

        base = GetPlotBase(
            xMin = xMin,
            xMax = xMax,
            yMin = yMin,
            yMax = yMax,
            xTitle = titles.get( xCoupling, xCoupling ),
            yTitle = titles.get( yCoupling, yCoupling ),
            )
        base.Draw('P')

        base.GetXaxis().SetTitleSize(0.06)
        base.GetXaxis().SetLabelSize(0.05)
        base.GetYaxis().SetTitleSize(0.06)
        base.GetYaxis().SetLabelSize(0.05)


        # leg = ROOT.TLegend(
        #     1 - c.GetRightMargin() - 0.22,
        #     c.GetBottomMargin() + 0.02,
        #     1 - c.GetRightMargin(),
        #     c.GetBottomMargin() + 0.21
        #     )

        leg = ROOT.TLegend(
            c.GetLeftMargin() + 0.01,
            c.GetBottomMargin() + 0.02,
            1 - c.GetRightMargin() - 0.01,
            c.GetBottomMargin() + 0.09
            )
        leg.SetNColumns(3)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)

        for container in containers:

            addToLegend = True
            for Tg in container.contours_1sigma:
                TheoryCommands.SetExtremaOfContour(Tg)
                if abs( Tg.xMax - Tg.xMin ) < 0.1*(xMax-xMin) or abs( Tg.yMax - Tg.yMin ) < 0.1*(yMax-yMin):
                    continue
                Tg.SetLineWidth(2)
                Tg.SetLineColor( container.color )
                Tg.SetLineStyle(1)
                Tg.Draw('LSAME')
                Tg.SetName( '{0}_contour_1sigma'.format(container.name) )
                if addToLegend:
                    leg.AddEntry( Tg.GetName(), titles.get( container.name, container.name ), 'l' )
                    addToLegend = False

            for Tg in container.contours_2sigma:
                TheoryCommands.SetExtremaOfContour(Tg)
                if abs( Tg.xMax - Tg.xMin ) < 0.1*(xMax-xMin) or abs( Tg.yMax - Tg.yMin ) < 0.1*(yMax-yMin):
                    continue
                Tg.SetLineWidth(2)
                Tg.SetLineColor( container.color )
                Tg.SetLineStyle(2)
                Tg.Draw('LSAME')

            Tpoint = ROOT.TGraph( 1, array( 'd', [container.xBestfit] ), array( 'd', [container.yBestfit] ) )
            ROOT.SetOwnership( Tpoint, False )
            Tpoint.SetMarkerSize(2)
            Tpoint.SetMarkerStyle(34)
            Tpoint.SetMarkerColor( container.color )
            Tpoint.Draw('PSAME')
            Tpoint.SetName( '{0}_bestfitpoint'.format( container.name ) )

        TpointSM = ROOT.TGraph( 1, array( 'd', [1.0] ), array( 'd', [1.0] ) )
        ROOT.SetOwnership( TpointSM, False )
        TpointSM.SetMarkerSize(2)
        TpointSM.SetMarkerStyle(21)
        TpointSM.SetMarkerColor( 12 )
        TpointSM.Draw('PSAME')

        leg.Draw()

        SaveC( 'contours_theoryUncs' )






########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'