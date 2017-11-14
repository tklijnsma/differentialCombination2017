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

    # Plotting
    parser.add_argument( '--OneKappaScanPlot_Yukawa',                             action=CustomAction )
    parser.add_argument( '--coupling2Dplot_Yukawa',                               action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa',                          action=CustomAction )
    parser.add_argument( '--couplingContourPlotAsimov_Yukawa',                    action=CustomAction )
    parser.add_argument( '--checkWSParametrization_Yukawa',                       action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa_TheoryCrossCheck',         action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa_TheoryUncertainties',      action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa_BRdependencyComparison',   action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa_highLumiStudy',            action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa_ProfiledTotalXS',          action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Yukawa_atfchi2',                  action=CustomAction )


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

        INCLUDE_BR_COUPLING_DEPENDENCY    = True
        # INCLUDE_BR_COUPLING_DEPENDENCY    = False

        # PROFILE_TOTAL_XS                  = True
        PROFILE_TOTAL_XS                  = False

        # DO_BR_UNCERTAINTIES               = True
        DO_BR_UNCERTAINTIES               = False


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

        # LUMISTUDY = True
        LUMISTUDY = False

        # UNCORRELATED_THEORY_UNCERTAINTIES = True
        UNCORRELATED_THEORY_UNCERTAINTIES = False

        # NO_THEORY_UNCERTAINTIES           = True
        NO_THEORY_UNCERTAINTIES           = False

        INCLUDE_BR_COUPLING_DEPENDENCY    = True
        # INCLUDE_BR_COUPLING_DEPENDENCY    = False

        # PROFILE_TOTAL_XS                  = True
        PROFILE_TOTAL_XS                  = False

        # FIX_KAPPAV                        = True
        FIX_KAPPAV                        = False

        # DO_BR_UNCERTAINTIES               = True
        DO_BR_UNCERTAINTIES               = False

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
        print 'INCLUDE_BR_COUPLING_DEPENDENCY    = ', INCLUDE_BR_COUPLING_DEPENDENCY
        print 'PROFILE_TOTAL_XS                  = ', PROFILE_TOTAL_XS
        print 'FIX_KAPPAV                        = ', FIX_KAPPAV
        print 'DO_BR_UNCERTAINTIES               = ', DO_BR_UNCERTAINTIES
        print 'DO_ONLY_ONE_KAPPA                 = ', DO_ONLY_ONE_KAPPA
        if DO_ONLY_ONE_KAPPA:
            print 'Chosen kappa                      = {0}'.format(theKappa )
            print 'Other kappa                       = {0}'.format(theOtherKappa )


        # ======================================
        # Set correct input

        # if INCLUDE_BR_COUPLING_DEPENDENCY:
        #     # combinedDatacard = LatestPaths.ws_combined_split_betterYukawa_couplingDependentBR
        #     combinedDatacard = LatestPaths.ws_combined_yukawa_couplingDependentBR

        #     if DO_BR_UNCERTAINTIES:
        #         combinedDatacard = LatestPaths.ws_combined_yukawa_couplingDependentBR_withBRunc

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

        if DO_ONLY_ONE_KAPPA:
            nPoints       = 39
            nPointsPerJob = 3
            queue         = 'short.q'


        # ======================================
        # Construct the fit command and process flags

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
            jobDirectory += '_oneKappa_{0}'.format(theKappa)
            extraOptions.append( '--floatNuisances {0}'.format(theOtherKappa) )
        if LUMISTUDY:
            jobDirectory += '_lumiStudy'
        if PROFILE_TOTAL_XS:
            jobDirectory += '_profiledTotalXS'
        if UNCORRELATED_THEORY_UNCERTAINTIES:
            jobDirectory += '_uncorrelatedTheoryUncertainties'
        if NO_THEORY_UNCERTAINTIES:
            jobDirectory += '_noTheoryUncertainties'
        if INCLUDE_BR_COUPLING_DEPENDENCY:
            jobDirectory += '_couplingDependentBR'
            if FIX_KAPPAV:
                jobDirectory += '_fixedKappaV'
        if ASIMOV:
            jobDirectory += '_asimov'

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
        extraOptions.append( '--saveSpecifiedFunc ' + ','.join(variablesToSave) )


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

            TheoryCommands.PlotMultipleScans(
                [ container ],
                xTitle   = kappa,
                yTitle   = '#Delta NLL',
                yMax     = 10.,
                plotname = 'oneKappaScan_{0}_{1}'.format( kappa, basename(scanDir).replace('/','') ),
                draw1sigmaline = True
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

        # ASIMOV = True
        ASIMOV = False
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


        TheoryCommands.BasicMixedContourPlot(
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
        combined.title = 'Regular'
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

        TheoryCommands.BasicMixedContourPlot(
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
        combined.title = 'Regular'
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

        TheoryCommands.BasicMixedContourPlot(
            containers,
            xMin = -35.,
            xMax = 35.,
            yMin = -15.,
            yMax = 15.,
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

        scalingBR_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_couplingDependentBR_asimov ) )
        scalingBR = TheoryCommands.GetTH2FromListOfRootFiles(
            scalingBR_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        scalingBR.color = 2
        scalingBR.name = 'scalingBR'
        scalingBR.title = 'BR(#kappa_{t}, #kappa_{V})'
        containers.append( scalingBR )

        scalingBRfixedKappaV_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_couplingDependentBR_fixedKappaV_asimov ) )
        scalingBRfixedKappaV = TheoryCommands.GetTH2FromListOfRootFiles(
            scalingBRfixedKappaV_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        scalingBRfixedKappaV.color = 4
        scalingBRfixedKappaV.name = 'scalingBRfixedKappaV'
        scalingBRfixedKappaV.title = 'BR(#kappa_{t}) (#kappa_{V} fixed)'
        containers.append( scalingBRfixedKappaV )

        # scalingBR_withBRunc_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_yukawa_combined_asimov_couplingDependentBR_withBRunc ) )
        # scalingBR_withBRunc = TheoryCommands.GetTH2FromListOfRootFiles(
        #     scalingBR_withBRunc_rootfiles,
        #     xCoupling,
        #     yCoupling,
        #     verbose   = False,
        #     )
        # scalingBR_withBRunc.color = 8
        # scalingBR_withBRunc.name = 'scalingBR_withBRunc'
        # scalingBR_withBRunc.title = 'BR(#kappa_{t}) (#kappa_{V} fixed)'


        TheoryCommands.BasicMixedContourPlot(
            containers,
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


        TheoryCommands.BasicMixedContourPlot(
            containers,
            xMin = -35.,
            xMax = 35.,
            yMin = -13.,
            yMax = 13.,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_highLumiStudy',
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

        combined_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Yukawa_asimov ) )
        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name  = 'regular'
        combined.title = 'Full likelihood'
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


        TheoryCommands.BasicMixedContourPlot(
            containers,
            xMin = -35.,
            xMax = 35.,
            yMin = -13.,
            yMax = 13.,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_atfchi2',
            x_SM      = 1.,
            y_SM      = 1.,
            plotIndividualH2s = True,
            )

    #____________________________________________________________________
    if args.checkWSParametrization_Yukawa:

        plotCombination = False

        # DrawExperimentalBinLines = False
        DrawExperimentalBinLines = True


        # # wsToCheck =  LatestPaths.ws_combined_unsplit_yukawa
        # wsToCheck =  LatestPaths.ws_combined_unsplit_yukawa_onlyGluonInduced

        # # theoryDir = LatestPaths.derivedTheoryFiles_YukawaSummed
        # theoryDir = LatestPaths.derivedTheoryFiles_YukawaGluonInduced

        # wsToCheck = LatestPaths.ws_onlyhzz_split_yukawa
        wsToCheck = LatestPaths.ws_combined_yukawa
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



########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'