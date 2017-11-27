#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, itertools, operator, re, argparse, sys, random, numpy
from math import isnan, isinf
from os.path import *
from glob import glob
from copy import deepcopy
from array import array

import LatestPaths

sys.path.append('src')
import Commands
import PhysicsCommands
import OneOfCommands
import TheoryCommands
import CorrelationMatrices
import MergeHGGWDatacards
import TheoryFileInterface
import OutputInterface
from Container import Container
from Parametrization import Parametrization, WSParametrization


from time import strftime
datestr = strftime( '%b%d' )



import ROOT
from TheoryCommands import c
from TheoryCommands import SaveC
from TheoryCommands import GetPlotBase
from TheoryCommands import SetCMargins


########################################
# Main
########################################

def AppendParserOptions( parser ):

    parser.add_argument( '--extrastudyCommands', default=False )
    class CustomAction(argparse.Action):
        def __init__(self, option_strings, dest, **kwargs):
            # Only 1 argument allowed (this is basically the "store_true" action)
            super(CustomAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            # Make sure the parser knows a combineCommand was entered
            setattr( namespace, 'extrastudyCommands', True )
            setattr( namespace, self.dest, True )
            # super(FooAction, self).__call__( parser, namespace, values, option_string=None )
            # setattr( namespace, self.dest, values )

    parser.add_argument( '--FitBR_t2ws',                              action=CustomAction )
    parser.add_argument( '--FitBR_bestfit',                           action=CustomAction )
    parser.add_argument( '--FitBR_scan',                              action=CustomAction )
    parser.add_argument( '--FitBR_plot',                              action=CustomAction )

    parser.add_argument( '--TotalXS_t2ws',                            action=CustomAction )
    parser.add_argument( '--TotalXS_bestfit',                         action=CustomAction )
    parser.add_argument( '--TotalXS_scan',                            action=CustomAction )
    parser.add_argument( '--TotalXS_plot',                            action=CustomAction )

    parser.add_argument( '--chi2fitToCombination_Yukawa',             action=CustomAction )
    parser.add_argument( '--RepeatTheoristFit',                       action=CustomAction )

    parser.add_argument( '--Make2DPlotOfVariableInWS',                action=CustomAction )
    parser.add_argument( '--PlotOfTotalXSInYukawaWS',                 action=CustomAction )
    parser.add_argument( '--PlotOfTotalXS_FromParametrization',       action=CustomAction )
    parser.add_argument( '--PlotBRsInOnePlot',                        action=CustomAction )

    parser.add_argument( '--CorrelationMatrixScaleDependence_Yukawa', action=CustomAction )

    parser.add_argument( '--PlotParametrizationShapes',                        action=CustomAction )


########################################
# Methods
########################################    

def main( args ):

    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )


    #____________________________________________________________________
    if args.FitBR_t2ws:

        extraOptions = [
            '--PO verbose=2',
            '--PO \'higgsMassRange=123,127\'',
            ]

        USE_GLOBAL_SCALES = True
        # USE_GLOBAL_SCALES = False

        if USE_GLOBAL_SCALES:
            Commands.BasicT2WSwithModel(
                LatestPaths.card_combined_smH_PTH,
                'FitBRModel.py',
                extraOptions = extraOptions,
                modelName    = 'fitGlobalBRModel',
                suffix       = 'globalScales',
                )

        else:
            Commands.BasicT2WSwithModel(
                LatestPaths.card_combined_smH_PTH,
                'FitBRModel.py',
                extraOptions = extraOptions,
                smartMaps    = [
                    ( r'.*/smH_PTH_([\d\_GT]+)', r'r_smH_PTH_\1[1.0,-1.0,4.0]' )
                    ],
                )

    #____________________________________________________________________
    if args.FitBR_bestfit:

        # USE_GLOBAL_SCALES = True
        USE_GLOBAL_SCALES = False

        if USE_GLOBAL_SCALES:
            ws = abspath( LatestPaths.ws_FitBR_combined_unsplit )
        else:
            ws = abspath( LatestPaths.ws_FitBR_combined_unsplit )

        cmd = [
            'combine',
            ws,
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            # '--algo=grid',
            '--floatOtherPOIs=1',
            '--saveNLL',
            '--saveInactivePOI 1',
            # '--fastScan',
            # '-P kappab',
            # '-P kappac',
            # '--setPhysicsModelParameters kappab=1.0,kappac=1.0',
            # '--saveSpecifiedFunc {0}'.format(','.join(
            #     Commands.ListSet( datacard, 'yieldParameters' ) + [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
            '--squareDistPoiStep',
            '-M MultiDimFit',
            '-m 125.00',
            # '--setPhysicsModelParameterRanges kappab=-20.0,20.0:kappac=-50.0,50.0',
            # '--setPhysicsModelParameterRanges kappab=0.5,1.0:kappac=1.0,2.0',
            # '--points 12800',
            # '--firstPoint 0',
            # '--lastPoint 799',
            '-n testjob',
            # '-v 3',
            ]

        Commands.BasicGenericCombineCommand(
            cmd,
            onBatch = False,
            )



    if args.FitBR_scan:

        doFastscan = False
        if args.fastscan: doFastscan = True

        ASIMOV = False
        if args.asimov: ASIMOV = True
        
        USE_GLOBAL_SCALES = True
        # USE_GLOBAL_SCALES = False

        jobDirectory = 'Scan_ratioOfBRs_{0}'.format( datestr )


        if USE_GLOBAL_SCALES:
            datacard = abspath( LatestPaths.ws_ratio_of_BRs_globalScales )
            jobDirectory += '_globalScales'
        else:
            datacard = abspath( LatestPaths.ws_ratio_of_BRs )

        ratio_BR_hgg_hzz_ranges = [ 0.03, 0.16 ]


        if ASIMOV: jobDirectory += '_asimov'
        jobDirectory = Commands.AppendNumberToDirNameUntilItDoesNotExistAnymore( jobDirectory )

        if doFastscan:
            nPoints = 42
            nPointsPerJob = 42
            queue = 'short.q'
        else:
            nPoints = 42
            nPointsPerJob = 3
            queue = 'short.q'

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
            extraOptions  = [
                # '--importanceSampling={0}:couplingScan'.format( abspath('scanTH2D_Jun01.root') ),
                '-P ratio_BR_hgg_hzz',
                '--setPhysicsModelParameterRanges ratio_BR_hgg_hzz={0},{1}'.format( ratio_BR_hgg_hzz_ranges[0], ratio_BR_hgg_hzz_ranges[1] ),
                '--setPhysicsModelParameters {0}'.format(
                    ','.join([ '{0}=1.0'.format(i) for i in Commands.ListSet( datacard, 'POI' ) if i.startswith('r_') ])
                    ),
                # '--setPhysicsModelParameterRanges kappab={0},{1}:kappac={2},{3}'.format(
                #     kappab_ranges[0], kappab_ranges[1], kappac_ranges[0], kappac_ranges[1] ),
                # '--saveSpecifiedFunc {0}'.format(','.join(
                #     Commands.ListSet( datacard, 'yieldParameters' ) + [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
                '--squareDistPoiStep',
                ]
            )


    if args.FitBR_plot:


        USE_GLOBAL_SCALES = True
        # USE_GLOBAL_SCALES = False

        if USE_GLOBAL_SCALES:
            scanRootFiles = glob( LatestPaths.scan_ratioOfBRs_globalScales + '/*.root' )
        else:
            scanRootFiles = glob( LatestPaths.scan_ratioOfBRs + '/*.root' )
        
        
        scanContainer = OutputInterface.OutputContainer()

        x_unfiltered, y_unfiltered = TheoryCommands.BasicReadScan(
            scanRootFiles,
            xAttr = 'ratio_BR_hgg_hzz',
            yAttr = 'deltaNLL',
            )

        scanContainer.x = []
        scanContainer.y = []
        for x, y in zip( x_unfiltered, y_unfiltered ):
            if (
                x > -10.0 and x < 10.
                and
                y > -10.0 and y < 10.
                ):
                scanContainer.x.append( x )
                scanContainer.y.append( y )


        # Do uncertainty determination before scaling
        # FindMinimaAndErrors expects deltaNLL, not 2*deltaNLL
        scanContainer.extrema = PhysicsCommands.FindMinimaAndErrors( scanContainer.x, scanContainer.y, returnContainer=True )

        print '[info] Multiplying by 2: deltaNLL --> chi^2'
        scanContainer.y = [ 2.*y for y in scanContainer.y ]

        scanContainer.GetTGraph( xAttr = 'x', yAttr = 'y', xAreBinBoundaries = False )
        scanContainer.Tg.SetMarkerStyle(8)
        scanContainer.Tg.SetMarkerSize(0.8)


        # ======================================
        # Make plot

        c.Clear()
        SetCMargins()

        yMinAbs = min( scanContainer.y )
        yMaxAbs = max( scanContainer.y )
        # yMin = yMinAbs - 0.1*(yMaxAbs-yMinAbs)
        yMin = 0.0
        # yMax = yMaxAbs + 0.5*(yMaxAbs-yMinAbs)
        yMax = 5.0

        xMin = min( scanContainer.x )
        xMax = max( scanContainer.x )

        base = GetPlotBase(
            xMin = xMin,
            xMax = xMax,
            yMin = yMin,
            yMax = yMax,
            xTitle = 'BR_{H #rightarrow #gamma#gamma} / BR_{H #rightarrow ZZ}',
            # yTitle = '#Delta NLL',
            yTitle = '2#DeltaNLL',
            )
        base.Draw('P')


        scanContainer.Tg.Draw('XPL')


        oneLine = ROOT.TLine( xMin, 1.0, xMax, 1.0 )
        oneLine.SetLineColor(12)
        oneLine.Draw()


        print '\n' + '-'*70
        print 'Found range: {0:.4f} < ratio < {1:.4f}'.format( scanContainer.extrema.leftBound, scanContainer.extrema.rightBound )

        # imin
        # min
        # leftError
        # leftBound
        # rightError
        # rightBound
        # wellDefinedRightBound
        # wellDefinedLeftBound

        # Check number more carefully
        SM_ratio = LatestPaths.SM_ratio_of_BRs
        SMLine = ROOT.TLine( SM_ratio, yMin, SM_ratio, yMax )
        SMLine.SetLineWidth(2)
        SMLine.SetLineColor(9)
        SMLine.Draw()


        xBestfit = scanContainer.x[ scanContainer.extrema.imin ]
        bestfitLine = ROOT.TLine( xBestfit, yMin, xBestfit, yMax )
        bestfitLine.SetLineWidth(2)
        bestfitLine.SetLineColor(2)
        bestfitLine.Draw()


        l = ROOT.TLatex()
        l.SetTextColor(2)
        l.SetTextSize(0.04)

        l.SetTextAlign(31)
        l.DrawLatex(
            scanContainer.extrema.leftBound, 1.0 + 0.013*(yMax-yMin),
            '-{0:.2f} ({1:d}%)'.format(
                abs(scanContainer.extrema.leftError),
                int( abs(scanContainer.extrema.leftError) / xBestfit * 100. )
                )
            )

        l.SetTextAlign(11)
        l.DrawLatex(
            scanContainer.extrema.rightBound, 1.0 + 0.013*(yMax-yMin),
            '+{0:.2f} ({1:d}%)'.format(
                abs(scanContainer.extrema.rightError),
                int( abs(scanContainer.extrema.rightError) / xBestfit * 100. )
                )
            )

        l.SetTextAlign(21)
        l.DrawLatex(
            xBestfit, 1.0 + 0.013*(yMax-yMin),
            '{0:.3f}'.format( xBestfit )
            )


        TgPoints = ROOT.TGraph( 2,
            array( 'f', [ scanContainer.extrema.leftBound, scanContainer.extrema.rightBound ] ),
            array( 'f', [ 1.0, 1.0 ] ),
            )
        TgPoints.SetMarkerSize(1.2)
        TgPoints.SetMarkerStyle(8)
        TgPoints.SetMarkerColor(2)
        TgPoints.Draw('PSAME')


        SaveC( 'BRscan' + ( '_globalScales' if USE_GLOBAL_SCALES else '' ) )




    #____________________________________________________________________
    if args.TotalXS_t2ws:

        extraOptions = [
            '--PO verbose=2',
            '--PO \'higgsMassRange=123,127\'',
            ]

        Commands.BasicT2WSwithModel(
            LatestPaths.card_combined_smH_PTH,
            pathToModel = 'FitBRModel.py',
            modelName   = 'fitTotalXSModel',
            suffix       = 'fitTotalXS',
            extraOptions = extraOptions
            )

    #____________________________________________________________________
    if args.TotalXS_bestfit:

        ws = LatestPaths.ws_totalXS

        cmd = [
            'combine',
            ws,
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            # '--algo=grid',
            '--floatOtherPOIs=1',
            '--saveNLL',
            '--saveInactivePOI 1',
            '-P r',
            # '--fastScan',
            # '-P kappab',
            # '-P kappac',
            # '--setPhysicsModelParameters kappab=1.0,kappac=1.0',
            # '--saveSpecifiedFunc {0}'.format(','.join(
            #     Commands.ListSet( datacard, 'yieldParameters' ) + [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
            '--squareDistPoiStep',
            '-M MultiDimFit',
            '-m 125.00',
            # '--setPhysicsModelParameterRanges kappab=-20.0,20.0:kappac=-50.0,50.0',
            # '--setPhysicsModelParameterRanges kappab=0.5,1.0:kappac=1.0,2.0',
            # '--points 12800',
            # '--firstPoint 0',
            # '--lastPoint 799',
            '-n testjob',
            # '-v 3',
            ]

        Commands.BasicGenericCombineCommand(
            cmd,
            onBatch = False,
            )


    #____________________________________________________________________
    if args.TotalXS_scan:

        ws = LatestPaths.ws_totalXS

        doFastscan = False
        if args.fastscan: doFastscan = True

        ASIMOV = False
        if args.asimov: ASIMOV = True

        totalXS_ranges = [ 0., 2. ]

        jobDirectory = 'Scan_TotalXS_{0}'.format( datestr )
        if ASIMOV: jobDirectory += '_asimov'
        jobDirectory = Commands.AppendNumberToDirNameUntilItDoesNotExistAnymore( jobDirectory )

        if doFastscan:
            nPoints = 42
            nPointsPerJob = 42
            queue = 'short.q'
            queue = '8nm'
        else:
            nPoints = 50
            nPointsPerJob = 5
            queue = 'short.q'

        Commands.MultiDimCombineTool(
            ws,
            nPoints       = nPoints,
            nPointsPerJob = nPointsPerJob,
            queue         = queue,
            notOnBatch    = False,
            jobDirectory  = jobDirectory,
            fastscan      = doFastscan,
            asimov        = ASIMOV,
            jobPriority   = 0,
            extraOptions  = [
                # '--importanceSampling={0}:couplingScan'.format( abspath('scanTH2D_Jun01.root') ),
                '-P r',
                '--setPhysicsModelParameterRanges r={0},{1}'.format( totalXS_ranges[0], totalXS_ranges[1] ),
                # '--setPhysicsModelParameters {0}'.format(
                #     ','.join([ '{0}=1.0'.format(i) for i in Commands.ListSet( datacard, 'POI' ) if i.startswith('r_') ])
                #     ),
                # '--setPhysicsModelParameterRanges kappab={0},{1}:kappac={2},{3}'.format(
                #     kappab_ranges[0], kappab_ranges[1], kappac_ranges[0], kappac_ranges[1] ),
                # '--saveSpecifiedFunc {0}'.format(','.join(
                #     Commands.ListSet( datacard, 'yieldParameters' ) + [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
                '--squareDistPoiStep',
                ]
            )



    if args.TotalXS_plot:

        scanRootFiles = glob( LatestPaths.scan_combined_totalXS + '/*.root' )

        scanContainer = OutputInterface.OutputContainer()


        x_unfiltered, y_unfiltered = TheoryCommands.BasicReadScan(
            scanRootFiles,
            xAttr = 'r',
            yAttr = 'deltaNLL',
            )

        scanContainer.x = []
        scanContainer.y = []
        for x, y in zip( x_unfiltered, y_unfiltered ):
            if (
                x > -10.0 and x < 10.
                and
                y > -10.0 and y < 10.
                ):
                scanContainer.x.append( x )
                scanContainer.y.append( y )


        # FindMinimaAndErrors assumes deltaNLL (not 2*deltaNLL), so compute it before scaling
        scanContainer.extrema = PhysicsCommands.FindMinimaAndErrors( scanContainer.x, scanContainer.y, returnContainer=True )

        print '[info] Multiplying by 2: deltaNLL --> chi^2'
        scanContainer.y = [ 2.*y for y in scanContainer.y ]

        scanContainer.GetTGraph( xAttr = 'x', yAttr = 'y', xAreBinBoundaries = False )
        scanContainer.Tg.SetMarkerStyle(8)
        scanContainer.Tg.SetMarkerSize(0.8)


        # ======================================
        # Make plot

        c.Clear()
        SetCMargins()

        yMinAbs = min( scanContainer.y )
        yMaxAbs = max( scanContainer.y )
        # yMin = yMinAbs - 0.1*(yMaxAbs-yMinAbs)
        yMin = 0.0
        # yMax = yMaxAbs + 0.5*(yMaxAbs-yMinAbs)
        yMax = 5.0

        # xMin = min( scanContainer.x )
        # xMax = max( scanContainer.x )
        xMin = 0.6
        xMax = 1.4

        base = GetPlotBase(
            xMin = xMin,
            xMax = xMax,
            yMin = yMin,
            yMax = yMax,
            xTitle = '#sigma/#sigma_{SM}',
            # yTitle = '#Delta NLL',
            yTitle = '2#DeltaNLL',
            )
        base.Draw('P')


        scanContainer.Tg.Draw('XPL')


        oneLine = ROOT.TLine( xMin, 1.0, xMax, 1.0 )
        oneLine.SetLineColor(12)
        oneLine.Draw()

        print '\n' + '-'*70
        print 'Found range: {0:.4f} < r_totalXS < {1:.4f}'.format( scanContainer.extrema.leftBound, scanContainer.extrema.rightBound )

        # imin
        # min
        # leftError
        # leftBound
        # rightError
        # rightBound
        # wellDefinedRightBound
        # wellDefinedLeftBound


        xBestfit = scanContainer.x[ scanContainer.extrema.imin ]
        bestfitLine = ROOT.TLine( xBestfit, yMin, xBestfit, yMax )
        bestfitLine.SetLineWidth(2)
        bestfitLine.SetLineColor(2)
        bestfitLine.Draw()

        # for x in [ scanContainer.extrema.leftBound, scanContainer.extrema.rightBound ]:
        #     uncLine = ROOT.TLine( x, 0.0, x, 1.0 )
        #     ROOT.SetOwnership( uncLine, False )
        #     uncLine.SetLineWidth(1)
        #     uncLine.SetLineColor(2)
        #     uncLine.Draw()

        smLine = ROOT.TLine( 1.0, yMin, 1.0, yMax )
        smLine.SetLineWidth(2)
        smLine.SetLineColor(9)
        smLine.Draw()




        l = ROOT.TLatex()
        l.SetTextColor(2)
        l.SetTextSize(0.04)

        l.SetTextAlign(31)
        l.DrawLatex(
            scanContainer.extrema.leftBound, 1.0 + 0.013*(yMax-yMin),
            '-{0:.2f} ({1:d}%)'.format(
                abs(scanContainer.extrema.leftError),
                int( abs(scanContainer.extrema.leftError) / xBestfit * 100. )
                )
            )

        l.SetTextAlign(11)
        l.DrawLatex(
            scanContainer.extrema.rightBound, 1.0 + 0.013*(yMax-yMin),
            '+{0:.2f} ({1:d}%)'.format(
                abs(scanContainer.extrema.rightError),
                int( abs(scanContainer.extrema.rightError) / xBestfit * 100. )
                )
            )

        l.SetTextAlign(21)
        l.DrawLatex(
            xBestfit, 1.0 + 0.013*(yMax-yMin),
            '{0:.3f}'.format( xBestfit )
            )


        TgPoints = ROOT.TGraph( 2,
            array( 'f', [ scanContainer.extrema.leftBound, scanContainer.extrema.rightBound ] ),
            array( 'f', [ 1.0, 1.0 ] ),
            )
        TgPoints.SetMarkerSize(1.2)
        TgPoints.SetMarkerStyle(8)
        TgPoints.SetMarkerColor(2)
        TgPoints.Draw('PSAME')



        SaveC( 'totalXSscan' )


    if args.chi2fitToCombination_Yukawa:
        TheoryCommands.SetPlotDir( 'plots_{0}_Yukawa'.format(datestr) )

        # Workspace to get the parametrization from
        # ws = LatestPaths.ws_combined_split_betterYukawa

        # Derived theory files to get the parametrization from

        # File to get the correlation matrix from
        # corrMatFile = 'corrMat_Oct17/higgsCombine_CORRMAT_combinedCard_Jul26.MultiDimFit.mH125.root'
        # corrMatFile = 'corrMat_Oct17/higgsCombine_CORRMAT_combinedCard_Aug21.MultiDimFit.mH125.root'
        # corrMatFile = 'corrMat_Oct19/higgsCombine_CORRMAT_combinedCard_Aug21_xHfixed.MultiDimFit.mH125.root'
        corrMatFile = 'corrMat_Nov08_combinedCard_Nov03_xHfixed/higgsCombine_CORRMAT_combinedCard_Nov03_xHfixed.MultiDimFit.mH125.root'

        # Scan to get the uncertainties from
        # scanDir = LatestPaths.scan_ptcombination_combined_profiled
        # scanDir = LatestPaths.scan_ptcombination_combined_profiled_xHfixed
        scanDir = LatestPaths.scan_combined_PTH_xHfixed


        expBinning = [ 0., 15., 30., 45., 85., 125. ]
        yieldParameterNames = []
        for left, right in zip( expBinning[:-1], expBinning[1:] ):
            yieldParameterNames.append( 'r_ggH_PTH_{0}_{1}'.format( int(left), int(right) ) )
        nBinsExp = len(expBinning) - 1


        # # ======================================
        # # Load parametrization from WS

        # rootFp = ROOT.TFile.Open(ws)
        # w = rootFp.Get('w')

        # yieldParameterNames = Commands.ListSet( ws, 'yieldParameters' )
        # yieldParameters = [ w.function(yPName) for yPName in yieldParameterNames ]
        # kappacRealVar = w.var('kappac')
        # kappabRealVar = w.var('kappab')

        # binBoundaries   = [ w.var(binBoundName).getVal() for binBoundName in Commands.ListSet( ws, 'expBinBoundaries' ) ]
        # nBins = len(binBoundaries)-1


        # ======================================
        # Load parametrization from derived theory files

        SM = TheoryFileInterface.FileFinder(
            kappab=1, kappac=1, muR=1, muF=1, Q=1,
            expectOneFile=True,
            directory = LatestPaths.derivedTheoryFiles_YukawaSummed,
            loadImmediately=True
            )
        derivedTheoryFileContainers = TheoryFileInterface.FileFinder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFiles_YukawaSummed,
            loadImmediately=True
            )
        parametrization = Parametrization()
        parametrization.SetSM(SM)
        parametrization.Parametrize( derivedTheoryFileContainers )

        theoryBinBoundaries = SM.binBoundaries


        # ======================================
        # Get the combination result (center values + uncertainties)

        combinationscans = PhysicsCommands.GetScanResults(
            yieldParameterNames,
            scanDir,
            # pattern = 'combinedCard'
            )
        TgCombination = PhysicsCommands.GetTGraphForSpectrum( yieldParameterNames, combinationscans, name='Combination' )

        # Get relevant bins and errors from scan
        bestfitCenters  = TgCombination.POICenters[:nBinsExp]
        bestfitDownErrs = TgCombination.POIErrsLeft[:nBinsExp]
        bestfitUpErrs   = TgCombination.POIErrsRight[:nBinsExp]
        bestfitSymmErrs = TgCombination.POIErrsSymm[:nBinsExp]


        # ======================================
        # Get the correlation matrix

        corrMat = []
        corrMatFp = ROOT.TFile.Open(corrMatFile)
        fit = corrMatFp.Get('fit')
        for poi1 in yieldParameterNames:
            corrMatRow = []
            for poi2 in yieldParameterNames:
                corrMatRow.append( fit.correlation( poi1, poi2 ) )
            corrMat.append( corrMatRow )

        print '\nFound corrMat from {0}:'.format(corrMatFile)
        print numpy.array(corrMat)

        # Compute covariance matrix
        covMat = [ [ 0. for i in xrange(nBinsExp) ] for j in xrange(nBinsExp) ]
        for i in xrange(nBinsExp):
            for j in xrange(nBinsExp):
                covMat[i][j] = bestfitSymmErrs[i] * bestfitSymmErrs[j] * corrMat[i][j]

        print '\nComputed covMat:'
        print numpy.array(covMat)

        covMat_npArray = numpy.array(covMat)
        covMat_inversed_npArray = numpy.linalg.inv( covMat_npArray )


        # ======================================
        # Doing a bestfit

        # Build simple chi2 function
        def chi2_function( inputTuple ):
            kappac, kappab = inputTuple

            mus_parametrization = parametrization.EvaluateForBinning(
                theoryBinBoundaries, expBinning,
                kappab = kappab, kappac = kappac,
                returnRatios=True
                )

            x_column = numpy.array( [ mu - 1.0 for mu in mus_parametrization ] ).T

            # x_column = numpy.array( [ y.getVal() for y in yieldParameters ] ).T

            # chi2Val = 0.
            # for i in xrange(nBinsExp):
            #     chi2Val += (yieldParameterValues[i]-bestfitCenters[i])**2 / bestfitSymmErrs[i]**2

            chi2 = x_column.T.dot(  covMat_inversed_npArray.dot( x_column )  )

            return chi2


        # Do a best fit
        from scipy.optimize import minimize
        print ''
        res = minimize( chi2_function, [ 1., 1. ], method='Nelder-Mead', tol=1e-6 )
        print 'End of minimization'
        print res

        chi2_bestfit = res.fun
        kappac_bestfit = res.x[0]
        kappab_bestfit = res.x[1]


        # Do a scan

        kappacMin = -35.
        kappacMax = 35.
        kappabMin = -13.
        kappabMax = 13.

        kappacNPoints = 100
        kappabNPoints = 100


        kappacBinBoundaries = [ kappacMin + i*(kappacMax-kappacMin)/float(kappacNPoints) for i in xrange(kappacNPoints+1) ]
        kappabBinBoundaries = [ kappabMin + i*(kappabMax-kappabMin)/float(kappabNPoints) for i in xrange(kappabNPoints+1) ]

        kappacPoints = [ 0.5*(kappacBinBoundaries[i]+kappacBinBoundaries[i+1]) for i in xrange(kappacNPoints) ]
        kappabPoints = [ 0.5*(kappabBinBoundaries[i]+kappabBinBoundaries[i+1]) for i in xrange(kappabNPoints) ]


        H2 = ROOT.TH2F(
            'H2', '',
            kappacNPoints, array( 'f', kappacBinBoundaries ),
            kappabNPoints, array( 'f', kappabBinBoundaries ),
            )

        for i_kappab, kappabVal in enumerate(kappabPoints):
            for i_kappac, kappacVal in enumerate(kappacPoints):
                H2.SetBinContent( i_kappac, i_kappab, chi2_function( (kappacVal, kappabVal) ) - chi2_bestfit )


        print ''
        contours_1sigma = TheoryCommands.GetContoursFromTH2( H2, 2.30 )
        contours_2sigma = TheoryCommands.GetContoursFromTH2( H2, 6.18 )


        # ======================================
        # Plotting

        H2.SetTitle('')
        H2.GetXaxis().SetTitle( '#kappa_{c}' )
        H2.GetYaxis().SetTitle( '#kappa_{b}' )
        
        H2.SetMaximum(7.0)

        c.Clear()
        SetCMargins(
            LeftMargin   = 0.12,
            RightMargin  = 0.10,
            BottomMargin = 0.12,
            TopMargin    = 0.09,
            )

        H2.Draw('COLZ')

        for Tg in contours_1sigma:
            Tg.SetLineWidth(2)
            Tg.Draw('LSAME')
        for Tg in contours_2sigma:
            Tg.SetLineWidth(2)
            Tg.SetLineStyle(2)
            Tg.Draw('LSAME')


        SMpoint = ROOT.TGraph( 1, array( 'f', [1.0] ), array( 'f', [1.0] ) )
        SMpoint.SetMarkerStyle(21)
        SMpoint.SetMarkerSize(2)
        SMpoint.Draw('PSAME')

        bestfitpoint = ROOT.TGraph( 1, array( 'f', [kappac_bestfit] ), array( 'f', [kappab_bestfit] ) )
        bestfitpoint.SetMarkerStyle(34)
        bestfitpoint.SetMarkerSize(1.1)
        bestfitpoint.SetMarkerColor(13)
        bestfitpoint.Draw('PSAME')
        bestfitpoint.SetName('bestfitpoint')

        SaveC( 'AfterTheFactChi2Fit', asROOT = True )

        # rootFp.Close()





    #____________________________________________________________________
    if args.PlotParametrizationShapes:
        newColorCycle = lambda: itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )

        TheoryCommands.SetPlotDir( 'plots_{0}_Yukawa'.format(datestr) )

        # scanDir     = LatestPaths.scan_combined_PTH_xHfixed
        scanDir     = LatestPaths.scan_combined_PTH_xHfixed_asimov
        corrMatFile = 'corrMat_Nov08_combinedCard_Nov03_xHfixed/higgsCombine_CORRMAT_combinedCard_Nov03_xHfixed.MultiDimFit.mH125.root'


        # ======================================
        # Some settings

        expBinBoundaries = [
            0., 15., 30., 45., 85., 125.,
            # 200., 350., 10000.
            ]

        nBins = len(expBinBoundaries)-1
        binWidths = [ expBinBoundaries[i+1] - expBinBoundaries[i] for i in xrange(nBins) ]
        binCenters = [ 0.5*( expBinBoundaries[i] + expBinBoundaries[i+1] ) for i in xrange(nBins) ]


        # ======================================
        # Load parametrization from derived theory files

        SM = TheoryFileInterface.FileFinder(
            kappab=1, kappac=1, muR=1, muF=1, Q=1,
            expectOneFile=True,
            directory = LatestPaths.derivedTheoryFiles_YukawaSummed,
            loadImmediately=True
            )
        derivedTheoryFileContainers = TheoryFileInterface.FileFinder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFiles_YukawaSummed,
            loadImmediately=True
            )
        parametrization = Parametrization()
        # parametrization.SetSM(SM)

        # parametrization.Parametrize( derivedTheoryFileContainers )
        parametrization.ParametrizeByFitting( derivedTheoryFileContainers, fitWithScipy=True )

        theoryBinBoundaries = SM.binBoundaries
        theoryBinWidths  = [ right - left for left, right in zip( theoryBinBoundaries[:-1], theoryBinBoundaries[1:] ) ]
        theoryBinCenters = [ 0.5*( left + right ) for left, right in zip( theoryBinBoundaries[:-1], theoryBinBoundaries[1:] ) ]

        SM_integralFunction = TheoryCommands.GetIntegral(
            SM.binBoundaries,
            SM.crosssection
            )


        # ======================================
        # Scan for some points

        points = [

            #  ( kappac, kappab )
            # ( 1.0,    1.0 ),
            # ( 10.,    1.0 ),
            # ( -5.,    2.0 ),
            # ( 0.,     200. ),
            # ( 505.0,  505.0 ),
            # ( -495.0, 505. ),
            # ( 0.,     100000. )

            # ( 0., 0. ),
            # ( 1., 1. ),
            # ( 2., 2. ),
            # ( 4., 4. ),
            # ( 10., 10. ),
            # ( 100., 100. ),
            # ( 1000., 1000. ),
            # ( 10000., 10000. ),

            ( 10., 0. ),
            ( 10., 1./10. ),
            ( 10., 1./3. ),
            ( 10., 1./5. ),
            ( 10., 1. ),
            ( 10., 3. ),
            ( 10., 5. ),
            ( 10., 10. ),

            ( 10., -1./10. ),
            ( 10., -1./3. ),
            ( 10., -1./5. ),
            ( 10., -1. ),
            ( 10., -3. ),
            ( 10., -5. ),
            ( 10., -10. ),

            ( -10., 0. ),            

            ]


        plotContainers = []

        for kappac, kappab in points:

            print '\nkappac = {0}, kappab = {1}'.format( kappac, kappab )

            XSs_perGeV_parametrization = parametrization.EvaluateForBinning(
                theoryBinBoundaries, expBinBoundaries,
                kappab = kappab, kappac = kappac,
                returnRatios=False,
                verbose = True
                )


            print 'XSs_perGeV_parametrization: ', XSs_perGeV_parametrization

            for derivedTheoryFileContainer in derivedTheoryFileContainers:
                if derivedTheoryFileContainer.kappab == kappab and derivedTheoryFileContainer.kappac == kappac:

                    integral = TheoryCommands.GetIntegral(
                        derivedTheoryFileContainer.binBoundaries,
                        derivedTheoryFileContainer.crosssection
                        )

                    exp_xs = []
                    for left, right in zip( expBinBoundaries[:-1], expBinBoundaries[1:] ):
                        exp_xs.append( integral( left, right ) / ( right - left ) )

                    print '  for matching container:   ', exp_xs
                    break


            XSs_parametrization = [ xs * binWidth for xs, binWidth in zip( XSs_perGeV_parametrization, binWidths ) ]

            # Calculate the shape for the exp binning
            S_parametrization = [ xs / sum(XSs_parametrization) for xs in XSs_parametrization ]

            print 'XSs_parametrization:        ', XSs_parametrization
            print 'S_parametrization:          ', S_parametrization


            # Also calculate shape for theory binning
            XSs_perGeV_parametrization_theoryBinning = parametrization.Evaluate(
                kappab = kappab, kappac = kappac,
                returnRatios=False,
                # verbose = True
                )

            XSs_parametrization_theoryBinning = [
                xs * binWidth for xs, binWidth in zip( XSs_perGeV_parametrization_theoryBinning, theoryBinWidths )
                ]

            S_parametrization_theoryBinning = [ xs / sum(XSs_parametrization_theoryBinning) for xs in XSs_parametrization_theoryBinning ]


            container = Container()

            container.kappab = kappab
            container.kappac = kappac
            container.S_parametrization = S_parametrization
            container.S_parametrization_theoryBinning = S_parametrization_theoryBinning

            plotContainers.append( container )



        # ======================================
        # Make plot

        c.Clear()
        SetCMargins()

        xMin = expBinBoundaries[0]
        xMax = expBinBoundaries[-1]

        yMin = 0.0
        yMax = 0.5


        base = GetPlotBase(
            xMin = xMin,
            xMax = xMax,
            yMin = yMin,
            yMax = yMax,
            xTitle = 'p_{T} [GeV]', yTitle = 'Shape [A.U.]'
            )
        base.Draw('P')


        leg = ROOT.TLegend(
            c.GetLeftMargin(),
            1 - c.GetTopMargin() - 0.25,
            1 - c.GetRightMargin(),
            1 - c.GetTopMargin() 
            )
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetNColumns(3)


        colorCycle = newColorCycle()
        for container in plotContainers:
            color = colorCycle.next()

            # First the theory
            Tg_theory = ROOT.TGraph(
                len(theoryBinCenters),
                array( 'f', theoryBinCenters ),
                array( 'f', container.S_parametrization_theoryBinning )
                )
            ROOT.SetOwnership( Tg_theory, False )

            Tg_theory.SetMarkerColor(color)
            Tg_theory.SetMarkerStyle(8)
            Tg_theory.SetMarkerSize(0.9)
            Tg_theory.Draw('SAMEP')

            Tg_theory.SetName( TheoryCommands.GetUniqueRootName() )


            # H_exp = ROOT.TH1F(
            #     TheoryCommands.GetUniqueRootName(), '',

            #     )


            leg.AddEntry(
                Tg_theory.GetName(),
                '#kappa_{{c}} = {0:d}, #kappa_{{b}} = {1:d}'.format( int(container.kappac), int(container.kappab) ),
                'p'
                )

        leg.Draw()

        SaveC( 'ParametrizationShapes' )



    #____________________________________________________________________
    if args.RepeatTheoristFit:
        TheoryCommands.SetPlotDir( 'plots_{0}_Yukawa'.format(datestr) )

        # scanDir     = LatestPaths.scan_combined_PTH_xHfixed
        scanDir     = LatestPaths.scan_combined_PTH_xHfixed_asimov
        corrMatFile = 'corrMat_Nov08_combinedCard_Nov03_xHfixed/higgsCombine_CORRMAT_combinedCard_Nov03_xHfixed.MultiDimFit.mH125.root'


        # NORMALIZE_BY_SMXS = True
        NORMALIZE_BY_SMXS = False



        # ======================================
        # Load parametrization from derived theory files

        SM = TheoryFileInterface.FileFinder(
            kappab=1, kappac=1, muR=1, muF=1, Q=1,
            expectOneFile=True,
            directory = LatestPaths.derivedTheoryFiles_YukawaSummed,
            loadImmediately=True
            )
        derivedTheoryFileContainers = TheoryFileInterface.FileFinder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFiles_YukawaSummed,
            loadImmediately=True
            )
        parametrization = Parametrization()
        parametrization.SetSM(SM)
        parametrization.Parametrize( derivedTheoryFileContainers )

        theoryBinBoundaries = SM.binBoundaries

        SM_integralFunction = TheoryCommands.GetIntegral(
            SM.binBoundaries,
            SM.crosssection
            )


        # ======================================
        # Calculate SM cross section

        allBinBoundaries = [
            0., 15., 30., 45., 85., 125.,
            200., 350., 10000.
            ]

        nBins = len(allBinBoundaries)-1
        binWidths = [ allBinBoundaries[i+1] - allBinBoundaries[i] for i in xrange(nBins) ]

        # shape = [
        #     0.208025,
        #     0.234770,
        #     0.165146,
        #     0.218345,
        #     0.087552,
        #     0.059154,
        #     0.022612,
        #     0.004398
        #     ]
        # # shape = [ s / sum(shape) for s in shape ]
        # SMXSs = [ s * LatestPaths.YR4_totalXS for s in shape ]


        # Set up exp binning used for the fit
        expBinning = [ 0., 15., 30., 45., 85., 125. ]
        yieldParameterNames = []
        for left, right in zip( expBinning[:-1], expBinning[1:] ):
            yieldParameterNames.append( 'r_ggH_PTH_{0}_{1}'.format( int(left), int(right) ) )
        nBinsExp = len(expBinning) - 1

        # Re-normalize the shape
        # shape = shape[:nBinsExp]
        # shape = [ s / sum(shape) for s in shape ]

        SMXStot = SM_integralFunction( expBinning[0], expBinning[-1] )
        # SMXSs = [ s * SMXStot for s in shape ]

        # Actually, just use the SM from the parametrization, to at least get expected results back
        SMXSs = []
        for i in xrange(nBinsExp):
            SMXSs.append(
                SM_integralFunction( expBinning[i], expBinning[i+1] )
                )


        print '\nDetermined SM:'
        print '  SMXStot    = ', SMXStot
        # print '  shape      = ', shape
        print '  SMXSs      = ', SMXSs


        # ======================================
        # Get the combination result (center values + uncertainties)

        combinationscans = PhysicsCommands.GetScanResults(
            yieldParameterNames,
            scanDir,
            # pattern = 'combinedCard'
            )
        TgCombination = PhysicsCommands.GetTGraphForSpectrum( yieldParameterNames, combinationscans, name='Combination' )

        # Get relevant bins and errors from scan
        bestfitCenters  = TgCombination.POICenters[:nBinsExp]
        bestfitDownErrs = TgCombination.POIErrsLeft[:nBinsExp]
        bestfitUpErrs   = TgCombination.POIErrsRight[:nBinsExp]
        bestfitSymmErrs = TgCombination.POIErrsSymm[:nBinsExp]


        XSs_data = [ mu * SMXS for mu, SMXS in zip( bestfitCenters, SMXSs ) ]
        XS_data_tot = sum(XSs_data)

        # Input to theorist fit:
        S_data = [ xs / XS_data_tot for xs in XSs_data ]
        S_data_errs = [ e * s for e, s in zip( bestfitSymmErrs, S_data ) ]


        print '\nRead data from', scanDir
        print '  Found total XS: ', XS_data_tot
        print '  Found mus: ', bestfitCenters
        print '  Found xss: ', XSs_data
        print '  xs/xstot:  ', S_data
        print ''
        print '  Uncertainties:'
        print '  deltaMu:   ', bestfitSymmErrs
        print '  dxs/xstot: ', S_data_errs
        print ''


        # ======================================
        # Get the correlation matrix

        corrMat = []
        corrMatFp = ROOT.TFile.Open(corrMatFile)
        fit = corrMatFp.Get('fit')
        for poi1 in yieldParameterNames:
            corrMatRow = []
            for poi2 in yieldParameterNames:
                corrMatRow.append( fit.correlation( poi1, poi2 ) )
            corrMat.append( corrMatRow )

        print '\nFound corrMat from {0}:'.format(corrMatFile)
        print numpy.array(corrMat)


        print '\nOverwriting with square matrix for testing purposes:'
        for i in xrange(len(yieldParameterNames)):
            for j in xrange(len(yieldParameterNames)):
                if i==j:
                    corrMat[i][j] = 1.
                else:
                    corrMat[i][j] = 0.
        print numpy.array(corrMat)


        # Compute covariance matrix
        covMat = [ [ 0. for i in xrange(nBinsExp) ] for j in xrange(nBinsExp) ]
        for i in xrange(nBinsExp):
            for j in xrange(nBinsExp):
                covMat[i][j] = S_data_errs[i] * S_data_errs[j] * corrMat[i][j]

        print '\nComputed covMat:'
        print numpy.array(covMat)

        covMat_npArray = numpy.array(covMat)
        covMat_inversed_npArray = numpy.linalg.inv( covMat_npArray )


        # ======================================
        # Build simple chi2 function

        VERBOSITY_IN_CHI2 = True
        # VERBOSITY_IN_CHI2 = False
        
        f = lambda number, width = 6: '{0:+{width}.{decimals}f}'.format( number, width=width, decimals=width-4 )
        def chi2_function( inputTuple ):
            kappac, kappab = inputTuple

            XSs_perGeV_parametrization = parametrization.EvaluateForBinning(
                theoryBinBoundaries, expBinning,
                kappab = kappab, kappac = kappac,
                returnRatios=False,
                verbose = True if VERBOSITY_IN_CHI2 else False
                )

            XSs_parametrization = [ xs * binWidth for xs, binWidth in zip( XSs_perGeV_parametrization, binWidths ) ]

            # Calculate the shape for the exp binning
            if NORMALIZE_BY_SMXS:
                S_parametrization = [ xs / SMXStot for xs in XSs_parametrization ]
            else:
                S_parametrization = [ xs / sum(XSs_parametrization) for xs in XSs_parametrization ]


            # Calculate the column
            x_column = numpy.array( [ s_parametrization - s_data for s_parametrization, s_data in zip( S_parametrization, S_data ) ] ).T

            chi2 = x_column.T.dot(  covMat_inversed_npArray.dot( x_column )  )

            if VERBOSITY_IN_CHI2:
                print ''
                print '    kappac = {0}, kappab = {1}'.format( f(kappac), f(kappab) )
                print '    XSs_parametrization: ' + ' | '.join([ f(number) for number in XSs_parametrization ])
                print '    Shape_param:         ' + ' | '.join([ f(number) for number in S_parametrization ])
                print '    Shape_data:          ' + ' | '.join([ f(number) for number in S_data ])
                print '    Chi2: {0:.8f}'.format(chi2)

            return chi2


        # ======================================
        # Do a best fit
        
        from scipy.optimize import minimize
        print ''
        res = minimize( chi2_function, [ 1., 1. ], method='Nelder-Mead', tol=1e-6 )
        print 'End of minimization'
        print res

        chi2_bestfit = res.fun
        kappac_bestfit = res.x[0]
        kappab_bestfit = res.x[1]


        # ======================================
        # Do a scan

        if NORMALIZE_BY_SMXS:
            kappacMin = -45.
            kappacMax = 45.
            kappabMin = -23.
            kappabMax = 23.
        else:
            # kappacMin = -1000.
            # kappacMax = 1000.
            # kappabMin = -1000.
            # kappabMax = 1000.

            kappacMin = -100000.
            kappacMax = 100000.
            kappabMin = -100000.
            kappabMax = 100000.



        kappacNPoints       = 200
        kappabNPoints       = 200
        kappacBinBoundaries = [ kappacMin + i*(kappacMax-kappacMin)/float(kappacNPoints) for i in xrange(kappacNPoints+1) ]
        kappabBinBoundaries = [ kappabMin + i*(kappabMax-kappabMin)/float(kappabNPoints) for i in xrange(kappabNPoints+1) ]
        kappacPoints        = [ 0.5*(kappacBinBoundaries[i]+kappacBinBoundaries[i+1]) for i in xrange(kappacNPoints) ]
        kappabPoints        = [ 0.5*(kappabBinBoundaries[i]+kappabBinBoundaries[i+1]) for i in xrange(kappabNPoints) ]

        H2 = ROOT.TH2F(
            'H2', '',
            kappacNPoints, array( 'f', kappacBinBoundaries ),
            kappabNPoints, array( 'f', kappabBinBoundaries ),
            )


        print '\n\nDoing scan'

        # VERBOSITY_IN_CHI2 = True
        VERBOSITY_IN_CHI2 = False
        iIteration = 0

        for i_kappab, kappabVal in enumerate(kappabPoints):
            for i_kappac, kappacVal in enumerate(kappacPoints):

                if i_kappab % 50 == 0 and i_kappac % 50 == 0:
                    # print '  ', i_kappab, i_kappac
                    VERBOSITY_IN_CHI2 = True

                H2.SetBinContent( i_kappac, i_kappab, chi2_function( (kappacVal, kappabVal) ) - chi2_bestfit )

                VERBOSITY_IN_CHI2 = False


        print ''
        contours_1sigma = TheoryCommands.GetContoursFromTH2( H2, 2.30 )
        contours_2sigma = TheoryCommands.GetContoursFromTH2( H2, 6.18 )


        # ======================================
        # Plotting

        H2.SetTitle('')
        H2.GetXaxis().SetTitle( '#kappa_{c}' )
        H2.GetYaxis().SetTitle( '#kappa_{b}' )
        H2.SetMaximum(7.0)

        c.Clear()
        SetCMargins(
            LeftMargin   = 0.12,
            RightMargin  = 0.10,
            BottomMargin = 0.12,
            TopMargin    = 0.09,
            )
        H2.Draw('COLZ')

        for Tg in contours_1sigma:
            Tg.SetLineWidth(2)
            Tg.Draw('LSAME')
        for Tg in contours_2sigma:
            Tg.SetLineWidth(2)
            Tg.SetLineStyle(2)
            Tg.Draw('LSAME')

        SMpoint = ROOT.TGraph( 1, array( 'f', [1.0] ), array( 'f', [1.0] ) )
        SMpoint.SetMarkerStyle(21)
        SMpoint.SetMarkerSize(2)
        SMpoint.Draw('PSAME')

        bestfitpoint = ROOT.TGraph( 1, array( 'f', [kappac_bestfit] ), array( 'f', [kappab_bestfit] ) )
        bestfitpoint.SetMarkerStyle(34)
        bestfitpoint.SetMarkerSize(1.1)
        bestfitpoint.SetMarkerColor(13)
        bestfitpoint.Draw('PSAME')
        bestfitpoint.SetName('bestfitpoint')

        SaveC( 'TheoristChi2Fit', asROOT = True )


    #____________________________________________________________________
    if args.Make2DPlotOfVariableInWS:
        TheoryCommands.SetPlotDir( 'plots_{0}_Yukawa'.format(datestr) )


        # ws = LatestPaths.ws_combined_split_betterYukawa_couplingDependentBR
        ws = LatestPaths.ws_combined_yukawa_couplingDependentBR

        # varToPlot = 'hggBRmodifier'
        # varToPlot = 'hzzBRmodifier'
        varToPlot = 'c7_BRscal_hgg'
        # varToPlot = 'c7_BRscal_hzz'
        # varToPlot = 'Scaling_hgluglu'


        xVar = 'kappac'
        xRange = [ -32., 32. ]
        xNBins = 200

        yVar = 'kappab'
        yRange = [ -15., 15. ]
        yNBins = 200


        # ======================================
        # 

        rootFp = ROOT.TFile.Open( ws )
        w = rootFp.Get('w')

        x = w.var(xVar)
        y = w.var(yVar)

        mH = w.var('MH')
        mH.setVal( 125. )

        for getter in [ 'var', 'function' ]:
            z = getattr( w, getter )(varToPlot)
            try:
                z.getVal()
                break
            except ReferenceError:
                pass
            except TypeError:
                pass

        else:
            print 'Variable \'{0}\' is probably not in the workspace'.format(varToPlot)
            sys.exit()



        def makeBinBoundaries( xMin, xMax, nBins ):
            return [ xMin + i*(xMax-xMin)/nBins for i in xrange(nBins+1) ]

        xBinBoundaries = makeBinBoundaries( xRange[0], xRange[1], xNBins )
        yBinBoundaries = makeBinBoundaries( yRange[0], yRange[1], yNBins )

        xBinCenters = [ 0.5*(xBinBoundaries[i]+xBinBoundaries[i+1]) for i in xrange(xNBins) ]
        yBinCenters = [ 0.5*(yBinBoundaries[i]+yBinBoundaries[i+1]) for i in xrange(yNBins) ]


        H2 = ROOT.TH2F(
            'H2', '',
            xNBins, array( 'f', xBinBoundaries ),
            yNBins, array( 'f', yBinBoundaries ),
            )

        H2.SetTitle(varToPlot)
        H2.GetXaxis().SetTitle( xVar )
        H2.GetYaxis().SetTitle( yVar )

        H2.SetMaximum(1.5)

        for ix in xrange(xNBins):
            for iy in xrange(yNBins):

                xCenter = xBinCenters[ix]
                yCenter = yBinCenters[iy]

                x.setVal(xCenter)
                y.setVal(yCenter)

                zValue = z.getVal()

                H2.SetBinContent( ix, iy, zValue )



        c.Clear()
        SetCMargins(
            LeftMargin   = 0.12,
            RightMargin  = 0.10,
            BottomMargin = 0.12,
            TopMargin    = 0.09,
            )

        H2.Draw('COLZ')

        SaveC( '2Dplot_{0}'.format(varToPlot) )


        rootFp.Close()




    #____________________________________________________________________
    if args.PlotBRsInOnePlot:
        TheoryCommands.SetPlotDir( 'plots_{0}_Yukawa'.format(datestr) )

        PLOT_SCALING = False

        LOGSCALE = True
        # LOGSCALE = False


        # ======================================
        # 

        # ws = LatestPaths.ws_combined_split_betterYukawa_couplingDependentBR_profiledTotalXS
        # ws = LatestPaths.ws_combined_split_betterYukawa_couplingDependentBR
        # ws = LatestPaths.ws_combined_yukawa_couplingDependentBR
        ws = LatestPaths.ws_combined_Yukawa_couplingDependentBR

        rootFp = ROOT.TFile.Open( ws )
        w = rootFp.Get('w')

        mH = w.var('MH')
        mH.setVal(125.)


        titleDict = {
            'kappab' : '#kappa_{b}',
            'kappac' : '#kappa_{c}',
            'kappa_V' : '#kappa_{V}',
            }


        SMBR = {
            'hww'     : 2.137E-01,
            'hzz'     : 2.619E-02,
            'htt'     : 6.272E-02,
            'hmm'     : 2.176E-04,
            'hbb'     : 5.824E-01,
            'hcc'     : 2.891E-02,
            'hgg'     : 2.270E-03,
            'hzg'     : 1.533E-03,
            'hgluglu' : 8.187E-02,
            }

        SMcouplings = {
            'kappab' : 1.0,
            'kappac' : 1.0,
            'kappa_V' : 1.0,
            }
        couplings = SMcouplings.keys()

        BRfunctions = {
            'hww'     : 'c7_BRscal_hww',
            'hzz'     : 'c7_BRscal_hzz',
            'htt'     : 'c7_BRscal_htt',
            'hmm'     : 'c7_BRscal_hmm',
            'hbb'     : 'c7_BRscal_hbb',
            'hcc'     : 'c7_BRscal_hcc',
            'hgg'     : 'c7_BRscal_hgg',
            'hzg'     : 'c7_BRscal_hzg',
            'hgluglu' : 'c7_BRscal_hgluglu',
            }

        decaysToPlot = BRfunctions.keys()


        for xVar in couplings:

            x = w.var(xVar)

            # Set other couplings to their SM value
            for otherCoupling in couplings:
                if otherCoupling == xVar: continue
                w.var(otherCoupling).setVal( SMcouplings[otherCoupling] )

            if xVar == 'kappab':
                kappab_min = -10
                kappab_max = 10
            elif xVar == 'kappac':
                kappab_min = -30
                kappab_max = 30
            elif xVar == 'kappa_V':
                kappab_min = -15
                kappab_max = 15

            # ======================================
            # From here on, 'kappab' is the xVar

            kappab = x

            nPoints = 100
            kappab_axis = [ kappab_min + i*(kappab_max-kappab_min)/(nPoints-1.) for i in xrange(nPoints) ]

            BRvalues = { d : [] for d in decaysToPlot }
            for decay in decaysToPlot:
                y = w.function( BRfunctions[decay] )
                for kappab_val in kappab_axis:
                    kappab.setVal(kappab_val)
                    BRvalues[decay].append(y.getVal())


            # ======================================
            # Plotting

            c.Clear()
            SetCMargins()

            yMin = -0.1
            yMax = 5.0 if PLOT_SCALING else 1.05

            if LOGSCALE:
                yMin = 0.0001
                c.SetLogy(True)

            base = GetPlotBase(
                xMin = kappab_min,
                xMax = kappab_max,
                yMin = yMin,
                yMax = yMax,
                xTitle = titleDict.get( xVar, xVar ),
                yTitle = 'BR scaling' if PLOT_SCALING else 'BR',
                )
            base.Draw('P')

            leg = ROOT.TLegend(
                1 - c.GetRightMargin() - 0.35,
                1 - c.GetTopMargin() - 0.50,
                1 - c.GetRightMargin(),
                1 - c.GetTopMargin() 
                )
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)


            colorCycle = itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )
            for decay in decaysToPlot:
                color = next(colorCycle)

                if PLOT_SCALING:
                    yValues = BRvalues[decay]
                else:
                    yValues = [ SMBR[decay] * val for val in BRvalues[decay] ]

                Tg = ROOT.TGraph(
                    nPoints,
                    array( 'd', kappab_axis ),
                    array( 'd', yValues ),
                    )
                ROOT.SetOwnership( Tg, False )
                Tg.SetLineWidth(2)
                Tg.SetLineColor(color)
                Tg.Draw('LSAME')
                Tg.SetName(decay)

                leg.AddEntry( Tg.GetName(), Tg.GetName(), 'l' )

            if not PLOT_SCALING:

                sumY = [ 0. for i in xrange(nPoints) ]
                for i in xrange(nPoints):
                    for d in decaysToPlot:
                        sumY[i] +=  SMBR[d] * BRvalues[d][i]

                sumTg = ROOT.TGraph(
                    nPoints,
                    array( 'd', kappab_axis ),
                    array( 'd', sumY ),
                    )
                ROOT.SetOwnership( sumTg, False )
                sumTg.SetLineWidth(4)
                sumTg.SetLineStyle(2)
                sumTg.SetLineColor(1)
                sumTg.Draw('LSAME')


            lineAtOne = ROOT.TLine( 1.0, yMin, 1.0, yMax )
            lineAtOne.Draw('L')

            lineAtKappabOne = ROOT.TLine( kappab_min, 1.0, kappab_max, 1.0 )
            lineAtKappabOne.Draw('L')

            leg.Draw()

            SaveC(
                ( 'BRscalings' if PLOT_SCALING else 'BRs' )
                + ( '_logscale' if LOGSCALE else '' )
                + '-' + xVar
                )
            rootFp.Close()





    #____________________________________________________________________
    if args.PlotOfTotalXSInYukawaWS:
        TheoryCommands.SetPlotDir( 'plots_{0}_Yukawa'.format(datestr) )

        # ws = LatestPaths.ws_combined_split_betterYukawa_couplingDependentBR_profiledTotalXS
        # ws = LatestPaths.ws_combined_split_betterYukawa_couplingDependentBR
        ws = LatestPaths.ws_combined_yukawa_couplingDependentBR

        rootFp = ROOT.TFile.Open( ws )
        w = rootFp.Get('w')


        for yVar in [ 'totalXS', 'Scaling_hgluglu' ]:

            y = w.function(yVar)

            kappab = w.var('kappab')
            kappac = w.var('kappac')
            mH     = w.var('MH')

            mH.setVal(125.)

            y_SM = y.getVal()


            minKappab = -5.
            maxKappab = 15.
            nPoints = 100

            kappabAxis = [ minKappab + i*(maxKappab-minKappab)/(nPoints-1.) for i in xrange(nPoints) ]

            yAxis = []
            for kappabVal in kappabAxis:
                kappab.setVal(kappabVal)

                yVal = y.getVal()

                if yVar == 'totalXS':
                    yAxis.append( yVal / y_SM )
                else:
                    yAxis.append( yVal )


            c.Clear()
            SetCMargins()


            # yMin = min(yAxis) - 0.1*(max(yAxis)-min(yAxis)),
            # yMax = max(yAxis) + 0.1*(max(yAxis)-min(yAxis)),

            yMin = 0.7
            yMax = 2.3

            base = GetPlotBase(
                xMin = minKappab,
                xMax = maxKappab,
                yMin = yMin,
                yMax = yMax,
                xTitle = '#kappa_{b}',
                yTitle = yVar,
                )
            base.Draw('P')

            base.SetTitle( '{0}_byEvaluatingFunctionInWorkspace'.format(yVar) )


            Tg = ROOT.TGraph(
                nPoints,
                array( 'd', kappabAxis ),
                array( 'd', yAxis ),
                )
            Tg.SetLineWidth(2)
            Tg.SetLineColor(2)
            Tg.Draw('LSAME')

            lineAtOne = ROOT.TLine( minKappab, 1.0, maxKappab, 1.0 )
            lineAtOne.Draw('L')


            SaveC( 'kappabDependency_{0}_byEvaluatingFunctionInWorkspace'.format(yVar) )

            rootFp.Close()


            if True:

                minY = min(yAxis)
                minX = kappabAxis[ yAxis.index(minY) ]

                print 'minimum: kappab = {0}, {1} = {2}'.format( minX, yVar, minY )

                kappab.setVal(0.)
                y_0 = y.getVal() / y_SM

                print 'kappab = 0.0, {0} = {1}'.format( yVar, y_0 )



    #____________________________________________________________________
    if args.PlotOfTotalXS_FromParametrization:
        TheoryCommands.SetPlotDir( 'plots_{0}_Yukawa'.format(datestr) )


        # ======================================
        # Gather all the theory shapes

        readList = lambda theoryFiles: [ TheoryFileInterface.ReadDerivedTheoryFile( theoryFile ) for theoryFile in theoryFiles ]

        theoryFiles_kappab_kappac = TheoryFileInterface.FileFinder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFilesDirectory_YukawaSummed
            )
        containers_kappab_kappac = readList( theoryFiles_kappab_kappac )


        theoryFiles_kappab_kappac_Gluon = TheoryFileInterface.FileFinder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFilesDirectory_YukawaGluonInduced
            )
        containers_kappab_kappac_Gluon = readList( theoryFiles_kappab_kappac_Gluon )


        theoryFiles_kappab_kappac_Quark = TheoryFileInterface.FileFinder(
            kappab='*', kappac='*', filter='muR',
            directory = LatestPaths.derivedTheoryFilesDirectory_YukawaQuarkInduced
            )
        containers_kappab_kappac_Quark = readList( theoryFiles_kappab_kappac_Quark )


        theoryFiles_kappab_kappac_QuarkScaled = TheoryFileInterface.FileFinder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFilesDirectory_YukawaQuarkInducedScaled
            )
        containers_kappab_kappac_QuarkScaled = readList( theoryFiles_kappab_kappac_QuarkScaled )


        # ======================================
        # 

        for key in [ 'Gluon', 'Quark', 'QuarkScaled', 'Summed' ]:

            if key == 'Gluon':
                containers = containers_kappab_kappac_Gluon
                parametrizeByFitting = False
            elif key == 'Quark':
                containers = containers_kappab_kappac_Quark
                parametrizeByFitting = True
            elif key == 'QuarkScaled':
                containers = containers_kappab_kappac_QuarkScaled
                parametrizeByFitting = False
            elif key == 'Summed':
                containers = containers_kappab_kappac
                parametrizeByFitting = False


            SM = [ container for container in containers if container.kappab == 1 and container.kappac == 1 and container.muR == 1 and container.muF ==1 and container.Q == 1 ][0]

            parametrization = Parametrization()
            if parametrizeByFitting:
                parametrization.ParametrizeByFitting( containers, fitWithScipy=True )
            else:
                parametrization.Parametrize( containers )
            parametrization.kappac = 1.0

            minKappab = -5.
            maxKappab = 15.
            nPoints = 100
            kappabAxis = [ minKappab + i*(maxKappab-minKappab)/(nPoints-1.) for i in xrange(nPoints) ]
            totalXSaxis = []

            for kappabVal in [ 0., 1. ] + kappabAxis:
                parametrization.kappab = kappabVal
                parametrizationResult = parametrization.GetOutputContainer()

                # Calculate integral
                binWidths = [ 0.5*(parametrizationResult.binBoundaries[i]+parametrizationResult.binBoundaries[i+1]) for i in xrange(len(parametrizationResult.binBoundaries)-1) ]
                integral = sum([
                    width * ratio * SMXS for width, ratio, SMXS in zip( binWidths, parametrizationResult.ratios, SM.crosssection )
                    ])

                totalXSaxis.append( integral )

            totalXS_at_kappab_zero = totalXSaxis.pop(0)
            totalXS_at_kappab_one  = totalXSaxis.pop(0)

            totalXSaxis = [ xs / totalXS_at_kappab_one for xs in totalXSaxis ]


            c.Clear()
            SetCMargins()


            # yMin = min(totalXSaxis) - 0.1*(max(totalXSaxis)-min(totalXSaxis)),
            # yMax = max(totalXSaxis) + 0.1*(max(totalXSaxis)-min(totalXSaxis)),

            yMin = 0.7
            yMax = 2.3

            base = GetPlotBase(
                xMin = minKappab,
                xMax = maxKappab,
                yMin = yMin,
                yMax = yMax,
                xTitle = '#kappa_{b}',
                yTitle = 'totalXS',
                )
            base.Draw('P')

            base.SetTitle('totalXS')

            Tg = ROOT.TGraph(
                nPoints,
                array( 'd', kappabAxis ),
                array( 'd', totalXSaxis ),
                )
            Tg.SetLineWidth(2)
            Tg.SetLineColor(2)
            Tg.Draw('LSAME')

            lineAtOne = ROOT.TLine( minKappab, 1.0, maxKappab, 1.0 )
            lineAtOne.Draw('L')


            SaveC( 'kappabDependency_totalXS_fromHistograms_{0}'.format(key) )


            if True:

                minY = min(totalXSaxis)
                minX = kappabAxis[ totalXSaxis.index(minY) ]

                print 'minimum: kappab = {0}, {1} = {2}'.format( minX, 'totalXS', minY )

                print 'kappab = 0.0, {0} = {1}'.format( 'totalXS', totalXS_at_kappab_zero / totalXS_at_kappab_one )


    #____________________________________________________________________
    if args.CorrelationMatrixScaleDependence_Yukawa:

        CorrelationMatrices.SetPlotDir( 'plots_CorrelationMatrixCrosscheck_{0}'.format(datestr) )
        TheoryCommands.SetPlotDir( 'plots_CorrelationMatrixCrosscheck_{0}'.format(datestr) )

        expCorrMatrices = []
        theoryCorrMatrices = []

        kappabs = [ -2, -1, 0, 1, 2 ]
        kappacs = [ -10, -5, 0, 1, 5, 10 ]
        # kappabs = [ -1, 0, 1 ]
        # kappacs = [ -5, 1, 5 ]

        for kappab, kappac in itertools.product( kappabs, kappacs ):
            print 'Processing kappab = {0}, kappac = {1}'.format( kappab, kappac )

            variationFiles = TheoryFileInterface.FileFinder(
                directory = LatestPaths.derivedTheoryFiles_YukawaGluonInduced,
                # directory = LatestPaths.derivedTheoryFilesDirectory_YukawaSummed,
                kappab = kappab, kappac = kappac
                )

            variations = [
                TheoryFileInterface.ReadDerivedTheoryFile( variationFile, returnContainer=True )
                    for variationFile in variationFiles ]

            theoryCorrMatrices.append( CorrelationMatrices.GetCorrelationMatrix(
                variations,
                makeScatterPlots          = False,
                makeCorrelationMatrixPlot = True,
                outname                   = 'corrMat_theory_kappab_{0}_kappac_{1}'.format( kappab, kappac ),
                verbose                   = True,
                ))

            variations_expbinning = deepcopy(variations)
            for variation in variations_expbinning:
                TheoryCommands.RebinDerivedTheoryContainer( variation, [ 0., 15., 30., 45., 85., 125. ] )

            expCorrMatrices.append( CorrelationMatrices.GetCorrelationMatrix(
                variations_expbinning,
                makeScatterPlots          = False,
                makeCorrelationMatrixPlot = True,
                outname                   = 'corrMat_exp_kappab_{0}_kappac_{1}'.format( kappab, kappac ),
                verbose                   = True,
                ))


        # ======================================
        # Min-max study

        print '\nCreating min-max matrix'

        toString = lambda number: str(int(number)) if number.is_integer() else '{0:.1f}'.format(number)
        def getBinLabels( variation ):
            binLabels = []
            for iBin in xrange(len(variation.binCenters)):
                if not hasattr( variation, 'binBoundaries' ):
                    binLabel = toString(variation.binCenters[iBin])
                else:
                    binLabel = '{0} - {1}'.format(
                        toString(variation.binBoundaries[iBin] ),
                        toString(variation.binBoundaries[iBin+1] )
                        )
                binLabels.append(binLabel)
            return binLabels

        expBinLabels    = getBinLabels( variations_expbinning[0] )
        theoryBinLabels = getBinLabels( variations[0] )


        for name, corrMatrices, binLabels in [
                ( 'exp', expCorrMatrices, expBinLabels ),
                ( 'theory', theoryCorrMatrices, theoryBinLabels )
                ]:

            nBins = len(corrMatrices[0])

            # Construct minimum and maximum corrMatrix
            minCorrMatrix = [ [ 999.  for j in xrange(nBins) ] for i in xrange(nBins) ]
            maxCorrMatrix = [ [ -999. for j in xrange(nBins) ] for i in xrange(nBins) ]
            degreeOfAsymmetryMatrix = [ [ 0. for j in xrange(nBins) ] for i in xrange(nBins) ]
            for iRow in xrange(nBins):
                for iCol in xrange(nBins):
                    minVal = 999.
                    maxVal = -999.
                    for corrMatrix in corrMatrices:
                        if corrMatrix[iRow][iCol] < minVal:
                            minVal = corrMatrix[iRow][iCol]
                        if corrMatrix[iRow][iCol] > maxVal:
                            maxVal = corrMatrix[iRow][iCol]
                    minCorrMatrix[iRow][iCol] = minVal
                    maxCorrMatrix[iRow][iCol] = maxVal
                    degreeOfAsymmetryMatrix[iRow][iCol] = maxVal - minVal


            # Make a plot

            c.Clear()
            SetCMargins(
                LeftMargin   = 0.18,
                RightMargin  = 0.12,
                TopMargin    = 0.06,
                BottomMargin = 0.16,
                )

            T = ROOT.TH2D(
                'corrMat', '#scale[0.85]{Min and max corr. for p_{T} bins}',
                nBins, 0., nBins,
                nBins, 0., nBins
                )
            ROOT.SetOwnership( T, False )

            for iRow in xrange(nBins):
                for iCol in xrange(nBins):
                    T.SetBinContent( iCol+1, iRow+1, degreeOfAsymmetryMatrix[iRow][iCol] )

            # Bin titles
            for iBin in xrange(nBins):
                if nBins < 20 or iBin % int(0.1*nBins) == 0:
                    T.GetXaxis().SetBinLabel( iBin+1, binLabels[iBin] )
                    T.GetYaxis().SetBinLabel( iBin+1, binLabels[iBin] )

            # Set color range

            # # With negative (to orange)
            # n_stops = 3
            # stops  = [ 0.0, 0.5, 1.0 ]
            # reds   = [ 0.0, 1.0, 1.0 ]
            # blues  = [ 1.0, 1.0, 0.0 ]
            # greens = [ 0.0, 1.0, 0.4 ]
            # zMin   = -1.2
            # zMax   = 1.2

            # With only positive
            n_stops = 2
            stops  = [ 0.0, 1.0 ]
            reds   = [ 1.0, 255./255. ]
            greens = [ 1.0, 191./255. ]
            blues  = [ 1.0, 128./255. ]
            zMin   = 0.
            zMax   = 1.05 * max([ max(row) for row in degreeOfAsymmetryMatrix ])


            ROOT.TColor.CreateGradientColorTable(
                n_stops,
                array('d', stops ),
                array('d', reds ),
                array('d', greens ),
                array('d', blues ),
                255 )
            T.GetZaxis().SetRangeUser( zMin, zMax )

            T.GetXaxis().SetTitle( 'p_{T}' )
            # T.GetXaxis().SetTitleOffset( 1.7 )
            T.GetXaxis().SetTitleOffset( 1.2 )
            T.GetXaxis().SetTitleSize(0.05)

            T.GetYaxis().SetTitle( 'p_{T}' )
            # T.GetYaxis().SetTitleOffset( 2.45 )
            T.GetYaxis().SetTitleOffset( 1.85 )
            T.GetYaxis().SetTitleSize(0.05)

            T.GetXaxis().SetLabelSize(0.045)
            T.GetYaxis().SetLabelSize(0.045)

            T.Draw('COLZ')

            c.cd()
            c.Update()


            # Draw in the min and max numbers

            labelMin = ROOT.TLatex()
            labelMin.SetTextSize(0.035)
            labelMin.SetTextColor(4)
            labelMin.SetTextAlign(23)

            labelMax = ROOT.TLatex()
            labelMax.SetTextSize(0.035)
            labelMax.SetTextColor(2)
            labelMax.SetTextAlign(21)

            xMin = T.GetXaxis().GetBinLowEdge(1)
            xMax = T.GetXaxis().GetBinUpEdge(nBins)
            yMin = T.GetYaxis().GetBinLowEdge(1)
            yMax = T.GetYaxis().GetBinUpEdge(nBins)
            dy = yMax - yMin

            significance = 2
            formatter = lambda number: '{0:+.{significance}f}'.format(
                number, significance=significance )

            if nBins >= 20:
                labelMin.SetTextSize(0.008)
                labelMax.SetTextSize(0.008)
                dy *= 0.0
                significance = 1
                formatter = lambda number: '{0:+.{significance}f}'.format(
                    number, significance=significance ).replace('+0','+').replace('-0','-').replace('1.0','1')


            for iRow in xrange(nBins):
                for iCol in xrange(nBins):

                    if (
                        nBins < 20
                        or ( iRow % int(0.1*nBins) == 0 and iCol % int(0.1*nBins) == 0 )
                        or ( degreeOfAsymmetryMatrix[iRow][iCol] > 0.5 )
                        ):

                        x = T.GetXaxis().GetBinCenter(iCol+1)
                        y = T.GetYaxis().GetBinCenter(iRow+1)

                        labelMin.DrawLatex(
                            x, y - 0.01*dy,
                            formatter( minCorrMatrix[iRow][iCol] )
                            )

                        labelMax.DrawLatex(
                            x, y + 0.01*dy,
                            formatter( maxCorrMatrix[iRow][iCol] )
                            )

            SaveC( '_minmax_corrMat_{0}'.format(name), asPNG=True, asROOT=True )



########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'