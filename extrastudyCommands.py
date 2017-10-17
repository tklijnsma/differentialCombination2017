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

    parser.add_argument( '--FitBR_t2ws',                             action=CustomAction )
    parser.add_argument( '--FitBR_bestfit',                          action=CustomAction )
    parser.add_argument( '--FitBR_scan',                             action=CustomAction )
    parser.add_argument( '--FitBR_plot',                             action=CustomAction )

    parser.add_argument( '--TotalXS_t2ws',                           action=CustomAction )
    parser.add_argument( '--TotalXS_bestfit',                        action=CustomAction )
    parser.add_argument( '--TotalXS_scan',                           action=CustomAction )
    parser.add_argument( '--TotalXS_plot',                           action=CustomAction )

    parser.add_argument( '--chi2fitToCombination_Yukawa',            action=CustomAction )

    parser.add_argument( '--Make2DPlotOfVariableInWS',               action=CustomAction )
    parser.add_argument( '--PlotOfTotalXSInYukawaWS',                action=CustomAction )
    parser.add_argument( '--PlotOfTotalXS_FromParametrization',      action=CustomAction )


########################################
# Methods
########################################    

def main( args ):

    TheoryCommands.SetPlotDir( 'plots_{0}_extraStudy'.format(datestr) )


    #____________________________________________________________________
    if args.FitBR_t2ws:

        extraOptions = [
            '--PO verbose=2',
            '--PO \'higgsMassRange=123,127\'',
            ]

        Commands.BasicT2WSwithModel(
            LatestPaths.card_combined_unsplit,
            'FitBRModel.py',
            # suffix       = '',
            extraOptions = extraOptions,
            smartMaps    = [
                ( r'.*/smH_PTH_([\d\_GT]+)', r'r_smH_PTH_\1[1.0,-1.0,4.0]' )
                ],
            )

    #____________________________________________________________________
    if args.FitBR_bestfit:

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

        doFastscan = True
        if args.notFastscan: doFastscan = False

        # ASIMOV = True
        ASIMOV = False
        
        datacard = LatestPaths.ws_FitBR_combined_unsplit

        ratio_BR_hgg_hzz_ranges = [ 0.03, 0.16 ]

        jobDirectory = 'Scan_BR_{0}'.format( datestr )
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
            jobPriority   = -5,
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

        scanRootFiles = glob( 'Scan_BR_Sep25/*.root' )

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


        scanContainer.GetTGraph( xAttr = 'x', yAttr = 'y', xAreBinBoundaries = False )
        scanContainer.Tg.SetMarkerStyle(8)
        scanContainer.Tg.SetMarkerSize(0.8)


        # ======================================
        # Make plot

        c.Clear()
        SetCMargins()

        yMinAbs = min( scanContainer.y )
        yMaxAbs = max( scanContainer.y )
        yMin = yMinAbs - 0.1*(yMaxAbs-yMinAbs)
        yMax = yMaxAbs + 0.5*(yMaxAbs-yMinAbs)

        xMin = min( scanContainer.x )
        xMax = max( scanContainer.x )

        base = GetPlotBase(
            xMin = xMin,
            xMax = xMax,
            yMin = yMin,
            yMax = yMax,
            xTitle = 'BR_{H #rightarrow #gamma#gamma} / BR_{H #rightarrow ZZ}',
            yTitle = '#Delta NLL',
            )
        base.Draw('P')


        scanContainer.Tg.Draw('XPL')


        oneLine = ROOT.TLine( xMin, 0.5, xMax, 0.5 )
        oneLine.SetLineColor(12)
        oneLine.Draw()

        scanContainer.extrema = PhysicsCommands.FindMinimaAndErrors( scanContainer.x, scanContainer.y, returnContainer=True )

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
        SM_ratio = 0.086
        SMLine = ROOT.TLine( SM_ratio, yMin, SM_ratio, yMax )
        SMLine.SetLineWidth(2)
        SMLine.SetLineColor(9)
        SMLine.Draw()


        xBestfit = scanContainer.x[ scanContainer.extrema.imin ]
        bestfitLine = ROOT.TLine( xBestfit, yMin, xBestfit, yMax )
        bestfitLine.SetLineWidth(2)
        bestfitLine.SetLineColor(2)
        bestfitLine.Draw()


        SaveC( 'BRscan' )




    #____________________________________________________________________
    if args.TotalXS_t2ws:

        extraOptions = [
            '--PO verbose=2',
            '--PO \'higgsMassRange=123,127\'',
            ]

        Commands.BasicT2WSwithModel(
            LatestPaths.card_combined_unsplit,
            pathToModel = 'FitBRModel.py',
            modelName   = 'fitTotalXSModel',
            suffix       = 'fitTotalXS',
            extraOptions = extraOptions,
            # smartMaps    = [
            #     ( r'.*/smH_PTH_([\d\_GT]+)', r'r_smH_PTH_\1[1.0,-1.0,4.0]' )
            #     ],
            )

    #____________________________________________________________________
    if args.TotalXS_bestfit:

        print 'Note: Move workspace to LatestPaths.py'
        ws = abspath( 'workspaces_Oct05/combinedCard_Jul26_FitBRModel_fitTotalXS.root' )

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

        print 'Note: Move workspace to LatestPaths.py'
        ws = abspath( 'workspaces_Oct05/combinedCard_Jul26_FitBRModel_fitTotalXS.root' )

        doFastscan = True
        if args.notFastscan: doFastscan = False

        # ASIMOV = True
        ASIMOV = False
        
        datacard = ws

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
            nPoints = 42
            nPointsPerJob = 3
            queue = 'short.q'

        Commands.MultiDimCombineTool(
            datacard,
            nPoints       = nPoints,
            nPointsPerJob = nPointsPerJob,
            queue         = queue,
            notOnBatch    = True,
            jobDirectory  = jobDirectory,
            fastscan      = doFastscan,
            asimov        = ASIMOV,
            jobPriority   = -5,
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

        scanRootFiles = glob( 'Scan_TotalXS_Oct05_0/*.root' )

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


        scanContainer.GetTGraph( xAttr = 'x', yAttr = 'y', xAreBinBoundaries = False )
        scanContainer.Tg.SetMarkerStyle(8)
        scanContainer.Tg.SetMarkerSize(0.8)


        # ======================================
        # Make plot

        c.Clear()
        SetCMargins()

        yMinAbs = min( scanContainer.y )
        yMaxAbs = max( scanContainer.y )
        yMin = yMinAbs - 0.1*(yMaxAbs-yMinAbs)
        yMax = yMaxAbs + 0.5*(yMaxAbs-yMinAbs)

        xMin = min( scanContainer.x )
        xMax = max( scanContainer.x )

        base = GetPlotBase(
            xMin = xMin,
            xMax = xMax,
            yMin = yMin,
            yMax = yMax,
            xTitle = '#sigma',
            yTitle = '#Delta NLL',
            )
        base.Draw('P')


        scanContainer.Tg.Draw('XPL')


        oneLine = ROOT.TLine( xMin, 0.5, xMax, 0.5 )
        oneLine.SetLineColor(12)
        oneLine.Draw()

        scanContainer.extrema = PhysicsCommands.FindMinimaAndErrors( scanContainer.x, scanContainer.y, returnContainer=True )

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

        # Check number more carefully
        SM_ratio = 0.086
        SMLine = ROOT.TLine( SM_ratio, yMin, SM_ratio, yMax )
        SMLine.SetLineWidth(2)
        SMLine.SetLineColor(9)
        SMLine.Draw()


        xBestfit = scanContainer.x[ scanContainer.extrema.imin ]
        bestfitLine = ROOT.TLine( xBestfit, yMin, xBestfit, yMax )
        bestfitLine.SetLineWidth(2)
        bestfitLine.SetLineColor(2)
        bestfitLine.Draw()


        SaveC( 'totalXSscan' )



    if args.chi2fitToCombination_Yukawa:

        # Workspace to get the parametrization from
        ws = LatestPaths.ws_combined_split_betterYukawa
        # corrMatFile = 'corrMat_Oct17/higgsCombine_CORRMAT_combinedCard_Jul26.MultiDimFit.mH125.root'
        corrMatFile = 'corrMat_Oct17/higgsCombine_CORRMAT_combinedCard_Aug21.MultiDimFit.mH125.root'

        rootFp = ROOT.TFile.Open(ws)
        w = rootFp.Get('w')

        yieldParameterNames = Commands.ListSet( ws, 'yieldParameters' )
        yieldParameters = [ w.var(yPName) for yPName in yieldParameterNames ]
        kappacRealVar = w.var('kappac')
        kappabRealVar = w.var('kappab')

        binBoundaries   = [ w.var(binBoundName).getVal() for binBoundName in Commands.ListSet( ws, 'expBinBoundaries' ) ]
        nBins = len(binBoundaries)-1


        # Get the combination result
        combinationPOIs = Commands.ListPOIs( LatestPaths.ws_combined_unsplit )
        combinationscans = PhysicsCommands.GetScanResults(
            combinationPOIs,
            LatestPaths.scan_ptcombination_combined_profiled,
            pattern = 'combinedCard'
            )
        TgCombination = PhysicsCommands.GetTGraphForSpectrum( combinationPOIs, combinationscans, name='Combination' )
        # binBoundaries, binCenters, binWidths, halfBinWidths, POICenters, POIErrsLeft, POIErrsRight, POIErrsSymm

        # Get relevant bins and errors from scan
        bestfitCenters  = TgCombination.POICenters[:nBins]
        bestfitDownErrs = TgCombination.POIErrsLeft[:nBins]
        bestfitUpErrs   = TgCombination.POIErrsRight[:nBins]
        bestfitSymmErrs = TgCombination.POIErrsSymm[:nBins]


        # Get the correlation matrix
        corrMat = []
        corrMatFp = ROOT.TFile.Open(corrMatFile)
        fit = corrMatFp.Get('fit')
        for poi1 in yieldParameterNames:
            corrMatRow = []
            for poi2 in yieldParameterNames:
                corrMatRow.append( fit.correlation( poi1, poi2 ) )
            corrMat.append( corrMatRow )

        print corrMat
        sys.exit()



        # Build simple chi2 function
        def chi2_function( kappac, kappab ):
            kappacRealVar.setVal(kappac)
            kappabRealVar.setVal(kappab)

            yP_column = numpy.array( [ y.getVal() for y in yieldParameters ] ).T

            # chi2Val = 0.
            # for i in xrange(nBins):
            #     chi2Val += (yieldParameterValues[i]-bestfitCenters[i])**2 / bestfitSymmErrs[i]**2

            chi2 = yP_column.T.dot(  corrMat.dot( yP_column )  )

            return chi2


        # Do a scan

        kappacMin = -35.,
        kappacMax = 35.,
        kappabMin = -13.,
        kappabMax = 13.,

        kappacNPoints = 100
        kappabNPoints = 100


        kappacBinBoundaries = [ kappacMin + i*(kappacMax-kappacMin)/float(kappacNPoints) for i in xrange(kappacNPoints+1) ]
        kappabBinBoundaries = [ kappabMin + i*(kappabMax-kappabMin)/float(kappabNPoints) for i in xrange(kappabNPoints+1) ]

        kappacPoints = [ 0.5*(kappacBinBoundaries[i]+kappacBinBoundaries[i+1]) for i in xrange(kappacNPoints+1) ]
        kappabPoints = [ 0.5*(kappabBinBoundaries[i]+kappabBinBoundaries[i+1]) for i in xrange(kappabNPoints+1) ]


        for kappacVal, kappabVal in zip( kappacPoints, kappabPoints ):

            print cchi2_function()


            break



        rootFp.Close()



    #____________________________________________________________________
    if args.Make2DPlotOfVariableInWS:
        TheoryCommands.SetPlotDir( 'plots_{0}_Yukawa'.format(datestr) )


        ws = LatestPaths.ws_combined_split_betterYukawa_couplingDependentBR

        # varToPlot = 'hggBRmodifier'
        # varToPlot = 'hzzBRmodifier'
        # varToPlot = 'c7_BRscal_hgg'
        # varToPlot = 'c7_BRscal_hzz'
        varToPlot = 'Scaling_hgluglu'


        xVar = 'kappac'
        xRange = [ -5., 5. ]
        xNBins = 200

        yVar = 'kappab'
        yRange = [ -5., 15. ]
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
    if args.PlotOfTotalXSInYukawaWS:
        TheoryCommands.SetPlotDir( 'plots_{0}_Yukawa'.format(datestr) )

        # ws = LatestPaths.ws_combined_split_betterYukawa_couplingDependentBR_profiledTotalXS
        ws = LatestPaths.ws_combined_split_betterYukawa_couplingDependentBR

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
            yMax = 6.3

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




########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'