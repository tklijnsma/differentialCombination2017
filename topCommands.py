#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, itertools, operator, re, argparse, sys, random
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

    parser.add_argument( '--topCommands', default=False )
    class CustomAction(argparse.Action):
        def __init__(self, option_strings, dest, **kwargs):
            # Only 1 argument allowed (this is basically the "store_true" action)
            super(CustomAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            # Make sure the parser knows a combineCommand was entered
            setattr( namespace, 'topCommands', True )
            setattr( namespace, self.dest, True )
            # super(FooAction, self).__call__( parser, namespace, values, option_string=None )
            # setattr( namespace, self.dest, values )

    parser.add_argument( '--createDerivedTheoryFiles_Top',                  action=CustomAction )
    parser.add_argument( '--CorrelationMatrices_Top',                       action=CustomAction )
    parser.add_argument( '--couplingT2WS_Top',                              action=CustomAction )
    parser.add_argument( '--couplingBestfit_Top',                           action=CustomAction )
    parser.add_argument( '--coupling2Dplot_Top',                            action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Top',                       action=CustomAction )
    parser.add_argument( '--checkWSParametrization_Top',                    action=CustomAction )



########################################
# Methods
########################################    

def main( args ):

    # expBinBoundaries = [ 0., 15., 30., 45., 85., 125. ]
    expBinBoundaries = [ 0., 15., 30., 45., 85., 125., 200., 350. ]
    print 'Hardcoded binBoundaries for Top:'
    print expBinBoundaries
    print ''

    LATESTDATACARD_Top            = 'suppliedInput/combinedCard_Jul26.txt'

    LATESTWORKSPACE_Top           = 'workspaces_Aug11/combinedCard_Jul26_CouplingModel_Top_withTheoryUncertainties.root'

    LATESTCORRELATIONMATRIX_Top   = 'plots_CorrelationMatrices_Aug11_Top/corrMat_exp.txt'
    LATESTTHEORYUNCERTAINTIES_Top = 'plots_CorrelationMatrices_Aug11_Top/errors_for_corrMat_exp.txt'

    TheoryCommands.SetPlotDir( 'plots_{0}_Top'.format(datestr) )


    #____________________________________________________________________
    if args.createDerivedTheoryFiles_Top:
        TheoryFileInterface.CreateDerivedTheoryFiles_Top(
            # verbose = True,
            )


    #____________________________________________________________________
    if args.CorrelationMatrices_Top:

        CorrelationMatrices.SetPlotDir( 'plots_CorrelationMatrices_{0}_Top'.format(datestr) )

        variationFiles = TheoryFileInterface.FileFinder(
            directory = LatestPaths.derivedTheoryFilesDirectory_Top,
            ct = 1.0, cg = 1.0, cb = 1.0
            )

        variations = [
            TheoryFileInterface.ReadDerivedTheoryFile( variationFile ) for variationFile in variationFiles ]

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

    #____________________________________________________________________
    if args.couplingT2WS_Top:

        INCLUDE_THEORY_UNCERTAINTIES = True
        # INCLUDE_THEORY_UNCERTAINTIES = False

        datacard = LatestPaths.card_combined_split
        if args.hgg:
            datacard = LatestPaths.card_onlyhgg_split_both_renamed
        if args.hzz:
            datacard = LatestPaths.card_onlyhzz_split_renamed

        TheoryFileInterface.SetFileFinderDir( LatestPaths.derivedTheoryFilesDirectory_Top )

        extraOptions = [
            '--PO verbose=2',
            '--PO \'higgsMassRange=123,127\'',
            '--PO linearTerms=False',
            '--PO splitggH=True',
            ]

        extraOptions.append(
            '--PO SM=[ct=1,cg=0,file={0}]'.format(
                TheoryFileInterface.FileFinder( ct=1, cg=1, cb=1, muR=1, muF=1, Q=1, expectOneFile=True )
                )
            )

        
        theoryFiles = TheoryFileInterface.FileFinder(
            ct='*', cg='*', muR=1, muF=1, Q=1, filter='cb'
            )
        possibleTheories = []
        for theoryFile in theoryFiles:
            ct = Commands.ConvertStrToFloat( re.search( r'ct_([\dmp]+)', theoryFile ).group(1) )
            cg = Commands.ConvertStrToFloat( re.search( r'cg_([\dmp]+)', theoryFile ).group(1) )
            possibleTheories.append(
                '--PO theory=[ct={0},cg={1},file={2}]'.format(
                    ct, cg, theoryFile
                    )                
                )
        extraOptions.extend(possibleTheories)


        suffix = 'Top'
        if INCLUDE_THEORY_UNCERTAINTIES:
            extraOptions.append(
                '--PO correlationMatrix={0}'.format(LATESTCORRELATIONMATRIX_Top) )
            extraOptions.append(
                '--PO theoryUncertainties={0}'.format(LATESTTHEORYUNCERTAINTIES_Top) )
            suffix += '_withTheoryUncertainties'


        # Scale these bins with 1.0 regardless of parametrization
        extraOptions.append(
            '--PO skipBins=GT350'
            )

        Commands.BasicT2WSwithModel(
            datacard,
            'CouplingModel.py',
            suffix = suffix,
            extraOptions = extraOptions,
            )


    #____________________________________________________________________
    if args.couplingBestfit_Top:

        doFastscan = True
        if args.notFastscan: doFastscan = False

        ASIMOV = False
        # ASIMOV = True

        datacard = LatestPaths.ws_combined_split_top
        if args.hgg:
            datacard = LatestPaths.ws_onlyhgg_split_top
        if args.hzz:
            datacard = LatestPaths.ws_onlyhzz_split_top


        ct_ranges = [ -1., 2. ]
        cg_ranges = [ -0.1, 0.2 ]


        jobDirectory = 'Scan_Top_{0}'.format( datestr )
        if ASIMOV: jobDirectory += '_asimov'
        jobDirectory = Commands.AppendNumberToDirNameUntilItDoesNotExistAnymore( jobDirectory )


        if doFastscan:
            nPoints = 12800
            nPointsPerJob = 800
            queue = 'short.q'
        else:
            nPoints = 6400
            nPointsPerJob = 8
            queue = 'all.q'
            if args.hzz:
                nPointsPerJob = 80
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
            jobPriority   = -2,
            extraOptions  = [
                # '--importanceSampling={0}:couplingScan'.format( abspath('scanTH2D_Jun01.root') ),
                '-P ct -P cg',
                '--setPhysicsModelParameters ct=1.0,cg=0.0',
                '--setPhysicsModelParameterRanges ct={0},{1}:cg={2},{3}'.format(
                    ct_ranges[0], ct_ranges[1], cg_ranges[0], cg_ranges[1] ),
                '--saveSpecifiedFunc {0}'.format(','.join(
                    Commands.ListSet( datacard, 'yieldParameters' ) + [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
                '--squareDistPoiStep',
                ]
            )


    #____________________________________________________________________
    if args.coupling2Dplot_Top:

        datacard = LatestPaths.ws_combined_split_top
        
        # scandir  = 'Scan_couplings_Aug11'
        # scandir  = 'Scan_couplings_Aug11_0'
        # scandir  = 'Scan_couplings_Aug11_1'
        # scandir  = 'Scan_Top_Aug21'
        scandir    = LatestPaths.scan_top_combined_profiled

        res = TheoryCommands.PlotCouplingScan2D(
            datacard,
            glob( '{0}/*.root'.format(scandir) ),
            xCoupling = 'ct',
            yCoupling = 'cg',
            SM = ( 1.0, 0.0 ),
            )

        print '\nBest fit:'
        print res.xCoupling, '=', res.xBestfit
        print res.yCoupling, '=', res.yBestfit


    #____________________________________________________________________
    if args.couplingContourPlot_Top:

        xCoupling = 'ct'
        yCoupling = 'cg'
        titles = {
            'ct': '#kappa_{t}', 'cg' : '#kappa_{g}',
            'hgg' : 'H #rightarrow #gamma#gamma',
            'hzz' : 'H #rightarrow 4l',
            'combined' : 'Combination',
            }

        hgg = TheoryCommands.GetTH2FromListOfRootFiles(
            glob( '{0}/*.root'.format(LatestPaths.scan_top_hgg_profiled) ),
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        hgg.color = 2
        hgg.name = 'hgg'

        hzz = TheoryCommands.GetTH2FromListOfRootFiles(
            glob( '{0}/*.root'.format(LatestPaths.scan_top_hzz_profiled) ),
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        hzz.color = 4
        hzz.name = 'hzz'

        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            glob( '{0}/*.root'.format( LatestPaths.scan_top_combined_profiled ) ),
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name = 'combined'

        containers = [
            hgg,
            hzz,
            combined
            ]


        for container in containers:
            container.contours_1sigma = TheoryCommands.GetContoursFromTH2( container.H2, 2.30 )
            container.contours_2sigma = TheoryCommands.GetContoursFromTH2( container.H2, 6.18 )

        c.cd()
        c.Clear()
        SetCMargins()

        xMin = -0.2
        xMax = 1.35
        yMin = -0.03
        yMax = 0.135

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
                Tg.SetLineWidth(2)
                Tg.SetLineColor( container.color )
                Tg.SetLineStyle(1)
                Tg.Draw('CSAME')
                if Tg == container.contours_1sigma[0]:
                    Tg.SetName( '{0}_contour_1sigma'.format(container.name) )
                    leg.AddEntry( Tg.GetName(), titles.get( container.name, container.name ), 'l' )

            for Tg in container.contours_2sigma:
                Tg.SetLineWidth(2)
                Tg.SetLineColor( container.color )
                Tg.SetLineStyle(2)
                Tg.Draw('CSAME')

            Tpoint = ROOT.TGraph( 1, array( 'd', [container.xBestfit] ), array( 'd', [container.yBestfit] ) )
            ROOT.SetOwnership( Tpoint, False )
            Tpoint.SetMarkerSize(2)
            Tpoint.SetMarkerStyle(34)
            Tpoint.SetMarkerColor( container.color )
            Tpoint.Draw('PSAME')
            Tpoint.SetName( '{0}_bestfitpoint'.format( container.name ) )

        TpointSM = ROOT.TGraph( 1, array( 'd', [1.0] ), array( 'd', [0.0] ) )
        ROOT.SetOwnership( TpointSM, False )
        TpointSM.SetMarkerSize(2)
        TpointSM.SetMarkerStyle(21)
        TpointSM.SetMarkerColor( 12 )
        TpointSM.Draw('PSAME')

        leg.Draw()

        SaveC( 'contours' )


    #____________________________________________________________________
    if args.checkWSParametrization_Top:

        plotCombination = False

        # DrawExperimentalBinLines = False
        DrawExperimentalBinLines = True

        # wsToCheck = 'workspaces_Aug11/combinedCard_Jul26_CouplingModel_Top_withTheoryUncertainties.root'
        wsToCheck = LATESTWORKSPACE_Top
        wsParametrization = WSParametrization( wsToCheck )

        theoryDir = LatestPaths.derivedTheoryFilesDirectory_Top
        TheoryFileInterface.SetFileFinderDir( theoryDir )
        topDerivedTheoryFiles = TheoryFileInterface.FileFinder(
            ct='*', cg='*', muR=1, muF=1, Q=1, filter='cb'
            )

        containers = []
        colorCycle = itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )
        for topDerivedTheoryFile in topDerivedTheoryFiles:
            color = next(colorCycle)

            print topDerivedTheoryFile

            container = TheoryFileInterface.ReadDerivedTheoryFile( topDerivedTheoryFile )
            container.name = 'ct_{0}_cg_{1}'.format(
                Commands.ConvertFloatToStr( container.ct ),
                Commands.ConvertFloatToStr( container.cg ),
                )

            container.Tg_theory = TheoryFileInterface.ReadDerivedTheoryFileToTGraph( topDerivedTheoryFile, name = container.name )
            container.Tg_theory.SetLineWidth(2)
            container.Tg_theory.SetLineColor(color)
            container.Tg_theory.SetMarkerColor(color)
            container.Tg_theory.SetLineStyle(2)
            container.Tg_theory.SetMarkerStyle(8)
            container.Tg_theory.SetMarkerSize(0.8)

            container.Tg_parametrization = wsParametrization.GetOutputContainer(
                ct = container.ct , cg = container.cg,
                returnWhat='theory' ).Tg
            container.Tg_parametrization.SetLineColor(color)
            container.Tg_parametrization.SetMarkerColor(color)
            container.Tg_parametrization.SetLineStyle(1)

            container.Tg_parametrization_expBinning = wsParametrization.GetOutputContainer(
                ct = container.ct , cg = container.cg,
                returnWhat='exp' ).Tg
            container.Tg_parametrization_expBinning.SetLineColor(color)
            container.Tg_parametrization_expBinning.SetMarkerColor(color)
            container.Tg_parametrization_expBinning.SetLineStyle(1)

            containers.append( container )


        # ======================================
        # Additional lines to check

        extraTestCouplings = [
            # { 'ct' : -0.565910458565, 'cg' : 9.62666416168 },
            ]

        extraTestContainers = []
        for testCouplings in extraTestCouplings:
            color = next(colorCycle)

            container = Container(
                name = '_'.join([ '{0}_{1:.2f}'.format(key, value) for key, value in sorted(testCouplings.iteritems()) ])
                )
            container.binBoundaries = containers[0].binBoundaries

            container.Tg_parametrization = wsParametrization.GetOutputContainer(
                ct = testCouplings['ct'] , cg = testCouplings['cg'],
                returnWhat='theory' ).Tg
            container.Tg_parametrization.SetLineColor(color)
            container.Tg_parametrization.SetMarkerColor(color)
            container.Tg_parametrization.SetLineStyle(1)

            container.Tg_parametrization_expBinning = wsParametrization.GetOutputContainer(
                ct = testCouplings['ct'] , cg = testCouplings['cg'],
                returnWhat='exp' ).Tg
            container.Tg_parametrization_expBinning.SetLineColor(color)
            container.Tg_parametrization_expBinning.SetMarkerColor(color)
            container.Tg_parametrization_expBinning.SetLineStyle(1)

            extraTestContainers.append( container )



        # ======================================
        # Make plot

        c.cd()
        c.Clear()
        SetCMargins( RightMargin=0.3 )

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
            xTitle = 'p_{T} [GeV]', yTitle = '#mu'
            )
        base.Draw('P')

        leg = ROOT.TLegend(
            1 - 0.3,
            c.GetBottomMargin(),
            1 - 0.02 ,
            1 - c.GetTopMargin() 
            )
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)


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
                    drawImmediately=True, legendObject=leg, noBoxes=True, xMaxExternal=xMax )

            leg.AddEntry( container.Tg_theory.GetName(), container.name, 'p' )


        for container in extraTestContainers:
            container.Tg_parametrization.Draw('XL')
            if DrawExperimentalBinLines:
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