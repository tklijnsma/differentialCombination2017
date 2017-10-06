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
    parser.add_argument( '--couplingContourPlot_Top_lumiStudy',             action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Top_BRdependencyComparison',  action=CustomAction )


########################################
# Methods
########################################    

def main( args ):

    # # expBinBoundaries = [ 0., 15., 30., 45., 85., 125. ]
    # expBinBoundaries = [ 0., 15., 30., 45., 85., 125., 200., 350. ]
    # print 'Hardcoded binBoundaries for Top:'
    # print expBinBoundaries
    # print ''

    # LATESTDATACARD_Top            = 'suppliedInput/combinedCard_Jul26.txt'

    # LATESTWORKSPACE_Top           = 'workspaces_Aug11/combinedCard_Jul26_CouplingModel_Top_withTheoryUncertainties.root'

    # LATESTCORRELATIONMATRIX_Top   = 'plots_CorrelationMatrices_Aug11_Top/corrMat_exp.txt'
    # LATESTTHEORYUNCERTAINTIES_Top = 'plots_CorrelationMatrices_Aug11_Top/errors_for_corrMat_exp.txt'

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

        # ======================================
        # Switches

        INCLUDE_THEORY_UNCERTAINTIES = True
        # INCLUDE_THEORY_UNCERTAINTIES = False

        # MAKELUMISCALABLE = True
        MAKELUMISCALABLE = False

        INCLUDE_BR_COUPLING_DEPENDENCY = True
        # INCLUDE_BR_COUPLING_DEPENDENCY = False


        # ======================================

        datacard = LatestPaths.card_combined_split
        if args.hgg:
            datacard = LatestPaths.card_onlyhgg_split_both_renamed
        if args.hzz:
            datacard = LatestPaths.card_onlyhzz_split_renamed

        TheoryFileInterface.SetFileFinderDir( LatestPaths.derivedTheoryFilesDirectory_Top )

        extraOptions = [
            '--PO verbose=2',
            # '--PO verbose=0',
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
                '--PO correlationMatrix={0}'.format(  LatestPaths.correlationMatrix_Top ) )
            extraOptions.append(
                '--PO theoryUncertainties={0}'.format(  LatestPaths.theoryUncertainties_Top ) )
            suffix += '_withTheoryUncertainties'

        if MAKELUMISCALABLE:
            extraOptions.append(
                '--PO lumiScale=True' )
            suffix += '_lumiScale'

        if INCLUDE_BR_COUPLING_DEPENDENCY:
            extraOptions.append( '--PO FitBR=True' )
            suffix += '_couplingDependentBR'


        # Scale these bins with 1.0 regardless of parametrization
        # extraOptions.append(
        #     '--PO skipBins=GT350'
        #     )
        extraOptions.append(
            '--PO binBoundaries=0,15,30,45,85,125,200,350,10000'
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

        # ASIMOV = False
        ASIMOV = True

        # LUMISTUDY = True
        LUMISTUDY = False

        INCLUDE_BR_COUPLING_DEPENDENCY = True
        # INCLUDE_BR_COUPLING_DEPENDENCY = False


        # datacard = LatestPaths.ws_combined_split_top
        datacard = LatestPaths.ws_combined_split_betterTop
        if INCLUDE_BR_COUPLING_DEPENDENCY:
            # datacard = LatestPaths.ws_combined_split_top_couplingDependentBR
            datacard = LatestPaths.ws_combined_split_betterTop_couplingDependentBR

        if LUMISTUDY:
            datacard = LatestPaths.ws_combined_split_top_lumiScalableWS
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
            # nPoints = 6400
            # nPointsPerJob = 8
            nPoints = 4900
            nPointsPerJob = 14
            queue = 'all.q'
            # if ASIMOV:
            #     # This does not converge in time
            #     queue = 'short.q'
            #     nPoints = 1600
            #     nPointsPerJob = 2
            if args.hzz:
                nPointsPerJob = 80
                queue = 'short.q'


        extraOptions = [
            '-P ct -P cg' + ( ' -P kappa_V' if INCLUDE_BR_COUPLING_DEPENDENCY else '' ),
            '--squareDistPoiStep',
            # 
            ( '--setPhysicsModelParameters ct=1.0,cg=0.0'
                + ( ',lumiScale=8.356546' if LUMISTUDY else '' )
                + ( ',kappa_V=1.0' if INCLUDE_BR_COUPLING_DEPENDENCY else '' )
                ),
            # 
            ( '--setPhysicsModelParameterRanges ct={0},{1}:cg={2},{3}'.format( ct_ranges[0], ct_ranges[1], cg_ranges[0], cg_ranges[1] )
                + ( ':kappa_V=-10.0,10.0' if INCLUDE_BR_COUPLING_DEPENDENCY else '' )
                ),
            ]

        variablesToSave = []
        variablesToSave.extend( Commands.ListSet( datacard, 'yieldParameters' ) )
        variablesToSave.extend( [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ] )
        if INCLUDE_BR_COUPLING_DEPENDENCY:
            variablesToSave.extend( Commands.ListSet( datacard, 'hgg_yieldParameters' ) )
            variablesToSave.extend( Commands.ListSet( datacard, 'hzz_yieldParameters' ) )
            variablesToSave.extend( Commands.ListSet( datacard, 'BRvariables' ) )

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


    #____________________________________________________________________
    if args.coupling2Dplot_Top:

        datacard = LatestPaths.ws_combined_split_top
        scandir    = LatestPaths.scan_top_combined_profiled
        if args.hgg:
            datacard = LatestPaths.ws_onlyhgg_split_top
            scandir = LatestPaths.scan_top_hgg_profiled
        if args.hzz:
            datacard = LatestPaths.ws_onlyhzz_split_top
            scandir = LatestPaths.scan_top_hzz_profiled

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
    if (
        args.couplingContourPlot_Top
        or args.couplingContourPlot_Top_lumiStudy
        or args.couplingContourPlot_Top_BRdependencyComparison
        ):

        xCoupling = 'ct'
        yCoupling = 'cg'
        titles = {
            'ct'          : '#kappa_{t}',
            'cg'          : '#kappa_{g}',
            'hgg'         : 'H #rightarrow #gamma#gamma',
            'hzz'         : 'H #rightarrow 4l',
            'combined'    : 'Combination',
            # For lumi study:
            'regularLumi' : 'Expected at 35.9 fb^{-1}',
            'highLumi'    : 'Expected at 300 fb^{-1}',
            }


    #____________________________________________________________________
    if args.couplingContourPlot_Top:


        combined_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_top_combined_profiled ) )
        hgg_rootfiles      = glob( '{0}/*.root'.format(LatestPaths.scan_top_hgg_profiled) )
        hzz_rootfiles      = glob( '{0}/*.root'.format(LatestPaths.scan_top_hzz_profiled) )

        hgg = TheoryCommands.GetTH2FromListOfRootFiles(
            hgg_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        hgg.color = 2
        hgg.name = 'hgg'

        hzz = TheoryCommands.GetTH2FromListOfRootFiles(
            hzz_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        hzz.color = 4
        hzz.name = 'hzz'

        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
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
        for container in containers: container.title = titles.get( container.name, container.name )

        TheoryCommands.BasicMixedContourPlot(
            containers,
            xMin = -0.2,
            xMax = 2.0,
            yMin = -0.08,
            yMax = 0.135,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours',
            x_SM      = 1.,
            y_SM      = 0.,
            )


    #____________________________________________________________________
    if args.couplingContourPlot_Top_lumiStudy:

        ASIMOV = True

        if ASIMOV:
            combined_rootfiles      = glob( '{0}/*.root'.format( LatestPaths.scan_top_combined_profiled_asimov ) )
            combined_lum8_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_top_combined_profiled_asimov_lum8 ) )
        else:
            print 'ERROR: Only implemented for Asimov now'
            sys.exit()

        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name = 'regularLumi'

        combined_lum8 = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_lum8_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined_lum8.color = 2
        combined_lum8.name = 'highLumi'


        containers = [
            combined,
            combined_lum8
            ]

        for container in containers: container.title = titles.get( container.name, container.name )

        TheoryCommands.BasicMixedContourPlot(
            containers,
            xMin = -0.2,
            xMax = 2.0,
            yMin = -0.08,
            yMax = 0.135,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_LumiStudy',
            x_SM      = 1.,
            y_SM      = 0.,
            )


    #____________________________________________________________________
    if args.couplingContourPlot_Top_BRdependencyComparison:

        # combined_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_top_combined_profiled_asimov ) )
        # scalingBR_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_top_combined_profiled_asimov_couplingDependentBR ) )

        combined_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_betterTop_combined_asimov ) )
        scalingBR_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_betterTop_combined_asimov_couplingDependentBR ) )


        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name = 'regularBR'
        combined.title = 'BR constant'

        scalingBR = TheoryCommands.GetTH2FromListOfRootFiles(
            scalingBR_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        scalingBR.color = 2
        scalingBR.name = 'scalingBR'
        scalingBR.title = 'BR(#kappa_{t})'


        TheoryCommands.BasicMixedContourPlot(
            [ combined, scalingBR ],
            xMin = -0.2,
            xMax = 2.0,
            yMin = -0.08,
            yMax = 0.135,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_BRcouplingDependency',
            x_SM      = 1.,
            y_SM      = 0.,
            plotIndividualH2s = True,
            )




    #____________________________________________________________________
    if args.checkWSParametrization_Top:

        # plotCombination = True
        plotCombination = False

        DrawExperimentalBinLines = False
        # DrawExperimentalBinLines = True

        # wsToCheck = 'workspaces_Aug11/combinedCard_Jul26_CouplingModel_Top_withTheoryUncertainties.root'
        # wsToCheck = LATESTWORKSPACE_Top
        wsToCheck = LatestPaths.ws_combined_split_top
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
            container.Tg_theory.SetLineStyle(1)
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
        yMaxAbs = 1.4*max([ container.Tg_theory.yMax for container in containers ])

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


        if doBigLegend:
            leg = ROOT.TLegend(
                1 - 0.3,
                c.GetBottomMargin(),
                1 - 0.02 ,
                1 - c.GetTopMargin() 
                )
        else:
            leg = ROOT.TLegend(
                c.GetLeftMargin(),
                1 - c.GetTopMargin() - 0.25,
                1 - c.GetRightMargin(),
                1 - c.GetTopMargin() 
                )
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetNColumns(3)


        Tg_dotDummy = ROOT.TGraph( 1, array( 'f', [-999.] ), array( 'f', [-999.] ) )
        Tg_dotDummy.SetMarkerStyle(8)
        Tg_dotDummy.SetMarkerSize(0.8)
        Tg_dotDummy.SetMarkerColor(49)
        Tg_dotDummy.SetName( 'dotDummy' )
        ROOT.SetOwnership( Tg_dotDummy, False )
        Tg_dotDummy.Draw('P')

        Tg_lineDummy = ROOT.TGraph( 1, array( 'f', [-999.] ), array( 'f', [-999.] ) )
        Tg_lineDummy.SetLineWidth(2)
        Tg_lineDummy.SetLineColor(49)
        Tg_lineDummy.SetName( 'lineDummy' )
        ROOT.SetOwnership( Tg_lineDummy, False )
        Tg_lineDummy.Draw('L')

        leg.AddEntry( Tg_dotDummy.GetName(),  'Theory calc.', 'P' )
        leg.AddEntry( Tg_lineDummy.GetName(), 'Parametrization', 'L' )


        if plotCombination:
            # Combination scan result
            combinationPOIs = Commands.ListPOIs( LatestPaths.ws_combined_unsplit )
            combinationscans = PhysicsCommands.GetScanResults(
                combinationPOIs,
                'Scan_May15',
                pattern = 'combinedCard'
                )
            TgCombination = PhysicsCommands.GetTGraphForSpectrum( combinationPOIs, combinationscans, name='Combination' )
            TgCombination.SetLineColor(1)
            TgCombination.SetFillColorAlpha( 1, 0.2 )
            TgCombination.SetName('Combination')

            CorrelationMatrices.ConvertTGraphToLinesAndBoxes(
                TgCombination,
                drawImmediately=True, legendObject=leg, noBoxes=False, xMaxExternal=xMax, yMinExternal=yMin )


        for container in containers:
            container.Tg_theory.Draw('XP')
            container.Tg_parametrization.Draw('XL')

            if DrawExperimentalBinLines:
                CorrelationMatrices.ConvertTGraphToLinesAndBoxes(
                    container.Tg_parametrization_expBinning,
                    drawImmediately=True,
                    # legendObject=leg,
                    noBoxes=True,
                    xMaxExternal=xMax
                    )


            if doBigLegend:
                leg.AddEntry( container.Tg_theory.GetName(), container.name, 'p' )
            else:
                leg.AddEntry(
                    container.Tg_theory.GetName(),
                    '#kappa_{{t}} = {0:.2f}, #kappa_{{g}} = {1:.2f}'.format( container.ct, container.cg ),
                    'lp'
                    )


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