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

    parser.add_argument( '--createDerivedTheoryFiles_YukawaQuarkInduced',   action=CustomAction )
    parser.add_argument( '--createDerivedTheoryFiles_YukawaGluonInduced',   action=CustomAction )
    parser.add_argument( '--mergeGluonInducedWithQuarkInduced',             action=CustomAction )
    parser.add_argument( '--CorrelationMatrices_Yukawa',                    action=CustomAction )
    parser.add_argument( '--couplingT2WS_Yukawa',                           action=CustomAction )
    parser.add_argument( '--couplingBestfit_Yukawa',                        action=CustomAction )
    parser.add_argument( '--coupling2Dplot_Yukawa',                         action=CustomAction )
    parser.add_argument( '--checkWSParametrization_Yukawa',                 action=CustomAction )


########################################
# Methods
########################################    

def main( args ):

    expBinBoundaries = [ 0., 15., 30., 45., 85., 125. ]
    # expBinBoundaries = [ 0., 15., 30., 45., 85., 125., 200., 350. ]
    print 'Hardcoded binBoundaries for Yukawa:'
    print expBinBoundaries
    print ''

    TheoryCommands.SetPlotDir( 'plots_{0}_Yukawa'.format(datestr) )

    LATESTTHEORYDIR_YukawaGluonInduced = 'derivedTheoryFiles_Aug09_YukawaGluonInduced'
    LATESTTHEORYDIR_YukawaQuarkInduced = 'derivedTheoryFiles_Aug09_YukawaQuarkInduced'
    LATESTTHEORYDIR_YukawaSummed       = 'derivedTheoryFiles_Aug09_YukawaSummed'

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
        TheoryFileInterface.MergeGluonAndQuarkInduced(
            gI_theoryFileDir = LATESTTHEORYDIR_YukawaGluonInduced,
            qI_theoryFileDir = LATESTTHEORYDIR_YukawaQuarkInduced,
            verbose = True
            )

    # ======================================
    # Dealing with theory uncertainties

    if args.CorrelationMatrices_Yukawa:

        variationFiles = TheoryFileInterface.FileFinder(
            directory = LATESTTHEORYDIR_YukawaGluonInduced,
            kappab = 1.0, kappac = 1.0
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

        INCLUDE_THEORY_UNCERTAINTIES = True
        # INCLUDE_THEORY_UNCERTAINTIES = False

        datacard = 'suppliedInput/combinedCard_Jul26.txt'

        TheoryFileInterface.SetFileFinderDir( LATESTTHEORYDIR_YukawaSummed )

        extraOptions = [
            '--PO verbose=2',
            '--PO \'higgsMassRange=123,127\'',
            '--PO linearTerms=True',
            ]

        extraOptions.append(
            '--PO SM=[kappab=1,kappac=1,file={0}]'.format(
                TheoryFileInterface.FileFinder( kappab=1, kappac=1, expectOneFile=True )
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
                            TheoryFileInterface.FileFinder( kappab=kappab, kappac=kappac, expectOneFile=True )
                            )
                        )


        suffix = 'Yukawa'
        if INCLUDE_THEORY_UNCERTAINTIES:
            extraOptions.append(
                '--PO correlationMatrix=plots_CorrelationMatrices_Aug09/corrMat_exp.txt' )
            extraOptions.append(
                '--PO theoryUncertainties=plots_CorrelationMatrices_Aug09/errors_for_corrMat_exp.txt' )
            suffix += '_withTheoryUncertainties'


        # Scale these bins with 1.0 regardless of parametrization
        extraOptions.append(
            '--PO skipBins=200_350,GT350'
            )

        random.seed(1002)
        extraOptions.extend( random.sample( possibleTheories, 6 ) )

        Commands.BasicT2WSwithModel(
            datacard,
            'CouplingModel.py',
            suffix = suffix,
            extraOptions = extraOptions,
            )


    if args.couplingBestfit_Yukawa:

        FORCE_FASTSCAN = False
        # FORCE_FASTSCAN = True
        if FORCE_FASTSCAN:
            args.doFastscan = True

        ASIMOV = False
        # ASIMOV = True

        # datacard = 'workspaces_Aug09/combinedCard_Jul26_CouplingModel_Yukawa.root'
        datacard = 'workspaces_Aug09/combinedCard_Jul26_CouplingModel_Yukawa_withTheoryUncertainties.root'

        kappab_ranges = [ -20., 20. ]
        kappac_ranges = [ -50., 100. ]

        jobDirectory = 'Scan_couplings_{0}'.format( datestr )
        if ASIMOV: jobDirectory += '_asimov'
        jobDirectory = Commands.AppendNumberToDirNameUntilItDoesNotExistAnymore( jobDirectory )


        if args.doFastscan:
            nPoints = 6400
            nPointsPerJob = 100
            queue = 'short.q'
        else:
            nPoints = 6400
            nPointsPerJob = 50
            queue = 'long.q'

        Commands.MultiDimCombineTool(
            datacard,
            nPoints       = nPoints,
            nPointsPerJob = nPointsPerJob,
            queue         = queue,
            notOnBatch    = False,
            jobDirectory  = jobDirectory,
            fastscan      = args.doFastscan,
            asimov        = ASIMOV,
            extraOptions  = [
                # '--importanceSampling={0}:couplingScan'.format( abspath('scanTH2D_Jun01.root') ),
                '-P kappab -P kappac',
                '--setPhysicsModelParameters kappab=1.0,kappac=1.0',
                '--setPhysicsModelParameterRanges kappab={0},{1}:kappac={2},{3}'.format(
                    kappab_ranges[0], kappab_ranges[1], kappac_ranges[0], kappac_ranges[1] ),
                '--saveSpecifiedFunc {0}'.format(','.join(
                    Commands.ListSet( datacard, 'yieldParameters' ) + [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
                '--squareDistPoiStep',
                ]
            )





    if args.coupling2Dplot_Yukawa:

        datacard = 'workspaces_Aug09/combinedCard_Jul26_CouplingModel_Yukawa_withTheoryUncertainties.root'
        
        # scandir  = 'Scan_couplings_Aug09_0'

        # Profiled scan:
        scandir  = 'Scan_couplings_Aug09_1'


        # datacard = 'workspaces_Aug09/combinedCard_Jul26_CouplingModel_Yukawa_withTheoryUncertainties.root'
        # scandir  = 'Scan_couplings_Aug09_asimov_0'

        res = TheoryCommands.PlotCouplingScan2D(
            datacard,
            glob( '{0}/*.root'.format(scandir) ),
            xCoupling = 'kappac',
            yCoupling = 'kappab',
            )

        print '\nBest fit:'
        print res.xCoupling, '=', res.xBestfit
        print res.yCoupling, '=', res.yBestfit




    if args.checkWSParametrization_Yukawa:

        plotCombination = False

        # DrawExperimentalBinLines = False
        DrawExperimentalBinLines = True

        wsToCheck = 'workspaces_Aug10/combinedCard_Jul26_CouplingModel_Yukawa_withTheoryUncertainties.root'
        wsParametrization = WSParametrization( wsToCheck )

        theoryDir = LATESTTHEORYDIR_YukawaSummed
        TheoryFileInterface.SetFileFinderDir( theoryDir )
        yukawaDerivedTheoryFiles = TheoryFileInterface.FileFinder( muF=1, muR=1 )


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
            container.Tg_theory.SetLineStyle(2)
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
            container.Tg_parametrization_expBinning.SetLineColor(color)
            container.Tg_parametrization_expBinning.SetMarkerColor(color)
            container.Tg_parametrization_expBinning.SetLineStyle(1)

            containers.append( container )


        # ======================================
        # Additional lines to check

        extraTestCouplings = [
            { 'kappab' : -0.565910458565, 'kappac' : 9.62666416168 },
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