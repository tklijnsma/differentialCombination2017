#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, itertools, operator, re, argparse, sys
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

    parser.add_argument( '--plotCommands', default=False )
    class CustomAction(argparse.Action):
        def __init__(self, option_strings, dest, **kwargs):
            # Only 1 argument allowed (this is basically the "store_true" action)
            super(CustomAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            # Make sure the parser knows a combineCommand was entered
            setattr( namespace, 'plotCommands', True )
            setattr( namespace, self.dest, True )
            # super(FooAction, self).__call__( parser, namespace, values, option_string=None )
            # setattr( namespace, self.dest, values )

    parser.add_argument( '--plot',                              action=CustomAction )
    parser.add_argument( '--plot2D',                            action=CustomAction )
    parser.add_argument( '--parametrize',                       action=CustomAction )
    parser.add_argument( '--rebinnedTheoryPlot',                action=CustomAction )
    parser.add_argument( '--CheckWSParametrization',            action=CustomAction )
    parser.add_argument( '--plotParametrizationsOnCombination', action=CustomAction )
    parser.add_argument( '--coupling2Dplot',                    action=CustomAction )
    parser.add_argument( '--ReproducePaperPlot',                action=CustomAction )
    parser.add_argument( '--TestParametrizationClassOnDerivedTheoryFiles', action=CustomAction )


########################################
# Methods
########################################    

def main( args ):

    if args.ReproducePaperPlot:

        theoryDir = 'derivedTheoryFiles_YukawaSummed_Aug08/'
        Commands.SetFileFinderDir( theoryDir )

        SMcontainer = TheoryCommands.ReadDerivedTheoryFile(
            Commands.FileFinder( muF=1, muR=1, kappab=1, kappac=1, expectOneFile=True ),
            returnContainer = True,
            verbose = False,
            )
        nBins = len( SMcontainer.binBoundaries ) - 1

        # SMxs_integrated = 0.
        # for iBin in xrange(nBins):
        #     SMxs_integrated += SMcontainer.crosssection[iBin] * ( SMcontainer.binBoundaries[iBin+1] - SMcontainer.binBoundaries[iBin] )

        # Integral to normalize shapes to 1.0
        SMcontainer.integralFn = TheoryCommands.GetIntegral(
            SMcontainer.binBoundaries,
            SMcontainer.crosssection
            )
        SMcontainer.integral = SMcontainer.integralFn( 0., SMcontainer.binBoundaries[-1] )


        yukawaDerivedTheoryFiles = Commands.FileFinder( muF = 1, muR = 1, kappab=1 )
        containers = []
        colorCycle = itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )
        for yukawaDerivedTheoryFile in yukawaDerivedTheoryFiles:
            color = next(colorCycle)

            container = TheoryCommands.ReadDerivedTheoryFile(
                yukawaDerivedTheoryFile,
                returnContainer = True,
                verbose = False,
                )

            container.name = 'kappab_{0}_kappac_{1}'.format(
                Commands.ConvertFloatToStr( container.kappab ),
                Commands.ConvertFloatToStr( container.kappac ),
                )

            container.binCenters = [
                0.5*(container.binBoundaries[i] + container.binBoundaries[i+1]) for i in xrange(nBins)
                ]


            # Integral to normalize shapes to 1.0
            container.integralFn = TheoryCommands.GetIntegral(
                container.binBoundaries,
                container.crosssection
                )
            container.integral = container.integralFn( 0., container.binBoundaries[-1] )


            container.Tg_theory = TheoryCommands.GetTheoryTGraph(
                container.name,
                container.binBoundaries,
                [ (xs/container.integral) / (SMxs/SMcontainer.integral) for xs, SMxs in zip( container.crosssection, SMcontainer.crosssection ) ],
                muBoundLeft   = None,
                muBoundRight  = None,
                boundaries    = True,
                )
            container.Tg_theory.SetLineWidth(2)
            container.Tg_theory.SetLineColor(color)
            container.Tg_theory.SetMarkerColor(color)
            container.Tg_theory.SetLineStyle(1)
            container.Tg_theory.SetMarkerStyle(8)
            container.Tg_theory.SetMarkerSize(0.8)

            containers.append( container )


            # ======================================
            # Make plot

            c.cd()
            c.Clear()
            SetCMargins( RightMargin=0.3 )

            xMinAbs = min([ container.Tg_theory.xMin for container in containers ])
            xMaxAbs = max([ container.Tg_theory.xMax for container in containers ])
            yMinAbs = min([ container.Tg_theory.yMin for container in containers ])
            yMaxAbs = max([ container.Tg_theory.yMax for container in containers ])

            # xMin = xMinAbs - 0.1*( xMaxAbs - xMinAbs )
            # xMax = xMaxAbs + 0.1*( xMaxAbs - xMinAbs )
            xMin = 0.
            xMax = 100.
            yMin = yMinAbs - 0.1*( yMaxAbs - yMinAbs )
            yMax = yMaxAbs + 0.1*( yMaxAbs - yMinAbs )
            # yMin = 0.
            # yMax = 0.14

            base = GetPlotBase(
                xMin = xMin,
                xMax = xMax,
                yMin = yMin,
                yMax = yMax,
                xTitle = 'p_{T} [GeV]', yTitle = 'd#sigma/dp_{T} [pb/GeV]'
                )
            base.Draw('P')

            leg = ROOT.TLegend(
                # 1 - c.GetRightMargin() - 0.3,
                # 1 - c.GetTopMargin() - 0.3,
                # 1 - c.GetRightMargin() ,
                # 1 - c.GetTopMargin() 
                1 - 0.3,
                c.GetBottomMargin(),
                1 - 0.02 ,
                1 - c.GetTopMargin() 

                )
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)

            for container in containers:
                CorrelationMatrices.ConvertTGraphToHistogram(
                    container.Tg_theory,
                    drawImmediately=True,
                    legendObject=leg,
                    verbose=True,
                    )

            leg.Draw()
            outname = 'ReproducePaperPlot'
            SaveC( outname )



    if args.CheckWSParametrization:

        if args.latest:

            # DrawExperimentalBinLines = False
            DrawExperimentalBinLines = True

            wsToCheck = 'workspaces_Aug10/combinedCard_Jul26_CouplingModel_Yukawa_withTheoryUncertainties.root'
            wsParametrization = WSParametrization( wsToCheck )

            theoryDir = 'derivedTheoryFiles_Aug09_YukawaSummed/'
            TheoryFileInterface.SetFileFinderDir( theoryDir )
            yukawaDerivedTheoryFiles = TheoryFileInterface.FileFinder( muF=1, muR=1 )
            # yukawaDerivedTheoryFiles = glob( 'derivedTheoryFiles_Jul25/muR_1_muF_1_Q_1_*.txt' )

            import random
            random.seed(1000)
            yukawaDerivedTheoryFiles = random.sample( yukawaDerivedTheoryFiles, 8 )

            nBins = len( TheoryFileInterface.ReadDerivedTheoryFile( yukawaDerivedTheoryFiles[0] ).binBoundaries ) - 1

            containers = []
            colorCycle = itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )
            for yukawaDerivedTheoryFile in yukawaDerivedTheoryFiles:
                color = next(colorCycle)

                container = TheoryFileInterface.ReadDerivedTheoryFile(
                    yukawaDerivedTheoryFile,
                    returnContainer = True,
                    verbose = True
                    )

                container.name = 'kappab_{0}_kappac_{1}'.format(
                    Commands.ConvertFloatToStr( container.kappab ),
                    Commands.ConvertFloatToStr( container.kappac ),
                    )

                container.binCenters = [
                    0.5*(container.binBoundaries[i] + container.binBoundaries[i+1]) for i in xrange(nBins)
                    ]

                # container.Tg = TheoryCommands.GetTheoryTGraph(
                #     container.name,
                #     container.binBoundaries,
                #     container.ratios,
                #     muBoundLeft   = None,
                #     muBoundRight  = None,
                #     boundaries    = True,
                #     )

                # parametrizationContainer = TheoryCommands.TestParametrizationsInWorkspace(
                #     wsToCheck,
                #     testcouplings = [
                #         { 'kappab' : container.kappab, 'kappac' : container.kappac }
                #         ]
                #     )[0]


                for iBin in xrange(nBins):
                    print 'Bin {0:3} | Value in file = {1:+9.5f},  parametrization = {2:+9.5f}'.format(
                        iBin, container.ratios[iBin], theoryParametrization[iBin]
                        )

                container.Tg_theory = TheoryCommands.GetTheoryTGraph(
                    container.name,
                    container.binBoundaries,
                    container.ratios,
                    muBoundLeft   = None,
                    muBoundRight  = None,
                    boundaries    = True,
                    )
                container.Tg_theory.SetLineWidth(2)
                container.Tg_theory.SetLineColor(color)
                container.Tg_theory.SetMarkerColor(color)
                container.Tg_theory.SetLineStyle(2)
                container.Tg_theory.SetMarkerStyle(8)
                container.Tg_theory.SetMarkerSize(0.8)


                # theoryParametrization = wsParametrization.Evaluate(
                #     kappab = container.kappab , kappac = container.kappac, returnWhat='theory'
                #     )[1:-1]

                # container.Tg_parametrization = TheoryCommands.GetTheoryTGraph(
                #     container.name,
                #     container.binBoundaries,
                #     theoryParametrization,
                #     muBoundLeft   = None,
                #     muBoundRight  = None,
                #     boundaries    = True,
                #     )


                container.Tg_parametrization = wsParametrization.GetOutputContainer( kappab = container.kappab , kappac = container.kappac, returnWhat='theory' ).Tg
                container.Tg_parametrization.SetLineColor(color)
                container.Tg_parametrization.SetMarkerColor(color)
                container.Tg_parametrization.SetLineStyle(1)



                expParametrization = wsParametrization.Evaluate(
                    kappab = container.kappab , kappac = container.kappac, returnWhat='exp'
                    )

                container.Tg_parametrization_expBinning = TheoryCommands.GetTheoryTGraph(
                    container.name,
                    expBinBoundaries,
                    expParametrization,
                    muBoundLeft   = None,
                    muBoundRight  = None,
                    boundaries    = True,
                    )
                container.Tg_parametrization_expBinning.SetLineColor(color)
                container.Tg_parametrization_expBinning.SetMarkerColor(color)
                container.Tg_parametrization_expBinning.SetLineStyle(1)

                containers.append( container )


            # ======================================
            # Additional lines to check

            extraTestCouplings = [
                # { 'kappab' : -6.843, 'kappac' : 21.0 },
                # { 'kappab' : -27.544008255, 'kappac' : 227.531097412 },
                ]

            extraTestContainers = []
            for testCouplings in extraTestCouplings:
                color = next(colorCycle)

                container = TheoryCommands.Container(
                    name = '_'.join([ '{0}_{1:.2f}'.format(key, value) for key, value in sorted(testCouplings.iteritems()) ])
                    )
                container.binBoundaries = containers[0].binBoundaries

                parametrizationContainer = TheoryCommands.TestParametrizationsInWorkspace(
                    wsToCheck,
                    testcouplings = [
                        testCouplings
                        ]
                    )[0]

                container.Tg_parametrization = TheoryCommands.GetTheoryTGraph(
                    container.name,
                    container.binBoundaries,
                    parametrizationContainer.mus_theoryBinning[1:-1],
                    muBoundLeft   = None,
                    muBoundRight  = None,
                    boundaries    = True,
                    )
                container.Tg_parametrization.SetLineColor(color)
                container.Tg_parametrization.SetMarkerColor(color)
                container.Tg_parametrization.SetLineStyle(1)

                # Remove overflow because kappab/kappac don't go very far
                parametrizationContainer.mus_expBinBoundaries = parametrizationContainer.mus_expBinBoundaries[:-1]
                parametrizationContainer.mus_expBinning = parametrizationContainer.mus_expBinning[:-1]

                # Set end of spectrum to whatever was the end for kappab/kappac
                parametrizationContainer.mus_expBinBoundaries[-1] = container.binBoundaries[-1]

                container.Tg_parametrization_expBinning = TheoryCommands.GetTheoryTGraph(
                    container.name,
                    parametrizationContainer.mus_expBinBoundaries,
                    parametrizationContainer.mus_expBinning,
                    muBoundLeft   = None,
                    muBoundRight  = None,
                    boundaries    = True,
                    )
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
                # 1 - c.GetRightMargin() - 0.3,
                # 1 - c.GetTopMargin() - 0.3,
                # 1 - c.GetRightMargin() ,
                # 1 - c.GetTopMargin() 
                1 - 0.3,
                c.GetBottomMargin(),
                1 - 0.02 ,
                1 - c.GetTopMargin() 

                )
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)

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



        else:
            print 'Not implemented'



    if args.plotParametrizationsOnCombination:

        # Combination scan result
        combinationPOIs = Commands.ListPOIs( 'workspaces_May15/combinedCard_May15.root' )
        combinationscans = PhysicsCommands.GetScanResults(
            combinationPOIs,
            'Scan_May15',
            pattern = 'combinedCard'
            )

        # Need theory binBoundaries; awkwardly read this from a derivedTheoryFile
        binBoundaries = TheoryFileInterface.ReadDerivedTheoryFile(
            'derivedTheoryFiles_Aug09_YukawaSummed/YukawaSummed_kappab_1_kappac_1_muR_1_muF_1_Q_1.txt'
            ).binBoundaries

        # Get parametrizations to plot
        couplingsToPlot = [
            { 'kappab' : -0.565910458565, 'kappac' : 9.62666416168 },
            ]


        wsToCheck = 'workspaces_Aug10/combinedCard_Jul26_CouplingModel_Yukawa_withTheoryUncertainties.root'
        wsParametrization = WSParametrization( wsToCheck )
        theoryBinBoundaries = Commands.ReadTheoryBinBoundariesFromWS(wsToCheck)
        expBinBoundaries    = Commands.ReadExpBinBoundariesFromWS(wsToCheck)

        Tgs = []
        colorCycle = itertools.cycle( range(2,9+1) + [ 30, 38, 40, 41, 42 ] + range( 45, 48+1 ) )
        for couplingDict in couplingsToPlot:            
            color = next(colorCycle)

            parametrizationName = '_'.join([ '{0}_{1:.2f}'.format(key, float(value)) for key, value in sorted(couplingDict.iteritems()) ])

            musTheory = wsParametrization.Evaluate(
                kappab=couplingDict['kappab'], kappac=couplingDict['kappac']
                )[1:-1]

            Tg_theory = TheoryCommands.GetTheoryTGraph(
                parametrizationName,
                theoryBinBoundaries,
                musTheory,
                muBoundLeft   = None,
                muBoundRight  = None,
                boundaries    = True,
                )
            Tg_theory.SetLineStyle(1)
            Tg_theory.SetLineColor(color)

            musExp = wsParametrization.Evaluate(
                returnWhat='exp',
                kappab=couplingDict['kappab'], kappac=couplingDict['kappac']
                )

            Tg_exp = TheoryCommands.GetTheoryTGraph(
                parametrizationName + '_exp',
                expBinBoundaries,
                musExp,
                muBoundLeft   = None,
                muBoundRight  = None,
                boundaries    = True,
                )
            Tg_exp.SetLineStyle(1)
            Tg_exp.SetLineColor(color)

            Tgs.append( Tg_theory )
            Tgs.append( Tg_exp )


        # Tgs = []
        # for couplingDict in couplingsToPlot:            

        #     parametrizationContainer = TheoryCommands.TestParametrizationsInWorkspace(
        #         wsToCheck,
        #         testcouplings = [ couplingDict ]
        #         )[0]

        #     parametrizationName = '_'.join([ '{0}_{1}'.format(key, value) for key, value in sorted(couplingDict.iteritems()) ])

        #     Tg_parametrization = TheoryCommands.GetTheoryTGraph(
        #         parametrizationName,
        #         binBoundaries,
        #         parametrizationContainer.mus_theoryBinning[1:-1],
        #         muBoundLeft   = None,
        #         muBoundRight  = None,
        #         boundaries    = True,
        #         )
        #     # Tg_parametrization.SetLineColor(color)
        #     # Tg_parametrization.SetMarkerColor(color)
        #     Tg_parametrization.SetLineStyle(1)

        #     Tgs.append( Tg_parametrization )


        # Draw both
        PhysicsCommands.BasicCombineSpectra(
            ( 'combination', combinationPOIs, combinationscans,
                ( 'SetLineColor', 1 ),
                ( 'SetMarkerStyle', 2 ),
                ( 'SetFillColorAlpha', 1, 0.2 ),
                # ( 'SetFillColor', 13 ),
                # ( 'SetFillStyle', 3544 ),
                # ( 'SetFillStyle', 3345 ),
                ),
            theoryCurves=Tgs,
            drawTheoryCurvesAsLines=True,
            autocolor=False,
            filenameSuffix='_withParametrizationsFromWS'
            )




    if args.plot:

        # res = Commands.ConvertTChainToArray(
        #     'limit',
        #     [
        #         'test001_May03/higgsCombineSCAN_May03_Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.POINTS.0.3.MultiDimFit.mH125.root',
        #         'test001_May03/higgsCombineSCAN_May03_Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.POINTS.12.15.MultiDimFit.mH125.root',
        #         'test001_May03/higgsCombineSCAN_May03_Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.POINTS.16.19.MultiDimFit.mH125.root',
        #         'test001_May03/higgsCombineSCAN_May03_Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.POINTS.4.7.MultiDimFit.mH125.root',
        #         'test001_May03/higgsCombineSCAN_May03_Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.POINTS.8.11.MultiDimFit.mH125.root',
        #         ],
        #     'r_',
        #     )
        # for i in res['r_smH_PTH_30_45']:
        #     print i
        # print
        # print len(res['r_smH_PTH_30_45'])



        hggPOIs = Commands.ListPOIs( 'workspaces_May15/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.root' )
        hggscans = PhysicsCommands.GetScanResults(
            hggPOIs,
            'Scan_May15',
            pattern = 'Datacard_13TeV_differential_pT'
            )
        PhysicsCommands.BasicDrawScanResults( hggPOIs, hggscans, name='hgg' )
        PhysicsCommands.BasicDrawSpectrum( hggPOIs, hggscans, name='hgg' )


        hzzPOIs = Commands.ListPOIs( 'workspaces_May15/usingDavidsCommand_OAshifted.root' )
        hzzscans = PhysicsCommands.GetScanResults(
            hzzPOIs,
            'Scan_May15',
            pattern = 'usingDavidsCommand_OAshifted'
            )
        PhysicsCommands.BasicDrawScanResults( hzzPOIs, hzzscans, name='hzz' )
        PhysicsCommands.BasicDrawSpectrum( hzzPOIs, hzzscans, name='hzz' )


        combinationPOIs = Commands.ListPOIs( 'workspaces_May15/combinedCard_May15.root' )
        combinationscans = PhysicsCommands.GetScanResults(
            combinationPOIs,
            'Scan_May15',
            pattern = 'combinedCard'
            )
        PhysicsCommands.BasicDrawScanResults( combinationPOIs, combinationscans, name='combination' )
        PhysicsCommands.BasicDrawSpectrum( combinationPOIs, combinationscans, name='combination' )


        PhysicsCommands.BasicCombineSpectra(
            ( 'combination', combinationPOIs, combinationscans,
                ( 'SetLineColor', 1 ),
                ( 'SetMarkerStyle', 2 ),
                # ( 'SetFillColorAlpha', 1, 0.2 ),
                ( 'SetFillColor', 13 ),
                # ( 'SetFillStyle', 3544 ),
                ( 'SetFillStyle', 3345 ),
                ),
            ( 'hgg', hggPOIs, hggscans,
                ( 'SetLineColor', 2 ),
                ( 'SetMarkerStyle', 5 ),
                ( 'SetFillColorAlpha', 2, 0.2 ),
                ),
            ( 'hzz', hzzPOIs, hzzscans,
                ( 'SetLineColor', 4 ),
                ( 'SetMarkerStyle', 8 ),
                ( 'SetFillColorAlpha', 4, 0.2 ),
                ),
            printTable=True,
            )


        PhysicsCommands.BasicDrawMultipleParabolas(
            ( 'combination', combinationPOIs, combinationscans,
                ( 'SetLineColor', 1 ),
                # ( 'SetMarkerStyle', 2 ),
                # ( 'SetFillColorAlpha', 1, 0.2 ),
                # ( 'SetFillColor', 13 ),
                # ( 'SetFillStyle', 3544 ),
                # ( 'SetFillStyle', 3345 ),
                ),
            ( 'hgg', hggPOIs, hggscans,
                ( 'SetLineColor', 2 ),
                # ( 'SetMarkerStyle', 5 ),
                # ( 'SetFillColorAlpha', 2, 0.2 ),
                ),
            ( 'hzz', hzzPOIs, hzzscans,
                ( 'SetLineColor', 4 ),
                # ( 'SetMarkerStyle', 8 ),
                # ( 'SetFillColorAlpha', 4, 0.2 ),
                ),
            )


        # Combination only with theory curves
        for pattern in [
            r'^cg_[mp\d]+$',
            r'^ct_[mp\d]+$',
            r'^cb_[mp\d]+$',
            r'ct_[mp\d]+_cb_[mp\d]+_cg_[mp\d]+',
            r'ct_[mp\d]+_cb_([mp\d]+)$', # ct_XX_cb_XX , but not ct_XX_cb_XX_cg_XX
            r'ct_[mp\d]+_cg_[mp\d]+',
            ]:

            theoryTgs = TheoryCommands.LoadTheoryCurves( pattern )

            PhysicsCommands.BasicCombineSpectra(
                ( 'combination', combinationPOIs, combinationscans,
                    ( 'SetLineColor', 1 ),
                    ( 'SetMarkerStyle', 2 ),
                    ( 'SetFillColorAlpha', 1, 0.2 ),
                    # ( 'SetFillColor', 13 ),
                    # ( 'SetFillStyle', 3544 ),
                    # ( 'SetFillStyle', 3345 ),
                    ),
                # ( 'hgg', hggPOIs, hggscans,
                #     ( 'SetLineColor', 2 ),
                #     ( 'SetMarkerStyle', 5 ),
                #     ( 'SetFillColorAlpha', 2, 0.2 ),
                #     ),
                # ( 'hzz', hzzPOIs, hzzscans,
                #     ( 'SetLineColor', 4 ),
                #     ( 'SetMarkerStyle', 8 ),
                #     ( 'SetFillColorAlpha', 4, 0.2 ),
                #     ),
                theoryCurves=theoryTgs
                )


    # ======================================
    # Try parametrization of theory curves

    if args.parametrize:

        # Only to get the pt axis
        pt_variation, ratio_ct_1p1 = TheoryCommands.ReadTheoryFile( 'suppliedInput/numbersFromAgnieszka/HRes_SMEFT_May16/ratio_ctup_new' )

        # Get the following shapes:
        # ct = 0.5, cg = 0.042
        # ct = 0.1, cg = 0.075
        # ct = 2.0, cg = -0.083
        # ct = 1.5, cg = -0.042
        Tgs = TheoryCommands.LoadTheoryCurves( r'ct_[mp\d]+_cg_[mp\d]+' )

        parametrization = TheoryCommands.GetParametrization(
            [ ( Tg.ct**2, Tg.ct*Tg.cg, Tg.cg**2 ) for Tg in Tgs ],
            [ Tg.binValues for Tg in Tgs ],
            # testMode=True
            testMode=False
            )


        # Test it out
        def GetParametrizedTg( ct, cg ):
            yVals = parametrization( ct**2, ct*cg, cg**2 )
            Tg = TheoryCommands.GetTheoryTGraph(
                'ct_{0}_cg_{1}_param'.format( str(ct).replace('.','p').replace('-','m'), str(cg).replace('.','p').replace('-','m') ),
                pt_variation, yVals,
                # prodlist( ratio_SM_NNLO_down, ratio ),
                # prodlist( ratio_SM_NNLO_up, ratio )
                )
            return Tg


        Tgs_copy = Tgs[:]
        Tgs_copy.append( GetParametrizedTg( ct = 1.5, cg = -0.042 ) )
        Tgs_copy.append( GetParametrizedTg( ct = 1.5, cg = -0.043 ) )
        Tgs_copy.append( GetParametrizedTg( ct = 2.0, cg = -0.083 ) )
        TheoryCommands.BasicTheoryPlot( Tgs_copy, drawErrors=False )



        # Scale a theory graph

        experimentalBinning = PhysicsCommands.FigureOutBinning( Commands.ListPOIs( 'workspaces_May15/combinedCard_May15.root' ) )

        ptAxisDummy, SM_NNLO_variationBinned, ratio_SM_NNLO_down_variationBinned, ratio_SM_NNLO_up_variationBinned = \
            TheoryCommands.LoadStandardModelCurves( ptAxis=Tgs[0].binCenters )



        # Make test plot only if other curves are loaded as well
        if args.plot:


            Tgs_expBinning = []
            for Tg in Tgs:

                Tg_expBinning = TheoryCommands.MapPredictionToExperimental(
                    Tg.binBoundaries,
                    map( operator.mul, Tg.binValues, SM_NNLO_variationBinned ),
                    experimentalBinning,
                    makeTGraph = Tg.name + '_rebin',
                    # verbose=True,
                    )

                Tgs_expBinning.append( Tg_expBinning )


            PhysicsCommands.BasicCombineSpectra(
                ( 'combination', combinationPOIs, combinationscans,
                    ( 'SetLineColor', 1 ),
                    ( 'SetMarkerStyle', 2 ),
                    ( 'SetFillColorAlpha', 1, 0.2 ),
                    # ( 'SetFillColor', 13 ),
                    # ( 'SetFillStyle', 3544 ),
                    # ( 'SetFillStyle', 3345 ),
                    ),
                theoryCurves=Tgs_expBinning
                )



    if args.coupling2Dplot:


        if args.latest:

            res = TheoryCommands.PlotCouplingScan2D(
                'workspaces_Aug08/combinedCard_Jul26_CouplingModel.root',
                glob( 'Scan_couplings_Aug09/*.root' ),
                xCoupling = 'kappab',
                yCoupling = 'kappac',
                )

            print '\nBest fit:'
            print res.xCoupling, '=', res.xBestfit
            print res.yCoupling, '=', res.yBestfit


        else:

            res = TheoryCommands.PlotCouplingScan2D(
                'workspaces_Jul27/combinedCard_Jul26_CouplingModel.root',
                glob( 'Scan_couplings_Jul28_1/*.root' ),
                xCoupling = 'kappab',
                yCoupling = 'kappac',
                )

            print '\nBest fit:'
            print res.xCoupling, '=', res.xBestfit
            print res.yCoupling, '=', res.yBestfit

                

    if args.plot2D:

        scan = TheoryCommands.PlotCouplingScan2D(
            # 'workspaces_May30/combinedCard_May15.root',
            # glob( 'Scan_couplings_Jun01_last/*.root' ),
            'workspaces_Jun22/combinedCard_May15_theoryUncertainties.root',
            glob( 'Scan_couplings_Jun22/*.root' ),
            )

        scan = [ s for s in scan if not ( isnan(s[2]) or isinf(s[2]) ) ]
        scan.sort( key = lambda i: i[2] )

        # Get first N unique results
        uniquePoints = []
        allNLLs    = [ s[2] for s in scan ]
        uniqueNLLs = list(set( allNLLs ))
        uniqueNLLs.sort()
        for nll in uniqueNLLs:
            iPoint = allNLLs.index(nll)
            uniquePoints.append( scan[iPoint] )

        testcouplings = []
        for cg, ct, deltaNLL in [ uniquePoints[0], uniquePoints[4], uniquePoints[8] ] :
            testcouplings.append( { 'ct' : ct, 'cg' : cg } )


        if args.plot:


            # ======================================
            # Theory curves from files

            SM = TheoryCommands.ReadDerivedTheoryFile( 'derivedTheoryFiles_May23/SM_NNLO.txt' )
            theories = []
            for derivedTheoryFile in glob( 'derivedTheoryFiles_May23/c*.txt' ):
                theory = TheoryCommands.ReadDerivedTheoryFile(derivedTheoryFile)
                theory['name'] = basename(derivedTheoryFile).replace('.txt','')
                theories.append( theory )

            expBinBoundaries = PhysicsCommands.FigureOutBinning( combinationPOIs )
            rebinnedSMXS = TheoryCommands.Rebin(
                SM['binBoundaries'], SM['crosssection'],
                expBinBoundaries,
                )

            Tgs = []
            colorCycle = itertools.cycle( range(2,5) + range(6,10) )
            for theory in theories:
                color = next(colorCycle)

                rebinnedXS = TheoryCommands.Rebin(
                    theory['binBoundaries'], theory['crosssection'],
                    expBinBoundaries,
                    verbose=False,
                    )
                rebinnedRatios = [ xs / smxs for xs, smxs in zip( rebinnedXS, rebinnedSMXS ) ]

                Tg_theoryBinning = TheoryCommands.GetTheoryTGraph(
                    theory['name'] + '_tb',
                    theory['binBoundaries'],
                    theory['ratios'],
                    boundaries = True,
                    )
                Tg_theoryBinning.SetLineColor( color )

                Tg_expBinning = TheoryCommands.GetTheoryTGraph(
                    theory['name'] + '_eb',
                    expBinBoundaries,
                    rebinnedRatios,
                    boundaries = True,
                    )
                Tg_expBinning.SetLineColor( color )

                # Tgs.append( Tg_theoryBinning )
                # Tgs.append( Tg_expBinning )


            # ======================================
            # Theory curves from parametrization in workspace

            # testcouplings = [
            #     { 'ct' : 3.00, 'cg' : -0.18 },
            #     { 'ct' : 2.00, 'cg' : -0.083 },
            #     { 'ct' : 1.40, 'cg' : -0.05 }
            #     # { 'ct' : 1.75, 'cg' : -0.0625 },
            #     # { 'ct' : 0.35, 'cg' : 0.0585 },
            #     ]

            yValuesPerTestCoupling = TheoryCommands.TestParametrizationsInWorkspace(
                # 'workspaces_May30/combinedCard_May15.root',
                'workspaces_Jun22/combinedCard_May15_theoryUncertainties.root',
                testcouplings
                )

            color = next(colorCycle)
            for couplings, ( y, yParametrization ) in zip( testcouplings, yValuesPerTestCoupling ):

                name = ''
                for key, value in couplings.iteritems():
                    name += '_' + key + '_' + '{0:.2f}'.format(value).replace('.','p').replace('-','m')
                name = name[1:]

                color = next(colorCycle)

                Tg = TheoryCommands.GetTheoryTGraph(
                    name,
                    expBinBoundaries,
                    y,
                    boundaries = True,
                    )
                Tg.SetLineColor( color )
                Tgs.append( Tg )

                yParametrization = yParametrization[1:-1]
                # yParametrization = [ xs / smxs for xs, smxs in zip( yParametrization, SM['crosssection'] ) ]

                Tg = TheoryCommands.GetTheoryTGraph(
                    name + '_param',
                    SM['binBoundaries'],
                    yParametrization,
                    boundaries = True,
                    )
                Tg.SetLineColor( color )
                Tgs.append( Tg )


            # ======================================
            # Plot

            PhysicsCommands.BasicCombineSpectra(
                ( 'combination', combinationPOIs, combinationscans,
                    ( 'SetLineColor', 1 ),
                    ( 'SetMarkerStyle', 2 ),
                    ( 'SetFillColorAlpha', 1, 0.2 ),
                    # ( 'SetFillColor', 13 ),
                    # ( 'SetFillStyle', 3544 ),
                    # ( 'SetFillStyle', 3345 ),
                    ),
                theoryCurves=Tgs[:1],
                drawTheoryCurvesAsLines=True,
                autocolor=False,
                filenameSuffix='_withParametrizationsFromWS'
                )


    if args.rebinnedTheoryPlot:

        container = TheoryCommands.ReadDerivedTheoryFile( 'derivedTheoryFiles_Jun22/SM_NNLO.txt', returnContainer=True )
        TheoryCommands.MakeRebinnedTheoryPlot( container )



    if args.TestParametrizationClassOnDerivedTheoryFiles:

        theoryFileDir = 'derivedTheoryFiles_YukawaQuarkInduced_Aug08'
        # theoryFileDir = 'derivedTheoryFiles_gluonInduced_Aug07/'

        Commands.SetFileFinderDir( theoryFileDir )

        if 'gluonInduced' in theoryFileDir:
            SMfile      = Commands.FileFinder( muF=1, muR=1, kappab=1, kappac=1, expectOneFile=True, Q=1 )
            theoryFiles = [ f for f in Commands.FileFinder( muF=1, muR=1, Q=1 ) ]
        else:
            SMfile      = Commands.FileFinder( muF=1, muR=1, kappab=1, kappac=1, expectOneFile=True )
            theoryFiles = [ f for f in Commands.FileFinder() if not 'mu' in f ]


        SM = TheoryCommands.ReadDerivedTheoryFile(
            SMfile,
            returnContainer = True,
            verbose = False,
            )

        containers = []
        for theoryFile in theoryFiles:
            containers.append(
                TheoryCommands.ReadDerivedTheoryFile(
                    theoryFile,
                    returnContainer = True,
                    verbose = False,
                    )
                )

        colorCycle = itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )
        parametrization = Parametrization()

        # parametrization.Parametrize( containers )
        parametrization.ParametrizeByFitting( containers )


        for container in containers:
            color = next(colorCycle)

            container.name = 'kappab_{0}_kappac_{1}'.format(
                Commands.ConvertFloatToStr( container.kappab ),
                Commands.ConvertFloatToStr( container.kappac )
                )

            container.crosssection_fromParametrization = parametrization.Evaluate(
                kappab = container.kappab, kappac = container.kappac
                )

            container.ratios_fromParametrization = \
                [ i/j if j!=0. else 0. for i, j in zip( container.crosssection_fromParametrization, SM.crosssection ) ]

            container.Tg = TheoryCommands.GetTheoryTGraph(
                container.name,
                container.binBoundaries,
                container.ratios,
                boundaries = True
                )
            container.Tg.SetMarkerStyle(8)
            container.Tg.SetMarkerColor(color)
            container.Tg.SetMarkerSize(0.75)

            container.Tg_fromParametrization = TheoryCommands.GetTheoryTGraph(
                container.name,
                container.binBoundaries,
                container.ratios_fromParametrization,
                boundaries = True
                )
            container.Tg_fromParametrization.SetLineColor(color)
            container.Tg_fromParametrization.SetLineWidth(2)


        # ======================================
        # Make a plot of it

        c.cd()
        c.Clear()
        SetCMargins( RightMargin=0.3 )

        xMinAbs = min([ container.Tg.xMin for container in containers ])
        xMaxAbs = max([ container.Tg.xMax for container in containers ])
        yMinAbs = min([ container.Tg.yMin for container in containers ])
        yMaxAbs = max([ container.Tg.yMax for container in containers ])

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
            # 1 - c.GetRightMargin() - 0.3,
            # 1 - c.GetTopMargin() - 0.3,
            # 1 - c.GetRightMargin() ,
            # 1 - c.GetTopMargin() 
            1 - 0.3,
            c.GetBottomMargin(),
            1 - 0.02 ,
            1 - c.GetTopMargin() 

            )
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)

        for container in containers:
            container.Tg.Draw('XP')
            container.Tg_fromParametrization.Draw('XL')

            # if DrawExperimentalBinLines:
            #     CorrelationMatrices.ConvertTGraphToLinesAndBoxes(
            #         container.Tg_parametrization_expBinning,
            #         drawImmediately=True, legendObject=leg, noBoxes=True )

            leg.AddEntry( container.Tg.GetName(), container.name, 'p' )

        leg.Draw()

        outname = 'ParametrizationClassCheck_{0}'.format( basename(abspath(theoryFileDir)).replace('/','') )
        SaveC( outname )





    





########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'