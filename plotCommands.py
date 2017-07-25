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

    parser.add_argument( '--plot',                            action=CustomAction )
    parser.add_argument( '--plot2D',                          action=CustomAction )
    parser.add_argument( '--parametrize',                     action=CustomAction )
    parser.add_argument( '--rebinnedTheoryPlot',              action=CustomAction )
    parser.add_argument( '--CheckWSParametrization',          action=CustomAction )



########################################
# Methods
########################################    

def main( args ):


    if args.CheckWSParametrization:

        if args.latest:

            wsToCheck = 'workspaces_Jul25/combinedCard_Jul21_CouplingModel.root'


            yukawaDerivedTheoryFiles = glob( 'derivedTheoryFiles_Jul24/muR_1_muF_1_Q_1_*.txt' )
            # yukawaDerivedTheoryFiles = yukawaDerivedTheoryFiles[:10]

            nBins = len( TheoryCommands.ReadDerivedTheoryFile( yukawaDerivedTheoryFiles[0], returnContainer=True ).binBoundaries ) - 1

            containers = []
            colorCycle = itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )
            for yukawaDerivedTheoryFile in yukawaDerivedTheoryFiles:
                color = next(colorCycle)

                container = TheoryCommands.ReadDerivedTheoryFile(
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

                yParametrization_expBinning, yParametrization = TheoryCommands.TestParametrizationsInWorkspace(
                    wsToCheck,
                    testcouplings = [
                        { 'kappab' : container.kappab, 'kappac' : container.kappac }
                        ]
                    )[0]


                for iBin in xrange(nBins):
                    print 'Bin {0:3} | Value in file = {1:+9.5f},  parametrization = {2:+9.5f}'.format(
                        iBin, container.ratios[iBin], yParametrization[iBin+1]
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

                container.Tg_parametrization = TheoryCommands.GetTheoryTGraph(
                    container.name,
                    container.binBoundaries,
                    yParametrization[1:-1],
                    muBoundLeft   = None,
                    muBoundRight  = None,
                    boundaries    = True,
                    )
                container.Tg_parametrization.SetLineColor(color)
                container.Tg_parametrization.SetMarkerColor(color)
                container.Tg_parametrization.SetLineStyle(1)

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

                leg.AddEntry( container.Tg_theory.GetName(), container.name, 'p' )

            leg.Draw()

            outname = '{0}_parametrizationCheck'.format( basename(wsToCheck) )
            SaveC( outname )



        else:
            print 'Not implemented'



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




########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'