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

    parser.add_argument( '--FitBR_t2ws',              action=CustomAction )
    parser.add_argument( '--FitBR_bestfit',           action=CustomAction )
    parser.add_argument( '--FitBR_scan',              action=CustomAction )
    parser.add_argument( '--FitBR_plot',              action=CustomAction )

    parser.add_argument( '--TotalXS_t2ws',            action=CustomAction )
    parser.add_argument( '--TotalXS_bestfit',         action=CustomAction )
    parser.add_argument( '--TotalXS_scan',            action=CustomAction )
    parser.add_argument( '--TotalXS_plot',            action=CustomAction )



########################################
# Methods
########################################    

def main( args ):

    TheoryCommands.SetPlotDir( 'plots_{0}_BR'.format(datestr) )


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






########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'