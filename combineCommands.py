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
import PlotCommands
from Container import Container

from time import strftime
datestr = strftime( '%b%d' )

import ROOT
from TheoryCommands import c
from TheoryCommands import SaveC
from TheoryCommands import GetPlotBase
from TheoryCommands import SetCMargins
from TheoryCommands import GetUniqueRootName


########################################
# Main
########################################

def AppendParserOptions( parser ):

    parser.add_argument( '--combineCommands', default=False )
    class CustomAction(argparse.Action):
        def __init__(self, option_strings, dest, **kwargs):
            # Only 1 argument allowed (this is basically the "store_true" action)
            super(CustomAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            # Make sure the parser knows a combineCommand was entered
            setattr( namespace, 'combineCommands', True )
            setattr( namespace, self.dest, True )
            # super(FooAction, self).__call__( parser, namespace, values, option_string=None )
            # setattr( namespace, self.dest, values )


    parser.add_argument( '--RenameHggProcesses_Aug21',      action=CustomAction )
    parser.add_argument( '--RenumberHzzProcesses_Aug21',    action=CustomAction )
    parser.add_argument( '--CombineCards_Aug21',            action=CustomAction )
    parser.add_argument( '--CombineCards_Dec15_hbb',        action=CustomAction )

    parser.add_argument( '--RenumberHzzProcesses_Jan24',    action=CustomAction )
    parser.add_argument( '--CombineCards_Jan24_hzz_hbb',    action=CustomAction )

    parser.add_argument( '--RenameHggProcesses_smHcard',    action=CustomAction )
    parser.add_argument( '--RenumberHzzProcesses_smHcard',  action=CustomAction )
    parser.add_argument( '--CombineCards_smHcard',          action=CustomAction )

    parser.add_argument( '--corrMat_PTH',                   action=CustomAction )
    parser.add_argument( '--corrMat_PTH_ggH',               action=CustomAction )
    parser.add_argument( '--corrMat_NJ',                    action=CustomAction )
    parser.add_argument( '--corrMat_YH',                    action=CustomAction )
    parser.add_argument( '--corrMat_PTJ',                   action=CustomAction )
    parser.add_argument( '--plotCorrelationMatrices',       action=CustomAction )

    parser.add_argument( '--t2ws_split',                      action=CustomAction )
    parser.add_argument( '--t2ws_withhbb',                    action=CustomAction )

    parser.add_argument( '--t2ws_combined_unsplit',           action=CustomAction )
    parser.add_argument( '--t2ws_hgg_unsplit',                action=CustomAction )
    parser.add_argument( '--t2ws_hzz_unsplit',                action=CustomAction )
    parser.add_argument( '--bestfit_combined_unsplit',        action=CustomAction )
    parser.add_argument( '--bestfit_hgg_unsplit',             action=CustomAction )
    parser.add_argument( '--bestfit_hzz_unsplit',             action=CustomAction )

    parser.add_argument( '--t2ws_combined_split',           action=CustomAction )
    parser.add_argument( '--t2ws_hgg_split',                action=CustomAction )
    parser.add_argument( '--t2ws_hzz_split',                action=CustomAction )
    parser.add_argument( '--bestfit_combined_split',        action=CustomAction )
    parser.add_argument( '--bestfit_hgg_split',             action=CustomAction )
    parser.add_argument( '--bestfit_hzz_split',             action=CustomAction )

    parser.add_argument( '--scan_combined_unsplit',         action=CustomAction )
    parser.add_argument( '--scan_hgg_unsplit',              action=CustomAction )
    parser.add_argument( '--scan_hzz_unsplit',              action=CustomAction )
    parser.add_argument( '--scan_split_xHfixed',            action=CustomAction )

    parser.add_argument( '--redoPostfit',                   action=CustomAction )
    parser.add_argument( '--corrMat_combined_unsplit',      action=CustomAction )


########################################
# Methods
########################################    

def main( args ):


    #____________________________________________________________________
    # Renaming and combining of split cards

    if args.RenameHggProcesses_Aug21:
        MergeHGGWDatacards.RenameProcesses_Aug21(
            LatestPaths.card_hgg_ggHxH_PTH_unprocessed,
            )

    if args.RenumberHzzProcesses_Aug21:
        MergeHGGWDatacards.RenumberProcessesHZZ_Aug21(
            LatestPaths.card_hzz_ggHxH_PTH_unprocessed,
            )

    if args.CombineCards_Aug21:
        Commands.BasicCombineCards(
            'suppliedInput/combinedCard_ggHxH_{0}.txt'.format(datestr),
            'hgg=' + LatestPaths.card_hgg_ggHxH_PTH,
            'hzz=' + LatestPaths.card_hzz_ggHxH_PTH
            )

    if args.CombineCards_Dec15_hbb:
        Commands.BasicCombineCards(
            'suppliedInput/combinedCard_hgg_hzz_hbb_ggHxH_{0}.txt'.format(datestr),
            'hgg=' + LatestPaths.card_hgg_ggHxH_PTH,
            'hzz=' + LatestPaths.card_hzz_ggHxH_PTH,
            'hbb=' + LatestPaths.card_hbb_ggHxH_PTH
            )


    # Test of new bins for hzz and hbb
    if args.RenumberHzzProcesses_Jan24:
        MergeHGGWDatacards.RenumberProcessesHZZ_Aug21(
            LatestPaths.card_hzz_ggHxH_PTH_newBins_unprocessed,
            )

    if args.CombineCards_Jan24_hzz_hbb:
        Commands.BasicCombineCards(
            'suppliedInput/combinedCard_newBins_hzz_hbb_ggHxH_{0}.txt'.format(datestr),
            # 'hgg=' + LatestPaths.card_hgg_ggHxH_PTH,
            'hzz=' + LatestPaths.card_hzz_ggHxH_PTH_newBins,
            'hbb=' + LatestPaths.card_hbb_ggHxH_PTH
            )

    #____________________________________________________________________
    # Renaming and combining of unsplit cards

    if args.RenameHggProcesses_smHcard:
        MergeHGGWDatacards.RenameProcesses_Aug21(
            LatestPaths.card_hgg_smH_PTH_unprocessed,
            )

    if args.RenumberHzzProcesses_smHcard:
        MergeHGGWDatacards.RenumberProcessesHZZ_Aug21(
            LatestPaths.card_hzz_smH_PTH_unprocessed,
            )

    if args.CombineCards_smHcard:
        Commands.BasicCombineCards(
            'suppliedInput/combinedCard_smH_{0}.txt'.format(datestr),
            'hgg=' + LatestPaths.card_hgg_smH_PTH,
            'hzz=' + LatestPaths.card_hzz_smH_PTH
            )

    #____________________________________________________________________
    if args.t2ws_combined_unsplit:
        Commands.BasicT2WS(
            LatestPaths.card_combined_smH_PTH,
            smartMaps = [
                ( r'.*/smH_PTH_([\d\_GT]+)', r'r_smH_PTH_\1[1.0,-1.0,4.0]' )
                ],
            )

    if args.t2ws_hgg_unsplit:
        Commands.BasicT2WS(
            LatestPaths.card_hgg_smH_PTH,
            smartMaps = [
                ( r'.*/smH_PTH_([\d\_GT]+)', r'r_smH_PTH_\1[1.0,-1.0,4.0]' )
                ],
            )

    if args.t2ws_hzz_unsplit:
        Commands.BasicT2WS(
            LatestPaths.card_hzz_smH_PTH,
            manualMaps = [
                '--PO \'map=.*/smH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_30_45:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_45_85:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_85_125:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_125_200:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_200_350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_GT350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
                ],
            )

    #____________________________________________________________________
    if args.t2ws_split:

        # No real use case for this
        # FIXXH = True
        # if not FIXXH:
        #     Commands.BasicT2WS(
        #         LatestPaths.card_combined_ggHxH_PTH,
        #         smartMaps = [
        #             ( r'.*/.*H_PTH_([\d\_GT]+)', r'r_smH_PTH_\1[1.0,-1.0,4.0]' )
        #             ],
        #         )
        # else:

        datacard = LatestPaths.card_combined_ggHxH_PTH
        if args.hgg: datacard = LatestPaths.card_hgg_ggHxH_PTH
        if args.hzz: datacard = LatestPaths.card_hzz_ggHxH_PTH

        if args.hzz:
            Commands.BasicT2WS(
                datacard,
                manualMaps=[
                    '--PO \'map=.*/ggH_PTH_0_15:r_ggH_PTH_0_15[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_15_30:r_ggH_PTH_15_30[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_30_45:r_ggH_PTH_30_85[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_45_85:r_ggH_PTH_30_85[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_85_125:r_ggH_PTH_85_200[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_125_200:r_ggH_PTH_85_200[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_200_350:r_ggH_PTH_GT200[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_GT350:r_ggH_PTH_GT200[1.0,0.0,3.0]\'',
                    ],
                outputWS = basename(datacard).replace( '.txt', '_xHfixed.root' )
                )

        else:
            Commands.BasicT2WS(
                datacard,
                smartMaps = [
                    ( r'.*/ggH_PTH_([\d\_GT]+)', r'r_ggH_PTH_\1[1.0,-1.0,4.0]' )
                    ],
                outputWS = basename(datacard).replace( '.txt', '_xHfixed.root' )
                )

    #____________________________________________________________________
    if args.t2ws_withhbb:

        if args.hbb:
            datacard = LatestPaths.card_hbb_ggHxH_PTH
            Commands.BasicT2WS(
                datacard,
                manualMaps=[
                    # '--PO \'map=.*/ggH_PTH_200_350:r_ggH_PTH_200_350[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_350_600:r_ggH_PTH_350_600[1.0,0.0,10.0]\'',
                    '--PO \'map=.*/ggH_PTH_GT600:r_ggH_PTH_GT600[1.0,0.0,10.0]\'',
                    ],
                outputWS = basename(datacard).replace( '.txt', '_xHfixed.root' )
                )
        else:
            datacard = LatestPaths.card_combinedWithHbb_ggHxH_PTH
            Commands.BasicT2WS(
                datacard,
                manualMaps=[
                    '--PO \'map=.*/ggH_PTH_GT350:r_ggH_PTH_350_600[1.0,0.0,10.0]\'',
                    '--PO \'map=.*/ggH_PTH_350_600:r_ggH_PTH_350_600[1.0,0.0,10.0]\'',
                    '--PO \'map=.*/ggH_PTH_GT600:r_ggH_PTH_GT600[1.0,0.0,10.0]\'',
                    ],
                smartMaps = [
                    ( r'.*/ggH_PTH_([\d\_GT]+)', r'r_ggH_PTH_\1[1.0,-1.0,4.0]' )
                    ],
                outputWS = basename(datacard).replace( '.txt', '_xHfixed.root' )
                )


    #____________________________________________________________________
    def bestfit( ws ):
        Commands.SetTempJobDir( 'plainWStests_{0}'.format(datestr) )
        Commands.BasicBestfit(
            ws,
            onBatch=False,
            batchJobSubDir = 'job_{0}'.format( basename(ws).replace('/','').replace('.root','') ),
            extraOptions = [
                '-m 125',
                '--floatOtherPOIs=1',
                ]
            )

    if args.bestfit_combined_unsplit:
        bestfit( ws_combined_unsplit )

    if args.bestfit_hgg_unsplit:
        bestfit( ws_onlyhgg_unsplit )
        
    if args.bestfit_hzz_unsplit:
        bestfit( ws_onlyhzz_unsplit )

    if args.bestfit_combined_split:
        bestfit( ws_combined_split )

    if args.bestfit_hgg_split:
        bestfit( ws_onlyhgg_split )
        
    if args.bestfit_hzz_split:
        bestfit( ws_onlyhzz_split )



    #____________________________________________________________________
    if args.corrMat_PTH:
        Commands.ComputeCorrMatrix( LatestPaths.ws_combined_smH, onBatch = False )

    #____________________________________________________________________
    if args.corrMat_PTH_ggH:
        Commands.ComputeCorrMatrix( LatestPaths.ws_combined_ggH_xHfixed, onBatch = False )

    #____________________________________________________________________
    if args.corrMat_NJ:
        Commands.ComputeCorrMatrix( LatestPaths.ws_combined_smH_NJ, onBatch = False )

    #____________________________________________________________________
    if args.corrMat_YH:
        Commands.Warning( 'Using asimov ws for YH (blinded)' )
        Commands.ComputeCorrMatrix( LatestPaths.ws_combined_smH_YH, asimov=True, onBatch = False )

    #____________________________________________________________________
    if args.corrMat_PTJ:
        Commands.Warning( 'Using asimov ws for PTJ (blinded)' )
        Commands.ComputeCorrMatrix( LatestPaths.ws_combined_smH_PTJ, asimov=True, onBatch = False )


    #____________________________________________________________________
    if args.plotCorrelationMatrices:

        pth = Container()
        pth.ws           = LatestPaths.ws_combined_smH
        pth.corrRootFile = LatestPaths.correlationMatrix_PTH
        pth.xTitle       = 'p_{T}^{H} (GeV)'
        PlotCommands.PlotCorrelationMatrix( pth )


        pth_ggH = Container()
        pth_ggH.ws           = LatestPaths.ws_combined_ggH_xHfixed
        pth_ggH.corrRootFile = LatestPaths.correlationMatrix_PTH_ggH
        pth_ggH.xTitle       = 'p_{T}^{H} (GeV) (non-ggH fixed to SM)'
        PlotCommands.PlotCorrelationMatrix( pth_ggH )

        njets = Container()
        njets.ws           = LatestPaths.ws_combined_smH_NJ
        njets.corrRootFile = LatestPaths.correlationMatrix_NJ
        njets.xTitle       = 'N_{jets}'
        PlotCommands.PlotCorrelationMatrix( njets )

        yh = Container()
        yh.ws           = LatestPaths.ws_combined_smH_YH
        yh.corrRootFile = LatestPaths.correlationMatrix_YH
        yh.xTitle       = '|y_{H}|'
        PlotCommands.PlotCorrelationMatrix( yh )

        ptjet = Container()
        ptjet.ws           = LatestPaths.ws_combined_smH_PTJ
        ptjet.corrRootFile = LatestPaths.correlationMatrix_PTJ
        ptjet.xTitle       = 'p_{T}^{j1}'
        PlotCommands.PlotCorrelationMatrix( ptjet )



    #____________________________________________________________________

    # ASIMOV = True
    ASIMOV = False
    if args.asimov: ASIMOV = True

    if args.scan_combined_unsplit:
        Commands.BasicCombineTool(
            LatestPaths.ws_combined_smH,
            POIpattern    = '*',
            nPoints       = 39,
            nPointsPerJob = 3,
            jobDirectory  = 'Scan_PTH_{0}'.format(datestr),
            queue         = 'short.q',
            asimov        = ASIMOV,
            )

    if args.scan_hgg_unsplit:
        Commands.BasicCombineTool(
            LatestPaths.ws_hgg_smH,
            POIpattern    = '*',
            nPoints       = 39,
            nPointsPerJob = 3,
            jobDirectory  = 'Scan_PTH_{0}_hgg'.format(datestr),
            queue         = 'short.q',
            asimov        = ASIMOV,
            )
        
    if args.scan_hzz_unsplit:
        Commands.BasicCombineTool(
            LatestPaths.ws_hzz_smH,
            POIpattern    = '*',
            nPoints       = 39,
            nPointsPerJob = 39,
            jobDirectory  = 'Scan_PTH_{0}_hzz'.format(datestr),
            queue         = 'short.q',
            asimov        = ASIMOV,
            POIRange      = [ 0.0, 3.0 ]
            )

    if args.scan_split_xHfixed:

        ws = LatestPaths.ws_combined_ggH_xHfixed
        if args.hgg: ws = LatestPaths.ws_hgg_ggH_xHfixed
        if args.hzz: ws = LatestPaths.ws_hzz_ggH_xHfixed
        if args.hbb: ws = LatestPaths.ws_hbb_ggH_xHfixed
        if args.combWithHbb: ws = LatestPaths.ws_combWithHbb_ggH_xHfixed
        # ws_combinedWithHbb_ggH_xHfixed

        jobDirectory  = 'out/Scan_PTH_{0}_xHfixed'.format(datestr)
        if args.hgg: jobDirectory += '_hgg'
        if args.hzz: jobDirectory += '_hzz'
        if args.hbb: jobDirectory += '_hbb'
        if args.combWithHbb: jobDirectory += '_combWithHbb'
        # if args.asimov: jobDirectory += '_asimov'

        nPoints = 42
        nPointsPerJob = 3
        if args.hzz: nPointsPerJob = nPoints

        extraOptions = None
        POIRange = None
        if args.hzz:
            POIRange = [ 0.0, 4.0 ]
        if args.hbb:
            nPointsPerJob = nPoints
            POIRange = [ -10.0, 10.0 ]
            extraOptions = [
                '--minimizerStrategy 2',
                '--minimizerTolerance 0.001',
                '--robustFit 1',
                '--minimizerAlgoForMinos Minuit2,Migrad',
                ]
        if args.combWithHbb:
            POIRange = [ 0.0, 8.5 ]
            nPoints = 150
            # nPoints = 2
            nPointsPerJob = 2
            extraOptions = [
                '--minimizerStrategy 2',
                '--minimizerTolerance 0.001',
                '--robustFit 1',
                '--minimizerAlgoForMinos Minuit2,Migrad',
                # '--freezeNuisances qcdeff,r1p0,r2p0,r3p0',
                # '--saveSpecifiedNuis r1p0,r2p0,r3p0,r0p1,r1p1,r2p1,r3p1,qcdeff',
                # '--saveSpecifiedNuis qcdeff',
                # '--saveSpecifiedFunc qcdeff',
                ]
            physicsModelParameterRanges = [
                [ 'qcdeff', 0.001, 8.0 ],
                [ 'r1p0', 0.0, 8.0 ],
                [ 'r2p0', 0.0, 8.0 ],
                [ 'r3p0', 0.0, 8.0 ],
                [ 'r0p1', 0.0, 8.0 ],
                [ 'r1p1', 0.0, 8.0 ],
                [ 'r2p1', 0.0, 8.0 ],
                [ 'r3p1', 0.0, 8.0 ],
                ]
        else:
            physicsModelParameterRanges = []

        Commands.BasicCombineTool(
            ws,
            POIpattern    = '*',
            # POIpattern    = '350_600',
            # POIpattern    = '45_85',
            # POIpattern    = 'GT600',
            POIRange      = POIRange,
            nPoints       = nPoints,
            nPointsPerJob = nPointsPerJob,
            jobDirectory  = jobDirectory,
            queue         = 'short.q',
            asimov        = ASIMOV,
            extraOptions  = extraOptions,
            # disableFloatOtherPOIs = ( True if args.combWithHbb else False ),
            disableFloatOtherPOIs = False,
            physicsModelParameterRanges = physicsModelParameterRanges
            )





########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'