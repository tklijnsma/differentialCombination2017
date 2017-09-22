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

from time import strftime
datestr = strftime( '%b%d' )


########################################
# Main
########################################

def AppendParserOptions( parser ):

    parser.add_argument( '--njetsCommands', default=False )
    class CustomAction(argparse.Action):
        def __init__(self, option_strings, dest, **kwargs):
            # Only 1 argument allowed (this is basically the "store_true" action)
            super(CustomAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            # Make sure the parser knows a combineCommand was entered
            setattr( namespace, 'njetsCommands', True )
            setattr( namespace, self.dest, True )
            # super(FooAction, self).__call__( parser, namespace, values, option_string=None )
            # setattr( namespace, self.dest, values )

    parser.add_argument( '--rename_hgg_njets',                       action=CustomAction )
    parser.add_argument( '--njets_combineCards',                     action=CustomAction )
    parser.add_argument( '--njets_t2ws_unsplit',                     action=CustomAction )
    parser.add_argument( '--njets_bestfit_unsplit',                  action=CustomAction )
    parser.add_argument( '--plot_nJetsSpectra',                  action=CustomAction )


########################################
# Methods
########################################    

def main( args ):

    TheoryCommands.SetPlotDir( 'plots_{0}_nJets'.format(datestr) )


    #____________________________________________________________________
    if args.rename_hgg_njets:
        MergeHGGWDatacards.RenameProcesses_Hgg_nJets(
            LatestPaths.card_njets_hgg_unsplit_unrenamed,
            )

    #____________________________________________________________________
    if args.njets_combineCards:

        Commands.BasicCombineCards(
            'suppliedInput/combinedCard_nJets_{0}.txt'.format(datestr),
            LatestPaths.card_njets_hgg_unsplit_renamed,
            LatestPaths.card_njets_hzz_unsplit
            )


    #____________________________________________________________________
    if args.njets_t2ws_unsplit:

        datacard = LatestPaths.card_njets_combined_unsplit
        if args.hgg:
            datacard = LatestPaths.card_njets_hgg_unsplit_renamed
        if args.hzz:
            datacard = LatestPaths.card_njets_hzz_unsplit

        if args.hzz:
            Commands.BasicT2WS(
                datacard,
                manualMaps=[
                    '--PO \'map=.*/smH_NJ_0:r_smH_NJ_0[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/smH_NJ_1:r_smH_NJ_1[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/smH_NJ_2:r_smH_NJ_2[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/smH_NJ_3:r_smH_NJ_GE3[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/smH_NJ_GE4:r_smH_NJ_GE3[1.0,-1.0,4.0]\'',
                    ],
                )
        else:
            Commands.BasicT2WS(
                datacard,
                smartMaps = [
                    ( r'.*/smH_NJ_([\d\_GE]+)', r'r_smH_NJ_\1[1.0,-1.0,4.0]' )
                    ],
                )

    #____________________________________________________________________
    if args.njets_bestfit_unsplit:

        ASIMOV = False

        nPoints       = 39
        nPointsPerJob = 3


        ws = LatestPaths.ws_njets_combined_unsplit
        if args.hgg:
            ws = LatestPaths.ws_njets_hgg_unsplit
        if args.hzz:
            ws = LatestPaths.ws_njets_hzz_unsplit
            nPoints       = 39
            nPointsPerJob = 39

        Commands.BasicCombineTool(
            ws,
            POIpattern    = '*',
            nPoints       = 39,
            nPointsPerJob = 3,
            jobDirectory  = 'Scan_nJets_{0}'.format(datestr),
            queue         = 'short.q',
            asimov        = ASIMOV,
            )


    #____________________________________________________________________
    if args.plot_nJetsSpectra:

        # ======================================
        # Specify SMXS

        # binBoundaries = [ 0, 1, 2, 3, 4, 5 ]

        # binWidths_hgg = [ binBoundaries_hgg[i+1] - binBoundaries_hgg[i] for i in xrange(len(binBoundaries_hgg)-1) ]
        # binWidths_hzz = [ binBoundaries_hzz[i+1] - binBoundaries_hzz[i] for i in xrange(len(binBoundaries_hzz)-1) ]

        # YR4_totalXS = 55.70628722 # pb

        # shape_hgg = [ 0.208025, 0.234770, 0.165146, 0.218345, 0.087552, 0.059154, 0.022612, 0.004398 ]
        # shape_hzz = [
        #     0.208025,
        #     0.234770,
        #     0.165146 + 0.218345,
        #     0.087552 + 0.059154,
        #     0.022612 + 0.004398,
        #     ]
        # hgg_crosssections = [ s * YR4_totalXS / binWidth for s, binWidth in zip( shape_hgg, binWidths_hgg ) ]
        # hzz_crosssections = [ s * YR4_totalXS / binWidth for s, binWidth in zip( shape_hzz, binWidths_hzz ) ]


        # ======================================
        # Start drawing

        def DrawScan(
                ws,
                scandir,
                pattern,
                name,
                ):

            POIs = Commands.ListPOIs( ws )
            scans = PhysicsCommands.GetScanResults(
                POIs,
                scandir,
                pattern = pattern
                )
            PhysicsCommands.BasicDrawScanResults( POIs, scans, name=name )
            # PhysicsCommands.BasicDrawSpectrum(    POIs, scans, name=name )

            return POIs, scans


        combinationPOIs, combinationscans = DrawScan(
            ws      = LatestPaths.ws_njets_combined_unsplit,
            scandir = 'Scan_nJets_Sep19',
            pattern = 'combinedCard',
            name    = 'combination',
            )

        hggPOIs, hggscans = DrawScan(
            ws      = LatestPaths.ws_njets_hgg_unsplit,
            scandir = 'Scan_nJets_Sep19_0',
            pattern = 'Datacard_13TeV_differential_Njets',
            name    = 'hgg',
            )

        hzzPOIs, hzzscans = DrawScan(
            ws      = LatestPaths.ws_njets_hzz_unsplit,
            scandir = 'Scan_nJets_Sep19_2',
            pattern = 'hzz4l_comb',
            name    = 'hzz',
            )

        fineBinning = [ 0, 1, 2, 3, 4, 5 ]
        hzzBinning  = [ 0, 1, 2, 3, 5 ]

        fineBinWidths = [ fineBinning[i+1]-fineBinning[i] for i in xrange(len(fineBinning)-1) ]
        hzzBinWidths  = [ hzzBinning[i+1]-hzzBinning[i]   for i in xrange(len(hzzBinning)-1) ]

        percentage_fineBinning = [ 0.620010, 0.260148, 0.082666, 0.023464, 0.013713 ]
        percentage_hzzBinning  = [ 0.620010, 0.260148, 0.082666, 0.023464 + 0.013713 ]

        xs_fineBinning = [ xs * LatestPaths.YR4_totalXS / binWidth for xs, binWidth in zip( percentage_fineBinning, fineBinWidths ) ]
        xs_hzzBinning  = [ xs * LatestPaths.YR4_totalXS / binWidth for xs, binWidth in zip( percentage_hzzBinning, hzzBinWidths ) ]


        # For big mu plot, comment out xs's
        PhysicsCommands.BasicCombineSpectra(
            ( 'combination', combinationPOIs, combinationscans,
                ( 'AsBlocks', False ),
                ( 'SetLineColor', 1 ),
                ( 'SetMarkerStyle', 8 ),
                # # Block settings
                # ( 'SetLineColor', 1 ),
                # ( 'SetMarkerStyle', 2 ),
                # # ( 'SetFillColorAlpha', 1, 0.2 ),
                # ( 'SetFillColor', 13 ),
                # # ( 'SetFillStyle', 3544 ),
                # ( 'SetFillStyle', 3345 ),
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
            hzz_SMXS         = xs_hzzBinning,
            hgg_SMXS         = xs_fineBinning,
            combination_SMXS = xs_fineBinning,
            # legendLeft       = True
            )



        PhysicsCommands.BasicCombineSpectra(
            ( 'combination', combinationPOIs, combinationscans,
                ( 'AsBlocks', False ),
                ( 'SetLineColor', 1 ),
                ( 'SetMarkerStyle', 8 ),
                # # Block settings
                # ( 'SetLineColor', 1 ),
                # ( 'SetMarkerStyle', 2 ),
                # # ( 'SetFillColorAlpha', 1, 0.2 ),
                # ( 'SetFillColor', 13 ),
                # # ( 'SetFillStyle', 3544 ),
                # ( 'SetFillStyle', 3345 ),
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
            bottomRatioPlot = True,
            )


        # hggPOIs = Commands.ListPOIs( ws_onlyhgg_unsplit )
        # hggscans = PhysicsCommands.GetScanResults(
        #     hggPOIs,
        #     # scan_ptcombination_hgg_profiled_asimov,
        #     scan_ptcombination_hgg_profiled,
        #     pattern = 'Datacard_13TeV_differential_pT'
        #     )
        # PhysicsCommands.BasicDrawScanResults( hggPOIs, hggscans, name='hgg' )
        # PhysicsCommands.BasicDrawSpectrum( hggPOIs, hggscans, name='hgg' )


        # hzzPOIs = Commands.ListPOIs( ws_onlyhzz_unsplit )
        # hzzscans = PhysicsCommands.GetScanResults(
        #     hzzPOIs,
        #     # scan_ptcombination_hzz_profiled_asimov,
        #     scan_ptcombination_hzz_profiled,
        #     pattern = 'usingDavidsCommand_OAshifted'
        #     )
        # PhysicsCommands.BasicDrawScanResults( hzzPOIs, hzzscans, name='hzz' )
        # PhysicsCommands.BasicDrawSpectrum( hzzPOIs, hzzscans, name='hzz' )


        # combinationPOIs = Commands.ListPOIs( ws_combined_unsplit )
        # combinationscans = PhysicsCommands.GetScanResults(
        #     combinationPOIs,
        #     # scan_ptcombination_combined_profiled_asimov,
        #     scan_ptcombination_combined_profiled,
        #     pattern = 'combinedCard'
        #     )
        # PhysicsCommands.BasicDrawScanResults( combinationPOIs, combinationscans, name='combination' )
        # PhysicsCommands.BasicDrawSpectrum( combinationPOIs, combinationscans, name='combination' )




########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'