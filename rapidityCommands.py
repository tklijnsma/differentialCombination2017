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
from Container import Container
import PlotCommands

from time import strftime
datestr = strftime( '%b%d' )


########################################
# Main
########################################

def AppendParserOptions( parser ):

    parser.add_argument( '--rapidityCommands', default=False )
    class CustomAction(argparse.Action):
        def __init__(self, option_strings, dest, **kwargs):
            # Only 1 argument allowed (this is basically the "store_true" action)
            super(CustomAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            # Make sure the parser knows a combineCommand was entered
            setattr( namespace, 'rapidityCommands', True )
            setattr( namespace, self.dest, True )
            # super(FooAction, self).__call__( parser, namespace, values, option_string=None )
            # setattr( namespace, self.dest, values )

    parser.add_argument( '--rename_hgg_rapidity',                action=CustomAction )
    parser.add_argument( '--rapidity_combineCards',              action=CustomAction )
    parser.add_argument( '--rapidity_t2ws_unsplit',              action=CustomAction )
    parser.add_argument( '--rapidity_bestfit_unsplit',           action=CustomAction )
    parser.add_argument( '--rapidity_plot',                      action=CustomAction )

    parser.add_argument( '--ptjet_rename_hgg',                action=CustomAction )
    parser.add_argument( '--ptjet_combineCards',              action=CustomAction )
    parser.add_argument( '--ptjet_t2ws_unsplit',              action=CustomAction )
    parser.add_argument( '--ptjet_bestfit_unsplit',           action=CustomAction )
    parser.add_argument( '--ptjet_plot',                      action=CustomAction )




########################################
# Methods
########################################    

def main( args ):

    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )


    ########################################
    # PTJET
    ########################################

    #____________________________________________________________________
    if args.ptjet_rename_hgg:
        MergeHGGWDatacards.RenameProcesses_Hgg_differentials(
            LatestPaths.card_hgg_smH_PTJ_unprocessed
            )


    #____________________________________________________________________
    if args.ptjet_combineCards:
        Commands.BasicCombineCards(
            'suppliedInput/combinedCard_PTJ_smH_{0}.txt'.format(datestr),
            'hgg=' + LatestPaths.card_hgg_smH_PTJ,
            'hzz=' + LatestPaths.card_hzz_smH_PTJ
            )


    #____________________________________________________________________
    if args.ptjet_t2ws_unsplit:

        datacard = LatestPaths.card_combined_smH_PTJ
        if args.hgg:
            datacard = LatestPaths.card_hgg_smH_PTJ
        if args.hzz:
            datacard = LatestPaths.card_hzz_smH_PTJ

        if args.hzz:

            # hzz_PTJ_LT30_cat4e
            # hzz_PTJ_30_55_cat4e
            # hzz_PTJ_55_95_cat4e
            # hzz_PTJ_GT95_cat4e

            # smH_PTJ_LT30
            # smH_PTJ_30_55
            # smH_PTJ_55_95
            # smH_PTJ_95_120
            # smH_PTJ_120_200
            # smH_PTJ_GT200

            # smH_PTJ_LT30
            # smH_PTJ_30_55
            # smH_PTJ_55_95
            # smH_PTJ_GT95

            Commands.BasicT2WS(
                datacard,
                manualMaps=[
                    '--PO \'map=.*/smH_PTJ_LT30:r_smH_PTJ_LT30[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/smH_PTJ_30_55:r_smH_PTJ_30_55[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/smH_PTJ_55_95:r_smH_PTJ_55_95[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/smH_PTJ_95_120:r_smH_PTJ_GT95[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/smH_PTJ_120_200:r_smH_PTJ_GT95[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/smH_PTJ_GT200:r_smH_PTJ_GT95[1.0,-1.0,4.0]\'',
                    ],
                suffix = '_ptjet'
                )

        else:
            Commands.BasicT2WS(
                datacard,
                smartMaps = [
                    ( r'.*/smH_PTJ_([pm\d\_GELT]+)', r'r_smH_PTJ_\1[1.0,-1.0,4.0]' )
                    ],
                # suffix = '_ptjet'
                )

    #____________________________________________________________________
    if args.ptjet_bestfit_unsplit:

        ASIMOV = True
        # ASIMOV = False

        nPoints       = 55
        nPointsPerJob = 5

        ws = LatestPaths.ws_combined_smH_PTJ
        if args.hgg:
            ws = LatestPaths.ws_hgg_smH_PTJ
        if args.hzz:
            ws = LatestPaths.ws_hzz_smH_PTJ
            nPointsPerJob = nPoints

        Commands.BasicCombineTool(
            ws,
            POIpattern    = '*',
            nPoints       = nPoints,
            nPointsPerJob = nPointsPerJob,
            jobDirectory  = 'Scan_PTJ_{0}'.format(datestr),
            queue         = 'short.q',
            asimov        = ASIMOV,
            )



    ########################################
    # AbsRapidity
    ########################################

    #____________________________________________________________________
    if args.rename_hgg_rapidity:
        MergeHGGWDatacards.RenameProcesses_Hgg_differentials(
            LatestPaths.card_hgg_smH_YH_unprocessed,
            )

    #____________________________________________________________________
    if args.rapidity_combineCards:

        Commands.BasicCombineCards(
            'suppliedInput/combinedCard_YH_smH_{0}.txt'.format(datestr),
            'hgg=' + LatestPaths.card_hgg_smH_YH,
            'hzz=' + LatestPaths.card_hzz_smH_YH
            )


    #____________________________________________________________________
    if args.rapidity_t2ws_unsplit:

        datacard = LatestPaths.card_combined_smH_YH
        if args.hgg:
            datacard = LatestPaths.card_hgg_smH_YH
        if args.hzz:
            datacard = LatestPaths.card_hzz_smH_YH

        # smH_YH_0p0_0p15
        # smH_YH_0p15_0p30
        # smH_YH_0p30_0p60
        # smH_YH_0p60_0p90
        # smH_YH_0p90_1p20
        # smH_YH_1p20_2p50

        Commands.BasicT2WS(
            datacard,
            smartMaps = [
                ( r'.*/smH_YH_([pm\d\_GE]+)', r'r_smH_YH_\1[1.0,-1.0,4.0]' )
                ],
            )

    #____________________________________________________________________
    if args.rapidity_bestfit_unsplit:

        # ASIMOV = False
        ASIMOV = True

        nPoints       = 55
        nPointsPerJob = 5

        ws = LatestPaths.ws_combined_smH_YH
        if args.hgg:
            ws = LatestPaths.ws_hgg_smH_YH
        if args.hzz:
            ws = LatestPaths.ws_hzz_smH_YH
            nPointsPerJob = nPoints

        Commands.BasicCombineTool(
            ws,
            POIpattern    = '*',
            nPoints       = nPoints,
            nPointsPerJob = nPointsPerJob,
            jobDirectory  = 'Scan_YH_{0}'.format(datestr),
            queue         = 'short.q',
            asimov        = ASIMOV,
            )


    #____________________________________________________________________
    def DrawScan(
            ws,
            scandir,
            name,
            pattern = '',
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

    # Use pt cross sections for now
    binBoundaries_hgg = [ 0., 15., 30., 45., 85., 125., 200., 350., 1000. ]
    binBoundaries_hzz = [ 0., 15., 30., 85., 200., 1000. ]
    binWidths_hgg = [ binBoundaries_hgg[i+1] - binBoundaries_hgg[i] for i in xrange(len(binBoundaries_hgg)-1) ]
    binWidths_hzz = [ binBoundaries_hzz[i+1] - binBoundaries_hzz[i] for i in xrange(len(binBoundaries_hzz)-1) ]

    YR4_totalXS = LatestPaths.YR4_totalXS
    shape_hgg = [ 0.208025, 0.234770, 0.165146, 0.218345, 0.087552, 0.059154, 0.022612, 0.004398 ]
    shape_hzz = [ 0.208025, 0.234770, 0.165146 + 0.218345, 0.087552 + 0.059154, 0.022612 + 0.004398, ]
    hgg_crosssections = [ s * YR4_totalXS / binWidth for s, binWidth in zip( shape_hgg, binWidths_hgg ) ]
    hzz_crosssections = [ s * YR4_totalXS / binWidth for s, binWidth in zip( shape_hzz, binWidths_hzz ) ]

    #____________________________________________________________________
    if args.rapidity_plot:

        combinationPOIs, combinationscans = DrawScan(
            ws      = LatestPaths.ws_combined_smH_YH,
            scandir = LatestPaths.scan_combined_YH,
            name    = 'combination',
            # pattern = 'combinedCard',
            )

        hggPOIs, hggscans = DrawScan(
            ws      = LatestPaths.ws_hgg_smH_YH,
            scandir = LatestPaths.scan_hgg_YH,
            name    = 'hgg',
            # pattern = 'Datacard_13TeV_differential_Njets',
            )

        hzzPOIs, hzzscans = DrawScan(
            ws      = LatestPaths.ws_hzz_smH_YH,
            scandir = LatestPaths.scan_hzz_YH,
            name    = 'hzz',
            # pattern = 'hzz4l_comb',
            )

        # ======================================
        # Load into suitable containers for two panel plotting

        Commands.Warning( 'No cross sections available yet; place holder!!' )
        xs_per_bin = hgg_crosssections[:len(combinationPOIs)]

        containers = []

        hgg               = Container()
        hgg.name          = 'hgg'
        hgg.title         = 'H#rightarrow#gamma#gamma'
        hgg.color         = 2
        hgg.POIs          = hggPOIs
        hgg.Scans         = hggscans
        hgg.SMcrosssections = xs_per_bin
        containers.append(hgg)

        hzz               = Container()
        hzz.name          = 'hzz'
        hzz.title         = 'H#rightarrowZZ'
        hzz.color         = 4
        hzz.POIs          = hzzPOIs
        hzz.Scans         = hzzscans
        hzz.SMcrosssections = xs_per_bin
        containers.append(hzz)

        combination               = Container()
        combination.name          = 'combination'
        combination.title         = 'Combination'
        combination.POIs          = combinationPOIs
        combination.Scans         = combinationscans
        combination.SMcrosssections = xs_per_bin
        containers.append(combination)

        PlotCommands.PlotSpectraOnTwoPanel(
            'rapiditySpectrum_twoPanel',
            containers,
            xTitle = '|y_{H}|',
            # yMinLimit = 0.07,
            # yMaxExternalTop = 500
            lastBinIsNotOverflow=True,
            )

    #____________________________________________________________________
    if args.ptjet_plot:

        combinationPOIs, combinationscans = DrawScan(
            ws      = LatestPaths.ws_combined_smH_PTJ,
            scandir = LatestPaths.scan_combined_PTJ,
            name    = 'combination',
            # pattern = 'combinedCard',
            )

        hggPOIs, hggscans = DrawScan(
            ws      = LatestPaths.ws_hgg_smH_PTJ,
            scandir = LatestPaths.scan_hgg_PTJ,
            name    = 'hgg',
            # pattern = 'Datacard_13TeV_differential_Njets',
            )

        Commands.Warning( 'Made mistake in t2ws, now correcting by renaming POIs...' )
        combinationPOIs = [ POI.replace('YH','PTJ') for POI in combinationPOIs ]
        hggPOIs = [ POI.replace('YH','PTJ') for POI in hggPOIs ]


        hzzPOIs, hzzscans = DrawScan(
            ws      = LatestPaths.ws_hzz_smH_PTJ,
            scandir = LatestPaths.scan_hzz_PTJ,
            name    = 'hzz',
            # pattern = 'hzz4l_comb',
            )

        # ======================================
        # Load into suitable containers for two panel plotting

        Commands.Warning( 'No cross sections available yet; place holder!!' )
        xs_per_bin = hgg_crosssections[:len(combinationPOIs)]

        containers = []

        hgg               = Container()
        hgg.name          = 'hgg'
        hgg.title         = 'H#rightarrow#gamma#gamma'
        hgg.color         = 2
        hgg.POIs          = hggPOIs
        hgg.Scans         = hggscans
        hgg.SMcrosssections = xs_per_bin
        containers.append(hgg)

        hzz               = Container()
        hzz.name          = 'hzz'
        hzz.title         = 'H#rightarrowZZ'
        hzz.color         = 4
        hzz.POIs          = hzzPOIs
        hzz.Scans         = hzzscans
        hzz.SMcrosssections = xs_per_bin
        containers.append(hzz)

        combination               = Container()
        combination.name          = 'combination'
        combination.title         = 'Combination'
        combination.POIs          = combinationPOIs
        combination.Scans         = combinationscans
        combination.SMcrosssections = xs_per_bin
        containers.append(combination)

        PlotCommands.PlotSpectraOnTwoPanel(
            'ptjetSpectrum_twoPanel',
            containers,
            xTitle = 'p_{T}^{jet}',
            # yMinLimit = 0.07,
            # yMaxExternalTop = 500
            xMinExternal = 0.0,
            )



########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'