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
import LatestBinning

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

    parser.add_argument( '--ptjet_rename_hgg',                   action=CustomAction )
    parser.add_argument( '--ptjet_combineCards',                 action=CustomAction )
    parser.add_argument( '--ptjet_t2ws_unsplit',                 action=CustomAction )
    parser.add_argument( '--ptjet_bestfit_unsplit',              action=CustomAction )
    parser.add_argument( '--ptjet_plot',                         action=CustomAction )

    parser.add_argument( '--njets_plot',                         action=CustomAction )
    parser.add_argument( '--pth_plot',                           action=CustomAction )
    parser.add_argument( '--pth_ggH_plot',                       action=CustomAction )
    parser.add_argument( '--pth_ggH_hbb_plot',                   action=CustomAction )
    parser.add_argument( '--plot_all',                           action=CustomAction )

########################################
# Methods
########################################    

def main( args ):

    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )

    if args.plot_all:
        args.rapidity_plot = True
        args.ptjet_plot    = True
        args.njets_plot    = True
        args.pth_plot      = True
        args.pth_ggH_plot  = True


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


    ########################################
    # Plotting
    ########################################

    #____________________________________________________________________
    if args.rapidity_plot:

        combinationPOIs, combinationscans = DrawScan(
            ws      = LatestPaths.ws_combined_smH_YH,
            scandir = LatestPaths.scan_combined_YH_asimov,
            name    = 'combination',
            # pattern = 'combinedCard',
            )

        hggPOIs, hggscans = DrawScan(
            ws      = LatestPaths.ws_hgg_smH_YH,
            scandir = LatestPaths.scan_hgg_YH_asimov,
            name    = 'hgg',
            # pattern = 'Datacard_13TeV_differential_Njets',
            )

        hzzPOIs, hzzscans = DrawScan(
            ws      = LatestPaths.ws_hzz_smH_YH,
            scandir = LatestPaths.scan_hzz_YH_asimov,
            name    = 'hzz',
            # pattern = 'hzz4l_comb',
            )

        # ======================================
        # Load into suitable containers for two panel plotting

        containers = []

        LatestBinning.obs_yh.Print()

        hgg               = Container()
        hgg.name          = 'hgg'
        hgg.title         = 'H#rightarrow#gamma#gamma'
        hgg.color         = 2
        hgg.POIs          = hggPOIs
        hgg.Scans         = hggscans
        hgg.SMcrosssections = LatestBinning.obs_yh.crosssection_over_binwidth()
        containers.append(hgg)

        hzz               = Container()
        hzz.name          = 'hzz'
        hzz.title         = 'H#rightarrowZZ'
        hzz.color         = 4
        hzz.POIs          = hzzPOIs
        hzz.Scans         = hzzscans
        hzz.SMcrosssections = LatestBinning.obs_yh.crosssection_over_binwidth()
        containers.append(hzz)

        combination               = Container()
        combination.name          = 'combination'
        combination.title         = 'Combination'
        combination.POIs          = combinationPOIs
        combination.Scans         = combinationscans
        combination.SMcrosssections = LatestBinning.obs_yh.crosssection_over_binwidth()
        containers.append(combination)

        PlotCommands.PlotSpectraOnTwoPanel(
            'twoPanel_rapiditySpectrum',
            containers,
            xTitle = '|y_{H}|',
            yTitleTop = '#Delta#sigma/#Delta|y_{H}| (pb)',
            # yMinLimit = 0.07,
            # yMaxExternalTop = 500
            lastBinIsNotOverflow=True,
            )

    #____________________________________________________________________
    if args.ptjet_plot:

        combinationPOIs, combinationscans = DrawScan(
            ws      = LatestPaths.ws_combined_smH_PTJ,
            scandir = LatestPaths.scan_combined_PTJ_asimov,
            name    = 'combination',
            # pattern = 'combinedCard',
            )

        hggPOIs, hggscans = DrawScan(
            ws      = LatestPaths.ws_hgg_smH_PTJ,
            scandir = LatestPaths.scan_hgg_PTJ_asimov,
            name    = 'hgg',
            # pattern = 'Datacard_13TeV_differential_Njets',
            )

        hzzPOIs, hzzscans = DrawScan(
            ws      = LatestPaths.ws_hzz_smH_PTJ,
            scandir = LatestPaths.scan_hzz_PTJ_asimov,
            name    = 'hzz',
            # pattern = 'hzz4l_comb',
            )


        # ======================================
        # Load into suitable containers for two panel plotting

        containers = []

        hgg                 = Container()
        hgg.name            = 'hgg'
        hgg.title           = 'H#rightarrow#gamma#gamma'
        hgg.color           = 2
        hgg.POIs            = hggPOIs
        hgg.Scans           = hggscans
        hgg.SMcrosssections = LatestBinning.obs_ptjet.crosssection_over_binwidth()
        containers.append(hgg)

        hzz                 = Container()
        hzz.name            = 'hzz'
        hzz.title           = 'H#rightarrowZZ'
        hzz.color           = 4
        hzz.POIs            = hzzPOIs
        hzz.Scans           = hzzscans
        hzz.SMcrosssections = LatestBinning.obs_ptjet_hzzBinning.crosssection_over_binwidth()
        containers.append(hzz)

        combination                 = Container()
        combination.name            = 'combination'
        combination.title           = 'Combination'
        combination.POIs            = combinationPOIs
        combination.Scans           = combinationscans
        combination.SMcrosssections = LatestBinning.obs_ptjet.crosssection_over_binwidth()
        containers.append(combination)

        Commands.Warning( 'Skipping first bin (should be the underflow)' )
        for container in containers:
            container.POIs = container.POIs[1:]
            container.Scans = container.Scans[1:]
            container.SMcrosssections = container.SMcrosssections[1:]

        PlotCommands.PlotSpectraOnTwoPanel(
            'twoPanel_ptjetSpectrum',
            containers,
            xTitle = 'p_{T}^{jet} (GeV)',
            yTitleTop = '#Delta#sigma/#Deltap_{T}^{jet} (pb/GeV)',
            # yMinLimit = 0.07,
            yMaxExternalTop = 10,
            xMinExternal = 30.0,
            # yMinLimit    = 0.1
            )


    #____________________________________________________________________
    if args.njets_plot:

        combinationPOIs, combinationscans = DrawScan(
            ws      = LatestPaths.ws_combined_smH_NJ,
            scandir = LatestPaths.scan_combined_NJ,
            name    = 'combination',
            )

        hggPOIs, hggscans = DrawScan(
            ws      = LatestPaths.ws_hgg_smH_NJ,
            scandir = LatestPaths.scan_hgg_NJ,
            name    = 'hgg',
            )

        hzzPOIs, hzzscans = DrawScan(
            ws      = LatestPaths.ws_hzz_smH_NJ,
            scandir = LatestPaths.scan_hzz_NJ,
            name    = 'hzz',
            )

        LatestBinning.obs_njets.Print()
        LatestBinning.obs_njets_hzzBinning.Print()

        containers = []

        hgg                 = Container()
        hgg.name            = 'hgg'
        hgg.title           = 'H#rightarrow#gamma#gamma'
        hgg.color           = 2
        hgg.POIs            = hggPOIs
        hgg.Scans           = hggscans
        hgg.SMcrosssections = LatestBinning.obs_njets.crosssection_over_binwidth()
        containers.append(hgg)

        hzz                 = Container()
        hzz.name            = 'hzz'
        hzz.title           = 'H#rightarrowZZ'
        hzz.color           = 4
        hzz.POIs            = hzzPOIs
        hzz.Scans           = hzzscans
        hzz.SMcrosssections = LatestBinning.obs_njets_hzzBinning.crosssection_over_binwidth()
        containers.append(hzz)

        combination                 = Container()
        combination.name            = 'combination'
        combination.title           = 'Combination'
        combination.POIs            = combinationPOIs
        combination.Scans           = combinationscans
        combination.SMcrosssections = LatestBinning.obs_njets.crosssection_over_binwidth()
        containers.append(combination)

        PlotCommands.PlotSpectraOnTwoPanel(
            'twoPanel_nJetsSpectrum',
            containers,
            xTitle = 'N_{jets}',
            # yMinLimit = 0.07,
            # yMaxExternalTop = 500,
            yTitleTop = '#Delta#sigma/#DeltaN_{jets} (pb)',
            )


    #____________________________________________________________________
    if args.pth_plot:

        combinationPOIs, combinationscans = DrawScan(
            ws      = LatestPaths.ws_combined_smH,
            scandir = LatestPaths.scan_combined_PTH,
            name    = 'combination',
            )

        hggPOIs, hggscans = DrawScan(
            ws      = LatestPaths.ws_hgg_smH,
            scandir = LatestPaths.scan_hgg_PTH,
            name    = 'hgg',
            )

        hzzPOIs, hzzscans = DrawScan(
            ws      = LatestPaths.ws_hzz_smH,
            scandir = LatestPaths.scan_hzz_PTH,
            name    = 'hzz',
            )

        LatestBinning.obs_pth.Print()
        LatestBinning.obs_pth_hzzBinning.Print()

        containers = []

        hgg                 = Container()
        hgg.name            = 'hgg'
        hgg.title           = 'H#rightarrow#gamma#gamma'
        hgg.color           = 2
        hgg.POIs            = hggPOIs
        hgg.Scans           = hggscans
        hgg.SMcrosssections = LatestBinning.obs_pth.crosssection_over_binwidth()
        containers.append(hgg)

        hzz                 = Container()
        hzz.name            = 'hzz'
        hzz.title           = 'H#rightarrowZZ'
        hzz.color           = 4
        hzz.POIs            = hzzPOIs
        hzz.Scans           = hzzscans
        hzz.SMcrosssections = LatestBinning.obs_pth_hzzBinning.crosssection_over_binwidth()
        containers.append(hzz)

        combination                 = Container()
        combination.name            = 'combination'
        combination.title           = 'Combination'
        combination.POIs            = combinationPOIs
        combination.Scans           = combinationscans
        combination.SMcrosssections = LatestBinning.obs_pth.crosssection_over_binwidth()
        containers.append(combination)

        PlotCommands.PlotSpectraOnTwoPanel(
            'twoPanel_pthSpectrum',
            containers,
            xTitle = 'p_{T}^{H} (GeV)',
            yTitleTop = '#Delta#sigma/#Deltap_{T}^{H} (pb/GeV)',
            # 
            # yMinExternalTop = 0.0005,
            # yMaxExternalTop = 110.,
            )

    #____________________________________________________________________
    if args.pth_ggH_plot:

        combinationPOIs, combinationscans = DrawScan(
            ws      = LatestPaths.ws_combined_ggH_xHfixed,
            scandir = LatestPaths.scan_combined_PTH_ggH,
            name    = 'combination',
            )

        hggPOIs, hggscans = DrawScan(
            ws      = LatestPaths.ws_hgg_ggH_xHfixed,
            scandir = LatestPaths.scan_hgg_PTH_ggH,
            name    = 'hgg',
            )

        hzzPOIs, hzzscans = DrawScan(
            ws      = LatestPaths.ws_hzz_ggH_xHfixed,
            scandir = LatestPaths.scan_hzz_PTH_ggH,
            name    = 'hzz',
            )

        Commands.Warning( 'No ggH-only cross sections known yet! Using ggH+xH (smH) cross sections for now' )
        LatestBinning.obs_pth.Print()
        LatestBinning.obs_pth_hzzBinning.Print()

        containers = []

        hgg                 = Container()
        hgg.name            = 'hgg'
        hgg.title           = 'H#rightarrow#gamma#gamma'
        hgg.color           = 2
        hgg.POIs            = hggPOIs
        hgg.Scans           = hggscans
        hgg.SMcrosssections = LatestBinning.obs_pth.crosssection_over_binwidth()
        containers.append(hgg)

        hzz                 = Container()
        hzz.name            = 'hzz'
        hzz.title           = 'H#rightarrowZZ'
        hzz.color           = 4
        hzz.POIs            = hzzPOIs
        hzz.Scans           = hzzscans
        hzz.SMcrosssections = LatestBinning.obs_pth_hzzBinning.crosssection_over_binwidth()
        containers.append(hzz)

        combination                 = Container()
        combination.name            = 'combination'
        combination.title           = 'Combination'
        combination.POIs            = combinationPOIs
        combination.Scans           = combinationscans
        combination.SMcrosssections = LatestBinning.obs_pth.crosssection_over_binwidth()
        containers.append(combination)

        l = PlotCommands.TLatexMultiPanel(
            lambda c: 1.0 - c.GetRightMargin() - 0.01,
            lambda c: 1.0 - c.GetTopMargin() - 0.14,
            '(non-ggH fixed to SM)'
            )
        l.SetTextSize(0.05)
        l.SetTextAlign(33)

        PlotCommands.PlotSpectraOnTwoPanel(
            'twoPanel_pth_ggH_Spectrum',
            containers,
            xTitle = 'p_{T}^{H} (GeV)',
            yTitleTop = '#Delta#sigma^{ggH}/#Deltap_{T}^{H} (pb/GeV)',
            topPanelObjects = [ ( l, '' ) ],
            )


    #____________________________________________________________________
    if args.pth_ggH_hbb_plot:

        if args.asimov:
            combined_scanDir    = LatestPaths.scan_combined_PTH_ggH_asimov
            hgg_scanDir         = LatestPaths.scan_hgg_PTH_ggH_asimov
            hzz_scanDir         = LatestPaths.scan_hzz_PTH_ggH_asimov
            hbb_scanDir         = LatestPaths.scan_hbb_PTH_ggH_asimov
            combWithHbb_scanDir = LatestPaths.scan_combWithHbb_PTH_ggH_asimov
        else:
            combined_scanDir    = LatestPaths.scan_combined_PTH_ggH
            hgg_scanDir         = LatestPaths.scan_hgg_PTH_ggH
            hzz_scanDir         = LatestPaths.scan_hzz_PTH_ggH
            hbb_scanDir         = LatestPaths.scan_hbb_PTH_ggH
            combWithHbb_scanDir = LatestPaths.scan_combWithHbb_PTH_ggH

        combinationPOIs, combinationscans = DrawScan(
            ws      = LatestPaths.ws_combined_ggH_xHfixed,
            scandir = combined_scanDir,
            name    = 'combination' + ('_asimov' if args.asimov else ''),
            )

        hggPOIs, hggscans = DrawScan(
            ws      = LatestPaths.ws_hgg_ggH_xHfixed,
            scandir = hgg_scanDir,
            name    = 'hgg' + ('_asimov' if args.asimov else ''),
            )

        hzzPOIs, hzzscans = DrawScan(
            ws      = LatestPaths.ws_hzz_ggH_xHfixed,
            scandir = hzz_scanDir,
            name    = 'hzz' + ('_asimov' if args.asimov else ''),
            )

        hbbPOIs, hbbscans = DrawScan(
            ws      = LatestPaths.ws_hbb_ggH_xHfixed,
            # scandir = 'out/Scan_PTH_Dec18_xHfixed_hbb_asimov_0', # With all Hbb options
            scandir = hbb_scanDir,
            name    = 'hbb' + ('_asimov' if args.asimov else ''),
            )

        combWithHbbPOIs, combWithHbbscans = DrawScan(
            ws      = LatestPaths.ws_combWithHbb_ggH_xHfixed,
            # scandir = 'out/Scan_PTH_Dec18_xHfixed_combWithHbb_asimov_1',
            scandir = combWithHbb_scanDir,
            name    = 'combWithHbb' + ('_asimov' if args.asimov else ''),
            )


        # Commands.Warning( 'No ggH-only cross sections known yet! Using ggH+xH (smH) cross sections for now' )
        # LatestBinning.obs_pth.Print()
        # LatestBinning.obs_pth_hzzBinning.Print()

        containers = []

        hgg                 = Container()
        hgg.name            = 'hgg'
        hgg.title           = 'H#rightarrow#gamma#gamma'
        hgg.color           = 2
        hgg.POIs            = hggPOIs
        hgg.Scans           = hggscans
        hgg.SMcrosssections = LatestBinning.obs_pth_ggH.crosssection_over_binwidth()
        containers.append(hgg)

        hzz                 = Container()
        hzz.name            = 'hzz'
        hzz.title           = 'H#rightarrowZZ'
        hzz.color           = 4
        hzz.POIs            = hzzPOIs
        hzz.Scans           = hzzscans
        hzz.SMcrosssections = LatestBinning.obs_pth_ggH_hzzBinning.crosssection_over_binwidth()
        containers.append(hzz)

        hbb                 = Container()
        hbb.name            = 'hbb'
        hbb.title           = 'H#rightarrowbb'
        hbb.color           = 8
        hbb.POIs            = hbbPOIs
        hbb.Scans           = hbbscans
        hbb.SMcrosssections = LatestBinning.obs_pth_ggH_hbbBinning.crosssection_over_binwidth()
        containers.append(hbb)

        combination                 = Container()
        combination.name            = 'combination'
        combination.title           = 'Combination'
        combination.color           = 1
        combination.POIs            = combinationPOIs
        combination.Scans           = combinationscans
        combination.SMcrosssections = LatestBinning.obs_pth_ggH.crosssection_over_binwidth()
        containers.append(combination)

        combWithHbb                 = Container()
        combWithHbb.name            = 'combWithHbb'
        combWithHbb.title           = 'Comb. with Hbb'
        combWithHbb.color           = 14
        combWithHbb.POIs            = combWithHbbPOIs
        combWithHbb.Scans           = combWithHbbscans
        combWithHbb.SMcrosssections = LatestBinning.obs_pth_ggH_combWithHbbBinning.crosssection_over_binwidth()
        containers.append(combWithHbb)


        for container in containers:
            if not len(container.POIs) == len(container.SMcrosssections):
                Commands.ThrowError(
                    'For container {2}, found {0} POIs, but {1} SM cross sections; something is misaligned.'.format(
                        len(container.POIs), len(container.SMcrosssections), container.name )
                    + '\n  POIs:  {0}'.format(container.POIs)
                    + '\n  SM xs: {0}'.format(container.SMcrosssections)
                    )



        l = PlotCommands.TLatexMultiPanel(
            lambda c: 1.0 - c.GetRightMargin() - 0.01,
            lambda c: 1.0 - c.GetTopMargin() - 0.14,
            '(non-ggH fixed to SM)'
            )
        l.SetTextSize(0.05)
        l.SetTextAlign(33)

        PlotCommands.PlotSpectraOnTwoPanel(
            'twoPanel_pth_ggH_hbb_Spectrum' + ( '_asimov' if args.asimov else '' ),
            containers,
            xTitle = 'p_{T}^{H} (GeV)',
            yTitleTop = '#Delta#sigma^{ggH}/#Deltap_{T}^{H} (pb/GeV)',
            # topPanelObjects = [ ( l, '' ) ],
            )





#____________________________________________________________________
def DrawScan(
        ws,
        scandir,
        name,
        pattern = '',
        verbose = True,
        ):

    POIs = Commands.ListPOIs( ws )
    POIs.sort( key=Commands.POIsorter )
    if verbose: print 'Sorted POIs:', POIs

    scans = PhysicsCommands.GetScanResults(
        POIs,
        scandir,
        pattern = pattern,
        filterNegatives = True
        )
    PhysicsCommands.BasicDrawScanResults( POIs, scans, name=name )
    # PhysicsCommands.BasicDrawSpectrum(    POIs, scans, name=name )

    return POIs, scans


########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'