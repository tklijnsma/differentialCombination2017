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

from LatestPaths import *

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

    parser.add_argument( '--renamehgg',                       action=CustomAction )
    parser.add_argument( '--mergehgg',                        action=CustomAction )
    parser.add_argument( '--t2ws',                            action=CustomAction )
    parser.add_argument( '--t2ws_OLD',                        action=CustomAction )
    parser.add_argument( '--bestfit',                         action=CustomAction )
    parser.add_argument( '--combineCards',                    action=CustomAction )
    parser.add_argument( '--couplingImportanceHistogram',     action=CustomAction )
    parser.add_argument( '--doFastscan',                      action=CustomAction )
    parser.add_argument( '--OutsideAcceptancetoSignal',       action=CustomAction )
    parser.add_argument( '--reparseT2WS',                     action=CustomAction )
    parser.add_argument( '--reparseBestfit',                  action=CustomAction )
    parser.add_argument( '--smartMapTest',                    action=CustomAction )


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

    parser.add_argument( '--RenameHggProcesses_Aug21',      action=CustomAction )
    parser.add_argument( '--RenumberHzzProcesses_Aug21',    action=CustomAction )
    parser.add_argument( '--CombineCards_Aug21',            action=CustomAction )

    parser.add_argument( '--plot_ptSpectra',                action=CustomAction )

    # group = parser.add_mutually_exclusive_group(required=False)
    # group.add_argument( '--latest', dest='latest', action='store_true', default=True )
    # group.add_argument( '--older',  dest='latest', action='store_false' )



########################################
# Methods
########################################    

def main( args ):


    #____________________________________________________________________
    if args.RenameHggProcesses_Aug21:
        MergeHGGWDatacards.RenameProcesses_Aug21(
            card_onlyhgg_split_both_unrenamed,
            )


    if args.RenumberHzzProcesses_Aug21:
        MergeHGGWDatacards.RenumberProcessesHZZ_Aug21(
            card_onlyhzz_split_unrenamed,
            )


    if args.CombineCards_Aug21:
        Commands.BasicCombineCards(
            'suppliedInput/combinedCard_{0}.txt'.format(datestr),
            card_onlyhgg_split_both_renamed,
            card_onlyhzz_split_renamed
            )


    #____________________________________________________________________
    if args.t2ws_combined_unsplit:
        Commands.BasicT2WS(
            card_combined_unsplit,
            smartMaps = [
                ( r'.*/smH_PTH_([\d\_GT]+)', r'r_smH_PTH_\1[1.0,-1.0,4.0]' )
                ],
            )

    if args.t2ws_hgg_unsplit:
        Commands.BasicT2WS(
            card_onlyhgg_unsplit_renamed,
            smartMaps = [
                ( r'.*/smH_PTH_([\d\_GT]+)', r'r_smH_PTH_\1[1.0,-1.0,4.0]' )
                ],
            )

    if args.t2ws_hzz_unsplit:
        Commands.BasicT2WS(
            card_onlyhzz_unsplit_OAsignal,
            # smartMaps = [
            #     ( r'.*/smH_PTH_([\d\_GT]+)', r'r_smH_PTH_\1[1.0,-1.0,4.0]' )
            #     ],
            manualMaps = [
                '--PO \'map=.*/smH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_30_45:r_smH_PTH_30_85[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTH_45_85:r_smH_PTH_30_85[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTH_85_125:r_smH_PTH_85_200[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTH_125_200:r_smH_PTH_85_200[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTH_200_350:r_smH_PTH_GT200[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTH_GT350:r_smH_PTH_GT200[1.0,-1.0,4.0]\'',
                ],
            )

    #____________________________________________________________________

    if args.t2ws_combined_split:

        Commands.BasicT2WS(
            card_combined_split,
            smartMaps = [
                ( r'.*/.*H_PTH_([\d\_GT]+)', r'r_smH_PTH_\1[1.0,-1.0,4.0]' )
                ],
            )


    if args.t2ws_hgg_split:

        # Specific for HZZ
        Commands.BasicT2WS(
            card_onlyhgg_split_both_renamed,
            smartMaps = [
                ( r'.*/.*H_PTH_([\d\_GT]+)', r'r_smH_PTH_\1[1.0,-1.0,4.0]' )
                ],
            )


    if args.t2ws_hzz_split:

        # Specific for HZZ
        Commands.BasicT2WS(
            # 'suppliedInput/hzz_ggH_xH_split_Jun26/hzz4l_all_13TeV_xs.txt',
            card_onlyhzz_split_renamed,
            manualMaps=[
                '--PO \'map=.*/ggH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
                '--PO \'map=.*/xH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
                '--PO \'map=.*/xH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_30_45:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/xH_PTH_30_45:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_45_85:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/xH_PTH_45_85:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_85_125:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/xH_PTH_85_125:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_125_200:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/xH_PTH_125_200:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_200_350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/xH_PTH_200_350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_GT350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/xH_PTH_GT350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
                ],
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



    # Commands.BasicCombineTool(
    #     'workspaces_May15/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.root',
    #     POIpattern    = '*',
    #     nPoints       = 40,
    #     nPointsPerJob = 4,
    #     # notOnBatch    = True,
    #     jobDirectory  = 'Scan_{0}'.format(datestr),
    #     queue         = 'short.q',
    #     )



    #____________________________________________________________________

    ASIMOV = True

    if args.scan_combined_unsplit:
        Commands.BasicCombineTool(
            ws_combined_unsplit,
            POIpattern    = '*',
            nPoints       = 39,
            nPointsPerJob = 3,
            jobDirectory  = 'Scan_{0}'.format(datestr),
            queue         = 'short.q',
            asimov        = ASIMOV,
            )

    if args.scan_hgg_unsplit:
        Commands.BasicCombineTool(
            ws_onlyhgg_unsplit,
            POIpattern    = '*',
            nPoints       = 39,
            nPointsPerJob = 3,
            jobDirectory  = 'Scan_{0}'.format(datestr),
            queue         = 'short.q',
            asimov        = ASIMOV,
            )
        
    if args.scan_hzz_unsplit:
        Commands.BasicCombineTool(
            ws_onlyhzz_unsplit,
            POIpattern    = '*',
            nPoints       = 39,
            nPointsPerJob = 13,
            jobDirectory  = 'Scan_{0}'.format(datestr),
            queue         = 'short.q',
            asimov        = ASIMOV,
            POIRange      = [ 0.0, 3.0 ]
            )



    if args.plot_ptSpectra:

        # ======================================
        # Specify SMXS

        binBoundaries_hgg = [ 0., 15., 30., 45., 85., 125., 200., 350., 1000. ]
        binBoundaries_hzz = [ 0., 15., 30., 85., 200., 1000. ]

        binWidths_hgg = [ binBoundaries_hgg[i+1] - binBoundaries_hgg[i] for i in xrange(len(binBoundaries_hgg)-1) ]
        binWidths_hzz = [ binBoundaries_hzz[i+1] - binBoundaries_hzz[i] for i in xrange(len(binBoundaries_hzz)-1) ]

        YR4_totalXS = 55.70628722 # pb

        shape_hgg = [ 0.208025, 0.234770, 0.165146, 0.218345, 0.087552, 0.059154, 0.022612, 0.004398 ]
        shape_hzz = [
            0.208025,
            0.234770,
            0.165146 + 0.218345,
            0.087552 + 0.059154,
            0.022612 + 0.004398,
            ]
        hgg_crosssections = [ s * YR4_totalXS / binWidth for s, binWidth in zip( shape_hgg, binWidths_hgg ) ]
        hzz_crosssections = [ s * YR4_totalXS / binWidth for s, binWidth in zip( shape_hzz, binWidths_hzz ) ]


        # ======================================
        # Start drawing

        hggPOIs = Commands.ListPOIs( ws_onlyhgg_unsplit )
        hggscans = PhysicsCommands.GetScanResults(
            hggPOIs,
            # scan_ptcombination_hgg_profiled_asimov,
            scan_ptcombination_hgg_profiled,
            pattern = 'Datacard_13TeV_differential_pT'
            )
        PhysicsCommands.BasicDrawScanResults( hggPOIs, hggscans, name='hgg' )
        PhysicsCommands.BasicDrawSpectrum( hggPOIs, hggscans, name='hgg' )


        hzzPOIs = Commands.ListPOIs( ws_onlyhzz_unsplit )
        hzzscans = PhysicsCommands.GetScanResults(
            hzzPOIs,
            # scan_ptcombination_hzz_profiled_asimov,
            scan_ptcombination_hzz_profiled,
            pattern = 'usingDavidsCommand_OAshifted'
            )
        PhysicsCommands.BasicDrawScanResults( hzzPOIs, hzzscans, name='hzz' )
        PhysicsCommands.BasicDrawSpectrum( hzzPOIs, hzzscans, name='hzz' )


        combinationPOIs = Commands.ListPOIs( ws_combined_unsplit )
        combinationscans = PhysicsCommands.GetScanResults(
            combinationPOIs,
            # scan_ptcombination_combined_profiled_asimov,
            scan_ptcombination_combined_profiled,
            pattern = 'combinedCard'
            )
        PhysicsCommands.BasicDrawScanResults( combinationPOIs, combinationscans, name='combination' )
        PhysicsCommands.BasicDrawSpectrum( combinationPOIs, combinationscans, name='combination' )


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
            hzz_SMXS         = hzz_crosssections,
            hgg_SMXS         = hgg_crosssections,
            combination_SMXS = hgg_crosssections,
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




    ########################################
    # Older
    ########################################

    #____________________________________________________________________
    if args.OutsideAcceptancetoSignal:
        OneOfCommands.ChangeOutsideAcceptanceToSignalProcess(
            'suppliedInput/fromDavid/PTH_May15/hzz4l_comb_13TeV_xs.txt',
            )

    #____________________________________________________________________
    if args.renamehgg:

        if args.latest:

            MergeHGGWDatacards.RenameProductionModeHgg(
                'ggH',
                'suppliedInput/fromVittorio/pT_ggHonly_Jun26/CMS-HGG_sigfit_differential_pT_moriond17_ggHonly_v2.root'
                )

            MergeHGGWDatacards.RenameProductionModeHgg(
                'xH',
                'suppliedInput/fromVittorio/differential_pT_moriond17_HxOnly_July19/CMS-HGG_sigfit_differential_pT_moriond17_HxOnly.root'
                )

        else:

            MergeHGGWDatacards.Rename_fea(
                'ggH',
                'suppliedInput/fromVittorio/pT_ggHonly_Jun26/CMS-HGG_sigfit_differential_pT_moriond17_ggHonly_v2.root'
                )

            MergeHGGWDatacards.Rename_fea(
                'xH',
                'suppliedInput/fromVittorio/differential_pT_moriond17_HxOnly_July19/CMS-HGG_sigfit_differential_pT_moriond17_HxOnly.root'
                )

    #____________________________________________________________________
    if args.mergehgg:

        if args.latest:

            xHCard  = 'suppliedInput/fromVittorio/differential_pT_moriond17_HxOnly_July19/Datacard_13TeV_differential_pT_moriond17_HxOnly.txt'
            ggHCard = 'suppliedInput/fromVittorio/pT_ggHonly_Jun26/Datacard_13TeV_differential_pT_moriond17_ggHonly_v2.txt'


            # Pre-process the text datacards

            xH_DCrenamed = MergeHGGWDatacards.RenameProcesses( 'xH',  xHCard,  outdatacard='auto',
                globalReplace = [
                    (
                        'CMS-HGG_sigfit_differential_pT_moriond17_HxOnly.root',
                        'CMS-HGG_sigfit_differential_pT_moriond17_HxOnly_FullyRenamed.root'
                        ),
                    (
                        'InsideAcceptance',
                        'xH_InsideAcceptance'
                        ),
                    (
                        'OutsideAcceptance',
                        'xH_OutsideAcceptance'
                        ),
                    ]
                )
            ggH_DCrenamed = MergeHGGWDatacards.RenameProcesses( 'ggH', ggHCard, outdatacard='auto',
                globalReplace = [
                    (
                        'CMS-HGG_sigfit_differential_pT_moriond17_ggHonly_v2.root',
                        'CMS-HGG_sigfit_differential_pT_moriond17_ggHonly_v2_FullyRenamed.root'
                        ),
                    (
                        'InsideAcceptance',
                        'ggH_InsideAcceptance'
                        ),
                    (
                        'OutsideAcceptance',
                        'ggH_OutsideAcceptance'
                        ),
                    ]
                )


            # Read into container
            xh  = MergeHGGWDatacards.GetDatacardContainer( xH_DCrenamed )
            ggh = MergeHGGWDatacards.GetDatacardContainer( ggH_DCrenamed )
            xh.process = 'xH'
            ggh.process = 'ggH'

            # Include one directory upwards in the path of the root files
            MergeHGGWDatacards.ExtendPathOfRootFiles( xh, dirname(xHCard) )
            MergeHGGWDatacards.ExtendPathOfRootFiles( ggh, dirname(ggHCard) )

            # Merge
            merged = MergeHGGWDatacards.MergeCards( ggh, xh )

            outpath = join( dirname(ggHCard), '..', 'splithggMerged_{0}.txt'.format(datestr) )
            MergeHGGWDatacards.ParseDataContainer( merged, writeToFile = outpath )



        else:

            xHCard  = 'suppliedInput/fromVittorio/differential_pT_moriond17_HxOnly_July19/Datacard_13TeV_differential_pT_moriond17_HxOnly.txt'
            ggHCard = 'suppliedInput/fromVittorio/pT_ggHonly_Jun26/Datacard_13TeV_differential_pT_moriond17_ggHonly_v2.txt'
            # xHCard  = 'suppliedInput/fromVittorio/differential_pT_moriond17_HxOnly_July19/CMS-HGG_sigfit_differential_pT_moriond17_HxOnly_feaRenamed.root'
            # ggHCard = 'suppliedInput/fromVittorio/pT_ggHonly_Jun26/CMS-HGG_sigfit_differential_pT_moriond17_ggHonly_v2_feaRenamed.root'

            MergeHGGWDatacards.Merge_xH_ggH_hgg_cards( xHCard, ggHCard )

            MergeHGGWDatacards.RenameProcesses( 'xH',  xHCard,  outdatacard='auto',
                globalReplace = [(
                    'CMS-HGG_sigfit_differential_pT_moriond17_HxOnly.root',
                    'CMS-HGG_sigfit_differential_pT_moriond17_HxOnly_feaRenamed.root'
                    )]
                )
            MergeHGGWDatacards.RenameProcesses( 'ggH', ggHCard, outdatacard='auto',
                globalReplace = [(
                    'CMS-HGG_sigfit_differential_pT_moriond17_ggHonly_v2.root',
                    'CMS-HGG_sigfit_differential_pT_moriond17_ggHonly_v2_feaRenamed.root'
                    )]
                )

            MergeHGGWDatacards.RenameProcesses(
                'smH',
                'suppliedInput/fromVittorio/pT_unsplit_Mar21/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5.txt',
                outdatacard='auto' )

    #____________________________________________________________________
    if args.combineCards:

        if args.latest:

            Commands.BasicCombineCards(
                'suppliedInput/combinedCard_{0}.txt'.format(datestr),
                'suppliedInput/fromVittorio/pT_unsplit_Mar21/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamedProcesses_Jul21.txt',
                'suppliedInput/fromDavid/PTH_May15/hzz4l_comb_13TeV_xs_processesShifted.txt'
                )

        else:
            print 'Nothing for args.combineCards'


    #____________________________________________________________________
    if args.smartMapTest:

        Commands.TestMode()

        card = 'suppliedInput/combinedCard_Jul26.txt'

        Commands.BasicT2WS(
            card,
            manualMaps=[
                '--PO \'map=.*/smH_PTH_0_15:r_smH_PTH_0_15[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTH_15_30:r_smH_PTH_15_30[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTH_30_45:r_smH_PTH_30_45[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTH_45_85:r_smH_PTH_45_85[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTH_85_125:r_smH_PTH_85_125[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTH_125_200:r_smH_PTH_125_200[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTH_200_350:r_smH_PTH_200_350[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTH_GT350:r_smH_PTH_GT350[1.0,-1.0,4.0]\'',
                ],
            )

        Commands.BasicT2WS(
            card,
            smartMaps = [
                ( r'.*/smH_PTH_([\d\_GT]+)', r'r_smH_PTH_\1[1.0,-1.0,4.0]' )
                ],
            )

    #____________________________________________________________________
    if args.reparseT2WS:

        GGH = False
        XH  = True

        if GGH:

            ggHCard = 'suppliedInput/fromVittorio/pT_ggHonly_Jun26/Datacard_13TeV_differential_pT_moriond17_ggHonly_v2_renamedProcesses_Jul31.txt'

            ggHCardContainer = MergeHGGWDatacards.GetDatacardContainer( ggHCard )

            ggHCardReparsed  = join( dirname(ggHCard), basename(ggHCard).replace( '.txt', '_reparsed.txt' ) )
            MergeHGGWDatacards.ParseDataContainer(
                ggHCardContainer,
                writeToFile = ggHCardReparsed
                )

            Commands.BasicT2WS(
                ggHCardReparsed,
                manualMaps=[
                    '--PO \'map=.*/ggH_PTH_0_15:r_ggH_PTH_0_15[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/ggH_PTH_15_30:r_ggH_PTH_15_30[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/ggH_PTH_30_45:r_ggH_PTH_30_45[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/ggH_PTH_45_85:r_ggH_PTH_45_85[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/ggH_PTH_85_125:r_ggH_PTH_85_125[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/ggH_PTH_125_200:r_ggH_PTH_125_200[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/ggH_PTH_200_350:r_ggH_PTH_200_350[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/ggH_PTH_GT350:r_ggH_PTH_GT350[1.0,-1.0,4.0]\'',
                    ],
                )

        if XH:

            xHCard = 'suppliedInput/fromVittorio/differential_pT_moriond17_HxOnly_July19/Datacard_13TeV_differential_pT_moriond17_HxOnly.txt'

            xHCardContainer = MergeHGGWDatacards.GetDatacardContainer( xHCard )

            xHCardReparsed  = join( dirname(xHCard), basename(xHCard).replace( '.txt', '_reparsed.txt' ) )
            MergeHGGWDatacards.ParseDataContainer(
                xHCardContainer,
                writeToFile = xHCardReparsed
                )

            Commands.BasicT2WS(
                xHCardReparsed,
                manualMaps=[
                    '--PO \'map=.*/InsideAcceptance_genPt_0p0_15p0:r_xH_PTH_0_15[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_15p0_30p0:r_xH_PTH_15_30[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_30p0_45p0:r_xH_PTH_30_45[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_45p0_85p0:r_xH_PTH_45_85[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_85p0_125p0:r_xH_PTH_85_125[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_125p0_200p0:r_xH_PTH_125_200[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_200p0_350p0:r_xH_PTH_200_350[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_350p0_10000p0:r_xH_PTH_GT350[1.0,-1.0,4.0]\'',
                    ],
                )



    if args.reparseBestfit:

        GGH = False
        XH  = True

        if GGH:
            ws = 'workspaces_Aug02/Datacard_13TeV_differential_pT_moriond17_ggHonly_v2_renamedProcesses_Jul31_reparsed.root'

        if XH:
            ws = 'workspaces_Aug02/Datacard_13TeV_differential_pT_moriond17_HxOnly_reparsed.root'


        Commands.BasicBestfit(
            ws,
            extraOptions = [
                '--minimizerStrategy 2',
                '-v 2',
                '-m 125',
                '--floatOtherPOIs=1',
                ]
            )




    #____________________________________________________________________
    if args.t2ws_OLD:

        # ======================================
        # July 24

        if args.latest:
            print 'Running latest'

            Commands.BasicT2WS(
                'suppliedInput/fromVittorio/splithggMerged_Jul31.txt',
                manualMaps=[
                    '--PO \'map=.*/ggH_PTH_0_15:r_smH_PTH_0_15[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/ggH_PTH_15_30:r_smH_PTH_15_30[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/ggH_PTH_30_45:r_smH_PTH_30_45[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/ggH_PTH_45_85:r_smH_PTH_45_85[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/ggH_PTH_85_125:r_smH_PTH_85_125[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/ggH_PTH_125_200:r_smH_PTH_125_200[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/ggH_PTH_200_350:r_smH_PTH_200_350[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/ggH_PTH_GT350:r_smH_PTH_GT350[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/xH_PTH_0_15:r_smH_PTH_0_15[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/xH_PTH_15_30:r_smH_PTH_15_30[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/xH_PTH_30_45:r_smH_PTH_30_45[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/xH_PTH_45_85:r_smH_PTH_45_85[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/xH_PTH_85_125:r_smH_PTH_85_125[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/xH_PTH_125_200:r_smH_PTH_125_200[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/xH_PTH_200_350:r_smH_PTH_200_350[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/xH_PTH_GT350:r_smH_PTH_GT350[1.0,-1.0,4.0]\'',
                    ],
                )


        else:
            print 'Running older'

            xHCard_directlyFromVittorio = 'suppliedInput/fromVittorio/differential_pT_moriond17_HxOnly_July19/Datacard_13TeV_differential_pT_moriond17_HxOnly.txt'
            Commands.BasicT2WS(
                xHCard_directlyFromVittorio,
                manualMaps=[
                    '--PO \'map=.*/InsideAcceptance_genPt_0p0_15p0:r_xH_PTH_0_15[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_15p0_30p0:r_xH_PTH_15_30[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_30p0_45p0:r_xH_PTH_30_45[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_45p0_85p0:r_xH_PTH_45_85[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_85p0_125p0:r_xH_PTH_85_125[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_125p0_200p0:r_xH_PTH_125_200[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_200p0_350p0:r_xH_PTH_200_350[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_350p0_10000p0:r_xH_PTH_GT350[1.0,-1.0,4.0]\'',
                    ],
                )

            sys.exit()


            ggHCard_directlyFromVittorio = 'suppliedInput/fromVittorio/pT_ggHonly_Jun26/Datacard_13TeV_differential_pT_moriond17_ggHonly_v2.txt'

            Commands.BasicT2WS(
                ggHCard_directlyFromVittorio,
                manualMaps=[
                    '--PO \'map=.*/InsideAcceptance_genPt_0p0_15p0:r_ggH_PTH_0_15[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_15p0_30p0:r_ggH_PTH_15_30[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_30p0_45p0:r_ggH_PTH_30_45[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_45p0_85p0:r_ggH_PTH_45_85[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_85p0_125p0:r_ggH_PTH_85_125[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_125p0_200p0:r_ggH_PTH_125_200[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_200p0_350p0:r_ggH_PTH_200_350[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/InsideAcceptance_genPt_350p0_10000p0:r_ggH_PTH_GT350[1.0,-1.0,4.0]\'',
                    ],
                )


            # ggHCard = 'suppliedInput/fromVittorio/pT_ggHonly_Jun26/Datacard_13TeV_differential_pT_moriond17_ggHonly_v2_renamedProcesses_Jul31.txt'
            # Commands.BasicT2WS(
            #     ggHCard,
            #     manualMaps=[
            #         '--PO \'map=.*/ggH_PTH_0_15:r_ggH_PTH_0_15[1.0,-1.0,4.0]\'',
            #         '--PO \'map=.*/ggH_PTH_15_30:r_ggH_PTH_15_30[1.0,-1.0,4.0]\'',
            #         '--PO \'map=.*/ggH_PTH_30_45:r_ggH_PTH_30_45[1.0,-1.0,4.0]\'',
            #         '--PO \'map=.*/ggH_PTH_45_85:r_ggH_PTH_45_85[1.0,-1.0,4.0]\'',
            #         '--PO \'map=.*/ggH_PTH_85_125:r_ggH_PTH_85_125[1.0,-1.0,4.0]\'',
            #         '--PO \'map=.*/ggH_PTH_125_200:r_ggH_PTH_125_200[1.0,-1.0,4.0]\'',
            #         '--PO \'map=.*/ggH_PTH_200_350:r_ggH_PTH_200_350[1.0,-1.0,4.0]\'',
            #         '--PO \'map=.*/ggH_PTH_GT350:r_ggH_PTH_GT350[1.0,-1.0,4.0]\'',
            #         ],
            #     )


            # unsplitCard = 'suppliedInput/fromVittorio/pT_unsplit_Mar21/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamedProcesses_Jul21.txt'
            # Commands.BasicT2WS(
            #     unsplitCard,
            #     manualMaps=[
            #         '--PO \'map=.*/smH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/smH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/smH_PTH_30_45:r_smH_PTH_30_45[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/smH_PTH_45_85:r_smH_PTH_45_85[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/smH_PTH_85_125:r_smH_PTH_85_125[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/smH_PTH_125_200:r_smH_PTH_125_200[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/smH_PTH_200_350:r_smH_PTH_200_350[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/smH_PTH_GT350:r_smH_PTH_GT350[1.0,0.0,3.0]\'',
            #         ],
            #     )



            # ======================================
            # Even older

            # # Specific for HZZ
            # Commands.BasicT2WS(
            #     'suppliedInput/hzz_ggH_xH_split_Jun26/hzz4l_all_13TeV_xs.txt',
            #     manualMaps=[
            #         '--PO \'map=.*/ggH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/xH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/ggH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/xH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/ggH_PTH_30_45:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/xH_PTH_30_45:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/ggH_PTH_45_85:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/xH_PTH_45_85:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/ggH_PTH_85_125:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/xH_PTH_85_125:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/ggH_PTH_125_200:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/xH_PTH_125_200:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/ggH_PTH_200_350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/xH_PTH_200_350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/ggH_PTH_GT350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/xH_PTH_GT350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
            #         ],
            #     )


            # Commands.BasicT2WS( 'suppliedInput/pT_v5_renamed/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.txt' )
            # Commands.BasicT2WS( 'suppliedInput/PTH/hzz4l_comb_13TeV_xs.txt' )
            # Commands.BasicT2WS( 'suppliedInput/fromDavid/PTH_May09/hzz4l_comb_13TeV_xs.txt' )
            # Commands.BasicT2WS( 'suppliedInput/combinedCard_{0}.txt'.format(datestr) )
            # Commands.BasicT2WS( 'suppliedInput/pT_v5_renamed/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.txt' )



    #____________________________________________________________________
    if args.bestfit:

        if args.latest:

            Commands.BasicBestfit(
                # 'workspaces_Jul24/hggMerged_Jul24.root',
                # 'workspaces_Jul21/combinedCard_Jul21.root',
                # 'workspaces_Jul27/combinedCard_Jul26.root',
                # 'workspaces_Jul31/splithggMerged_Jul31.root',
                'workspaces_Aug02/splithggMerged_Jul31.root',
                # 
                # onBatch = True
                extraOptions = [
                    '--minimizerStrategy 2',
                    '-v 2',
                    '-m 125',
                    '--floatOtherPOIs=1',
                    ]
                )


        else:

            # Commands.BasicBestfit(
            #     'workspaces_Aug02/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamedProcesses_Jul21.root',
            #     # onBatch = True
            #     extraOptions = [
            #         '--minimizerStrategy 2',
            #         '-v 2',
            #         '-m 125',
            #         '--floatOtherPOIs=1',
            #         ]
            #     )


            Commands.BasicBestfit(
                # 'workspaces_Aug02/Datacard_13TeV_differential_pT_moriond17_ggHonly_v2_renamedProcesses_Jul31.root',
                # 'workspaces_Aug02/Datacard_13TeV_differential_pT_moriond17_ggHonly_v2.root',
                'workspaces_Aug02/Datacard_13TeV_differential_pT_moriond17_HxOnly.root',
                extraOptions = [
                    '--minimizerStrategy 2',
                    '-v 2',
                    '-m 125',
                    '--floatOtherPOIs=1',
                    ]
                )

            sys.exit()


            # hgg
            # Commands.BasicBestfit(
            #     'workspaces_May15/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.root',
            #     # onBatch = True
            #     )
            
            # hzz
            # Commands.BasicBestfit(
            #     'workspaces_May10/hzz4l_comb_13TeV_xs.root',
            #     # setPOIs = False,
            #     # onBatch = True
            #     )



            # May 15: workspaces with included but faulty _norm

            # extraVars = []
            # nProcesses = len(Commands.ListProcesses( 'suppliedInput/PTH/hzz4l_comb_13TeV_xs.txt' )[0])
            # extraVars.extend( [ 'K1Bin{0}'.format(i) for i in xrange(nProcesses) ] )
            # extraVars.extend( [ 'K2Bin{0}'.format(i) for i in xrange(nProcesses) ] )

            # # hzz
            # Commands.BasicBestfit(
            #     'workspaces_May15/hzz4l_comb_13TeV_xs.root',
            #     setPOIs = False,
            #     # onBatch = True,
            #     extraOptions = (
            #         '--saveSpecifiedFunc {0}'.format( ','.join(extraVars) )
            #         + ' -m 125.09 --floatOtherPOIs=1 --freezeNuisances MH'
            #         )
            #     )



            # Commands.BasicBestfit(
            #     'workspaces_May15/combinedCard_May15.root',
            #     # setPOIs = False,
            #     # onBatch = True
            #     extraOptions = (
            #         ' --floatOtherPOIs=1'
            #         # + ' -m 125.09 --freezeNuisances MH'
            #         )
            #     )



    #____________________________________________________________________
    if args.couplingImportanceHistogram:

        if args.doFastscan:

            datacard = 'workspaces_May30/combinedCard_May15.root'
            Commands.MultiDimCombineTool(
                datacard,
                nPoints       = 10000,
                nPointsPerJob = 10,
                queue         = 'all.q',
                # notOnBatch    = False,
                notOnBatch    = True,
                jobDirectory  = 'Fastscan_couplings_{0}'.format( datestr ),
                extraOptions  = [
                    '-P ct -P cg',
                    '--setPhysicsModelParameters ct=1.0,cg=0.0',
                    '--setPhysicsModelParameterRanges ct=-0.5,4.2:cg=-0.32,0.1',
                    '--fastScan',
                    '--saveSpecifiedFunc {0}'.format( ','.join(Commands.ListSet( datacard, 'yieldParameters' )) ),
                    ]
                )
            TheoryCommands.WriteTH2DToFile(
                glob( 'Fastscan_couplings_{0}/*.root'.format( datestr ) )
                )

        else:

            TheoryCommands.WriteTH2DToFile(
                # glob( 'Scan_couplings_May31_3/*.root' )
                glob( 'Fastscan_couplings_{0}/*.root'.format( datestr ) )
                )



########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'