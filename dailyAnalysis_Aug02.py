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

import combineCommands
import plotCommands

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

def main():

    # ======================================
    # Parser

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument( '--test',                            action='store_true' )

    parser.add_argument( '--newCommand',                      action='store_true' )
    parser.add_argument( '--testPathExtension',               action='store_true' )
    parser.add_argument( '--testPathExtension_bestfit',       action='store_true' )
    parser.add_argument( '--recofitT2WS',                     action='store_true' )
    parser.add_argument( '--recofit_bestfit',                 action='store_true' )

    parser.add_argument( '--unsplit',                 action='store_true' )
    parser.add_argument( '--unsplit_recofit',                 action='store_true' )
    parser.add_argument( '--unsplit_recofit_bestfit',         action='store_true' )



    combineCommands.AppendParserOptions(parser)
    plotCommands.AppendParserOptions(parser)

    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument( '--latest', dest='latest', action='store_true', default=True )
    group.add_argument( '--older',  dest='latest', action='store_false' )

    parser.add_argument( '--ggH',                      action='store_true' )
    parser.add_argument( '--xH',                       action='store_true' )
    parser.add_argument( '--both',                     action='store_true' )


    parser.add_argument( '--merge',                    action='store_true' )
    parser.add_argument( '--merge_recofit',            action='store_true' )
    parser.add_argument( '--merge_bestfit',            action='store_true' )
    parser.add_argument( '--merge_recofit_bestfit',    action='store_true' )

    parser.add_argument( '--doWSrenaming',             action='store_true' )
    
    args = parser.parse_args()

    if not args.ggH and not args.xH:
        args.ggH = True
    if args.both:
        args.ggH = True
        args.xH  = True

    print args
    print ''



    if args.test:
        Commands.TestMode()

    Commands.SetTempJobDir( 'combineOutput_{0}'.format(datestr) )


    ########################################
    # Stuff dealing with combine (datacard merging/combining, t2ws, bestfits, scans, etc.)
    ########################################

    # Moved to separate file
    if args.combineCommands:
        combineCommands.main(args)


    ########################################
    # Result and Test Plotting
    ########################################

    # Moved to separate file
    if args.plotCommands:
        plotCommands.main(args)



    ########################################
    # New commands
    ########################################


    if args.newCommand:
        print 'test'


    if args.testPathExtension:


        if args.ggH:

            ggHCardReparsed, ggh = FullyProcess_ggH()

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

        if args.xH:

            xHCardReparsed, xh = FullyProcess_xH()

            Commands.BasicT2WS(
                xHCardReparsed,
                manualMaps=[
                    '--PO \'map=.*/xH_PTH_0_15:r_xH_PTH_0_15[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/xH_PTH_15_30:r_xH_PTH_15_30[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/xH_PTH_30_45:r_xH_PTH_30_45[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/xH_PTH_45_85:r_xH_PTH_45_85[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/xH_PTH_85_125:r_xH_PTH_85_125[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/xH_PTH_125_200:r_xH_PTH_125_200[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/xH_PTH_200_350:r_xH_PTH_200_350[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/xH_PTH_GT350:r_xH_PTH_GT350[1.0,-1.0,4.0]\'',
                    ],
                )


    if args.recofitT2WS:

        if args.ggH:

            ggHCardReparsed, ggh = FullyProcess_ggH()

            Commands.BasicT2WS(
                ggHCardReparsed,
                outputWS = basename(ggHCardReparsed).replace( '.txt', '_recofit.root' ),
                manualMaps=[
                    '--PO \'map=.*_SigmaMpTTag_0.*/ggH_PTH_0_15:r_ggH_PTH_0_15_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/ggH_PTH_0_15:r_ggH_PTH_0_15_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/ggH_PTH_0_15:r_ggH_PTH_0_15_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/ggH_PTH_15_30:r_ggH_PTH_15_30_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/ggH_PTH_15_30:r_ggH_PTH_15_30_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/ggH_PTH_15_30:r_ggH_PTH_15_30_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/ggH_PTH_30_45:r_ggH_PTH_30_45_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/ggH_PTH_30_45:r_ggH_PTH_30_45_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/ggH_PTH_30_45:r_ggH_PTH_30_45_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/ggH_PTH_45_85:r_ggH_PTH_45_85_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/ggH_PTH_45_85:r_ggH_PTH_45_85_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/ggH_PTH_45_85:r_ggH_PTH_45_85_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/ggH_PTH_85_125:r_ggH_PTH_85_125_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/ggH_PTH_85_125:r_ggH_PTH_85_125_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/ggH_PTH_85_125:r_ggH_PTH_85_125_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/ggH_PTH_125_200:r_ggH_PTH_125_200_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/ggH_PTH_125_200:r_ggH_PTH_125_200_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/ggH_PTH_125_200:r_ggH_PTH_125_200_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/ggH_PTH_200_350:r_ggH_PTH_200_350_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/ggH_PTH_200_350:r_ggH_PTH_200_350_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/ggH_PTH_200_350:r_ggH_PTH_200_350_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/ggH_PTH_GT350:r_ggH_PTH_GT350_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/ggH_PTH_GT350:r_ggH_PTH_GT350_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/ggH_PTH_GT350:r_ggH_PTH_GT350_cat2[1.0,-1.0,20.0]\'',
                    ],
                )

        if args.xH:

            xHCardReparsed, xh = FullyProcess_xH()

            Commands.BasicT2WS(
                xHCardReparsed,
                outputWS = basename(xHCardReparsed).replace( '.txt', '_recofit.root' ),
                manualMaps=[
                    '--PO \'map=.*_SigmaMpTTag_0.*/xH_PTH_0_15:r_xH_PTH_0_15_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/xH_PTH_0_15:r_xH_PTH_0_15_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/xH_PTH_0_15:r_xH_PTH_0_15_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/xH_PTH_15_30:r_xH_PTH_15_30_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/xH_PTH_15_30:r_xH_PTH_15_30_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/xH_PTH_15_30:r_xH_PTH_15_30_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/xH_PTH_30_45:r_xH_PTH_30_45_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/xH_PTH_30_45:r_xH_PTH_30_45_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/xH_PTH_30_45:r_xH_PTH_30_45_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/xH_PTH_45_85:r_xH_PTH_45_85_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/xH_PTH_45_85:r_xH_PTH_45_85_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/xH_PTH_45_85:r_xH_PTH_45_85_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/xH_PTH_85_125:r_xH_PTH_85_125_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/xH_PTH_85_125:r_xH_PTH_85_125_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/xH_PTH_85_125:r_xH_PTH_85_125_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/xH_PTH_125_200:r_xH_PTH_125_200_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/xH_PTH_125_200:r_xH_PTH_125_200_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/xH_PTH_125_200:r_xH_PTH_125_200_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/xH_PTH_200_350:r_xH_PTH_200_350_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/xH_PTH_200_350:r_xH_PTH_200_350_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/xH_PTH_200_350:r_xH_PTH_200_350_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/xH_PTH_GT350:r_xH_PTH_GT350_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/xH_PTH_GT350:r_xH_PTH_GT350_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/xH_PTH_GT350:r_xH_PTH_GT350_cat2[1.0,-1.0,20.0]\'',
                    ],
                )



    if args.testPathExtension_bestfit:

        if args.ggH:
            ws = 'workspaces_Aug03/Datacard_13TeV_differential_pT_moriond17_ggHonly_v2_reparsed.root'
            Commands.BasicBestfit(
                ws,
                onBatch=True,
                batchJobSubDir = 'ggH_only_reparsedTest',
                extraOptions = [
                    '--minimizerStrategy 2',
                    '-v 2',
                    '-m 125',
                    '--floatOtherPOIs=1',
                    '--saveWorkspace',
                    ]
                )

        if args.xH:
            ws = 'workspaces_Aug03/Datacard_13TeV_differential_pT_moriond17_HxOnly_reparsed.root'
            Commands.BasicBestfit(
                ws,
                onBatch=True,
                batchJobSubDir = 'xH_only_reparsedTest',
                extraOptions = [
                    '--minimizerStrategy 2',
                    '-v 2',
                    '-m 125',
                    '--floatOtherPOIs=1',
                    '--saveWorkspace',
                    ]
                )



    if args.recofit_bestfit:

        if args.ggH:
            ws = 'workspaces_Aug03/Datacard_13TeV_differential_pT_moriond17_ggHonly_v2_reparsed_recofit.root'
            Commands.BasicBestfit(
                ws,
                onBatch=True,
                batchJobSubDir = 'ggH_only_recofit',
                extraOptions = [
                    '--minimizerStrategy 2',
                    '-v 2',
                    '-m 125',
                    '--floatOtherPOIs=1',
                    '--saveWorkspace',
                    ]
                )

        if args.xH:
            ws = 'workspaces_Aug03/Datacard_13TeV_differential_pT_moriond17_HxOnly_reparsed_recofit.root'
            Commands.BasicBestfit(
                ws,
                onBatch=True,
                batchJobSubDir = 'xH_only_recofit',
                extraOptions = [
                    '--minimizerStrategy 2',
                    '-v 2',
                    '-m 125',
                    '--floatOtherPOIs=1',
                    '--saveWorkspace',
                    ]
                )


    # ======================================
    # Redo merge

    if args.merge or args.merge_recofit:

        ggHCardReparsed, ggHContainer = FullyProcess_ggH()
        xHCardReparsed, xHContainer   = FullyProcess_xH()

        mergedContainer = MergeHGGWDatacards.MergeCards( ggHContainer, xHContainer )

        mergedCard = join( dirname(ggHCardReparsed), 'Merged_{0}.txt'.format(datestr) )
        MergeHGGWDatacards.ParseDataContainer(
            mergedContainer,
            writeToFile = mergedCard
            )

        if args.merge:
            Commands.BasicT2WS(
                mergedCard,
                manualMaps=[
                    '--PO \'map=.*/.*_PTH_0_15:r_smH_PTH_0_15[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/.*_PTH_15_30:r_smH_PTH_15_30[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/.*_PTH_30_45:r_smH_PTH_30_45[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/.*_PTH_45_85:r_smH_PTH_45_85[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/.*_PTH_85_125:r_smH_PTH_85_125[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/.*_PTH_125_200:r_smH_PTH_125_200[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/.*_PTH_200_350:r_smH_PTH_200_350[1.0,-1.0,4.0]\'',
                    '--PO \'map=.*/.*_PTH_GT350:r_smH_PTH_GT350[1.0,-1.0,4.0]\'',
                    # 
                    # '--PO \'map=.*/ggH_PTH_0_15:r_ggH_PTH_0_15[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/ggH_PTH_15_30:r_ggH_PTH_15_30[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/ggH_PTH_30_45:r_ggH_PTH_30_45[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/ggH_PTH_45_85:r_ggH_PTH_45_85[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/ggH_PTH_85_125:r_ggH_PTH_85_125[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/ggH_PTH_125_200:r_ggH_PTH_125_200[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/ggH_PTH_200_350:r_ggH_PTH_200_350[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/ggH_PTH_GT350:r_ggH_PTH_GT350[1.0,-1.0,4.0]\'',
                    # 
                    # '--PO \'map=.*/xH_PTH_0_15:r_xH_PTH_0_15[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/xH_PTH_15_30:r_xH_PTH_15_30[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/xH_PTH_30_45:r_xH_PTH_30_45[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/xH_PTH_45_85:r_xH_PTH_45_85[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/xH_PTH_85_125:r_xH_PTH_85_125[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/xH_PTH_125_200:r_xH_PTH_125_200[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/xH_PTH_200_350:r_xH_PTH_200_350[1.0,-1.0,4.0]\'',
                    # '--PO \'map=.*/xH_PTH_GT350:r_xH_PTH_GT350[1.0,-1.0,4.0]\'',
                    ],
                )

        if args.merge_recofit:

            Commands.BasicT2WS(
                mergedCard,
                outputWS = basename(mergedCard).replace( '.txt', '_recofit.root' ),
                manualMaps=[
                    '--PO \'map=.*_SigmaMpTTag_0.*/.*_PTH_0_15:r_smH_PTH_0_15_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/.*_PTH_0_15:r_smH_PTH_0_15_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/.*_PTH_0_15:r_smH_PTH_0_15_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/.*_PTH_15_30:r_smH_PTH_15_30_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/.*_PTH_15_30:r_smH_PTH_15_30_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/.*_PTH_15_30:r_smH_PTH_15_30_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/.*_PTH_30_45:r_smH_PTH_30_45_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/.*_PTH_30_45:r_smH_PTH_30_45_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/.*_PTH_30_45:r_smH_PTH_30_45_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/.*_PTH_45_85:r_smH_PTH_45_85_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/.*_PTH_45_85:r_smH_PTH_45_85_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/.*_PTH_45_85:r_smH_PTH_45_85_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/.*_PTH_85_125:r_smH_PTH_85_125_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/.*_PTH_85_125:r_smH_PTH_85_125_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/.*_PTH_85_125:r_smH_PTH_85_125_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/.*_PTH_125_200:r_smH_PTH_125_200_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/.*_PTH_125_200:r_smH_PTH_125_200_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/.*_PTH_125_200:r_smH_PTH_125_200_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/.*_PTH_200_350:r_smH_PTH_200_350_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/.*_PTH_200_350:r_smH_PTH_200_350_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/.*_PTH_200_350:r_smH_PTH_200_350_cat2[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_0.*/.*_PTH_GT350:r_smH_PTH_GT350_cat0[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_1.*/.*_PTH_GT350:r_smH_PTH_GT350_cat1[1.0,-1.0,20.0]\'',
                    '--PO \'map=.*_SigmaMpTTag_2.*/.*_PTH_GT350:r_smH_PTH_GT350_cat2[1.0,-1.0,20.0]\'',
                    ],
                )



    if args.unsplit:

        unsplitCard = 'suppliedInput/fromVittorio/pT_unsplit_Mar21/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamedProcesses_Jul24.txt'

        Commands.BasicT2WS(
            unsplitCard,
            manualMaps=[
                '--PO \'map=.*/.*_PTH_0_15:r_smH_PTH_0_15[1.0,-1.0,20.0]\'',
                '--PO \'map=.*/.*_PTH_15_30:r_smH_PTH_15_30[1.0,-1.0,20.0]\'',
                '--PO \'map=.*/.*_PTH_30_45:r_smH_PTH_30_45[1.0,-1.0,20.0]\'',
                '--PO \'map=.*/.*_PTH_45_85:r_smH_PTH_45_85[1.0,-1.0,20.0]\'',
                '--PO \'map=.*/.*_PTH_85_125:r_smH_PTH_85_125[1.0,-1.0,20.0]\'',
                '--PO \'map=.*/.*_PTH_125_200:r_smH_PTH_125_200[1.0,-1.0,20.0]\'',
                '--PO \'map=.*/.*_PTH_200_350:r_smH_PTH_200_350[1.0,-1.0,20.0]\'',
                '--PO \'map=.*/.*_PTH_GT350:r_smH_PTH_GT350[1.0,-1.0,20.0]\'',
                ],
            )



    if args.unsplit_recofit:

        unsplitCard = 'suppliedInput/fromVittorio/pT_unsplit_Mar21/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamedProcesses_Jul24.txt'

        Commands.BasicT2WS(
            unsplitCard,
            outputWS = basename(unsplitCard).replace('.txt','_recofit.root'),
            manualMaps=[
                '--PO \'map=.*SigmaMpTTag_0.*/.*_PTH_0_15:r_smH_PTH_0_15_cat0[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_1.*/.*_PTH_0_15:r_smH_PTH_0_15_cat1[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_2.*/.*_PTH_0_15:r_smH_PTH_0_15_cat2[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_0.*/.*_PTH_15_30:r_smH_PTH_15_30_cat0[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_1.*/.*_PTH_15_30:r_smH_PTH_15_30_cat1[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_2.*/.*_PTH_15_30:r_smH_PTH_15_30_cat2[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_0.*/.*_PTH_30_45:r_smH_PTH_30_45_cat0[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_1.*/.*_PTH_30_45:r_smH_PTH_30_45_cat1[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_2.*/.*_PTH_30_45:r_smH_PTH_30_45_cat2[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_0.*/.*_PTH_45_85:r_smH_PTH_45_85_cat0[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_1.*/.*_PTH_45_85:r_smH_PTH_45_85_cat1[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_2.*/.*_PTH_45_85:r_smH_PTH_45_85_cat2[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_0.*/.*_PTH_85_125:r_smH_PTH_85_125_cat0[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_1.*/.*_PTH_85_125:r_smH_PTH_85_125_cat1[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_2.*/.*_PTH_85_125:r_smH_PTH_85_125_cat2[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_0.*/.*_PTH_125_200:r_smH_PTH_125_200_cat0[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_1.*/.*_PTH_125_200:r_smH_PTH_125_200_cat1[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_2.*/.*_PTH_125_200:r_smH_PTH_125_200_cat2[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_0.*/.*_PTH_200_350:r_smH_PTH_200_350_cat0[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_1.*/.*_PTH_200_350:r_smH_PTH_200_350_cat1[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_2.*/.*_PTH_200_350:r_smH_PTH_200_350_cat2[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_0.*/.*_PTH_GT350:r_smH_PTH_GT350_cat0[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_1.*/.*_PTH_GT350:r_smH_PTH_GT350_cat1[1.0,-1.0,20.0]\'',
                '--PO \'map=.*SigmaMpTTag_2.*/.*_PTH_GT350:r_smH_PTH_GT350_cat2[1.0,-1.0,20.0]\'',
                ],
            )


    if args.merge_bestfit:

        Commands.BasicBestfit(
            'workspaces_Aug03/Merged_Aug03.root',
            onBatch=True,
            batchJobSubDir = 'merged',
            extraOptions = [
                '--minimizerStrategy 2',
                '-v 2',
                '-m 125',
                '--floatOtherPOIs=1',
                '--saveWorkspace',
                ]
            )


    if args.merge_recofit_bestfit:

        Commands.BasicBestfit(
            'workspaces_Aug03/Merged_Aug03_recofit.root',
            onBatch=True,
            batchJobSubDir = 'merged_recofit',
            extraOptions = [
                '--minimizerStrategy 2',
                '-v 2',
                '-m 125',
                '--floatOtherPOIs=1',
                '--saveWorkspace',
                ]
            )

    if args.unsplit_recofit_bestfit:

        Commands.BasicBestfit(
            'workspaces_Aug03/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamedProcesses_Jul24_recofit.root',
            onBatch=True,
            batchJobSubDir = 'unsplit_recofit',
            extraOptions = [
                '--minimizerStrategy 2',
                '-v 2',
                '-m 125',
                '--floatOtherPOIs=1',
                '--saveWorkspace',
                ]
            )




    if args.doWSrenaming:

        if args.ggH:
            MergeHGGWDatacards.RenameProductionModeHgg(
                'ggH',
                'suppliedInput/fromVittorio/pT_ggHonly_Jun26/CMS-HGG_sigfit_differential_pT_moriond17_ggHonly_v2.root'
                )

        if args.xH:
            MergeHGGWDatacards.RenameProductionModeHgg(
                'xH',
                'suppliedInput/fromVittorio/differential_pT_moriond17_HxOnly_July19/CMS-HGG_sigfit_differential_pT_moriond17_HxOnly.root'
                )




########################################
# Helper functions
########################################

def FullyProcess_ggH():

    ggHCard = 'suppliedInput/fromVittorio/pT_ggHonly_Jun26/Datacard_13TeV_differential_pT_moriond17_ggHonly_v2.txt'

    ggH_DCrenamed = MergeHGGWDatacards.RenameProcesses(
        'ggH',
        ggHCard,
        # outdatacard='auto',
        renameOutsideAcceptance=False,
        globalReplace = [
            (
                'CMS-HGG_sigfit_differential_pT_moriond17_ggHonly_v2.root',
                'CMS-HGG_sigfit_differential_pT_moriond17_ggHonly_v2_FullyRenamed_Aug03.root'
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

    ggh = MergeHGGWDatacards.GetDatacardContainer( ggH_DCrenamed )
    ggh.process = 'ggH'

    MergeHGGWDatacards.ExtendPathOfRootFiles( ggh, dirname(ggHCard) )


    ggHCardReparsed  = join( dirname(ggHCard), '..', basename(ggHCard).replace( '.txt', '_reparsed.txt' ) )
    MergeHGGWDatacards.ParseDataContainer(
        ggh,
        writeToFile = ggHCardReparsed
        )

    return ggHCardReparsed, ggh


def FullyProcess_xH():

    xHCard = 'suppliedInput/fromVittorio/differential_pT_moriond17_HxOnly_July19/Datacard_13TeV_differential_pT_moriond17_HxOnly.txt'

    xH_DCrenamed = MergeHGGWDatacards.RenameProcesses(
        'xH',
        xHCard,
        # outdatacard='auto',
        renameOutsideAcceptance=False,
        globalReplace = [
            (
                'CMS-HGG_sigfit_differential_pT_moriond17_HxOnly.root',
                'CMS-HGG_sigfit_differential_pT_moriond17_HxOnly_FullyRenamed_Aug03.root'
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


    xh = MergeHGGWDatacards.GetDatacardContainer( xH_DCrenamed )
    xh.process = 'xH'

    MergeHGGWDatacards.ExtendPathOfRootFiles( xh, dirname(xHCard) )


    xHCardReparsed  = join( dirname(xHCard), '..', basename(xHCard).replace( '.txt', '_reparsed.txt' ) )
    MergeHGGWDatacards.ParseDataContainer(
        xh,
        writeToFile = xHCardReparsed
        )

    return xHCardReparsed, xh


########################################
# End of Main
########################################
if __name__ == "__main__":
    main()
