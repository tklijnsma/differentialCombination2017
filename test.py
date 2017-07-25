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
    parser.add_argument( '--makeStewartTackmannDatacard',     action='store_true' )
    parser.add_argument( '--createDerivedTheoryFiles',        action='store_true' )
    parser.add_argument( '--createDerivedTheoryFiles_Yukawa', action='store_true' )
    parser.add_argument( '--CorrelationMatrices',             action='store_true' )

    combineCommands.AppendParserOptions(parser)
    plotCommands.AppendParserOptions(parser)

    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument( '--latest', dest='latest', action='store_true', default=True )
    group.add_argument( '--older',  dest='latest', action='store_false' )

    args = parser.parse_args()

    print args



    if args.test:
        Commands.TestMode()



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
    # Stuff dealing with theory spectra (creation of derivedTheoryFiles, rebinning, correlation matrices, etc.)
    ########################################


    if args.CorrelationMatrices:

        variationFiles = glob( 'suppliedInput/fromAgnieszka/ScaleVarNNLO_Jul17/*.top' )

        # ======================================
        # Correlation matrix in theory binning

        variations = [ CorrelationMatrices.ReadVariationFile( variationFile, fromAgnieszka=True ) for variationFile in variationFiles ]
        CorrelationMatrices.GetCorrelationMatrix(
            variations,
            makeScatterPlots          = False,
            makeCorrelationMatrixPlot = True,
            outname                   = 'corrMat_theory',
            verbose                   = True,
            )

        # ======================================
        # Correlation matrix in exp binning
        
        print '[fixme] Exp bin boundaries hardcoded'
        expBinBoundaries    = [ 0., 15., 30., 45., 85., 125., 200., 350., 800. ]
        expBinCenters       = [ 0.5*(expBinBoundaries[iBin]+expBinBoundaries[iBin+1]) for iBin in xrange(len(expBinBoundaries)-1) ]

        variations_expbinning = []
        for variation in variations:
            theoryBinCenters, theoryBinBoundaries, theoryBinWidths = TheoryCommands.BinningHeuristic( variation.binCenters, manualSwitchAt50=True )
            expBinValues = TheoryCommands.MapFineToCoarse(
                theoryBinBoundaries = theoryBinBoundaries,
                theoryBinValues     = variation.binValues,
                expBinBoundaries    = expBinBoundaries,
                lastBinIsOverflow   = True,
                )

            variation_exp = deepcopy( variation )
            variation_exp.binCenters = expBinCenters
            variation_exp.binValues  = expBinValues
            variation_exp.binBoundaries = expBinBoundaries

            variations_expbinning.append( variation_exp )


        CorrelationMatrices.GetCorrelationMatrix(
            variations_expbinning,
            # makeScatterPlots          = False,
            makeScatterPlots          = True,
            makeCorrelationMatrixPlot = True,
            outname                   = 'corrMat_exp',
            verbose                   = True,
            )

        CorrelationMatrices.PlotVariationSpectra( variations_expbinning, 'exp' )
        CorrelationMatrices.PlotVariationSpectra( variations, 'theory' )

        CorrelationMatrices.PlotRelativeUncertainty( variations_expbinning, 'exp' )
        CorrelationMatrices.PlotRelativeUncertainty( variations, 'theory' )


        # ======================================
        # Read HqT variations; These are split up over 2 files and need to be merged

        LOHqTvariationFiles_lowPt = glob( 'HqT/infiles_Jul19/LO*minPt_1_*.out' )
        LOHqTvariations_lowPt = [ CorrelationMatrices.ReadVariationFile( variationFile, fromHqT=True ) for variationFile in LOHqTvariationFiles_lowPt ]
        LOHqTvariationFiles_highPt = glob( 'HqT/infiles_Jul19/LO*minPt_50_*.out' )
        LOHqTvariations_highPt = [ CorrelationMatrices.ReadVariationFile( variationFile, fromHqT=True ) for variationFile in LOHqTvariationFiles_highPt ]
        LOHqTvariations = CorrelationMatrices.MergeLowHighPtFilesFromHqT( LOHqTvariations_lowPt, LOHqTvariations_highPt )

        CorrelationMatrices.PlotVariationSpectra( LOHqTvariations, 'hqtLO' )
        CorrelationMatrices.PlotRelativeUncertainty( LOHqTvariations, 'hqtLO' )

        # [July 19 14:40] No working NLO files yet
        NLOHqTvariationFiles_lowPt = glob( 'HqT/infiles_Jul19/NLO*minPt_1_*.out' )
        NLOHqTvariations_lowPt = [ CorrelationMatrices.ReadVariationFile( variationFile, fromHqT=True ) for variationFile in NLOHqTvariationFiles_lowPt ]
        NLOHqTvariationFiles_highPt = glob( 'HqT/infiles_Jul19/NLO*minPt_50_*.out' )
        NLOHqTvariations_highPt = [ CorrelationMatrices.ReadVariationFile( variationFile, fromHqT=True ) for variationFile in NLOHqTvariationFiles_highPt ]
        NLOHqTvariations = CorrelationMatrices.MergeLowHighPtFilesFromHqT( NLOHqTvariations_lowPt, NLOHqTvariations_highPt )

        CorrelationMatrices.PlotVariationSpectra( NLOHqTvariations, 'hqtNLO' )
        CorrelationMatrices.PlotRelativeUncertainty( NLOHqTvariations, 'hqtNLO' )

        CorrelationMatrices.PlotWithEnvelop(
            ( 'hqtLO', LOHqTvariations ),
            ( 'hqtNLO', NLOHqTvariations ),
            ( 'exp', variations_expbinning ),
            ( 'hresNNLO', variations ),
            ptMax = 150.
            )


    if args.makeStewartTackmannDatacard:

        container = TheoryCommands.ReadDerivedTheoryFile( 'derivedTheoryFiles_Jun22/SM_NNLO.txt', returnContainer=True )
        
        covMat = TheoryCommands.GetStewartTackmannCovarianceMatrix( container )

        # No longer works! Was anyway incorrectly implemented
        # TheoryCommands.AddCovarianceMatrixAsNuisanceParameters(
        #     'suppliedInput/combinedCard_May15.txt',
        #     covMat
        #     )


    if args.createDerivedTheoryFiles:
        TheoryCommands.CreateDerivedTheoryFiles( pattern=r'ct_[mp\d]+_cg_[mp\d]+' )

    if args.createDerivedTheoryFiles_Yukawa:
        TheoryCommands.CreateDerivedTheoryFiles_Yukawa()




########################################
# End of Main
########################################
if __name__ == "__main__":
    main()