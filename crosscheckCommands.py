#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, itertools, operator, re, argparse, sys, random, numpy
from math import isnan, isinf, sqrt
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
import TheoryFileInterface
import OutputInterface
from Container import Container
from Parametrization import Parametrization, WSParametrization
import PlotCommands
import LatestBinning

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

    parser.add_argument( '--crosscheckCommands', default=False )
    class CustomAction(argparse.Action):
        def __init__(self, option_strings, dest, **kwargs):
            # Only 1 argument allowed (this is basically the "store_true" action)
            super(CustomAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            # Make sure the parser knows a combineCommand was entered
            setattr( namespace, 'crosscheckCommands', True )
            setattr( namespace, self.dest, True )
            # super(FooAction, self).__call__( parser, namespace, values, option_string=None )
            # setattr( namespace, self.dest, values )

    parser.add_argument( '--normalization_Top',                  action=CustomAction )
    parser.add_argument( '--normalization_Yukawa',                  action=CustomAction )


########################################
# Methods
########################################    

def main( args ):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )

    #____________________________________________________________________
    if args.normalization_Top:

        LOGSCALE = True

        from LatestBinning import obs_pth_combWithHbbBinning as obs_pth
        obs_pth.YR4_totalXS = LatestBinning.YR4_ggF_n3lo
        SMXS = obs_pth.crosssection_over_binwidth()

        H_SMXS = obs_pth.BasicHistogram()
        H_SMXS.SetLineColor(1)

        xMin = obs_pth.binning[0]
        if obs_pth.lastBinIsOverflow:
            xMax = obs_pth.binning[-2] + (obs_pth.binning[-2]-obs_pth.binning[-3])
        else:
            xMax = obs_pth.binning[-1]

        if LOGSCALE:
            yMinAbs = min([ y for y in SMXS if y >= 0.0 ])
            yMaxAbs = max(SMXS)
            yMin = 0.01 * yMinAbs
            yMax = 2.0 * yMaxAbs
        else:
            yMinAbs = min(SMXS)
            yMaxAbs = max(SMXS)
            yMin = yMinAbs - 0.1*(yMaxAbs-yMinAbs)
            yMax = yMaxAbs + 0.1*(yMaxAbs-yMinAbs)

        c.Clear()
        SetCMargins()

        base = GetPlotBase(
            xMin = xMin,
            xMax = xMax,
            yMin = yMin,
            yMax = yMax,
            xTitle = 'p_{T} (GeV)',
            # yTitle = '#Delta NLL',
            yTitle = '#sigma (pb/GeV)',
            )
        base.Draw('P')

        H_SMXS.Draw('HISTSAME')


        # ======================================
        # Load theory predictions on top

        SM_Top = TheoryFileInterface.FileFinder(
            ct=1, cg=0, cb=1, muR=1, muF=1, Q=1,
            expectOneFile=True, loadImmediately=True,
            # directory=LatestPaths.derivedTheoryFiles_TopHighPt
            directory=LatestPaths.derivedTheoryFiles_Top
            )
        SM_Top_normalization = BasicIntegral( SM_Top.binBoundaries, SM_Top.crosssection, lastBinIsOverflow=False )

        H_Top = BasicHistogramFromContainer( SM_Top, lastBinIsOverflow=False )
        # H_Top.Draw('HISTSAME')



        SM_Top_Rebinned = deepcopy(SM_Top)
        TheoryCommands.RebinDerivedTheoryContainer( SM_Top_Rebinned, obs_pth.binning, lastBinIsOverflow=False )

        H_Top_Rebinned = BasicHistogramFromContainer( SM_Top_Rebinned, lastBinIsOverflow=False )
        H_Top_Rebinned.SetLineColor(9)
        H_Top_Rebinned.Draw('HISTSAME')

        # Also draw a line for the last bin that is not /GeV
        lastBinIntegratedLine = GetIntegratedLineForLastBin( SM_Top_Rebinned )
        lastBinIntegratedLine.SetLineColor(9)
        lastBinIntegratedLine.Draw()


        # ======================================
        # Now rescale and do again

        SM_Top_Rescaled = deepcopy( SM_Top )
        SM_Top_Rescaled.crosssection = [ xs * ( LatestBinning.YR4_ggF_n3lo / SM_Top_normalization ) for xs in SM_Top_Rescaled.crosssection ]


        H_Top_Rescaled = BasicHistogramFromContainer( SM_Top_Rescaled, lastBinIsOverflow=False )
        H_Top_Rescaled.SetLineColor(2)
        # H_Top_Rescaled.Draw('HISTSAME')


        SM_Top_Rescaled_Rebinned = deepcopy(SM_Top_Rescaled)
        TheoryCommands.RebinDerivedTheoryContainer( SM_Top_Rescaled_Rebinned, obs_pth.binning, lastBinIsOverflow=False )

        H_Top_Rescaled_Rebinned = BasicHistogramFromContainer( SM_Top_Rescaled_Rebinned, lastBinIsOverflow=False )
        H_Top_Rescaled_Rebinned.SetLineColor(46)
        H_Top_Rescaled_Rebinned.Draw('HISTSAME')

        # Also draw a line for the last bin that is not /GeV
        lastBinIntegratedLine_Rescaled = GetIntegratedLineForLastBin( SM_Top_Rescaled_Rebinned )
        lastBinIntegratedLine_Rescaled.SetLineColor(46)
        lastBinIntegratedLine_Rescaled.Draw()



        # ======================================

        print '\nSM:'
        BasicIntegral( obs_pth.binning, obs_pth.crosssection_over_binwidth(), lastBinIsOverflow=True, verbose=True )
        # print '\nTop:'
        # BasicIntegral( SM_Top.binBoundaries, SM_Top.crosssection, lastBinIsOverflow=False, verbose=True )
        print '\nTop_Rebinned:'
        BasicIntegral( SM_Top_Rebinned.binBoundaries, SM_Top_Rebinned.crosssection, lastBinIsOverflow=False, verbose=True )
        # print '\nTop_Rescaled:'
        # BasicIntegral( SM_Top_Rescaled.binBoundaries, SM_Top_Rescaled.crosssection, lastBinIsOverflow=False, verbose=True )
        print '\nTop_Rescaled_Rebinned:'
        BasicIntegral( SM_Top_Rescaled_Rebinned.binBoundaries, SM_Top_Rescaled_Rebinned.crosssection, lastBinIsOverflow=False, verbose=True )


        print '\n'
        print 'YR4:          {1:.4f}  (ggH+xH={0:.4f}, xH={2:.4f})'.format( LatestBinning.YR4_totalXS, LatestBinning.YR4_ggF_n3lo, LatestBinning.YR4_xH )
        print 'SM:           {0:.4f}'.format( BasicIntegral( obs_pth.binning, obs_pth.crosssection_over_binwidth(), lastBinIsOverflow=True ) )
        print 'Top:          {0:.4f}'.format( BasicIntegral( SM_Top.binBoundaries, SM_Top.crosssection, lastBinIsOverflow=False ) )
        print 'Top_Rescaled: {0:.4f}'.format( BasicIntegral( SM_Top_Rescaled.binBoundaries, SM_Top_Rescaled.crosssection, lastBinIsOverflow=False ) )
        # print 'Top rebinned: ', BasicIntegral( SM_Top_Rebinned.binBoundaries, SM_Top_Rebinned.crosssection )

        if LOGSCALE: c.SetLogy()
        SaveC( 'normalizationCrosscheck_Top' )


    #____________________________________________________________________
    if args.normalization_Yukawa:

        LOGSCALE = True

        from LatestBinning import obs_pth_combWithHbbBinning as obs_pth
        obs_pth.YR4_totalXS = LatestBinning.YR4_ggF_n3lo
        obs_pth.keepOnlyFirstNBins(5)

        H_SMXS = obs_pth.BasicHistogram()
        H_SMXS.SetLineColor(1)

        xMin = obs_pth.binning[0]
        if obs_pth.lastBinIsOverflow:
            xMax = obs_pth.binning[-2] + (obs_pth.binning[-2]-obs_pth.binning[-3])
        else:
            xMax = obs_pth.binning[-1]

        if LOGSCALE:
            yMinAbs = min([ y for y in obs_pth.crosssection_over_binwidth() if y >= 0.0 ])
            yMaxAbs = max(obs_pth.crosssection_over_binwidth())
            yMin = 0.01 * yMinAbs
            yMax = 2.0 * yMaxAbs
        else:
            yMinAbs = min(obs_pth.crosssection_over_binwidth())
            yMaxAbs = max(obs_pth.crosssection_over_binwidth())
            yMin = yMinAbs - 0.1*(yMaxAbs-yMinAbs)
            yMax = yMaxAbs + 0.1*(yMaxAbs-yMinAbs)

        c.Clear()
        SetCMargins()

        base = GetPlotBase(
            xMin = xMin,
            xMax = xMax,
            yMin = yMin,
            yMax = yMax,
            xTitle = 'p_{T} (GeV)',
            # yTitle = '#Delta NLL',
            yTitle = '#sigma (pb/GeV)',
            )
        base.Draw('P')

        H_SMXS.Draw('HISTSAME')


        #____________________________________________________________________

        SM_Yukawa = TheoryFileInterface.FileFinder(
            kappab=1, kappac=1, muR=1, muF=1, Q=1,
            expectOneFile=True, loadImmediately=True,
            directory=LatestPaths.derivedTheoryFiles_YukawaSummed
            )
        SM_Yukawa_normalization = BasicIntegral( SM_Yukawa.binBoundaries, SM_Yukawa.crosssection, lastBinIsOverflow=False )

        H_Yukawa = BasicHistogramFromContainer( SM_Yukawa, lastBinIsOverflow=False )
        H_Yukawa.Draw('HISTSAME')



        SM_Yukawa_Rebinned = deepcopy(SM_Yukawa)
        TheoryCommands.RebinDerivedTheoryContainer( SM_Yukawa_Rebinned, obs_pth.binning, lastBinIsOverflow=False )

        H_Yukawa_Rebinned = BasicHistogramFromContainer( SM_Yukawa_Rebinned, lastBinIsOverflow=False )
        H_Yukawa_Rebinned.SetLineColor(9)
        H_Yukawa_Rebinned.Draw('HISTSAME')


        print '\nSM:'
        BasicIntegral( obs_pth.binning, obs_pth.crosssection_over_binwidth(), lastBinIsOverflow=False, verbose=True )
        print '\nYukawa_Rebinned:'
        BasicIntegral( SM_Yukawa_Rebinned.binBoundaries, SM_Yukawa_Rebinned.crosssection, lastBinIsOverflow=False, verbose=True )
        # print '\nYukawa_Rescaled_Rebinned:'
        # BasicIntegral( SM_Yukawa_Rescaled_Rebinned.binBoundaries, SM_Yukawa_Rescaled_Rebinned.crosssection, lastBinIsOverflow=False, verbose=True )



        if LOGSCALE: c.SetLogy()
        SaveC( 'normalizationCrosscheck_Yukawa' )




#____________________________________________________________________
def BasicHistogramFromContainer( container, lastBinIsOverflow=False ):

    nBins = len(container.crosssection)

    H = ROOT.TH1F(
        TheoryCommands.GetUniqueRootName(), TheoryCommands.GetUniqueRootName(),
        nBins, array( 'f', container.binBoundaries )
        )
    ROOT.SetOwnership( H, False )

    for iBin in xrange(nBins):
        H.SetBinContent( iBin+1, container.crosssection[iBin] )

    # if lastBinIsOverflow:
    #     H.SetBinContent( nBins, container.crosssection[-1] / ( container.binBoundaries[-2] - container.binBoundaries[-3] ) )

    H.SetLineWidth(2)
    H.SetLineColor(4)

    return H
    

#____________________________________________________________________
def GetIntegratedLineForLastBin( container ):
    left  = container.binBoundaries[-2]
    right = container.binBoundaries[-1]
    lastBinIntegratedLine = ROOT.TLine(
        left,  container.crosssection[-1] * (right-left),
        right, container.crosssection[-1] * (right-left),
        )
    lastBinIntegratedLine.SetLineWidth(2)
    lastBinIntegratedLine.SetLineStyle(2)
    ROOT.SetOwnership( lastBinIntegratedLine, False )
    return lastBinIntegratedLine
 

#____________________________________________________________________
def BasicIntegral( binning, xss_per_binWidth, lastBinIsOverflow=False, verbose=False ):

    if not len(binning)-1 == len(xss_per_binWidth):
        Commands.ThrowError( 'len(binning) = {0}, len(xss) = {1}; Cannot calculate integral.'.format( len(binning), len(xss_per_binWidth) ) )

    nBins = len(xss_per_binWidth)

    s = 0.0
    # if verbose: print ''
    for iBin in xrange(nBins-1):
        s += xss_per_binWidth[iBin] * ( binning[iBin+1] - binning[iBin] )
        if verbose: print 'Bin {0:3} | xs = {1} | left = {2} | right = {3} | integral = {4}'.format(
            iBin,
            formatNumber( xss_per_binWidth[iBin] ),
            formatNumber( binning[iBin] ),
            formatNumber( binning[iBin+1] ),
            formatNumber( xss_per_binWidth[iBin] * ( binning[iBin+1] - binning[iBin]  ) ),
            )

    if lastBinIsOverflow:
        if verbose: print 'Last Bin| xs = {0} | unintegrated'.format( formatNumber(xss_per_binWidth[-1]) )
        s += xss_per_binWidth[-1]
    else:
        if verbose: print 'Bin {0:3} | xs = {1} | left = {2} | right = {3} | integral = {4}'.format(
            nBins-1,
            formatNumber( xss_per_binWidth[-1] ),
            formatNumber( binning[-2] ),
            formatNumber( binning[-1] ),
            formatNumber( xss_per_binWidth[-1] * ( binning[-1] - binning[-2] ) ),
            )
        s += xss_per_binWidth[-1] * ( binning[-1] - binning[-2] )

    if verbose: print 'Total xs: {0}'.format(formatNumber(s))
    return s

#____________________________________________________________________
def formatNumber( number, width=10 ):
    s = '{0:.4f}'.format(number)
    s = '{0:{width}}'.format( s, width=width )
    return s


########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'