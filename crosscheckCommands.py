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
    parser.add_argument( '--normalization_Top_reimpl',           action=CustomAction )
    parser.add_argument( '--normalization_Yukawa',               action=CustomAction )
    parser.add_argument( '--normalization_Yukawa_reimpl',        action=CustomAction )

########################################
# Methods
########################################    

class Spectrum(object):
    """docstring for Spectrum"""
    def __init__(self, name, values, isSM, title=None, color=4, line_style=1):
        self.values = values
        self.name = name
        self.isSM = isSM
        self.color = color
        self.line_style = line_style
        if title is None:
            self.title = name
        else:
            self.title = title

class NormalizationCrossCheck(object):
    """docstring for NormalizationCrossCheck"""
    colorCycle = itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )
 
    def __init__(self):
        self.logscale = True
        self.last_bin_is_overflow = False
        self.spectra = []

    def add_spectrum(self, name, values, isSM=False, color=None, line_style=1):
        if len(values) < self.n_bins:
            raise ValueError(
                'The passed spectrum \'{0}\' has {1} bins, but the number of passed bins is {2}'.format(
                    name, len(values), self.n_bins) +
                '\n    values: {0}'.format(values) +
                '\n    bin_boundaries: {0}'.format(self.bin_boundaries)
                )
        elif len(values) > self.n_bins:
            Commands.warning(
                'The passed spectrum \'{0}\' has {1} bins, but the number of passed bins is {2}'.format(
                    name, len(values), self.n_bins) +
                '\n    Will keep only the first {0} bins'.format(self.n_bins)
                )
            values = values[:self.n_bins]

        if color is None:
            color = next(self.colorCycle)
        self.spectra.append(Spectrum(name, values, isSM, color=color, line_style=line_style))

    def set_bin_boundaries(self, bin_boundaries, add_overflow=False):
        self.bin_boundaries = bin_boundaries
        if add_overflow:
            self.n_bins = len(self.bin_boundaries)
            self.last_bin_is_overflow = True
            self.bin_boundaries.append(10000.)
        else:
            self.n_bins = len(self.bin_boundaries)-1
    
    def get_SM(self):
        return [s for s in self.spectra if s.isSM][0]

    def plot_spectra(self, tag=None):
        SM = self.get_SM()

        x_min = self.bin_boundaries[0]
        if self.last_bin_is_overflow:
            x_max = self.bin_boundaries[-2] + (self.bin_boundaries[-2]-self.bin_boundaries[-3])
        else:
            x_max = self.bin_boundaries[-1]

        if self.logscale:
            y_minAbs = min([ y for y in SM.values if y >= 0.0 ])
            y_maxAbs = max(SM.values)
            y_min = 0.01 * y_minAbs
            y_max = 2.0 * y_maxAbs
        else:
            y_minAbs = min(SM.values)
            y_maxAbs = max(SM.values)
            y_min = y_minAbs - 0.1*(y_maxAbs-y_minAbs)
            y_max = y_maxAbs + 0.1*(y_maxAbs-y_minAbs)

        c.Clear()
        set_cmargins()
        base = get_plot_base(
            xMin = x_min,
            xMax = x_max,
            yMin = y_min,
            yMax = y_max,
            xTitle = 'p_{T} (GeV)',
            # yTitle = '#Delta NLL',
            yTitle = '#sigma (pb/GeV)',
            )
        base.Draw('P')

        leg = ROOT.TLegend(
            1 - c.GetRightMargin() - 0.3,
            1 - c.GetTopMargin()   - 0.25,
            1 - c.GetRightMargin(),
            1 - c.GetTopMargin() 
            )
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)

        for spectrum in self.spectra:
            H = self.spectrum_to_hist(spectrum)
            H.Draw('HISTSAME')
            leg.AddEntry( H.GetName(), spectrum.title, 'l' )

        leg.Draw()

        if self.logscale: c.SetLogy()

        outname = 'normalizationCrosscheck_reimpl'
        if not(tag is None):
            outname += '_' + tag
        save_c(outname)


    def spectrum_to_hist(self, s):
        H = ROOT.TH1F(
            s.name+'_'+TheoryCommands.get_unique_root_name(), s.title,
            self.n_bins, array('f', self.bin_boundaries)
            )
        ROOT.SetOwnership( H, False )
        for i_bin in xrange(self.n_bins):
            H.SetBinContent(i_bin+1, s.values[i_bin])
        H.SetLineWidth(2)
        H.SetLineColor(s.color)
        H.SetLineStyle(s.line_style)
        return H




def main( args ):
    TheoryCommands.set_plot_dir( 'plots_{0}'.format(datestr) )


    #____________________________________________________________________
    if args.normalization_Top_reimpl:

        normalizationCrossCheck = NormalizationCrossCheck()

        from LatestBinning import obs_pth_combWithHbbBinning as obs_pth
        obs_pth.YR4_totalXS = LatestBinning.YR4_ggF_n3lo

        Commands.warning('Merging some SM bin boundaries; fix when hgg comes in (BAD FIX!)')
        obs_pth.binning[obs_pth.binning.index(85.)] = 80.
        obs_pth.binning[obs_pth.binning.index(125.)] = 120.
        obs_pth = obs_pth.merge_bins([ 0, 1, [2,3], [4,5], 6, 7, 8 ])

        normalizationCrossCheck.set_bin_boundaries(obs_pth.binning[:-1], add_overflow=True)

        # Add SM
        SMXSs = obs_pth.crosssection_over_binwidth()
        normalizationCrossCheck.add_spectrum('YR4, shape Vitt.', SMXSs, isSM=True, color=4)

        # Add raw from Agnieszka
        SM_Top = TheoryFileInterface.file_finder(
            ct=1, cg=0, cb=1, muR=1, muF=1, Q=1,
            expectOneFile=True, loadImmediately=True,
            # directory=LatestPaths.derivedTheoryFiles_TopHighPt
            directory=LatestPaths.derivedTheoryFiles_Top
            )
        xs_Top = TheoryCommands.rebin(SM_Top.binBoundaries, SM_Top.crosssection, obs_pth.binning, lastBinIsOverflow=True)
        normalizationCrossCheck.add_spectrum('Agnieszka', xs_Top, color=2)

        # Add SM from the ws (what is actually used)
        ws = 'out/workspaces_Jan29/combinedCard_newBins_hzz_hbb_ggHxH_Jan24_CouplingModel_TopHighPt_noTheoryUncertainties_hzz_hbb_newBins.root'
        with Commands.OpenRootFile(ws) as wsFp:
            w = wsFp.Get('w')
            yieldParameters = Commands.list_set(w, 'all_ggH_yieldParameters')
            yieldParameters.sort(key=Commands.rangeSorter)
            reweightors = [ w.function(y).getVal() for y in yieldParameters ]
            SMXS_ws = [ w.var(y).getVal() for y in Commands.list_set(w, 'SMXS') ]
        binWidths = [ right-left for left, right in zip(normalizationCrossCheck.bin_boundaries[:-1], normalizationCrossCheck.bin_boundaries[1:]) ]
        binWidths[-1] = 1.0 # Don't divide by bin width in overflow bin
        ws_xs = [ reweightor * SMXS / binWidth for reweightor, SMXS, binWidth in zip(reweightors, SMXS_ws, binWidths) ]

        normalizationCrossCheck.add_spectrum('Agn. reweighted', ws_xs, line_style=2)

        normalizationCrossCheck.plot_spectra('Top')


    #____________________________________________________________________
    if args.normalization_Yukawa_reimpl:

        normalizationCrossCheck = NormalizationCrossCheck()

        from LatestBinning import obs_pth_ggH as obs_pth
        obs_pth.YR4_totalXS = LatestBinning.YR4_ggF_n3lo
        obs_pth.keep_only_first_nbins(5)

        normalizationCrossCheck.set_bin_boundaries(obs_pth.binning, add_overflow=False)

        # Add SM
        SMXSs = obs_pth.crosssection_over_binwidth()
        normalizationCrossCheck.add_spectrum('YR4, shape Vitt.', SMXSs, isSM=True, color=4)

        # Add raw from Pier
        SM_Yukawa = TheoryFileInterface.file_finder(
            kappab=1, kappac=1, muR=1, muF=1, Q=1,
            expectOneFile=True, loadImmediately=True,
            directory=LatestPaths.derivedTheoryFiles_YukawaSummed
            )
        xs_Yukawa = TheoryCommands.rebin(SM_Yukawa.binBoundaries, SM_Yukawa.crosssection, obs_pth.binning, lastBinIsOverflow=False)
        normalizationCrossCheck.add_spectrum('Pier', xs_Yukawa, color=2)

        # Add SM from the ws (what is actually used)
        ws = 'out/workspaces_Jan29/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_couplingDependentBR_reweighted.root'
        with Commands.OpenRootFile(ws) as wsFp:
            w = wsFp.Get('w')
            reweightors = []
            for left, right in zip(obs_pth.binning[:-1], obs_pth.binning[1:]):
                name = 'reweightor_ggH_PTH_{0:d}_{1:d}'.format(int(left), int(right))
                try:
                    reweightors.append( w.function(name).getVal() )
                except TypeError:
                    print 'Error: left = {0}, right = {1}'.format(left, right)
                    print name
                    raise
            SMXS_ws = [ w.var(y).getVal() for y in Commands.list_set(w, 'SMXS') ]
        binWidths = [ right-left for left, right in zip(normalizationCrossCheck.bin_boundaries[:-1], normalizationCrossCheck.bin_boundaries[1:]) ]
        ws_xs = [ reweightor * SMXS / binWidth for reweightor, SMXS, binWidth in zip(reweightors, SMXS_ws, binWidths) ]
        normalizationCrossCheck.add_spectrum('Pier reweighted', ws_xs, line_style=2)

        normalizationCrossCheck.plot_spectra('Yukawa')


    #____________________________________________________________________
    if args.normalization_Top:

        LOGSCALE = True

        from LatestBinning import obs_pth_combWithHbbBinning as obs_pth
        obs_pth.YR4_totalXS = LatestBinning.YR4_ggF_n3lo
        SMXS = obs_pth.crosssection_over_binwidth()

        H_SMXS = obs_pth.basic_histogram()
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
        set_cmargins()

        base = get_plot_base(
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

        SM_Top = TheoryFileInterface.file_finder(
            ct=1, cg=0, cb=1, muR=1, muF=1, Q=1,
            expectOneFile=True, loadImmediately=True,
            # directory=LatestPaths.derivedTheoryFiles_TopHighPt
            directory=LatestPaths.derivedTheoryFiles_Top
            )
        SM_Top_normalization = BasicIntegral( SM_Top.binBoundaries, SM_Top.crosssection, lastBinIsOverflow=False )

        H_Top = BasicHistogramFromContainer( SM_Top, lastBinIsOverflow=False )
        # H_Top.Draw('HISTSAME')



        SM_Top_Rebinned = deepcopy(SM_Top)
        TheoryCommands.rebin_derived_theory_container( SM_Top_Rebinned, obs_pth.binning, lastBinIsOverflow=False )

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
        TheoryCommands.rebin_derived_theory_container( SM_Top_Rescaled_Rebinned, obs_pth.binning, lastBinIsOverflow=False )

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
        save_c( 'normalizationCrosscheck_Top' )


    #____________________________________________________________________
    if args.normalization_Yukawa:

        LOGSCALE = True

        from LatestBinning import obs_pth_combWithHbbBinning as obs_pth
        obs_pth.YR4_totalXS = LatestBinning.YR4_ggF_n3lo
        obs_pth.keep_only_first_nbins(5)

        H_SMXS = obs_pth.basic_histogram()
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
        set_cmargins()

        base = get_plot_base(
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

        SM_Yukawa = TheoryFileInterface.file_finder(
            kappab=1, kappac=1, muR=1, muF=1, Q=1,
            expectOneFile=True, loadImmediately=True,
            directory=LatestPaths.derivedTheoryFiles_YukawaSummed
            )
        SM_Yukawa_normalization = BasicIntegral( SM_Yukawa.binBoundaries, SM_Yukawa.crosssection, lastBinIsOverflow=False )

        H_Yukawa = BasicHistogramFromContainer( SM_Yukawa, lastBinIsOverflow=False )
        H_Yukawa.Draw('HISTSAME')



        SM_Yukawa_Rebinned = deepcopy(SM_Yukawa)
        TheoryCommands.rebin_derived_theory_container( SM_Yukawa_Rebinned, obs_pth.binning, lastBinIsOverflow=False )

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
        save_c( 'normalizationCrosscheck_Yukawa' )




#____________________________________________________________________
def BasicHistogramFromContainer( container, lastBinIsOverflow=False ):

    nBins = len(container.crosssection)

    H = ROOT.TH1F(
        TheoryCommands.get_unique_root_name(), TheoryCommands.get_unique_root_name(),
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
        Commands.throw_error( 'len(binning) = {0}, len(xss) = {1}; Cannot calculate integral.'.format( len(binning), len(xss_per_binWidth) ) )

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