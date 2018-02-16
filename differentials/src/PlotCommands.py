#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import Commands

import os, itertools, operator, re, argparse, sys, random, numpy
from math import isnan, isinf, log10
from os.path import *
from glob import glob
from copy import deepcopy
from array import array


from time import strftime
datestr = strftime( '%b%d' )

import ROOT


import TheoryCommands
from TheoryCommands import c
from TheoryCommands import SetCMargins
from TheoryCommands import SetPlotDir
from TheoryCommands import SaveAsRoot
from TheoryCommands import SaveC
from TheoryCommands import GetUniqueRootName
from TheoryCommands import GetPlotBase

from Container import Container

import PhysicsCommands
import CorrelationMatrices
from Parametrization import Parametrization, WSParametrization

from time import strftime
datestr = strftime( '%b%d' )


########################################
# Functions
########################################

#____________________________________________________________________
def set_color_palette(option=None):

    # n_stops = 3
    # stops  = [ 0.0, 0.5, 1.0 ]
    # reds   = [ 0.0, 1.0, 1.0 ]
    # blues  = [ 1.0, 1.0, 0.0 ]
    # greens = [ 0.0, 1.0, 0.0 ]

    # n_stops = 2
    # stops  = [ 0.0, 1.0 ]
    # reds   = [ 55./255.,  1.0 ]
    # greens = [ 138./255., 1.0 ]
    # blues  = [ 221./255., 1.0 ]

    if option == 'twocolor':
        n_stops = 3
        stops  = [ 0.0, 0.3, 1.0 ]
        reds   = [ 55./255.,  1.0, 1.0 ]
        greens = [ 138./255., 1.0, 26./255. ]
        blues  = [ 221./255., 1.0, 26./255. ]
    else:
        n_stops = 3
        stops  = [ 0.0, 0.3, 1.0 ]
        reds   = [ 55./255.,  166./255., 1.0 ]
        greens = [ 138./255., 203./255., 1.0 ]
        blues  = [ 221./255., 238./255., 1.0 ]

    ROOT.TColor.CreateGradientColorTable(
        n_stops,
        array('d', stops ),
        array('d', reds ),
        array('d', greens ),
        array('d', blues ),
        255 )


#____________________________________________________________________
class TLegendMultiPanel(object):
    """docstring for TLegendMultiPanel"""
    def __init__(
            self,
            x1, y1, x2, y2
            ):
        self._x1 = x1
        self._y1 = y1
        self._x2 = x2
        self._y2 = y2
        self._entries = []

        self._SetNColumns = None
        self._SetBorderSize = None
        self._SetFillStyle = None

    def AddEntry( self, *args ):
        self._entries.append( args )

    def set_n_columns( self, val ):
        self._SetNColumns = val
    def SetBorderSize( self, val ):
        self._SetBorderSize = val
    def SetFillStyle( self, val ):
        self._SetFillStyle = val

    def Draw( self, drawStr = '' ):

        if callable(self._x1): self._x1 = self._x1( ROOT.gPad )
        if callable(self._y1): self._y1 = self._y1( ROOT.gPad )
        if callable(self._x2): self._x2 = self._x2( ROOT.gPad )
        if callable(self._y2): self._y2 = self._y2( ROOT.gPad )

        leg = ROOT.TLegend( self._x1, self._y1, self._x2, self._y2 )
        ROOT.SetOwnership( leg, False )

        for args in self._entries:
            leg.AddEntry( *args )

        if not( self._SetNColumns is None ):
            leg.set_n_columns(   self._SetNColumns )
        if not( self._SetBorderSize is None ):
            leg.SetBorderSize( self._SetBorderSize )
        if not( self._SetFillStyle is None ):
            leg.SetFillStyle(  self._SetFillStyle )

        leg.Draw( drawStr )


#____________________________________________________________________
class TLatexMultiPanel(object):
    """docstring"""
    def __init__(
            self,
            x, y, text,
            ):
        self.x = x
        self.y = y
        self.text = text

        self._SetTextSize = None
        self._SetTextAlign = 31
        self._SetTextColor = None
        self._SetNDC = True
        self._SetTextFont = 42

    def SetTextSize( self, val ):
        self._SetTextSize = val
    def SetTextAlign( self, val ):
        self._SetTextAlign = val
    def SetTextColor( self, val ):
        self._SetTextColor = val
    def SetNDC( self, val ):
        self._SetNDC = val
    def SetTextFont( self, val ):
        self._SetTextFont = val

    def Draw( self, drawStr = '' ):

        if callable(self.x): self.x = self.x( ROOT.gPad )
        if callable(self.y): self.y = self.y( ROOT.gPad )

        l = ROOT.TLatex()
        ROOT.SetOwnership( l, False )

        if not self._SetTextSize is None:
            l.SetTextSize( self._SetTextSize )
        if not self._SetTextAlign is None:
            l.SetTextAlign( self._SetTextAlign )
        if not self._SetTextColor is None:
            l.SetTextColor( self._SetTextColor )
        if not self._SetNDC is None:
            l.SetNDC( self._SetNDC )
        if not self._SetTextFont is None:
            l.SetTextFont( self._SetTextFont )

        print 'In Draw():'
        print '  x    = ', self.x
        print '  y    = ', self.y
        print '  text = ', self.text

        l.DrawLatex( self.x, self.y, self.text )



#____________________________________________________________________
def plot_with_bottom_panel(
        plotName,
        topPanelObjects,
        bottomPanelObjects,
        centralPadObjects = None,
        xTitle        = '',
        yTitleTop     = '',
        yTitleBottom  = '',
        padSplitPoint = 0.33,
        helpLines     = False,
        # 
        topPadBottomMargin = 0.02,
        topPadTopMargin    = 0.10,
        topPadLeftMargin   = 0.12,
        topPadRightMargin  = 0.02,
        # 
        bottomPadTopMargin    = 0.00,
        bottomPadBottomMargin = 0.30,
        bottomPadLeftMargin   = 0.12,
        bottomPadRightMargin  = 0.02,
        # 
        SetTopPanelLogScale = False,
        # 
        disableCMSText = False
        ):

    c.Clear()
    # c.SetRightMargin( 0.03 )
    # c.SetLeftMargin(  0.12 )

    _width = c.GetWindowWidth()
    _height = c.GetWindowHeight()
    c.SetCanvasSize( 800, 900 )


    topPad_bottom    = padSplitPoint
    bottomPad_top    = padSplitPoint

    heightRatio = float( 1.0 - topPad_bottom ) / float( bottomPad_top )


    topPad = ROOT.TPad(
        get_unique_root_name(), '',
        # c.GetLeftMargin(), topPad_bottom, 1-c.GetRightMargin(), topPad_top
        0.0, topPad_bottom, 1.0, 1.0
        )
    topPad.SetBottomMargin( topPadBottomMargin) # Distance to the bottom panel
    topPad.SetTopMargin(    topPadTopMargin)     # Space for labels
    topPad.SetLeftMargin(   topPadLeftMargin)
    topPad.SetRightMargin(  topPadRightMargin)
    topPad.Draw()

    bottomPad = ROOT.TPad(
        get_unique_root_name(), '',
        # c.GetLeftMargin(), bottomPad_bottom, 1-c.GetRightMargin(), bottomPad_top
        0.0, 0.0, 1.0, bottomPad_top
        )
    bottomPad.SetTopMargin(    bottomPadTopMargin)    # Distance to the bottom panel
    bottomPad.SetBottomMargin( bottomPadBottomMargin) # Space for labels
    bottomPad.SetLeftMargin(   bottomPadLeftMargin)
    bottomPad.SetRightMargin(  bottomPadRightMargin)
    bottomPad.Draw()


    # ======================================
    # Help lines

    if helpLines:
        for pad in [ topPad, bottomPad ]:
            pad.cd()
            for x1, y1, x2, y2 in [
                    ( 0.0, 0.0, 1.0, 0.0 ),
                    ( 1.0, 0.0, 1.0, 1.0 ),
                    ( 0.0, 0.0, 0.0, 1.0 ),
                    ( 0.0, 1.0, 1.0, 1.0 ),
                    ]:
                line = ROOT.TLine( x1, y1, x2, y2 )
                ROOT.SetOwnership( line, False )
                line.Draw()

        c.cd()
        line = ROOT.TLine( 0.0, bottomPad_top, 1.0, topPad_bottom )
        ROOT.SetOwnership( line, False )
        line.Draw()


    # ======================================
    # Draw objects

    topPad.cd()
    if SetTopPanelLogScale: topPad.SetLogy()

    # print '\nPrinting ROOT.gPad stats'

    # print ROOT.gPad.GetX1()
    # print ROOT.gPad.GetY1()
    # print ROOT.gPad.GetX2()
    # print ROOT.gPad.GetY2()

    # print ROOT.gPad.GetXlowNDC()
    # print ROOT.gPad.GetYlowNDC()
    # print ROOT.gPad.GetXHNDC()
    # print ROOT.gPad.GetYHNDC()


    for obj, drawStr in topPanelObjects:
        obj.Draw(drawStr)

    bottomPad.cd()
    for obj, drawStr in bottomPanelObjects:
        obj.Draw(drawStr)



    topPad.cd()
    axisHolderTop = topPanelObjects[0][0]
    axisHolderTop.GetXaxis().SetLabelOffset(999.)
    axisHolderTop.GetYaxis().SetTitle(yTitleTop)


    bottomPad.cd()
    axisHolderBottom = bottomPanelObjects[0][0]

    axisHolderBottom.GetXaxis().SetLabelSize( axisHolderTop.GetXaxis().GetLabelSize() * heightRatio )
    axisHolderBottom.GetYaxis().SetLabelSize( axisHolderTop.GetYaxis().GetLabelSize() * heightRatio )

    axisHolderBottom.GetXaxis().SetTickLength( axisHolderTop.GetXaxis().GetTickLength() * heightRatio )
    # axisHolderBottom.GetYaxis().SetTickLength( axisHolderTop.GetYaxis().GetTickLength() * heightRatio )

    axisHolderBottom.GetYaxis().SetTitle(yTitleBottom)
    axisHolderBottom.GetXaxis().SetTitle(xTitle)
    axisHolderBottom.GetXaxis().SetTitleSize( axisHolderTop.GetXaxis().GetTitleSize() * heightRatio )
    axisHolderBottom.GetYaxis().SetTitleSize( axisHolderTop.GetYaxis().GetTitleSize() * heightRatio )
    axisHolderBottom.GetYaxis().SetTitleOffset( 1./heightRatio)

    if not disableCMSText:
        topPad.cd()
        Commands.get_cms_label( textSize=0.08 )
        Commands.get_cms_lumi( textSize=0.07 )

    save_c(plotName)
    c.SetCanvasSize( _width, _height )

#____________________________________________________________________
def get_uncs_from_scans(container):
    container.binBoundaries = PhysicsCommands.figure_out_binning( container.POIs )
    container.uncs = []
    for POI, scan in zip( container.POIs, container.Scans ):
        POIvals, deltaNLLs = PhysicsCommands.filter_scan( scan )
        unc = PhysicsCommands.find_minima_and_errors( POIvals, deltaNLLs, returnContainer=True )
        container.uncs.append(unc)

    container.nBins         = len(container.binBoundaries)-1
    container.binCenters    = [ 0.5*(left+right) for left, right in zip( container.binBoundaries[:-1], container.binBoundaries[1:] ) ]
    container.binWidths     = [ right-left for left, right in zip( container.binBoundaries[:-1], container.binBoundaries[1:] ) ]
    container.halfBinWidths = [ 0.5*(right-left) for left, right in zip( container.binBoundaries[:-1], container.binBoundaries[1:] ) ]





#____________________________________________________________________
class GenericHistogram(object):
    """docstring for GenericHistogram"""

    colorCycle = Commands.new_color_cycle()
    fillStyleCycle = itertools.cycle([ 3245, 3254, 3205 ])

    def __init__(self, name, title, bin_boundaries, bin_values, color=None):
        super(GenericHistogram, self).__init__()
        self.name = name
        self.title = title
        self.bin_boundaries = bin_boundaries
        self.n_bins = len(bin_boundaries)-1
        self.bin_values = bin_values
        self.bins = []
        self._expecting_asymm_unc = False
        if color is None:
            self.color = self.colorCycle.next()
        else:
            self.color = color
        
    def set_err_up(self, errs_up):
        self.errs_up = [ abs(i) for i in errs_up ]
        self.bounds_up = [ c+e for c, e in zip(self.bin_values, self.errs_up) ]
        self._expecting_asymm_err = True

    def set_err_down(self, errs_down):
        self.errs_down = [ abs(i) for i in errs_down ]
        self.bounds_down = [ c+e for c, e in zip(self.bin_values, self.errs_down) ]
        self._expecting_asymm_err = True

    def get_bin_centers(self):
        return [ 0.5*(left+right) for left, right in zip(self.bin_boundaries[:-1], self.bin_boundaries[1:]) ]

    def get_bin_widths(self):
        return [ right-left for left, right in zip(self.bin_boundaries[:-1], self.bin_boundaries[1:]) ]

    def get_half_bin_widths(self):
        return [ 0.5*i for i in self.get_bin_widths() ]

    def get_zeroes(self):
        return [0.0 for i in xrange(self.n_bins)]


    def repr_basic_histogram(self, leg=None):
        H = ROOT.TH1F(
            get_unique_root_name(), '',
            len(self.bin_boundaries)-1, array( 'f', self.bin_boundaries)
            )
        ROOT.SetOwnership( H, False )
        H.SetLineColor(self.color)
        H.SetLineWidth(2)
        for i_bin in xrange(self.n_bins):
            H.SetBinContent( i_bin+1, self.bin_values[i_bin] )

        if not(leg is None):
            leg.AddEntry(H.GetName(), self.title, 'l')

        return [ (H, 'HISTSAME') ]


    def repr_horizontal_bars(self):
        Tg = ROOT.TGraphErrors(
            self.n_bins,
            array( 'f', self.get_bin_centers() ),
            array( 'f', self.bin_values ),
            array( 'f', self.get_half_bin_widths() ),
            array( 'f', self.get_zeroes() ),
            )
        ROOT.SetOwnership( Tg, False )
        Tg.SetName( get_unique_root_name() )

        Tg.SetLineColor(   getattr(self, 'setLineColor',   self.color ) )
        # Tg.SetLineWidth(   getattr(self, 'setLineWidth',   2 ) )

        return [ (Tg, 'EPSAME') ]


    def repr_uncertainties_filled_area(self, leg=None):
        Tg = ROOT.TGraphAsymmErrors(
            self.n_bins,
            array( 'f', self.get_bin_centers() ),
            array( 'f', self.bin_values ),
            array( 'f', [ 0.45*w for w in self.get_bin_widths() ] ),
            array( 'f', [ 0.45*w for w in self.get_bin_widths() ] ),
            array( 'f', self.errs_down ),
            array( 'f', self.errs_up ),
            )
        ROOT.SetOwnership( Tg, False )
        Tg.SetName( get_unique_root_name() )
    
        Tg.SetFillStyle(   getattr(self, 'setFillStyle',   3245 ) )
        Tg.SetMarkerStyle( getattr(self, 'setMarkerStyle', 8 ) )
        Tg.SetMarkerSize(  getattr(self, 'setMarkerSize',  0 ) )
        Tg.SetFillColor(   getattr(self, 'setFillColor',   self.color ) )
        Tg.SetMarkerColor( getattr(self, 'setMarkerColor', self.color ) )
        Tg.SetLineColor(   getattr(self, 'setLineColor',   self.color ) )

        if not(leg is None):
            leg.AddEntry( Tg.GetName(), self.title, 'LF' )

        return [ (Tg, 'E2PSAME') ]


    def repr_uncertainties_narrow_filled_area(self, leg=None):
        Tg = ROOT.TGraphAsymmErrors(
            self.n_bins,
            array( 'f', self.get_bin_centers() ),
            array( 'f', self.bin_values ),
            array( 'f', [ 0.1*w for w in self.get_bin_widths() ] ),
            array( 'f', [ 0.1*w for w in self.get_bin_widths() ] ),
            array( 'f', self.errs_down ),
            array( 'f', self.errs_up ),
            )
        ROOT.SetOwnership( Tg, False )
        Tg.SetName( get_unique_root_name() )
    
        Tg.SetMarkerStyle( getattr(self, 'setMarkerStyle', 8 ) )
        Tg.SetMarkerSize(  getattr(self, 'setMarkerSize',  0 ) )
        Tg.SetMarkerColor( getattr(self, 'setMarkerColor', self.color ) )
        Tg.SetLineColor(   getattr(self, 'setLineColor',   self.color ) )

        # Tg.SetFillStyle(   getattr(self, 'setFillStyle',   3245 ) )
        # Tg.SetFillColor(   getattr(self, 'setFillColor',   self.color ) )
        Tg.SetFillColorAlpha(self.color, 0.30)

        if not(leg is None):
            leg.AddEntry( Tg.GetName(), self.title, 'LF' )

        return [ (Tg, 'E2PSAME') ]


    def repr_uncertainties_fully_filled_area(self, leg=None):
        Tg = ROOT.TGraphAsymmErrors(
            self.n_bins,
            array( 'f', self.get_bin_centers() ),
            array( 'f', self.bin_values ),
            array( 'f', self.get_half_bin_widths() ),
            array( 'f', self.get_half_bin_widths() ),
            array( 'f', self.errs_down ),
            array( 'f', self.errs_up ),
            )
        ROOT.SetOwnership( Tg, False )
        Tg.SetName( get_unique_root_name() )
    
        Tg.SetMarkerStyle( getattr(self, 'setMarkerStyle', 8 ) )
        Tg.SetMarkerSize(  getattr(self, 'setMarkerSize',  0 ) )
        Tg.SetMarkerColor( getattr(self, 'setMarkerColor', self.color ) )
        Tg.SetLineColor(   getattr(self, 'setLineColor',   self.color ) )

        # Tg.SetFillStyle(   getattr(self, 'setFillStyle',   3245 ) )
        # Tg.SetFillColor(   getattr(self, 'setFillColor',   self.color ) )
        Tg.SetFillColorAlpha(self.color, 0.30)

        if not(leg is None):
            leg.AddEntry( Tg.GetName(), self.title, 'LF' )

        return [ (Tg, 'E2PSAME') ]


    def repr_point_with_vertical_bar(self, leg=None):

        Tg = ROOT.TGraphAsymmErrors(
            self.n_bins,
            array( 'f', self.get_bin_centers() ),
            array( 'f', self.bin_values ),
            array( 'f', [ 0.0 for i in xrange(self.n_bins) ] ),
            array( 'f', [ 0.0 for i in xrange(self.n_bins) ] ),
            array( 'f', self.errs_down ),
            array( 'f', self.errs_up ),
            )
        ROOT.SetOwnership( Tg, False )
        Tg.SetName( get_unique_root_name() )

        Tg.SetMarkerStyle( getattr(self, 'setMarkerStyle', 8 ) )
        Tg.SetFillColor(   getattr(self, 'setFillColor',   self.color ) )
        Tg.SetMarkerColor( getattr(self, 'setMarkerColor', self.color ) )
        Tg.SetLineColor(   getattr(self, 'setLineColor',   self.color ) )

        if not(leg is None):
            leg.AddEntry( Tg.GetName(), self.title, 'PE' )

        return [(Tg, 'PSAME')]


    def Draw(self, draw_style):
        for obj, draw_str in getattr(self, draw_style):
            obj.Draw(draw_str)


#____________________________________________________________________
def plot_spectra_on_two_panel(
        plotname,
        containers,
        xTitle       = 'p_{T}^{H}',
        yTitleTop    = '#Delta#sigma (pb/GeV)',
        yMinLimit = 0.001,
        yMinExternalTop = None,
        yMaxExternalTop = None,
        lastBinIsNotOverflow = False,
        xMinExternal = None,
        xMaxExternal = None,
        verbose = True,
        # 
        scaleLastBin = 'previousBin',
        # 
        topPanelObjects = None,
        bottomPanelObjects = None,
        ):

    if verbose: print 'Plotting \'{0}\''.format(plotname)


    TOP_PANEL_LOGSCALE = True

    addLaterToTopPanelObjects = []
    addLaterToBottomPanelObjects = []
    if not topPanelObjects is None:
        addLaterToTopPanelObjects = topPanelObjects
    if not bottomPanelObjects is None:
        addLaterToBottomPanelObjects = bottomPanelObjects
    topPanelObjects    = []
    bottomPanelObjects = []

    # ---------------------
    # Determine extrema of plot

    # Filter out passed containers that have hardcoded cross sections
    hardcoded_numbers_containers = []
    _containers = []
    for container in containers:
        if hasattr(container, 'crosssection'):
            hardcoded_numbers_containers.append(container)
        else:
            _containers.append(container)
    containers = _containers

    for container in containers:
        get_uncs_from_scans(container)

        # If last bin is indeed an overflow, and scaleLastBin is not set to False, scale the last bin
        if not lastBinIsNotOverflow and not scaleLastBin == False:
            if scaleLastBin == 'previousBin':
                container.scale = container.binBoundaries[-1] - container.binBoundaries[-2]
            else:
                container.scale = scaleLastBin
            container.SMcrosssections[-1] /= container.scale


        container.xMin = container.binBoundaries[0]
        if lastBinIsNotOverflow:
            container.xMax = container.binBoundaries[-1]
        else:
            container.xMax = container.binBoundaries[-2] + ( container.binBoundaries[-2] - container.binBoundaries[-3] )
        container.yMinRatio = min([ unc.leftBound for unc in container.uncs ])
        container.yMaxRatio = max([ unc.rightBound for unc in container.uncs ])
        if TOP_PANEL_LOGSCALE:
            container.yMinCrosssection = min([
                min([ unc.leftBound  * xs for unc, xs in zip( container.uncs, container.SMcrosssections ) if unc.leftBound > 0.0 ] + [0.01] ),
                min([ unc.min        * xs for unc, xs in zip( container.uncs, container.SMcrosssections ) if unc.min > 0.0 ]),
                min([ unc.rightBound * xs for unc, xs in zip( container.uncs, container.SMcrosssections ) if unc.rightBound > 0.0 ]),
                ]) 
        else:
            container.yMinCrosssection = min([ unc.leftBound * xs for unc, xs in zip( container.uncs, container.SMcrosssections ) ])
        container.yMaxCrosssection = max([ unc.rightBound * xs for unc, xs in zip( container.uncs, container.SMcrosssections ) ])

    xMin = min([ container.xMin for container in containers ])
    xMax = max([ container.xMax for container in containers ])
    if not xMinExternal is None:
        xMin = xMinExternal
    if not xMaxExternal is None:
        xMax = xMaxExternal

    for container in containers:
        if container.binBoundaries[-1] < xMax: container.binBoundaries[-1] = xMax
        container.nBins         = len(container.binBoundaries)-1
        container.binCenters    = [ 0.5*(left+right) for left, right in zip( container.binBoundaries[:-1], container.binBoundaries[1:] ) ]
        container.binWidths     = [ right-left for left, right in zip( container.binBoundaries[:-1], container.binBoundaries[1:] ) ]
        container.halfBinWidths = [ 0.5*(right-left) for left, right in zip( container.binBoundaries[:-1], container.binBoundaries[1:] ) ]


    # Bottom plot base
    yMinAbsBottom = min([ container.yMinRatio for container in containers if not abs(container.yMinRatio) == 999. ])
    yMaxAbsBottom = max([ container.yMaxRatio for container in containers if not abs(container.yMaxRatio) == 999. ])
    yMinBottom = yMinAbsBottom - 0.1*(yMaxAbsBottom-yMinAbsBottom)
    yMaxBottom = yMaxAbsBottom + 0.1*(yMaxAbsBottom-yMinAbsBottom)

    baseBottom = get_plot_base(
        xMin = xMin,
        xMax = xMax,
        yMin = yMinBottom,
        yMax = yMaxBottom,
        )
    bottomPanelObjects.append( ( baseBottom, 'P' ) )

    lineAtOne = ROOT.TLine( xMin, 1.0, xMax, 1.0 )
    lineAtOne.SetLineColor(14)
    bottomPanelObjects.append( ( lineAtOne, '' ) )


    # Top plot base
    yMinAbsTop = min([ container.yMinCrosssection for container in containers ])
    yMaxAbsTop = max([ container.yMaxCrosssection for container in containers ])

    if TOP_PANEL_LOGSCALE:
        yMaxTop = yMaxAbsTop + ( yMaxAbsTop / yMinAbsTop )**0.2 * yMaxAbsTop
        # yMinTop = yMinAbsTop - yMinAbsTop * ( yMaxAbsTop / yMinAbsTop )**-0.2
        yMinTop = 0.5*yMinAbsTop
    else:
        yMinTop = yMinAbsTop - 0.1*(yMaxAbsTop-yMinAbsTop)
        yMaxTop = yMaxAbsTop + 0.1*(yMaxAbsTop-yMinAbsTop)

    if yMaxExternalTop: 
        yMaxTop = yMaxExternalTop
    if yMinExternalTop: 
        yMinTop = yMinExternalTop

    baseTop = get_plot_base(
        xMin = xMin,
        xMax = xMax,
        yMin = yMinTop,
        yMax = yMaxTop,
        )
    topPanelObjects.append( ( baseTop, 'P' ) )


    # ======================================
    # Load results into plot

    leg = TLegendMultiPanel(
        lambda c: c.GetLeftMargin() + 0.01,
        lambda c: 1 - c.GetTopMargin() - 0.10,
        lambda c: 1 - c.GetRightMargin() - 0.01,
        lambda c: 1 - c.GetTopMargin()
        )
    leg.set_n_columns( min( len(containers)+len(hardcoded_numbers_containers), 4 ) )
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    colorCycle = Commands.new_color_cycle()
    fillStyleCycle = itertools.cycle([ 3245, 3254, 3205 ])


    for container in hardcoded_numbers_containers:
        if not hasattr( container, 'color' ): container.color = next(colorCycle)
        container.nBins = getattr(container, 'nBins', len(container.binBoundaries)-1)

        Hratio = ROOT.TH1F(
            get_unique_root_name(), '',
            len(container.binBoundaries)-1, array( 'f', container.binBoundaries )
            )
        ROOT.SetOwnership( Hratio, False )
        Hratio.SetLineColor(container.color)
        Hratio.SetLineWidth(2)
        for iBin in xrange(container.nBins):
            Hratio.SetBinContent( iBin+1, container.ratios[iBin] )
        bottomPanelObjects.append( ( Hratio, 'HISTSAME' ) )

        Hcrosssection = ROOT.TH1F(
            get_unique_root_name(), '',
            len(container.binBoundaries)-1, array( 'f', container.binBoundaries )
            )
        ROOT.SetOwnership( Hcrosssection, False )
        Hcrosssection.SetLineColor(container.color)
        Hcrosssection.SetLineWidth(2)
        for iBin in xrange(container.nBins):
            Hcrosssection.SetBinContent( iBin+1, container.crosssection[iBin] )
        topPanelObjects.append( ( Hcrosssection, 'HISTSAME' ) )

        leg.AddEntry(Hcrosssection.GetName(), container.title, 'l')


    for container in containers:
        if not hasattr( container, 'color' ): container.color = next(colorCycle)
        container.fillStyle = fillStyleCycle.next()

        # ---------------------
        # Construct objects for bottom panel

        ratios = [container.uncs[i].min for i in xrange(container.nBins)]
        H_ratio_generic = GenericHistogram(container.name, container.title, container.binBoundaries, ratios, color=container.color)
        H_ratio_generic.set_err_up([ container.uncs[i].rightError for i in xrange(container.nBins) ])
        H_ratio_generic.set_err_down([ container.uncs[i].leftError for i in xrange(container.nBins) ])

        crosssections = [container.uncs[i].min * container.SMcrosssections[i] for i in xrange(container.nBins)]
        H_crosssection_generic = GenericHistogram(container.name, container.title, container.binBoundaries, crosssections, color=container.color)
        H_crosssection_generic.set_err_up([ container.uncs[i].rightError * container.SMcrosssections[i] for i in xrange(container.nBins) ])
        H_crosssection_generic.set_err_down([ container.uncs[i].leftError * container.SMcrosssections[i] for i in xrange(container.nBins) ])

        if not(container.name == 'combination' or container.name == 'combWithHbb') and not(getattr(container, 'error_line', False)):
            bottomPanelObjects.extend(H_ratio_generic.repr_horizontal_bars())
            bottomPanelObjects.extend(H_ratio_generic.repr_uncertainties_narrow_filled_area())
            # bottomPanelObjects.extend(H_ratio_generic.repr_uncertainties_fully_filled_area())
            topPanelObjects.extend(H_crosssection_generic.repr_horizontal_bars())
            topPanelObjects.extend(H_crosssection_generic.repr_uncertainties_narrow_filled_area(leg if not getattr(container, 'suppress_text', False) else None))
            # topPanelObjects.extend(H_crosssection_generic.repr_uncertainties_fully_filled_area(leg if not getattr(container, 'suppress_text', False) else None))
        else:
            bottomPanelObjects.extend(H_ratio_generic.repr_point_with_vertical_bar())
            topPanelObjects.extend(H_crosssection_generic.repr_point_with_vertical_bar(leg if not getattr(container, 'suppress_text', False) else None))


        # ======================================
        # Print table of uncertainties and cross sections

        strList = lambda L: ', '.join([ '{0:<9.3f}'.format(f) for f in L ])
        print '\nSignal strengths and cross sections in {0}, {1}:'.format( plotname, container.name )
        print 'binBoundaries : ' + strList( container.binBoundaries )
        print 'mu            : ' + strList( [ container.uncs[i].min for i in xrange(container.nBins) ] )
        print 'mu_error_down : ' + strList( [ container.uncs[i].leftError for i in xrange(container.nBins) ] )
        print 'mu_error_up   : ' + strList( [ container.uncs[i].rightError for i in xrange(container.nBins) ] )
        print 'mu_bound_down : ' + strList( [ container.uncs[i].leftBound for i in xrange(container.nBins) ] )
        print 'mu_bound_up   : ' + strList( [ container.uncs[i].rightBound for i in xrange(container.nBins) ] )
        print 'xs            : ' + strList( [ container.SMcrosssections[i] * container.uncs[i].min for i in xrange(container.nBins) ] )
        print 'xs_error_down : ' + strList( [ container.SMcrosssections[i] * container.uncs[i].leftError for i in xrange(container.nBins) ] )
        print 'xs_error_up   : ' + strList( [ container.SMcrosssections[i] * container.uncs[i].rightError for i in xrange(container.nBins) ] )
        print 'xs_bound_down : ' + strList( [ container.SMcrosssections[i] * container.uncs[i].leftBound for i in xrange(container.nBins) ] )
        print 'xs_bound_up   : ' + strList( [ container.SMcrosssections[i] * container.uncs[i].rightBound for i in xrange(container.nBins) ] )


    # ======================================
    # Add label for last bin regarding scaling

    if not lastBinIsNotOverflow:
        textSize = 0.03

        # Pick highest y to start (including combination, this only matters for y height)
        yOverflow = max([ container.SMcrosssections[-1] * container.uncs[-1].rightBound for container in containers ])

        # Add a label per scaled spectrum
        for i, container in enumerate([ con for con in containers if not con.name == 'combination']):
            if getattr(container, 'suppress_text', False): continue

            yOverflowInNDC = lambda c, i_lambda=i: (
                c.GetBottomMargin()
                + log10(yOverflow/yMinTop)/log10(yMaxTop/yMinTop) * ( 1.0 - c.GetTopMargin() - c.GetBottomMargin() )
                + i_lambda * 4.5*textSize
                + 2.5*textSize
                )

            scaleText = str(int(container.scale)) if float(container.scale).is_integer() else '{0:.2f}'.format(container.scale)
            l = TLatexMultiPanel(
                    lambda c: 1. - c.GetRightMargin() - 0.01,
                    yOverflowInNDC,
                    '#int_{{#lower[0.3]{{{0}}}}}^{{#lower[0.0]{{#infty}}}}#sigma({1}) d{1} / {2}'.format(
                        container.binBoundaries[-2],
                        xTitle.replace(' (GeV)',''),
                        scaleText
                        )
                    )
            l.SetTextSize(textSize)
            l.SetTextColor(container.color)
            topPanelObjects.append( ( l, '' ) )



    topPanelObjects.append( ( leg, '' ) )

    # Add objects that were given in the input
    for obj in addLaterToTopPanelObjects:
        topPanelObjects.append(obj)
    for obj in addLaterToBottomPanelObjects:
        bottomPanelObjects.append(obj)

    ROOT.gStyle.SetEndErrorSize(3)
    plot_with_bottom_panel(
        plotname,
        topPanelObjects,
        bottomPanelObjects,
        xTitle = xTitle,
        yTitleTop    = yTitleTop,
        yTitleBottom = 'ratio w.r.t. SM',
        SetTopPanelLogScale = TOP_PANEL_LOGSCALE,
        topPadLeftMargin = 0.14,
        bottomPadLeftMargin = 0.14,
        )
    ROOT.gStyle.SetEndErrorSize(1)


#____________________________________________________________________
def write_scans_to_table(
        container,
        tag,
        xTitle       = 'p_{T}^{H}',
        yTitle       = '#Delta#sigma (pb/GeV)',
        lastBinIsOverflow = True,
        verbose = True,
        scaleLastBin = False,
        ):

    container.binBoundaries = PhysicsCommands.figure_out_binning( container.POIs )
    container.nBins = len(container.binBoundaries)-1
    
    # Determine uncertainties from scan
    container.uncs = []
    for POI, scan in zip( container.POIs, container.Scans ):
        POIvals, deltaNLLs = PhysicsCommands.filter_scan( scan )
        unc = PhysicsCommands.find_minima_and_errors( POIvals, deltaNLLs, returnContainer=True )
        container.uncs.append(unc)

    # If last bin is indeed an overflow, and scaleLastBin is not set to False, scale the last bin
    if lastBinIsOverflow and not scaleLastBin == False:
        if scaleLastBin == 'previousBin':
            container.scale = container.binBoundaries[-1] - container.binBoundaries[-2]
        else:
            container.scale = scaleLastBin
        container.SMcrosssections[-1] /= container.scale

    # ======================================
    # Print table of uncertainties and cross sections

    strList = lambda L: ', '.join([ '{0:<9.3f}'.format(f) for f in L ])
    print '\nSignal strengths and cross sections in {0}:'.format( container.name )
    print 'binBoundaries : ' + strList( container.binBoundaries )
    print 'mu            : ' + strList( [ container.uncs[i].min for i in xrange(container.nBins) ] )
    print 'mu_error_down : ' + strList( [ container.uncs[i].leftError for i in xrange(container.nBins) ] )
    print 'mu_error_up   : ' + strList( [ container.uncs[i].rightError for i in xrange(container.nBins) ] )
    print 'mu_bound_down : ' + strList( [ container.uncs[i].leftBound for i in xrange(container.nBins) ] )
    print 'mu_bound_up   : ' + strList( [ container.uncs[i].rightBound for i in xrange(container.nBins) ] )
    print 'xs            : ' + strList( [ container.SMcrosssections[i] * container.uncs[i].min for i in xrange(container.nBins) ] )
    print 'xs_error_down : ' + strList( [ container.SMcrosssections[i] * container.uncs[i].leftError for i in xrange(container.nBins) ] )
    print 'xs_error_up   : ' + strList( [ container.SMcrosssections[i] * container.uncs[i].rightError for i in xrange(container.nBins) ] )
    print 'xs_bound_down : ' + strList( [ container.SMcrosssections[i] * container.uncs[i].leftBound for i in xrange(container.nBins) ] )
    print 'xs_bound_up   : ' + strList( [ container.SMcrosssections[i] * container.uncs[i].rightBound for i in xrange(container.nBins) ] )

    # ======================================
    # Construct table

    def f(number, decimals=2, force_sign=False):
        sign = '+' if force_sign else ''
        if isinstance(number, int) or number.is_integer():
            ret =  str(int(number))
        else:
            basic_float_str = '{0:{sign}.{decimals}f}'.format(number, decimals=decimals, sign=sign)
            float_str = "{0:{sign}.{decimals}E}".format(number, decimals=decimals, sign=sign)
            if "E" in float_str:
                base, exponent = float_str.split("E")
                exponent = int(exponent)
                if exponent <= -3:
                    ret = "{{{0} \\times 10^{{{1}}}}}".format(base, exponent)
                else:
                    ret = basic_float_str
            else:
                ret = basic_float_str
        return ret

    table = []

    binBoundariesLine = ['$' + xTitle.replace('#','\\') + '$']
    crosssectionLine  = ['$' + yTitle.replace('#','\\') + '$']
    ratiosLine        = ['$' + '\\mu' + '$']

    for iBin in xrange(container.nBins):

        if iBin < container.nBins-1 or not(lastBinIsOverflow):
            left  = f(container.binBoundaries[iBin])
            right = f(container.binBoundaries[iBin+1])
            binBoundariesLine.append('${0}-{1}$'.format(left, right))
        else:
            left  = f(container.binBoundaries[iBin])
            binBoundariesLine.append('${0}-\\infty$'.format(left))

        xs      = container.SMcrosssections[iBin] * container.uncs[iBin].min
        xs_down = -abs(container.SMcrosssections[iBin] * container.uncs[iBin].leftError)
        xs_up   = container.SMcrosssections[iBin] * container.uncs[iBin].rightError

        # xs_down_perc = xs_down / xs * 100. if not xs==0.0 else 0.0
        # xs_up_perc   = xs_up   / xs * 100. if not xs==0.0 else 0.0
        # xsStr = '{0:.1f}_{{{1:+.1f} ({2:.1f}\\%)}}^{{{3:+.1f} ({4:.1f}\\%)}}'.format(
        #     xs, xs_down, xs_down_perc, xs_up, xs_up_perc
        #     )
        xsStr = '${0}_{{{1}}}^{{{2}}}$'.format(
            f(xs), f(xs_down), f(xs_up, force_sign=True)
            )

        crosssectionLine.append(xsStr)

        mu      = container.uncs[iBin].min
        mu_down = -abs(container.uncs[iBin].leftError)
        mu_up   = container.uncs[iBin].rightError
        muStr = '${0}_{{{1}}}^{{{2}}}$'.format(
            f(mu), f(mu_down), f(mu_up, force_sign=True)
            )
        ratiosLine.append(muStr)

    table.append(binBoundariesLine)
    table.append(crosssectionLine)
    table.append(ratiosLine)

    tableText = '% ' + Commands.tag_git_commit_and_module() + '\n'
    tableText += Commands.print_table(table, maxColWidth=100, sep=' & ', newline_sep=' \\\\\n' )

    outname = join(TheoryCommands.PLOTDIR, 'mutable_{0}_{1}.tex'.format(tag, container.name))
    with open(outname, 'w') as outFp:
        outFp.write(tableText)


#____________________________________________________________________
def plot_parametrizations_on_combination(
        container,
        OnOneCanvas = False
        ):

    def debug(txt):
        print '[debug]', txt
        sys.exit()

    expBinBoundaries    = container.expBinBoundaries
    ws_combination      = container.ws_combination
    scanDir_combination = container.scanDir_combination
    ws_coupling         = container.ws_coupling
    scanDir_coupling    = container.scanDir_coupling
    xCoupling           = container.xCoupling
    yCoupling           = container.yCoupling
    xCouplingTitle      = container.xCouplingTitle
    yCouplingTitle      = container.yCouplingTitle
    plotTitle           = container.plotTitle


    DO_X_MAXIMA = False
    if hasattr( container, 'OnlyXMaximaOfContour' ) and container.OnlyXMaximaOfContour:
        DO_X_MAXIMA = True

    DO_Y_MAXIMA = False
    if hasattr( container, 'OnlyYMaximaOfContour' ) and container.OnlyYMaximaOfContour:
        DO_Y_MAXIMA = True

    if hasattr( container, 'MaximaOfContour' ) and container.MaximaOfContour:
        DO_X_MAXIMA = True
        DO_Y_MAXIMA = True


    DO_BESTFIT = False
    if hasattr( container, 'BestFit' ) and container.BestFit:
        DO_BESTFIT = True

    STRAIGHT_LINE_TO_SM = False
    if hasattr( container, 'StraightLineToSM' ) and container.StraightLineToSM:
        STRAIGHT_LINE_TO_SM = True
        if not hasattr( container, 'xSM' ):
            Commands.throw_error( 'Requested mode \'StraightLineToSM\', but xSM is not specified' )
        if not hasattr( container, 'ySM' ):
            Commands.throw_error( 'Requested mode \'StraightLineToSM\', but ySM is not specified' )


    # ======================================
    # Begin plotting

    newColorCycle = lambda: itertools.cycle( [ 2, 4, 6, 41, 46, 30, 43, 3, 5, 8, 9 ] )

    c.Clear()
    set_cmargins( TopMargin = 0.09 )

    xMin = expBinBoundaries[0]
    xMax = expBinBoundaries[-2] + ( expBinBoundaries[-2] - expBinBoundaries[-3] ) # Overflow will screw up plot
    yMin = 0.0
    yMax = 1.0

    base = get_plot_base(
        xMin = xMin,
        xMax = xMax,
        yMin = yMin,
        yMax = yMax,
        xTitle = 'p_{T} [GeV]', yTitle = '#mu'
        )
    base.Draw('P')

    if OnOneCanvas:
        leg = ROOT.TLegend(
            c.GetLeftMargin(),
            1 - c.GetTopMargin() - 0.25,
            c.GetLeftMargin() + 0.5*( 1. - c.GetRightMargin() - c.GetLeftMargin() ),
            1 - c.GetTopMargin() 
            )
        leg.set_n_columns(1)
    else:
        leg = ROOT.TLegend(
            c.GetLeftMargin(),
            1 - c.GetTopMargin() - 0.17,
            1 - c.GetRightMargin(),
            1 - c.GetTopMargin() 
            )
        leg.set_n_columns(2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)


    # ======================================
    # Draw the combination as blocks

    combinationPOIs = Commands.list_pois( ws_combination )
    combinationscans = PhysicsCommands.get_scan_results( combinationPOIs, scanDir_combination, pattern = '' )

    TgCombination = PhysicsCommands.get_TGraph_for_spectrum(
        combinationPOIs,
        combinationscans,
        name = get_unique_root_name()
        )

    TgCombination.SetLineColor( 1 )
    # TgCombination.SetMarkerStyle( 2 )
    TgCombination.SetFillColorAlpha( 1, 0.2 )
    # TgCombination.SetFillColor( 13 )
    # TgCombination.SetFillStyle( 3544 )
    # TgCombination.SetFillStyle( 3345 )

    TgCombination.title = 'Combination'

    # CorrelationMatrices.convert_TGraph_to_lines_and_boxes(
    #     TgCombination,
    #     drawImmediately=True,
    #     legendObject=leg,
    #     verbose=False,
    #     noBoxes=False,
    #     # xMinExternal=xMin,
    #     xMaxExternal=xMax,
    #     yMinExternal=yMin,
    #     yMaxExternal=yMax,
    #     )

    TgCombination.SetLineWidth(1)
    TgCombination.SetMarkerStyle(8)
    # TgCombination.Draw('PSAME') # Postpone to later

    yMinAbs = min( TgCombination.POIBoundsLeft[:len(expBinBoundaries)-1] )
    yMaxAbs = max( TgCombination.POIBoundsRight[:len(expBinBoundaries)-1] )


    # ======================================
    # Get a contour

    combined_rootfiles  = glob( '{0}/*.root'.format( scanDir_coupling ) )
    combined = TheoryCommands.get_TH2_from_list_of_root_files(
        combined_rootfiles, xCoupling, yCoupling, verbose = False,
        )
    combined.color = 1
    combined.name = 'regular'
    combined.title = 'Nominal'

    allcontours = TheoryCommands.get_contours_from_TH2( combined.H2, 2.30 )
    # debug('Got contours')
    candidatecontours = []

    xBestfit = combined.xBestfit
    yBestfit = combined.yBestfit
    for Tg in allcontours:
        Tg.x, Tg.y = TheoryCommands.get_xyfrom_TGraph(Tg)
        Tg.xMin = min(Tg.x)
        Tg.xMax = max(Tg.x)
        Tg.yMin = min(Tg.y)
        Tg.yMax = max(Tg.y)

        # Check if bestfit is at least inside minima and maxima
        if not (
                Tg.xMin < xBestfit
                and Tg.xMax > xBestfit
                and Tg.yMin < yBestfit
                and Tg.yMax > yBestfit
                ):
            continue

        # Compute some numbers that may help in selection: minimum distance to bestfit
        Tg.minDist = min([ ( Tg.x[i] - xBestfit )**2 + ( Tg.y[i] - yBestfit )**2 for i in xrange(Tg.GetN()) ])

        # Distance to bestfit should be close to half the difference between xMax-xMin (or yMax-yMin)
        # Compute ratio of this, minus 1 - result *should* be close to zero
        Tg.distRatio = Tg.minDist / ( 0.5 * min( Tg.yMax-Tg.yMin, Tg.xMax-Tg.xMin ) )  -  1.0
        # if abs(Tg.distRatio) > 0.8: continue

        candidatecontours.append( Tg )
    # debug('Filtered contours')

    if len(candidatecontours) == 0:
        Commands.throw_error( 'Can\'t find contour' )
    elif len(candidatecontours) > 1:
        candidatecontours.sort( key = lambda Tg: Tg.minDist )
        candidatecontours = candidatecontours[:2]

    # # Pick the contour with the highest ratio (more likely to be the 'outer shell')
    # candidatecontours.sort( key = lambda Tg: Tg.distRatio, reverse=True )

    # Actually, pick the 'inner' shell, outer shell is too likely to be a misfit
    candidatecontours.sort( key = lambda Tg: Tg.distRatio )    

    contour = candidatecontours[0]

    # ======================================
    # Determine the points to analyze:

    points = []

    # debug('Sorted contours')

    if DO_BESTFIT:
        points.append( ( xBestfit, yBestfit ) )

    if DO_X_MAXIMA:
        points.extend([
            ( contour.xMin, contour.y[ contour.x.index(contour.xMin) ] ),
            ( contour.xMax, contour.y[ contour.x.index(contour.xMax) ] ),
            ])

    if DO_Y_MAXIMA:
        points.extend([
            ( contour.x[ contour.y.index(contour.yMin) ], contour.yMin ),
            ( contour.x[ contour.y.index(contour.yMax) ], contour.yMax ),
            ])

        if getattr( container, 'OnlyYMaximaOfContour', False ):
            x1, y1 = points[-2]
            x2, y2 = points[-1]
            points.extend([
                ( xBestfit + 0.5*(x1-xBestfit), yBestfit + 0.5*(y1-yBestfit) ),
                ( xBestfit + 0.5*(x2-xBestfit), yBestfit + 0.5*(y2-yBestfit) )
                ])

    if STRAIGHT_LINE_TO_SM:
        dxSM = xBestfit - container.xSM
        dySM = yBestfit - container.ySM

        points.extend([
            ( xBestfit - 2.0*dxSM, yBestfit - 2.0*dySM ),
            ( xBestfit - 1.5*dxSM, yBestfit - 1.5*dySM ),
            ( xBestfit - 0.5*dxSM, yBestfit - 0.5*dySM ),
            ( xBestfit + 0.5*dxSM, yBestfit + 0.5*dySM ),
            # ( xBestfit - 1.5*dxSM, yBestfit - 1.5*dySM ),
            ])


    if hasattr( container, 'ManualPoints' ):
        for xRaw, yRaw in container.ManualPoints:

            if callable(xRaw):
                try:
                    x = xRaw( xBestfit )
                except TypeError:
                    x = xRaw( xBestfit, container.xSM )
            else:
                x = xRaw

            if callable(yRaw):
                try:
                    y = yRaw( yBestfit )
                except TypeError:
                    y = yRaw( yBestfit, container.ySM )
            else:
                y = yRaw

            points.append( ( x, y ) )


    # ======================================
    # Get the parametrization

    wsParametrization = WSParametrization( ws_coupling, verbose=True, newStyle=container.newStyleCoupling )

    colorCycle = new_color_cycle()
    for xPoint, yPoint in points:
        color = next(colorCycle)

        kwargs = { xCoupling : xPoint, yCoupling : yPoint }
        Tg_param = wsParametrization.get_output_container( returnWhat='exp', xMax=xMax, **kwargs ).Tg

        Tg_param.SetLineColor(color)
        Tg_param.SetMarkerColor(color)
        Tg_param.SetLineStyle(1)
        Tg_param.SetLineWidth(3)

        Tg_param.title = '{0} = {1:.1f}, {2} = {3:.1f}'.format(
            xCouplingTitle, xPoint, yCouplingTitle, yPoint
            )

        CorrelationMatrices.convert_TGraph_to_lines_and_boxes(
            Tg_param,
            drawImmediately=True,
            legendObject=leg,
            noBoxes=True,
            xMaxExternal=xMax
            )

        # For the connecting line
        # Tg_param.SetLineWidth(1)
        # Tg_param.Draw('SAMEL')

        if max(Tg_param.yValues) > yMaxAbs:
            yMaxAbs = max(Tg_param.yValues)
        if min(Tg_param.yValues) < yMinAbs:
            print 'Replacing existing yMin ({0}) with {1}, due to'.format( yMinAbs, min(Tg_param.yValues) ), ( xPoint, yPoint )
            yMinAbs = min(Tg_param.yValues)



    print '\nFound yMin = {0}, yMax = {1} after looping over TGraphs'.format( yMinAbs, yMaxAbs )

    # Draw actual combination last, so it's on top
    TgCombination.Draw('PSAME')

    yMin = yMinAbs - 0.1*(yMaxAbs-yMinAbs)
    yMax = yMaxAbs + 0.4*(yMaxAbs-yMinAbs)
    if OnOneCanvas:
        yMax = yMaxAbs + 1.2*(yMaxAbs-yMinAbs)
    base.SetMinimum(yMin)
    base.SetMaximum(yMax)


    leg.Draw()

    Commands.get_cms_label()
    Commands.get_cms_lumi()

    if not OnOneCanvas:
        save_c( '{0}_onCombination'.format(plotTitle) )


    # ======================================
    # Second plot: The contour with the points on it

    if not OnOneCanvas:
        c.Clear()
        set_cmargins( for2Dhist = True )
    else:
        cw = 1.0 - c.GetLeftMargin() - c.GetRightMargin()
        ch = 1.0 - c.GetBottomMargin() - c.GetTopMargin()
        smallPad = ROOT.TPad(
            get_unique_root_name(), '',
            c.GetLeftMargin() + 0.50*cw, c.GetBottomMargin() + 0.50*ch,
            c.GetLeftMargin() + 0.99*cw, c.GetBottomMargin() + 0.99*ch,
            )
        ROOT.SetOwnership( smallPad, False )
        smallPad.SetBottomMargin( 0.14 )
        smallPad.SetTopMargin(    0.03 )
        smallPad.SetLeftMargin(   0.12 )
        smallPad.SetRightMargin(  0.10 )
        smallPad.Draw()
        smallPad.cd()


    set_color_palette()

    combined.H2.GetXaxis().SetLabelSize(0.045)
    combined.H2.GetYaxis().SetLabelSize(0.045)

    combined.H2.GetXaxis().SetTitleSize(0.055)
    combined.H2.GetYaxis().SetTitleSize(0.055)

    # combined.H2.GetXaxis().SetTitleOffset(1.2)
    # combined.H2.GetYaxis().SetTitleOffset(1.2)

    combined.H2.GetXaxis().SetTitle( xCouplingTitle )
    combined.H2.GetYaxis().SetTitle( yCouplingTitle )

    combined.H2.Draw('COLZ')
    combined.H2.SetMaximum(7.0)

    contour.SetLineColor(1)
    contour.SetLineWidth(3)
    if OnOneCanvas: contour.SetLineWidth(2)
    contour.Draw('SAMEL')


    colorCycle = new_color_cycle()
    for xPoint, yPoint in points:
        color = next(colorCycle)

        point = ROOT.TGraph( 1, array( 'f', [xPoint] ), array( 'f', [yPoint] ) )
        ROOT.SetOwnership( point, False )

        point.SetMarkerStyle(33)
        point.SetMarkerColor(color)
        point.SetMarkerSize( 4.5 if not hasattr( container, 'MarkerSize' ) else container.MarkerSize )
        if OnOneCanvas: point.SetMarkerSize( 3.0 )
        point.Draw('SAMEP')

        point2 = ROOT.TGraph( point )
        ROOT.SetOwnership( point2, False )
        point2.SetMarkerStyle( 27 )
        point2.SetMarkerColor(1)
        point2.Draw('SAMEP')

    xMinAbs = min([ x for x,y in points ])
    xMaxAbs = max([ x for x,y in points ])
    yMinAbs = min([ y for x,y in points ])
    yMaxAbs = max([ y for x,y in points ])

    xMin = xMinAbs - 0.15*(xMaxAbs-xMinAbs)
    xMax = xMaxAbs + 0.15*(xMaxAbs-xMinAbs)
    yMin = yMinAbs - 0.15*(yMaxAbs-yMinAbs)
    yMax = yMaxAbs + 0.15*(yMaxAbs-yMinAbs)

    # combined.H2.SetMinimum( yMin )
    # combined.H2.SetMaximum( yMax )
    combined.H2.GetXaxis().SetRangeUser( xMin, xMax )
    combined.H2.GetYaxis().SetRangeUser( yMin, yMax )

    c.Update()

    if not OnOneCanvas:
        Commands.get_cms_label()
        Commands.get_cms_lumi()
        save_c( '{0}_pointsOnContour'.format(plotTitle) )
    else:
        save_c( '{0}_pointsOnContourOnePlot'.format(plotTitle) )


#____________________________________________________________________
def plot_multiple_scans(
        containers,
        xMin      = None,
        xMax      = None,
        yMin      = None,
        yMax      = None,
        xTitle    = 'x',
        yTitle    = 'y',
        plotname  = 'unnamedscans',
        draw1sigmaline = True,
        draw2sigmaline = True,
        translateToChi2 = False,
        printUncertainties = False,
        nonFancyPrint = True,
        ):

    if xMin is None: xMin = min([ min(container.x) for container in containers if not hasattr( container, 'line' ) ])
    if xMax is None: xMax = max([ max(container.x) for container in containers if not hasattr( container, 'line' ) ])
    if yMin is None: yMin = min([ min(container.y) for container in containers if not hasattr( container, 'line' ) ])
    if yMax is None: yMax = max([ max(container.y) for container in containers if not hasattr( container, 'line' ) ])

    titles = {
        'kappab' : '#kappa_{b}',
        'kappac' : '#kappa_{c}',
        'ct'     : '#kappa_{t}',
        'cg'     : '#kappa_{g}',
        'cb'     : '#kappa_{b}',
        }


    # ======================================
    # Make plot

    c.cd()
    c.Clear()
    set_cmargins( TopMargin = 0.09 )

    base = get_plot_base(
        xMin = xMin,
        xMax = xMax,
        yMin = yMin,
        yMax = yMax,
        xTitle = titles.get( xTitle, xTitle ),
        yTitle = titles.get( yTitle, yTitle ),
        )
    base.Draw('P')

    base.GetXaxis().SetTitleSize(0.06)
    base.GetXaxis().SetLabelSize(0.05)
    base.GetYaxis().SetTitleSize(0.06)
    base.GetYaxis().SetLabelSize(0.05)


    leg = ROOT.TLegend(
        c.GetLeftMargin() + 0.01,
        c.GetBottomMargin() + 0.02,
        1 - c.GetRightMargin() - 0.01,
        c.GetBottomMargin() + 0.09
        )
    leg.set_n_columns(3)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    colorCycle = itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )
    for i_container, container in enumerate(containers):

        if hasattr( container, 'line' ):
            if container.line.GetX1() == -999: container.line.SetX1( xMin )
            if container.line.GetX2() ==  999: container.line.SetX2( xMax )
            if container.line.GetY1() == -999: container.line.SetY1( yMin )
            if container.line.GetY2() ==  999: container.line.SetY2( yMax )
            container.line.Draw()
            continue

        if printUncertainties:
            container.unc = PhysicsCommands.find_minima_and_errors( container.x, container.y, returnContainer=True )
            print 'uncertainties in container',i_container
            print container.unc
            for attr in container.unc.list_attributes():
                print '  ',attr,' : ',getattr(container.unc, attr)

        if translateToChi2:
            container.y = [ 2.*y for y in container.y ]

        if not hasattr( container, 'Tg' ):
            container.Tg = ROOT.TGraph(
                len(container.x), array( 'f', container.x ), array( 'f', container.y )
                )


        if not hasattr( container, 'color' ):
            container.color = next(colorCycle)

        container.Tg.SetLineColor(container.color)
        container.Tg.SetMarkerColor(container.color)
        container.Tg.SetLineWidth(2)
        container.Tg.SetMarkerSize(0.9)
        container.Tg.SetMarkerStyle(5)

        container.Tg.Draw('PLSAME')


        if printUncertainties:

            withPercentages = False

            bestfitLine = ROOT.TLine( container.unc.min, yMin, container.unc.min, yMax )
            bestfitLine.SetLineWidth(2)
            bestfitLine.SetLineColor(container.color)
            bestfitLine.Draw()

            l = ROOT.TLatex()
            # l.SetTextColor(container.color)
            l.SetTextColor(1)
            l.SetTextSize(0.04)

            if nonFancyPrint:

                if withPercentages:
                    text = '{0} = {0:.3f}_{{{1:+.2f} ({2:+d}%)}}^{{{3:+.2f} ({4:+d}%)}}'.format(
                        titles.get( xTitle, xTitle ),
                        container.unc.min, 
                        abs(container.unc.leftError),
                        int( abs(container.unc.leftError) / container.unc.min * 100. ),
                        abs(container.unc.rightError),
                        int( abs(container.unc.rightError) / container.unc.min * 100. ),
                        )
                else:
                    # text = '{0} = {1:.3f}_{{{2:+.2f}}}^{{{3:+.2f}}}'.format(
                    #     titles.get( xTitle, xTitle ),
                    #     container.unc.min, 
                    #     abs(container.unc.leftError),
                    #     abs(container.unc.rightError),
                    #     )

                    # Decent but try 'official'
                    # text = '{1:+.2f} < {0} < {2:+.2f}'.format(
                    #     titles.get( xTitle, xTitle ),
                    #     container.unc.leftBound,
                    #     container.unc.rightBound
                    #     )
                    text = '{0} #in ({1:+.2f}, {2:+.2f}) (68%)'.format(
                        titles.get( xTitle, xTitle ),
                        container.unc.leftBound,
                        container.unc.rightBound
                        )


                l.SetTextColor(container.color)
                l.SetTextAlign(11)
                l.SetNDC()
                l.DrawLatex( 0.2, 0.8-(i_container*0.06), text )

            else:
                if withPercentages:
                    l.SetTextAlign(31)
                    l.DrawLatex(
                        container.unc.leftBound, 1.0 + 0.013*(yMax-yMin),
                        '-{0:.2f} ({1:d}%)'.format(
                            abs(container.unc.leftError),
                            int( abs(container.unc.leftError) / container.unc.min * 100. )
                            ))

                    l.SetTextAlign(11)
                    l.DrawLatex(
                        container.unc.rightBound, 1.0 + 0.013*(yMax-yMin),
                        '+{0:.2f} ({1:d}%)'.format(
                            abs(container.unc.rightError),
                            int( abs(container.unc.rightError) / container.unc.min * 100. )
                            ))
                else:
                    l.SetTextAlign(31)
                    l.DrawLatex(
                        container.unc.leftBound - 0.013*(xMax-xMin), 1.0 + 0.013*(yMax-yMin),
                        '{0:+.2f}'.format(
                            container.unc.leftBound,
                            ))
                    l.SetTextAlign(11)
                    l.DrawLatex(
                        container.unc.rightBound + 0.013*(xMax-xMin), 1.0 + 0.013*(yMax-yMin),
                        '{0:+.2f}'.format(
                            container.unc.rightBound,
                            ))

                l.SetTextAlign(21)
                l.DrawLatex(
                    container.unc.min, 1.0 + 0.013*(yMax-yMin),
                    '{0:.3f}'.format( container.unc.min )
                    )

            TgPoints = ROOT.TGraph( 2,
                array( 'f', [ container.unc.leftBound, container.unc.rightBound ] ),
                array( 'f', [ 1.0, 1.0 ] ),
                )
            TgPoints.SetMarkerSize(1.2)
            TgPoints.SetMarkerStyle(8)
            TgPoints.SetMarkerColor(container.color)
            TgPoints.Draw('PSAME')


    if draw1sigmaline:
        y1sigma = 0.5
        if translateToChi2: y1sigma = 1.0
        line1sigma = ROOT.TLine( xMin, y1sigma, xMax, y1sigma )
        line1sigma.SetLineColor(14)
        line1sigma.Draw()

    if draw2sigmaline:
        y2sigma = 1.0
        if translateToChi2: y2sigma = 2.0
        line2sigma = ROOT.TLine( xMin, y2sigma, xMax, y2sigma )
        line2sigma.SetLineColor(14)
        line2sigma.Draw()

    Commands.get_cms_label()
    Commands.get_cms_lumi()

    save_c( plotname )



#____________________________________________________________________
def filter_contour_heuristic(
        contours,
        xBestfit,
        yBestfit,
        ):

    filteredContours = []

    for Tg in contours:
        Tg.x, Tg.y = TheoryCommands.get_xyfrom_TGraph(Tg)
        Tg.xMin = min(Tg.x)
        Tg.xMax = max(Tg.x)
        Tg.yMin = min(Tg.y)
        Tg.yMax = max(Tg.y)

        # Check if bestfit is at least inside minima and maxima
        if not (
                Tg.xMin < xBestfit
                and Tg.xMax > xBestfit
                and Tg.yMin < yBestfit
                and Tg.yMax > yBestfit
                ):
            continue

        # Compute some numbers that may help in selection: minimum distance to bestfit
        Tg.minDist = min([ ( Tg.x[i] - xBestfit )**2 + ( Tg.y[i] - yBestfit )**2 for i in xrange(Tg.GetN()) ])

        # Distance to bestfit should be close to half the difference between xMax-xMin (or yMax-yMin)
        # Compute ratio of this, minus 1 - result *should* be close to zero
        Tg.distRatio = Tg.minDist / ( 0.5 * min( Tg.yMax-Tg.yMin, Tg.xMax-Tg.xMin ) )  -  1.0
        # if abs(Tg.distRatio) > 0.8: continue

        filteredContours.append( Tg )


    if len(filteredContours) == 0:
        # Heuristic failed
        filteredContours = contours
        return filteredContours
    elif len(filteredContours) > 1:
        filteredContours.sort( key = lambda Tg: Tg.minDist )
        filteredContours = filteredContours[:2]

    # Pick the contour with the highest ratio (more likely to be the 'outer shell')
    filteredContours.sort( key = lambda Tg: Tg.distRatio, reverse=True )

    return filteredContours


#____________________________________________________________________
def basic_mixed_contour_plot(
        containers,
        xMin      = 0.,
        xMax      = 1.,
        yMin      = 0.,
        yMax      = 1.,
        xTitle    = 'x',
        yTitle    = 'y',
        plotname  = 'contours',
        x_SM      = 1.,
        y_SM      = 1.,
        plotIndividualH2s = False,
        filterContours = True,
        only1sigmaContours = False,
        nLegendColumns = None,
        ):

    print '\nRunning BasicMixedContourPlot for \'{0}\''.format( plotname )

    # ======================================
    # Check whether the passed containers fulfill requirements

    for container in containers:
        attrs = container.list_attributes()


        # for expectedAttr in [ 'H2', 'name' ]:
        #     if not expectedAttr in attrs:
        #         Commands.throw_error(
        #             'Container misses mandatory attribute \'{0}\' (defined attributes: {1})'.format( expectedAttr, ', '.join(attrs) ),
        #             throwException = True
        #             )

        if not hasattr( container, 'color' ):
            container.color = 1


    # ======================================
    # Calculate contours

    for container in containers:

        if not hasattr( container, 'contours_1sigma' ):
            print '\nGetting contours for {0}'.format( container.name )
            container.contours_1sigma = TheoryCommands.get_contours_from_TH2( container.H2, 2.30 )
            container.contours_2sigma = TheoryCommands.get_contours_from_TH2( container.H2, 6.18 )

        if filterContours:
            container.contours_1sigma = filter_contour_heuristic( container.contours_1sigma, container.xBestfit, container.yBestfit )
            container.contours_2sigma = filter_contour_heuristic( container.contours_2sigma, container.xBestfit, container.yBestfit )

        if hasattr( container, 'bestfitPoint' ):
            x_root = ROOT.Double(0.)
            y_root = ROOT.Double(0.)
            container.bestfitPoint.GetPoint( 0, x_root, y_root )
            container.xBestfit = float(x_root)
            container.yBestfit = float(y_root)

        print '  Bestfit:'
        print '  {0} = {1}'.format( xTitle, container.xBestfit )
        print '  {0} = {1}'.format( yTitle, container.yBestfit )


    # ======================================
    # Make plot

    c.cd()
    c.Clear()
    set_cmargins( TopMargin=0.08 )


    base = get_plot_base(
        xMin = xMin,
        xMax = xMax,
        yMin = yMin,
        yMax = yMax,
        xTitle = xTitle,
        yTitle = yTitle,
        )
    base.Draw('P')

    base.GetXaxis().SetTitleSize(0.06)
    base.GetXaxis().SetLabelSize(0.05)
    base.GetYaxis().SetTitleSize(0.06)
    base.GetYaxis().SetLabelSize(0.05)


    leg = ROOT.TLegend(
        c.GetLeftMargin() + 0.01,
        c.GetBottomMargin() + 0.02,
        1 - c.GetRightMargin() - 0.01,
        c.GetBottomMargin() + 0.09
        )
    ROOT.SetOwnership( leg, False )
    leg.set_n_columns( min( 3, len(containers) ) )
    if not nLegendColumns is None:
        leg.set_n_columns( nLegendColumns )
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)


    TpointSM = ROOT.TGraph( 1, array( 'd', [x_SM] ), array( 'd', [y_SM] ) )
    ROOT.SetOwnership( TpointSM, False )
    TpointSM.SetMarkerSize(2)
    TpointSM.SetMarkerStyle(21)
    TpointSM.SetMarkerColor( 12 )
    TpointSM.Draw('PSAME')

    for container in containers:
        for Tg_original in container.contours_1sigma:
            Tg = Tg_original.Clone( 'contour_1sigma_' + get_unique_root_name() )
            Tg.SetLineWidth(2)
            Tg.SetLineColor( container.color )
            Tg.SetLineStyle(1)
            # Tg.SetName( 'contour_1sigma_' + get_unique_root_name() )
            Tg.Draw('CSAME')
            if Tg_original == container.contours_1sigma[0]:
                Tg.SetName( '{0}_contour_1sigma'.format(container.name) )
                leg.AddEntry(
                    Tg.GetName(),
                    container.name if not hasattr( container, 'title' ) else container.title,
                    'l' )
            ROOT.SetOwnership( Tg, False )

        if not only1sigmaContours:
            for Tg in container.contours_2sigma:
                Tg.SetLineWidth(2)
                Tg.SetLineColor( container.color )
                Tg.SetLineStyle(2)
                Tg.Draw('CSAME')

        Tpoint = ROOT.TGraph( 1, array( 'd', [container.xBestfit] ), array( 'd', [container.yBestfit] ) )
        ROOT.SetOwnership( Tpoint, False )
        Tpoint.SetMarkerSize(2)
        Tpoint.SetMarkerStyle(34)
        Tpoint.SetMarkerColor( container.color )
        Tpoint.Draw('PSAME')
        Tpoint.SetName( '{0}_bestfitpoint'.format( container.name ) )
        container.bestfitPoint = Tpoint

    leg.Draw()

    Commands.get_cms_label()
    Commands.get_cms_lumi()

    contour_dummy_legend(
        c.GetLeftMargin() + 0.01,
        1. - c.GetTopMargin() - 0.1,
        1. - c.GetRightMargin() - 0.01,
        1. - c.GetTopMargin() - 0.01,
        )

    save_c( plotname )


    if plotIndividualH2s:
        for container in containers:
            plot_single_TH2(
                container,
                xMin, xMax,
                yMin, yMax,
                xTitle, yTitle,
                plotname,                
                )


#____________________________________________________________________
def plot_single_TH2(
        container,
        xMin, xMax,
        yMin, yMax,
        xTitle, yTitle,
        plotname,
        doPNG=False, doROOT=False,
        zMax = 7.0,
        palette=None,
        getCustomContour=None,
        ):

    set_color_palette(palette)

    c.Clear()
    set_cmargins(
        # LeftMargin   = 0.12,
        RightMargin  = 0.10,
        # BottomMargin = 0.12,
        TopMargin    = 0.08,
        )

    container.H2.SetMaximum(zMax)
    container.H2.Draw('COLZ')

    # container.H2.GetXaxis().SetLimits( xMin, xMax )
    container.H2.GetXaxis().SetRangeUser( xMin, xMax )
    container.H2.GetYaxis().SetRangeUser( yMin, yMax )

    container.H2.GetXaxis().SetTitle( xTitle )
    container.H2.GetYaxis().SetTitle( yTitle )
    container.H2.GetXaxis().SetTitleSize(0.06)
    container.H2.GetXaxis().SetLabelSize(0.05)
    container.H2.GetYaxis().SetTitleSize(0.06)
    container.H2.GetYaxis().SetLabelSize(0.05)

    # if hasattr(container, 'contours_1sigma'):
    #     for Tg in container.contours_1sigma:
    #         Tg.SetName( 'contour_1sigma_' + get_unique_root_name() )
    #         Tg.Draw('CSAME')
    # if hasattr(container, 'contours_2sigma'):
    #     for Tg in container.contours_2sigma:
    #         Tg.SetName( 'contour_2sigma_' + get_unique_root_name() )
    #         Tg.Draw('CSAME')
    # if not getCustomContour is None:
    #     container.contours_custom = TheoryCommands.get_contours_from_TH2( container.H2, getCustomContour )
    #     for Tg in container.contours_custom:
    #         Tg.SetName( 'contour_custom_' + get_unique_root_name() )
    #         Tg.Draw('CSAME')

    # if hasattr( container, 'color' ):
    #     container.bestfitPoint.SetMarkerColor(container.color)

    # if hasattr(container, 'bestfitPoint'):
    #     container.bestfitPoint.Draw('PSAME')

    c.Update()
    c.RedrawAxis()

    contour_dummy_legend(
        c.GetLeftMargin() + 0.01,
        1. - c.GetTopMargin() - 0.1,
        1. - c.GetRightMargin() - 0.01,
        1. - c.GetTopMargin() - 0.01,
        )

    Commands.get_cms_label()
    Commands.get_cms_lumi()
    save_c( plotname + '_' + container.name, asPNG=doPNG, asROOT=doROOT )


#____________________________________________________________________
ContourDummyObjectsCreated = False
globalDummies = []
def contour_dummy_legend(
        x1, y1, x2, y2
        ):

    global globalDummies
    global ContourDummyObjectsCreated

    if not ContourDummyObjectsCreated:
        dummy1sigma = ROOT.TGraph( 1, array( 'f' , [-999.] ), array( 'f' , [-999.] )  )
        dummy1sigma.SetLineWidth(2)
        dummy1sigma.SetName('dummy1sigma')
        ROOT.SetOwnership( dummy1sigma, False )

        dummy2sigma = ROOT.TGraph( 1, array( 'f' , [-999.] ), array( 'f' , [-999.] )  )
        dummy2sigma.SetLineWidth(2)
        dummy2sigma.SetLineStyle(2)
        dummy2sigma.SetName('dummy2sigma')
        ROOT.SetOwnership( dummy2sigma, False )

        dummySM = ROOT.TGraph( 1, array( 'f' , [-999.] ), array( 'f' , [-999.] )  )
        dummySM.SetMarkerSize(2)
        dummySM.SetMarkerStyle(21)
        dummySM.SetMarkerColor( 12 )
        dummySM.SetName('dummySM')
        ROOT.SetOwnership( dummySM, False )

        dummybestfit = ROOT.TGraph( 1, array( 'f' , [-999.] ), array( 'f' , [-999.] )  )
        dummybestfit.SetMarkerSize(2)
        dummybestfit.SetMarkerStyle(34)
        dummybestfit.SetName('dummybestfit')
        ROOT.SetOwnership( dummybestfit, False )

        ContourDummyObjectsCreated = True

        globalDummies.extend([
            dummy1sigma,
            dummy2sigma,
            dummySM,
            dummybestfit,
            ])

    for dummy in globalDummies:
        dummy.Draw('PSAME')

    leg2 = ROOT.TLegend( x1, y1, x2, y2 )
    ROOT.SetOwnership( leg2, False )
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.set_n_columns(4)
    leg2.AddEntry( 'dummy1sigma',  '1 #sigma', 'l' )
    leg2.AddEntry( 'dummy2sigma',  '2 #sigma', 'l' )
    leg2.AddEntry( 'dummybestfit', 'Bestfit', 'p' )
    leg2.AddEntry( 'dummySM',      'SM', 'p' )
    leg2.Draw()


#____________________________________________________________________
def plot_correlation_matrix(
        container
        ):
    numpy.set_printoptions( precision=2, linewidth=100 )

    POIs = Commands.list_pois(container.ws)
    Commands.sort_pois(POIs)
    nBins = len(POIs)

    # Obtain correlation matrix
    corrMat = [ [ 999 for j in xrange(nBins) ] for i in xrange(nBins) ]
    with Commands.OpenRootFile( container.corrRootFile ) as rootFp:
        fit = rootFp.Get('fit')
        for iBin1, POI1 in enumerate( POIs ):
            for iBin2, POI2 in enumerate( POIs ):
                corrMat[iBin1][iBin2] = fit.correlation( POI1, POI2 )

    print 'Found the following corrMat from', container.corrRootFile
    print numpy.array(corrMat)
    

    # ======================================
    # Make plot

    c.Clear()
    set_cmargins(
        TopMargin   = 0.08,
        RightMargin = 0.14,
        BottomMargin = 0.17,
        )

    titleDict = {
        'PTH' : 'p_{T}^{H} (GeV)',
        }
    productionMode, observableName, _ = Commands.interpret_poi(POIs[0])
    observableName = titleDict.get( observableName, observableName )
    xTitle = getattr( container, 'xTitle', observableName )

    # Construct the binning labels
    def to_str( number ):
        if number == '-INF':
            string = '-#infty'
        elif number == 'INF':
            string = '#infty'
        elif number.is_integer():
            string = '{0:d}'.format(int(number))
        else:
            string = '{0:0.2f}'.format(number)
        return string

    binningLabels = []
    for POI in POIs:
        _1, _2, binBoundaries = Commands.interpret_poi(POI)
        binBoundariesAsStrs = [ to_str(i) for i in binBoundaries ]
        if len(binBoundaries) == 1:
            binningLabels.append( binBoundariesAsStrs[0] )
        elif len(binBoundaries) == 2:
            binningLabels.append( '(' + ', '.join(binBoundariesAsStrs) + ')'  )


    H = ROOT.TH2D(
        get_unique_root_name(),
        # '#scale[0.85]{{Bin-to-bin correlation matrix for {0}}}'.format(observableName),
        '',
        nBins, 0., nBins,
        nBins, 0., nBins
        )
    ROOT.SetOwnership( H, False )
    H.SetContour(100)

    for iRow in xrange(nBins):
        for iCol in xrange(nBins):
            H.SetBinContent( iCol+1, iRow+1, corrMat[iRow][iCol] )
    H.GetZaxis().SetRangeUser(-1.0,1.0)

    for iBin in xrange(nBins):
        H.GetXaxis().SetBinLabel( iBin+1, binningLabels[iBin] )
        H.GetYaxis().SetBinLabel( iBin+1, binningLabels[iBin] )
    H.GetXaxis().SetTitle( xTitle )
    H.GetXaxis().SetTitleSize(0.05)
    H.GetXaxis().SetTitleOffset( 1.6 )
    H.GetXaxis().SetLabelSize(0.045)
    H.GetYaxis().SetLabelSize(0.045)

    H.Draw('COLZ TEXT')


    # ======================================
    # Set some style

    ROOT.gStyle.SetHistMinimumZero() # To draw the "0", otherwise ROOT leaves it empty
    ROOT.gStyle.SetPaintTextFormat('1.2g')

    n_stops = 3
    stops  = [ 0.0, 0.5, 1.0 ]
    reds   = [ 0.0, 1.0, 1.0 ]
    blues  = [ 1.0, 1.0, 0.0 ]
    greens = [ 0.0, 1.0, 0.0 ]

    ROOT.TColor.CreateGradientColorTable(
        n_stops,
        array('d', stops ),
        array('d', reds ),
        array('d', greens ),
        array('d', blues ),
        255 )

    Commands.get_cms_label()
    Commands.get_cms_lumi()

    plotname = 'corrMat_' + basename(container.corrRootFile).replace('/','').replace('higgsCombine_','').replace('higgsCombine','').replace('.root','').replace('.','_')
    save_c( plotname )

    # Set back to default
    numpy.set_printoptions( precision=8, linewidth=75 )


########################################
# Execute this if called directly
########################################

if __name__ == '__main__':
    print 'Import functions from this file'
    print 'Available functions:'
    functions = [ name for (name, thing) in locals().items() if callable(thing) ]
    print '\n'.join( functions )