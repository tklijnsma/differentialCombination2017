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
def SetColorPalette(option=None):

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

    def SetNColumns( self, val ):
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
            leg.SetNColumns(   self._SetNColumns )
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
def PlotWithBottomPanel(
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
        GetUniqueRootName(), '',
        # c.GetLeftMargin(), topPad_bottom, 1-c.GetRightMargin(), topPad_top
        0.0, topPad_bottom, 1.0, 1.0
        )
    topPad.SetBottomMargin( topPadBottomMargin) # Distance to the bottom panel
    topPad.SetTopMargin(    topPadTopMargin)     # Space for labels
    topPad.SetLeftMargin(   topPadLeftMargin)
    topPad.SetRightMargin(  topPadRightMargin)
    topPad.Draw()

    bottomPad = ROOT.TPad(
        GetUniqueRootName(), '',
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
        Commands.GetCMSLabel( textSize=0.08 )
        Commands.GetCMSLumi( textSize=0.07 )

    SaveC(plotName)
    c.SetCanvasSize( _width, _height )


#____________________________________________________________________
def PlotSpectraOnTwoPanel(
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
        container.binBoundaries = PhysicsCommands.FigureOutBinning( container.POIs )
        
        # Determine uncertainties from scan
        container.uncs = []
        for POI, scan in zip( container.POIs, container.Scans ):
            POIvals, deltaNLLs = PhysicsCommands.FilterScan( scan )
            unc = PhysicsCommands.FindMinimaAndErrors( POIvals, deltaNLLs, returnContainer=True )
            container.uncs.append(unc)

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

    baseBottom = GetPlotBase(
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

    baseTop = GetPlotBase(
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
    leg.SetNColumns( min( len(containers)+len(hardcoded_numbers_containers), 4 ) )
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    colorCycle = Commands.newColorCycle()
    fillStyleCycle = itertools.cycle([ 3245, 3254, 3205 ])


    for container in hardcoded_numbers_containers:
        if not hasattr( container, 'color' ): container.color = next(colorCycle)
        container.nBins = getattr(container, 'nBins', len(container.binBoundaries)-1)

        Hratio = ROOT.TH1F(
            GetUniqueRootName(), '',
            len(container.binBoundaries)-1, array( 'f', container.binBoundaries )
            )
        ROOT.SetOwnership( Hratio, False )
        Hratio.SetLineColor(container.color)
        Hratio.SetLineWidth(2)
        for iBin in xrange(container.nBins):
            Hratio.SetBinContent( iBin+1, container.ratios[iBin] )
        bottomPanelObjects.append( ( Hratio, 'HISTSAME' ) )

        Hcrosssection = ROOT.TH1F(
            GetUniqueRootName(), '',
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

        if not ( container.name == 'combination' or container.name == 'combWithHbb' ):

            Hratio = ROOT.TH1F(
                GetUniqueRootName(), '',
                len(container.binBoundaries)-1, array( 'f', container.binBoundaries )
                )
            ROOT.SetOwnership( Hratio, False )
            Hratio.SetLineColor(container.color)
            Hratio.SetLineWidth(2)
            for iBin in xrange(container.nBins):
                Hratio.SetBinContent( iBin+1, container.uncs[iBin].min )
            bottomPanelObjects.append( ( Hratio, 'HISTSAME' ) )

            Hcrosssection = ROOT.TH1F(
                GetUniqueRootName(), '',
                len(container.binBoundaries)-1, array( 'f', container.binBoundaries )
                )
            ROOT.SetOwnership( Hcrosssection, False )
            Hcrosssection.SetLineColor(container.color)
            Hcrosssection.SetLineWidth(2)
            for iBin in xrange(container.nBins):
                Hcrosssection.SetBinContent( iBin+1, container.uncs[iBin].min * container.SMcrosssections[iBin] )
            topPanelObjects.append( ( Hcrosssection, 'HISTSAME' ) )


            # ---------------------
            # Construct separate TGraphAsymmErrors for the uncertainties

            Tgratio = ROOT.TGraphAsymmErrors(
                container.nBins,
                array( 'f', container.binCenters ),
                array( 'f', [ container.uncs[i].min for i in xrange(container.nBins) ] ),
                array( 'f', [ 0.45*binWidth for binWidth in container.binWidths ] ),
                array( 'f', [ 0.45*binWidth for binWidth in container.binWidths ] ),
                array( 'f', [ container.uncs[i].leftError for i in xrange(container.nBins) ] ),
                array( 'f', [ container.uncs[i].rightError for i in xrange(container.nBins) ] ),
                )
            ROOT.SetOwnership( Tgratio, False )
            Tgratio.SetName( GetUniqueRootName() )
            Tgratio.SetFillStyle(container.fillStyle)
            Tgratio.SetFillColor(container.color)
            Tgratio.SetMarkerStyle(1)
            Tgratio.SetMarkerSize(0)
            Tgratio.SetMarkerColor(container.color)
            bottomPanelObjects.append( ( Tgratio, 'E2PSAME' ) )

            Tgcrosssection = ROOT.TGraphAsymmErrors(
                container.nBins,
                array( 'f', container.binCenters ),
                array( 'f', [ container.SMcrosssections[i] * container.uncs[i].min for i in xrange(container.nBins) ] ),
                array( 'f', [ 0.45*binWidth for binWidth in container.binWidths ] ),
                array( 'f', [ 0.45*binWidth for binWidth in container.binWidths ] ),
                array( 'f', [ container.SMcrosssections[i] * container.uncs[i].leftError for i in xrange(container.nBins) ] ),
                array( 'f', [ container.SMcrosssections[i] * container.uncs[i].rightError for i in xrange(container.nBins) ] ),
                )
            ROOT.SetOwnership( Tgcrosssection, False )
            Tgcrosssection.SetName( GetUniqueRootName() )
            Tgcrosssection.SetFillStyle(container.fillStyle)
            Tgcrosssection.SetFillColor(container.color)
            Tgcrosssection.SetMarkerStyle(1)
            Tgcrosssection.SetMarkerSize(0)
            Tgcrosssection.SetMarkerColor(container.color)
            Tgcrosssection.SetLineColor(container.color)
            Tgcrosssection.SetLineWidth(2)
            topPanelObjects.append( ( Tgcrosssection, 'E2PSAME' ) )

            leg.AddEntry( Tgcrosssection.GetName(), container.title, 'LF' )


        else:


            Tgratio = ROOT.TGraphAsymmErrors(
                container.nBins,
                array( 'f', container.binCenters ),
                array( 'f', [ container.uncs[i].min for i in xrange(container.nBins) ] ),
                array( 'f', [ 0.0 for i in xrange(container.nBins) ] ),
                array( 'f', [ 0.0 for i in xrange(container.nBins) ] ),
                array( 'f', [ container.uncs[i].leftError for i in xrange(container.nBins) ] ),
                array( 'f', [ container.uncs[i].rightError for i in xrange(container.nBins) ] ),
                )
            ROOT.SetOwnership( Tgratio, False )
            Tgratio.SetName( GetUniqueRootName() )
            Tgratio.SetMarkerStyle(8)
            Tgratio.SetMarkerSize(1.2)
            # Tgratio.SetLineWidth(2)
            Tgratio.SetMarkerColor(container.color)
            Tgratio.SetLineColor(container.color)
            bottomPanelObjects.append( ( Tgratio, 'PSAME' ) )

            Tgcrosssection = ROOT.TGraphAsymmErrors(
                container.nBins,
                array( 'f', container.binCenters ),
                array( 'f', [ container.SMcrosssections[i] * container.uncs[i].min for i in xrange(container.nBins) ] ),
                array( 'f', [ 0.0 for i in xrange(container.nBins) ] ),
                array( 'f', [ 0.0 for i in xrange(container.nBins) ] ),
                array( 'f', [ container.SMcrosssections[i] * container.uncs[i].leftError for i in xrange(container.nBins) ] ),
                array( 'f', [ container.SMcrosssections[i] * container.uncs[i].rightError for i in xrange(container.nBins) ] ),
                )
            ROOT.SetOwnership( Tgcrosssection, False )
            Tgcrosssection.SetName( GetUniqueRootName() )
            Tgcrosssection.SetMarkerStyle(8)
            Tgcrosssection.SetMarkerSize(1.2)
            Tgcrosssection.SetMarkerColor(container.color)
            Tgcrosssection.SetLineColor(container.color)
            # Tgcrosssection.SetLineWidth(2)
            topPanelObjects.append( ( Tgcrosssection, 'PSAME' ) )

            leg.AddEntry( Tgcrosssection.GetName(), container.title, 'PE' )


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
    PlotWithBottomPanel(
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
def WriteScansToTable(
        container,
        tag,
        xTitle       = 'p_{T}^{H}',
        yTitle       = '#Delta#sigma (pb/GeV)',
        lastBinIsOverflow = True,
        verbose = True,
        scaleLastBin = False,
        ):

    container.binBoundaries = PhysicsCommands.FigureOutBinning( container.POIs )
    container.nBins = len(container.binBoundaries)-1
    
    # Determine uncertainties from scan
    container.uncs = []
    for POI, scan in zip( container.POIs, container.Scans ):
        POIvals, deltaNLLs = PhysicsCommands.FilterScan( scan )
        unc = PhysicsCommands.FindMinimaAndErrors( POIvals, deltaNLLs, returnContainer=True )
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

    tableText = '% ' + Commands.TagGitCommitAndModule() + '\n'
    tableText += Commands.PrintTable(table, maxColWidth=100, sep=' & ', newline_sep=' \\\\\n' )

    outname = join(TheoryCommands.PLOTDIR, 'mutable_{0}_{1}.tex'.format(tag, container.name))
    with open(outname, 'w') as outFp:
        outFp.write(tableText)


#____________________________________________________________________
def PlotParametrizationsOnCombination(
        container,
        OnOneCanvas = False
        ):

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
            Commands.ThrowError( 'Requested mode \'StraightLineToSM\', but xSM is not specified' )
        if not hasattr( container, 'ySM' ):
            Commands.ThrowError( 'Requested mode \'StraightLineToSM\', but ySM is not specified' )


    # ======================================
    # Begin plotting

    newColorCycle = lambda: itertools.cycle( [ 2, 4, 6, 41, 46, 30, 43, 3, 5, 8, 9 ] )

    c.Clear()
    SetCMargins( TopMargin = 0.09 )

    xMin = expBinBoundaries[0]
    xMax = expBinBoundaries[-2] + ( expBinBoundaries[-2] - expBinBoundaries[-3] ) # Overflow will screw up plot
    yMin = 0.0
    yMax = 1.0

    base = GetPlotBase(
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
        leg.SetNColumns(1)
    else:
        leg = ROOT.TLegend(
            c.GetLeftMargin(),
            1 - c.GetTopMargin() - 0.17,
            1 - c.GetRightMargin(),
            1 - c.GetTopMargin() 
            )
        leg.SetNColumns(2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)


    # ======================================
    # Draw the combination as blocks

    combinationPOIs = Commands.ListPOIs( ws_combination )
    combinationscans = PhysicsCommands.GetScanResults( combinationPOIs, scanDir_combination, pattern = '' )

    TgCombination = PhysicsCommands.GetTGraphForSpectrum(
        combinationPOIs,
        combinationscans,
        name = GetUniqueRootName()
        )

    TgCombination.SetLineColor( 1 )
    # TgCombination.SetMarkerStyle( 2 )
    TgCombination.SetFillColorAlpha( 1, 0.2 )
    # TgCombination.SetFillColor( 13 )
    # TgCombination.SetFillStyle( 3544 )
    # TgCombination.SetFillStyle( 3345 )

    TgCombination.title = 'Combination'

    # CorrelationMatrices.ConvertTGraphToLinesAndBoxes(
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
    combined = TheoryCommands.GetTH2FromListOfRootFiles(
        combined_rootfiles, xCoupling, yCoupling, verbose = False,
        )
    combined.color = 1
    combined.name = 'regular'
    combined.title = 'Nominal'

    allcontours = TheoryCommands.GetContoursFromTH2( combined.H2, 2.30 )
    candidatecontours = []

    xBestfit = combined.xBestfit
    yBestfit = combined.yBestfit
    for Tg in allcontours:
        Tg.x, Tg.y = TheoryCommands.GetXYfromTGraph(Tg)
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


    if len(candidatecontours) == 0:
        Commands.ThrowError( 'Can\'t find contour' )
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
                y = xRaw


            points.append( ( x, y ) )


    # ======================================
    # Get the parametrization

    wsParametrization = WSParametrization( ws_coupling )

    colorCycle = newColorCycle()
    for xPoint, yPoint in points:
        color = next(colorCycle)

        kwargs = { xCoupling : xPoint, yCoupling : yPoint }
        Tg_param = wsParametrization.GetOutputContainer( returnWhat='exp', xMax=xMax, **kwargs ).Tg

        Tg_param.SetLineColor(color)
        Tg_param.SetMarkerColor(color)
        Tg_param.SetLineStyle(1)
        Tg_param.SetLineWidth(3)

        Tg_param.title = '{0} = {1:.1f}, {2} = {3:.1f}'.format(
            xCouplingTitle, xPoint, yCouplingTitle, yPoint
            )

        CorrelationMatrices.ConvertTGraphToLinesAndBoxes(
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

    Commands.GetCMSLabel()
    Commands.GetCMSLumi()

    if not OnOneCanvas:
        SaveC( '{0}_onCombination'.format(plotTitle) )


    # ======================================
    # Second plot: The contour with the points on it

    if not OnOneCanvas:
        c.Clear()
        SetCMargins( for2Dhist = True )
    else:
        cw = 1.0 - c.GetLeftMargin() - c.GetRightMargin()
        ch = 1.0 - c.GetBottomMargin() - c.GetTopMargin()
        smallPad = ROOT.TPad(
            GetUniqueRootName(), '',
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


    SetColorPalette()

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


    colorCycle = newColorCycle()
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
        Commands.GetCMSLabel()
        Commands.GetCMSLumi()
        SaveC( '{0}_pointsOnContour'.format(plotTitle) )
    else:
        SaveC( '{0}_pointsOnContourOnePlot'.format(plotTitle) )


#____________________________________________________________________
def PlotMultipleScans(
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
    SetCMargins( TopMargin = 0.09 )

    base = GetPlotBase(
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
    leg.SetNColumns(3)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    colorCycle = itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )
    for container in containers:

        if hasattr( container, 'line' ):
            if container.line.GetX1() == -999: container.line.SetX1( xMin )
            if container.line.GetX2() ==  999: container.line.SetX2( xMax )
            if container.line.GetY1() == -999: container.line.SetY1( yMin )
            if container.line.GetY2() ==  999: container.line.SetY2( yMax )
            container.line.Draw()
            continue

        if printUncertainties:
            container.unc = PhysicsCommands.FindMinimaAndErrors( container.x, container.y, returnContainer=True )

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


                l.SetTextAlign(11)
                l.SetNDC()
                l.DrawLatex( 0.2, 0.8, text )

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

    Commands.GetCMSLabel()
    Commands.GetCMSLumi()

    SaveC( plotname )



#____________________________________________________________________
def FilterContourHeuristic(
        contours,
        xBestfit,
        yBestfit,
        ):

    filteredContours = []

    for Tg in contours:
        Tg.x, Tg.y = TheoryCommands.GetXYfromTGraph(Tg)
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
def BasicMixedContourPlot(
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
        attrs = container.ListAttributes()


        # for expectedAttr in [ 'H2', 'name' ]:
        #     if not expectedAttr in attrs:
        #         Commands.ThrowError(
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
            container.contours_1sigma = TheoryCommands.GetContoursFromTH2( container.H2, 2.30 )
            container.contours_2sigma = TheoryCommands.GetContoursFromTH2( container.H2, 6.18 )

        if filterContours:
            container.contours_1sigma = FilterContourHeuristic( container.contours_1sigma, container.xBestfit, container.yBestfit )
            container.contours_2sigma = FilterContourHeuristic( container.contours_2sigma, container.xBestfit, container.yBestfit )

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
    SetCMargins( TopMargin=0.08 )


    base = GetPlotBase(
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
    leg.SetNColumns( min( 3, len(containers) ) )
    if not nLegendColumns is None:
        leg.SetNColumns( nLegendColumns )
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
            Tg = Tg_original.Clone( 'contour_1sigma_' + GetUniqueRootName() )
            Tg.SetLineWidth(2)
            Tg.SetLineColor( container.color )
            Tg.SetLineStyle(1)
            # Tg.SetName( 'contour_1sigma_' + GetUniqueRootName() )
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

    Commands.GetCMSLabel()
    Commands.GetCMSLumi()

    ContourDummyLegend(
        c.GetLeftMargin() + 0.01,
        1. - c.GetTopMargin() - 0.1,
        1. - c.GetRightMargin() - 0.01,
        1. - c.GetTopMargin() - 0.01,
        )

    SaveC( plotname )


    if plotIndividualH2s:
        for container in containers:
            PlotSingle2DHistogram(
                container,
                xMin, xMax,
                yMin, yMax,
                xTitle, yTitle,
                plotname,                
                )


#____________________________________________________________________
def PlotSingle2DHistogram(
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

    SetColorPalette(palette)

    c.Clear()
    SetCMargins(
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

    if hasattr(container, 'contours_1sigma'):
        for Tg in container.contours_1sigma:
            Tg.SetName( 'contour_1sigma_' + GetUniqueRootName() )
            Tg.Draw('CSAME')
    if hasattr(container, 'contours_2sigma'):
        for Tg in container.contours_2sigma:
            Tg.SetName( 'contour_2sigma_' + GetUniqueRootName() )
            Tg.Draw('CSAME')
    if not getCustomContour is None:
        container.contours_custom = TheoryCommands.GetContoursFromTH2( container.H2, getCustomContour )
        for Tg in container.contours_custom:
            Tg.SetName( 'contour_custom_' + GetUniqueRootName() )
            Tg.Draw('CSAME')

    if hasattr( container, 'color' ):
        container.bestfitPoint.SetMarkerColor(container.color)

    if hasattr(container, 'bestfitPoint'):
        container.bestfitPoint.Draw('PSAME')

    c.Update()
    c.RedrawAxis()

    ContourDummyLegend(
        c.GetLeftMargin() + 0.01,
        1. - c.GetTopMargin() - 0.1,
        1. - c.GetRightMargin() - 0.01,
        1. - c.GetTopMargin() - 0.01,
        )

    Commands.GetCMSLabel()
    Commands.GetCMSLumi()
    SaveC( plotname + '_' + container.name, asPNG=doPNG, asROOT=doROOT )


#____________________________________________________________________
ContourDummyObjectsCreated = False
globalDummies = []
def ContourDummyLegend(
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
    leg2.SetNColumns(4)
    leg2.AddEntry( 'dummy1sigma',  '1 #sigma', 'l' )
    leg2.AddEntry( 'dummy2sigma',  '2 #sigma', 'l' )
    leg2.AddEntry( 'dummybestfit', 'Bestfit', 'p' )
    leg2.AddEntry( 'dummySM',      'SM', 'p' )
    leg2.Draw()


#____________________________________________________________________
def PlotCorrelationMatrix(
        container
        ):
    numpy.set_printoptions( precision=2, linewidth=100 )

    POIs = Commands.ListPOIs(container.ws)
    Commands.SortPOIs(POIs)
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
    SetCMargins(
        TopMargin   = 0.08,
        RightMargin = 0.14,
        BottomMargin = 0.17,
        )

    titleDict = {
        'PTH' : 'p_{T}^{H} (GeV)',
        }
    productionMode, observableName, _ = Commands.InterpretPOI(POIs[0])
    observableName = titleDict.get( observableName, observableName )
    xTitle = getattr( container, 'xTitle', observableName )

    # Construct the binning labels
    def toStr( number ):
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
        _1, _2, binBoundaries = Commands.InterpretPOI(POI)
        binBoundariesAsStrs = [ toStr(i) for i in binBoundaries ]
        if len(binBoundaries) == 1:
            binningLabels.append( binBoundariesAsStrs[0] )
        elif len(binBoundaries) == 2:
            binningLabels.append( '(' + ', '.join(binBoundariesAsStrs) + ')'  )


    H = ROOT.TH2D(
        GetUniqueRootName(),
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

    Commands.GetCMSLabel()
    Commands.GetCMSLumi()

    plotname = 'corrMat_' + basename(container.corrRootFile).replace('/','').replace('higgsCombine_','').replace('higgsCombine','').replace('.root','').replace('.','_')
    SaveC( plotname )

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