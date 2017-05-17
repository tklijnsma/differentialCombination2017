#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import Commands

import os, tempfile, shutil, re, glob, itertools, sys
from os.path import *
from operator import itemgetter
from array import array

from time import strftime
datestr = strftime( '%b%d' )

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
ROOT.gStyle.SetOptStat(0)
c = ROOT.TCanvas( 'cpc', 'cpc', 1000, 800 )

LeftMargin   = 0.15
RightMargin  = 0.03
BottomMargin = 0.15
TopMargin    = 0.03
c.SetLeftMargin( LeftMargin )
c.SetRightMargin( RightMargin )
c.SetBottomMargin( BottomMargin )
c.SetTopMargin( TopMargin )


########################################
# Main
########################################


PLOTDIR = 'plots_{0}'.format(datestr)
def SaveC( outname, PNG=False ):
    global PLOTDIR
    if not isdir(PLOTDIR): os.makedirs(PLOTDIR)
    outname = join( PLOTDIR, basename(outname).replace('.pdf','').replace('.png','') )
    c.SaveAs( outname + '.pdf' )
    if PNG:
        c.SaveAs( outname + '.png' )




def GetScanResults(
    POIs,
    scanDirectory,
    pattern='*',
    verbose=False,
    ):

    POIs.sort( key=lambda POI: Commands.InterpretPOI( POI )[2][0] )

    scans = []
    for POI in POIs:
        prodMode, obsName, obsRange = Commands.InterpretPOI( POI )

        if pattern == '*':
            rootfilepattern = '*{0}*.root'.format(POI)
        else:
            rootfilepattern = '*{0}*{1}*.root'.format( POI, pattern )

        scanFiles = glob.glob( join( scanDirectory, rootfilepattern ) )
        if len(scanFiles) == 0:
            Commands.ThrowError( 'The pattern \'{0}\' yielded no results in directory {1}'.format( rootfilepattern, scanDirectory ) )
            return

        scanResult = Commands.ConvertTChainToArray(
            'limit',
            scanFiles,
            r'{0}|deltaNLL'.format( POI )
            )

        scanResult = [ (POIval, deltaNLL) for (POIval, deltaNLL) in sorted(zip( scanResult[POI], scanResult['deltaNLL'] ))]


        if verbose:
            print '\nPOI: {0}'.format(POI)
            print '{0:10}  |  {1:10}'.format( 'POIval', 'deltaNLL' )
            for POIval, deltaNLL in scanResult:
                print '{0:+10.4f}  |  {1:+10.4f}'.format( POIval, deltaNLL )

        scans.append(scanResult)

    return scans



def FilterScan( scan ):
    scan = list(set( scan ))
    scan.sort()
    POIvals   = []
    deltaNLLs = []
    for POIval, deltaNLL in scan:
        if (
            POIval > -1000. and POIval < 1000. and
            deltaNLL > -1000. and deltaNLL < 1000. ):
            POIvals.append( POIval )
            deltaNLLs.append( deltaNLL )
    return POIvals, deltaNLLs




def BasicDrawScanResults(
    POIs,
    scans,
    POIRange = '*',
    deltaNLLRange = '*',
    name = 'unnamed'
    ):

    # Use first POI to get production mode and observable name
    productionMode, observableName, dummyRange = Commands.InterpretPOI( POIs[0] )
    
    if POIRange == '*':
        allPOIvals = map( itemgetter(0), itertools.chain.from_iterable( scans ) )
        POIRange = [
            min( filter( lambda x: x < 1000., allPOIvals ) ),
            max( filter( lambda x: x < 1000., allPOIvals ) ),
            ]
    # print POIRange

    if deltaNLLRange == '*':
        allDeltaNLLs = map( itemgetter(1), itertools.chain.from_iterable( scans ) )
        deltaNLLRange = [
            min( filter( lambda x: x < 1000., allDeltaNLLs ) ),
            max( filter( lambda x: x < 1000., allDeltaNLLs ) ),
            ]
    # print deltaNLLRange

    c.Clear()

    base = ROOT.TH1F()
    base.Draw('P')
    base.GetXaxis().SetLimits( POIRange[0], POIRange[1] )
    base.SetMinimum( deltaNLLRange[0] )
    base.SetMaximum( deltaNLLRange[1] )
    base.SetMarkerColor(0)
    base.GetXaxis().SetTitle( '#mu' )
    base.GetYaxis().SetTitle( '#Delta NLL' )

    leg = ROOT.TLegend( 0.65, 0.5, 1-RightMargin, 1-TopMargin )
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)


    # Create predefined list of colors
    colors = range(2,9+1) + [ 30, 38, 40, 41, 42 ] + range( 45, 48+1 )
    colors = colors + colors + colors + colors
    iColor = 0

    for POI, scan in zip( POIs, scans ):

        POIvals, deltaNLLs = FilterScan( scan )

        nPoints = len(POIvals)

        # print '\n\n'
        # for POIval, deltaNLL in zip(POIvals, deltaNLLs):
        #     print POIval, deltaNLL


        Tg = ROOT.TGraph(
            nPoints,
            array( 'd', POIvals ),
            array( 'd', deltaNLLs ),
            )
        ROOT.SetOwnership( Tg, False )
        Tg.SetName( 'Tg_' + POI )

        color = colors[iColor]
        iColor += 1

        Tg.SetLineColor(color)
        Tg.SetMarkerColor(color)
        Tg.SetLineWidth(2)
        Tg.SetMarkerStyle(5)
        Tg.Draw('SAME L')

        # Small line at minimum
        minimaAndUncertainties = FindMinimaAndErrors( POIvals, deltaNLLs )
        lMin = ROOT.TLine( minimaAndUncertainties['min'], 0., minimaAndUncertainties['min'], 2.0 )
        ROOT.SetOwnership( lMin, False )
        lMin.SetLineColor( color )
        lMin.SetLineWidth(2)
        lMin.Draw()

        leg.AddEntry( Tg.GetName(), POI, 'l' )

    leg.Draw()

    SaveC( 'parabolas_{0}_{1}_{2}'.format( productionMode, observableName, name ) )



def FigureOutBinning( POIs ):
    POIs.sort( key=lambda POI: Commands.InterpretPOI( POI )[2][0] )
    Ranges = [ Commands.InterpretPOI( POI )[2] for POI in POIs ]

    rangeLengthList = list(set(map( len, Ranges )))
    if len(rangeLengthList) != 1:
        print 'ERROR: Ranges are inconsistent:'
        print Ranges
        return
    rangeLength = rangeLengthList[0]

    binBoundaries = []

    if rangeLength == 2:

        for iRange, Range in enumerate(Ranges):
            leftBound, rightBound = Range

            if leftBound == '-INF':
                nextLeftBound, nextRightBound = map( float, Ranges[iRange+1] )
                leftBound = rightBound - ( nextRightBound - nextLeftBound )

            if rightBound == 'INF':
                previousLeftBound, previousRightBound = map( float, Ranges[iRange-1] )
                rightBound = leftBound + ( previousRightBound - previousLeftBound )

            binBoundaries.append( leftBound )
            if iRange == len(Ranges)-1:
                binBoundaries.append( rightBound )

    elif rangeLength == 1:
        print 'ERROR: Bin observable range has length 1, not yet implemented'
        return
    else:
        print 'ERROR: Bin observable range has length {0}, which does not make sense'.format( rangeLength )
        return

    return binBoundaries



def GetTGraphForSpectrum(
    POIs,
    scans,
    name = 'unnamed'
    ):

    # Sort by first number of observable range
    POIs.sort( key=lambda POI: Commands.InterpretPOI( POI )[2][0] )
    binBoundaries = FigureOutBinning( POIs )

    binCenters    = [ 0.5*(left+right) for left, right in zip( binBoundaries[:-1], binBoundaries[1:] ) ]
    binWidths     = [ right-left for left, right in zip( binBoundaries[:-1], binBoundaries[1:] ) ]
    halfBinWidths = [ 0.5*(right-left) for left, right in zip( binBoundaries[:-1], binBoundaries[1:] ) ]

    POICenters  = []
    POIErrsLeft  = []
    POIErrsRight = []
    for POI, scan in zip( POIs, scans ):

        POIvals, deltaNLLs = FilterScan( scan )
        nPoints = len(POIvals)
        minimaAndUncertainties = FindMinimaAndErrors( POIvals, deltaNLLs )

        POICenters.append( minimaAndUncertainties['min'] )
        POIErrsLeft.append( minimaAndUncertainties['leftError'] )
        POIErrsRight.append( minimaAndUncertainties['rightError'] )


    Tg = ROOT.TGraphAsymmErrors(
        len(POIs),
        array( 'd', binCenters ),
        array( 'd', POICenters ),
        array( 'd', halfBinWidths ),
        array( 'd', halfBinWidths ),
        array( 'd', POIErrsLeft ),
        array( 'd', POIErrsRight ),
        )
    ROOT.SetOwnership( Tg, False )

    Tg.SetName( 'Tg_' + name )
    Tg.SetLineWidth(2)
    Tg.SetFillStyle(0)


    # Add some more information to the object

    Tg.name = name
    Tg.xMin = binBoundaries[0]
    Tg.xMax = binBoundaries[-1]
    Tg.yMin = min( [ POICenter - POIErrLeft for POICenter, POIErrLeft in zip( POICenters, POIErrsLeft ) ] )
    Tg.yMax = max( [ POICenter + POIErrRight for POICenter, POIErrRight in zip( POICenters, POIErrsRight ) ] )

    Tg.binBoundaries = binBoundaries
    Tg.binCenters    = binCenters
    Tg.binWidths     = binWidths
    Tg.halfBinWidths = halfBinWidths
    Tg.POICenters    = POICenters
    Tg.POIErrsLeft   = POIErrsLeft
    Tg.POIErrsRight  = POIErrsRight

    return Tg



def BasicDrawSpectrum(
    POIs,
    scans,
    name = 'unnamed'
    ):

    # Use first POI to get production mode and observable name
    productionMode, observableName, dummyRange = Commands.InterpretPOI( POIs[0] )

    c.Clear()

    Tg = GetTGraphForSpectrum( POIs, scans, name = 'unnamed' )

    base = ROOT.TH1F()
    base.Draw('P')
    base.GetXaxis().SetLimits( Tg.xMin, Tg.xMax )
    base.SetMinimum( Tg.yMin )
    base.SetMaximum( Tg.yMax )
    base.SetMarkerColor(0)
    base.GetXaxis().SetTitle( observableName )
    base.GetYaxis().SetTitle( '#mu' )

    # Create predefined list of colors
    colors = range(2,9+1) + [ 30, 38, 40, 41, 42 ] + range( 45, 48+1 )
    colors = colors + colors + colors + colors
    iColor = 0
    color = colors[iColor]

    Tg.Draw('SAME P')
    Tg.SetLineColor(color)
    Tg.SetMarkerColor(color)

    SaveC( 'spectrum_{0}_{1}_{2}'.format( productionMode, observableName, name ) )




def BasicCombineSpectra(
    *args,
    **kwargs
    ):

    names       = []
    POIs        = []
    scans       = []
    drawoptions = []

    for arg in args:
        names.append( arg[0] )
        POIs.append(  arg[1] )
        scans.append( arg[2] )
        if len(arg) > 3:
            drawoptions.append( arg[3:] )
        else:
            drawoptions.append( [] )


    Tgs = []
    for name, POIList, scanList, drawoptionList in zip( names, POIs, scans, drawoptions ):
        Tg = GetTGraphForSpectrum( POIList, scanList, name=name )

        for drawoption in drawoptionList:
            try:
                getattr( Tg, drawoption[0] )( *drawoption[1:] )
                # print 'Successfully called \'Tg.{0}({1})\''.format( drawoption[0], drawoption[1:] )
            except:
                Commands.ThrowError( 'Problem calling \'Tg.{0}({1})\', skipping'.format( drawoption[0], drawoption[1:] ) )

        Tgs.append( Tg )



    # Use first POI to get production mode and observable name
    productionMode, observableName, dummyRange = Commands.InterpretPOI( POIs[0][0] )


    xMin = min([ Tg.xMin for Tg in Tgs ])
    xMax = max([ Tg.xMax for Tg in Tgs ])
    dx = xMax-xMin
    yMin = min([ Tg.yMin for Tg in Tgs ])
    yMax = max([ Tg.yMax for Tg in Tgs ])
    dy = yMax-yMin
    
    c.Clear()
    base = ROOT.TH1F()
    base.Draw('P')
    base.GetXaxis().SetLimits( xMin, xMax )
    base.SetMinimum( yMin - 0.05*dy )
    base.SetMaximum( yMax + 0.05*dy )
    base.SetMarkerColor(0)
    base.GetXaxis().SetTitle( observableName )
    base.GetYaxis().SetTitle( '#mu' )

    lineAtOne = ROOT.TLine( xMin, 1.0, xMax, 1.0 )
    lineAtOne.SetLineWidth(3)
    lineAtOne.SetLineColor(17)
    # lineAtOne.SetLineStyle(2)
    lineAtOne.Draw()

    leg = ROOT.TLegend( 1-RightMargin-0.25, 1-TopMargin-0.3, 1-RightMargin, 1-TopMargin )
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)

    optionTLines = True
    TgDrawStr = 'SAME E2 P'

    for Tg in Tgs:

        # Make overflow bin extend until max of spectrum
        if Tg.binBoundaries[-1] < xMax:
            Tg.binBoundaries[-1] = xMax
            newBinCenter    = 0.5*( xMax + Tg.binBoundaries[-2] )
            newHalfBinWidth = 0.5*( xMax - Tg.binBoundaries[-2] )
            lastPoint = len(Tg.binCenters)-1
            Tg.SetPoint( lastPoint, newBinCenter, Tg.POICenters[lastPoint] )
            Tg.SetPointEXlow( lastPoint, newHalfBinWidth )
            Tg.SetPointEXhigh( lastPoint, newHalfBinWidth )


        if optionTLines:

            # This is already transparent!
            # Regular color 4 becomes something like 972 when transparent
            fillColor = Tg.GetFillColor()
            fillStyle = Tg.GetFillStyle()
            lineColor = Tg.GetLineColor()

            for iBin in xrange(len(Tg.binBoundaries)-1):

                box = ROOT.TBox(
                    Tg.binBoundaries[iBin], Tg.POICenters[iBin]-Tg.POIErrsLeft[iBin],
                    Tg.binBoundaries[iBin+1], Tg.POICenters[iBin]+Tg.POIErrsRight[iBin]
                    )
                ROOT.SetOwnership( box, False )

                box.SetLineWidth(0)
                box.SetFillColor( fillColor )
                if not fillStyle == 0: box.SetFillStyle( fillStyle )
                box.Draw()

                line = ROOT.TLine(
                    Tg.binBoundaries[iBin], Tg.POICenters[iBin],
                    Tg.binBoundaries[iBin+1], Tg.POICenters[iBin],
                    )
                ROOT.SetOwnership( line, False )

                line.SetLineColor( lineColor )
                line.SetLineWidth( Tg.GetLineWidth() )
                line.Draw()

            # Dummy for the legend
            TgDummy = ROOT.TGraph( 1, array( 'd', [-999.]), array( 'd', [-999.]) )
            ROOT.SetOwnership( TgDummy, False )
            TgDummy.SetLineColor( lineColor )
            TgDummy.SetLineWidth( Tg.GetLineWidth() )
            TgDummy.SetFillColor( fillColor )
            if not fillStyle == 0: TgDummy.SetFillStyle( fillStyle )
            TgDummy.SetMarkerStyle(8)
            TgDummy.SetMarkerColor( lineColor )
            TgDummy.SetMarkerSize(0)
            TgDummy.SetName( Tg.name + '_dummy' )
            TgDummy.Draw('SAME')

            leg.AddEntry( TgDummy.GetName(), Tg.name, 'lf' )


        else:
            Tg.Draw( TgDrawStr )

    leg.Draw()


    SaveC( 'mspectrum_{0}_{1}_{2}'.format( productionMode, observableName, '_'.join(names) ) )



def BasicDrawMultipleParabolas(
    *args,
    **kwargs
    ):

    # ======================================
    # Handle input

    names       = []
    POIs        = []
    scans       = []
    drawoptions = []

    for arg in args:
        names.append( arg[0] )
        POIs.append(  arg[1] )
        scans.append( arg[2] )
        if len(arg) > 3:
            drawoptions.append( arg[3:] )
        else:
            drawoptions.append( [] )


    # ======================================
    # Determine the ranges of the plot

    POImin      = 999.
    POImax      = -999.
    deltaNLLmin = 999.
    deltaNLLmax = -999.
    for scanList in scans:
        allPOIvals = map( itemgetter(0), itertools.chain.from_iterable( scanList ) )
        this_POImin = min( filter( lambda x: x < 1000., allPOIvals ) )
        this_POImax = max( filter( lambda x: x < 1000., allPOIvals ) )
        allDeltaNLLs = map( itemgetter(1), itertools.chain.from_iterable( scanList ) )
        this_deltaNLLmin = min( filter( lambda x: x < 1000., allDeltaNLLs ) )
        this_deltaNLLmax = max( filter( lambda x: x < 1000., allDeltaNLLs ) )
        if this_POImin < POImin: POImin = this_POImin
        if this_POImax > POImax: POImax = this_POImax
        if this_deltaNLLmin < deltaNLLmin: deltaNLLmin = this_deltaNLLmin
        if this_deltaNLLmax > deltaNLLmax: deltaNLLmax = this_deltaNLLmax


    # ======================================
    # Make TGraph objects for ALL parabolas

    Tgs = { name : {} for name in names }
    for name, POIList, scanList, drawoptionList in zip( names, POIs, scans, drawoptions ):

        for POI, scan in zip( POIList, scanList ):

            POIvals, deltaNLLs = FilterScan( scan )
            nPoints = len(POIvals)

            Tg = ROOT.TGraph(
                nPoints,
                array( 'd', POIvals ),
                array( 'd', deltaNLLs ),
                )
            ROOT.SetOwnership( Tg, False )
            Tg.SetName( name + '_' + POI )
            Tg.SetLineWidth(2)

            for drawoption in drawoptionList:
                try:
                    getattr( Tg, drawoption[0] )( *drawoption[1:] )
                    # print 'Successfully called \'Tg.{0}({1})\''.format( drawoption[0], drawoption[1:] )
                except:
                    Commands.ThrowError( 'Problem calling \'Tg.{0}({1})\', skipping'.format( drawoption[0], drawoption[1:] ) )

            # For small line at minimum
            minimaAndUncertainties = FindMinimaAndErrors( POIvals, deltaNLLs )
            Tg.imin        = minimaAndUncertainties['imin']
            Tg.min         = minimaAndUncertainties['min']
            Tg.leftError   = minimaAndUncertainties['leftError']
            Tg.leftBound   = minimaAndUncertainties['leftBound']
            Tg.rightError  = minimaAndUncertainties['rightError']
            Tg.rightBound  = minimaAndUncertainties['rightBound']

            Tg.channel = name
            Tg.POI     = POI

            Tgs[name][POI] = Tg


    # ======================================
    # Figure out the binning from the POIs, and determine which scans to plot in where

    nBins = max( map( len, POIs ) )
    finestPOIs = [ POIList for POIList in POIs if len(POIList)==nBins ][0]
    finestBoundaries = FigureOutBinning( finestPOIs )
    binBoundaries = { name : FigureOutBinning( POIList ) for name, POIList in zip( names, POIs ) }

    parabolasPerBin = []

    for iBin in xrange(nBins):

        POI = finestPOIs[iBin]
        observableMin = finestBoundaries[iBin]
        observableMax = finestBoundaries[iBin+1]

        parabolas = []

        # Loop over the channels, and find the parabola that overlaps with these ranges
        for name, POIList in zip( names, POIs ):

            # If channel is already in finest binning, it's easy
            if len(POIList) == nBins and POI == POIList[iBin]:
                parabolas.append( Tgs[name][POI] )
            # Otherwise, have to find what partial bin matches
            else:
                for iBinPartial in xrange( len(POIList) ):
                    leftBound  = binBoundaries[name][iBinPartial]
                    rightBound = binBoundaries[name][iBinPartial+1]
                    POIPartial = POIList[iBinPartial]
                    if observableMin == leftBound or observableMax == rightBound:
                        parabolas.append( Tgs[name][POIPartial] )
                        break

        parabolasPerBin.append( parabolas )


    # ======================================
    # Create plots

    for iBin in xrange(nBins):

        c.Clear()

        base = ROOT.TH1F()
        base.Draw('P')
        base.GetXaxis().SetLimits( POImin, POImax )
        base.SetMinimum( deltaNLLmin )
        base.SetMaximum( deltaNLLmax )
        base.SetMarkerColor(0)
        base.GetXaxis().SetTitle( '#mu' )
        base.GetYaxis().SetTitle( '#Delta NLL' )

        leg = ROOT.TLegend( 0.65, 0.5, 1-RightMargin, 1-TopMargin )
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.Draw()

        for Tg in parabolasPerBin[iBin]:
            Tg.Draw('SAME')
            leg.AddEntry( Tg.GetName(), '{0}_{1}'.format( Tg.channel, Tg.POI ), 'l' )

        SaveC( 'cparabolas_{0}'.format( finestPOIs[iBin] ) )








        
def rindex( someList, val ):
    # Regular list.index() finds first instance in list, this function finds the last
    return len(someList) - someList[::-1].index(val) - 1

def FindMinimaAndErrors( POIvals, deltaNLLs ):

    minDeltaNLL   = min(deltaNLLs)
    iMin          = rindex( deltaNLLs, minDeltaNLL )
    minimumPOIval = POIvals[iMin]

    # Dict that is returned in case of a total error
    errReturn = {
        'imin'       : iMin,
        'min'        : minimumPOIval,
        'leftError'  : 999,
        'leftBound'  : 999,
        'rightError' : 999,
        'rightBound' : 999,
        }

    if iMin > 2:
        # Find left minimum
        Tg_left = ROOT.TGraph(
            iMin+1,
            array( 'd', deltaNLLs[:iMin+1] ),
            array( 'd', POIvals[:iMin+1] )
            )
        ROOT.SetOwnership( Tg_left, False )

        leftBound = Tg_left.Eval( 0.5 )
        if leftBound <= POIvals[0]:
            wellDefinedLeftBound = False
        else:
            wellDefinedLeftBound = True
    else:
        wellDefinedLeftBound = False

    if iMin < len(POIvals)-2:
        # Find right minimum
        Tg_right = ROOT.TGraph(
            len(POIvals)-iMin+1,
            array( 'd', deltaNLLs[iMin:] ),
            array( 'd', POIvals[iMin:] )
            )
        ROOT.SetOwnership( Tg_right, False )

        rightBound = Tg_right.Eval( 0.5 )
        if rightBound >= POIvals[-1]:
            wellDefinedRightBound = False
        else:
            wellDefinedRightBound = True
    else:
        wellDefinedRightBound = False


    if wellDefinedLeftBound and wellDefinedRightBound:
        pass
    # Symmetrize if one of the bounds was poorly defined
    elif wellDefinedLeftBound and not wellDefinedRightBound:
        rightBound = minimumPOIval + ( minimumPOIval - leftBound )
    elif wellDefinedRightBound and not wellDefinedLeftBound:
        leftBound  = minimumPOIval - ( rightBound - minimumPOIval )
    else:
        print 'Hopeless interpolation case; unable to determine uncertainties.'
        return errReturn

    leftError = abs(minimumPOIval - leftBound)
    rightError = abs(minimumPOIval - rightBound)

    return {
        'imin'       : iMin,
        'min'        : minimumPOIval,
        'leftError'  : leftError,
        'leftBound'  : leftBound,
        'rightError' : rightError,
        'rightBound' : rightBound,
        }








########################################
# Execute this if called directly
########################################

if __name__ == '__main__':
    print 'Import functions from this file'
    print 'Available functions:'
    functions = [ name for (name, thing) in locals().items() if callable(thing) ]
    print '\n'.join( functions )