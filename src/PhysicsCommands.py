#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import Commands
import TheoryCommands
import CorrelationMatrices

import os, tempfile, shutil, re, glob, itertools, sys, numpy, operator
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

def SaveC( outname, PNG=False, asROOT=False ):
    global PLOTDIR
    if not isdir(PLOTDIR): os.makedirs(PLOTDIR)

    subdir = ''
    if len(outname.split('/')) == 2:
        subdir = outname.split('/')[0]
        if not isdir( join( PLOTDIR, subdir ) ): os.makedirs( join( PLOTDIR, subdir ) )

    outname = join( PLOTDIR, subdir, basename(outname).replace('.pdf','').replace('.png','') )
    c.SaveAs( outname + '.pdf' )
    if PNG:
        c.SaveAs( outname + '.png' )
    if asROOT or TheoryCommands.SAVEROOT:
        c.SaveAs( outname + '.root' )


########################################
# Common functions
########################################

ROOTCOUNTER = 1000
def GetUniqueRootName():
    global ROOTCOUNTER
    name = 'root{0}'.format(ROOTCOUNTER)
    ROOTCOUNTER += 1
    return name


def GetPlotBase(
    xMin = 0, xMax = 1,
    yMin = 0, yMax = 1,
    xTitle = 'x', yTitle = 'y',
    ):

    base = ROOT.TH1F()
    ROOT.SetOwnership( base, False )
    base.SetName( GetUniqueRootName() )
    base.GetXaxis().SetLimits( xMin, xMax )
    base.SetMinimum( yMin )
    base.SetMaximum( yMax )
    base.SetMarkerColor(0)
    base.GetXaxis().SetTitle( xTitle )
    base.GetYaxis().SetTitle( yTitle )

    return base


########################################
# Dealing with scans in combine
########################################

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
            scanFiles,
            'limit',
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
    base.SetMaximum( min( 6., deltaNLLRange[1] ) )
    base.SetMarkerColor(0)
    base.GetXaxis().SetTitle( '#mu' )
    base.GetYaxis().SetTitle( '2#DeltaNLL' )

    base.GetXaxis().SetTitleSize( 0.06 )
    base.GetYaxis().SetTitleSize( 0.06 )

    leg = ROOT.TLegend( 0.5, 0.5, 1-RightMargin, 1-TopMargin )
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)


    # Create predefined list of colors
    colors = range(2,9+1) + [ 30, 38, 40, 41, 42 ] + range( 45, 48+1 )
    colorCycle = itertools.cycle( colors )

    def GetColor( POI ):
        # knownPOIs = {
        #     'r_smH_PTH_0_15'    : ROOT.TColor.GetColor( colors[0] ),
        #     'r_smH_PTH_15_30'   : ROOT.TColor.GetColor( colors[1] ),
        #     'r_smH_PTH_30_45'   : ROOT.TColor.GetColor( colors[2] ) - 1,
        #     'r_smH_PTH_45_85'   : ROOT.TColor.GetColor( colors[2] ) + 1,
        #     'r_smH_PTH_85_125'  : ROOT.TColor.GetColor( colors[3] ) - 1,
        #     'r_smH_PTH_125_200' : ROOT.TColor.GetColor( colors[3] ) + 1,
        #     'r_smH_PTH_200_350' : ROOT.TColor.GetColor( colors[4] ),
        #     'r_smH_PTH_GT350'   : ROOT.TColor.GetColor( colors[5] ),
        #     'r_smH_PTH_30_85'   : ROOT.TColor.GetColor( colors[2] ),
        #     'r_smH_PTH_85_200'  : ROOT.TColor.GetColor( colors[3] ),
        #     'r_smH_PTH_GT200'   : ROOT.TColor.GetColor( colors[5] ),
        #     }

        knownPOIs = {
            'r_smH_PTH_0_15'    : colorCycle.next(),
            'r_smH_PTH_15_30'   : colorCycle.next(),
            'r_smH_PTH_30_45'   : colorCycle.next(),
            'r_smH_PTH_45_85'   : colorCycle.next(),
            'r_smH_PTH_85_125'  : colorCycle.next(),
            'r_smH_PTH_125_200' : colorCycle.next(),
            'r_smH_PTH_200_350' : colorCycle.next(),
            'r_smH_PTH_GT350'   : colorCycle.next(),
            'r_smH_PTH_30_85'   : colorCycle.next(),
            'r_smH_PTH_85_200'  : colorCycle.next(),
            # 'r_smH_PTH_GT200'   : colorCycle.next(),
            'dummy'             : [ colorCycle.next() for i in xrange(5) ],
            # 
            'r_smH_NJ_0'        : colorCycle.next(),
            'r_smH_NJ_1'        : colorCycle.next(),
            'r_smH_NJ_2'        : colorCycle.next(),
            'r_smH_NJ_3'        : colorCycle.next(),
            'r_smH_NJ_GE4'      : colorCycle.next(),
            # 
            'r_smH_PTH_GT200'   : 29,
            }

        if POI in knownPOIs:
            return knownPOIs[POI]
        elif '_GT' in POI:
            return knownPOIs['r_smH_PTH_GT350'] 
        else:
            return colorCycle.next()





    for POI, scan in zip( POIs, scans ):

        POIvals, deltaNLLs = FilterScan( scan )

        nPoints = len(POIvals)

        # print '\n\n'
        # for POIval, deltaNLL in zip(POIvals, deltaNLLs):
        #     print POIval, deltaNLL


        Tg = ROOT.TGraph(
            nPoints,
            array( 'd', POIvals ),
            array( 'd', [ 2.*val for val in deltaNLLs ] ),
            )
        ROOT.SetOwnership( Tg, False )
        Tg.SetName( 'Tg_' + POI )

        # color = colors[iColor]
        # iColor += 1
        color = GetColor(POI)

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


        # Dots at the error bounds

        rightDot = ROOT.TGraph(
            1, array( 'f', [ minimaAndUncertainties['rightBound'] ] ), array( 'f', [ 1.0 ] )
            )
        ROOT.SetOwnership( rightDot, False )
        if minimaAndUncertainties['wellDefinedRightBound']:
            rightDot.SetMarkerStyle( 8 )
        else:
            rightDot.SetMarkerStyle( 5 )
        rightDot.SetMarkerSize(0.8)
        rightDot.SetMarkerColor(color)
        rightDot.Draw('PSAME')

        leftDot = ROOT.TGraph(
            1, array( 'f', [ minimaAndUncertainties['leftBound'] ] ), array( 'f', [ 1.0 ] )
            )
        ROOT.SetOwnership( leftDot, False )
        if minimaAndUncertainties['wellDefinedLeftBound']:
            leftDot.SetMarkerStyle( 8 )
        else:
            leftDot.SetMarkerStyle( 5 )
        leftDot.SetMarkerSize(0.8)
        leftDot.SetMarkerColor(color)
        leftDot.Draw('PSAME')

        leg.AddEntry( Tg.GetName(), POI, 'l' )

    leg.Draw()

    SaveC( 'parabolas_{0}_{1}_{2}'.format( productionMode, observableName, name ) )



def FigureOutBinning( POIs ):
    POIs.sort( key=lambda POI: Commands.InterpretPOI( POI )[2][0] )

    # Very ugly hack to make sure that LT30 is the first element of the sorted list
    if '_LT' in POIs[-1] or '_LE' in POIs[-1]:
        POIs = [ POIs[-1] ] + POIs[:-1]

    Ranges = [ Commands.InterpretPOI( POI )[2] for POI in POIs ]
    
    rangeLengthList = list(set(map( len, Ranges[:-1] )))
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
        binBoundaries = [ element[0] for element in Ranges ]
        lastBinWidth = binBoundaries[-1] - binBoundaries[-2]
        binBoundaries.append( binBoundaries[-1] + lastBinWidth )

    else:
        print 'ERROR: Bin observable range has length {0}, which does not make sense'.format( rangeLength )
        return

    return binBoundaries



def GetTGraphForSpectrum(
    POIs,
    scans,
    name = 'unnamed',
    scalePOIs = None,
    ):

    # Sort by first number of observable range
    POIs.sort( key=lambda POI: Commands.InterpretPOI( POI )[2][0] )
    binBoundaries = FigureOutBinning( POIs )

    binCenters    = [ 0.5*(left+right) for left, right in zip( binBoundaries[:-1], binBoundaries[1:] ) ]
    binWidths     = [ right-left for left, right in zip( binBoundaries[:-1], binBoundaries[1:] ) ]
    halfBinWidths = [ 0.5*(right-left) for left, right in zip( binBoundaries[:-1], binBoundaries[1:] ) ]

    if scalePOIs and len(scalePOIs) != len(binCenters):
        Commands.ThrowError( 'Length of scalePOIs {0} is not equal to length of binCenters {1}'.format( scalePOIs, binCenters ) )
        sys.exit()        

    POICenters  = []
    POIErrsLeft  = []
    POIErrsRight = []
    for POI, scan in zip( POIs, scans ):

        POIvals, deltaNLLs = FilterScan( scan )
        nPoints = len(POIvals)
        minimaAndUncertainties = FindMinimaAndErrors( POIvals, deltaNLLs )

        POICenter   = minimaAndUncertainties['min']
        POIErrLeft  = minimaAndUncertainties['leftError']
        POIErrRight = minimaAndUncertainties['rightError']

        if scalePOIs:
            POICenter   *= scalePOIs[ POIs.index(POI) ]
            POIErrLeft  *= scalePOIs[ POIs.index(POI) ]
            POIErrRight *= scalePOIs[ POIs.index(POI) ]

        POICenters.append(   POICenter )
        POIErrsLeft.append(  POIErrLeft )
        POIErrsRight.append( POIErrRight )


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

    Tg.POIErrsSymm   = [ 0.5*(abs(i)+abs(j)) for i, j in zip( POIErrsLeft, POIErrsRight ) ]
    Tg.POIErrsSymmPerc = [ 100.*err/(val+0.0000000001) for err, val in zip( Tg.POIErrsSymm, Tg.POICenters ) ]

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

    isRatioPlot = True
    for kwarg in kwargs:
        if '_SMXS' in kwarg:
            isRatioPlot = False

    IsBottomRatioPlot = kwargs.get( 'bottomRatioPlot', False )

    Tgs = []
    for name, POIList, scanList, drawoptionList in zip( names, POIs, scans, drawoptions ):
        Tg = GetTGraphForSpectrum( POIList, scanList, name=name, scalePOIs=kwargs.get( '{0}_SMXS'.format(name), None ) )

        Tg.DrawAsBlocks = True
        for drawoption in drawoptionList:
            if drawoption[0] == 'AsBlocks':
                Tg.DrawAsBlocks = drawoption[1]
                continue

            try:
                getattr( Tg, drawoption[0] )( *drawoption[1:] )
                # print 'Successfully called \'Tg.{0}({1})\''.format( drawoption[0], drawoption[1:] )
            except:
                Commands.ThrowError( 'Problem calling \'Tg.{0}({1})\', skipping'.format( drawoption[0], drawoption[1:] ) )

        Tgs.append( Tg )


    if kwargs.get( 'printTable', False ):
        for Tg in Tgs:

            print 'Binning, mus, err_symm, err_down, err_up for {0}'.format(Tg.name)

            # line1 = '{0:<15s}'.format( Tg.name )
            line1 = []
            for iBin in xrange(len(Tg.binBoundaries)-2):
                line1.append( '{0:<9s}'.format( '{0}-{1}'.format( int(Tg.binBoundaries[iBin]), int(Tg.binBoundaries[iBin+1]) ) ) )
            line1.append( '>{0:<8d}'.format( int(Tg.binBoundaries[-2]) ) )
            print ' , '.join( line1 )

            # line3 = '{0:<15s}'.format( Tg.name )
            line3 = []
            for iBin in xrange(len(Tg.binBoundaries)-1):
                line3.append( '{0:<+9.2f}'.format( Tg.POICenters[iBin] ) )
            print ' , '.join( line3 )

            # line2 = '{0:<15s}'.format( Tg.name )
            line2 = []
            for iBin in xrange(len(Tg.binBoundaries)-1):
                line2.append( '{0:<+9.2f}'.format( Tg.POIErrsSymm[iBin] ) )
            print ' , '.join( line2 )

            # line2a = '{0:<15s}'.format( 'err down' )
            line2a = []
            for iBin in xrange(len(Tg.binBoundaries)-1):
                line2a.append( '{0:<+9.2f}'.format( Tg.POIErrsLeft[iBin] ) )
            print ' , '.join( line2a )

            # line2b = '{0:<15s}'.format( 'err up' )
            line2b = []
            for iBin in xrange(len(Tg.binBoundaries)-1):
                line2b.append( '{0:<+9.2f}'.format( Tg.POIErrsRight[iBin] ) )
            print ' , '.join( line2b )





    # Use first POI to get production mode and observable name
    productionMode, observableName, dummyRange = Commands.InterpretPOI( POIs[0][0] )


    xMin = min([ Tg.xMin for Tg in Tgs ])
    xMax = max([ Tg.xMax for Tg in Tgs ])
    dx = xMax-xMin
    yMin = min([ Tg.yMin for Tg in Tgs ])
    yMax = max([ Tg.yMax for Tg in Tgs ])
    dy = yMax-yMin
    
    yMin -= 0.05*dy
    yMax += 0.05*dy

    titledict = {
        'combination' : 'Combination',
        'hgg'         : 'H #rightarrow #gamma#gamma',
        'hzz'         : 'H #rightarrow ZZ',
        # 'PTH'         : 'p_{T}^{H} [GeV]',
        'PTH'         : 'p_{T}^{H}',
        'NJ'          : 'N_{jets}'
        }

    c.Clear()
    c.SetLeftMargin( 0.18 )
    if not isRatioPlot:
        c.SetLogy()
        yMax = yMax * 2.
        yMinExternal = kwargs.get( 'yMin', 0.01 )
        yMin = max( yMin, yMinExternal )
        # xMax = 600.

    if IsBottomRatioPlot:
        c.SetTopMargin( 0.65 )
        c.SetBottomMargin( 0.01 )

    base = ROOT.TH1F()
    base.Draw('P')
    base.GetXaxis().SetLimits( xMin, xMax )
    base.SetMinimum( yMin )
    base.SetMaximum( yMax )
    base.SetMarkerColor(0)
    base.GetXaxis().SetTitle( titledict.get( observableName, observableName ) )
    
    if isRatioPlot:
        base.GetYaxis().SetTitle( '#mu' )
        base.GetYaxis().SetTitleSize( 0.06 )
    else:
        # 'd#sigma/dp_{T} [pb/GeV]'
        base.GetYaxis().SetTitle( '#frac{{#Delta#sigma({0})}}{{#Delta {0}}} [pb/GeV]'.format( titledict.get( observableName, observableName ) ) )
        base.GetYaxis().SetTitleSize( 0.05 )

    if IsBottomRatioPlot:
        base.GetXaxis().SetTitle('')
        base.GetXaxis().SetLabelOffset(999.)
        base.GetYaxis().SetNdivisions(505)

    base.GetXaxis().SetTitleSize( 0.06 )
    base.GetXaxis().SetLabelSize( 0.05 )
    base.GetYaxis().SetLabelSize( 0.05 )
    base.GetYaxis().SetTitleOffset( 1.5 )

    if isRatioPlot:
        lineAtOne = ROOT.TLine( xMin, 1.0, xMax, 1.0 )
        lineAtOne.SetLineWidth(3)
        lineAtOne.SetLineColor(17)
        # lineAtOne.SetLineStyle(2)
        lineAtOne.Draw()

    if kwargs.get( 'legendLeft', False ):
        leg = ROOT.TLegend( c.GetLeftMargin()+0.02, 1-TopMargin-0.3, c.GetLeftMargin()+0.25, 1-TopMargin )
    else:
        leg = ROOT.TLegend( 1-RightMargin-0.23, 1-TopMargin-0.3, 1-RightMargin, 1-TopMargin )
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)

    TgDrawStr = 'SAME P'

    for Tg in Tgs:

        Tg.title = titledict.get( Tg.name, Tg.name )

        # Make overflow bin extend until max of spectrum
        if Tg.binBoundaries[-1] < xMax:
            Tg.binBoundaries[-1] = xMax
            newBinCenter    = 0.5*( xMax + Tg.binBoundaries[-2] )
            newHalfBinWidth = 0.5*( xMax - Tg.binBoundaries[-2] )
            lastPoint = len(Tg.binCenters)-1
            Tg.SetPoint( lastPoint, newBinCenter, Tg.POICenters[lastPoint] )
            Tg.SetPointEXlow( lastPoint, newHalfBinWidth )
            Tg.SetPointEXhigh( lastPoint, newHalfBinWidth )


        if Tg.DrawAsBlocks:

            CorrelationMatrices.ConvertTGraphToLinesAndBoxes(
                Tg,
                drawImmediately=True,
                legendObject=leg,
                verbose=False,
                noBoxes=False,
                xMaxExternal=xMax,
                yMinExternal=yMin,
                yMaxExternal=yMax,
                )

            # # This is already transparent!
            # # Regular color 4 becomes something like 972 when transparent
            # fillColor = Tg.GetFillColor()
            # fillStyle = Tg.GetFillStyle()
            # lineColor = Tg.GetLineColor()

            # for iBin in xrange(len(Tg.binBoundaries)-1):

            #     if Tg.binBoundaries[iBin] > xMax:
            #         continue

            #     revertChange = False
            #     if Tg.binBoundaries[iBin+1] > xMax:
            #         # Change this back at end of loop
            #         actualBound = Tg.binBoundaries[iBin+1]
            #         Tg.binBoundaries[iBin+1] = xMax
            #         revertChange = True

            #     box = ROOT.TBox(
            #         Tg.binBoundaries[iBin],
            #         Tg.POICenters[iBin]-Tg.POIErrsLeft[iBin] if Tg.POICenters[iBin]-Tg.POIErrsLeft[iBin] > yMin else yMin,
            #         Tg.binBoundaries[iBin+1]                 if Tg.binBoundaries[iBin+1] < xMax else xMax,
            #         Tg.POICenters[iBin]+Tg.POIErrsRight[iBin]
            #         )
            #     ROOT.SetOwnership( box, False )

            #     box.SetLineWidth(0)
            #     box.SetFillColor( fillColor )
            #     if not fillStyle == 0: box.SetFillStyle( fillStyle )
            #     box.Draw()

            #     line = ROOT.TLine(
            #         Tg.binBoundaries[iBin], Tg.POICenters[iBin],
            #         Tg.binBoundaries[iBin+1], Tg.POICenters[iBin],
            #         )
            #     ROOT.SetOwnership( line, False )

            #     line.SetLineColor( lineColor )
            #     line.SetLineWidth( Tg.GetLineWidth() )
            #     line.Draw()

            #     if revertChange:
            #         Tg.binBoundaries[iBin+1] = actualBound
            #         revertChange = False


            # # Dummy for the legend
            # TgDummy = ROOT.TGraph( 1, array( 'd', [-999.]), array( 'd', [-999.]) )
            # ROOT.SetOwnership( TgDummy, False )
            # TgDummy.SetLineColor( lineColor )
            # TgDummy.SetLineWidth( Tg.GetLineWidth() )
            # TgDummy.SetFillColor( fillColor )
            # if not fillStyle == 0: TgDummy.SetFillStyle( fillStyle )
            # TgDummy.SetMarkerStyle(8)
            # TgDummy.SetMarkerColor( lineColor )
            # TgDummy.SetMarkerSize(0)
            # TgDummy.SetName( Tg.name + '_dummy' )
            # TgDummy.Draw('SAME')

            # leg.AddEntry( TgDummy.GetName(), Tg.title, 'lf' )


        else:
            Tg.Draw( TgDrawStr )
            leg.AddEntry( Tg.GetName(), Tg.title, 'lp' )




    if 'theoryCurves' in kwargs:
        theoryTgs = kwargs['theoryCurves']
        colorCycle = itertools.cycle( [ 30, 38, 40, 41, 42 ] + range( 45, 48+1 ) )
        for Tg in theoryTgs:

            if kwargs.get( 'autocolor', True ):
                color = next(colorCycle)
                Tg.SetLineColor(color)

            if 'drawTheoryCurvesAsLines' in kwargs and kwargs['drawTheoryCurvesAsLines']:
                lines, legendDummy = ConvertToLinesAndBoxes( Tg )
                for line in lines:
                    line.Draw()
                legendDummy.Draw('SAME')
                leg.AddEntry( legendDummy.GetName(), Tg.name, 'l' )
            else:
                Tg.Draw('L X')
                leg.AddEntry( Tg.GetName(), Tg.name, 'l' )

    if not IsBottomRatioPlot:
        leg.Draw()


    outname = 'mspectrum_{0}_{1}_{2}'.format( productionMode, observableName, '_'.join(names) )
    if 'theoryCurves' in kwargs: outname += '_withTheoryCurves_{0}'.format( TheoryCommands.GetShortTheoryName([Tg.name for Tg in theoryTgs]) )
    outname += kwargs.get( 'filenameSuffix', '' )
    if IsBottomRatioPlot: outname += 'bottomRatioPlot'
    SaveC( outname, asROOT = kwargs.get( 'asROOT', False ) )

    c.SetLogy(False)
    c.Clear()


def ConvertToLinesAndBoxes(
    Tg,
    noYerrors = True,
    ):
    
    fillColor = Tg.GetFillColor()
    fillStyle = Tg.GetFillStyle()
    lineColor = Tg.GetLineColor()


    if not hasattr( Tg, 'binValues' ):
        print 'ERROR: No attribute binValues found, and automatic computation not yet implemented'
        return


    lines = []
    for iBin in xrange(len(Tg.binBoundaries)-1):

        # Box not necessary as long as no Yerrors are done
        # box = ROOT.TBox(
        #     Tg.binBoundaries[iBin], Tg.POICenters[iBin]-Tg.POIErrsLeft[iBin],
        #     Tg.binBoundaries[iBin+1], Tg.POICenters[iBin]+Tg.POIErrsRight[iBin]
        #     )
        # ROOT.SetOwnership( box, False )

        line = ROOT.TLine(
            Tg.binBoundaries[iBin], Tg.binValues[iBin],
            Tg.binBoundaries[iBin+1], Tg.binValues[iBin],
            )
        ROOT.SetOwnership( line, False )

        line.SetLineColor( lineColor )
        line.SetLineWidth( Tg.GetLineWidth() )

        lines.append(line)

    # Create a dummy for the legend
    legendDummy = ROOT.TGraph( 1, array( 'd', [-999.]), array( 'd', [-999.]) )
    ROOT.SetOwnership( legendDummy, False )
    legendDummy.SetLineColor( lineColor )
    legendDummy.SetLineWidth( Tg.GetLineWidth() )
    legendDummy.SetFillColor( fillColor )
    if not fillStyle == 0: legendDummy.SetFillStyle( fillStyle )
    legendDummy.SetMarkerStyle(8)
    legendDummy.SetMarkerColor( lineColor )
    legendDummy.SetMarkerSize(0)
    legendDummy.SetName( Tg.name + '_dummy' )

    return lines, legendDummy







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

        SaveC( 'compareScans/cparabolas_{0}'.format( finestPOIs[iBin] ) )








        
def rindex( someList, val ):
    # Regular list.index() finds first instance in list, this function finds the last
    return len(someList) - someList[::-1].index(val) - 1

def FindMinimaAndErrors( POIvals, deltaNLLs, returnContainer=False ):

    minDeltaNLL   = min(deltaNLLs)
    iMin          = rindex( deltaNLLs, minDeltaNLL )
    minimumPOIval = POIvals[iMin]

    # Dict that is returned in case of a total error
    errReturn = {
        'imin'       : iMin,
        'min'        : minimumPOIval,
        'leftError'  : -999,
        'leftBound'  : -999,
        'rightError' : 999,
        'rightBound' : 999,
        'wellDefinedRightBound' : False,
        'wellDefinedLeftBound'  : False,
        }
    if returnContainer:
        errReturnContainer = TheoryCommands.Container()
        for key, value in errReturn.iteritems():
            setattr( errReturnContainer, key, value )
        errReturn = errReturnContainer

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


    returnDict = {
        'imin'       : iMin,
        'min'        : minimumPOIval,
        'leftError'  : leftError,
        'leftBound'  : leftBound,
        'rightError' : rightError,
        'rightBound' : rightBound,
        'wellDefinedRightBound' : wellDefinedRightBound,
        'wellDefinedLeftBound'  : wellDefinedLeftBound,
        }

    if returnContainer:
        container = TheoryCommands.Container()
        for key, value in returnDict.iteritems():
            setattr( container, key, value )
        return container
    else:
        return returnDict










########################################
# Execute this if called directly
########################################

if __name__ == '__main__':
    print 'Import functions from this file'
    print 'Available functions:'
    functions = [ name for (name, thing) in locals().items() if callable(thing) ]
    print '\n'.join( functions )