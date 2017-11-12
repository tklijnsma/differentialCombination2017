#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import Commands

import os, tempfile, shutil, re, glob, itertools, sys, numpy, operator, pprint, re
from os.path import *
from operator import itemgetter
from array import array
from math import log, exp, sqrt, copysign
from copy import deepcopy

from time import strftime
datestr = strftime( '%b%d' )

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
ROOT.gStyle.SetOptStat(0)
c = ROOT.TCanvas( 'ctc', 'ctc', 1000, 800 )

LeftMargin   = 0.15
RightMargin  = 0.03
BottomMargin = 0.15
TopMargin    = 0.03
def SetCMargins(
    LeftMargin   = 0.15,
    RightMargin  = 0.03,
    BottomMargin = 0.15,
    TopMargin    = 0.03,
    ):
    c.SetLeftMargin( LeftMargin )
    c.SetRightMargin( RightMargin )
    c.SetBottomMargin( BottomMargin )
    c.SetTopMargin( TopMargin )


PLOTDIR = 'plots_{0}'.format(datestr)
def SetPlotDir( newdir ):
    global PLOTDIR
    PLOTDIR = newdir
def SaveC( outname, asPNG=False, asROOT=False ):
    global PLOTDIR
    if not isdir(PLOTDIR): os.makedirs(PLOTDIR)

    subdir = ''
    if len(outname.split('/')) == 2:
        subdir = outname.split('/')[0]
        if not isdir( join( PLOTDIR, subdir ) ): os.makedirs( join( PLOTDIR, subdir ) )

    outname = join( PLOTDIR, subdir, basename(outname).replace('.pdf','').replace('.png','') )
    c.SaveAs( outname + '.pdf' )
    if asPNG:
        c.SaveAs( outname + '.png' )

    if asROOT:
        c.SaveAs( outname + '.root' )


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
    SetTitleSizes = True,
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

    if SetTitleSizes:
        base.GetXaxis().SetTitleSize( 0.06 )
        base.GetYaxis().SetTitleSize( 0.06 )

    return base



# class Container:
#     def __init__(self, **kwds):
#         print 'Warning: Should replace this class with the more common Container.Container'
#         self.__dict__.update(kwds)
from Container import Container


########################################
# Derived Theory Containers
########################################

#____________________________________________________________________
def RebinDerivedTheoryContainer( container, newBinBoundaries ):

    newNBins = len(newBinBoundaries) - 1

    newCrosssection = Rebin(
        theoryBinBoundaries = container.binBoundaries,
        theoryBinValues     = container.crosssection,
        expBinBoundaries    = newBinBoundaries,
        lastBinIsOverflow   = True,
        )

    newRatios = Rebin(
        theoryBinBoundaries = container.binBoundaries,
        theoryBinValues     = container.ratios,
        expBinBoundaries    = newBinBoundaries,
        lastBinIsOverflow   = True,
        )

    # newBinCenters = [ 0.5*(newBinBoundaries[i]+newBinBoundaries[i+1]) for i in xrange(newNBins) ]

    container.binBoundaries = newBinBoundaries
    container.crosssection  = newCrosssection
    container.ratios        = newRatios


#____________________________________________________________________
# Theory files contain a somewhat non-consistent binning, turn it into a well-defined binning
def BinningHeuristic(
    binCenters,
    manualSwitchAt50 = True,
    manualSwitchAt5  = False,
    ):

    nBins = len(binCenters)
    binBoundaries = []

    # Assume first bin center is in the middle; use second bin to determine first boundary
    binBoundaries.append( binCenters[0] - 0.5*(binCenters[1]-binCenters[0]) )

    # Bin boundaries are defined by being the middle of bin centers
    for iBin in xrange(0,nBins-1):

        # Manually overwrite for the cross-over at pt = 50 (which is irregular)
        if manualSwitchAt50 and binCenters[iBin] == 48.5:
            binBoundaries.append( 51 )

        # For the quark induced histograms there is an irregularity at pt = 5.0
        elif manualSwitchAt5 and ( binCenters[iBin] == 3.0 or binCenters[iBin] == 2.75 ):
            binBoundaries.append( 5.0 )

        else:
            binBoundaries.append(
            binCenters[iBin] + 0.5*(binCenters[iBin+1]-binCenters[iBin])
            )

    binBoundaries.sort()

    # Last point gets the same width as previous bin
    binBoundaries.append( binCenters[-1] + 0.5*(binCenters[-1]-binCenters[-2]) )

    # Bin centers may have changed because of this definition
    newBinCenters = [ 0.5*(binBoundaries[i+1]+binBoundaries[i]) for i in xrange(nBins-1) ]
    binWidths = [ binBoundaries[i+1]-binBoundaries[i] for i in xrange(nBins-1) ]

    return newBinCenters, binBoundaries, binWidths



#____________________________________________________________________
def GetTheoryTGraph(
    name,
    ptPoints,
    muPoints,
    muBoundLeft   = None,
    muBoundRight  = None,
    boundaries    = False,
    ):
    
    print 'WARNING: Better to use OutputInterface.OutputContainer.GetTGraph()'

    if not boundaries:
        # Default setting; supplied pt points are bin centers

        if not len(ptPoints) == len(muPoints):
            Commands.ThrowError(
                'Length of input lists are not set right\n'
                '    len(ptPoints) = {0} is not equal to len(muPoints) = {1}'.format(
                    len(ptPoints), len(muPoints) )
                )
            return

        # ======================================
        # Binning heuristic for pt points

        nBins = len(ptPoints)

        binCenters, binBoundaries, binWidths = BinningHeuristic( ptPoints )
        halfBinWidths = [ 0.5*w for w in binWidths ]

    else:
        # Supplied pt points are bin boundaries

        if not len(ptPoints)-1 == len(muPoints):
            Commands.ThrowError(
                'Length of input lists are not set right\n'
                '    len(ptPoints)-1 = {0} is not equal to len(muPoints) = {1}'.format(
                    len(ptPoints)-1, len(muPoints) )
                )
            return

        # ======================================
        # Standard calculation for centers and widths

        nBins = len(ptPoints)-1

        binBoundaries = ptPoints
        binCenters    = [ 0.5*(binBoundaries[i]+binBoundaries[i+1]) for i in xrange(nBins) ]
        binWidths     = [ (binBoundaries[i+1]-binBoundaries[i]) for i in xrange(nBins) ]
        halfBinWidths = [ 0.5*(binBoundaries[i+1]-binBoundaries[i]) for i in xrange(nBins) ]



    # ======================================
    # Make TGraph

    if muBoundRight == None or muBoundLeft == None:
        Tg = ROOT.TGraphAsymmErrors(
            len(binCenters),
            array( 'd', binCenters ),
            array( 'd', muPoints ),
            array( 'd', halfBinWidths ),
            array( 'd', halfBinWidths ),
            array( 'd', [ 0 for i in xrange(nBins) ] ),
            array( 'd', [ 0 for i in xrange(nBins) ] ),
            )
        # muBoundLeft  = [ 0 for i in xrange(nBins) ]
        # muBoundRight = [ 0 for i in xrange(nBins) ]

    else:
        Tg = ROOT.TGraphAsymmErrors(
            len(binCenters),
            array( 'd', binCenters ),
            array( 'd', muPoints ),
            array( 'd', halfBinWidths ),
            array( 'd', halfBinWidths ),
            array( 'd', [ c-abs(l) for l, c in zip( muBoundLeft, muPoints ) ] ),
            array( 'd', [ abs(r)-c for r, c in zip( muBoundRight, muPoints ) ] ),
            )

    ROOT.SetOwnership( Tg, False )
    Tg.SetLineWidth(2)
    Tg.SetName( name )


    # ======================================
    # Give some extra attributes

    Tg.name          = name
    Tg.binValues     = muPoints
    Tg.binCenters    = binCenters
    Tg.binWidths     = binWidths
    Tg.binBoundaries = binBoundaries
    # Tg.binErrLeft    = binErrLeft
    # Tg.binErrRight   = binErrRight

    Tg.xMin = min( binBoundaries )
    Tg.xMax = max( binBoundaries )
    Tg.fourth_xMin = sorted( binBoundaries )[4]
    Tg.fourth_xMax = sorted( binBoundaries, reverse=True )[4]

    if muBoundRight == None or muBoundLeft == None:
        Tg.yMin = min( muPoints )
        Tg.yMax = max( muPoints )
        Tg.fourth_yMin = sorted( muPoints )[4]
        Tg.fourth_yMax = sorted( muPoints, reverse=True )[4]
    else:
        Tg.yMin = min( muBoundLeft )
        Tg.yMax = max( muBoundRight )
        Tg.fourth_yMin = sorted( muBoundLeft )[4]
        Tg.fourth_yMax = sorted( muBoundRight, reverse=True )[4]

    # Try to set a coupling
    for coupling in [ 'ct', 'cb', 'cg' ]:
        match = re.search( r'{0}_([mp\d]+)'.format(coupling), name )
        if match:
            setattr( Tg, coupling, float(match.group(1).replace('p','.').replace('m','-')) )
        else:
            if coupling == 'cg':
                setattr( Tg, coupling, 0.0 )
            else:
                setattr( Tg, coupling, 1.0 )

    return Tg



#____________________________________________________________________
# Generic function the returns a function integral( a, b, verbose=False ), which is evaluated for binBoundaries and binValues
def GetIntegral( binBoundaries, binValues ):
    nBins = len(binValues)
    binCenters = [ 0.5*( binBoundaries[i+1] + binBoundaries[i] ) for i in xrange(nBins-1) ]

    # linearInterpolate = lambda a, x1, x2, y1, y2: \
    #     y1 + ( a - x1 ) / ( x2 - x1 ) * ( y2 - y1 )

    def integralfunction( a, b, verbose=False ):

        if verbose: print '\n\nInterpolation function called with a = {0} and b = {1}'.format( a, b )

        aInterpolated = False
        bInterpolated = False

        # Extrapolating
        if a < binBoundaries[0]:
            if verbose: print 'Warning: Have to extrapolate for integral (a={0}, min x={1})'.format( a, binBoundaries[0] )
            ya = 0.0
            ia = -1
            aInterpolated = True
        if b > binBoundaries[-1]:
            if verbose: print 'Warning: Have to extrapolate for integral (b={0}, max x={1})'.format( b, binBoundaries[-1] )
            # yb = linearInterpolate( b, binCenters[-2], binCenters[-1], binValues[-2], binValues[-1] )
            yb = 0.0
            ib = nBins
            bInterpolated = True

        # Find binValues at x = a and x = b
        for iBin in xrange(nBins):

            if aInterpolated and bInterpolated: break

            if not aInterpolated and ( a >= binBoundaries[iBin] and a < binBoundaries[iBin+1] ):
                ya = binValues[iBin]
                ia = iBin
                aInterpolated = True

            if not bInterpolated and ( b > binBoundaries[iBin] and b <= binBoundaries[iBin+1] ):
                yb = binValues[iBin]
                ib = iBin
                bInterpolated = True


        if not aInterpolated or not bInterpolated:
            Commands.ThrowError( 'Somehow this case was not interpolated; this is wrong' )
            return 0.

        binBoundaries_integration = [ a  ] + binBoundaries[ia+1:ib+1] + [ b  ]
        
        if ia == ib:
            binValues_integration = [ ya ]
        else:
            binValues_integration = [ ya ] + binValues[ia+1:ib] + [ yb ]



        if verbose:
            print 'Integrating over bins'
            print 'binBoundaries:'
            print '    ' + '\n    '.join( [ str(f) for f in  binBoundaries_integration ] )
            print 'binValues:'
            print '    ' + '\n    '.join( [ str(f) for f in  binValues_integration ] )

        integral = 0.
        if verbose: contributions = []
        for iBin in xrange(len(binValues_integration)):
            integral += binValues_integration[iBin] * ( binBoundaries_integration[iBin+1] - binBoundaries_integration[iBin] )

            if verbose:
                contribution = binValues_integration[iBin] * ( binBoundaries_integration[iBin+1] - binBoundaries_integration[iBin] )
                print 'Contribution of bin {0:4} = {1:8.3f}  ( binValue = {2:8.4f}, binBoundaries = {3:8.2f} to {4:8.2f}'.format(
                    iBin,
                    contribution,
                    binValues_integration[iBin],
                    binBoundaries_integration[iBin],
                    binBoundaries_integration[iBin+1],
                    )
                contributions.append(contribution)


        if verbose:
            print 'Sum of contributions: {0}'.format( sum(contributions) )

        return integral

    return integralfunction

#____________________________________________________________________
# Rebin from any binning scheme to another using the integral function
def Rebin(
        theoryBinBoundaries,
        theoryBinValues,
        expBinBoundaries,
        lastBinIsOverflow = False,
        verbose = False,
        ):
    theoryIntegralFunction = GetIntegral( theoryBinBoundaries, theoryBinValues )
    expBinValues = []
    for iBinExp in xrange(len(expBinBoundaries)-1):
        expBinValues.append(
            theoryIntegralFunction( expBinBoundaries[iBinExp], expBinBoundaries[iBinExp+1] ) / ( expBinBoundaries[iBinExp+1] - expBinBoundaries[iBinExp] )
            )
        if verbose: print 'Integral for {0:7.2f} to {1:7.2f}: {2}'.format( expBinValues[iBinExp], expBinValues[iBinExp+1], expBinValues[iBinExp] )

    if lastBinIsOverflow:
        expBinValues[-1] = theoryIntegralFunction( expBinBoundaries[-2], theoryBinBoundaries[-1] ) / ( theoryBinBoundaries[-1] - expBinBoundaries[-2] )

    return expBinValues


#____________________________________________________________________
# Not sure which funcion uses this
def GetShortTheoryName( names ):

    print 'Notification: TheoryCommands.GetShortTheoryName() was called'

    keyStrings = [ 'ct', 'cb', 'cg' ]

    uniqueCombinations = set()
    for name in names:
        presentKeyStrings = ''
        for keyString in keyStrings:
            if keyString in name: presentKeyStrings += keyString
        uniqueCombinations.add( presentKeyStrings )
    return '_'.join( list(uniqueCombinations) )



########################################
# Helper functions for plotting
########################################

#____________________________________________________________________
def BasicReadScan(
        rootfiles,
        xAttr,
        yAttr = 'deltaNLL',
        ):
    
    containers = Commands.ConvertTChainToArray(
        rootfiles,
        'limit',
        r'{0}|{1}'.format( xAttr, yAttr ),
        returnStyle = 'containerPerPoint',
        )

    # Filter out duplicate x values
    unique_xs = list(set([ getattr( c, xAttr ) for c in containers ]))

    filteredContainers = []
    for container in containers:
        x = getattr( container, xAttr )
        if x in unique_xs:
            filteredContainers.append( container )
            unique_xs.pop( unique_xs.index(x) )

        if len(unique_xs) == 0:
            break


    # Sort along x
    containers.sort( key = lambda c: getattr( c, xAttr ) )

    xs = [ getattr( c, xAttr ) for c in containers ]
    ys = [ getattr( c, yAttr ) for c in containers ]

    return xs, ys


#____________________________________________________________________
def WriteTH2DToFile(
        rootfiles,
        xCoupling = 'ct',
        yCoupling = 'cg',
        verbose = True,
        ):

    scan = Commands.ConvertTChainToArray(
        rootfiles
        )
    nPoints = len(scan[xCoupling])


    keys = [ yCoupling, xCoupling, 'deltaNLL' ]
    ret = []
    for iPoint in xrange(nPoints):
        ret.append( [ scan[key][iPoint] for key in keys ] )
    if verbose:
        pprint.pprint( [ keys ] + ret )


    def inferBinBoundaries( binCenters ):
        binBoundaries = []
        for iBin in xrange(len(binCenters)-1):
            binBoundaries.append( 0.5*(binCenters[iBin]+binCenters[iBin]) )
        binBoundaries = (
            [ binCenters[0] - (binBoundaries[0]-binCenters[0]) ] +
            binBoundaries +
            [ binCenters[-1] + (binCenters[-1]-binBoundaries[-1]) ]
            )
        return binBoundaries


    iBestfit = scan['deltaNLL'].index( 0.0 )
    xBestfit = scan[xCoupling][iBestfit]
    yBestfit = scan[yCoupling][iBestfit]

    xBinCenters = list(set(scan[xCoupling]))
    xBinCenters.pop( xBinCenters.index(xBestfit) )
    xBinCenters.sort()
    xNBins = len(xBinCenters)
    xBinBoundaries = inferBinBoundaries( xBinCenters )

    yBinCenters = list(set(scan[yCoupling]))
    yBinCenters.pop( yBinCenters.index(yBestfit) )
    yBinCenters.sort()
    yNBins = len(yBinCenters)
    yBinBoundaries = inferBinBoundaries( yBinCenters )


    H2 = ROOT.TH2D(
        'couplingScan', 'couplingScan',
        xNBins, array( 'd', xBinBoundaries ),
        yNBins, array( 'd', yBinBoundaries ),
        )

    for iPoint in xrange(nPoints):

        if scan[xCoupling][iPoint] == xBestfit and scan[yCoupling][iPoint] == yBestfit:
            continue

        iBinX = xBinCenters.index( scan[xCoupling][iPoint] )
        iBinY = yBinCenters.index( scan[yCoupling][iPoint] )

        nll = scan['deltaNLL'][iPoint]
        gauss = exp(-nll)

        H2.SetBinContent( iBinX+1, iBinY+1, gauss )


    H2outFile = 'scanTH2D_{0}.root'.format(datestr)
    H2outFp = ROOT.TFile.Open( H2outFile, 'RECREATE' )
    H2.Write( 'couplingScan' )
    H2outFp.Close()


#____________________________________________________________________
def GetContoursFromTH2( TH2_original, threshold, verbose=True ):

    # Open a temporary canvas so the contour business does not screw up other plots
    ctemp = ROOT.TCanvas( 'ctemp', 'ctemp', 1000, 800 )
    ctemp.cd()

    TH2 = TH2_original.Clone()
    TH2.SetName( GetUniqueRootName() )
    TH2.SetContour( 1, array( 'd', [threshold] ) )

    if verbose:
        print 'Trying to get contours from \'{0}\''.format( TH2.GetName() ) + ( ' ({0})'.format( TH2_original.name ) if hasattr( TH2_original, 'name' ) else '' )

    TH2.Draw( 'CONT Z LIST' )
    ROOT.gPad.Update()

    if verbose:
        print '    Contours found'

    contours_genObj = ROOT.gROOT.GetListOfSpecials().FindObject('contours')

    Tgs = []
    for iContour in xrange( contours_genObj.GetSize() ):
        contour_TList = contours_genObj.At(iContour)
        for iVal in xrange( contour_TList.GetSize() ):
            Tg = contour_TList.At(iVal)

            TgClone = Tg.Clone()
            ROOT.SetOwnership( TgClone, False )
            Tgs.append( TgClone )

    if verbose:
        for Tg in Tgs:
            N = Tg.GetN()
            xList = [ Tg.GetX()[iPoint] for iPoint in xrange(N) ]
            yList = [ Tg.GetY()[iPoint] for iPoint in xrange(N) ]
            xMin = min( xList )
            xMax = max( xList )
            yMin = min( yList )
            yMax = max( yList )
            print '    xMin = {0:+9.4f}, xMax = {1:+9.4f}, yMin = {2:+9.4f}, yMax = {3:+9.4f}'.format( xMin, xMax, yMin, yMax )


    del TH2
    del contours_genObj

    ctemp.Close()
    del ctemp
    c.cd()
    ROOT.gPad.Update()

    return Tgs


#____________________________________________________________________
def SetExtremaOfContour( Tg ):

    xBuffer = Tg.GetX()
    xs = [ xBuffer[i] for i in xrange( Tg.GetN() ) ]

    yBuffer = Tg.GetY()
    ys = [ yBuffer[i] for i in xrange( Tg.GetN() ) ]

    Tg.xMin = min(xs)
    Tg.xMax = max(xs)
    Tg.yMin = min(ys)
    Tg.yMax = max(ys)


#____________________________________________________________________
def GetTH2FromListOfRootFiles(
        rootfiles,
        xCoupling = 'ct',
        yCoupling = 'cg',
        verbose   = False,
        xMin = None, xMax = None, yMin = None, yMax = None,
        multiplyByTwo = True,
        zVariable = 'deltaNLL'
        ):

    # Read values from specified rootfiles
    scan = Commands.ConvertTChainToArray(
        rootfiles,
        returnStyle = 'dictPerPoint'
        )
    if xMin: scan = [ s for s in scan if s[xCoupling] >= xMin ]
    if xMax: scan = [ s for s in scan if s[xCoupling] <= xMax ]
    if yMin: scan = [ s for s in scan if s[yCoupling] >= yMin ]
    if yMax: scan = [ s for s in scan if s[yCoupling] <= yMax ]
    nPoints = len(scan)

    GetListFromScan = lambda key: [ s[key] for s in scan ]

    # keys = [ yCoupling, xCoupling, 'deltaNLL' ]
    # ret = Container( scanPoints=[] )
    # for iPoint in xrange(nPoints):
    #     ret.scanPoints.append( [ scan[key][iPoint] for key in keys ] )
    # if verbose:
    #     pprint.pprint( [ keys ] + ret.scanPoints )


    def inferBinBoundaries( binCenters ):
        binBoundaries = []
        for iBin in xrange(len(binCenters)-1):
            binBoundaries.append( 0.5*(binCenters[iBin]+binCenters[iBin]) )
        binBoundaries = (
            [ binCenters[0] - (binBoundaries[0]-binCenters[0]) ] +
            binBoundaries +
            [ binCenters[-1] + (binCenters[-1]-binBoundaries[-1]) ]
            )
        return binBoundaries


    iBestfit = GetListFromScan('deltaNLL').index( 0.0 )
    xBestfit = GetListFromScan(xCoupling)[iBestfit]
    yBestfit = GetListFromScan(yCoupling)[iBestfit]

    xBinCenters = list(set( GetListFromScan(xCoupling) ))
    xBinCenters.pop( xBinCenters.index(xBestfit) )
    xBinCenters.sort()
    xNBins = len(xBinCenters)
    xBinBoundaries = inferBinBoundaries( xBinCenters )

    yBinCenters = list(set( GetListFromScan(yCoupling) ))
    yBinCenters.pop( yBinCenters.index(yBestfit) )
    yBinCenters.sort()
    yNBins = len(yBinCenters)
    yBinBoundaries = inferBinBoundaries( yBinCenters )


    H2name = GetUniqueRootName()
    H2 = ROOT.TH2D(
        H2name, '',
        xNBins, array( 'd', xBinBoundaries ),
        yNBins, array( 'd', yBinBoundaries ),
        )
    ROOT.SetOwnership( H2, False )


    for iPoint in xrange(nPoints):

        if scan[iPoint][xCoupling] == xBestfit and scan[iPoint][yCoupling] == yBestfit:
            continue

        try:
            iBinX = xBinCenters.index( scan[iPoint][xCoupling] )
        except ValueError:
            print '[ERROR] Point {0} ({1}) not in list'.format( iPoint, scan[iPoint][xCoupling] )
            continue

        try:
            iBinY = yBinCenters.index( scan[iPoint][yCoupling] )
        except ValueError:
            print '[ERROR] Point {0} ({1}) not in list'.format( iPoint, scan[iPoint][yCoupling] )
            continue

        if multiplyByTwo:
            H2.SetBinContent( iBinX+1, iBinY+1, 2.*scan[iPoint][zVariable] )
        else:
            H2.SetBinContent( iBinX+1, iBinY+1, scan[iPoint][zVariable] )


    # Open return object
    ret = Container()

    ret.H2             = H2
    ret.zVariable      = zVariable
    ret.xCoupling      = xCoupling
    ret.yCoupling      = yCoupling
    ret.iBestfit       = iBestfit
    ret.xBestfit       = xBestfit
    ret.yBestfit       = yBestfit
    ret.xBinCenters    = xBinCenters
    ret.xNBins         = xNBins
    ret.xBinBoundaries = xBinBoundaries
    ret.yBinCenters    = yBinCenters
    ret.yNBins         = yNBins
    ret.yBinBoundaries = yBinBoundaries

    return ret
    


########################################
# Plotting
########################################

#____________________________________________________________________
def PlotCouplingScan2D(
        # datacard,
        rootfiles,
        xCoupling    = 'ct',
        yCoupling    = 'cg',
        SM           = None,
        verbose      = True,
        drawContours = True,
        xMin = None, xMax = None, yMin = None, yMax = None,
        multiplyBinContents = None,
        # 
        zVariable = 'deltaNLL',
        multiplyByTwo = True,
        zMin = None, zMax = None,
        ):

    res = GetTH2FromListOfRootFiles(
        rootfiles,
        xCoupling,
        yCoupling,
        verbose,
        xMin = xMin, xMax = xMax, yMin = yMin, yMax = yMax,
        zVariable = zVariable,
        multiplyByTwo = multiplyByTwo,
        )
    H2 = res.H2

    if not zVariable == 'deltaNLL':
        H2.SetTitle(zVariable)

    if verbose or multiplyBinContents:
        for iX in xrange( H2.GetNbinsX() ):
            for iY in xrange( H2.GetNbinsY() ):

                if verbose:
                    if abs(iX-(H2.GetNbinsX()/2)) < 5 and abs(iY-(H2.GetNbinsY()/2)) < 5:
                        print 'iX = {0:2d}, iY = {1:2d}, X = {2:+5.2f}, Y = {3:+5.2f}, {5} = {4}'.format(
                            iX, iY,
                            H2.GetXaxis().GetBinCenter(iX+1), H2.GetYaxis().GetBinCenter(iY+1),
                            H2.GetBinContent( iX+1, iY+1 ),
                            res.zVariable
                            )
                if multiplyBinContents:
                    H2.SetBinContent(
                        iX+1, iY+1,
                        H2.GetBinContent( iX+1, iY+1 ) * multiplyBinContents
                        )

    c.Clear()
    SetCMargins(
        LeftMargin   = 0.12,
        RightMargin  = 0.10,
        BottomMargin = 0.12,
        TopMargin    = 0.09,
        )


    H2.Draw('COLZ')


    if drawContours:

        contourTgs_1sigma = GetContoursFromTH2( H2, 2.30 )
        for Tg_1sigma in contourTgs_1sigma:
            Tg_1sigma.SetLineColor(1)
            Tg_1sigma.SetLineWidth(3)
            Tg_1sigma.Draw('LSAME')
            Tg_1sigma.SetName('contour1sigma')

        contourTgs_2sigma = GetContoursFromTH2( H2, 6.18 )
        for Tg_2sigma in contourTgs_2sigma:
            Tg_2sigma.SetLineColor(1)
            Tg_2sigma.SetLineWidth(3)
            Tg_2sigma.SetLineStyle(2)
            Tg_2sigma.Draw('LSAME')
            Tg_2sigma.SetName('contour2sigma')


    axisTitleDict = {
        'kappab' : '#kappa_{b}',
        'kappac' : '#kappa_{c}',
        'ct'     : '#kappa_{t}',
        'cg'     : '#kappa_{g}',
        }

    H2.GetXaxis().SetTitle( axisTitleDict.get( xCoupling, xCoupling ) )
    H2.GetYaxis().SetTitle( axisTitleDict.get( yCoupling, yCoupling ) )
    H2.GetXaxis().SetTitleSize(0.06)
    H2.GetYaxis().SetTitleSize(0.06)
    H2.GetXaxis().SetTitleOffset(0.9)
    H2.GetYaxis().SetTitleOffset(0.9)

    H2.GetXaxis().SetLabelSize(0.045)
    H2.GetYaxis().SetLabelSize(0.045)
    H2.GetZaxis().SetLabelSize(0.045)


    # H2.GetZaxis().SetLimits( 0., 50. )
    # H2.GetZaxis().SetRange( 0, 10 )

    if zMax:
        H2.SetMaximum(zMax)
    elif res.zVariable == 'deltaNLL':
        H2.SetMaximum( 7. )

    if zMin: H2.SetMinimum(zMin)

    c.Update()


    Tpoint = ROOT.TGraph( 1, array( 'd', [res.xBestfit] ), array( 'd', [res.yBestfit] ) )
    Tpoint.SetMarkerSize(2)
    Tpoint.SetMarkerStyle(34)
    Tpoint.Draw('P')
    Tpoint.SetName('BestfitPoint')

    if SM is None:
        xSM = 1.0
        ySM = 1.0
    else:
        xSM, ySM = SM

    Tpoint_SM = ROOT.TGraph( 1, array( 'd', [xSM] ), array( 'd', [ySM] ) )
    Tpoint_SM.SetMarkerSize(2)
    Tpoint_SM.SetMarkerStyle(21)
    Tpoint_SM.Draw('P')
    Tpoint_SM.SetName('SMPoint')

    if drawContours:
        leg = ROOT.TLegend(
            c.GetRightMargin() + 0.08,
            0.99,
            c.GetRightMargin() + 0.08 + 0.7,
            0.99 - c.GetTopMargin() + 0.02
            )
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetNColumns(4)
        leg.AddEntry( Tg_1sigma.GetName(), '1 #sigma', 'l' )
        leg.AddEntry( Tg_2sigma.GetName(), '2 #sigma', 'l' )
        leg.AddEntry( Tpoint.GetName(),    'Best fit', 'p' )
        leg.AddEntry( Tpoint_SM.GetName(), 'SM',       'p' )
        leg.Draw()

    outname = 'couplingscan2D_{0}'.format( basename(dirname(rootfiles[0])).replace('/','') )
    if not zVariable == 'deltaNLL': outname = outname.replace( 'couplingscan2D_', 'couplingscan2D_{0}_'.format(zVariable) )
    SaveC( outname )
    SetCMargins()

    return res
    

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
        ):


    # ======================================
    # Check whether the passed containers fulfill requirements

    for container in containers:
        attrs = container.ListAttributes()

        for expectedAttr in [ 'H2', 'name' ]:
            if not expectedAttr in attrs:
                Commands.ThrowError(
                    'Container misses mandatory attribute \'{0}\' (defined attributes: {1})'.format( expectedAttr, ', '.join(attrs) ),
                    throwException = True
                    )

        if not hasattr( container, 'color' ):
            container.color = 1


    # ======================================
    # Calculate contours

    for container in containers:
        print 'Getting contours for {0}'.format( container.name )
        container.contours_1sigma = GetContoursFromTH2( container.H2, 2.30 )
        container.contours_2sigma = GetContoursFromTH2( container.H2, 6.18 )



    # ======================================
    # Make plot

    c.cd()
    c.Clear()
    SetCMargins()


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
    leg.SetNColumns(3)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)


    for container in containers:

        for Tg in container.contours_1sigma:
            Tg.SetLineWidth(2)
            Tg.SetLineColor( container.color )
            Tg.SetLineStyle(1)
            Tg.Draw('CSAME')
            if Tg == container.contours_1sigma[0]:
                Tg.SetName( '{0}_contour_1sigma'.format(container.name) )
                leg.AddEntry(
                    Tg.GetName(),
                    container.name if not hasattr( container, 'title' ) else container.title,
                    'l' )

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

    TpointSM = ROOT.TGraph( 1, array( 'd', [x_SM] ), array( 'd', [y_SM] ) )
    ROOT.SetOwnership( TpointSM, False )
    TpointSM.SetMarkerSize(2)
    TpointSM.SetMarkerStyle(21)
    TpointSM.SetMarkerColor( 12 )
    TpointSM.Draw('PSAME')

    leg.Draw()

    SaveC( plotname )


    if plotIndividualH2s:

        for container in containers:

            c.Clear()
            SetCMargins(
                LeftMargin   = 0.12,
                RightMargin  = 0.10,
                BottomMargin = 0.12,
                TopMargin    = 0.09,
                )

            container.H2.Draw('COLZ')

            container.H2.GetXaxis().SetRangeUser( xMin, xMax )
            container.H2.GetYaxis().SetRangeUser( yMin, yMax )
            container.H2.SetMaximum( 7. )

            container.H2.GetXaxis().SetTitle( xTitle )
            container.H2.GetYaxis().SetTitle( yTitle )
            container.H2.GetXaxis().SetTitleSize(0.06)
            container.H2.GetXaxis().SetLabelSize(0.05)
            container.H2.GetYaxis().SetTitleSize(0.06)
            container.H2.GetYaxis().SetLabelSize(0.05)

            for Tg in container.contours_1sigma: Tg.Draw('CSAME')
            for Tg in container.contours_2sigma: Tg.Draw('CSAME')
            container.bestfitPoint.Draw('PSAME')

            c.Update()

            SaveC( plotname + '_' + container.name )


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
        ):

    if xMin is None: xMin = min([ min(container.x) for container in containers ])
    if xMax is None: xMax = max([ max(container.x) for container in containers ])
    if yMin is None: yMin = min([ min(container.y) for container in containers ])
    if yMax is None: yMax = max([ max(container.y) for container in containers ])


    # ======================================
    # Make plot

    c.cd()
    c.Clear()
    SetCMargins()


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
    leg.SetNColumns(3)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    colorCycle = itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )
    for container in containers:

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

            container.Tg.Draw('PCSAME')

    if draw1sigmaline:
        line1sigma = ROOT.TLine( xMin, 0.5, xMax, 0.5 )
        line1sigma.SetLineColor(14)
        line1sigma.Draw()

    SaveC( plotname )




########################################
# Execute this if called directly
########################################

if __name__ == '__main__':
    print 'Import functions from this file'
    print 'Available functions:'
    functions = [ name for (name, thing) in locals().items() if callable(thing) ]
    print '\n'.join( functions )