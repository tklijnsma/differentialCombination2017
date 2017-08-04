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


########################################
# Main
########################################

class Container:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


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

    if asROOT:
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



########################################
# Dealing with derived theory files
########################################

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


def ReadDerivedTheoryFile(
    derivedTheoryFile,
    returnContainer = True,
    verbose = False,
    ):

    with open( derivedTheoryFile, 'r' ) as theoryFp:
        lines = [ l.strip() for l in theoryFp.readlines() if len(l.strip()) > 0 and not l.strip().startswith('#') ]


    if returnContainer:
        container = Container()

        container.file = derivedTheoryFile

        if verbose:
            print '\nCreating container for \'{0}\''.format( derivedTheoryFile )

        for line in lines:
            key, value = line.split('=',1)

            try:
                value = [ float(v) for v in value.split(',') ]
                if len(value) == 1:
                    value = value[0]
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    pass

            if verbose:
                print 'Attribute \'{0}\' is set to:'.format( key )
                print value

            setattr( container, key, value )


        return container


    else:

        ret  = {}
        couplings = {}
        for line in lines:
            key, value = line.split('=',1)
            if key == 'binBoundaries':
                binBoundaries = [ float(v) for v in value.split(',') ]
                ret['binBoundaries'] = binBoundaries
            elif key == 'crosssection':
                crosssection = [ float(v) for v in value.split(',') ]
                ret['crosssection'] = crosssection
            elif key == 'ratios':
                ratios = [ float(v) for v in value.split(',') ]
                ret['ratios'] = ratios
            else:
                couplings[key] = float(value)

        ret['couplings'] = couplings

        # return binBoundaries, crosssection
        return ret



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



def GetTheoryTGraph(
    name,
    ptPoints,
    muPoints,
    muBoundLeft   = None,
    muBoundRight  = None,
    boundaries    = False,
    ):


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


def MapFineToCoarse(
        theoryBinBoundaries,
        theoryBinValues,
        expBinBoundaries,
        lastBinIsOverflow = False,
        verbose = False,
        ):

    print '[warning] \'MapFineToCoarse\' is now just an alias of \'Rebin\', and should no longer be used' 
    
    Rebin(
        theoryBinBoundaries,
        theoryBinValues,
        expBinBoundaries,
        lastBinIsOverflow = False,
        verbose = False,
        )

# def Rebin(
#     ptFine, sigmasFine,
#     ptCoarse,
#     verbose=False,
#     ):

#     integralfunction = GetIntegral( ptFine, sigmasFine )
    
#     sigmasCoarse = []
#     for iBinCoarse in xrange( len(ptCoarse)-1 ):
#         integral       = integralfunction( ptCoarse[iBinCoarse], ptCoarse[iBinCoarse+1], verbose=verbose )
#         integralPerGeV = integral / ( ptCoarse[iBinCoarse+1] - ptCoarse[iBinCoarse] )
#         if verbose: print 'Integral for {0:7.2f} to {1:7.2f}: {2}'.format( ptCoarse[iBinCoarse], ptCoarse[iBinCoarse+1], integralPerGeV )
#         sigmasCoarse.append( integralPerGeV )

#     return sigmasCoarse




def BasicTheoryPlot( Tgs, drawErrors=True ):

    c.Clear()

    xMin = min( [ Tg.xMin for Tg in Tgs ] )
    xMax = max( [ Tg.xMax for Tg in Tgs ] )

    # Actually take 1.1 * 4th maximum to kill some spikes
    # yMin = min( [ Tg.yMin for Tg in Tgs ] )
    # yMax = max( [ Tg.yMax for Tg in Tgs ] )
    yMin = min( [ Tg.fourth_yMin for Tg in Tgs ] )
    yMax = max( [ Tg.fourth_yMax for Tg in Tgs ] )


    base = GetPlotBase(
        xMin = xMin, xMax = xMax, yMin = yMin, yMax = yMax,
        xTitle = 'PTH', yTitle = '#sigma'
        )
    base.Draw('P')

    leg = ROOT.TLegend( 1-RightMargin-0.3, 1-TopMargin-0.3, 1-RightMargin, 1-TopMargin )
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    colorCycle = itertools.cycle( range(2,9+1) + [ 30, 38, 40, 41, 42 ] + range( 45, 48+1 ) )
    for Tg in Tgs:
        color = next(colorCycle)

        if drawErrors:
            Tg.Draw('L')
        else:
            Tg.Draw('L X')
        
        Tg.SetLineColor(color)
        leg.AddEntry( Tg.GetName(), Tg.name, 'l' )

    leg.Draw()

    SaveC( 'theory_{0}'.format( GetShortTheoryName([ Tg.name for Tg in Tgs ]) ) )



def GetShortTheoryName( names ):
    keyStrings = [ 'ct', 'cb', 'cg' ]

    uniqueCombinations = set()
    for name in names:
        presentKeyStrings = ''
        for keyString in keyStrings:
            if keyString in name: presentKeyStrings += keyString
        uniqueCombinations.add( presentKeyStrings )
    return '_'.join( list(uniqueCombinations) )





def GetParametrization(
    points,
    # = [
    #    ( 0.1**2 , 0.1*0.075 , 0.075 ),
    #    ( 0.5**2 , 0.5*0.042 , 0.042 ),
    #    ( 1.5**2 , 1.5*-0.042 , -0.042 ),
    #    ( 2.0**2 , 2.0*-0.083 , -0.083 ),
    #    ],
    yValues,
    testMode = False,
    ):

    # function = 'y = A*x1 + B*x2 + C*x3'

    nPoints = len(points)
    nPars = len(points[0])

    if nPoints > nPars:
        print 'More points than parameters given, system is overconstrained; Taking only the first {0} points'.format( nPars )
        nPoints = nPars
        points = points[:nPars]
        yValues = yValues[:nPars]
    elif nPars < nPoints:
        print 'Less points than parameters given, system is underconstrained'
        return


    if testMode:
        print 'Used points:'
        for point in points:
            print 'y = fn( {0} )'.format( ', '.join([ str(p) for p in point]) )


    xMatrix = numpy.array(points)
    xInv = numpy.linalg.inv(xMatrix)


    pointFunctions = []
    parameterValuesPerPoint = []
    for i in xrange(len(yValues[0])):

        yValsPerPoint = numpy.array( [ [ys[i]] for ys in yValues ])
        parameterValues = list(itertools.chain.from_iterable( xInv.dot( yValsPerPoint ) ))
        parameterValuesPerPoint.append( parameterValues )

        # NOTE THE i=i! otherwise i always points to the last element of the loop
        pointFunction = lambda point, i=i: sum([ parValue * arg for parValue, arg in zip( point, parameterValuesPerPoint[i] ) ])
        pointFunctions.append( pointFunction )


        if testMode:

            print '\n' + '-'*70 + '\nPoint {0}'.format(i)
            print '\nxInv:'
            print xInv
            print '\nyValsPerPoint:'
            print yValsPerPoint
            print '\nparameterValues:'
            print parameterValues

            print '\nPoint function tests:'
            print 'Real value = {1},  pointFunction = {0}'.format( pointFunction( points[0] ), yValsPerPoint[0] )
            copy = list(points[0][:])

            copy[0] = 1.1*points[0][0]
            print 'small up variation of pointFunction   = {0}'.format( pointFunction( copy ) )
            copy[0] = 0.9*points[0][0]
            print 'small down variation of pointFunction = {0}'.format( pointFunction( copy ) )


    functionForList = lambda *args: [ function(args) for function in pointFunctions ]
    return functionForList




def MapPredictionToExperimental(
    ptTheory, sigmaTheory,
    binning,
    verbose = False,
    makeTGraph = None,
    ):

    Commands.ThrowError( '\'MapPredictionToExperimental\' should really not be used anymore; use \'Rebin\' instead' )
    sys.exit()


    sigmas        = sigmaTheory
    binBoundaries = ptTheory
    nBins         = len(binBoundaries)-1
    binCenters    = [ 0.5*( binBoundaries[i] + binBoundaries[i+1] ) for i in xrange(nBins) ]
    binWidths     = [ ( binBoundaries[i+1] - binBoundaries[i] ) for i in xrange(nBins) ]
    halfBinWidths = [ 0.5*( binBoundaries[i+1] - binBoundaries[i] ) for i in xrange(nBins) ]

    nBinsExp      = len(binning)-1
    binCentersExp    = [ 0.5*( binning[i] + binning[i+1] ) for i in xrange(nBinsExp) ]
    binWidthsExp     = [ ( binning[i+1] - binning[i] ) for i in xrange(nBinsExp) ]
    halfBinWidthsExp = [ 0.5*( binning[i+1] - binning[i] ) for i in xrange(nBinsExp) ]


    if verbose:
        print 'Theory curve:'
        print '{0:9}  |  {1:9}  |  {2:9}  |  {3:9}'.format( 'pt left', 'pt right', 'pt center', 'sigma' )
        for iBin in xrange(nBins):
            print '{0:+9.2f}  |  {1:+9.2f}  |  {2:+9.2f}  |  {3:+9.5f}'.format(
                binBoundaries[iBin], binBoundaries[iBin+1], binCenters[iBin], sigmas[iBin]
                )

    if verbose: print '\nInterpolating and rebinning'
    # integralfunction = GetIntegral( binCenters, sigmas )
    integralfunction = GetIntegral( binBoundaries, sigmas )

    sigmaExpBinning = []
    for iBinExp in xrange(nBinsExp):
        leftBound      = binning[iBinExp]
        rightBound     = binning[iBinExp+1]
        integral       = integralfunction( leftBound, rightBound, verbose=verbose )
        integralPerGeV = integral / ( rightBound - leftBound )
        if verbose: print 'Integral for {0:7.2f} to {1:7.2f}: {2}'.format( leftBound, rightBound, integralPerGeV )
        sigmaExpBinning.append( integralPerGeV )


    if not makeTGraph == None:

        Tg = ROOT.TGraphAsymmErrors(
            nBinsExp,
            array( 'd', binCentersExp ),
            array( 'd', sigmaExpBinning, ),
            array( 'd', halfBinWidthsExp ),
            array( 'd', halfBinWidthsExp ),
            array( 'd', [ 0 for i in xrange(nBinsExp) ] ),
            array( 'd', [ 0 for i in xrange(nBinsExp) ] ),
            )
        ROOT.SetOwnership( Tg, False )
        Tg.SetName( makeTGraph )
        Tg.name = makeTGraph
        return Tg

    else:
        return sigmaExpBinning



########################################
# Making plots with output workspaces
########################################

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



def PlotCouplingScan2D(
        datacard,
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
    ret = Container( scanPoints=[] )
    for iPoint in xrange(nPoints):
        ret.scanPoints.append( [ scan[key][iPoint] for key in keys ] )
    if verbose:
        pprint.pprint( [ keys ] + ret.scanPoints )


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



    c.Clear()
    SetCMargins(
        LeftMargin   = 0.10,
        RightMargin  = 0.10,
        BottomMargin = 0.10,
        TopMargin    = 0.10,
        )

    # n_stops = 3
    # stops  = [ 0.0, 0.5, 1.0 ]
    # reds   = [ 0.0, 1.0, 1.0 ]
    # blues  = [ 1.0, 1.0, 0.0 ]
    # greens = [ 0.0, 1.0, 0.0 ]

    # n_stops = 2
    # stops  = [ 0.0, 1.0 ]
    # reds   = [ 1.0, 1.0 ]
    # blues  = [ 0.0, 1.0 ]
    # greens = [ 0.0, 1.0 ]


    # ROOT.TColor.CreateGradientColorTable(
    #     n_stops,
    #     array('d', stops ),
    #     array('d', reds ),
    #     array('d', greens ),
    #     array('d', blues ),
    #     255 )

    H2 = ROOT.TH2D(
        'couplingScan', '',
        xNBins, array( 'd', xBinBoundaries ),
        yNBins, array( 'd', yBinBoundaries ),
        )

    for iPoint in xrange(nPoints):

        if scan[xCoupling][iPoint] == xBestfit and scan[yCoupling][iPoint] == yBestfit:
            continue

        try:
            iBinX = xBinCenters.index( scan[xCoupling][iPoint] )
        except ValueError:
            print '[ERROR] Point {0} ({1}) not in list'.format( iPoint, scan[xCoupling][iPoint] )
            continue

        try:
            iBinY = yBinCenters.index( scan[yCoupling][iPoint] )
        except ValueError:
            print '[ERROR] Point {0} ({1}) not in list'.format( iPoint, scan[yCoupling][iPoint] )
            continue

        H2.SetBinContent( iBinX+1, iBinY+1, scan['deltaNLL'][iPoint] )


    H2.Draw('COLZ')

    H2.GetXaxis().SetTitle(xCoupling)
    H2.GetYaxis().SetTitle(yCoupling)
    H2.GetXaxis().SetTitleSize(0.05)
    H2.GetYaxis().SetTitleSize(0.05)
    H2.GetXaxis().SetTitleOffset(0.9)
    H2.GetYaxis().SetTitleOffset(0.9)

    # H2.GetZaxis().SetLimits( 0., 50. )
    # H2.GetZaxis().SetRange( 0, 10 )
    H2.SetMaximum( 7. )
    c.Update()


    Tpoint = ROOT.TGraph( 1, array( 'd', [xBestfit] ), array( 'd', [yBestfit] ) )
    Tpoint.SetMarkerSize(2)
    Tpoint.SetMarkerStyle(34)
    Tpoint.Draw('P')

    Tpoint_SM = ROOT.TGraph( 1, array( 'd', [1.] ), array( 'd', [1.] ) )
    Tpoint_SM.SetMarkerSize(2)
    Tpoint_SM.SetMarkerStyle(21)
    Tpoint_SM.Draw('P')


    SaveC( 'couplingscan2D', asROOT=True )
    SetCMargins()


    ret.xCoupling = xCoupling
    ret.yCoupling = yCoupling
    ret.xBestfit  = xBestfit
    ret.yBestfit  = yBestfit

    return ret
    


def TestParametrizationsInWorkspace(
    datacard,
    testcouplings = [
        { 'ct' : 1.75, 'cg' : -0.0625 },
        ],
    ):
    
    datacardFp = ROOT.TFile.Open( datacard )
    w = datacardFp.Get('w')

    couplings = []
    couplingsList = ROOT.RooArgList( w.set('POI') )
    for i in xrange( couplingsList.getSize() ):
        couplings.append( couplingsList[i] )

    yieldParameters = []
    yieldParameterList = ROOT.RooArgList( w.set('yieldParameters') )
    for i in xrange( yieldParameterList.getSize() ):
        yieldParameters.append( yieldParameterList[i] )

    parametrizations = []
    parametrizationList = ROOT.RooArgList( w.set('parametrizations') )
    for i in xrange( parametrizationList.getSize() ):
        parametrizations.append( parametrizationList[i] )

    datacardFp.Close()

    yieldParameters.sort( key = lambda i: Commands.InterpretPOI( i.GetName() )[2][0] )


    expBinBoundaries = [ Commands.InterpretPOI( yp.GetName() )[2][0] for yp in yieldParameters ]
    expBinBoundaries = [ float(b) for b in expBinBoundaries if isinstance(b, float) ]
    expBinBoundaries.append( 999. )

    print '\nSM Couplings:'
    for coupling in couplings:
        print '    {0:5}: {1}'.format( coupling.GetName(), coupling.getVal() )

    print 'SM yieldParameters:'
    for yieldParameter in yieldParameters:
        print '    {0:20}: {1}'.format( yieldParameter.GetName(), yieldParameter.getVal() )



    containers = []
    # yPerCoupling = []
    for newcouplings in testcouplings:

        # print '\nSetting ct = {0}, cg = {1}'.format( newcouplings['ct'], newcouplings['cg'] )
        print '\nSetting:'
        for coupling in couplings:
            couplingName = coupling.GetName()
            print '    {0} = {1}'.format( couplingName, newcouplings[couplingName] )

        print '\nCouplings:'
        for coupling in couplings:
            coupling.setVal( newcouplings[coupling.GetName()] )
            print '    {0:5}: {1}'.format( coupling.GetName(), coupling.getVal() )

        print 'yieldParameters:'
        y = []
        for yieldParameter in yieldParameters:
            print '    {0:20}: {1}'.format( yieldParameter.GetName(), yieldParameter.getVal() )
            y.append( yieldParameter.getVal() )

        print 'parametrizations:'
        yParametrization = []
        for parametrization in parametrizations:
            print '    {0:20}: {1}'.format( parametrization.GetName(), parametrization.getVal() )
            yParametrization.append( parametrization.getVal() )

        # yPerCoupling.append( ( y, yParametrization ) )

        container = Container()
        container.mus_expBinning       = y
        container.mus_expBinBoundaries = expBinBoundaries
        container.mus_theoryBinning    = yParametrization

        containers.append( container )


    # return yPerCoupling
    return containers





########################################
# Stewart-Tackmann business
########################################

# Takes either:
# - 1 argument, which is a TGraph object
# - 1 argument, which is a Container object created by ReadDerivedTheoryFile
def GetStewartTackmannCovarianceMatrix(
    theoryTg,
    scaleToRate = True,
    ):

    if isinstance( theoryTg, Container ):
        # Rename attributes so that it's consistent with TGraphs
        oldContainer = theoryTg
        theoryTg = Container(
            name          = 'container',
            binBoundaries = oldContainer.binBoundaries,
            SM            = oldContainer.crosssection,
            SM_up         = oldContainer.crosssection_up,
            SM_down       = oldContainer.crosssection_down,
            )


    # From Yellow Report 4, N3LO ggH inclusive
    # Errors are ONLY SCALE errors
    ggHinclusive_xs           = 48.58
    ggHinclusive_xs_errup     = 0.10
    ggHinclusive_xs_errdown   = 1.15
    ggHinclusive_xs_errsymm   = 0.5*( ggHinclusive_xs_errup + ggHinclusive_xs_errdown )
    ggHinclusive_xs_ratioup   = 0.0021
    ggHinclusive_xs_ratiodown = 0.0237

    # patterns = [
    #     r'^cg_[mp\d]+$',
    #     r'^ct_[mp\d]+$',
    #     r'^cb_[mp\d]+$',
    #     r'ct_[mp\d]+_cb_[mp\d]+_cg_[mp\d]+',
    #     r'ct_[mp\d]+_cb_([mp\d]+)$', # ct_XX_cb_XX , but not ct_XX_cb_XX_cg_XX
    #     r'ct_[mp\d]+_cg_[mp\d]+',
    #     ]
    # theoryTgs = TheoryCommands.LoadTheoryCurves( r'ct_[mp\d]+_cb_([mp\d]+)$' )
    # theoryTg = theoryTgs[0]

    print '[info] Processing {0}'.format( theoryTg.name )


    print ''
    print '[info] Be aware there is no 1/2.27 scaling here'
    print '[fixme] Forcing SM now'
    xs_theoryBinning = theoryTg.SM

    if not hasattr( theoryTg, 'SM_up' ):
        xs_up_theoryBinning   = [ i*j for i, j in zip( xs_theoryBinning, theoryTg.SM_up_ratio ) ]
        xs_down_theoryBinning = [ i*j for i, j in zip( xs_theoryBinning, theoryTg.SM_down_ratio ) ]
    else:
        xs_up_theoryBinning   = theoryTg.SM_up
        xs_down_theoryBinning = theoryTg.SM_down


    print '\n[info] Mapping fine theory binning to the experimental bins'
    print '[fixme] expBinBoundaries now hardcoded here'
    expBinBoundaries = [ 0., 15., 30., 45., 85., 125., 200., 350., 1000. ]
    expNBins = len(expBinBoundaries)-1
    expBinWidths = [ expBinBoundaries[i+1]-expBinBoundaries[i] for i in xrange(expNBins-1) ] + [ 1. ]

    def mapTheoryToExp( theoryBinBoundaries, theoryBinValues, expBinBoundaries ):
        theoryIntegralFunction = GetIntegral( theoryBinBoundaries, theoryBinValues )
        expBinValues = []
        for iBinExp in xrange(len(expBinBoundaries)-1):
            expBinValues.append(
                theoryIntegralFunction( expBinBoundaries[iBinExp], expBinBoundaries[iBinExp+1] ) / ( expBinBoundaries[iBinExp+1] - expBinBoundaries[iBinExp] )
                )
        return expBinValues

    xs       = mapTheoryToExp( theoryTg.binBoundaries, xs_theoryBinning, expBinBoundaries )
    xs_up    = mapTheoryToExp( theoryTg.binBoundaries, xs_up_theoryBinning, expBinBoundaries )
    xs_down  = mapTheoryToExp( theoryTg.binBoundaries, xs_down_theoryBinning, expBinBoundaries )

    xs_symmErrs = [ 0.5*( abs(xs[iBinExp]-xs_up[iBinExp]) + abs(xs[iBinExp]-xs_down[iBinExp]) ) for iBinExp in xrange(expNBins) ]
    xs_symmErrsInclusive = [ err*binWidth for err, binWidth in zip( xs_symmErrs, expBinWidths ) ]
    xs_symmErrsRelative  = [ err/rate for err, rate in zip( xs_symmErrs, xs ) ]

    print '[info] Found the following theory spectrum in experimental binning:'
    for iBinExp in xrange(len(expBinBoundaries)-1):
        print '    Bin {0} ( {1:6.1f} to {2:6.1f} ):  xs = {3:8.3f}  +/- {4:8.3f}   ({5:8.3f} inc.)'.format(
            iBinExp,
            expBinBoundaries[iBinExp],
            expBinBoundaries[iBinExp+1],
            xs[iBinExp],
            xs_symmErrs[iBinExp],
            xs_symmErrsInclusive[iBinExp]
            )


    # Build functions that evaluate the integral of spectrum in exp binning
    xs_integralFunction       = GetIntegral( expBinBoundaries, xs )
    xs_up_integralFunction    = GetIntegral( expBinBoundaries, xs_up )
    xs_down_integralFunction  = GetIntegral( expBinBoundaries, xs_down )





    xs_total = xs_integralFunction( 0, 1000 )
    xs_err_total = 0.5*(
        abs( xs_up_integralFunction( 0, 1000 ) - xs_total )
        + abs( xs_down_integralFunction( 0, 1000 ) - xs_total )
        )

    print '\n[info] Found total inclusive cross section: {0:9.4f}   +/- {1}'.format( xs_total, xs_err_total )
    print '[info] N3LO inclusive prediction is:        {0:9.4f}   +/- {1}'.format( ggHinclusive_xs, ggHinclusive_xs_errsymm )

    print ''


    offDiagonal = []
    onDiagonal  = []

    previous_sigma_GEcut_err = None
    for iCut, ptcut in enumerate(expBinBoundaries[1:-1]):

        sigma_GEcut = xs_integralFunction( ptcut, 1000 )
        sigma_GEcut_err = 0.5*(
            abs( xs_up_integralFunction( ptcut, 1000 ) - sigma_GEcut )
            + abs( xs_down_integralFunction( ptcut, 1000 ) - sigma_GEcut )
            )

        print '    offDiagonal_{0:<2}: XS( pt >= {1:7.2f} ) = {2:8.4f}  +/- {3:8.4f}'.format(
            iCut, ptcut, sigma_GEcut, sigma_GEcut_err )


        if previous_sigma_GEcut_err == None:
            Delta_iCut     = sqrt( xs_err_total**2 + sigma_GEcut_err**2 )
        else:
            Delta_iCut     = sqrt( previous_sigma_GEcut_err**2 + sigma_GEcut_err**2 )       

        print '     onDiagonal_{0:<2}: Delta_{1:<3} = {2:8.4f}                 (for pt = {3:6.2f} to {4:6.2f} )'.format(
            iCut, iCut, Delta_iCut, expBinBoundaries[iCut], expBinBoundaries[iCut+1] )

        previous_sigma_GEcut_err = sigma_GEcut_err

        offDiagonal.append( -sigma_GEcut_err**2 )
        onDiagonal.append(  Delta_iCut**2 )

    # >LastCut is the overflow
    onDiagonal.append( sigma_GEcut_err**2 )


    print '\n[info] Building covariance matrix:'
    covMat = [ [] for i in xrange(expNBins) ]

    for iBinExp1 in xrange(expNBins):
        for iBinExp2 in xrange(expNBins):
            if iBinExp1 == iBinExp2:
                covMat[iBinExp1].append( onDiagonal[iBinExp1] )
            elif abs(iBinExp1 - iBinExp2) == 1:
                covMat[iBinExp1].append( offDiagonal[ min(iBinExp1,iBinExp2) ] )
            else:
                covMat[iBinExp1].append( 0 )

    # Remove first row and column
    # covMat = [ [ covMat[iBinExp1][iBinExp2] for iBinExp2 in xrange(1,expNBins) ] for iBinExp1 in xrange(1,expNBins) ]

    numpy.set_printoptions( precision=3, linewidth=120 )
    print numpy.matrix(covMat)


    if scaleToRate:
        print '\n[info] Dividing by uncertainties (turns into a correlation matrix(:'

        for iBinExp1 in xrange(expNBins):
            for iBinExp2 in xrange(expNBins):
                if abs(iBinExp1-iBinExp2) > 1:
                    covMat[iBinExp1][iBinExp2] = 0
                else:
                    covMat[iBinExp1][iBinExp2] /= ( xs_symmErrsInclusive[iBinExp1] * xs_symmErrsInclusive[iBinExp2] )
        print numpy.matrix(covMat)

        print '\n[info] Scaling each row to it\'s diagonal element:'
        for iBinExp1 in xrange(expNBins):
            diagonalElement = covMat[iBinExp1][iBinExp1]
            for iBinExp2 in xrange(expNBins):
                covMat[iBinExp1][iBinExp2] /= diagonalElement
        print numpy.matrix(covMat)

        print '\n[info] Turning each row into 1+ratio of uncertainty:'
        for iBinExp1 in xrange(expNBins):
            for iBinExp2 in xrange(expNBins):
                if abs(iBinExp1-iBinExp2) > 1:
                    continue
                covMat[iBinExp1][iBinExp2] *= xs_symmErrsRelative[iBinExp1]
                covMat[iBinExp1][iBinExp2] += 1
        print numpy.matrix(covMat)

        print '\n[info] Killing spikes'
        for iBinExp1 in xrange(expNBins):
            for iBinExp2 in xrange(expNBins):
                if abs(iBinExp1-iBinExp2) > 1:
                    continue
                if covMat[iBinExp1][iBinExp2] < 0.:
                    covMat[iBinExp1][iBinExp2] = 0.
                elif covMat[iBinExp1][iBinExp2] > 10.:
                    covMat[iBinExp1][iBinExp2] = 10.
        print numpy.matrix(covMat)



    return covMat



# def AddCovarianceMatrixAsNuisanceParameters( datacard, covMat, verbose=True ):

#     signalprocesses, processes, bins = Commands.ListProcesses( datacard )


#     # Get the actual process line to which the nuisance parameter line should correspond
#     with open( datacard, 'r' ) as datacardFp:
#         for processLine in datacardFp.readlines():
#             if processLine.startswith('process'):
#                 break
#     datacardProcessLine = [ i.strip() for i in processLine.split() ][1:]


#     nuisPars = []    
#     for iProcess, process in enumerate(signalprocesses):

#         nuisPar = [ 'theoryUnc_' + process, 'lnN' ]

#         for inlineProcess in datacardProcessLine:

#             if inlineProcess in signalprocesses:
#                 iInlineProcess = signalprocesses.index(inlineProcess)
#                 if abs( iInlineProcess - iProcess ) <= 1:
#                     nuisPar.append( '{0:.3f}'.format(covMat[iProcess][iInlineProcess]) )
#                 else:
#                     nuisPar.append( '-' )    
#             else:
#                 nuisPar.append( '-' )

#         nuisPar = ' '.join( nuisPar )
#         nuisPars.append( nuisPar )
#         if verbose: print nuisPar


#     datacardAppended = datacard.replace( '.txt', '_theoryUncertainties.txt' )
#     with open( datacardAppended, 'w' ) as datacardAppendedFp:

#         with open( datacard, 'r' ) as datacardFp:
#             datacardText = datacardFp.read()

#         datacardText = re.sub( r'kmax ([0-9]+)', 'kmax *', datacardText )

#         datacardAppendedFp.write( datacardText + '\n\n' )

#         for nuisPar in nuisPars:
#             datacardAppendedFp.write( nuisPar + '\n' )



def AddCovarianceMatrixAsNuisanceParameters( datacard, covMat, totalProcessRates, verbose=True ):

    signalprocesses, processes, bins = Commands.ListProcesses( datacard )


    ########################################
    # Perform whitening of the covariance matrix
    ########################################

    covMat = numpy.array( covMat )
    nPars = covMat.shape[0]

    eigenValues, eigenVectors = numpy.linalg.eig(
        covMat
        )

    # Convert to ordinary lists
    eigenValues = list( eigenValues )
    eigenVectors = [ list(eigenVectors[:,i]) for i in xrange(eigenVectors.shape[1]) ]

    # Make eigenvectors unit length (should already be the case but to be sure)
    eigenVectors = [ normalize(eigenVector) for eigenVector in eigenVectors ]

    # Sort according to eigenvalues
    eigenValues, eigenVectors = ( list(t) for t in zip(*sorted( zip(eigenValues, eigenVectors), reverse=True )) )


    # ======================================
    # Calculate whitening matrix:

    oneOverSqrtD = numpy.diag( [ 1./sqrt(val) for val in eigenValues ] )
    E            = numpy.hstack( [ numpy_column(vec) for vec in eigenVectors ] )
    ET           = E.T

    #   E * D^-0.5 * E^T
    whiteningMatrix = E.dot( oneOverSqrtD.dot( ET ) )

    # whiteningMatrix maps correlated noise to white noise;
    # Need the inverse operation:
    inverseWhiteningMatrix = numpy.linalg( whiteningMatrix )

    print 'Map from uncorrelated to correlated:'
    print inverseWhiteningMatrix


    # ======================================
    # Divide by process normalization
    # This has to be done column-wise

    for iPar in xrange(nPars):
        totalProcessRate = totalProcessRates[iPar]
        for iCol in xrange(nPars):
            inverseWhiteningMatrix[iPar][iCol] /= totalProcessRate

    print inverseWhiteningMatrix



def normalize( l ):
    norm = sqrt(sum([ val**2 for val in l ]))
    l = [ val/norm for val in l ]
    return l

def numpy_column( l ):
    return numpy.array(l)[:,numpy.newaxis]

def numpy_row( l ):
    return numpy.array(l)[:,numpy.newaxis]    




# Takes either:
# - 1 argument, which is a TGraph object
# - 1 argument, which is a Container object created by ReadDerivedTheoryFile
def PrepareRebinnedTheory(
    theoryTg,
    ):

    if isinstance( theoryTg, Container ):
        # Rename attributes so that it's consistent with TGraphs
        oldContainer = theoryTg
        theoryTg = Container(
            name          = 'container',
            binBoundaries = oldContainer.binBoundaries,
            SM            = oldContainer.crosssection,
            SM_up         = oldContainer.crosssection_up,
            SM_down       = oldContainer.crosssection_down,
            )

    output = Container( name = 'rebinnedTheory' )


    # From Yellow Report 4, N3LO ggH inclusive
    # Errors are ONLY SCALE errors
    output.ggHinclusive_xs           = 48.58
    output.ggHinclusive_xs_errup     = 0.10
    output.ggHinclusive_xs_errdown   = 1.15
    output.ggHinclusive_xs_errsymm   = 0.5*( output.ggHinclusive_xs_errup + output.ggHinclusive_xs_errdown )
    output.ggHinclusive_xs_ratioup   = 0.0021
    output.ggHinclusive_xs_ratiodown = 0.0237

    # patterns = [
    #     r'^cg_[mp\d]+$',
    #     r'^ct_[mp\d]+$',
    #     r'^cb_[mp\d]+$',
    #     r'ct_[mp\d]+_cb_[mp\d]+_cg_[mp\d]+',
    #     r'ct_[mp\d]+_cb_([mp\d]+)$', # ct_XX_cb_XX , but not ct_XX_cb_XX_cg_XX
    #     r'ct_[mp\d]+_cg_[mp\d]+',
    #     ]
    # theoryTgs = TheoryCommands.LoadTheoryCurves( r'ct_[mp\d]+_cb_([mp\d]+)$' )
    # theoryTg = theoryTgs[0]

    print '[info] Processing {0}'.format( theoryTg.name )


    print ''
    print '[info] Be aware there is no 1/2.27 scaling here'
    print '[fixme] Forcing SM now'
    xs_theoryBinning = theoryTg.SM

    if not hasattr( theoryTg, 'SM_up' ):
        xs_up_theoryBinning   = [ i*j for i, j in zip( xs_theoryBinning, theoryTg.SM_up_ratio ) ]
        xs_down_theoryBinning = [ i*j for i, j in zip( xs_theoryBinning, theoryTg.SM_down_ratio ) ]
    else:
        xs_up_theoryBinning   = theoryTg.SM_up
        xs_down_theoryBinning = theoryTg.SM_down


    print '\n[info] Mapping fine theory binning to the experimental bins'
    print '[fixme] expBinBoundaries now hardcoded here'
    expBinBoundaries = [ 0., 15., 30., 45., 85., 125., 200., 350., 1000. ]
    expNBins = len(expBinBoundaries)-1
    expBinWidths = [ expBinBoundaries[i+1]-expBinBoundaries[i] for i in xrange(expNBins-1) ] + [ 1. ]

    def mapTheoryToExp( theoryBinBoundaries, theoryBinValues, expBinBoundaries ):
        theoryIntegralFunction = GetIntegral( theoryBinBoundaries, theoryBinValues )
        expBinValues = []
        for iBinExp in xrange(len(expBinBoundaries)-1):
            expBinValues.append(
                theoryIntegralFunction( expBinBoundaries[iBinExp], expBinBoundaries[iBinExp+1] ) / ( expBinBoundaries[iBinExp+1] - expBinBoundaries[iBinExp] )
                )
        return expBinValues

    output.xs       = MapFineToCoarse( theoryTg.binBoundaries, xs_theoryBinning, expBinBoundaries )
    output.xs_up    = MapFineToCoarse( theoryTg.binBoundaries, xs_up_theoryBinning, expBinBoundaries )
    output.xs_down  = MapFineToCoarse( theoryTg.binBoundaries, xs_down_theoryBinning, expBinBoundaries )
    output.xs_symmErrs          = [ 0.5*( abs(output.xs[iBinExp]-output.xs_up[iBinExp]) + abs(output.xs[iBinExp]-output.xs_down[iBinExp]) ) for iBinExp in xrange(expNBins) ]
    output.xs_symmErrsInclusive = [ err*binWidth for err, binWidth in zip( output.xs_symmErrs, expBinWidths ) ]
    output.xs_symmErrsRelative  = [ err/rate for err, rate in zip( output.xs_symmErrs, output.xs ) ]

    output.binBoundaries    = expBinBoundaries
    output.expBinBoundaries = expBinBoundaries
    output.expNBins         = expNBins
    output.expBinWidths     = expBinWidths

    output.originalTheory = Container(
        name          = 'originalTheory',
        binBoundaries = theoryTg.binBoundaries,
        SM            = theoryTg.SM,
        SM_up         = theoryTg.SM_up,
        SM_down       = theoryTg.SM_down,
        xs            = theoryTg.SM,
        xs_up         = theoryTg.SM_up,
        xs_down       = theoryTg.SM_down,
        )

    return output


def MakeRebinnedTheoryPlot(
        theoryTg
        ):

    c.Clear()
    SetCMargins()

    rebinnedContainer = PrepareRebinnedTheory( theoryTg )

    fineTheoryTg = GetTheoryTGraph(
        'fine',
        rebinnedContainer.originalTheory.binBoundaries,
        rebinnedContainer.originalTheory.xs,
        rebinnedContainer.originalTheory.xs_down,
        rebinnedContainer.originalTheory.xs_up,
        boundaries = True
        )
    fineTheoryTg.SetLineColor(2)

    coarseTheoryTg = GetTheoryTGraph(
        'coarse',
        rebinnedContainer.binBoundaries,
        rebinnedContainer.xs,
        rebinnedContainer.xs_down,
        rebinnedContainer.xs_up,
        boundaries = True
        )
    coarseTheoryTg.SetLineColor(4)

    base = GetPlotBase(
        xMax = max( max(rebinnedContainer.originalTheory.binBoundaries[:-1]), max(rebinnedContainer.binBoundaries[:-1]) ),
        yMax = max( max(rebinnedContainer.originalTheory.xs_up), max(rebinnedContainer.xs_up) ),
        xTitle = 'pT [GeV]', yTitle='d#sigma/dp_{T} [pb/GeV]'
        )
    base.Draw()

    coarseTheoryTg.Draw( 'P' )

    leg = ROOT.TLegend( 1-RightMargin-0.5, 1-TopMargin-0.25, 1-RightMargin, 1-TopMargin )
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    leg.AddEntry( coarseTheoryTg.GetName(), 'Theory spectrum in exp. bins', 'l' )
    leg.Draw()

    SaveC( 'RebinnedTheoryPlot_onlyCoarse' )

    fineTheoryTg.Draw(   'P' )
    leg.AddEntry( fineTheoryTg.GetName(),   'Given theory spectrum', 'l' )
    leg.Draw()

    SaveC( 'RebinnedTheoryPlot' )



########################################
# Execute this if called directly
########################################

if __name__ == '__main__':
    print 'Import functions from this file'
    print 'Available functions:'
    functions = [ name for (name, thing) in locals().items() if callable(thing) ]
    print '\n'.join( functions )