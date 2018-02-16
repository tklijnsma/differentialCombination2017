#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import Commands
from Container import Container

import os, tempfile, shutil, re, glob, itertools, sys, numpy, operator, pprint, re
from os.path import *
from operator import itemgetter
from array import array
from math import log, exp, sqrt, copysign
from copy import deepcopy

# Needed for the observable class
sys.path.append('..')
import LatestPaths

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
def set_cmargins(
    LeftMargin   = 0.15,
    RightMargin  = 0.03,
    BottomMargin = 0.15,
    TopMargin    = 0.03,
    for2Dhist    = False
    ):
    if for2Dhist:
        c.SetLeftMargin(   0.12 )
        c.SetRightMargin(  0.10 )
        c.SetBottomMargin( 0.12 )
        c.SetTopMargin(    0.09 )
    else:
        c.SetLeftMargin( LeftMargin )
        c.SetRightMargin( RightMargin )
        c.SetBottomMargin( BottomMargin )
        c.SetTopMargin( TopMargin )


PLOTDIR = 'plots_{0}'.format(datestr)
def set_plot_dir( newdir ):
    global PLOTDIR
    PLOTDIR = newdir

SAVEROOT = False
def save_as_root( newvalue=True ):
    global SAVEROOT
    SAVEROOT = newvalue
SAVEPNG = False
def save_as_png( newvalue=True ):
    global SAVEPNG
    SAVEPNG = newvalue
SAVEPNG_THROUGH_CONVERT = False
def save_as_png_through_convert( newvalue=True ):
    global SAVEPNG_THROUGH_CONVERT
    SAVEPNG_THROUGH_CONVERT = newvalue

def save_c( outname, asPNG=False, asROOT=False ):
    global PLOTDIR
    if not isdir(PLOTDIR): os.makedirs(PLOTDIR)

    subdir = ''
    if len(outname.split('/')) == 2:
        subdir = outname.split('/')[0]
        if not isdir( join( PLOTDIR, subdir ) ): os.makedirs( join( PLOTDIR, subdir ) )

    outname = join( PLOTDIR, subdir, basename(outname).replace('.pdf','').replace('.png','') )
    c.SaveAs( outname + '.pdf' )
    if asPNG or SAVEPNG:
        c.SaveAs( outname + '.png' )
    if asROOT or SAVEROOT:
        c.SaveAs( outname + '.root' )

    if SAVEPNG_THROUGH_CONVERT:
        # See: https://stackoverflow.com/a/6605085/9209944
        cmd = 'convert -density 300 -quality 100 {0}.pdf {0}.png'.format(outname)
        Commands.execute_command(cmd)


ROOTCOUNTER = 1000
def get_unique_root_name():
    global ROOTCOUNTER
    name = 'root{0}_{1}'.format( ROOTCOUNTER, Commands.__uniqueid__().next() )
    ROOTCOUNTER += 1
    return name


def get_plot_base(
    xMin = 0, xMax = 1,
    yMin = 0, yMax = 1,
    xTitle = 'x', yTitle = 'y',
    SetTitleSizes = True,
    ):

    base = ROOT.TH1F()
    ROOT.SetOwnership( base, False )
    base.SetName( get_unique_root_name() )
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
# Derived Theory Containers
########################################

#____________________________________________________________________
def rebin_derived_theory_container( container, newBinBoundaries, lastBinIsOverflow=True, verbose=False ):

    newNBins = len(newBinBoundaries) - 1

    newCrosssection = rebin(
        theoryBinBoundaries = container.binBoundaries,
        theoryBinValues     = container.crosssection,
        expBinBoundaries    = newBinBoundaries,
        lastBinIsOverflow   = lastBinIsOverflow,
        verbose             = verbose,
        )

    newRatios = rebin(
        theoryBinBoundaries = container.binBoundaries,
        theoryBinValues     = container.ratios,
        expBinBoundaries    = newBinBoundaries,
        lastBinIsOverflow   = lastBinIsOverflow,
        verbose             = verbose,
        )

    # newBinCenters = [ 0.5*(newBinBoundaries[i]+newBinBoundaries[i+1]) for i in xrange(newNBins) ]

    container.binBoundaries = newBinBoundaries
    container.crosssection  = newCrosssection
    container.ratios        = newRatios


#____________________________________________________________________
# Theory files contain a somewhat non-consistent binning, turn it into a well-defined binning
def binning_heuristic(
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
def get_theory_TGraph(
    name,
    ptPoints,
    muPoints,
    muBoundLeft   = None,
    muBoundRight  = None,
    boundaries    = False,
    ):
    
    print 'WARNING: Better to use OutputInterface.OutputContainer.get_TGraph()'

    if not boundaries:
        # Default setting; supplied pt points are bin centers

        if not len(ptPoints) == len(muPoints):
            Commands.throw_error(
                'Length of input lists are not set right\n'
                '    len(ptPoints) = {0} is not equal to len(muPoints) = {1}'.format(
                    len(ptPoints), len(muPoints) )
                )
            return

        # ======================================
        # Binning heuristic for pt points

        nBins = len(ptPoints)

        binCenters, binBoundaries, binWidths = binning_heuristic( ptPoints )
        halfBinWidths = [ 0.5*w for w in binWidths ]

    else:
        # Supplied pt points are bin boundaries

        if not len(ptPoints)-1 == len(muPoints):
            Commands.throw_error(
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
def get_integral( binBoundaries, binValues ):

    if len(binBoundaries)-1 != len(binValues):
        # more_or_less = 'more' if len(binValues) > len(binBoundaries)-1 else 'less'
        # Commands.warning( 'Found {0} binValues[{1}] than bins[{2}]; will limit to {3} bins' )
        # Actually, just throw an error; the chance that this is intended is really slim
        Commands.throw_error( 'Found len(binBoundaries) = {0}, but len(binValues) = {1}'.format( len(binBoundaries), len(binValues) ) )

    nBins = len(binValues)
    binCenters = [ 0.5*( binBoundaries[i+1] + binBoundaries[i] ) for i in xrange(nBins-1) ]

    # linearInterpolate = lambda a, x1, x2, y1, y2: \
    #     y1 + ( a - x1 ) / ( x2 - x1 ) * ( y2 - y1 )

    def integralfunction( a, b, verbose=False ):

        if verbose: print '\n\nInterpolation function called with a = {0} and b = {1} (defined range: {2} to {3})'.format( a, b, binBoundaries[0], binBoundaries[-1] )

        aInterpolated = False
        bInterpolated = False

        if a < binBoundaries[0] and b < binBoundaries[0]:
            return 0.0
        elif a > binBoundaries[-1] and b > binBoundaries[-1]:
            return 0.0

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
        if verbose: print 'Looking for left and right bound in binBoundaries'
        for iBin in xrange(nBins):

            if aInterpolated and bInterpolated: break

            if verbose: print '\niBin = {0}'.format(iBin)
            if verbose: print '( a[{0}] >= binBoundaries[iBin][{1}] and a[{0}] < binBoundaries[iBin+1][{2}] ) == {3}'.format(
                a, binBoundaries[iBin], binBoundaries[iBin+1], ( a >= binBoundaries[iBin] and a < binBoundaries[iBin+1] )
                )
            if not aInterpolated and ( a >= binBoundaries[iBin] and a < binBoundaries[iBin+1] ):
                ya = binValues[iBin]
                ia = iBin
                aInterpolated = True

            if verbose: print '( b[{0}] > binBoundaries[iBin][{1}] and b[{0}] <= binBoundaries[iBin+1][{2}] ) == {3}'.format(
                b, binBoundaries[iBin], binBoundaries[iBin+1], ( b > binBoundaries[iBin] and b <= binBoundaries[iBin+1] )
                )
            if not bInterpolated and ( b > binBoundaries[iBin] and b <= binBoundaries[iBin+1] ):
                yb = binValues[iBin]
                ib = iBin
                bInterpolated = True

        if not aInterpolated:
            Commands.throw_error( 'Problem determining integration range for left bound' )

        if not bInterpolated:
            Commands.throw_error( 'Problem determining integration range for right bound' )

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
def rebin(
        theoryBinBoundaries,
        theoryBinValues,
        expBinBoundaries,
        lastBinIsOverflow = False,
        verbose = False,
        ):
    # Commands.warning( 'Calling TheoryCommands.rebin(); This code needs to be checked!!' )

    nBinsExp = len(expBinBoundaries)-1

    if verbose: print 'Rebinning to {0}; lastBinIsOverflow = {1}'.format( expBinBoundaries, lastBinIsOverflow )

    theoryIntegralFunction = get_integral( theoryBinBoundaries, theoryBinValues )
    expBinValues = []
    for iBinExp in xrange(nBinsExp-1):
        expBinValues.append(
            theoryIntegralFunction( expBinBoundaries[iBinExp], expBinBoundaries[iBinExp+1] ) / ( expBinBoundaries[iBinExp+1] - expBinBoundaries[iBinExp] )
            )
        if verbose: print '  Doing integral({0:7.2f}, {1:7.2f}) / ({1:7.2f}-{0:7.2f}): {2} (/GeV)'.format( expBinBoundaries[iBinExp], expBinBoundaries[iBinExp+1], expBinValues[iBinExp] )


    # ======================================
    # Last bin needs some special care

    if lastBinIsOverflow:
        expBinValues.append( theoryIntegralFunction( expBinBoundaries[-2], theoryBinBoundaries[-1] ) )
        Commands.warning( 'Last bin uses the rightmost theory bin bound {0} instead of requested exp bin bound {1}'.format( theoryBinBoundaries[-1], expBinBoundaries[-1] ) )
        if verbose: print '  Doing integral({0:7.2f}, {1:7.2f})                      : {2} (not per GeV!)'.format( expBinBoundaries[-2], theoryBinBoundaries[-1], expBinValues[-1] )

    else:
        if expBinBoundaries[-1] > theoryBinBoundaries[-1]:
            right = theoryBinBoundaries[-1]
            Commands.warning( 'Last bin uses the rightmost theory bin bound {0} instead of requested exp bin bound {1}'.format( theoryBinBoundaries[-1], expBinBoundaries[-1] ) )
        else:
            right = expBinBoundaries[-1]
        # expBinBoundaries[-1] = right
        expBinValues.append( theoryIntegralFunction( expBinBoundaries[-2], right ) / ( right - expBinBoundaries[-2] ) )
        if verbose: print '  Doing integral({0:7.2f}, {1:7.2f}) / ({1:7.2f}-{0:7.2f}): {2} (/GeV)'.format( expBinBoundaries[-2], right, expBinValues[-1] )

    return expBinValues

#____________________________________________________________________
def merge_bins( some_binned_data, binMerging, binning=None ):
    """
    Use e.g. "[ 1, 2, [3,4], [5,6] ]" for the binMerging. The order matters!
    If binning is given, the binned data is assumed to be *density*, rather than an inclusive number
    """
    ret = []
    for bins in binMerging:
        if isinstance( bins, int ):
            # Case for only 1 bin; simply copy
            ret.append( some_binned_data[bins] )
        else:
            sum = 0.
            if binning is None:
                for iBin in bins:
                    sum += some_binned_data[iBin]
            else:
                for iBin in bins:
                    sum += some_binned_data[iBin] * (binning[iBin+1]-binning[iBin])
                sum /= ( binning[ bins[-1]+1 ] - binning[0] )
            ret.append(sum)                
    return ret

#____________________________________________________________________
# Not sure which funcion uses this
def get_short_theory_name( names ):

    print 'Notification: TheoryCommands.get_short_theory_name() was called'

    keyStrings = [ 'ct', 'cb', 'cg' ]

    uniqueCombinations = set()
    for name in names:
        presentKeyStrings = ''
        for keyString in keyStrings:
            if keyString in name: presentKeyStrings += keyString
        uniqueCombinations.add( presentKeyStrings )
    return '_'.join( list(uniqueCombinations) )


#____________________________________________________________________
def get_xyfrom_TGraph( Tg ):
    N = Tg.GetN()
    xs = []
    ys = []
    x_Double = ROOT.Double(0)
    y_Double = ROOT.Double(0)

    for i in xrange(N):
        Tg.GetPoint( i, x_Double, y_Double )
        xs.append( float(x_Double) )
        ys.append( float(y_Double) )

    return xs, ys


########################################
# Helper functions for plotting
########################################

#____________________________________________________________________
def basic_read_scan(
        rootfiles,
        xAttr,
        yAttr = 'deltaNLL',
        ):
    
    containers = Commands.convert_TChain_to_array(
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
def write_TH2_to_file(
        rootfiles,
        xCoupling = 'ct',
        yCoupling = 'cg',
        verbose = True,
        ):

    scan = Commands.convert_TChain_to_array(
        rootfiles
        )
    nPoints = len(scan[xCoupling])


    keys = [ yCoupling, xCoupling, 'deltaNLL' ]
    ret = []
    for iPoint in xrange(nPoints):
        ret.append( [ scan[key][iPoint] for key in keys ] )
    if verbose:
        pprint.pprint( [ keys ] + ret )


    def infer_bin_boundaries( binCenters ):
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
    xBinBoundaries = infer_bin_boundaries( xBinCenters )

    yBinCenters = list(set(scan[yCoupling]))
    yBinCenters.pop( yBinCenters.index(yBestfit) )
    yBinCenters.sort()
    yNBins = len(yBinCenters)
    yBinBoundaries = infer_bin_boundaries( yBinCenters )


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
def get_contours_from_TH2( TH2_original, threshold, verbose=True ):

    # Open a temporary canvas so the contour business does not screw up other plots
    ctemp = ROOT.TCanvas( 'ctemp', 'ctemp', 1000, 800 )
    ctemp.cd()

    TH2 = TH2_original.Clone()
    TH2.SetName( get_unique_root_name() )
    TH2.SetContour( 1, array( 'd', [threshold] ) )

    if verbose:
        print 'Trying to get contours from \'{0}\''.format( TH2.GetName() ) + ( ' ({0})'.format( TH2_original.name ) if hasattr( TH2_original, 'name' ) else '' )

    TH2.Draw( 'CONT Z LIST' )
    ROOT.gPad.Update()

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
        print '    {0} contours found'.format(len(Tgs))
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
def set_extrema_of_contour( Tg ):

    xBuffer = Tg.GetX()
    xs = [ xBuffer[i] for i in xrange( Tg.GetN() ) ]

    yBuffer = Tg.GetY()
    ys = [ yBuffer[i] for i in xrange( Tg.GetN() ) ]

    Tg.xMin = min(xs)
    Tg.xMax = max(xs)
    Tg.yMin = min(ys)
    Tg.yMax = max(ys)


#____________________________________________________________________
def get_TH2_from_list_of_root_files(
        rootfiles,
        xCoupling = 'ct',
        yCoupling = 'cg',
        verbose   = False,
        xMin = None, xMax = None, yMin = None, yMax = None,
        multiplyByTwo = True,
        zVariable = 'deltaNLL',
        # defaultHValue = None,
        defaultHValue = 999.,
        refind_minimum_if_dnll_negative = False,
        ignore_deltaNLL_negativity = False,
        multiplier = None
        ):

    # Read values from specified rootfiles
    scan = Commands.convert_TChain_to_array(
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


    def infer_bin_boundaries( binCenters ):
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
    xBinBoundaries = infer_bin_boundaries( xBinCenters )

    yBinCenters = list(set( GetListFromScan(yCoupling) ))
    yBinCenters.pop( yBinCenters.index(yBestfit) )
    yBinCenters.sort()
    yNBins = len(yBinCenters)
    yBinBoundaries = infer_bin_boundaries( yBinCenters )


    # When looping over the points, don't consider the bestfit point a part of the histogram
    skipBestfit = True
    old_xBestfit = -99999.
    old_yBestfit = -99999.

    # Find actual minimum in all deltaNLL values
    minDeltaNLL = min(GetListFromScan('deltaNLL'))
    if minDeltaNLL < 0.0:
        scandir_heuristicly = basename(dirname(abspath(rootfiles[0])))
        Commands.warning( '[dir={0}] min deltaNLL = {1} - the start fit was NOT the best fit!!'.format(scandir_heuristicly, minDeltaNLL) )
        if refind_minimum_if_dnll_negative:
            # Save old bestfit values so they can be skipped when making the histogram
            old_xBestfit = xBestfit
            old_yBestfit = yBestfit
            iBestfit = GetListFromScan('deltaNLL').index( minDeltaNLL )
            xBestfit = GetListFromScan(xCoupling)[iBestfit]
            yBestfit = GetListFromScan(yCoupling)[iBestfit]
            Commands.warning( 'Will take the ACTUAL minimum to be deltaNLL={0} ({1}={2}, {3}={4}), and will shift histogram UP by {0}'.format(
                minDeltaNLL, xCoupling, xBestfit, yCoupling, yBestfit
                ))
            for point in scan:
                point['deltaNLL'] -= minDeltaNLL
        elif minDeltaNLL < -0.01:
            if ignore_deltaNLL_negativity:
                Commands.warning('Ignoring the fact that some points are negative!')
            else:
                Commands.throw_error('The minimum deltaNLL threshold is -0.01, but minDeltaNLL = {0}'.format(minDeltaNLL))

    H2name = get_unique_root_name()
    H2 = ROOT.TH2D(
        H2name, '',
        xNBins, array( 'd', xBinBoundaries ),
        yNBins, array( 'd', yBinBoundaries ),
        )
    ROOT.SetOwnership( H2, False )

    if not defaultHValue is None:
        for iX in xrange(xNBins):
            for iY in xrange(yNBins):
                H2.SetBinContent( iX+1, iY+1, defaultHValue )


    for iPoint in xrange(nPoints):

        x = scan[iPoint][xCoupling]
        y = scan[iPoint][yCoupling]

        if skipBestfit:
            if x == xBestfit and y == yBestfit:
                continue
            if x == old_xBestfit and y == old_yBestfit:
                continue

        try:
            iBinX = xBinCenters.index(x)
        except ValueError:
            print '[ERROR] Point {0} ({1}) not in list'.format( iPoint, x )
            continue

        try:
            iBinY = yBinCenters.index(y)
        except ValueError:
            print '[ERROR] Point {0} ({1}) not in list'.format( iPoint, y )
            continue

        val = scan[iPoint][zVariable]
        if multiplyByTwo:
            val *= 2.
        if not(multiplier is None):
            val *= multiplier
        H2.SetBinContent( iBinX+1, iBinY+1, val )


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
def plot_coupling_scan_2d(
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

    res = get_TH2_from_list_of_root_files(
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
    set_cmargins(
        LeftMargin   = 0.12,
        RightMargin  = 0.10,
        BottomMargin = 0.12,
        TopMargin    = 0.09,
        )


    H2.Draw('COLZ')


    if drawContours:

        contourTgs_1sigma = get_contours_from_TH2( H2, 2.30 )
        for Tg_1sigma in contourTgs_1sigma:
            Tg_1sigma.SetLineColor(1)
            Tg_1sigma.SetLineWidth(3)
            Tg_1sigma.Draw('LSAME')
            Tg_1sigma.SetName('contour1sigma')

        contourTgs_2sigma = get_contours_from_TH2( H2, 6.18 )
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

    Tpoint = ROOT.TGraph( 1, array( 'd', [res.xBestfit] ), array( 'd', [res.yBestfit] ) )
    ROOT.SetOwnership( Tpoint, False )
    Tpoint.SetMarkerColor(1)
    Tpoint.SetMarkerSize(2)
    Tpoint.SetMarkerStyle(34)
    Tpoint.Draw('P')
    Tpoint.SetName('BestfitPoint')


    if drawContours:
        leg = ROOT.TLegend(
            c.GetRightMargin() + 0.08,
            0.99,
            c.GetRightMargin() + 0.08 + 0.7,
            0.99 - c.GetTopMargin() + 0.02
            )
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.set_n_columns(4)
        leg.AddEntry( Tg_1sigma.GetName(), '1 #sigma', 'l' )
        leg.AddEntry( Tg_2sigma.GetName(), '2 #sigma', 'l' )
        leg.AddEntry( Tpoint.GetName(),    'Best fit', 'p' )
        leg.AddEntry( Tpoint_SM.GetName(), 'SM',       'p' )
        leg.Draw()

    outname = 'couplingscan2D_{0}'.format( basename(dirname(rootfiles[0])).replace('/','') )
    if not zVariable == 'deltaNLL': outname = outname.replace( 'couplingscan2D_', 'couplingscan2D_{0}_'.format(zVariable) )
    save_c( outname )
    set_cmargins()

    return res
    

########################################
# Execute this if called directly
########################################

if __name__ == '__main__':
    print 'Import functions from this file'
    print 'Available functions:'
    functions = [ name for (name, thing) in locals().items() if callable(thing) ]
    print '\n'.join( functions )