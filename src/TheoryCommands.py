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
# Dealing with theory spectra from Agnieszka
########################################

def TheoryPath( filename ):
    return join( 'suppliedInput/theory', filename )

physicsFileDict = {
    'ct_1p1'                   : TheoryPath( 'ratio_ctup_new' ),
    'ct_0p9'                   : TheoryPath( 'ratio_ctdown_new' ),
    'cg_0p008'                 : TheoryPath( 'ratio_cgup_new' ),
    'cg_m0p008'                : TheoryPath( 'ratio_cgdown_new' ),
    'cb_4p0'                   : TheoryPath( 'ratio_cbup_new' ),
    'cb_m2p0'                  : TheoryPath( 'ratio_cbdown_new' ),
    'ct_1p2_cb_m2p98_cg_m0p03' : TheoryPath( 'ratio_cg003ct12_new' ),
    'ct_1p3_cb_m0p85_cg_m0p03' : TheoryPath( 'ratio_cg003ct13sw_new' ),
    'ct_1p4_cb_3p31_cg_m0p03'  : TheoryPath( 'ratio_cg003ct14sw_new' ),
    'ct_1p2_cb_m4p89_cg_m0p04' : TheoryPath( 'ratio_cg004ct12sw_new' ),
    'ct_1p3_cb_m3p34_cg_m0p04' : TheoryPath( 'ratio_cg004ct13sw_new' ),
    'ct_1p5_cb_1p88_cg_m0p04'  : TheoryPath( 'ratio_cg004ct15sw_new' ),
    'ct_1p4_cb_m3p67_cg_m0p05' : TheoryPath( 'ratio_cg005ct14sw_new' ),
    'ct_1p5_cb_m1p79_cg_m0p05' : TheoryPath( 'ratio_cg005ct15sw_new' ),
    'ct_0p5_cb_m7p46'          : TheoryPath( 'ratio_ctcb05_new' ),
    'ct_0p8_cb_m3p67'          : TheoryPath( 'ratio_ctcb08_new' ),
    'ct_0p9_cb_m1p79'          : TheoryPath( 'ratio_ctcb09_new' ),
    'ct_1p1_cb_3p79'           : TheoryPath( 'ratio_ctcb11_new' ),
    'ct_1p2_cb_4p67'           : TheoryPath( 'ratio_ctcb12_new' ),
    'ct_0p1_cg_0p075'          : TheoryPath( 'ratio_ctcg01_new' ),
    'ct_0p5_cg_0p042'          : TheoryPath( 'ratio_ctcg05_new' ),
    'ct_1p5_cg_m0p042'         : TheoryPath( 'ratio_ctcg15_new' ),
    'ct_2p0_cg_m0p083'         : TheoryPath( 'ratio_ctcg2_new' ),

    'SM_NLO'                   : TheoryPath( 'SM_NLO' ),
    'SM_NLO_downRatio'         : TheoryPath( 'SMmin_NLO' ),
    'SM_NLO_upRatio'           : TheoryPath( 'SMmax_NLO' ),
    'SM_NNLO'                  : TheoryPath( 'SM_NNLO' ),
    'SM_NNLO_downRatio'        : TheoryPath( 'SMmin_NNLO' ),
    'SM_NNLO_upRatio'          : TheoryPath( 'SMmax_NNLO' ),
    }



def ReadLinesOfYukawaTheoryFile( theoryFile, verbose=False ):
    # Read lines
    with open( theoryFile, 'r' ) as theoryFp:
        lines = [ line.strip() for line in theoryFp.readlines() ]
    commentlines = [ line.strip() for line in lines if line.startswith('#') ]
    lines        = [ line for line in lines if not line.startswith('#') and len(line) > 0 ]

    pts          = []
    matched_xss  = []
    resummed_xss = []

    for line in lines:
        components = line.split()
        pt          = float(components[0])
        matched_xs  = float(components[1])
        resummed_xs = float(components[2])

        if verbose:
            print '    pt = {0:<8.3f} |  matched_xs = {1:<10.6f} |  resummed_xs = {2:<10.6f}'.format(
                pt, matched_xs, resummed_xs
                )

        pts.append( pt )
        matched_xss.append( matched_xs )
        resummed_xss.append( resummed_xs )

    return pts, matched_xss, resummed_xss


def CreateDerivedTheoryFiles_Yukawa(
        verbose = True,
        ):

    outdir = 'derivedTheoryFiles_{0}'.format( datestr )
    if not isdir( outdir ): os.makedirs( outdir )

    theoryDir = abspath( 'suppliedInput/theory_Yukawa' )
    theoryFiles = glob.glob( join( theoryDir, '*.res' ) )


    # ======================================
    # Try to find a 'SM' file first so it's possible to calculate ratios

    for theoryFile in theoryFiles:
        theoryFilename = basename( theoryFile )

        pat = r'H125-LHC13-R04-MSbar-xmur(?P<muR>\d+)-xmuf(?P<muF>\d+)_(?P<kappab>[\-\d]+)_(?P<kappac>[\-\d]+)-xQ(?P<Q>\d+)'
        match = re.search( pat, theoryFilename )

        if not match: continue

        muR    = float(match.group('muR')) / 50.
        muF    = float(match.group('muF')) / 50.
        Q      = float(match.group('Q')) / 50.
        kappab = float(match.group('kappab'))
        kappac = float(match.group('kappac'))

        if (
                muR == 1.0 and
                muF == 1.0 and
                Q == 1.0 and
                kappac == 1.0 and
                kappab == 1.0
                ):

            smFound = True
            pts, matched_xss, resummed_xss = ReadLinesOfYukawaTheoryFile( theoryFile, verbose )
            newBinCenters, binBoundaries, binWidths = BinningHeuristic( pts, manualSwitchAt50=False )

            SM = Container()
            SM.pts           = deepcopy(pts)
            SM.matched_xss   = deepcopy(matched_xss)
            SM.resummed_xss  = deepcopy(resummed_xss)
            SM.newBinCenters = deepcopy(newBinCenters)
            SM.binBoundaries = deepcopy(binBoundaries)
            SM.binWidths     = deepcopy(binWidths)
            del pts
            del matched_xss
            del resummed_xss
            del newBinCenters
            del binBoundaries
            del binWidths

            break
    else:
        print '[info] Did not find a SM file; ratios will be unavailable in the derivedTheoryFile'
        smFound = False


    # ======================================
    # Process all other files
    
    for theoryFile in theoryFiles:
        theoryFilename = basename( theoryFile )

        if verbose:
            print 'Processing theory file \'{0}\''.format( theoryFile )

        # H125-LHC13-R04-MSbar-xmur025-xmuf050_-2_-5-xQ050-NNLO+NNLLmult.-2_-5.res

        pat = r'H125-LHC13-R04-MSbar-xmur(?P<muR>\d+)-xmuf(?P<muF>\d+)_(?P<kappab>[\-\d]+)_(?P<kappac>[\-\d]+)-xQ(?P<Q>\d+)'

        match = re.search( pat, theoryFilename )

        if not match:
            print '    No match found for theory file \'{0}\''.format( theoryFilename )
            print '    The pattern was:'
            print '    ' + pat
            continue

        muR    = float(match.group('muR')) / 50.
        muF    = float(match.group('muF')) / 50.
        Q      = float(match.group('Q')) / 50.
        kappab = float(match.group('kappab'))
        kappac = float(match.group('kappac'))

        if verbose:
            print '    muR = {0:<6.1f}, muF = {1:<6.1f}, Q = {2:<6.1f}, kappab = {3:<6.1f}, kappac = {4:<6.1f}'.format(
                muR, muF, Q, kappab, kappac
                )

        pts, matched_xss, resummed_xss = ReadLinesOfYukawaTheoryFile( theoryFile, verbose )
        newBinCenters, binBoundaries, binWidths = BinningHeuristic( pts, manualSwitchAt50=False )


        def numberStr(number):
            if number.is_integer():
                return ('{0:d}'.format(int(number))).replace('.','p').replace('-','m')
            else:
                return ('{0:.2f}'.format(number)).replace('.','p').replace('-','m')
        outname = join( outdir, 'muR_{0}_muF_{1}_Q_{2}_kappab_{3}_kappac_{4}.txt'.format(
            numberStr(muR),
            numberStr(muF),
            numberStr(Q),
            numberStr(kappab),
            numberStr(kappac),
            ))
        with open( outname, 'w' ) as outFp:
            w = lambda text: outFp.write( text + '\n' )

            w( 'muR={0}'.format(    muR ) )
            w( 'muF={0}'.format(    muF ) )
            w( 'Q={0}'.format(      Q ) )
            w( 'kappab={0}'.format( kappab ) )
            w( 'kappac={0}'.format( kappac ) )

            w( 'binBoundaries={0}'.format(         ','.join( map(str, binBoundaries ) ) ) )
            w( 'crosssection={0}'.format(          ','.join( map(str, resummed_xss ) ) ) )

            w( 'resummed_crosssection={0}'.format( ','.join( map(str, resummed_xss ) ) ) )
            w( 'matched_crosssection={0}'.format(  ','.join( map(str, matched_xss ) ) ) )

            if smFound:
                # print '[FIXME] Also add ratio to derived theory file'
                ratios = [ xs / SMxs for xs, SMxs in zip( resummed_xss, SM.resummed_xss ) ]
                w( 'ratios={0}'.format( ','.join( map(str, ratios ) ) ) )


                print ratios
                print





def CreateDerivedTheoryFiles(
    pattern = None
    ):

    outdir = 'derivedTheoryFiles_{0}'.format( datestr )
    if not isdir( outdir ): os.makedirs( outdir )


    # ======================================
    # Determine for which files to create the derived files

    theories = []
    for ratio in physicsFileDict.keys():
        # print 'Trying pattern \'{0}\' on \'{1}\''.format( pattern, ratio )
        if re.search( pattern, ratio ):
            theories.append( ratio  )
            # print re.search( pattern, ratio ).group(1)
        # else:
        #     print '    No match'

    if len(theories) == 0:
        Commands.ThrowError( 'Pattern \'{0}\' does not match any know keys'.format(pattern) )
        return


    # ======================================
    # First get the SM predictions; need to read at least one theory file to know the binning

    dummymu, binBoundaries, binWidths, binCenters = ReadTheoryFile( physicsFileDict[theories[0]], applyHeuristic=True )

    SM_NNLO = ReadTheoryFile( physicsFileDict['SM_NNLO'], applyHeuristic=True )[0]
    print '[info] Multiplying list by 1/2.27 (1/(HggBR*1000, from Agnieszka\'s plotting script)'
    SM_NNLO = [ 1/2.27 * i for i in SM_NNLO ]
    print '[fixme] Cutting point >400GeV away for now'
    SM_NNLO = SM_NNLO[:len(binCenters)+1]  # Results go up to 800

    SM_NNLO_upRatio   = ReadTheoryFile( physicsFileDict['SM_NNLO_upRatio'], applyHeuristic=True )[0]
    SM_NNLO_upRatio   = SM_NNLO_upRatio[:len(binCenters)+1]  # Results go up to 800
    SM_NNLO_downRatio = ReadTheoryFile( physicsFileDict['SM_NNLO_downRatio'], applyHeuristic=True )[0]
    SM_NNLO_downRatio = SM_NNLO_downRatio[:len(binCenters)+1]  # Results go up to 800

    SM_NNLO_up   = [ i*j for i,j in zip( SM_NNLO_upRatio, SM_NNLO ) ]
    SM_NNLO_down = [ i*j for i,j in zip( SM_NNLO_downRatio, SM_NNLO ) ]

    # SM_NNLO has the same bins as coupling variation files
    # No rebinning necessary now, but the code below does the trick

    # # Read standard model values in original binning
    # SM_NNLO_originalValues, SM_NNLO_originalBinBoundaries = ReadTheoryFile(
    #     physicsFileDict['SM_NNLO'], applyHeuristic=True )[0:2]
    
    # # Get the rebinned SM values w.r.t. the binning that was used for the scaling of couplings
    # SM_NNLO = Rebin(
    #     SM_NNLO_originalBinBoundaries, SM_NNLO_originalValues,
    #     binBoundaries,
    #     verbose=True
    #     )

    outname = join( outdir, 'SM_NNLO.txt' )
    with open( outname, 'w' ) as outFp:
        w = lambda text: outFp.write( text + '\n' )
        w( 'binBoundaries={0}'.format( ','.join( map(str, binBoundaries ) ) ) )
        w( 'crosssection={0}'.format( ','.join( map(str, SM_NNLO ) ) ) )
        w( 'crosssection_up={0}'.format( ','.join( map(str, SM_NNLO_up ) ) ) )
        w( 'crosssection_down={0}'.format( ','.join( map(str, SM_NNLO_down ) ) ) )


    # ======================================
    # Read files and determine binning

    knownCouplings = [ 'ct', 'cb', 'cg' ]

    for theory in theories:
        mus = ReadTheoryFile( physicsFileDict[theory], applyHeuristic=True )[0]
        crosssection = [ mu * SMXS for mu, SMXS in zip( mus, SM_NNLO ) ]

        couplings = {}
        for coupling in knownCouplings:
            match = re.search( r'{0}_([mp\d]+)'.format(coupling), theory )
            if match:
                couplings[coupling] = float(match.group(1).replace('p','.').replace('m','-'))


        outname = join( outdir, theory + '.txt' )
        with open( outname, 'w' ) as outFp:
            w = lambda text: outFp.write( text + '\n' )

            for coupling in knownCouplings:
                match = re.search( r'{0}_([mp\d]+)'.format(coupling), theory )
                if match:
                    couplingValue = float(match.group(1).replace('p','.').replace('m','-'))
                    w( '{0}={1}'.format( coupling, couplingValue ) )

            w( 'binBoundaries={0}'.format( ','.join( map(str, binBoundaries ) ) ) )
            w( 'crosssection={0}'.format( ','.join( map(str, crosssection ) ) ) )
            w( 'ratios={0}'.format( ','.join( map(str, mus ) ) ) )


def ReadDerivedTheoryFile(
    derivedTheoryFile,
    returnContainer = False,
    ):

    with open( derivedTheoryFile, 'r' ) as theoryFp:
        lines = [ l.strip() for l in theoryFp.readlines() if len(l.strip()) > 0 and not l.strip().startswith('#') ]


    if returnContainer:
        container = Container()

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





def LoadTheoryCurves( pattern=None ):

    # theoryPath = lambda filename: join( 'suppliedInput/numbersFromAgnieszka', filename )

    divlist  = lambda l1, l2: map( operator.truediv, l1, l2 )
    prodlist = lambda l1, l2: map( operator.mul, l1, l2 )

    pt_SM_NNLO, SM_NNLO, ratio_SM_NNLO_down, ratio_SM_NNLO_up = LoadStandardModelCurves(
        ptAxis=ReadTheoryFile( TheoryPath('ratio_ctup_new') )[0] )


    if not pattern:
        doRatios = [
            # 'ct_1p2_cb_m2p98_cg_m0p03',
            # 'ct_1p3_cb_m0p85_cg_m0p03',
            # 'ct_1p4_cb_3p31_cg_m0p03',
            # 'ct_1p2_cb_m4p89_cg_m0p04',
            # 'ct_1p3_cb_m3p34_cg_m0p04',
            # 'ct_1p5_cb_1p88_cg_m0p04',
            # 'ct_1p4_cb_m3p67_cg_m0p05',
            # 'ct_1p5_cb_m1p79_cg_m0p05',
            'ct_0p1_cg_0p075',
            'ct_0p5_cg_0p042',
            'ct_1p5_cg_m0p042',
            'ct_2p0_cg_m0p083',
            ]
    
    elif pattern == 'SM':
        return pt_SM_NNLO, SM_NNLO, ratio_SM_NNLO_down, ratio_SM_NNLO_up

    else:
        doRatios = []
        for ratio in physicsFileDict.keys():
            # print 'Trying pattern \'{0}\' on \'{1}\''.format( pattern, ratio )
            if re.search( pattern, ratio ):
                doRatios.append( ratio )
                # print re.search( pattern, ratio ).group(1)
            # else:
            #     print '    No match'

        if len(doRatios) == 0:
            Commands.ThrowError( 'Pattern \'{0}\' does not match any know keys'.format(pattern) )
            return


    Tgs = []
    for ratioKey in doRatios:

        pt_variation, ratio = ReadTheoryFile( physicsFileDict[ratioKey] )

        Tg = GetTheoryTGraph(
            ratioKey,
            pt_variation, ratio,
            prodlist( ratio_SM_NNLO_down, ratio ),
            prodlist( ratio_SM_NNLO_up, ratio )
            )

        Tg.SM            = SM_NNLO
        Tg.SM_up_ratio   = ratio_SM_NNLO_up
        Tg.SM_down_ratio = ratio_SM_NNLO_down

        Tgs.append( Tg )


    BasicTheoryPlot(Tgs)
    return Tgs



def LoadStandardModelCurves(
    ptAxis = None
    ):

    # if not ptAxis:
    #     # Read 1 file to get the pt axis of the variation
    #     ptAxis, ratio_ct_1p1    = ReadTheoryFile( TheoryPath('ratio_ctup_new') )

    # Read the values of the SM cross sections
    pt_SM_NLO,  SM_NLO             = ReadTheoryFile( TheoryPath('SM_NLO') )
    pt_SM_NLO,  ratio_SM_NLO_down  = ReadTheoryFile( TheoryPath('SMmin_NLO') )
    pt_SM_NLO,  ratio_SM_NLO_up    = ReadTheoryFile( TheoryPath('SMmax_NLO') )
    pt_SM_NNLO, SM_NNLO            = ReadTheoryFile( TheoryPath('SM_NNLO') )
    pt_SM_NNLO, ratio_SM_NNLO_down = ReadTheoryFile( TheoryPath('SMmin_NNLO') )
    pt_SM_NNLO, ratio_SM_NNLO_up   = ReadTheoryFile( TheoryPath('SMmax_NNLO') )

    if not ptAxis == None:
        SM_NNLO            = MapFineToCoarse( pt_SM_NNLO, SM_NNLO, ptAxis )
        ratio_SM_NNLO_down = MapFineToCoarse( pt_SM_NNLO, ratio_SM_NNLO_down, ptAxis )
        ratio_SM_NNLO_up   = MapFineToCoarse( pt_SM_NNLO, ratio_SM_NNLO_up, ptAxis )
        pt_SM_NNLO         = ptAxis

    return pt_SM_NNLO, SM_NNLO, ratio_SM_NNLO_down, ratio_SM_NNLO_up



def ReadTheoryFile(
    theoryFile,
    applyHeuristic=False,
    ):

    with open( theoryFile, 'r' ) as theoryFp:
        lines = [ l.strip() for l in theoryFp.readlines() if not len(l.strip())==0 and not l.strip().startswith('#') ]

    binCenters = []
    mu = []
    for line in lines:
        components = line.split()
        binCenters.append( float(components[0]) )
        mu.append( float(components[1]) )

    if applyHeuristic:
        binCenters, binBoundaries, binWidths = BinningHeuristic( binCenters )
        return mu, binBoundaries, binWidths, binCenters
    else:
        return binCenters, mu


# Theory files contain a somewhat non-consistent binning, turn it into a well-defined binning
def BinningHeuristic(
    binCenters,
    manualSwitchAt50 = True,
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
        muBoundLeft  = [ 0 for i in xrange(nBins) ]
        muBoundRight = [ 0 for i in xrange(nBins) ]

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
    Tg.yMin = min( muBoundLeft )
    Tg.yMax = max( muBoundRight )

    Tg.fourth_xMin = sorted( binBoundaries )[4]
    Tg.fourth_xMax = sorted( binBoundaries, reverse=True )[4]
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


# def MapFineToCoarse(
#     xsFine, ysFine, xsCoarse, ysCoarse=None,
#     ):

#     nFine = len(xsFine)
#     nCoarse = len(xsCoarse)

#     xFineMin = min(xsFine)
#     xFineMax = max(xsFine)
#     yInterpolatedFine = []

#     for iCoarse in xrange(nCoarse):
#         xCoarse = xsCoarse[iCoarse]

#         if xCoarse in xsFine:
#             # Simple case, just copy value
#             yInterpolatedFine.append( ysFine[ xsFine.index(xCoarse) ] )

#         elif xCoarse < xFineMin and xCoarse > xFineMax:
#             print 'Interpolation error'

#         else:

#             # Loop over fine bins, get 2 enclosing boundaries
#             for iFine in xrange(nFine):

#                 # break condition:
#                 if xCoarse < xsFine[iFine]:
#                     leftx  = xsFine[iFine-1]
#                     rightx = xsFine[iFine]
#                     lefty  = ysFine[iFine-1]
#                     righty = ysFine[iFine]

#                     print 'Interpolating for x = {0} between x = {1} and x = {2}'.format( xCoarse, leftx, rightx )

#                     yInterpolatedFine.append(
#                         lefty + 
#                             (xCoarse-leftx) / (rightx-leftx)
#                             *
#                             (righty-lefty)
#                         )

#                     print '    Found y_interpolated = {0} (between y = {1} and y = {2})'.format( yInterpolatedFine[-1], lefty, righty )

#                     break

#     return yInterpolatedFine

def MapFineToCoarse(
        theoryBinBoundaries,
        theoryBinValues,
        expBinBoundaries,
        lastBinIsOverflow = False,
        ):
    theoryIntegralFunction = GetIntegral( theoryBinBoundaries, theoryBinValues )
    expBinValues = []
    for iBinExp in xrange(len(expBinBoundaries)-1):
        expBinValues.append(
            theoryIntegralFunction( expBinBoundaries[iBinExp], expBinBoundaries[iBinExp+1] ) / ( expBinBoundaries[iBinExp+1] - expBinBoundaries[iBinExp] )
            )

    if lastBinIsOverflow:
        expBinValues[-1] = theoryIntegralFunction( expBinBoundaries[-2], theoryBinBoundaries[-1] ) / ( theoryBinBoundaries[-1] - expBinBoundaries[-2] )

    return expBinValues



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






def Rebin(
    ptFine, sigmasFine,
    ptCoarse,
    verbose=False,
    ):

    integralfunction = GetIntegral( ptFine, sigmasFine )
    
    sigmasCoarse = []
    for iBinCoarse in xrange( len(ptCoarse)-1 ):
        integral       = integralfunction( ptCoarse[iBinCoarse], ptCoarse[iBinCoarse+1], verbose=verbose )
        integralPerGeV = integral / ( ptCoarse[iBinCoarse+1] - ptCoarse[iBinCoarse] )
        if verbose: print 'Integral for {0:7.2f} to {1:7.2f}: {2}'.format( ptCoarse[iBinCoarse], ptCoarse[iBinCoarse+1], integralPerGeV )
        sigmasCoarse.append( integralPerGeV )

    return sigmasCoarse



def MapPredictionToExperimental(
    ptTheory, sigmaTheory,
    binning,
    verbose = False,
    makeTGraph = None,
    ):

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
    nPoints = len(scan['ct'])


    keys = [ 'cg', 'ct', 'deltaNLL' ]
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

        iBinX = xBinCenters.index( scan[xCoupling][iPoint] )
        iBinY = yBinCenters.index( scan[yCoupling][iPoint] )

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
    H2.SetMaximum( 50. )
    c.Update()


    Tpoint = ROOT.TGraph( 1, array( 'd', [xBestfit] ), array( 'd', [yBestfit] ) )
    Tpoint.SetMarkerSize(2)
    Tpoint.SetMarkerStyle(34)
    Tpoint.Draw('P')

    Tpoint_SM = ROOT.TGraph( 1, array( 'd', [1.] ), array( 'd', [0.] ) )
    Tpoint_SM.SetMarkerSize(2)
    Tpoint_SM.SetMarkerStyle(21)
    Tpoint_SM.Draw('P')


    SaveC( 'couplingscan2D', asROOT=True )
    SetCMargins()

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



    print '\nSM Couplings:'
    for coupling in couplings:
        print '    {0:5}: {1}'.format( coupling.GetName(), coupling.getVal() )

    print 'SM yieldParameters:'
    for yieldParameter in yieldParameters:
        print '    {0:20}: {1}'.format( yieldParameter.GetName(), yieldParameter.getVal() )



    yPerCoupling = []
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

        yPerCoupling.append( ( y, yParametrization ) )

    return yPerCoupling





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