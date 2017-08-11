#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, tempfile, shutil, re, glob, itertools, sys, numpy, operator, pprint, re, random
from os.path import *
from operator import itemgetter
from array import array
from math import log, exp, sqrt, copysign
from copy import deepcopy
from glob import glob

from time import strftime
datestr = strftime( '%b%d' )

from Container import Container
from OutputInterface import OutputContainer
import Commands
import TheoryCommands
from Parametrization import Parametrization


########################################
# General functions
########################################

#____________________________________________________________________
def Div( l1, l2 ):
    ret = [ i / j if j!=0. else 0. for i, j in zip( l1, l2 ) ]
    return ret


#____________________________________________________________________
FILEFINDERDIR = '.'
def SetFileFinderDir( directory ):
    global FILEFINDERDIR
    FILEFINDERDIR = abspath( directory )


#____________________________________________________________________
def FileFinder( **kwargs ):

    # ======================================
    # Some options
    
    if 'verbose' in kwargs:
        verbose = kwargs['verbose']
        del kwargs['verbose']
    else:
        verbose = False

    if 'directory' in kwargs:
        allFilesDir = kwargs['directory']
        del kwargs['directory']
    else:
        allFilesDir = FILEFINDERDIR

    if 'expectOneFile' in kwargs:
        expectOne = kwargs['expectOneFile']
        del kwargs['expectOneFile']
    else:
        expectOne = False

    if 'filter' in kwargs:
        if isinstance( kwargs['filter'], basestring ):
            filefilters = [ kwargs['filter'] ]
        else:
            filefilters = kwargs['filter']
        del kwargs['filter']
    else:
        filefilters = []


    # ======================================
    # Find files

    allFiles = [ f for f in glob( join( allFilesDir, '*' ) ) if isfile(f) ]

    acceptedFiles = []
    for theoryFile in allFiles:

        passedFilter = True
        for filefilter in filefilters:
            if filefilter in theoryFile: passedFilter = False
        if not passedFilter: continue

        acceptancePerKey = []
        for key, value in kwargs.iteritems():
            if value == '*':
                if key in theoryFile:
                    acceptedByThisKey = True
                else:
                    acceptedByThisKey = False
            elif not(
                    re.search( r'{0}_{1}\D'.format( key, value ), theoryFile )
                    or re.search( r'{0}_{1}\D'.format( key, Commands.ConvertFloatToStr(value) ), theoryFile )
                    ):
                acceptedByThisKey = False
            else:
                acceptedByThisKey = True
            acceptancePerKey.append(acceptedByThisKey)

        if all(acceptancePerKey):
            acceptedFiles.append( theoryFile )


    if len(acceptedFiles) == 0:
        Commands.ThrowError(
            'The following call to FileFinder failed to retrieve any files:\n      '
            + '\n      '.join( [ '{0} = {1}'.format(key,value) for key, value in kwargs.iteritems() ] )
            )
        sys.exit()

    if expectOne and len(acceptedFiles) > 1:
        Commands.ThrowError(
            'Found more than 1 file for the following keys:'
            + '\n    ' + '\n    '.join( [ '{0} = {1}'.format(key,value) for key, value in kwargs.iteritems() ] )
            + '\n  File list:'
            + '\n    ' + '\n    '.join([ relpath( f, '.' ) for f in acceptedFiles ])
            )
        sys.exit()

    if verbose:
        print '[info] FileFinder: Using the following keys:'
        print '    ' + '\n    '.join( [ '{0} = {1}'.format(key,value) for key, value in kwargs.iteritems() ] )
        print '  Found the following file(s):'
        print '    ' + '\n    '.join([ relpath( f, '.' ) for f in acceptedFiles ])

    if expectOne:
        return acceptedFiles[0]
    else:
        return acceptedFiles


#____________________________________________________________________
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

        # Assume central scales if not specified
        if not hasattr( container, 'muR' ):
            container.muR = 1.0
        if not hasattr( container, 'muF' ):
            container.muF = 1.0
        if not hasattr( container, 'Q' ):
            container.Q = 1.0

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


def ReadDerivedTheoryFileToTGraph(
        derivedTheoryFile,
        name=None,
        yAttr=None,
        ):

    container = ReadDerivedTheoryFile( derivedTheoryFile )
    outputcontainer = OutputContainer( container )
    outputcontainer.GetTGraph( name=name, yAttr=yAttr )
    return outputcontainer.Tg


########################################
# Interface fromPier
########################################

#____________________________________________________________________
def DumpContainerToFile( container, prefix, outdir ):

    outname = '{0}_kappab_{1}_kappac_{2}'.format(
        prefix,
        Commands.ConvertFloatToStr( container.kappab ),
        Commands.ConvertFloatToStr( container.kappac )
        )
    
    if hasattr( container, 'muR' ) and hasattr( container, 'muF' ):
        outname += '_muR_{0}_muF_{1}'.format(
            Commands.ConvertFloatToStr( container.muR ),
            Commands.ConvertFloatToStr( container.muF )
            )

    if hasattr( container, 'Q' ):
        outname += '_Q_{0}'.format( Commands.ConvertFloatToStr(container.Q) )

    outname += '.txt'

    container.derivedTheoryFilePath = abspath( join( outdir, outname ) )

    with open( container.derivedTheoryFilePath, 'w' ) as outFp:
        w = lambda text: outFp.write( text + ' \n' )

        w( 'file={0}'.format(container.file) )
        w( 'kappab={0}'.format(container.kappab) )
        w( 'kappac={0}'.format(container.kappac) )

        if hasattr( container, 'muR' ) and hasattr( container, 'muF' ):
            w( 'muR={0}'.format(container.muR) )
            w( 'muF={0}'.format(container.muF) )

        if hasattr( container, 'Q' ):
            w( 'Q={0}'.format(container.Q) )

        w( 'binBoundaries={0}'.format( ','.join(map( str, container.binBoundaries )) ) )
        w( 'crosssection={0}'.format( ','.join(map( str, container.crosssection )) ) )
        if hasattr( container, 'binCenters' ): w( 'binCenters={0}'.format( ','.join(map( str, container.binCenters )) ) )
        w( 'ratios={0}'.format( ','.join(map( str, container.ratios )) ) )


#____________________________________________________________________
def addContainer( self, c2 ):
    nBinsC2 = len(c2.binBoundaries)-1
    for iBin in xrange(nBinsC2):
        self.crosssection[iBin] += c2.crosssection[iBin]
Container.addContainer = addContainer


#____________________________________________________________________
def MergeGluonAndQuarkInduced(
        gI_theoryFileDir,
        qI_theoryFileDir,
        verbose = True,
        ):

    # ======================================
    # IO

    outdir = 'derivedTheoryFiles_{0}_YukawaSummed'.format( datestr )
    if not isdir( outdir ): os.makedirs( outdir )

    qI_theoryFiles = FileFinder( directory=qI_theoryFileDir, filter='muR' )
    gI_theoryFiles = FileFinder( directory=gI_theoryFileDir, muR=1.0, muF=1.0, Q=1.0 )
    
    qI_containers = [ ReadDerivedTheoryFile(tF) for tF in qI_theoryFiles ]
    gI_containers = [ ReadDerivedTheoryFile(tF) for tF in gI_theoryFiles ]

    # Pointers to SM containers
    qI_SM = [ c for c in qI_containers if c.kappab == 1.0 and c.kappac == 1.0 ][0]
    gI_SM = [ c for c in gI_containers if c.kappab == 1.0 and c.kappac == 1.0 ][0]


    # ======================================
    # Rebin qI containers so the boundaries match with gI

    # gI has much more bins - only rebin up to where qI is defined
    lastBinBoundary = qI_SM.binBoundaries[-1]
    if not lastBinBoundary in gI_SM.binBoundaries:
        Commands.ThrowError( 'qI bin boundary {0} not found in gI; need to implement something more robust'.format( lastBinBoundary ) )
        sys.exit()
    newBinBoundaries = gI_SM.binBoundaries[: gI_SM.binBoundaries.index(lastBinBoundary)+1 ]

    for container in qI_containers:
        TheoryCommands.RebinDerivedTheoryContainer( container, newBinBoundaries )
    TheoryCommands.RebinDerivedTheoryContainer( qI_SM, newBinBoundaries )


    # ======================================
    # Create parametrizations for qI
    # (needed because mass effects need to be scaled)

    parametrization = Parametrization()

    # parametrization.Parametrize( qI_containers )
    parametrization.ParametrizeByFitting( qI_containers )

    mb_old = 4.65
    mb_new = 2.963
    mc_old = 1.275
    mc_new = 0.655


    qI_SM.crosssection = parametrization.Evaluate(
        kappab = qI_SM.kappab * ( mb_new/mb_old ),
        kappac = qI_SM.kappac * ( mc_new/mc_old )
        )

    # Calculate the effect of the different masses
    for container in qI_containers:
        container.crosssection = parametrization.Evaluate(
            kappab = container.kappab * ( mb_new/mb_old ),
            kappac = container.kappac * ( mc_new/mc_old )
            )
        container.ratios = [ xs / smxs if not smxs == 0. else 0. for xs, smxs in zip( container.crosssection, qI_SM.crosssection ) ]


    # Sum up; start with gluon induced
    summed_containers = deepcopy( gI_containers )
    summed_SM = deepcopy( gI_SM )
    DumpContainerToFile( summed_SM, prefix='YukawaSummed', outdir=outdir )

    for container in summed_containers:

        # Find the matching qI container
        for qI in qI_containers:

            if (
                    container.kappab == qI.kappab
                    and container.kappac == qI.kappac
                    ):

                container.addContainer( qI )
                container.ratios = [ xs / smxs for xs, smxs in zip( container.crosssection, summed_SM.crosssection ) ]
                DumpContainerToFile( container, prefix='YukawaSummed', outdir=outdir )
                break

        else:
            print '[warning] Could not find a matching container for kappab = {0}, kappac = {1}'.format( container.kappab, container.kappac )


#____________________________________________________________________
def ReadLinesOfTheoryFile_YukawaQuarkInduced( theoryFile, verbose=False, SM=None ):

    # Read lines
    with open( theoryFile, 'r' ) as theoryFp:
        lines = [ line.strip() for line in theoryFp.readlines() ]
    commentlines = [ line.strip() for line in lines if line.startswith('#') ]
    lines        = [ line for line in lines if not line.startswith('#') and len(line) > 0 ]

    pts          = []
    xss          = []

    for line in lines:
        components = line.split()
        pt         = float(components[0])

        # Factor 1000 to scale [nb] to [pb]
        xs         = 1000. * float(components[1])

        if verbose:
            if SM is None:
                print '    pt = {0:<8.3f} |  xs = {1:<10.6f}'.format(
                    pt, xs
                    )
            else:
                i = SM.pt.index(pt)
                print '    pt = {0:<8.3f} | xs = {1:<10.6f} | SM xs = {2:<10.6f} | ratio = {3:<10.6f}'.format(
                    pt, xs, SM.xs[i], xs / SM.xs[i] if SM.xs[i] != 0. else 0.
                    )


        pts.append( pt )
        xss.append( xs )

    return pts, xss


#____________________________________________________________________
def CreateDerivedTheoryFiles_YukawaQuarkInduced(
        theoryDir = 'suppliedInput/fromPier/13tev-pth_quarkInduced_Aug04/',
        verbose = True,
        ):

    outdir = 'derivedTheoryFiles_{0}_YukawaQuarkInduced'.format( datestr )
    if not isdir( outdir ): os.makedirs( outdir )

    theoryDir = abspath( theoryDir )
    theoryFiles = glob( join( theoryDir, '*.pth' ) )

    # ======================================
    # Try to find a 'SM' file first so it's possible to calculate ratios

    SM = Container()
    SM.file = join( theoryDir, 'higgs_plus_jet_13tev_1_1_mur050_muf050.pth' )

    if verbose:
        print '\n\n' + '='*80
        print '='*80
        print 'Processing SM \'{0}\'\n'.format( SM.file )

    SM.pt, SM.xs = ReadLinesOfTheoryFile_YukawaQuarkInduced( SM.file, verbose )
    SM.kappab = 1
    SM.kappac = 1

    if verbose:
        print '\n\nFound SM.xs:'
        print SM.xs

        print '\n\nContents of \'{0}\''.format( SM.file )
        with open( SM.file) as fp:
            print fp.read()


    # ======================================
    # Rest of theory files

    pat = r'higgs_plus_jet_13tev_([0-9mp]+)_([0-9mp]+)'
    scaleVarPat = r'mur([0-9]+)_muf([0-9]+)'

    for theoryFile in theoryFiles:
        container = Container()
        container.file = theoryFile

        if verbose:
            print '\n\n' + '='*80
            print '='*80
            print 'Processing \'{0}\'\n'.format( container.file )

        match = re.search( pat, container.file )
        if not match:
            print '[warning] Could not deduce couplings from file \'{0}\''.format( container.file )
        else:
            container.kappab = Commands.ConvertStrToFloat( match.group(1) )
            container.kappac = Commands.ConvertStrToFloat( match.group(2) )

        isScaleVariation = False
        match = re.search( scaleVarPat, container.file )
        if match:
            isScaleVariation = True
            container.muR = float( match.group(1) ) / 50
            container.muF = float( match.group(2) ) / 50

        container.pt, container.xs = ReadLinesOfTheoryFile_YukawaQuarkInduced( container.file, verbose, SM=SM )
        container.ratios = [ xs / smxs if not smxs == 0. else 0. for xs, smxs in zip( container.xs, SM.xs ) ]
        container.crosssection = container.xs

        container.binCenters, container.binBoundaries, container.binWidths = TheoryCommands.BinningHeuristic(
            container.pt, manualSwitchAt50 = False, manualSwitchAt5 = True )


        # ======================================
        # Dump container to file

        DumpContainerToFile( container, prefix='YukawaQuarkInduced', outdir=outdir )

        if verbose:
            print '\n\nContents of \'{0}\':\n'.format( relpath( container.file, '.' ) )
            with open( container.file) as fp:
                print fp.read()

            print '\n\nContents of \'{0}\':\n'.format( relpath( container.derivedTheoryFilePath, '.' ) )
            with open(container.derivedTheoryFilePath) as fp:
                print fp.read()


#____________________________________________________________________
def ReadLinesOfTheoryFile_YukawaGluonInduced( theoryFile, verbose=False ):
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
        matched_xs  = float(components[1]) * 1000.
        resummed_xs = float(components[2]) * 1000.

        if verbose:
            print '    pt = {0:<8.3f} |  matched_xs = {1:<10.6f} |  resummed_xs = {2:<10.6f}'.format(
                pt, matched_xs, resummed_xs
                )

        pts.append( pt )
        matched_xss.append( matched_xs )
        resummed_xss.append( resummed_xs )

    return pts, matched_xss, resummed_xss


#____________________________________________________________________
def CreateDerivedTheoryFiles_YukawaGluonInduced(
        theoryDir = 'suppliedInput/fromPier/histograms_ggH_May17/',
        mainCrossSection = 'matched',
        verbose = True,
        ):

    outdir = 'derivedTheoryFiles_{0}_YukawaGluonInduced'.format( datestr )
    if not isdir( outdir ): os.makedirs( outdir )

    theoryDir = abspath( theoryDir )
    theoryFiles = glob( join( theoryDir, '*.res' ) )


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
            pts, matched_xss, resummed_xss = ReadLinesOfTheoryFile_YukawaGluonInduced( theoryFile, verbose )
            newBinCenters, binBoundaries, binWidths = TheoryCommands.BinningHeuristic( pts, manualSwitchAt50=False, manualSwitchAt5=True )

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

        pts, matched_xss, resummed_xss = ReadLinesOfTheoryFile_YukawaGluonInduced( theoryFile, verbose )
        newBinCenters, binBoundaries, binWidths = TheoryCommands.BinningHeuristic( pts, manualSwitchAt50=False, manualSwitchAt5=True )


        container = Container()

        container.file = theoryFile

        if mainCrossSection == 'matched':
            container.crosssection = matched_xss
            container.ratios = Div( container.crosssection, SM.matched_xss )
        elif mainCrossSection == 'resummed':
            container.crosssection = resummed_xss
            container.ratios = Div( container.crosssection, SM.resummed_xss )

        container.binBoundaries = binBoundaries

        container.muR    = muR
        container.muF    = muF
        container.Q      = Q
        container.kappab = kappab
        container.kappac = kappac

        DumpContainerToFile( container, prefix='YukawaGluonInduced', outdir=outdir )



########################################
# Interface fromAgnieszka
########################################

AgnieszkasFilenameDecoder = {
    'ratio_ctup_new'          : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.1 },
    'ratio_ctdown_new'        : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 0.9 },
    'ratio_cgup_new'          : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cg' : 0.008 },
    'ratio_cgdown_new'        : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cg' : -0.008 },
    'ratio_cbup_new'          : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cb' : 4.0 },
    'ratio_cbdown_new'        : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cb' : -2.0 },
    'ratio_cg003ct12_new'     : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.2, 'cb' : -2.98 ,  'cg' : -0.03 },
    'ratio_cg003ct13sw_new'   : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.3, 'cb' : -0.85 ,  'cg' : -0.03 },
    'ratio_cg003ct14sw_new'   : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.4, 'cb' : 3.31 ,   'cg' : -0.03 },
    'ratio_cg004ct12sw_new'   : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.2, 'cb' : -4.89 ,  'cg' : -0.04 },
    'ratio_cg004ct13sw_new'   : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.3, 'cb' : -3.34 ,  'cg' : -0.04 },
    'ratio_cg004ct15sw_new'   : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.5, 'cb' : 1.88 ,   'cg' : -0.04 },
    'ratio_cg005ct14sw_new'   : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.4, 'cb' : -3.67 ,  'cg' : -0.05 },
    'ratio_cg005ct15sw_new'   : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.5, 'cb' : -1.79 ,  'cg' : -0.05 },
    'ratio_ctcb05_new'        : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 0.5, 'cb' : -7.46 },
    'ratio_ctcb08_new'        : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 0.8, 'cb' : -3.67 },
    'ratio_ctcb09_new'        : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 0.9, 'cb' : -1.79 },
    'ratio_ctcb11_new'        : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.1, 'cb' : 3.79 },
    'ratio_ctcb12_new'        : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.2, 'cb' : 4.67 },
    'ratio_ctcg01_new'        : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 0.1, 'cg' : 0.075 },
    'ratio_ctcg05_new'        : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 0.5, 'cg' : 0.042 },
    'ratio_ctcg15_new'        : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.5, 'cg' : -0.042 },
    'ratio_ctcg2_new'         : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 2.0, 'cg' : -0.083 },
    # 'SM_NLO'                  : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 1.0 },
    # 'SMmin_NLO'               : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 1.0 },
    # 'SMmax_NLO'               : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 1.0 },
    'SM_NNLO'                 : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 1.0 },
    # 'SMmin_NNLO'              : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 1.0 },
    # 'SMmax_NNLO'              : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 1.0 },
    'HRes_mR1mF1.top'         : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 1.0 },
    'HRes_mR1mF1Q2.top'       : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 2.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 1.0 },
    'HRes_mR1mF1Qh.top'       : { 'muR' : 1.0, 'muF' : 1.0, 'Q' : 0.5, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 1.0 },
    'HRes_mR1mF2.top'         : { 'muR' : 1.0, 'muF' : 2.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 1.0 },
    'HRes_mR1mFh.top'         : { 'muR' : 1.0, 'muF' : 0.5, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 1.0 },
    'HRes_mR2mF1.top'         : { 'muR' : 2.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 1.0 },
    'HRes_mR2mF2.top'         : { 'muR' : 2.0, 'muF' : 2.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 1.0 },
    'HRes_mRhmF1.top'         : { 'muR' : 0.5, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 1.0 },
    'HRes_mRhmFh.top'         : { 'muR' : 0.5, 'muF' : 0.5, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 1.0 },
    }


#____________________________________________________________________
def DumpContainerToFile_Top( container, prefix, outdir ):

    outname = prefix
    for coupling in [ 'ct', 'cg', 'cb' ]:
        if hasattr( container, coupling ):
            outname += '_{0}_{1}'.format(
                coupling, Commands.ConvertFloatToStr( getattr(container, coupling) )
                )
    
    if hasattr( container, 'muR' ) and hasattr( container, 'muF' ):
        outname += '_muR_{0}_muF_{1}'.format(
            Commands.ConvertFloatToStr( container.muR ),
            Commands.ConvertFloatToStr( container.muF )
            )

    if hasattr( container, 'Q' ):
        outname += '_Q_{0}'.format( Commands.ConvertFloatToStr(container.Q) )

    outname += '.txt'

    container.derivedTheoryFilePath = abspath( join( outdir, outname ) )

    with open( container.derivedTheoryFilePath, 'w' ) as outFp:
        w = lambda text: outFp.write( text + ' \n' )

        w( 'file={0}'.format(container.file) )

        if hasattr( container, 'ct' ): w( 'ct={0}'.format(container.ct) )
        if hasattr( container, 'cg' ): w( 'cg={0}'.format(container.cg) )
        if hasattr( container, 'cb' ): w( 'cb={0}'.format(container.cb) )

        if hasattr( container, 'muR' ) and hasattr( container, 'muF' ):
            w( 'muR={0}'.format(container.muR) )
            w( 'muF={0}'.format(container.muF) )

        if hasattr( container, 'Q' ):  w( 'Q={0}'.format(container.Q) )

        w( 'binBoundaries={0}'.format( ','.join(map( str, container.binBoundaries )) ) )
        w( 'crosssection={0}'.format( ','.join(map( str, container.crosssection )) ) )
        if hasattr( container, 'binCenters' ): w( 'binCenters={0}'.format( ','.join(map( str, container.binCenters )) ) )
        w( 'ratios={0}'.format( ','.join(map( str, container.ratios )) ) )


#____________________________________________________________________
def ReadLinesOfTheoryFile_Top( theoryFile, verbose=False, SM=None, isSMFile=False ):
    # Read lines
    with open( theoryFile, 'r' ) as theoryFp:
        lines = [ line.strip() for line in theoryFp.readlines() ]
    commentlines = [ line.strip() for line in lines if line.startswith('#') ]
    lines        = [ line for line in lines if not line.startswith('#') and len(line) > 0 ]

    if isSMFile:

        pts  = []
        xss  = []

        for line in lines:
            components = line.split()
            pt         = float(components[0])
            xs         = float(components[1]) / 2.27 # 1 over Hgg BR (peculiarity from Agnieszka)

            if verbose:
                if SM is None:
                    print '    pt = {0:<8.3f} |  xs = {1:<10.6f}'.format( pt, xs )
                else:
                    SMxs = SM.crosssection[ SM.pts.index(pt) ]
                    print '    pt = {0:<8.3f} |  xs = {1:<10.6f} |  SMxs = {2:<10.6f} |  ratio = {3:<10.6f}'.format(
                        pt, xs, SMxs, xs/SMxs if not SMxs == 0. else 0. )

            pts.append( pt )
            xss.append( xs )

        return pts, xss

    else:

        if SM is None:
            Commands.ThrowError( 'SM should be specified when reading a non-SM file' )
            sys.exit()

        pts    = []
        ratios = []
        xss    = []

        for line in lines:

            components = line.split()
            pt         = float(components[0])
            ratio      = float(components[1])

            SMxs = SM.crosssection[ SM.pts.index(pt) ]
            xs   = SMxs * ratio

            if verbose:
                print '    pt = {0:<8.3f} |  xs = {1:<10.6f} |  SMxs = {2:<10.6f} |  ratio = {3:<10.6f}'.format(
                    pt, xs, SMxs, ratio )

            pts.append( pt )
            ratios.append( ratio )
            xss.append( xs )

        return pts, xss, ratios

    


#____________________________________________________________________
def CreateDerivedTheoryFiles_Top(
        theoryDirs = [
            'suppliedInput/fromAgnieszka/HRes_SMEFT_May16',
            'suppliedInput/fromAgnieszka/SMEFTscaling_May16',
            'suppliedInput/fromAgnieszka/ScaleVarNNLO_Jul17',
            ],
        verbose = False,
        ):

    outdir = 'derivedTheoryFiles_{0}_Top'.format( datestr )
    if not isdir( outdir ): os.makedirs( outdir )

    theoryDirs = [ abspath( theoryDir ) for theoryDir in theoryDirs ]
    theoryFiles = []
    for theoryDir in theoryDirs:
        theoryFiles.extend( glob( join( theoryDir, '*' ) ) )


    # ======================================
    # First look for SM file

    for theoryFile in theoryFiles:
        if basename(theoryFile).replace('/','') == 'SM_NNLO':
            smFile = theoryFile
            break
    else:
        Commands.ThrowError( 'Could not find a SM file' )
        sys.exit()

    SM = Container()
    SM.file = smFile

    binCenters, crosssection = ReadLinesOfTheoryFile_Top( SM.file, isSMFile=True, verbose=verbose )
    newBinCenters, binBoundaries, binWidths = TheoryCommands.BinningHeuristic(
        binCenters, manualSwitchAt50=True, manualSwitchAt5=False )

    SM.pts           = binCenters
    SM.binCenters    = deepcopy( newBinCenters )
    SM.binBoundaries = deepcopy( binBoundaries )
    SM.crosssection  = deepcopy( crosssection )
    SM.ratios = [ 1.0 for xs in SM.crosssection ]

    for key, value in AgnieszkasFilenameDecoder[basename(SM.file).replace('/','')].iteritems():
        setattr( SM, key, value )


    # ======================================
    # Process the other theory files

    for theoryFile in theoryFiles:

        if verbose:
            print '\n\n' + '-'*80
            print 'Processing {0}'.format( theoryFile )
            print '\n'

        container = Container()
        container.file = theoryFile

        shortFileName = basename(container.file).replace('/','')
        if not shortFileName in AgnieszkasFilenameDecoder:
            print 'Skipping \'{0}\''.format(shortFileName)
            continue

        binCenters, crosssection, ratios = ReadLinesOfTheoryFile_Top( container.file, verbose, SM=SM )
        newBinCenters, binBoundaries, binWidths = TheoryCommands.BinningHeuristic(
            binCenters, manualSwitchAt50=True, manualSwitchAt5=False )

        container.binCenters    = newBinCenters
        container.binBoundaries = binBoundaries
        container.crosssection  = crosssection
        container.ratios        = ratios

        if verbose:
            print '\nWriting the following to a file:'
            for i in xrange(len(container.crosssection)):
                print '    {0:6.1f} - {1:6.1f}  |  xs = {2:10.6f}  |  ratio = {3:10.6f}'.format(
                    container.binBoundaries[i], container.binBoundaries[i+1],
                    container.crosssection[i], container.ratios[i]
                    )

        for key, value in AgnieszkasFilenameDecoder[shortFileName].iteritems():
            setattr( container, key, value )

        DumpContainerToFile_Top( container, prefix='Top', outdir=outdir )


########################################
# Older
########################################

#____________________________________________________________________
def CreateDerivedTheoryFiles_Agnieszka_OLD(
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



########################################
# Execute this if called directly
########################################

if __name__ == '__main__':
    print 'Import functions from this file'
    print 'Available functions:'
    functions = [ name for (name, thing) in locals().items() if callable(thing) ]
    print '\n'.join( functions )