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

    loadImmediately = False
    if 'loadImmediately' in kwargs:
        loadImmediately = kwargs['loadImmediately']
        del kwargs['loadImmediately']


    # ======================================
    # Find files

    allFiles = [ f for f in glob( join( allFilesDir, '*' ) ) if isfile(f) ]

    acceptedFiles = []
    for theoryFile in allFiles:
        fullPathTheoryFile = theoryFile
        theoryFile = basename(theoryFile)

        if verbose: print '\nChecking acceptance of file \'{0}\''.format(theoryFile)

        passedFilter = True
        for filefilter in filefilters:
            if filefilter in theoryFile: passedFilter = False
        if not passedFilter:
            if verbose: print '  Not accepted by filters {0}'.format(filefilters)
            continue

        acceptancePerKey = []
        for key, value in kwargs.iteritems():

            if not key in theoryFile:
                acceptedByThisKey = False
            elif value == '*':
                acceptedByThisKey = True
            else:
                valueInFile = Commands.ConvertStrToFloat(
                    re.search( r'{0}_([\dpm]+)'.format( key), theoryFile ).group(1)
                    )
                if value == valueInFile or Commands.ConvertFloatToStr(value) == valueInFile:
                    acceptedByThisKey = True
                else:
                    acceptedByThisKey = False
            # elif not(
            #         re.search( r'{0}_{1}\D'.format( key, value ), theoryFile )
            #         or re.search( r'{0}_{1}\D'.format( key, Commands.ConvertFloatToStr(value) ), theoryFile )
            #         ):
            #     acceptedByThisKey = False
            # else:
            #     acceptedByThisKey = True

            if verbose:
                print '  Key = {0:10}, Value = {1:10}, Accepted = {2}'.format( key, value, acceptedByThisKey )

            acceptancePerKey.append(acceptedByThisKey)

        if all(acceptancePerKey):
            acceptedFiles.append( fullPathTheoryFile )
            if verbose:
                print '  File was accepted'


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
        print '\n[info] FileFinder: Using the following keys:'
        print '    ' + '\n    '.join( [ '{0} = {1}'.format(key,value) for key, value in kwargs.iteritems() ] )
        print '  Found the following file(s):'
        print '    ' + '\n    '.join([ relpath( f, '.' ) for f in acceptedFiles ])


    if expectOne:
        if loadImmediately:
            return ReadDerivedTheoryFile( acceptedFiles[0] )
        else:
            return acceptedFiles[0]
    else:
        if loadImmediately:
            return [ ReadDerivedTheoryFile(f) for f in acceptedFiles ]
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


def ReadDerivedTheoryContainerToTGraph(
        derivedTheoryContainer,
        name=None,
        yAttr=None,
        ):
    outputcontainer = OutputContainer( derivedTheoryContainer )
    outputcontainer.GetTGraph( name=name, yAttr=yAttr )
    return outputcontainer.Tg


########################################
# General operations
########################################

def NormalizeToSMCrossSection(
        container,
        SMXS = 55.70628722,
        extraCheck = False,
        ):

    integralFn = TheoryCommands.GetIntegral( container.binBoundaries, container.crosssection )

    rescale = SMXS / integralFn( 0., container.binBoundaries[-1]*10 )

    container.crosssection = [ rescale * xs for xs in container.crosssection ]
    container.ratios       = [ rescale * r  for r  in container.ratios ]

    if extraCheck:
        print '\nIntegral before scaling:'
        print integralFn( 0., container.binBoundaries[-1]*10 )
        print 'Integral after scaling (should be {0}):'.format( SMXS )
        newintegralFn = TheoryCommands.GetIntegral( container.binBoundaries, container.crosssection )
        print newintegralFn( 0., container.binBoundaries[-1]*10 )


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
def ScaleQuarkInduced(
        qI_theoryFileDir,
        verbose = True,
        ):

    # ======================================
    # Create the scaled non-scale-variations files

    # Get non-scale-variation containers
    qI_theoryFiles = FileFinder( directory=qI_theoryFileDir, filter='muR' )
    qI_containers = [ ReadDerivedTheoryFile(tF) for tF in qI_theoryFiles ]

    # Pointer to SM containers
    qI_SM = [ c for c in qI_containers if c.kappab == 1.0 and c.kappac == 1.0 ][0]

    # Create parametrizations for qI
    # (needed because mass effects need to be scaled)

    parametrization = Parametrization()


    parametrization.SetSM( qI_SM )
    parametrization.ParametrizeByFitting( qI_containers, fitWithScipy=True )

    # parametrization.ParametrizeByFitting( qI_containers )


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


    outdir = 'derivedTheoryFiles_{0}_YukawaQuarkInducedScaled'.format( datestr )
    if not isdir( outdir ): os.makedirs( outdir )
    for container in qI_containers:
        DumpContainerToFile( container, prefix='YukawaQuarkInducedScaled', outdir=outdir )


    # ======================================
    # Now do the scale variations
    # Assume simply that the scale variations have the same proportionality w.r.t.
    # the SM *after* the quark mass scaling as before

    scaleVar_theoryFiles = [
        FileFinder( directory=qI_theoryFileDir, kappab=1.0, kappac=1.0, muR=0.5, muF=0.5, expectOneFile=True ),
        FileFinder( directory=qI_theoryFileDir, kappab=1.0, kappac=1.0, muR=0.5, muF=1.0, expectOneFile=True ),
        FileFinder( directory=qI_theoryFileDir, kappab=1.0, kappac=1.0, muR=1.0, muF=0.5, expectOneFile=True ),
        FileFinder( directory=qI_theoryFileDir, kappab=1.0, kappac=1.0, muR=1.0, muF=1.0, expectOneFile=True ),
        FileFinder( directory=qI_theoryFileDir, kappab=1.0, kappac=1.0, muR=2.0, muF=1.0, expectOneFile=True ),
        FileFinder( directory=qI_theoryFileDir, kappab=1.0, kappac=1.0, muR=1.0, muF=2.0, expectOneFile=True ),
        FileFinder( directory=qI_theoryFileDir, kappab=1.0, kappac=1.0, muR=2.0, muF=2.0, expectOneFile=True ),
        ]
    scaleVar_containers = [ ReadDerivedTheoryFile(tF) for tF in scaleVar_theoryFiles ]

    for scaleVar_container in scaleVar_containers:
        # Take the scaled SMXS, and multiply it by the same ratio as before the mass scaling
        scaleVar_container.crosssection = [
            ratio * SMXS for ratio, SMXS in zip( scaleVar_container.ratios, qI_SM.crosssection )
            ]
        DumpContainerToFile( scaleVar_container, prefix='YukawaQuarkInducedScaled', outdir=outdir )


def SumGluonAndQuarkFiles(
        gI_theoryFileDir,
        qI_theoryFileDir,
        verbose = True,
        ):

    # ======================================
    # IO

    qI_theoryFiles = FileFinder( directory=qI_theoryFileDir )
    gI_theoryFiles = FileFinder( directory=gI_theoryFileDir )
    
    qI_containers = [ ReadDerivedTheoryFile(tF) for tF in qI_theoryFiles ]
    gI_containers = [ ReadDerivedTheoryFile(tF) for tF in gI_theoryFiles ]

    qI_SMFile     = FileFinder( directory=qI_theoryFileDir, kappab=1.0, kappac=1.0, Q=1.0, muF=1.0, muR=1.0, expectOneFile=True )
    qI_SM         = ReadDerivedTheoryFile( qI_SMFile )
    gI_SMFile     = FileFinder( directory=gI_theoryFileDir, kappab=1.0, kappac=1.0, Q=1.0, muF=1.0, muR=1.0, expectOneFile=True )
    gI_SM         = ReadDerivedTheoryFile( gI_SMFile )
    TheoryCommands.RebinDerivedTheoryContainer( gI_SM, qI_SM.binBoundaries )

    outdir = 'derivedTheoryFiles_{0}_YukawaSummed'.format( datestr )
    if not isdir( outdir ): os.makedirs( outdir )

    summed_SM = deepcopy( qI_SM )
    summed_SM.crosssection = [ xs_qI + xs_gI for xs_qI, xs_gI in zip( qI_SM.crosssection, gI_SM.crosssection ) ]
    DumpContainerToFile( summed_SM, prefix='YukawaSummed', outdir=outdir )

    for qI in qI_containers:

        # Find matching gI container
        for gI in gI_containers:

            if ( True
                and qI.kappac == gI.kappac
                and qI.kappab == gI.kappab
                and qI.muR    == gI.muR
                and qI.muF    == gI.muF
                and qI.Q      == gI.Q
                ):

                break
        else:
            print 'Could not find a matching gluon-induced container for:'
            print '    kappac = ', qI.kappac
            print '    kappab = ', qI.kappab
            print '    muR    = ', qI.muR
            print '    muF    = ', qI.muF
            print '    Q      = ', qI.Q
            print '    Continuing'
            continue

        TheoryCommands.RebinDerivedTheoryContainer( gI, qI_SM.binBoundaries )

        summed = deepcopy( qI )
        summed.crosssection = [ xs_qI + xs_gI for xs_qI, xs_gI in zip( qI.crosssection, gI.crosssection ) ]
        summed.ratios = [ xs / SMxs for xs, SMxs in zip( summed.crosssection, summed_SM.crosssection ) ]

        DumpContainerToFile( summed, prefix='YukawaSummed', outdir=outdir )



#____________________________________________________________________
def MergeGluonAndQuarkInduced(
        gI_theoryFileDir,
        qI_theoryFileDir,
        verbose = True,
        ):

    # ======================================
    # IO

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


    outdir = 'derivedTheoryFiles_{0}_YukawaSummed'.format( datestr )
    if not isdir( outdir ): os.makedirs( outdir )

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
    'ratio_ctup_new'          : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.1 },
    'ratio_ctdown_new'        : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 0.9 },
    'ratio_cgup_new'          : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cg' : 0.008 },
    'ratio_cgdown_new'        : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cg' : -0.008 },
    'ratio_cbup_new'          : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cb' : 4.0 },
    'ratio_cbdown_new'        : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cb' : -2.0 },
    'ratio_cg003ct12_new'     : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.2, 'cb' : -2.98 ,  'cg' : -0.03 },
    'ratio_cg003ct13sw_new'   : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.3, 'cb' : -0.85 ,  'cg' : -0.03 },
    'ratio_cg003ct14sw_new'   : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.4, 'cb' : 3.31 ,   'cg' : -0.03 },
    'ratio_cg004ct12sw_new'   : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.2, 'cb' : -4.89 ,  'cg' : -0.04 },
    'ratio_cg004ct13sw_new'   : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.3, 'cb' : -3.34 ,  'cg' : -0.04 },
    'ratio_cg004ct15sw_new'   : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.5, 'cb' : 1.88 ,   'cg' : -0.04 },
    'ratio_cg005ct14sw_new'   : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.4, 'cb' : -3.67 ,  'cg' : -0.05 },
    'ratio_cg005ct15sw_new'   : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.5, 'cb' : -1.79 ,  'cg' : -0.05 },
    'ratio_ctcb05_new'        : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cg' : 0.0, 'ct' : 0.5, 'cb' : -7.46 },
    'ratio_ctcb08_new'        : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cg' : 0.0, 'ct' : 0.8, 'cb' : -3.67 },
    'ratio_ctcb09_new'        : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cg' : 0.0, 'ct' : 0.9, 'cb' : -1.79 },
    'ratio_ctcb11_new'        : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cg' : 0.0, 'ct' : 1.1, 'cb' : 3.79 },
    'ratio_ctcb12_new'        : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cg' : 0.0, 'ct' : 1.2, 'cb' : 4.67 },
    'ratio_ctcg01_new'        : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cb' : 1.0, 'ct' : 0.1, 'cg' : 0.075 },
    'ratio_ctcg05_new'        : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cb' : 1.0, 'ct' : 0.5, 'cg' : 0.042 },
    'ratio_ctcg15_new'        : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cb' : 1.0, 'ct' : 1.5, 'cg' : -0.042 },
    'ratio_ctcg2_new'         : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cb' : 1.0, 'ct' : 2.0, 'cg' : -0.083 },
    # 'SM_NLO'                  : { 'firstColumnIsRatio' : False, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    # 'SMmin_NLO'               : { 'firstColumnIsRatio' : False, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    # 'SMmax_NLO'               : { 'firstColumnIsRatio' : False, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    'SM_NNLO'                 : { 'firstColumnIsRatio' : False, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    # 'SMmin_NNLO'              : { 'firstColumnIsRatio' : False, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    # 'SMmax_NNLO'              : { 'firstColumnIsRatio' : False, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    'HRes_mR1mF1.top'         : { 'firstColumnIsRatio' : False, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    'HRes_mR1mF1Q2.top'       : { 'firstColumnIsRatio' : False, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 2.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    'HRes_mR1mF1Qh.top'       : { 'firstColumnIsRatio' : False, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 0.5, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    'HRes_mR1mF2.top'         : { 'firstColumnIsRatio' : False, 'muR' : 1.0, 'muF' : 2.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    'HRes_mR1mFh.top'         : { 'firstColumnIsRatio' : False, 'muR' : 1.0, 'muF' : 0.5, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    'HRes_mR2mF1.top'         : { 'firstColumnIsRatio' : False, 'muR' : 2.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    'HRes_mR2mF2.top'         : { 'firstColumnIsRatio' : False, 'muR' : 2.0, 'muF' : 2.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    'HRes_mRhmF1.top'         : { 'firstColumnIsRatio' : False, 'muR' : 0.5, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    'HRes_mRhmFh.top'         : { 'firstColumnIsRatio' : False, 'muR' : 0.5, 'muF' : 0.5, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    # 
    # Dec 05: High range files
    # 
    'SMNNLOhpt'               : { 'firstColumnIsRatio' : False, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    'SMNNLOmaxhpt'            : { 'firstColumnIsRatio' : False, 'muR' : 99.0, 'muF' : 99.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    'SMNNLOminhpt'            : { 'firstColumnIsRatio' : False, 'muR' : 99.0, 'muF' : 99.0, 'Q' : 1.0, 'ct' : 1.0, 'cb' : 1.0, 'cg' : 0.0 },
    # 
    'ratiovh_cbct_05_new'     : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cg' : 0.0, 'ct' : 0.5, 'cb' : -7.46 },
    'ratiovh_cbct_08_new'     : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cg' : 0.0, 'ct' : 0.8, 'cb' : -3.67 },
    'ratiovh_cbct_09_new'     : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cg' : 0.0, 'ct' : 0.9, 'cb' : -1.79 },
    'ratiovh_cbct_11_new'     : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cg' : 0.0, 'ct' : 1.1, 'cb' : 3.79 },
    'ratiovh_cbct_12_new'     : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cg' : 0.0, 'ct' : 1.2, 'cb' : 4.67 },
    # 
    'ratiovh_ctcg_01_new'     : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cb' : 1.0, 'ct' : 0.1, 'cg' : 0.075 },
    'ratiovh_ctcg_05_new'     : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cb' : 1.0, 'ct' : 0.5, 'cg' : 0.042 },
    'ratiovh_ctcg_15_new'     : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cb' : 1.0, 'ct' : 1.5, 'cg' : -0.042 },
    'ratiovh_ctcg_2_new'      : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cb' : 1.0, 'ct' : 2.0, 'cg' : -0.083 },
    # 
    'ratiovh_cg003ct12_new'   : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.2, 'cb' : -2.98 ,  'cg' : -0.03 },
    'ratiovh_cg003ct13_new'   : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.3, 'cb' : -0.85 ,  'cg' : -0.03 },
    'ratiovh_cg003ct14_new'   : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.4, 'cb' : 3.31 ,   'cg' : -0.03 },
    'ratiovh_cg004ct12_new'   : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.2, 'cb' : -4.89 ,  'cg' : -0.04 },
    'ratiovh_cg004ct13_new'   : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.3, 'cb' : -3.34 ,  'cg' : -0.04 },
    'ratiovh_cg004ct15_new'   : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.5, 'cb' : 1.88 ,   'cg' : -0.04 },
    'ratiovh_cg005ct14_new'   : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.4, 'cb' : -3.67 ,  'cg' : -0.05 },
    'ratiovh_cg005ct15_new'   : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.5, 'cb' : -1.79 ,  'cg' : -0.05 },
    # 
    'ratiovh_ctdown_new'      : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 1.1 },
    'ratiovh_ctup_new'        : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'ct' : 0.9 },
    'ratiovh_cgdown_new'      : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cg' : 0.008 },
    'ratiovh_cgup_new'        : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cg' : -0.008 },
    'ratiovh_cbdown_new'      : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cb' : 4.0 },
    'ratiovh_cbup_new'        : { 'firstColumnIsRatio' : True, 'muR' : 1.0, 'muF' : 1.0, 'Q' : 1.0, 'cb' : -2.0 },
    }

#____________________________________________________________________
def DumpContainerToFile_Top( container, prefix, outdir ):

    # Check if ratios, crosssection and binBoundaries have expected lenghts
    if not(
        len(container.crosssection) == len(container.binBoundaries)-1
        and len(container.ratios) == len(container.binBoundaries)-1
        ):
        Commands.ThrowError((
            'Lists have unexpected lenghts. Found:'
            '\n  file = {3}'
            '\n  len(binBoundaries) = {0}'
            '\n  len(crosssection)  = {1}'
            '\n  len(ratios)        = {2}'
            ).format( len(container.binBoundaries), len(container.crosssection), len(container.ratios), container.file )
            )

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

    if not isdir( outdir ): os.makedirs( outdir )

    with open( container.derivedTheoryFilePath, 'w' ) as outFp:
        w = lambda text: outFp.write( text + ' \n' )

        w( 'file={0}'.format(container.file) )
        if hasattr( container, 'secondfile' ): w( 'secondfile={0}'.format(container.secondfile) )

        if hasattr( container, 'ct' ): w( 'ct={0}'.format(container.ct) )
        if hasattr( container, 'cg' ): w( 'cg={0}'.format(container.cg) )
        if hasattr( container, 'cb' ): w( 'cb={0}'.format(container.cb) )

        if hasattr( container, 'muR' ) and hasattr( container, 'muF' ):
            w( 'muR={0}'.format(container.muR) )
            w( 'muF={0}'.format(container.muF) )

        if hasattr( container, 'Q' ):  w( 'Q={0}'.format(container.Q) )

        w( 'n_binBoundaries={0}'.format(len(container.binBoundaries)) )
        w( 'n_crosssection={0}'.format(len(container.crosssection)) )
        w( 'n_ratios={0}'.format(len(container.ratios)) )

        w( 'binBoundaries={0}'.format( ','.join(map( str, container.binBoundaries )) ) )
        w( 'crosssection={0}'.format( ','.join(map( str, container.crosssection )) ) )
        if hasattr( container, 'binCenters' ): w( 'binCenters={0}'.format( ','.join(map( str, container.binCenters )) ) )
        w( 'ratios={0}'.format( ','.join(map( str, container.ratios )) ) )


#____________________________________________________________________
faultyBinCenters = [ 401., 410., 420., 430., 440., 450., 460., 470., 480., 490., 500., 510., 520., 530., 540., 550., 560., 570., 580., 590., 600., 610., 620., 630., 640., 650., 660., 670., 680., 690., 700., 710., 720., 730., 740., 750., 760., 770., 780., 790. ]
correctedBinCenters = [ 405., 415., 425., 435., 445., 455., 465., 475., 485., 495., 505., 515., 525., 535., 545., 555., 565., 575., 585., 595., 605., 615., 625., 635., 645., 655., 665., 675., 685., 695., 705., 715., 725., 735., 745., 755., 765., 775., 785., 795. ]
def ReadLinesOfTheoryFile_Top( theoryFile, verbose=False, SM=None, isSMFile=False ):
    # Read lines
    with open( theoryFile, 'r' ) as theoryFp:
        lines = [ line.strip() for line in theoryFp.readlines() ]
    commentlines = [ line.strip() for line in lines if line.startswith('#') ]
    lines        = [ line for line in lines if not line.startswith('#') and len(line) > 0 ]

    if isSMFile:

        print '\nFound SM file', theoryFile

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
        elif theoryFile == SM.file:
            print 'This file was found to be the SM file; simply return SM values'
            return SM.pts, SM.crosssection, [ 1.0 for xs in SM.crosssection ]
        else:
            print '\nReading', theoryFile


        pts    = []
        ratios = []
        xss    = []

        # Check for the faulty bin centers
        correctingForFaultyBinCenters = False
        binCenters = [ float(line.split()[0]) for line in lines ]
        if binCenters == faultyBinCenters:
            Commands.Warning( 'Correcting faulty bin centers to corrected bin centers' )
            getcorrectedpt = lambda pt: correctedBinCenters[faultyBinCenters.index(pt)]
            correctingForFaultyBinCenters = True

        for line in lines:

            components = line.split()

            if AgnieszkasFilenameDecoder[ basename(theoryFile).replace('/','') ]['firstColumnIsRatio']:
                print '    Assuming the first column is a RATIO'
                pt         = float(components[0])
                ratio      = float(components[1])
                if correctingForFaultyBinCenters: pt = getcorrectedpt(pt)
                SMxs = SM.crosssection[ SM.pts.index(pt) ]
                xs   = SMxs * ratio
            else:
                print '    Assuming the first column is a CROSS SECTION - ALSO DIVIDING BY 2.27!!'
                pt         = float(components[0])
                xs         = float(components[1]) / 2.27
                if correctingForFaultyBinCenters: pt = getcorrectedpt(pt)
                SMxs  = SM.crosssection[ SM.pts.index(pt) ]
                ratio = xs / SMxs

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
    SM.ratios        = [ 1.0 for xs in SM.crosssection ]

    for key, value in AgnieszkasFilenameDecoder[basename(SM.file).replace('/','')].iteritems():
        setattr( SM, key, value )

    # DumpContainerToFile_Top( SM, prefix='PureSM', outdir=outdir )


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





#____________________________________________________________________
def CreateDerivedTheoryFiles_Top_highPt(
        verbose = True,
        ):

    SMfile = 'suppliedInput/fromAgnieszka/HRes_SMEFT_May16/SM_NNLO'

    lowRangeFiles = glob( 'suppliedInput/fromAgnieszka/SMEFTscaling_May16/*' )
    highRangeFiles = (
        glob( 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_cbct_*_new' )
        + glob( 'suppliedInput/fromAgnieszka/ratio400800_Dec05/ratiovh_ctcg_*_new' )
        )


    # ======================================
    # Make SM variation

    SM = Container()
    SM.file = SMfile

    binCenters, crosssection = ReadLinesOfTheoryFile_Top( SM.file, isSMFile=True, verbose=verbose )
    newBinCenters, binBoundaries, binWidths = TheoryCommands.BinningHeuristic(
        binCenters, manualSwitchAt50=True, manualSwitchAt5=False )

    SM.pts           = binCenters
    SM.binCenters    = deepcopy( newBinCenters )
    SM.binBoundaries = deepcopy( binBoundaries )
    SM.crosssection  = deepcopy( crosssection )
    SM.ratios        = [ 1.0 for xs in SM.crosssection ]

    for key, value in AgnieszkasFilenameDecoder[basename(SM.file).replace('/','')].iteritems():
        setattr( SM, key, value )


    # ======================================
    # Finding pairs of low and high range files

    def f(theoryFile):
        if verbose:
            print '\n\n' + '-'*80
            print 'Processing {0}'.format( theoryFile )
            print '\n'

        shortFileName = basename(theoryFile).replace('/','')
        if not shortFileName in AgnieszkasFilenameDecoder:
            print 'Skipping \'{0}\''.format(shortFileName)
            return

        container = Container()
        container.file = theoryFile

        binCenters, crosssection, ratios = ReadLinesOfTheoryFile_Top( container.file, verbose, SM=SM )

        newBinCenters, binBoundaries, binWidths = TheoryCommands.BinningHeuristic(
            binCenters, manualSwitchAt50=True, manualSwitchAt5=False )

        container.binCenters    = newBinCenters
        container.binBoundaries = binBoundaries
        container.crosssection  = crosssection
        container.ratios        = ratios

        # Load values of coupling into container
        for key, value in AgnieszkasFilenameDecoder[shortFileName].iteritems():
            setattr( container, key, value )

        return container


    lowRangeContainers = [ f(theoryFile) for theoryFile in lowRangeFiles ]
    highRangeContainers = [ f(theoryFile) for theoryFile in highRangeFiles ]

    lowhighPairs = []
    for lowContainer in lowRangeContainers:
        for highContainer in highRangeContainers:
            try:
                if (
                    lowContainer.ct == highContainer.ct
                    and lowContainer.cg == highContainer.cg
                    and lowContainer.cb == highContainer.cb
                    ):
                    print 'Matched {0} with {1}'.format( lowContainer.file, highContainer.file )
                    print '   ( ct = {0}, cg = {1}, cb = {2} )'.format( lowContainer.ct, lowContainer.cg, lowContainer.cb )
                    break
            except AttributeError:
                print 'File {0} or {1} has a coupling problem'.format( lowContainer.file, highContainer.file )
                raise
        else:
            print 'Could not find matching high range container for {0}'.format( lowContainer.file )
            continue

        lowhighPairs.append( [ lowContainer, highContainer ] )


    outdir = 'derivedTheoryFiles_{0}_TopHighPt'.format( datestr )
    if not isdir( outdir ): os.makedirs( outdir )

    print '\n'
    print 'SM.file:'
    print SM.file
    print 'SM.ratios ({0}):'.format( len(SM.ratios) )
    print SM.ratios
    print 'SM.crosssection ({0}):'.format( len(SM.crosssection) )
    print SM.crosssection

    print '\n' + '='*80
    for low, high in lowhighPairs:

        print '\nMerging low {0} with high {1}'.format( low.file, high.file )

        merged = deepcopy( high )
        merged.secondfile = low.file

        assert low.binBoundaries[-1] == high.binBoundaries[0]
        merged.binBoundaries = low.binBoundaries[:-1] + high.binBoundaries
        merged.crosssection  = low.crosssection + high.crosssection
        merged.ratios         = [ xs / SMxs if not SMxs == 0. else 0. for xs, SMxs in zip( merged.crosssection, SM.crosssection ) ]

        # print len(merged.binBoundaries)
        # print merged.binBoundaries
        # print len(merged.crosssection)
        # print merged.crosssection

        # print
        # print len(SM.binBoundaries)
        # print SM.binBoundaries
        # print len(SM.crosssection)
        # print SM.crosssection

        DumpContainerToFile_Top( merged, prefix='TopHighPt', outdir=outdir )
    DumpContainerToFile_Top( SM, prefix='TopHighPt', outdir=outdir )



########################################
# Execute this if called directly
########################################

if __name__ == '__main__':
    print 'Import functions from this file'
    print 'Available functions:'
    functions = [ name for (name, thing) in locals().items() if callable(thing) ]
    print '\n'.join( functions )