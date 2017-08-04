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

import Commands
import TheoryCommands



########################################
# Interface fromPier
########################################


class Parametrization():
    def __init__(self):
        self.verbose = True


    def Parametrize(
            self,
            containers,
            couplingsToParametrize = [ 'kappab', 'kappac' ],
            includeLinearTerms = True,
            ):

        self.couplings = couplingsToParametrize

        couplingCombinations = []
        couplingCombinations += [ list(couplingTuple) for couplingTuple in itertools.combinations( self.couplings, 2 ) ]
        couplingCombinations += [ [ coupling, coupling ] for coupling in self.couplings ]
        if includeLinearTerms:
            couplingCombinations += [ [coupling] for coupling in self.couplings ] + [[]]
        self.couplingCombinations = couplingCombinations

        nComponents = len(couplingCombinations)
        self.nComponents = nComponents

        if len(containers) < nComponents:
            Commands.ThrowError( 'Need at least as much input ({0}) as desired components ({1})'.format( len(containers), nComponents ) )
            sys.exit()
        elif len(containers) > nComponents:
            print '[info] Need only {0} containers; limiting number of containers'.format( nComponents )
            print '[fixme] Now taking the first {0} (may lead to degenerate matrix...)'.format( nComponents )
            # containers = containers[:nComponents]
            random.seed(1001)
            containers = random.sample( containers, nComponents )


        couplingMatrix = []
        for container in containers:
            row = []
            for couplingCombination in couplingCombinations:
                product = 1.0
                for coupling in couplingCombination:
                    product *= getattr( container, coupling )
                row.append( product )
            couplingMatrix.append( row )

        couplingMatrix = numpy.array( couplingMatrix )
        if self.verbose:
            print '\nSquared coupling terms:'
            print '  ', couplingCombinations
            print 'Found the following coupling matrix:'
            print couplingMatrix
        couplingMatrixInv = numpy.linalg.inv( couplingMatrix )

        nBins = len(containers[0].binBoundaries) - 1
        self.nBins = nBins

        # ======================================
        # Create parametrization per bin

        # def make_evaluateParametrization( components ):
        #     if not len(self.couplingCombinations) == len(components):
        #         Commands.ThrowError('Number of coupling combinations ({0}) does not agree with number of components ({1})'.format( len(couplingCombinations), len(components) ) )
        #         sys.exit()

        #     def evaluateParametrization():
        #         res = 0.
        #         for couplingCombination, component in zip( self.couplingCombinations, components ):
        #             product = 1.0
        #             for couping in couplingCombination:
        #                 product *= getattr( self, coupling )
        #             res += product * component
        #         return res

        #     return 

        componentsPerBin = []
        for iBin in xrange(nBins):
            yValues = numpy.array([ [container.crosssection[iBin]] for container in containers ])
            components = list(itertools.chain.from_iterable( couplingMatrixInv.dot( yValues ) ))
            componentsPerBin.append( components )

        self.componentsPerBin = componentsPerBin


    def Evaluate( self, **kwargs ):
        for key, value in kwargs.iteritems():
            setattr( self, key, value )

        xsFromParametrization = []
        for iBin in xrange(self.nBins):
            components = self.componentsPerBin[iBin]

            res = 0.
            for couplingCombination, component in zip( self.couplingCombinations, components ):
                product = 1.0
                for coupling in couplingCombination:
                    product *= getattr( self, coupling )
                res += product * component

            xsFromParametrization.append( res )

        return xsFromParametrization




def MergeGluonAndQuarkInduced(
        verbose = True,
        ):

    # ======================================
    # IO

    gI_theoryFileDir = 'derivedTheoryFiles_gluonInduced_Aug04'
    qI_theoryFileDir = 'derivedTheoryFiles_YukawaQuarkInduced_Aug04'

    qI_theoryFiles = [ theoryFile for theoryFile in glob( join( qI_theoryFileDir, '*.txt' ) ) if not 'muR' in theoryFile ]
    gI_theoryFiles = [ theoryFile for theoryFile in glob( join( gI_theoryFileDir, 'muR_1_muF_1_Q_1_*.txt' ) ) ]

    qI_containers = [ TheoryCommands.ReadDerivedTheoryFile(tF) for tF in qI_theoryFiles ]
    gI_containers = [ TheoryCommands.ReadDerivedTheoryFile(tF) for tF in gI_theoryFiles ]

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
    parametrization.Parametrize( qI_containers )

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


    # HIER VERDER

    # qI moet opgeteld worden bij gI, en dan naar een derivedTheoryFile opslaan
    # Nog niet zeker wat te doen met de range, qI gaat maar tot 120 GeV
    # Kan gewoon de laatste waarde kopieren (gluon induced is niet flat > 120 GeV, alleen de ratio is flat)



def ReadLinesOfYukawaTheoryFile_quarkInduced( theoryFile, verbose=False ):

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
            print '    pt = {0:<8.3f} |  xs = {1:<10.6f}'.format(
                pt, xs
                )

        pts.append( pt )
        xss.append( xs )

    return pts, xss


def CreateDerivedTheoryFiles_YukawaQuarkInduced(
        theoryDir = 'suppliedInput/fromPier/13tev-pth_quarkInduced_Aug04/',
        verbose = True,
        ):

    outdir = 'derivedTheoryFiles_YukawaQuarkInduced_{0}'.format( datestr )
    if not isdir( outdir ): os.makedirs( outdir )

    theoryDir = abspath( theoryDir )
    theoryFiles = glob( join( theoryDir, '*.pth' ) )

    # ======================================
    # Try to find a 'SM' file first so it's possible to calculate ratios

    SM = TheoryCommands.Container()
    SM.file = join( theoryDir, 'higgs_plus_jet_13tev_1_1_mur050_muf050.pth' )
    SM.pt, SM.xs = ReadLinesOfYukawaTheoryFile_quarkInduced( SM.file )
    SM.kappab = 1
    SM.kappac = 1



    pat = r'higgs_plus_jet_13tev_([0-9mp]+)_([0-9mp]+)'
    scaleVarPat = r'mur([0-9]+)_muf([0-9]+)'

    for theoryFile in theoryFiles:
        container = TheoryCommands.Container()
        container.file = theoryFile

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
            container.muR = float( match.group(1) ) / 100
            container.muF = float( match.group(2) ) / 100

        container.pt, container.xs = ReadLinesOfYukawaTheoryFile_quarkInduced( container.file )
        container.ratios = [ xs / smxs if not smxs == 0. else 0. for xs, smxs in zip( container.xs, SM.xs ) ]

        container.binCenters, container.binBoundaries, container.binWidths = TheoryCommands.BinningHeuristic(
            container.pt, manualSwitchAt50 = False, manualSwitchAt5 = True )


        # ======================================
        # Dump container to file

        outname = 'YukawaQuarkInduced_kappab_{0}_kappac_{1}'.format(
            Commands.ConvertFloatToStr( container.kappab ),
            Commands.ConvertFloatToStr( container.kappac )
            )
        if isScaleVariation:
            outname += '_muR_{0}_muF_{1}'.format(
                Commands.ConvertFloatToStr( container.muR ),
                Commands.ConvertFloatToStr( container.muF )
                )
        outname += '.txt'

        with open( join( outdir, outname ), 'w' ) as outFp:
            w = lambda text: outFp.write( text + ' \n' )

            w( 'file={0}'.format(container.file) )
            w( 'kappab={0}'.format(container.kappab) )
            w( 'kappac={0}'.format(container.kappac) )

            if isScaleVariation:
                w( 'muR={0}'.format(container.muR) )
                w( 'muF={0}'.format(container.muF) )

            w( 'binBoundaries={0}'.format( ','.join(map( str, container.binBoundaries )) ) )
            w( 'crosssection={0}'.format( ','.join(map( str, container.xs )) ) )
            w( 'binCenters={0}'.format( ','.join(map( str, container.binCenters )) ) )
            w( 'ratios={0}'.format( ','.join(map( str, container.ratios )) ) )



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
        theoryDir = 'suppliedInput/fromPier/histograms_ggH_May17/',
        mainCrossSection = 'resummed',
        verbose = True,
        ):

    outdir = 'derivedTheoryFiles_gluonInduced_{0}'.format( datestr )
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
            pts, matched_xss, resummed_xss = ReadLinesOfYukawaTheoryFile( theoryFile, verbose )
            newBinCenters, binBoundaries, binWidths = TheoryCommands.BinningHeuristic( pts, manualSwitchAt50=False, manualSwitchAt5=True )

            SM = TheoryCommands.Container()
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
        newBinCenters, binBoundaries, binWidths = TheoryCommands.BinningHeuristic( pts, manualSwitchAt50=False, manualSwitchAt5=True )


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
            
            if mainCrossSection == 'resummed':
                w( 'crosssection={0}'.format(          ','.join( map(str, resummed_xss ) ) ) )
            elif mainCrossSection == 'matched':
                w( 'crosssection={0}'.format(          ','.join( map(str, matched_xss ) ) ) )

            w( 'resummed_crosssection={0}'.format( ','.join( map(str, resummed_xss ) ) ) )
            w( 'matched_crosssection={0}'.format(  ','.join( map(str, matched_xss ) ) ) )

            if smFound:
                if mainCrossSection == 'resummed':
                    ratios = [ xs / SMxs for xs, SMxs in zip( resummed_xss, SM.resummed_xss ) ]
                elif mainCrossSection == 'matched':
                    ratios = [ xs / SMxs for xs, SMxs in zip( matched_xss, SM.matched_xss ) ]
                
                w( 'ratios={0}'.format( ','.join( map(str, ratios ) ) ) )


                print ratios
                print


########################################
# Interface fromAgnieszka
########################################

def CreateDerivedTheoryFiles_Agnieszka(
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