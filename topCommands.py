#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, itertools, operator, re, argparse, sys, random
from math import isnan, isinf
from os.path import *
from glob import glob
from copy import deepcopy
from array import array

import LatestPaths

sys.path.append('src')
import Commands
import PhysicsCommands
import OneOfCommands
import TheoryCommands
import CorrelationMatrices
import MergeHGGWDatacards
import TheoryFileInterface
from Container import Container
from Parametrization import Parametrization, WSParametrization
import PlotCommands

from time import strftime
datestr = strftime( '%b%d' )

import ROOT
from TheoryCommands import c
from TheoryCommands import SaveC
from TheoryCommands import GetPlotBase
from TheoryCommands import SetCMargins


########################################
# Main
########################################

def AppendParserOptions( parser ):

    parser.add_argument( '--topCommands', default=False )
    class CustomAction(argparse.Action):
        def __init__(self, option_strings, dest, **kwargs):
            # Only 1 argument allowed (this is basically the "store_true" action)
            super(CustomAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            # Make sure the parser knows a combineCommand was entered
            setattr( namespace, 'topCommands', True )
            setattr( namespace, self.dest, True )
            # super(FooAction, self).__call__( parser, namespace, values, option_string=None )
            # setattr( namespace, self.dest, values )

    parser.add_argument( '--createDerivedTheoryFiles_Top',                  action=CustomAction )
    parser.add_argument( '--CorrelationMatrices_Top',                       action=CustomAction )
    parser.add_argument( '--couplingT2WS_Top',                              action=CustomAction )
    parser.add_argument( '--couplingBestfit_Top',                           action=CustomAction )
    parser.add_argument( '--coupling2Dplot_Top',                            action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Top',                       action=CustomAction )
    parser.add_argument( '--checkWSParametrization_Top',                    action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Top_lumiStudy',             action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Top_BRdependencyComparison',  action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Top_BRdependencyComparison_bigRange',  action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Top_ProfiledTotalXS',       action=CustomAction )

    parser.add_argument( '--couplingContourPlot_Top_onlyNormalization',     action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Top_skippedLastBin',        action=CustomAction )
    parser.add_argument( '--couplingContourPlot_Top_bigRange',              action=CustomAction )

    parser.add_argument( '--couplingT2WS_TopCtCb',                          action=CustomAction )

    parser.add_argument( '--couplingContourPlot_TopCtCb',                   action=CustomAction )

    parser.add_argument( '--combinationAndContour_TopCtCb',                 action=CustomAction )
    parser.add_argument( '--combinationAndContour_Top',                     action=CustomAction )


########################################
# Methods
########################################    

def main( args ):

    # # expBinBoundaries = [ 0., 15., 30., 45., 85., 125. ]
    expBinBoundaries = [ 0., 15., 30., 45., 85., 125., 200., 350., 10000. ]
    print 'Hardcoded binBoundaries for Top:'
    print expBinBoundaries
    print ''

    # LATESTDATACARD_Top            = 'suppliedInput/combinedCard_Jul26.txt'

    # LATESTWORKSPACE_Top           = 'workspaces_Aug11/combinedCard_Jul26_CouplingModel_Top_withTheoryUncertainties.root'

    # LATESTCORRELATIONMATRIX_Top   = 'plots_CorrelationMatrices_Aug11_Top/corrMat_exp.txt'
    # LATESTTHEORYUNCERTAINTIES_Top = 'plots_CorrelationMatrices_Aug11_Top/errors_for_corrMat_exp.txt'

    TheoryCommands.SetPlotDir( 'plots_{0}_Top'.format(datestr) )


    #____________________________________________________________________
    if args.createDerivedTheoryFiles_Top:
        TheoryFileInterface.CreateDerivedTheoryFiles_Top(
            verbose = True,
            )


    #____________________________________________________________________
    if args.CorrelationMatrices_Top:

        CorrelationMatrices.SetPlotDir( 'correlationMatrices_{0}_Top'.format(datestr) )

        variationFiles = TheoryFileInterface.FileFinder(
            directory = LatestPaths.derivedTheoryFiles_Top,
            ct = 1.0, cg = 1.0, cb = 1.0
            )

        variations = [
            TheoryFileInterface.ReadDerivedTheoryFile( variationFile ) for variationFile in variationFiles ]

        CorrelationMatrices.GetCorrelationMatrix(
            variations,
            makeScatterPlots          = False,
            makeCorrelationMatrixPlot = True,
            outname                   = 'corrMat_theory',
            verbose                   = True,
            )

        variations_expbinning = deepcopy(variations)
        for variation in variations_expbinning:
            TheoryCommands.RebinDerivedTheoryContainer( variation, expBinBoundaries )

        CorrelationMatrices.GetCorrelationMatrix(
            variations_expbinning,
            makeScatterPlots          = True,
            makeCorrelationMatrixPlot = True,
            outname                   = 'corrMat_exp',
            verbose                   = True,
            )

    #____________________________________________________________________
    if args.couplingT2WS_Top:

        # ======================================
        # Switches

        INCLUDE_THEORY_UNCERTAINTIES   = True
        # INCLUDE_THEORY_UNCERTAINTIES   = False

        # MAKELUMISCALABLE               = True
        MAKELUMISCALABLE               = False

        # INCLUDE_BR_COUPLING_DEPENDENCY = True
        INCLUDE_BR_COUPLING_DEPENDENCY = False

        # PROFILE_TOTAL_XS               = True
        PROFILE_TOTAL_XS               = False

        # FIT_ONLY_NORMALIZATION         = True
        FIT_ONLY_NORMALIZATION         = False


        EXCLUDE_LAST_BIN               = True
        # EXCLUDE_LAST_BIN               = False


        # ======================================

        datacard = LatestPaths.card_combined_ggHxH_PTH
        if args.hgg:
            datacard = LatestPaths.card_hgg_ggHxH_PTH
        if args.hzz:
            datacard = LatestPaths.card_hzz_ggHxH_PTH

        TheoryFileInterface.SetFileFinderDir( LatestPaths.derivedTheoryFiles_Top )

        if INCLUDE_THEORY_UNCERTAINTIES:
            correlationMatrix   = LatestPaths.correlationMatrix_Top
            theoryUncertainties = LatestPaths.theoryUncertainties_Top
            # if UNCORRELATED_THEORY_UNCERTAINTIES:
            #     correlationMatrix = LatestPaths.correlationMatrix_Top_Uncorrelated


        # ======================================
        # 

        extraOptions = [
            '--PO verbose=2',
            # '--PO verbose=0',
            '--PO \'higgsMassRange=123,127\'',
            '--PO linearTerms=False',
            '--PO splitggH=True',
            ]

        if args.hzz:
            extraOptions.append( '--PO isOnlyHZZ=True' )
        if args.hgg:
            extraOptions.append( '--PO isOnlyHgg=True' )

        extraOptions.append(
            '--PO SM=[ct=1,cg=0,file={0}]'.format(
                TheoryFileInterface.FileFinder( ct=1, cg=1, cb=1, muR=1, muF=1, Q=1, expectOneFile=True )
                )
            )
        
        theoryFiles = TheoryFileInterface.FileFinder(
            ct='*', cg='*', muR=1, muF=1, Q=1, filter='cb'
            )
        possibleTheories = []
        for theoryFile in theoryFiles:
            ct = Commands.ConvertStrToFloat( re.search( r'ct_([\dmp]+)', theoryFile ).group(1) )
            cg = Commands.ConvertStrToFloat( re.search( r'cg_([\dmp]+)', theoryFile ).group(1) )
            possibleTheories.append(
                '--PO theory=[ct={0},cg={1},file={2}]'.format(
                    ct, cg, theoryFile
                    )                
                )
        extraOptions.extend(possibleTheories)


        suffix = 'Top'
        if INCLUDE_THEORY_UNCERTAINTIES:
            extraOptions.extend([
                '--PO correlationMatrix={0}'.format(  LatestPaths.correlationMatrix_Top ),
                '--PO theoryUncertainties={0}'.format(  LatestPaths.theoryUncertainties_Top ),
                '--PO skipOverflowBinTheoryUncertainty=True'
                ])
            suffix += '_withTheoryUncertainties'
        else:
            suffix += '_noTheoryUncertainties'

        if MAKELUMISCALABLE:
            extraOptions.append( '--PO lumiScale=True'  )
            suffix += '_lumiScale'

        if PROFILE_TOTAL_XS:
            extraOptions.append( '--PO ProfileTotalXS=True' )
            suffix += '_profiledTotalXS'

        if INCLUDE_BR_COUPLING_DEPENDENCY:
            extraOptions.append( '--PO FitBR=True' )
            suffix += '_couplingDependentBR'

        if FIT_ONLY_NORMALIZATION:
            extraOptions.append( '--PO FitOnlyNormalization=True' )
            suffix += '_fitOnlyNormalization'


        if EXCLUDE_LAST_BIN:
            expBinBoundaries = expBinBoundaries[:-1]
            suffix += '_skippedLastBin'

        extraOptions.append(
            '--PO binBoundaries={0}'.format( ','.join([ str(b) for b in expBinBoundaries ]) )
            )

        Commands.BasicT2WSwithModel(
            datacard,
            'CouplingModel.py',
            suffix = suffix,
            extraOptions = extraOptions,
            )


    #____________________________________________________________________
    if args.couplingT2WS_TopCtCb:


        # ======================================
        # Switches

        INCLUDE_THEORY_UNCERTAINTIES   = True
        # INCLUDE_THEORY_UNCERTAINTIES   = False

        # MAKELUMISCALABLE               = True
        MAKELUMISCALABLE               = False

        # INCLUDE_BR_COUPLING_DEPENDENCY = True
        INCLUDE_BR_COUPLING_DEPENDENCY = False

        # PROFILE_TOTAL_XS               = True
        PROFILE_TOTAL_XS               = False


        # ======================================

        datacard = LatestPaths.card_combined_ggHxH_PTH
        if args.hgg:
            datacard = LatestPaths.card_hgg_ggHxH_PTH
        if args.hzz:
            datacard = LatestPaths.card_hzz_ggHxH_PTH

        TheoryFileInterface.SetFileFinderDir( LatestPaths.derivedTheoryFiles_Top )

        if INCLUDE_THEORY_UNCERTAINTIES:
            correlationMatrix   = LatestPaths.correlationMatrix_Top
            theoryUncertainties = LatestPaths.theoryUncertainties_Top


        # ======================================
        # 

        extraOptions = [
            '--PO verbose=2',
            '--PO \'higgsMassRange=123,127\'',
            '--PO linearTerms=False',  # Check this, no other worthwhile contributions?
            '--PO splitggH=True',
            ]

        if args.hzz:
            extraOptions.append( '--PO isOnlyHZZ=True' )
        if args.hgg:
            extraOptions.append( '--PO isOnlyHgg=True' )

        extraOptions.append(
            '--PO SM=[ct=1,cb=1,file={0}]'.format(
                TheoryFileInterface.FileFinder( ct=1, cb=1, cg=1, muR=1, muF=1, Q=1, expectOneFile=True )
                )
            )
        
        theoryFiles = TheoryFileInterface.FileFinder(
            ct='*', cb='*', muR=1, muF=1, Q=1, filter='cg'
            )
        possibleTheories = []
        for theoryFile in theoryFiles:
            ct = Commands.ConvertStrToFloat( re.search( r'ct_([\dmp]+)', theoryFile ).group(1) )
            cb = Commands.ConvertStrToFloat( re.search( r'cb_([\dmp]+)', theoryFile ).group(1) )
            possibleTheories.append(
                '--PO theory=[ct={0},cb={1},file={2}]'.format(
                    ct, cb, theoryFile
                    )                
                )
        extraOptions.extend(possibleTheories)


        suffix = 'TopCtCb'
        if INCLUDE_THEORY_UNCERTAINTIES:
            extraOptions.extend([
                '--PO correlationMatrix={0}'.format(  LatestPaths.correlationMatrix_Top ),
                '--PO theoryUncertainties={0}'.format(  LatestPaths.theoryUncertainties_Top ),
                '--PO skipOverflowBinTheoryUncertainty=True'
                ])
            suffix += '_withTheoryUncertainties'
        else:
            suffix += '_noTheoryUncertainties'

        if MAKELUMISCALABLE:
            extraOptions.append( '--PO lumiScale=True'  )
            suffix += '_lumiScale'

        if PROFILE_TOTAL_XS:
            extraOptions.append( '--PO ProfileTotalXS=True' )
            suffix += '_profiledTotalXS'

        if INCLUDE_BR_COUPLING_DEPENDENCY:
            extraOptions.append( '--PO FitBR=True' )
            suffix += '_couplingDependentBR'


        # Scale these bins with 1.0 regardless of parametrization
        # extraOptions.append(
        #     '--PO skipBins=GT350'
        #     )
        extraOptions.append(
            '--PO binBoundaries={0}'.format( ','.join([ str(b) for b in expBinBoundaries ]) )
            )

        Commands.BasicT2WSwithModel(
            datacard,
            'CouplingModel.py',
            suffix = suffix,
            extraOptions = extraOptions,
            )




    #____________________________________________________________________
    if args.couplingBestfit_Top:

        # ======================================
        # Manage flags

        doFastscan = False
        if args.notFastscan: doFastscan = False
        if args.fastscan:    doFastscan = True

        ASIMOV = False
        if args.asimov: ASIMOV = True

        # LUMISTUDY                      = True
        LUMISTUDY                      = False

        # PROFILE_TOTAL_XS               = True
        PROFILE_TOTAL_XS               = False

        # FIT_ONLY_NORMALIZATION         = True
        FIT_ONLY_NORMALIZATION         = False

        # INCLUDE_BR_COUPLING_DEPENDENCY = True
        INCLUDE_BR_COUPLING_DEPENDENCY = False

        # FIX_KAPPAV                     = True
        FIX_KAPPAV                     = False

        # MAX_KAPPAV_ONE                 = True
        MAX_KAPPAV_ONE                 = False


        # DO_KAPPAT_KAPPAB               = True
        DO_KAPPAT_KAPPAB               = False

        # DO_ONLY_ONE_KAPPA              = True
        DO_ONLY_ONE_KAPPA              = False


        # EXCLUDE_LAST_BIN               = True
        EXCLUDE_LAST_BIN               = False


        if not INCLUDE_BR_COUPLING_DEPENDENCY and FIX_KAPPAV:
            Commands.ThrowError( 'INCLUDE_BR_COUPLING_DEPENDENCY == False and FIX_KAPPAV == True is not allowed' )
        if FIX_KAPPAV and MAX_KAPPAV_ONE:
            Commands.ThrowError( 'FIX_KAPPAV == True and MAX_KAPPAV_ONE == True is not allowed' )


        # ======================================
        # Select right datacard

        datacard = LatestPaths.ws_combined_Top
        if args.hgg:
            datacard = LatestPaths.ws_hgg_Top
        if args.hzz:
            datacard = LatestPaths.ws_hzz_Top

        if ( LUMISTUDY ) and ( args.hzz or args.hgg ):
            print '[fixme] These studies not implemented for hgg or hzz'
            sys.exit()
        if LUMISTUDY:
            datacard = LatestPaths.ws_combined_Top_lumiScalable
        if PROFILE_TOTAL_XS:
            datacard = LatestPaths.ws_combined_Top_profiledTotalXS
        if INCLUDE_BR_COUPLING_DEPENDENCY:
            datacard = LatestPaths.ws_combined_Top_couplingDependentBR
        if FIT_ONLY_NORMALIZATION:
            datacard = LatestPaths.ws_combined_Top_profiledTotalXS_fitOnlyNormalization
        if EXCLUDE_LAST_BIN:
            datacard = LatestPaths.ws_combined_Top_skippedLastBin

        if DO_KAPPAT_KAPPAB:
            datacard = LatestPaths.ws_combined_TopCtCb
            if args.hgg:
                datacard = LatestPaths.ws_hgg_TopCtCb
            if args.hzz:
                datacard = LatestPaths.ws_hzz_TopCtCb


        # ======================================
        # Set some job specifics (ranges, number of points)

        jobDirectory = 'Scan_Top_{0}'.format( datestr )
        if DO_KAPPAT_KAPPAB: jobDirectory = 'Scan_TopCtCb_{0}'.format( datestr )
        if args.hgg: jobDirectory += '_hgg'
        if args.hzz: jobDirectory += '_hzz'

        ct_ranges = [ -1., 2. ]
        cg_ranges = [ -0.1, 0.2 ]

        if DO_KAPPAT_KAPPAB:
            ct_ranges = [ -0.1, 2. ]
            cb_ranges = [ -10.0, 14.0 ]

        if doFastscan:
            nPoints = 6400
            nPointsPerJob = 1600
            queue = 'short.q'
        else:
            nPoints = 6400
            nPointsPerJob = 20
            # nPoints = 4900
            # nPointsPerJob = 14
            queue = 'all.q'
            # if ASIMOV:
            #     # This does not converge in time
            #     queue = 'short.q'
            #     nPoints = 1600
            #     nPointsPerJob = 2
            if args.hzz:
                nPointsPerJob = 320
                queue = 'short.q'


        print 'WARNING ' * 7
        print 'TEMPORARY CODE - REMOVE THIS'
        print 'ALSO TURN BACK ON squareDistPoiStep !!'
        # nPoints = 6400 * 4
        # nPoints = 10000
        nPoints   = 60*60
        ct_ranges = [ -3.5, 7.5 ]
        cg_ranges = [ -0.65, 0.5 ]


        # ======================================
        # Construct the fit command and process flags

        # --------------------
        # Setting POIs

        POIs = [ 'ct', 'cg' ]
        if DO_KAPPAT_KAPPAB:
            POIs = [ 'ct', 'cb' ]


        # --------------------
        # Setting physicsModelParameters

        physicsModelParameters = [ 'ct=1.0', 'cg=0.0' ]
        if DO_KAPPAT_KAPPAB:
            physicsModelParameters = [ 'ct=1.0', 'cb=1.0' ]
        if LUMISTUDY:
            physicsModelParameters.append( 'lumiScale=8.356546' )
        if INCLUDE_BR_COUPLING_DEPENDENCY:
            physicsModelParameters.append( 'kappa_V=0.99' )


        # --------------------
        # Setting physicsModelParameterRanges

        physicsModelParameterRanges = [
            'ct={0},{1}'.format( ct_ranges[0], ct_ranges[1] ),
            'cg={0},{1}'.format( cg_ranges[0], cg_ranges[1] )
            ]
        if DO_KAPPAT_KAPPAB:
            # Actually overwrite list contents
            physicsModelParameterRanges = [
                'ct={0},{1}'.format( ct_ranges[0], ct_ranges[1] ),
                'cb={0},{1}'.format( cb_ranges[0], cb_ranges[1] )
                ]
        if INCLUDE_BR_COUPLING_DEPENDENCY and not FIX_KAPPAV:
            if MAX_KAPPAV_ONE:
                physicsModelParameterRanges.append( 'kappa_V=-100.0,1.0' )
            else:
                physicsModelParameterRanges.append( 'kappa_V=-100.0,100.0' )


        # --------------------
        # Specify floating and frozen nuisances

        floatNuisances  = []
        freezeNuisances = []

        if INCLUDE_BR_COUPLING_DEPENDENCY:
            if FIX_KAPPAV:
                freezeNuisances.append( 'kappa_V' )
            else:
                floatNuisances.append( 'kappa_V' )

        if DO_ONLY_ONE_KAPPA:
            floatNuisances.append( theOtherKappa )


        # --------------------
        # Construct extraOptions

        extraOptions = [
            # '--squareDistPoiStep',
            '-P ' + ' -P '.join(POIs),
            '--setPhysicsModelParameters '      + ','.join(physicsModelParameters),
            '--setPhysicsModelParameterRanges ' + ':'.join(physicsModelParameterRanges),
            ]

        if len(floatNuisances) > 0:
            extraOptions.append( '--floatNuisances ' + ','.join(floatNuisances) )
        if len(freezeNuisances) > 0:
            extraOptions.append( '--freezeNuisances ' + ','.join(freezeNuisances) )


        # --------------------
        # Compile list of variables to save

        variablesToSave = []
        variablesToSave.extend( Commands.ListSet( datacard, 'yieldParameters' ) )
        variablesToSave.extend( [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ] )
        if INCLUDE_BR_COUPLING_DEPENDENCY:
            variablesToSave.extend( Commands.ListSet( datacard, 'hgg_yieldParameters' ) )
            variablesToSave.extend( Commands.ListSet( datacard, 'hzz_yieldParameters' ) )
            variablesToSave.extend( Commands.ListSet( datacard, 'BRvariables' ) )
        if PROFILE_TOTAL_XS:
            variablesToSave.extend([ 'r_totalXS', 'totalXSmodifier', 'totalXS_SM', 'totalXS' ])
            if FIT_ONLY_NORMALIZATION:
                variablesToSave.append( 'globalTotalXSmodifier' )
        extraOptions.append( '--saveSpecifiedFunc ' + ','.join(variablesToSave) )


        # ======================================
        # Appropriately name scan, create jobDirectoy and submit command

        if PROFILE_TOTAL_XS:
            jobDirectory += '_profiledTotalXS'
            if FIT_ONLY_NORMALIZATION:
                jobDirectory += '_fitOnlyNormalization'
        if LUMISTUDY:
            jobDirectory += '_lumiStudy'
        if INCLUDE_BR_COUPLING_DEPENDENCY:
            jobDirectory += '_couplingDependentBR'
            if FIX_KAPPAV:
                jobDirectory += '_fixedKappaV'
            elif MAX_KAPPAV_ONE:
                jobDirectory += '_kappaVMaxOne'
        if EXCLUDE_LAST_BIN:
            jobDirectory += '_skippedLastBin'
        if ASIMOV:                          jobDirectory += '_asimov'
        if doFastscan:                      jobDirectory += '_fastscan'

        jobDirectory = Commands.AppendNumberToDirNameUntilItDoesNotExistAnymore( jobDirectory )
        if Commands.IsTestMode(): print '\nWould now create new directory: {0}'.format( basename(jobDirectory) )

        Commands.MultiDimCombineTool(
            datacard,
            nPoints       = nPoints,
            nPointsPerJob = nPointsPerJob,
            queue         = queue,
            notOnBatch    = False,
            jobDirectory  = jobDirectory,
            fastscan      = doFastscan,
            asimov        = ASIMOV,
            jobPriority   = 0,
            extraOptions  = extraOptions
            )


    #____________________________________________________________________
    if args.coupling2Dplot_Top:

        # rootfiles = glob( LatestPaths.scan_combined_Top_couplingDependentBR_kappaVMaxOne_asimov + '/*.root' )

        # rootfiles = glob( 'Scan_Top_Nov23_hzz_0' + '/*.root' )
        rootfiles = glob( 'Scan_Top_Nov23_hzz_asimov_2' + '/*.root' )

        res = TheoryCommands.PlotCouplingScan2D(
            rootfiles,
            xCoupling = 'ct',
            yCoupling = 'cg',
            SM = ( 1.0, 0.0 ),
            )

        print '\nBest fit:'
        print res.xCoupling, '=', res.xBestfit
        print res.yCoupling, '=', res.yBestfit



    #____________________________________________________________________

    xCoupling = 'ct'
    yCoupling = 'cg'
    titles = {
        'ct'          : '#kappa_{t}',
        'cg'          : '#kappa_{g}',
        'cb'          : '#kappa_{b}',
        'hgg'         : 'H #rightarrow #gamma#gamma',
        'hzz'         : 'H #rightarrow 4l',
        'combined'    : 'Combination',
        # For lumi study:
        'regularLumi' : 'Expected at 35.9 fb^{-1}',
        'highLumi'    : 'Expected at 300 fb^{-1}',
        }


    #____________________________________________________________________
    if args.couplingContourPlot_Top:

        ASIMOV = False
        if args.asimov: ASIMOV = True

        # EXTENDED_RANGE = True
        EXTENDED_RANGE = False

        if not ASIMOV:
            combined_rootfiles = glob( '{0}/*.root'.format(LatestPaths.scan_combined_Top) )
            hzz_rootfiles      = glob( '{0}/*.root'.format(LatestPaths.scan_hzz_Top) )
            hgg_rootfiles      = glob( '{0}/*.root'.format(LatestPaths.scan_hgg_Top) )
        else:
            combined_rootfiles = glob( '{0}/*.root'.format(LatestPaths.scan_combined_Top_asimov) )
            hzz_rootfiles      = glob( '{0}/*.root'.format(LatestPaths.scan_hzz_Top_asimov) )
            hgg_rootfiles      = glob( '{0}/*.root'.format(LatestPaths.scan_hgg_Top_asimov) )


        if EXTENDED_RANGE:
            combined_rootfiles = glob( '{0}/*.root'.format(LatestPaths.scan_combined_Top_extendedRange_asimov) )
            hzz_rootfiles      = glob( '{0}/*.root'.format(LatestPaths.scan_hzz_Top_extendedRange_asimov) )
            hgg_rootfiles      = glob( '{0}/*.root'.format(LatestPaths.scan_hgg_Top_extendedRange_asimov) )



        hgg = TheoryCommands.GetTH2FromListOfRootFiles(
            hgg_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        hgg.color = 2
        hgg.name = 'hgg'

        hzz = TheoryCommands.GetTH2FromListOfRootFiles(
            hzz_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        hzz.color = 4
        hzz.name = 'hzz'

        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name = 'combined'

        containers = [
            hgg,
            hzz,
            combined
            ]
        for container in containers: container.title = titles.get( container.name, container.name )

        PlotCommands.BasicMixedContourPlot(
            containers,
            # xMin = -0.2,
            # xMax = 2.0,
            # yMin = -0.08,
            # yMax = 0.135,
            xMin = -1.0,  xMax = 2.0,
            yMin = -0.1,  yMax = 0.2,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_Top' + ( '_asimov' if ASIMOV else '' ) + ( '_extendedRange' if ASIMOV else '' ),
            x_SM      = 1.,
            y_SM      = 0.,
            plotIndividualH2s = True,
            )


    #____________________________________________________________________
    if args.couplingContourPlot_Top_skippedLastBin:

        ASIMOV = False
        if args.asimov: ASIMOV = True

        containers = []

        combined_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top ) )
        if ASIMOV: combined_rootfiles = glob( '{0}/*.root'.format(LatestPaths.scan_combined_Top_asimov) )
        combined = TheoryCommands.GetTH2FromListOfRootFiles( combined_rootfiles, xCoupling, yCoupling, verbose   = False, )
        combined.color = 1
        combined.name  = 'combined'
        combined.title = 'Combination'
        containers.append(combined)

        skippedLastBin_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top_skippedLastBin ) )
        if ASIMOV: skippedLastBin_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top_skippedLastBin_asimov ) )
        skippedLastBin = TheoryCommands.GetTH2FromListOfRootFiles( skippedLastBin_rootfiles, xCoupling, yCoupling, verbose   = False, )
        skippedLastBin.color = 2
        skippedLastBin.name  = 'skippedLastBin'
        skippedLastBin.title = 'Last bin not fitted'
        containers.append(skippedLastBin)

        PlotCommands.BasicMixedContourPlot(
            containers,
            # xMin = -0.2,
            # xMax = 2.0,
            # yMin = -0.08,
            # yMax = 0.135,
            xMin = -1.0,
            xMax = 2.0,
            yMin = -0.1,
            yMax = 0.2,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_skippedLastBin' + ( '_asimov' if ASIMOV else '' ),
            x_SM      = 1.,
            y_SM      = 0.,
            plotIndividualH2s = True,
            )


    #____________________________________________________________________
    if args.couplingContourPlot_Top_bigRange:

        containers = []

        # combined_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top_bigRange ) )
        # ctMin = -2.0,
        # ctMax = 2.5,
        # cgMin = -0.5,
        # cgMax = 0.25,

        combined_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top_bigRange2 ) )
        ctMin = -2.0,
        ctMax = 6.5,
        cgMin = -0.45,
        cgMax = 0.25,

        combined = TheoryCommands.GetTH2FromListOfRootFiles( combined_rootfiles, xCoupling, yCoupling, verbose   = False, )
        combined.color = 1
        combined.name  = 'combined'
        combined.title = 'Combination'
        containers.append(combined)

        PlotCommands.BasicMixedContourPlot(
            containers,
            xMin = -2.0,
            xMax = 6.5,
            yMin = -0.45,
            yMax = 0.25,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_bigRange',
            x_SM      = 1.,
            y_SM      = 0.,
            plotIndividualH2s = True,
            )

    #____________________________________________________________________
    if args.couplingContourPlot_Top_onlyNormalization:

        containers = []

        # ctMin = -2.0
        # ctMax = 6.5
        # cgMin = -0.45
        # cgMax = 0.25

        ctMin = -1.0
        ctMax = 2.0
        cgMin = -0.1
        cgMax = 0.2

        combined_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top_asimov ) )
        combined = TheoryCommands.GetTH2FromListOfRootFiles( combined_rootfiles, xCoupling, yCoupling, verbose   = False, )
        combined.color = 1
        combined.name  = 'combined'
        combined.title = 'Nominal'
        containers.append(combined)

        onlyNormalization_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top_fitOnlyNormalization_asimov ) )
        onlyNormalization = TheoryCommands.GetTH2FromListOfRootFiles( onlyNormalization_rootfiles, xCoupling, yCoupling, verbose   = False, )
        onlyNormalization.color = 2
        onlyNormalization.name  = 'onlyNormalization'
        onlyNormalization.title = 'Only normalization'
        containers.append(onlyNormalization)

        PlotCommands.BasicMixedContourPlot(
            containers,
            xMin = ctMin,
            xMax = ctMax,
            yMin = cgMin,
            yMax = cgMax,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_onlyNormalization',
            x_SM      = 1.,
            y_SM      = 0.,
            plotIndividualH2s = True,
            )

    #____________________________________________________________________
    if args.couplingContourPlot_TopCtCb:
        containers = []

        xCoupling = 'ct'
        yCoupling = 'cb'


        ASIMOV = False
        if args.asimov: ASIMOV = True

        if not ASIMOV:
            combined_rootfiles = glob( '{0}/*.root'.format(LatestPaths.scan_combined_TopCtCb) )
            hzz_rootfiles      = glob( '{0}/*.root'.format(LatestPaths.scan_hzz_TopCtCb) )
            hgg_rootfiles      = glob( '{0}/*.root'.format(LatestPaths.scan_hgg_TopCtCb) )
        else:
            combined_rootfiles = glob( '{0}/*.root'.format(LatestPaths.scan_combined_TopCtCb_asimov) )
            hzz_rootfiles      = glob( '{0}/*.root'.format(LatestPaths.scan_hzz_TopCtCb_asimov) )
            hgg_rootfiles      = glob( '{0}/*.root'.format(LatestPaths.scan_hgg_TopCtCb_asimov) )


        hgg = TheoryCommands.GetTH2FromListOfRootFiles(
            hgg_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        hgg.color = 2
        hgg.name = 'hgg'
        hgg.title = 'H #rightarrow #gamma#gamma'
        containers.append(hgg)

        hzz = TheoryCommands.GetTH2FromListOfRootFiles(
            hzz_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        hzz.color = 4
        hzz.name  = 'hzz'
        hzz.title = 'H #rightarrow ZZ'
        containers.append(hzz)

        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name  = 'combined'
        combined.title = 'Combination'
        containers.append(combined)


        PlotCommands.BasicMixedContourPlot(
            containers,
            # xMin = -0.2,
            # xMax = 2.0,
            # yMin = -0.08,
            # yMax = 0.135,
            xMin = -1.0,  xMax = 2.0,
            yMin = -10.0,  yMax = 15.0,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_TopCtCb' + ( '_asimov' if ASIMOV else '' ),
            x_SM      = 1.,
            y_SM      = 1.,
            plotIndividualH2s = True,
            )


    #____________________________________________________________________
    if args.couplingContourPlot_Top_lumiStudy:

        containers = []

        ctMin = -1.0
        ctMax = 2.0
        cgMin = -0.1
        cgMax = 0.2

        combined_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top_asimov ) )
        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name  = 'regularLumi'
        combined.title = '35.9 fb^{-1}'
        containers.append(combined)


        combined_lum8_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top_lumiStudy_asimov ) )
        combined_lum8 = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_lum8_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined_lum8.color = 2
        combined_lum8.name  = 'highLumi'
        combined_lum8.title = '300 fb^{-1}'
        containers.append(combined_lum8)


        PlotCommands.BasicMixedContourPlot(
            containers,
            # xMin = -0.2,
            # xMax = 2.0,
            # yMin = -0.08,
            # yMax = 0.135,
            xMin = ctMin,
            xMax = ctMax,
            yMin = cgMin,
            yMax = cgMax,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_LumiStudy',
            x_SM      = 1.,
            y_SM      = 0.,
            plotIndividualH2s = True,
            )

    #____________________________________________________________________
    if args.couplingContourPlot_Top_ProfiledTotalXS:

        containers = []

        combined_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top_asimov ) )
        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name = 'regular'
        combined.title = 'Nominal'
        containers.append(combined)

        profiledTotalXS_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top_profiledTotalXS_asimov ) )
        profiledTotalXS = TheoryCommands.GetTH2FromListOfRootFiles(
            profiledTotalXS_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        profiledTotalXS.color = 2
        profiledTotalXS.name = 'profiledTotalXS'
        profiledTotalXS.title = '#sigma profiled'
        containers.append(profiledTotalXS)

        PlotCommands.BasicMixedContourPlot(
            containers,
            xMin = -0.2,
            xMax = 2.0,
            yMin = -0.08,
            yMax = 0.135,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_profiledTotalXS',
            x_SM      = 1.,
            y_SM      = 0.,
            plotIndividualH2s = True,
            )


    #____________________________________________________________________
    if args.couplingContourPlot_Top_BRdependencyComparison:

        containers = []

        xCoupling = 'ct'
        yCoupling = 'cg'
        titles = { 'ct': '#kappa_{t}', 'cg' : '#kappa_{g}' }

        combined_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top_asimov) )
        combined = TheoryCommands.GetTH2FromListOfRootFiles(
            combined_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        combined.color = 1
        combined.name = 'regularBR'
        combined.title = 'BR constant'
        containers.append( combined )

        # scalingBR_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top_couplingDependentBR_asimov ) )
        # scalingBR = TheoryCommands.GetTH2FromListOfRootFiles(
        #     scalingBR_rootfiles,
        #     xCoupling,
        #     yCoupling,
        #     verbose   = False,
        #     )
        # scalingBR.color = 2
        # scalingBR.name = 'scalingBR'
        # scalingBR.title = 'BR(#kappa_{t}, #kappa_{V})'
        # containers.append( scalingBR )

        scalingBRfixedKappaV_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top_couplingDependentBR_fixedKappaV_asimov ) )
        scalingBRfixedKappaV = TheoryCommands.GetTH2FromListOfRootFiles(
            scalingBRfixedKappaV_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        scalingBRfixedKappaV.color = 4
        scalingBRfixedKappaV.name = 'scalingBRfixedKappaV'
        scalingBRfixedKappaV.title = 'BR(#kappa_{t}) (#kappa_{V} fixed)'
        containers.append( scalingBRfixedKappaV )

        # scalingBRkappaVMaxOne_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top_couplingDependentBR_kappaVMaxOne_asimov ) )
        # scalingBRkappaVMaxOne = TheoryCommands.GetTH2FromListOfRootFiles(
        #     scalingBRkappaVMaxOne_rootfiles,
        #     xCoupling,
        #     yCoupling,
        #     verbose   = False,
        #     )
        # scalingBRkappaVMaxOne.color = 8
        # scalingBRkappaVMaxOne.name = 'scalingBRkappaVMaxOne'
        # scalingBRkappaVMaxOne.title = 'BR(#kappa_{t}, #kappa_{V}#leq1)'
        # containers.append( scalingBRkappaVMaxOne )


        PlotCommands.BasicMixedContourPlot(
            containers,
            xMin = -0.2,
            xMax = 2.0,
            yMin = -0.08,
            yMax = 0.2,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_BRcouplingDependency',
            x_SM      = 1.,
            y_SM      = 0.,
            plotIndividualH2s = True,
            )



    #____________________________________________________________________
    if args.couplingContourPlot_Top_BRdependencyComparison_bigRange:

        containers = []

        xCoupling = 'ct'
        yCoupling = 'cg'
        titles = { 'ct': '#kappa_{t}', 'cg' : '#kappa_{g}' }

        # combined_rootfiles  = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top_asimov) )
        # combined = TheoryCommands.GetTH2FromListOfRootFiles(
        #     combined_rootfiles,
        #     xCoupling,
        #     yCoupling,
        #     verbose   = False,
        #     )
        # combined.color = 1
        # combined.name = 'regularBR'
        # combined.title = 'BR constant'
        # containers.append( combined )

        scalingBR_rootfiles = glob( '{0}/*.root'.format( LatestPaths.scan_combined_Top_couplingDependentBR_bigRange_asimov ) )
        scalingBR = TheoryCommands.GetTH2FromListOfRootFiles(
            scalingBR_rootfiles,
            xCoupling,
            yCoupling,
            verbose   = False,
            )
        scalingBR.color = 1
        scalingBR.name = 'scalingBR'
        scalingBR.title = 'BR(#kappa_{t}, #kappa_{V})'
        containers.append( scalingBR )

        PlotCommands.BasicMixedContourPlot(
            containers,
            xMin = -0.5,
            xMax = 3.0,
            yMin = -0.25,
            yMax = 0.45,
            xTitle    = titles.get( xCoupling, xCoupling ),
            yTitle    = titles.get( yCoupling, yCoupling ),
            plotname  = 'contours_BRcouplingDependency_bigRange',
            x_SM      = 1.,
            y_SM      = 0.,
            plotIndividualH2s = True,
            )


    #____________________________________________________________________
    if args.checkWSParametrization_Top:

        # ======================================
        # Flags

        plotCombination = True
        # plotCombination = False

        DrawExperimentalBinLines = False
        # DrawExperimentalBinLines = True


        # ======================================
        # Set input

        wsToCheck = LatestPaths.ws_combined_Top
        if args.hgg:
            wsToCheck = LatestPaths.ws_hgg_Top
        if args.hzz:
            wsToCheck = LatestPaths.ws_hzz_Top

        TheoryFileInterface.SetFileFinderDir( LatestPaths.derivedTheoryFiles_Top )


        # ======================================
        # Run

        wsParametrization = WSParametrization( wsToCheck )

        topDerivedTheoryFiles = TheoryFileInterface.FileFinder(
            ct='*', cg='*', muR=1, muF=1, Q=1, filter='cb'
            )

        containers = []
        colorCycle = itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )
        for topDerivedTheoryFile in topDerivedTheoryFiles:
            color = next(colorCycle)

            print topDerivedTheoryFile

            container = TheoryFileInterface.ReadDerivedTheoryFile( topDerivedTheoryFile )
            container.name = 'ct_{0}_cg_{1}'.format(
                Commands.ConvertFloatToStr( container.ct ),
                Commands.ConvertFloatToStr( container.cg ),
                )

            container.Tg_theory = TheoryFileInterface.ReadDerivedTheoryFileToTGraph( topDerivedTheoryFile, name = container.name )
            container.Tg_theory.SetLineWidth(2)
            container.Tg_theory.SetLineColor(color)
            container.Tg_theory.SetMarkerColor(color)
            container.Tg_theory.SetLineStyle(1)
            container.Tg_theory.SetMarkerStyle(8)
            container.Tg_theory.SetMarkerSize(0.8)

            container.Tg_parametrization = wsParametrization.GetOutputContainer(
                ct = container.ct , cg = container.cg,
                returnWhat='theory' ).Tg
            container.Tg_parametrization.SetLineColor(color)
            container.Tg_parametrization.SetMarkerColor(color)
            container.Tg_parametrization.SetLineStyle(1)

            container.Tg_parametrization_expBinning = wsParametrization.GetOutputContainer(
                ct = container.ct , cg = container.cg,
                returnWhat='exp' ).Tg
            container.Tg_parametrization_expBinning.SetLineColor(color)
            container.Tg_parametrization_expBinning.SetMarkerColor(color)
            container.Tg_parametrization_expBinning.SetLineStyle(1)

            containers.append( container )


        # ======================================
        # Additional lines to check

        extraTestCouplings = [
            # { 'ct' : -0.565910458565, 'cg' : 9.62666416168 },
            ]

        extraTestContainers = []
        for testCouplings in extraTestCouplings:
            color = next(colorCycle)

            container = Container(
                name = '_'.join([ '{0}_{1:.2f}'.format(key, value) for key, value in sorted(testCouplings.iteritems()) ])
                )
            container.binBoundaries = containers[0].binBoundaries

            container.Tg_parametrization = wsParametrization.GetOutputContainer(
                ct = testCouplings['ct'] , cg = testCouplings['cg'],
                returnWhat='theory' ).Tg
            container.Tg_parametrization.SetLineColor(color)
            container.Tg_parametrization.SetMarkerColor(color)
            container.Tg_parametrization.SetLineStyle(1)

            container.Tg_parametrization_expBinning = wsParametrization.GetOutputContainer(
                ct = testCouplings['ct'] , cg = testCouplings['cg'],
                returnWhat='exp' ).Tg
            container.Tg_parametrization_expBinning.SetLineColor(color)
            container.Tg_parametrization_expBinning.SetMarkerColor(color)
            container.Tg_parametrization_expBinning.SetLineStyle(1)

            extraTestContainers.append( container )



        # ======================================
        # Make plot


        doBigLegend = False

        c.cd()
        c.Clear()
        if doBigLegend:
            SetCMargins( RightMargin=0.3 )
        else:
            SetCMargins()

        xMinAbs = min([ container.Tg_theory.xMin for container in containers ])
        xMaxAbs = max([ container.Tg_theory.xMax for container in containers ])
        yMinAbs = min([ container.Tg_theory.yMin for container in containers ])
        yMaxAbs = 1.4*max([ container.Tg_theory.yMax for container in containers ])

        xMin = xMinAbs - 0.1*( xMaxAbs - xMinAbs )
        xMax = xMaxAbs + 0.1*( xMaxAbs - xMinAbs )
        yMin = yMinAbs - 0.1*( yMaxAbs - yMinAbs )
        yMax = yMaxAbs + 0.1*( yMaxAbs - yMinAbs )

        base = GetPlotBase(
            xMin = xMin,
            xMax = xMax,
            yMin = yMin,
            yMax = yMax,
            xTitle = 'p_{T} [GeV]', yTitle = '#mu'
            )
        base.Draw('P')


        if doBigLegend:
            leg = ROOT.TLegend(
                1 - 0.3,
                c.GetBottomMargin(),
                1 - 0.02 ,
                1 - c.GetTopMargin() 
                )
        else:
            leg = ROOT.TLegend(
                c.GetLeftMargin(),
                1 - c.GetTopMargin() - 0.25,
                1 - c.GetRightMargin(),
                1 - c.GetTopMargin() 
                )
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetNColumns(3)


        Tg_dotDummy = ROOT.TGraph( 1, array( 'f', [-999.] ), array( 'f', [-999.] ) )
        Tg_dotDummy.SetMarkerStyle(8)
        Tg_dotDummy.SetMarkerSize(0.8)
        Tg_dotDummy.SetMarkerColor(49)
        Tg_dotDummy.SetName( 'dotDummy' )
        ROOT.SetOwnership( Tg_dotDummy, False )
        Tg_dotDummy.Draw('P')

        Tg_lineDummy = ROOT.TGraph( 1, array( 'f', [-999.] ), array( 'f', [-999.] ) )
        Tg_lineDummy.SetLineWidth(2)
        Tg_lineDummy.SetLineColor(49)
        Tg_lineDummy.SetName( 'lineDummy' )
        ROOT.SetOwnership( Tg_lineDummy, False )
        Tg_lineDummy.Draw('L')

        leg.AddEntry( Tg_dotDummy.GetName(),  'Theory calc.', 'P' )
        leg.AddEntry( Tg_lineDummy.GetName(), 'Parametrization', 'L' )


        if plotCombination:
            # Combination scan result
            combinationPOIs = Commands.ListPOIs( LatestPaths.ws_combined_smH )
            combinationscans = PhysicsCommands.GetScanResults(
                combinationPOIs,
                LatestPaths.scan_combined_PTH,
                pattern = 'combinedCard'
                )
            TgCombination = PhysicsCommands.GetTGraphForSpectrum( combinationPOIs, combinationscans, name='Combination' )
            TgCombination.SetLineColor(1)
            TgCombination.SetFillColorAlpha( 1, 0.2 )
            TgCombination.SetName('Combination')

            CorrelationMatrices.ConvertTGraphToLinesAndBoxes(
                TgCombination,
                drawImmediately=True, legendObject=leg, noBoxes=False,
                xMaxExternal=xMax,
                yMinExternal=yMin,
                yMaxExternal=yMax,
                )


        for container in containers:
            container.Tg_theory.Draw('XP')
            container.Tg_parametrization.Draw('XL')

            if DrawExperimentalBinLines:
                CorrelationMatrices.ConvertTGraphToLinesAndBoxes(
                    container.Tg_parametrization_expBinning,
                    drawImmediately=True,
                    # legendObject=leg,
                    noBoxes=True,
                    xMaxExternal=xMax,
                    yMinExternal=yMin,
                    yMaxExternal=yMax,
                    )


            if doBigLegend:
                leg.AddEntry( container.Tg_theory.GetName(), container.name, 'p' )
            else:
                leg.AddEntry(
                    container.Tg_theory.GetName(),
                    '#kappa_{{t}} = {0:.2f}, #kappa_{{g}} = {1:.2f}'.format( container.ct, container.cg ),
                    'lp'
                    )


        for container in extraTestContainers:
            container.Tg_parametrization.Draw('XL')
            if DrawExperimentalBinLines:
                CorrelationMatrices.ConvertTGraphToLinesAndBoxes(
                    container.Tg_parametrization_expBinning,
                    drawImmediately=True, legendObject=leg, noBoxes=True, xMaxExternal=xMax )


        leg.Draw()

        outname = '{0}_parametrizationCheck'.format( basename(wsToCheck) )
        SaveC( outname )


    #____________________________________________________________________
    if args.combinationAndContour_TopCtCb:

        container = Container()

        container.expBinBoundaries    = expBinBoundaries
        container.ws_combination      = LatestPaths.ws_combined_smH
        container.scanDir_combination = LatestPaths.scan_combined_PTH
        container.ws_coupling         = LatestPaths.ws_combined_TopCtCb
        container.scanDir_coupling    = LatestPaths.scan_combined_TopCtCb

        container.xCoupling           = 'ct'
        container.yCoupling           = 'cb'
        container.xCouplingTitle      = '#kappa_{t}'
        container.yCouplingTitle      = '#kappa_{b}'

        container.plotTitle           = 'comparison_TopCtCb'

        container.MaximaOfContour     = True
        container.BestFit             = True

        PlotCommands.PlotParametrizationsOnCombination( container )


    #____________________________________________________________________
    if args.combinationAndContour_Top:

        container = Container()

        container.expBinBoundaries    = expBinBoundaries
        container.ws_combination      = LatestPaths.ws_combined_smH
        container.scanDir_combination = LatestPaths.scan_combined_PTH

        suffix = ''

        container.ws_coupling         = LatestPaths.ws_combined_Top
        # container.scanDir_coupling    = LatestPaths.scan_combined_Top
        container.scanDir_coupling    = LatestPaths.scan_combined_Top_bigRange2

        # container.ws_coupling         = LatestPaths.ws_combined_Top_skippedLastBin
        # container.scanDir_coupling    = LatestPaths.scan_combined_Top_skippedLastBin        
        # suffix += '_skippedLastBin'

        container.xCoupling           = 'ct'
        container.yCoupling           = 'cg'
        container.xCouplingTitle      = '#kappa_{t}'
        container.yCouplingTitle      = '#kappa_{g}'

        container.plotTitle           = 'comparison_Top' + suffix

        # container.MaximaOfContour     = True
        # container.StraightLineToSM      = True
        container.BestFit               = True
        container.OnlyYMaximaOfContour  = True

        container.xSM                 = 1.0
        container.ySM                 = 0.0

        container.ManualPoints = [
            ( lambda xBestfit, xSM: xSM, lambda yBestfit, ySM: ySM ),
            ( lambda xBestfit, xSM: xBestfit + 1.0*(xBestfit-xSM) -0.05 , lambda yBestfit, ySM: yBestfit  + 1.0*(yBestfit-ySM) ),
            ( lambda xBestfit, xSM: xBestfit + 2.0*(xBestfit-xSM) -0.10 , lambda yBestfit, ySM: yBestfit  + 2.0*(yBestfit-ySM) ),
            # 
            # ( lambda xBestfit, xSM: xSM + 0.5*(xBestfit-xSM) , lambda yBestfit, ySM: ySM ),
            # ( lambda xBestfit, xSM: xSM + 0.5*(xBestfit-xSM) , lambda yBestfit, ySM: ySM  + 0.01 ),
            # ( lambda xBestfit, xSM: xBestfit + 0.5 , lambda yBestfit, ySM: yBestfit ),
            # ( lambda xBestfit, xSM: xBestfit - 0.5 , lambda yBestfit, ySM: yBestfit )
            ]


        container.MarkerSize = 3.0

        PlotCommands.PlotParametrizationsOnCombination( container )




########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'