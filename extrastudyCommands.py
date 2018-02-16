#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, itertools, operator, re, argparse, sys, random, numpy
from math import isnan, isinf, sqrt
from os.path import *
from glob import glob
from copy import deepcopy
from array import array

import LatestPaths
import LatestBinning

sys.path.append('src')
import Commands
import PhysicsCommands
import OneOfCommands
import TheoryCommands
import CorrelationMatrices
import MergeHGGWDatacards
import TheoryFileInterface
import OutputInterface
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

    parser.add_argument( '--extrastudyCommands', default=False )
    class CustomAction(argparse.Action):
        def __init__(self, option_strings, dest, **kwargs):
            # Only 1 argument allowed (this is basically the "store_true" action)
            super(CustomAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            # Make sure the parser knows a combineCommand was entered
            setattr( namespace, 'extrastudyCommands', True )
            setattr( namespace, self.dest, True )
            # super(FooAction, self).__call__( parser, namespace, values, option_string=None )
            # setattr( namespace, self.dest, values )

    parser.add_argument( '--inclusive_combineCards',                  action=CustomAction )
    parser.add_argument( '--inclusive_t2ws',                          action=CustomAction )
    parser.add_argument( '--inclusive_scan',                          action=CustomAction )

    parser.add_argument( '--onlyNormalization_Top',                   action=CustomAction )
    parser.add_argument( '--onlyNormalization_Yukawa',                action=CustomAction )

    parser.add_argument( '--FitBR_t2ws',                              action=CustomAction )
    parser.add_argument( '--FitBR_bestfit',                           action=CustomAction )
    parser.add_argument( '--FitBR_scan',                              action=CustomAction )
    parser.add_argument( '--FitBR_scan_crossCheck',                   action=CustomAction )
    parser.add_argument( '--FitBR_plot',                              action=CustomAction )

    parser.add_argument( '--TotalXS_t2ws',                            action=CustomAction )
    parser.add_argument( '--TotalXS_bestfit',                         action=CustomAction )
    parser.add_argument( '--TotalXS_scan',                            action=CustomAction )
    parser.add_argument( '--TotalXS_plot',                            action=CustomAction )

    parser.add_argument( '--chi2fitToCombination_Yukawa',             action=CustomAction )
    parser.add_argument( '--RepeatTheoristFit',                       action=CustomAction )

    parser.add_argument( '--Make2DPlotOfVariableInWS',                action=CustomAction )
    parser.add_argument( '--PlotOfTotalXSInYukawaWS',                 action=CustomAction )
    parser.add_argument( '--PlotOfTotalXS_FromParametrization',       action=CustomAction )
    parser.add_argument( '--PlotBRsInOnePlot',                        action=CustomAction )

    parser.add_argument( '--CorrelationMatrixScaleDependence_Yukawa', action=CustomAction )

    parser.add_argument( '--PlotParametrizationShapes',               action=CustomAction )

    parser.add_argument( '--CreateParametrizationTables',              action=CustomAction )    
    
    parser.add_argument( '--MiniCombineTest',                          action=CustomAction )    
    parser.add_argument( '--HbbInclusionPlots',                        action=CustomAction )


########################################
# Methods
########################################    

def main( args ):

    TheoryCommands.set_plot_dir( 'plots_{0}'.format(datestr) )


    #____________________________________________________________________
    if args.MiniCombineTest:
        import MiniCombine

        ASIMOV = True

        DO_TOP = True
        DO_YUKAWA = not(DO_TOP)

        DO_HIGH_PT = True
        # DO_HIGH_PT = False


        # ======================================
        # 

        if ASIMOV:
            combinationResult = LatestPaths.scan_combined_PTH_xHfixed_asimov
        else:
            combinationResult = LatestPaths.scan_combined_PTH_ggH

        Commands.warning( 'The correlation matrix is not calculated on Asimov!!' )
        correlationMatrixFile = LatestPaths.correlationMatrix_PTH_ggH


        # ======================================
        # Set histograms to take parametrization from

        if DO_YUKAWA:
            couplings = [ 'kappac', 'kappab' ]
            SM = TheoryFileInterface.file_finder(
                kappab=1, kappac=1, muR=1, muF=1, Q=1,
                expectOneFile=True,
                directory = LatestPaths.derivedTheoryFiles_YukawaSummed,
                loadImmediately=True
                )
            derivedTheoryFileContainers = TheoryFileInterface.file_finder(
                kappab='*', kappac='*', muR=1, muF=1, Q=1,
                directory = LatestPaths.derivedTheoryFiles_YukawaSummed,
                loadImmediately=True
                )
            includeLinearTerms = True

            # It's not convenient to have the >400 SM bins for <400 variations
            nBinsVariations = len(derivedTheoryFileContainers[0].binBoundaries)-1
            SM.binBoundaries = SM.binBoundaries[:nBinsVariations+1]
            SM.crosssection  = SM.crosssection[:nBinsVariations]
            SM.ratios        = SM.ratios[:nBinsVariations]

        elif DO_TOP:
            couplings = [ 'ct', 'cg' ]
            if DO_HIGH_PT:
                SM = TheoryFileInterface.file_finder(
                    ct=1, cg=0, cb=1, muR=1, muF=1, Q=1,
                    expectOneFile=True,
                    directory = LatestPaths.derivedTheoryFiles_TopHighPt,
                    loadImmediately=True
                    )
                derivedTheoryFileContainers = TheoryFileInterface.file_finder(
                    ct='*', cg='*', muR=1, muF=1, Q=1, cb=1,
                    directory = LatestPaths.derivedTheoryFiles_TopHighPt,
                    loadImmediately=True
                    )
                derivedTheoryFileContainers = [ d for d in derivedTheoryFileContainers if not (d.ct == 1. and d.cg == 0.) ]
            else:
                SM = TheoryFileInterface.file_finder(
                    ct=1, cg=0, cb=1, muR=1, muF=1, Q=1,
                    expectOneFile=True,
                    directory = LatestPaths.derivedTheoryFiles_Top,
                    loadImmediately=True
                    )
                derivedTheoryFileContainers = TheoryFileInterface.file_finder(
                    ct='*', cg='*', muR=1, muF=1, Q=1, filter='cb',
                    directory = LatestPaths.derivedTheoryFiles_Top,
                    loadImmediately=True
                    )
                # It's not convenient to have the >400 SM bins for <400 variations
                nBinsVariations = len(derivedTheoryFileContainers[0].binBoundaries)-1
                SM.binBoundaries = SM.binBoundaries[:nBinsVariations+1]
                SM.crosssection  = SM.crosssection[:nBinsVariations]
                SM.ratios        = SM.ratios[:nBinsVariations]

            includeLinearTerms = False
            xMin = -6.8
            xMax = 6.8
            yMin = -0.6
            yMax = 0.6
            NPointsX = 200
            NPointsY = 200


        # ======================================
        # 

        miniCombine_allBins = MiniCombine.MiniCombine()
        miniCombine_allBins.lastBinIsOverflow = [ True if DO_TOP else False ]

        Commands.warning( 'Using the scan result of the last bin' )
        miniCombine_allBins.set_combination_result( combinationResult, excludeLastBin=False )

        miniCombine_allBins.set_correlation_matrix(
            correlationMatrixFile
            # forceDiagonalMatrix = True
            )

        miniCombine_allBins.couplings = couplings
        miniCombine_allBins.set_sm(SM)
        miniCombine_allBins.parametrize( derivedTheoryFileContainers, includeLinearTerms )

        miniCombine_allBins.link()
        miniCombine_allBins.build_chi2_function()

        miniCombine_allBins.name = 'MiniCombine_allBins'
        miniCombine_allBins.scan_grid(
            NPointsX, xMin, xMax,
            NPointsY, yMin, yMax,
            )


        # ======================================
        # 

        miniCombine = MiniCombine.MiniCombine()
        miniCombine.lastBinIsOverflow = [ True if DO_TOP else False ]

        Commands.warning( 'Not using the scan result of the last bin' )
        miniCombine.set_combination_result( combinationResult, excludeLastBin=True )

        miniCombine.set_correlation_matrix(
            correlationMatrixFile
            # forceDiagonalMatrix = True
            )

        miniCombine.couplings = couplings
        miniCombine.set_sm(SM)
        miniCombine.parametrize( derivedTheoryFileContainers, includeLinearTerms )


        # ======================================
        # 

        miniCombine.link()

        # Plug in the Hgg sensitivity manually
        # miniCombine.add_mu_constraint(         mu = 1.0, unc = 0.5*(0.98+0.76)/0.57, left = 350., right = 450. )
        # miniCombine.add_overflow_mu_constraint( mu = 1.0, unc = 0.5*(1.54+0.97)/0.43, left = 450. )

        miniCombine.add_asym_mu_constraint( 
            mu = 1.0,
            unc_up = 0.82/0.53, unc_down = -0.66/0.53,
            left = 350, right = 600
            )
        miniCombine.add_asym_mu_overflow_constraint( 
            mu = 1.0,
            unc_up = 3.63/0.38, unc_down = -1.65/0.38,
            left = 600
            )

        miniCombine.build_chi2_function()
        miniCombine.scan_grid(
            NPointsX, xMin, xMax,
            NPointsY, yMin, yMax,
            )


        # ======================================
        # Now add Hbb

        # Hbb covmat:
        #        r1                       r2
        # r1     7.637            0.3554
        # r2     0.3554           6.52


        miniCombine.add_mu_constraint(         mu = 1.0, unc = sqrt(7.637), left = 350., right = 575. )
        miniCombine.add_overflow_mu_constraint( mu = 1.0, unc = sqrt(6.520), left = 575. )

        miniCombine.build_chi2_function()
        miniCombine.scan_grid(
            NPointsX, xMin, xMax,
            NPointsY, yMin, yMax,
            )


        # Also again, using only Hbb
        miniCombine.clear_constraints()

        miniCombine.add_mu_constraint(         mu = 1.0, unc = sqrt(7.637), left = 350., right = 575. )
        miniCombine.add_overflow_mu_constraint( mu = 1.0, unc = sqrt(6.520), left = 575. )
        
        miniCombine.build_chi2_function()
        miniCombine.scan_grid(
            NPointsX, xMin, xMax,
            NPointsY, yMin, yMax,
            )



        # for unc in [ 1.0, 2.0, 2.5, 3.0, 4.0 ]:

        #     print '\n' + '-'*70
        #     print 'Scanning with constraint uncertainty', unc

        #     miniCombine.clear_constraints()
        #     miniCombine.add_mu_constraint(
        #         mu = 1.0, unc = unc, left = 350., right = 600.
        #         )
        #     miniCombine.add_overflow_mu_constraint(
        #         mu = 1.0, unc = unc, left = 600.
        #         )
        #     miniCombine.build_chi2_function()

        #     miniCombine.scan_grid(
        #         NPointsX, xMin, xMax,
        #         NPointsY, yMin, yMax,
        #         )


    #____________________________________________________________________
    if args.HbbInclusionPlots:

        def GetObjectsFromRootFile( rootFile ):
            with Commands.OpenRootFile( rootFile ) as rootFp:
                canvas = rootFp.Get('ctc')
                L = canvas.GetListOfPrimitives()

                ret = {
                    'contours_1sigma' : [],
                    'bestfitPoint'    : None,
                    'H2'              : None,
                    }

                # print '\nPrinting content of', rootFile
                for i in xrange(L.GetEntries()):
                    obj = L.At(i)
                    # print obj.GetName()

                    if 'contour_1sigma' in obj.GetName():
                        ret['contours_1sigma'].append( obj )
                    elif 'bestfitpoint' in obj.GetName() or 'bestfitPoint' in obj.GetName():
                        ret['bestfitPoint'] = obj
                    elif obj.GetName() == 'H2':
                        ret['H2'] = obj
                    else:
                        continue
                    
                    ROOT.SetOwnership( obj, False )    

            if ret['contours_1sigma'] == 0:
                Commands.throw_error( 'No objects with substr \'contour_1sigma\' found in {0}'.format(rootFile) )

            return ret


        containers = []

        # plots_Dec13/MiniCombine_allBins_ct_cg.root
        # plots_Dec13/MiniCombine_ct_cg_350to600_uncup1p55_uncdownm1p25_600to800_uncup9p55_uncdownm4p34.root
        # plots_Dec13/MiniCombine_ct_cg_350to575_unc2p76_GT575_unc2p55.root
        # plots_Dec13/MiniCombine_ct_cg_350to600_uncup1p55_uncdownm1p25_600to800_uncup9p55_uncdownm4p34_350to575_unc2p76_GT575_unc2p55.root

        nominal = Container()
        nominal.name     = 'nominal'
        nominal.title    = 'nominal'
        nominal.rootFile = 'plots_Dec13/MiniCombine_allBins_ct_cg.root'
        nominal.color    = 1
        containers.append(nominal)

        miniCombine_hgg = Container()
        miniCombine_hgg.name     = 'hgg'
        miniCombine_hgg.title    = 'H#rightarrow#gamma#gamma p_{T}^{H} (350, 600, #infty)'
        miniCombine_hgg.rootFile = 'plots_Dec13/MiniCombine_ct_cg_350to600_uncup1p55_uncdownm1p25_600to800_uncup9p55_uncdownm4p34.root'
        miniCombine_hgg.color    = 2
        containers.append(miniCombine_hgg)

        miniCombine_hbb = Container()
        miniCombine_hbb.name     = 'hbb'
        miniCombine_hbb.title    = 'H#rightarrowbb p_{T}^{H} (350, 575, #infty)'
        miniCombine_hbb.rootFile = 'plots_Dec13/MiniCombine_ct_cg_350to575_unc2p76_GT575_unc2p55.root'
        miniCombine_hbb.color    = 4
        containers.append(miniCombine_hbb)

        miniCombine_hgg_hbb = Container()
        miniCombine_hgg_hbb.name     = 'hgghbb'
        miniCombine_hgg_hbb.title    = 'H#rightarrow#gamma#gamma (split bins) and H#rightarrowbb'
        miniCombine_hgg_hbb.rootFile = 'plots_Dec13/MiniCombine_ct_cg_350to600_uncup1p55_uncdownm1p25_600to800_uncup9p55_uncdownm4p34_350to575_unc2p76_GT575_unc2p55.root'
        miniCombine_hgg_hbb.color    = 8
        containers.append(miniCombine_hgg_hbb)

        Persistence = []
        for container in containers:
            objs = GetObjectsFromRootFile( container.rootFile )
            container.contours_1sigma = objs['contours_1sigma']
            container.contours_2sigma = []
            container.bestfitPoint    = objs['bestfitPoint']
            container.H2              = objs['H2']
            Persistence.append(objs)

            print '\n1sigma contour extrema of {0}'.format( container.name )
            for Tg in container.contours_1sigma:
                xs, ys = TheoryCommands.get_xyfrom_TGraph(Tg)
                print ' ', Tg.GetName()
                print '    xMin =', min(xs)
                print '    xMax =', max(xs)
                print '    yMin =', min(ys)
                print '    yMax =', max(ys)


        PlotCommands.basic_mixed_contour_plot(
            [ nominal, miniCombine_hgg, miniCombine_hbb ],
            xMin = -5.0,
            xMax = 5.0,
            yMin = -0.40,
            yMax = 0.40,
            # 
            xTitle    = '#kappa_{t}',
            yTitle    = '#kappa_{g}',
            plotname  = 'contours_MiniCombine_1',
            x_SM      = 1.,
            y_SM      = 0.,
            plotIndividualH2s = False,
            filterContours = False,
            only1sigmaContours = True,
            nLegendColumns = 3
            )

        PlotCommands.basic_mixed_contour_plot(
            [ nominal, miniCombine_hgg_hbb ],
            xMin = -5.0,
            xMax = 5.0,
            yMin = -0.40,
            yMax = 0.40,
            # 
            xTitle    = '#kappa_{t}',
            yTitle    = '#kappa_{g}',
            plotname  = 'contours_MiniCombine_2',
            x_SM      = 1.,
            y_SM      = 0.,
            plotIndividualH2s = False,
            filterContours = False,
            only1sigmaContours = True,
            nLegendColumns = 2
            )

        miniCombine_hgg_hbb.color = 1
        PlotCommands.plot_single_TH2(
            miniCombine_hgg_hbb,
            xMin = -5.0,
            xMax = 5.0,
            yMin = -0.40,
            yMax = 0.40,
            xTitle    = '#kappa_{t}',
            yTitle    = '#kappa_{g}',
            plotname = 'MiniCombine_0_hgg_hbb',
            )


    #____________________________________________________________________
    if args.inclusive_combineCards:
        Commands.basic_combine_cards(
            'suppliedInput/combinedCard_smH_{0}_INCLUSIVE.txt'.format(datestr),
            'hgg=' + LatestPaths.card_hgg_INC_unprocessed,
            'hzz=' + LatestPaths.card_hzz_INC_unprocessed
            )

    #____________________________________________________________________
    if args.inclusive_t2ws:

        extraOptions = [
            '--PO verbose=2',
            '--PO \'higgsMassRange=123,127\'',
            ]

        Commands.basic_t2ws_with_model(
            LatestPaths.card_combined_INC,
            pathToModel = 'FitBRModel.py',
            modelName   = 'fitTotalXSModel',
            suffix       = 'fitTotalXS',
            extraOptions = extraOptions,
            )

    #____________________________________________________________________
    if args.inclusive_scan:

        # STATONLY = True
        STATONLY = False
        
        ws = LatestPaths.ws_combined_totalXS

        ASIMOV = False
        if args.asimov: ASIMOV = True

        totalXS_ranges = [ 0., 2. ]

        jobDirectory = 'Scan_TotalXS_{0}'.format( datestr )
        if STATONLY: jobDirectory += '_statonly'
        if ASIMOV: jobDirectory += '_asimov'
        jobDirectory = Commands.make_unique_directory( jobDirectory )

        nPoints = 55
        nPointsPerJob = 5
        queue = 'short.q'

        extraOptions  = [
            # '--importanceSampling={0}:couplingScan'.format( abspath('scanTH2D_Jun01.root') ),
            '-P r',
            '--setPhysicsModelParameterRanges r={0},{1}'.format( totalXS_ranges[0], totalXS_ranges[1] ),
            # '--setPhysicsModelParameters {0}'.format(
            #     ','.join([ '{0}=1.0'.format(i) for i in Commands.list_set( datacard, 'POI' ) if i.startswith('r_') ])
            #     ),
            # '--setPhysicsModelParameterRanges kappab={0},{1}:kappac={2},{3}'.format(
            #     kappab_ranges[0], kappab_ranges[1], kappac_ranges[0], kappac_ranges[1] ),
            # '--saveSpecifiedFunc {0}'.format(','.join(
            #     Commands.list_set( datacard, 'yieldParameters' ) + [ i for i in Commands.list_set( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
            '--squareDistPoiStep',
            ]

        if STATONLY:
            extraOptions.append('--freezeNuisances rgx{.*}')

        Commands.multidim_combine_tool(
            ws,
            nPoints       = nPoints,
            nPointsPerJob = nPointsPerJob,
            queue         = queue,
            notOnBatch    = False,
            jobDirectory  = jobDirectory,
            fastscan      = False,
            asimov        = ASIMOV,
            jobPriority   = 0,
            extraOptions  = extraOptions
            )


    #____________________________________________________________________
    if args.onlyNormalization_Top:

        USE_HIGHPT_PARAMETRIZATION = True
        
        expBinBoundaries = [ 0., 15., 30., 45., 85., 125., 200., 350., 10000. ]

        extraOptions = [
            '--PO verbose=2',
            '--PO \'higgsMassRange=123,127\'',
            '--PO linearTerms=False',
            '--PO splitggH=True',
            ]

        if args.hzz:
            extraOptions.append( '--PO isOnlyHZZ=True' )
        if args.hgg:
            extraOptions.append( '--PO isOnlyHgg=True' )


        if USE_HIGHPT_PARAMETRIZATION:
            TheoryFileInterface.set_file_finder_dir( LatestPaths.derivedTheoryFiles_TopHighPt )
            Commands.warning( 'Correlation matrix was made using low pt spectra!!' )
        else:
            TheoryFileInterface.set_file_finder_dir( LatestPaths.derivedTheoryFiles_Top )

        extraOptions.append(
            '--PO SM=[ct=1,cg=0,file={0}]'.format(
                TheoryFileInterface.file_finder( ct=1, cg=0, cb=1, muR=1, muF=1, Q=1, expectOneFile=True )
                )
            )

        if USE_HIGHPT_PARAMETRIZATION:
            theoryFiles = TheoryFileInterface.file_finder(
                ct='*', cg='*', cb=1, muR=1, muF=1, Q=1, filter='ct_1_cg_0',
                )            
        else:
            theoryFiles = TheoryFileInterface.file_finder(
                ct='*', cg='*', muR=1, muF=1, Q=1, filter='cb'
                )

        possibleTheories = []
        for theoryFile in theoryFiles:
            ct = Commands.convert_str_to_float( re.search( r'ct_([\dmp]+)', theoryFile ).group(1) )
            cg = Commands.convert_str_to_float( re.search( r'cg_([\dmp]+)', theoryFile ).group(1) )
            possibleTheories.append(
                '--PO theory=[ct={0},cg={1},file={2}]'.format(
                    ct, cg, theoryFile
                    )                
                )
        extraOptions.extend(possibleTheories)

        suffix = 'Top'
        extraOptions.append( '--PO ProfileTotalXS=True' )
        suffix += '_profiledTotalXS'
        extraOptions.append( '--PO FitOnlyNormalization=True' )
        suffix += '_fitOnlyNormalization'

        extraOptions.append(
            '--PO binBoundaries={0}'.format( ','.join([ str(b) for b in expBinBoundaries ]) )
            )

        Commands.basic_t2ws_with_model(
            LatestPaths.card_combined_INC,
            'CouplingModel.py',
            suffix = suffix,
            extraOptions = extraOptions,
            )


    #____________________________________________________________________
    if args.onlyNormalization_Yukawa:

        expBinBoundaries = [ 0., 15., 30., 45., 85., 125. ]

        TheoryFileInterface.set_file_finder_dir( LatestPaths.derivedTheoryFiles_YukawaSummed )

        extraOptions = [
            '--PO verbose=2',
            '--PO \'higgsMassRange=123,127\'',
            '--PO linearTerms=True',
            '--PO splitggH=True',
            ]

        extraOptions.append(
            '--PO binBoundaries={0}'.format( ','.join([ str(b) for b in expBinBoundaries ]) )
            )

        extraOptions.append(
            '--PO SM=[kappab=1,kappac=1,file={0}]'.format(
                TheoryFileInterface.file_finder( kappab=1, kappac=1, muR=1, muF=1, Q=1, expectOneFile=True )
                )
            )

        possibleTheories = []
        for kappab in [ -2, -1, 0, 1, 2 ]:
            for kappac in [ -10, -5, 0, 1, 5, 10 ]:
                if ( kappab == 1 and kappac == 1 ) or ( kappab == 0 and kappac == 0 ): continue
                else:
                    possibleTheories.append(
                        '--PO theory=[kappab={0},kappac={1},file={2}]'.format(
                            kappab, kappac,
                            TheoryFileInterface.file_finder( kappab=kappab, kappac=kappac, muR=1, muF=1, Q=1, expectOneFile=True )
                            )
                        )
        # Sample only a few needed theories
        import random
        random.seed(1002)
        extraOptions.extend( random.sample( possibleTheories, 6 ) )


        suffix = 'Yukawa'

        extraOptions.append( '--PO ProfileTotalXS=True' )
        suffix += '_profiledTotalXS'

        extraOptions.append( '--PO FitOnlyNormalization=True' )
        suffix += '_fitOnlyNormalization'


        Commands.basic_t2ws_with_model(
            LatestPaths.card_combined_INC,
            'CouplingModel.py',
            suffix = suffix,
            extraOptions = extraOptions,
            )



    #____________________________________________________________________
    if args.FitBR_t2ws:

        extraOptions = [
            '--PO verbose=2',
            '--PO \'higgsMassRange=123,127\'',
            ]

        # USE_GLOBAL_SCALES = True
        USE_GLOBAL_SCALES = False

        if USE_GLOBAL_SCALES:
            Commands.basic_t2ws_with_model(
                LatestPaths.card_combined_INC,
                'FitBRModel.py',
                extraOptions = extraOptions,
                modelName    = 'fitGlobalBRModel',
                suffix       = 'globalScales',
                )

        else:

            Commands.basic_t2ws_with_model(
                LatestPaths.card_combined_smH_PTH,
                'FitBRModel.py',
                extraOptions = extraOptions,
                smartMaps    = [
                    ( r'.*/smH_PTH_([\d\_GT]+)', r'r_smH_PTH_\1[1.0,-1.0,4.0]' )
                    ],
                )

            Commands.basic_t2ws_with_model(
                LatestPaths.card_combined_smH_NJ,
                'FitBRModel.py',
                extraOptions = extraOptions,
                smartMaps    = [
                    ( r'.*/smH_NJ_([\d\_GTLEpm]+)', r'r_smH_NJ_\1[1.0,-1.0,4.0]' )
                    ],
                )

            Commands.basic_t2ws_with_model(
                LatestPaths.card_combined_smH_YH,
                'FitBRModel.py',
                extraOptions = extraOptions,
                smartMaps    = [
                    ( r'.*/smH_YH_([\d\_GTLEpm]+)', r'r_smH_YH_\1[1.0,-1.0,4.0]' )
                    ],
                )

            Commands.basic_t2ws_with_model(
                LatestPaths.card_combined_smH_PTJ,
                'FitBRModel.py',
                extraOptions = extraOptions,
                smartMaps    = [
                    ( r'.*/smH_PTJ_([\d\_GT]+)', r'r_smH_PTJ_\1[1.0,-1.0,4.0]' )
                    ],
                )


    #____________________________________________________________________
    if args.FitBR_bestfit:

        # USE_GLOBAL_SCALES = True
        USE_GLOBAL_SCALES = False

        if USE_GLOBAL_SCALES:
            ws = abspath( LatestPaths.ws_FitBR_combined_unsplit )
        else:
            ws = abspath( LatestPaths.ws_FitBR_combined_unsplit )

        cmd = [
            'combine',
            ws,
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            # '--algo=grid',
            '--floatOtherPOIs=1',
            '--saveNLL',
            '--saveInactivePOI 1',
            # '--fastScan',
            # '-P kappab',
            # '-P kappac',
            # '--setPhysicsModelParameters kappab=1.0,kappac=1.0',
            # '--saveSpecifiedFunc {0}'.format(','.join(
            #     Commands.list_set( datacard, 'yieldParameters' ) + [ i for i in Commands.list_set( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
            '--squareDistPoiStep',
            '-M MultiDimFit',
            '-m 125.00',
            # '--setPhysicsModelParameterRanges kappab=-20.0,20.0:kappac=-50.0,50.0',
            # '--setPhysicsModelParameterRanges kappab=0.5,1.0:kappac=1.0,2.0',
            # '--points 12800',
            # '--firstPoint 0',
            # '--lastPoint 799',
            '-n testjob',
            # '-v 3',
            ]

        Commands.basic_generic_combine_command(
            cmd,
            onBatch = False,
            )



    if args.FitBR_scan:

        doFastscan = False
        if args.fastscan: doFastscan = True

        ASIMOV = False
        if args.asimov: ASIMOV = True
        
        # USE_GLOBAL_SCALES = True
        USE_GLOBAL_SCALES = False

        jobDirectory = 'Scan_ratioOfBRs_{0}'.format( datestr )


        if USE_GLOBAL_SCALES:
            datacard = abspath( LatestPaths.ws_ratio_of_BRs_globalScales )
            jobDirectory += '_globalScales'
        else:
            datacard = abspath( LatestPaths.ws_ratio_of_BRs )

        ratio_BR_hgg_hzz_ranges = [ 0.03, 0.16 ]


        if ASIMOV: jobDirectory += '_asimov'
        jobDirectory = Commands.make_unique_directory( jobDirectory )

        if doFastscan:
            nPoints = 42
            nPointsPerJob = 42
            queue = 'short.q'
        else:
            nPoints = 42
            nPointsPerJob = 3
            queue = 'short.q'

        Commands.multidim_combine_tool(
            datacard,
            nPoints       = nPoints,
            nPointsPerJob = nPointsPerJob,
            queue         = queue,
            notOnBatch    = False,
            jobDirectory  = jobDirectory,
            fastscan      = doFastscan,
            asimov        = ASIMOV,
            jobPriority   = 0,
            extraOptions  = [
                # '--importanceSampling={0}:couplingScan'.format( abspath('scanTH2D_Jun01.root') ),
                '-P ratio_BR_hgg_hzz',
                '--setPhysicsModelParameterRanges ratio_BR_hgg_hzz={0},{1}'.format( ratio_BR_hgg_hzz_ranges[0], ratio_BR_hgg_hzz_ranges[1] ),
                '--setPhysicsModelParameters {0}'.format(
                    ','.join([ '{0}=1.0'.format(i) for i in Commands.list_set( datacard, 'POI' ) if i.startswith('r_') ])
                    ),
                # '--setPhysicsModelParameterRanges kappab={0},{1}:kappac={2},{3}'.format(
                #     kappab_ranges[0], kappab_ranges[1], kappac_ranges[0], kappac_ranges[1] ),
                # '--saveSpecifiedFunc {0}'.format(','.join(
                #     Commands.list_set( datacard, 'yieldParameters' ) + [ i for i in Commands.list_set( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
                '--squareDistPoiStep',
                ]
            )


    if args.FitBR_scan_crossCheck:

        datacards = [
            LatestPaths.ws_ratio_of_BRs_NJ,
            LatestPaths.ws_ratio_of_BRs_PTJ,
            LatestPaths.ws_ratio_of_BRs_YH,
            LatestPaths.ws_ratio_of_BRs_PTH,
            LatestPaths.ws_ratio_of_BRs_globalScales,
            ]

        ASIMOV = ( True if args.asimov else False )

        ratio_BR_hgg_hzz_ranges = [ 0.03, 0.16 ]

        jobDirectory = 'out/Scan_ratioOfBRs_{0}'.format( datestr )
        if ASIMOV: jobDirectory += '_asimov'

        nPoints = 42
        nPointsPerJob = 3
        queue = 'short.q'

        for datacard in datacards:

            actualJobDirectory = Commands.make_unique_directory( jobDirectory )

            extraOptions = [
                '-P ratio_BR_hgg_hzz',
                '--setPhysicsModelParameterRanges ratio_BR_hgg_hzz={0},{1}'.format( ratio_BR_hgg_hzz_ranges[0], ratio_BR_hgg_hzz_ranges[1] ),
                # '--saveSpecifiedFunc {0}'.format(','.join(
                #     Commands.list_set( datacard, 'yieldParameters' ) + [ i for i in Commands.list_set( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
                '--squareDistPoiStep',
                ]

            yieldParameterSets = [ '{0}=1.0'.format(i) for i in Commands.list_set( datacard, 'POI' ) if i.startswith('r_') ]
            if len(yieldParameterSets) > 0:
                extraOptions.append( '--setPhysicsModelParameters {0}'.format( ','.join(yieldParameterSets) ) )

            Commands.multidim_combine_tool(
                datacard,
                nPoints       = nPoints,
                nPointsPerJob = nPointsPerJob,
                queue         = queue,
                notOnBatch    = False,
                jobDirectory  = actualJobDirectory,
                fastscan      = False,
                asimov        = ASIMOV,
                jobPriority   = 0,
                extraOptions  = extraOptions
                )


    if args.FitBR_plot:

        if args.asimov:
            scanDirs = [
                LatestPaths.scan_ratioOfBRs_NJ_asimov,
                LatestPaths.scan_ratioOfBRs_PTJ_asimov,
                LatestPaths.scan_ratioOfBRs_YH_asimov,
                LatestPaths.scan_ratioOfBRs_PTH_asimov,
                LatestPaths.scan_ratioOfBRs_INC_asimov,
                ]
        else:
            scanDirs = [
                LatestPaths.scan_ratioOfBRs_NJ,
                LatestPaths.scan_ratioOfBRs_PTJ,
                LatestPaths.scan_ratioOfBRs_YH,
                LatestPaths.scan_ratioOfBRs_PTH,
                LatestPaths.scan_ratioOfBRs_INC,
                ]

        varNames = [
            'NJ',
            'PTJ',
            'YH',
            'PTH',
            'INC',
            ]

        for scanDir, varName in zip( scanDirs, varNames ):

            # if USE_GLOBAL_SCALES:
            #     scanRootFiles = glob( LatestPaths.scan_ratioOfBRs_globalScales + '/*.root' )
            # else:
            #     scanRootFiles = glob( LatestPaths.scan_ratioOfBRs + '/*.root' )
            
            scanRootFiles = glob( scanDir + '/*.root' )
            
            scanContainer = OutputInterface.OutputContainer()

            x_unfiltered, y_unfiltered = TheoryCommands.basic_read_scan(
                scanRootFiles,
                xAttr = 'ratio_BR_hgg_hzz',
                yAttr = 'deltaNLL',
                )

            scanContainer.x = []
            scanContainer.y = []
            for x, y in zip( x_unfiltered, y_unfiltered ):
                if (
                    x > -10.0 and x < 10.
                    and
                    y > -10.0 and y < 10.
                    ):
                    scanContainer.x.append( x )
                    scanContainer.y.append( y )


            # Do uncertainty determination before scaling
            # FindMinimaAndErrors expects deltaNLL, not 2*deltaNLL
            scanContainer.extrema = PhysicsCommands.find_minima_and_errors( scanContainer.x, scanContainer.y, returnContainer=True )

            print '[info] Multiplying by 2: deltaNLL --> chi^2'
            scanContainer.y = [ 2.*y for y in scanContainer.y ]

            scanContainer.get_TGraph( xAttr = 'x', yAttr = 'y', xAreBinBoundaries = False )
            scanContainer.Tg.SetMarkerStyle(8)
            scanContainer.Tg.SetMarkerSize(0.8)


            # ======================================
            # Make plot

            c.Clear()
            set_cmargins( TopMargin=0.09 )

            yMinAbs = min( scanContainer.y )
            yMaxAbs = max( scanContainer.y )
            # yMin = yMinAbs - 0.1*(yMaxAbs-yMinAbs)
            yMin = 0.0
            # yMax = yMaxAbs + 0.5*(yMaxAbs-yMinAbs)
            yMax = 5.0

            # xMin = min( scanContainer.x )
            # xMax = max( scanContainer.x )
            xMin = 0.04
            xMax = 0.15

            base = get_plot_base(
                xMin = xMin,
                xMax = xMax,
                yMin = yMin,
                yMax = yMax,
                xTitle = 'BR_{H #rightarrow #gamma#gamma} / BR_{H #rightarrow ZZ}',
                # yTitle = '#Delta NLL',
                yTitle = '2#DeltaNLL',
                )
            base.Draw('P')


            scanContainer.Tg.Draw('XPL')


            oneLine = ROOT.TLine( xMin, 1.0, xMax, 1.0 )
            oneLine.SetLineColor(12)
            oneLine.Draw()


            print '\n' + '-'*70
            print 'Found range: {0:.4f} < ratio < {1:.4f}'.format( scanContainer.extrema.leftBound, scanContainer.extrema.rightBound )

            # imin
            # min
            # leftError
            # leftBound
            # rightError
            # rightBound
            # wellDefinedRightBound
            # wellDefinedLeftBound

            # Check number more carefully
            SM_ratio = LatestBinning.SM_ratio_of_BRs
            SMLine = ROOT.TLine( SM_ratio, yMin, SM_ratio, yMax )
            SMLine.SetLineWidth(2)
            SMLine.SetLineColor(9)
            SMLine.Draw()


            xBestfit = scanContainer.x[ scanContainer.extrema.imin ]
            bestfitLine = ROOT.TLine( xBestfit, yMin, xBestfit, yMax )
            bestfitLine.SetLineWidth(2)
            bestfitLine.SetLineColor(2)
            bestfitLine.Draw()


            l = ROOT.TLatex()
            l.SetNDC()
            l.SetTextColor(1)
            l.SetTextSize(0.05)

            # l.SetTextAlign(31)
            # l.DrawLatex(
            #     scanContainer.extrema.leftBound, 1.0 + 0.013*(yMax-yMin),
            #     '-{0:.2f} ({1:d}%)'.format(
            #         abs(scanContainer.extrema.leftError),
            #         int( abs(scanContainer.extrema.leftError) / xBestfit * 100. )
            #         )
            #     )
            # l.SetTextAlign(11)
            # l.DrawLatex(
            #     scanContainer.extrema.rightBound, 1.0 + 0.013*(yMax-yMin),
            #     '+{0:.2f} ({1:d}%)'.format(
            #         abs(scanContainer.extrema.rightError),
            #         int( abs(scanContainer.extrema.rightError) / xBestfit * 100. )
            #         )
            #     )
            # l.SetTextAlign(21)
            # l.DrawLatex(
            #     xBestfit, 1.0 + 0.013*(yMax-yMin),
            #     '{0:.3f}'.format( xBestfit )
            #     )

            l.DrawLatex( 0.55, 0.8, 'R = {0:.3f} _{{{1:+.3f}}}^{{{2:+.3f}}}'.format(
                xBestfit, -abs(scanContainer.extrema.leftError), scanContainer.extrema.rightError
                ))
            l.DrawLatex( 0.55, 0.6, varName )


            TgPoints = ROOT.TGraph( 2,
                array( 'f', [ scanContainer.extrema.leftBound, scanContainer.extrema.rightBound ] ),
                array( 'f', [ 1.0, 1.0 ] ),
                )
            TgPoints.SetMarkerSize(1.2)
            TgPoints.SetMarkerStyle(8)
            TgPoints.SetMarkerColor(2)
            TgPoints.Draw('PSAME')

            Commands.get_cms_label()
            Commands.get_cms_lumi()

            # save_c( 'BRscan' + ( '_globalScales' if USE_GLOBAL_SCALES else '' ) )
            save_c( 'BRscan_' + basename(scanDir).replace('/','') )




    #____________________________________________________________________
    if args.TotalXS_t2ws:

        extraOptions = [
            '--PO verbose=2',
            '--PO \'higgsMassRange=123,127\'',
            ]

        Commands.basic_t2ws_with_model(
            LatestPaths.card_combined_smH_PTH,
            pathToModel = 'FitBRModel.py',
            modelName   = 'fitTotalXSModel',
            suffix       = 'fitTotalXS',
            extraOptions = extraOptions
            )

    #____________________________________________________________________
    if args.TotalXS_bestfit:

        ws = LatestPaths.ws_totalXS

        cmd = [
            'combine',
            ws,
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            # '--algo=grid',
            '--floatOtherPOIs=1',
            '--saveNLL',
            '--saveInactivePOI 1',
            '-P r',
            # '--fastScan',
            # '-P kappab',
            # '-P kappac',
            # '--setPhysicsModelParameters kappab=1.0,kappac=1.0',
            # '--saveSpecifiedFunc {0}'.format(','.join(
            #     Commands.list_set( datacard, 'yieldParameters' ) + [ i for i in Commands.list_set( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
            '--squareDistPoiStep',
            '-M MultiDimFit',
            '-m 125.00',
            # '--setPhysicsModelParameterRanges kappab=-20.0,20.0:kappac=-50.0,50.0',
            # '--setPhysicsModelParameterRanges kappab=0.5,1.0:kappac=1.0,2.0',
            # '--points 12800',
            # '--firstPoint 0',
            # '--lastPoint 799',
            '-n testjob',
            # '-v 3',
            ]

        Commands.basic_generic_combine_command(
            cmd,
            onBatch = False,
            )


    #____________________________________________________________________
    if args.TotalXS_scan:

        # STATONLY = True
        STATONLY = False

        ws = LatestPaths.ws_totalXS

        doFastscan = False
        if args.fastscan: doFastscan = True

        ASIMOV = False
        if args.asimov: ASIMOV = True

        totalXS_ranges = [ 0., 2. ]

        jobDirectory = 'Scan_TotalXS_{0}'.format( datestr )
        if STATONLY: jobDirectory += '_statonly'
        if ASIMOV: jobDirectory += '_asimov'
        jobDirectory = Commands.make_unique_directory( jobDirectory )

        if doFastscan:
            nPoints = 42
            nPointsPerJob = 42
            queue = 'short.q'
            queue = '8nm'
        else:
            nPoints = 50
            nPointsPerJob = 5
            queue = 'short.q'

        extraOptions  = [
            # '--importanceSampling={0}:couplingScan'.format( abspath('scanTH2D_Jun01.root') ),
            '-P r',
            '--setPhysicsModelParameterRanges r={0},{1}'.format( totalXS_ranges[0], totalXS_ranges[1] ),
            # '--setPhysicsModelParameters {0}'.format(
            #     ','.join([ '{0}=1.0'.format(i) for i in Commands.list_set( datacard, 'POI' ) if i.startswith('r_') ])
            #     ),
            # '--setPhysicsModelParameterRanges kappab={0},{1}:kappac={2},{3}'.format(
            #     kappab_ranges[0], kappab_ranges[1], kappac_ranges[0], kappac_ranges[1] ),
            # '--saveSpecifiedFunc {0}'.format(','.join(
            #     Commands.list_set( datacard, 'yieldParameters' ) + [ i for i in Commands.list_set( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
            '--squareDistPoiStep',
            ]

        if STATONLY:
            extraOptions.append('--freezeNuisances rgx{.*}')

        Commands.multidim_combine_tool(
            ws,
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



    def filter_scan( xs, ys, xMin=-10., xMax=10., yMin=-10., yMax=10.):
        xs_filtered = []
        ys_filtered = []
        for x, y in zip(xs, ys):
            if (
                x > xMin and x < xMax
                and
                y > yMin and y < yMax
                ):
                xs_filtered.append(x)
                ys_filtered.append(y)
        return xs_filtered, ys_filtered


    if args.TotalXS_plot:

        yMaxObjects = 3.75

        def make_container(scanDir):

            x_unfiltered, y_unfiltered = TheoryCommands.basic_read_scan(
                Commands.glob_root_files(scanDir),
                xAttr = 'r',
                yAttr = 'deltaNLL',
                )
            if len(x_unfiltered)==0:
                raise ValueError('x_unfiltered is an empty list, scanDir = {0}'.format(scanDir))
            if len(y_unfiltered)==0:
                raise ValueError('y_unfiltered is an empty list, scanDir = {0}'.format(scanDir))

            scanContainer = OutputInterface.OutputContainer()
            scanContainer.x, scanContainer.y = filter_scan(x_unfiltered, y_unfiltered, yMax=yMaxObjects*.5)

            print '[info] Multiplying all x by {0}'.format( LatestBinning.YR4_totalXS )
            scanContainer.x = [ LatestBinning.YR4_totalXS*x for x in scanContainer.x ]

            # FindMinimaAndErrors assumes deltaNLL (not 2*deltaNLL), so compute it before scaling
            scanContainer.extrema = PhysicsCommands.find_minima_and_errors( scanContainer.x, scanContainer.y, returnContainer=True )

            print '[info] Multiplying by 2: deltaNLL --> chi^2'
            scanContainer.y = [ 2.*y for y in scanContainer.y ]

            scanContainer.get_TGraph( xAttr = 'x', yAttr = 'y', xAreBinBoundaries = False )
            scanContainer.Tg.SetMarkerStyle(8)
            scanContainer.Tg.SetMarkerSize(0.8)

            scanContainer.Tg.SetName(TheoryCommands.get_unique_root_name())

            return scanContainer


        scanDir_statsyst = 'out/Scan_TotalXS_Feb09_0'
        scanContainer = make_container(scanDir_statsyst)

        scanDir_statonly = 'out/Scan_TotalXS_Feb09_0_statonly'
        scanContainer_statonly = make_container(scanDir_statonly)
        scanContainer_statonly.Tg.SetLineColor(14)
        scanContainer_statonly.Tg.SetMarkerColor(14)


        # ======================================
        # Make plot

        c.Clear()
        set_cmargins( TopMargin = 0.09 )

        yMinAbs = min( scanContainer.y )
        yMaxAbs = max( scanContainer.y )
        # yMin = yMinAbs - 0.1*(yMaxAbs-yMinAbs)
        yMin = 0.0
        # yMax = yMaxAbs + 0.5*(yMaxAbs-yMinAbs)
        yMax = 5.0

        # xMin = min( scanContainer.x )
        # xMax = max( scanContainer.x )
        xMin = 0.8 * LatestBinning.YR4_totalXS
        xMax = 1.4 * LatestBinning.YR4_totalXS

        base = get_plot_base(
            xMin = xMin,
            xMax = xMax,
            yMin = yMin,
            yMax = yMax,
            xTitle = '#sigma_{tot} (pb)',
            # yTitle = '#Delta NLL',
            yTitle = '2#DeltaNLL',
            )
        base.Draw('P')


        scanContainer.Tg.Draw('XC')
        scanContainer_statonly.Tg.Draw('XC')


        oneLine = ROOT.TLine( xMin, 1.0, xMax, 1.0 )
        oneLine.SetLineColor(12)
        oneLine.Draw()

        print '\n' + '-'*70
        print 'Found range: {0:.4f} < r_totalXS < {1:.4f}'.format( scanContainer.extrema.leftBound, scanContainer.extrema.rightBound )

        xBestfit = scanContainer.x[ scanContainer.extrema.imin ]
        left_full  = scanContainer.extrema.leftError
        right_full = scanContainer.extrema.rightError
        left_stat  = scanContainer_statonly.extrema.leftError
        right_stat = scanContainer_statonly.extrema.rightError
        left_syst  = -sqrt(left_full**2 - left_stat**2)
        right_syst = sqrt(right_full**2 - right_stat**2)
        print '\nSeparated syst-stat:'
        print 'r_totalXS = {0:.4f}  {1:.4f}/{2:.4f} (stat) {3:.4f}/{4:.4f} (syst)'.format(
            xBestfit, right_stat, left_stat, right_syst, left_syst)

        # TLines can't be added to the legend
        # bestfitLine = ROOT.TLine( xBestfit, yMin, xBestfit, yMaxObjects )
        bestfitLine = ROOT.TGraph(2, array('f', [xBestfit, xBestfit]), array('f', [yMin, yMaxObjects]))
        bestfitLine.SetName('bestfitLine')
        bestfitLine.SetLineWidth(2)
        bestfitLine.SetLineColor(2)
        bestfitLine.Draw('LSAME')

        # for x in [ scanContainer.extrema.leftBound, scanContainer.extrema.rightBound ]:
        #     uncLine = ROOT.TLine( x, 0.0, x, 1.0 )
        #     ROOT.SetOwnership( uncLine, False )
        #     uncLine.SetLineWidth(1)
        #     uncLine.SetLineColor(2)
        #     uncLine.Draw()

        smLine = ROOT.TGraph(2, array('f', [LatestBinning.YR4_totalXS, LatestBinning.YR4_totalXS]), array('f', [yMin, yMaxObjects]))
        smLine.SetName('smLine')
        smLine.SetLineWidth(2)
        smLine.SetLineColor(9)
        smLine.Draw('LSAME')


        l = ROOT.TLatex()
        l.SetNDC()
        l.SetTextColor(1)
        l.SetTextSize(0.040)

        # l.SetTextAlign(31)
        # l.DrawLatex(
        #     scanContainer.extrema.leftBound, 1.0 + 0.013*(yMax-yMin),
        #     '-{0:.2f} ({1:d}%)'.format(
        #         abs(scanContainer.extrema.leftError),
        #         int( abs(scanContainer.extrema.leftError) / xBestfit * 100. )
        #         )
        #     )
        # l.SetTextAlign(11)
        # l.DrawLatex(
        #     scanContainer.extrema.rightBound, 1.0 + 0.013*(yMax-yMin),
        #     '+{0:.2f} ({1:d}%)'.format(
        #         abs(scanContainer.extrema.rightError),
        #         int( abs(scanContainer.extrema.rightError) / xBestfit * 100. )
        #         )
        #     )
        # l.SetTextAlign(21)
        # l.DrawLatex(
        #     xBestfit, 1.0 + 0.013*(yMax-yMin),
        #     '{0:.3f}'.format( xBestfit )
        #     )


        # l.DrawLatex( 0.65, 0.8, '#sigma_{{tot}} = {0:.1f} _{{{1:+.1f}}}^{{{2:+.1f}}}'.format(
        #     xBestfit, -abs(scanContainer.extrema.leftError), scanContainer.extrema.rightError
        #     ))
        # l.DrawLatex( c.GetLeftMargin()+0.03, 0.78, '#sigma_{{tot}} = {0:.1f} _{{{1:+.1f} (stat.)  {2:+.1f} (syst.)}}^{{{3:+.1f} (stat.)  {4:+.1f} (syst.)}} pb'.format(
        #     xBestfit, -abs(left_stat), -abs(left_syst), right_stat, right_syst
        #     ))
        l.DrawLatex( c.GetLeftMargin()+0.03, 0.78, '#sigma_{{tot}} = {0:.1f}  #pm{1:.1f} (stat.) #pm{2:.1f} (syst.)  pb'.format(
            xBestfit, 0.5*(abs(left_stat)+abs(right_stat)), 0.5*(abs(left_syst)+abs(right_syst))
            ))


        TgPoints = ROOT.TGraph( 2,
            array( 'f', [ scanContainer.extrema.leftBound, scanContainer.extrema.rightBound ] ),
            array( 'f', [ 1.0, 1.0 ] ),
            )
        TgPoints.SetMarkerSize(1.2)
        TgPoints.SetMarkerStyle(8)
        TgPoints.SetMarkerColor(2)
        TgPoints.Draw('PSAME')

        leg = ROOT.TLegend(
            c.GetLeftMargin()+0.02,
            1-c.GetTopMargin()-0.07, 
            c.GetLeftMargin()+0.7,
            1-c.GetTopMargin(),
            )
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.set_n_columns(4)

        leg.AddEntry(scanContainer.Tg.GetName(), 'Stat.+syst.', 'l')
        leg.AddEntry(scanContainer_statonly.Tg.GetName(), 'Stat. only', 'l')
        leg.AddEntry(bestfitLine.GetName(), 'Best fit', 'l')
        leg.AddEntry(smLine.GetName(), 'SM', 'l')

        leg.Draw()

        Commands.get_cms_label()
        Commands.get_cms_lumi()

        save_c( 'totalXSscan' )


    if args.chi2fitToCombination_Yukawa:
        TheoryCommands.set_plot_dir( 'plots_{0}_Yukawa'.format(datestr) )

        # Workspace to get the parametrization from
        # ws = LatestPaths.ws_combined_split_betterYukawa

        # Derived theory files to get the parametrization from

        # File to get the correlation matrix from
        # corrMatFile = 'corrMat_Oct17/higgsCombine_CORRMAT_combinedCard_Jul26.MultiDimFit.mH125.root'
        # corrMatFile = 'corrMat_Oct17/higgsCombine_CORRMAT_combinedCard_Aug21.MultiDimFit.mH125.root'
        # corrMatFile = 'corrMat_Oct19/higgsCombine_CORRMAT_combinedCard_Aug21_xHfixed.MultiDimFit.mH125.root'
        corrMatFile = 'corrMat_Nov08_combinedCard_Nov03_xHfixed/higgsCombine_CORRMAT_combinedCard_Nov03_xHfixed.MultiDimFit.mH125.root'

        # Scan to get the uncertainties from
        # scanDir = LatestPaths.scan_ptcombination_combined_profiled
        # scanDir = LatestPaths.scan_ptcombination_combined_profiled_xHfixed
        scanDir = LatestPaths.scan_combined_PTH_xHfixed


        expBinning = [ 0., 15., 30., 45., 85., 125. ]
        yieldParameterNames = []
        for left, right in zip( expBinning[:-1], expBinning[1:] ):
            yieldParameterNames.append( 'r_ggH_PTH_{0}_{1}'.format( int(left), int(right) ) )
        nBinsExp = len(expBinning) - 1


        # # ======================================
        # # Load parametrization from WS

        # rootFp = ROOT.TFile.Open(ws)
        # w = rootFp.Get('w')

        # yieldParameterNames = Commands.list_set( ws, 'yieldParameters' )
        # yieldParameters = [ w.function(yPName) for yPName in yieldParameterNames ]
        # kappacRealVar = w.var('kappac')
        # kappabRealVar = w.var('kappab')

        # binBoundaries   = [ w.var(binBoundName).getVal() for binBoundName in Commands.list_set( ws, 'expBinBoundaries' ) ]
        # nBins = len(binBoundaries)-1


        # ======================================
        # Load parametrization from derived theory files

        SM = TheoryFileInterface.file_finder(
            kappab=1, kappac=1, muR=1, muF=1, Q=1,
            expectOneFile=True,
            directory = LatestPaths.derivedTheoryFiles_YukawaSummed,
            loadImmediately=True
            )
        derivedTheoryFileContainers = TheoryFileInterface.file_finder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFiles_YukawaSummed,
            loadImmediately=True
            )
        parametrization = Parametrization()
        parametrization.set_sm(SM)
        parametrization.parametrize( derivedTheoryFileContainers )

        theoryBinBoundaries = SM.binBoundaries


        # ======================================
        # Get the combination result (center values + uncertainties)

        combinationscans = PhysicsCommands.get_scan_results(
            yieldParameterNames,
            scanDir,
            # pattern = 'combinedCard'
            )
        TgCombination = PhysicsCommands.get_TGraph_for_spectrum( yieldParameterNames, combinationscans, name='Combination' )

        # Get relevant bins and errors from scan
        bestfitCenters  = TgCombination.POICenters[:nBinsExp]
        bestfitDownErrs = TgCombination.POIErrsLeft[:nBinsExp]
        bestfitUpErrs   = TgCombination.POIErrsRight[:nBinsExp]
        bestfitSymmErrs = TgCombination.POIErrsSymm[:nBinsExp]


        # ======================================
        # Get the correlation matrix

        corrMat = []
        corrMatFp = ROOT.TFile.Open(corrMatFile)
        fit = corrMatFp.Get('fit')
        for poi1 in yieldParameterNames:
            corrMatRow = []
            for poi2 in yieldParameterNames:
                corrMatRow.append( fit.correlation( poi1, poi2 ) )
            corrMat.append( corrMatRow )

        print '\nFound corrMat from {0}:'.format(corrMatFile)
        print numpy.array(corrMat)

        # Compute covariance matrix
        covMat = [ [ 0. for i in xrange(nBinsExp) ] for j in xrange(nBinsExp) ]
        for i in xrange(nBinsExp):
            for j in xrange(nBinsExp):
                covMat[i][j] = bestfitSymmErrs[i] * bestfitSymmErrs[j] * corrMat[i][j]

        print '\nComputed covMat:'
        print numpy.array(covMat)

        covMat_npArray = numpy.array(covMat)
        covMat_inversed_npArray = numpy.linalg.inv( covMat_npArray )


        # ======================================
        # Doing a bestfit

        # Build simple chi2 function
        def chi2_function( inputTuple ):
            kappac, kappab = inputTuple

            mus_parametrization = parametrization.evaluate_for_binning(
                theoryBinBoundaries, expBinning,
                kappab = kappab, kappac = kappac,
                returnRatios=True
                )

            x_column = numpy.array( [ mu - 1.0 for mu in mus_parametrization ] ).T

            # x_column = numpy.array( [ y.getVal() for y in yieldParameters ] ).T

            # chi2Val = 0.
            # for i in xrange(nBinsExp):
            #     chi2Val += (yieldParameterValues[i]-bestfitCenters[i])**2 / bestfitSymmErrs[i]**2

            chi2 = x_column.T.dot(  covMat_inversed_npArray.dot( x_column )  )

            return chi2


        # Do a best fit
        from scipy.optimize import minimize
        print ''
        res = minimize( chi2_function, [ 1., 1. ], method='Nelder-Mead', tol=1e-6 )
        print 'End of minimization'
        print res

        chi2_bestfit = res.fun
        kappac_bestfit = res.x[0]
        kappab_bestfit = res.x[1]


        # Do a scan

        kappacMin = -35.
        kappacMax = 35.
        kappabMin = -13.
        kappabMax = 13.

        kappacNPoints = 100
        kappabNPoints = 100


        kappacBinBoundaries = [ kappacMin + i*(kappacMax-kappacMin)/float(kappacNPoints) for i in xrange(kappacNPoints+1) ]
        kappabBinBoundaries = [ kappabMin + i*(kappabMax-kappabMin)/float(kappabNPoints) for i in xrange(kappabNPoints+1) ]

        kappacPoints = [ 0.5*(kappacBinBoundaries[i]+kappacBinBoundaries[i+1]) for i in xrange(kappacNPoints) ]
        kappabPoints = [ 0.5*(kappabBinBoundaries[i]+kappabBinBoundaries[i+1]) for i in xrange(kappabNPoints) ]


        H2 = ROOT.TH2F(
            'H2', '',
            kappacNPoints, array( 'f', kappacBinBoundaries ),
            kappabNPoints, array( 'f', kappabBinBoundaries ),
            )

        for i_kappab, kappabVal in enumerate(kappabPoints):
            for i_kappac, kappacVal in enumerate(kappacPoints):
                H2.SetBinContent( i_kappac, i_kappab, chi2_function( (kappacVal, kappabVal) ) - chi2_bestfit )


        print ''
        contours_1sigma = TheoryCommands.get_contours_from_TH2( H2, 2.30 )
        contours_2sigma = TheoryCommands.get_contours_from_TH2( H2, 6.18 )


        # ======================================
        # Plotting

        H2.SetTitle('')
        H2.GetXaxis().SetTitle( '#kappa_{c}' )
        H2.GetYaxis().SetTitle( '#kappa_{b}' )
        
        H2.SetMaximum(7.0)

        c.Clear()
        set_cmargins(
            LeftMargin   = 0.12,
            RightMargin  = 0.10,
            BottomMargin = 0.12,
            TopMargin    = 0.09,
            )

        H2.Draw('COLZ')

        for Tg in contours_1sigma:
            Tg.SetLineWidth(2)
            Tg.Draw('LSAME')
        for Tg in contours_2sigma:
            Tg.SetLineWidth(2)
            Tg.SetLineStyle(2)
            Tg.Draw('LSAME')


        SMpoint = ROOT.TGraph( 1, array( 'f', [1.0] ), array( 'f', [1.0] ) )
        SMpoint.SetMarkerStyle(21)
        SMpoint.SetMarkerSize(2)
        SMpoint.Draw('PSAME')

        bestfitpoint = ROOT.TGraph( 1, array( 'f', [kappac_bestfit] ), array( 'f', [kappab_bestfit] ) )
        bestfitpoint.SetMarkerStyle(34)
        bestfitpoint.SetMarkerSize(1.1)
        bestfitpoint.SetMarkerColor(13)
        bestfitpoint.Draw('PSAME')
        bestfitpoint.SetName('bestfitpoint')

        save_c( 'AfterTheFactChi2Fit', asROOT = True )

        # rootFp.Close()





    #____________________________________________________________________
    if args.PlotParametrizationShapes:
        newColorCycle = lambda: itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )

        TheoryCommands.set_plot_dir( 'plots_{0}_Yukawa'.format(datestr) )

        # scanDir     = LatestPaths.scan_combined_PTH_xHfixed
        scanDir     = LatestPaths.scan_combined_PTH_xHfixed_asimov
        corrMatFile = 'corrMat_Nov08_combinedCard_Nov03_xHfixed/higgsCombine_CORRMAT_combinedCard_Nov03_xHfixed.MultiDimFit.mH125.root'


        # ======================================
        # Some settings

        expBinBoundaries = [
            0., 15., 30., 45., 85., 125.,
            # 200., 350., 10000.
            ]

        nBins = len(expBinBoundaries)-1
        binWidths = [ expBinBoundaries[i+1] - expBinBoundaries[i] for i in xrange(nBins) ]
        binCenters = [ 0.5*( expBinBoundaries[i] + expBinBoundaries[i+1] ) for i in xrange(nBins) ]


        # ======================================
        # Load parametrization from derived theory files

        SM = TheoryFileInterface.file_finder(
            kappab=1, kappac=1, muR=1, muF=1, Q=1,
            expectOneFile=True,
            directory = LatestPaths.derivedTheoryFiles_YukawaSummed,
            loadImmediately=True
            )
        derivedTheoryFileContainers = TheoryFileInterface.file_finder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFiles_YukawaSummed,
            loadImmediately=True
            )
        parametrization = Parametrization()
        # parametrization.set_sm(SM)

        # parametrization.parametrize( derivedTheoryFileContainers )
        parametrization.parametrize_by_fitting( derivedTheoryFileContainers, fitWithScipy=True )

        theoryBinBoundaries = SM.binBoundaries
        theoryBinWidths  = [ right - left for left, right in zip( theoryBinBoundaries[:-1], theoryBinBoundaries[1:] ) ]
        theoryBinCenters = [ 0.5*( left + right ) for left, right in zip( theoryBinBoundaries[:-1], theoryBinBoundaries[1:] ) ]

        SM_integralFunction = TheoryCommands.get_integral(
            SM.binBoundaries,
            SM.crosssection
            )


        # ======================================
        # Scan for some points

        points = [

            #  ( kappac, kappab )
            # ( 1.0,    1.0 ),
            # ( 10.,    1.0 ),
            # ( -5.,    2.0 ),
            # ( 0.,     200. ),
            # ( 505.0,  505.0 ),
            # ( -495.0, 505. ),
            # ( 0.,     100000. )

            # ( 0., 0. ),
            # ( 1., 1. ),
            # ( 2., 2. ),
            # ( 4., 4. ),
            # ( 10., 10. ),
            # ( 100., 100. ),
            # ( 1000., 1000. ),
            # ( 10000., 10000. ),

            ( 10., 0. ),
            ( 10., 1./10. ),
            ( 10., 1./3. ),
            ( 10., 1./5. ),
            ( 10., 1. ),
            ( 10., 3. ),
            ( 10., 5. ),
            ( 10., 10. ),

            ( 10., -1./10. ),
            ( 10., -1./3. ),
            ( 10., -1./5. ),
            ( 10., -1. ),
            ( 10., -3. ),
            ( 10., -5. ),
            ( 10., -10. ),

            ( -10., 0. ),            

            ]


        plotContainers = []

        for kappac, kappab in points:

            print '\nkappac = {0}, kappab = {1}'.format( kappac, kappab )

            XSs_perGeV_parametrization = parametrization.evaluate_for_binning(
                theoryBinBoundaries, expBinBoundaries,
                kappab = kappab, kappac = kappac,
                returnRatios=False,
                verbose = True
                )


            print 'XSs_perGeV_parametrization: ', XSs_perGeV_parametrization

            for derivedTheoryFileContainer in derivedTheoryFileContainers:
                if derivedTheoryFileContainer.kappab == kappab and derivedTheoryFileContainer.kappac == kappac:

                    integral = TheoryCommands.get_integral(
                        derivedTheoryFileContainer.binBoundaries,
                        derivedTheoryFileContainer.crosssection
                        )

                    exp_xs = []
                    for left, right in zip( expBinBoundaries[:-1], expBinBoundaries[1:] ):
                        exp_xs.append( integral( left, right ) / ( right - left ) )

                    print '  for matching container:   ', exp_xs
                    break


            XSs_parametrization = [ xs * binWidth for xs, binWidth in zip( XSs_perGeV_parametrization, binWidths ) ]

            # Calculate the shape for the exp binning
            S_parametrization = [ xs / sum(XSs_parametrization) for xs in XSs_parametrization ]

            print 'XSs_parametrization:        ', XSs_parametrization
            print 'S_parametrization:          ', S_parametrization


            # Also calculate shape for theory binning
            XSs_perGeV_parametrization_theoryBinning = parametrization.evaluate(
                kappab = kappab, kappac = kappac,
                returnRatios=False,
                # verbose = True
                )

            XSs_parametrization_theoryBinning = [
                xs * binWidth for xs, binWidth in zip( XSs_perGeV_parametrization_theoryBinning, theoryBinWidths )
                ]

            S_parametrization_theoryBinning = [ xs / sum(XSs_parametrization_theoryBinning) for xs in XSs_parametrization_theoryBinning ]


            container = Container()

            container.kappab = kappab
            container.kappac = kappac
            container.S_parametrization = S_parametrization
            container.S_parametrization_theoryBinning = S_parametrization_theoryBinning

            plotContainers.append( container )



        # ======================================
        # Make plot

        c.Clear()
        set_cmargins()

        xMin = expBinBoundaries[0]
        xMax = expBinBoundaries[-1]

        yMin = 0.0
        yMax = 0.5


        base = get_plot_base(
            xMin = xMin,
            xMax = xMax,
            yMin = yMin,
            yMax = yMax,
            xTitle = 'p_{T} [GeV]', yTitle = 'Shape [A.U.]'
            )
        base.Draw('P')


        leg = ROOT.TLegend(
            c.GetLeftMargin(),
            1 - c.GetTopMargin() - 0.25,
            1 - c.GetRightMargin(),
            1 - c.GetTopMargin() 
            )
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.set_n_columns(3)


        colorCycle = new_color_cycle()
        for container in plotContainers:
            color = colorCycle.next()

            # First the theory
            Tg_theory = ROOT.TGraph(
                len(theoryBinCenters),
                array( 'f', theoryBinCenters ),
                array( 'f', container.S_parametrization_theoryBinning )
                )
            ROOT.SetOwnership( Tg_theory, False )

            Tg_theory.SetMarkerColor(color)
            Tg_theory.SetMarkerStyle(8)
            Tg_theory.SetMarkerSize(0.9)
            Tg_theory.Draw('SAMEP')

            Tg_theory.SetName( TheoryCommands.get_unique_root_name() )


            # H_exp = ROOT.TH1F(
            #     TheoryCommands.get_unique_root_name(), '',

            #     )


            leg.AddEntry(
                Tg_theory.GetName(),
                '#kappa_{{c}} = {0:d}, #kappa_{{b}} = {1:d}'.format( int(container.kappac), int(container.kappab) ),
                'p'
                )

        leg.Draw()

        save_c( 'ParametrizationShapes' )



    #____________________________________________________________________
    if args.RepeatTheoristFit:
        TheoryCommands.set_plot_dir( 'plots_{0}_Yukawa'.format(datestr) )

        # scanDir     = LatestPaths.scan_combined_PTH_xHfixed
        scanDir     = LatestPaths.scan_combined_PTH_xHfixed_asimov
        corrMatFile = 'corrMat_Nov08_combinedCard_Nov03_xHfixed/higgsCombine_CORRMAT_combinedCard_Nov03_xHfixed.MultiDimFit.mH125.root'


        NORMALIZE_BY_SMXS = True
        # NORMALIZE_BY_SMXS = False

        DO_TOP            = True
        # DO_TOP            = False

        DO_YUKAWA = not(DO_TOP)

        PUT_CONSTRAINT_ON_LAST_BIN = True
        # PUT_CONSTRAINT_ON_LAST_BIN = False

        # VARYING_HBB_UNCERTAINTY = True
        VARYING_HBB_UNCERTAINTY = False

        DO_HIGH_PT = True
        # DO_HIGH_PT = False


        # ======================================
        # Load parametrization from derived theory files

        if VARYING_HBB_UNCERTAINTY and not PUT_CONSTRAINT_ON_LAST_BIN:
            Commands.throw_error( 'Don\'t put VARYING_HBB_UNCERTAINTY to True and PUT_CONSTRAINT_ON_LAST_BIN to False' )

        if DO_YUKAWA:
            # Set up exp binning used for the fit
            expBinning = [ 0., 15., 30., 45., 85., 125. ]

            SM = TheoryFileInterface.file_finder(
                kappab=1, kappac=1, muR=1, muF=1, Q=1,
                expectOneFile=True,
                directory = LatestPaths.derivedTheoryFiles_YukawaSummed,
                loadImmediately=True
                )
            derivedTheoryFileContainers = TheoryFileInterface.file_finder(
                kappab='*', kappac='*', muR=1, muF=1, Q=1,
                directory = LatestPaths.derivedTheoryFiles_YukawaSummed,
                loadImmediately=True
                )

            lastBinIsOverflow = False

        elif DO_TOP:

            # Set up exp binning used for the fit
            expBinning = [ 0., 15., 30., 45., 85., 125., 200., 350., 10000. ]

            if DO_HIGH_PT:

                SM = TheoryFileInterface.file_finder(
                    ct=1, cg=0, cb=1, muR=1, muF=1, Q=1,
                    expectOneFile=True,
                    directory = LatestPaths.derivedTheoryFiles_TopHighPt,
                    loadImmediately=True
                    )
                derivedTheoryFileContainers = TheoryFileInterface.file_finder(
                    ct='*', cg='*', muR=1, muF=1, Q=1, cb=1,
                    directory = LatestPaths.derivedTheoryFiles_TopHighPt,
                    loadImmediately=True
                    )
                derivedTheoryFileContainers = [ d for d in derivedTheoryFileContainers if not (d.ct == 1. and d.cg == 0.) ]

            else:
                SM = TheoryFileInterface.file_finder(
                    ct=1, cg=0, cb=1, muR=1, muF=1, Q=1,
                    expectOneFile=True,
                    directory = LatestPaths.derivedTheoryFiles_Top,
                    loadImmediately=True
                    )
                derivedTheoryFileContainers = TheoryFileInterface.file_finder(
                    ct='*', cg='*', muR=1, muF=1, Q=1, filter='cb',
                    directory = LatestPaths.derivedTheoryFiles_Top,
                    loadImmediately=True
                    )

            lastBinIsOverflow = True


        yieldParameterNames = []
        for left, right in zip( expBinning[:-1], expBinning[1:] ):
            if lastBinIsOverflow and right == expBinning[-1]:
                yieldParameterNames.append( 'r_ggH_PTH_GT{0}'.format( int(left) ) )
            else:
                yieldParameterNames.append( 'r_ggH_PTH_{0}_{1}'.format( int(left), int(right) ) )
        nBinsExp = len(expBinning) - 1


        parametrization = Parametrization()
        parametrization.set_sm(SM)
        parametrization.parametrize( derivedTheoryFileContainers,
            includeLinearTerms = False,
            couplingsToParametrize = [ 'ct', 'cg' ]
            )

        theoryBinBoundaries = SM.binBoundaries

        SM_integralFunction = TheoryCommands.get_integral(
            SM.binBoundaries,
            SM.crosssection
            )


        # ======================================
        # Calculate SM cross section

        allBinBoundaries = [
            0., 15., 30., 45., 85., 125.,
            200., 350., 10000.
            ]

        nBins = len(allBinBoundaries)-1
        binWidths = [ allBinBoundaries[i+1] - allBinBoundaries[i] for i in xrange(nBins) ]

        # shape = [
        #     0.208025,
        #     0.234770,
        #     0.165146,
        #     0.218345,
        #     0.087552,
        #     0.059154,
        #     0.022612,
        #     0.004398
        #     ]
        # # shape = [ s / sum(shape) for s in shape ]
        # SMXSs = [ s * LatestBinning.YR4_totalXS for s in shape ]


        # Re-normalize the shape
        # shape = shape[:nBinsExp]
        # shape = [ s / sum(shape) for s in shape ]

        SMXStot = SM_integralFunction( expBinning[0], expBinning[-1] )
        # SMXSs = [ s * SMXStot for s in shape ]

        # Actually, just use the SM from the parametrization, to at least get expected results back
        SMXSs = []
        for i in xrange(nBinsExp):
            SMXSs.append(
                SM_integralFunction( expBinning[i], expBinning[i+1] )
                )

        print '\nDetermined SM:'
        print '  SMXStot    = ', SMXStot
        # print '  shape      = ', shape
        print '  SMXSs      = ', SMXSs

        if DO_HIGH_PT:
            SMXS_350_600 = SM_integralFunction( 350., 600. )
            SMXS_GT600   = SM_integralFunction( 600., 10000. )
            print 'SMXS_350_600 = ', SMXS_350_600
            print 'SMXS_GT600   = ', SMXS_GT600


        # ======================================
        # Get the combination result (center values + uncertainties)

        combinationscans = PhysicsCommands.get_scan_results(
            yieldParameterNames,
            scanDir,
            # pattern = 'combinedCard'
            )
        TgCombination = PhysicsCommands.get_TGraph_for_spectrum( yieldParameterNames, combinationscans, name='Combination' )

        # Get relevant bins and errors from scan
        bestfitCenters  = TgCombination.POICenters[:nBinsExp]
        bestfitDownErrs = TgCombination.POIErrsLeft[:nBinsExp]
        bestfitUpErrs   = TgCombination.POIErrsRight[:nBinsExp]
        bestfitSymmErrs = TgCombination.POIErrsSymm[:nBinsExp]


        XSs_data = [ mu * SMXS for mu, SMXS in zip( bestfitCenters, SMXSs ) ]
        XS_data_tot = sum(XSs_data)

        # Input to theorist fit:
        S_data = [ xs / XS_data_tot for xs in XSs_data ]
        S_data_errs = [ e * s for e, s in zip( bestfitSymmErrs, S_data ) ]


        print '\nRead data from', scanDir
        print '  Found total XS: ', XS_data_tot
        print '  Found mus: ', bestfitCenters
        print '  Found xss: ', XSs_data
        print '  xs/xstot:  ', S_data
        print ''
        print '  Uncertainties:'
        print '  deltaMu:   ', bestfitSymmErrs
        print '  dxs/xstot: ', S_data_errs
        print ''


        # ======================================
        # Get the correlation matrix

        corrMat = []
        corrMatFp = ROOT.TFile.Open(corrMatFile)
        fit = corrMatFp.Get('fit')
        for poi1 in yieldParameterNames:
            corrMatRow = []
            for poi2 in yieldParameterNames:
                corrMatRow.append( fit.correlation( poi1, poi2 ) )
            corrMat.append( corrMatRow )

        print '\nFound corrMat from {0}:'.format(corrMatFile)
        print numpy.array(corrMat)


        print '\nOverwriting with square matrix for testing purposes:'
        for i in xrange(len(yieldParameterNames)):
            for j in xrange(len(yieldParameterNames)):
                if i==j:
                    corrMat[i][j] = 1.
                else:
                    corrMat[i][j] = 0.
        print numpy.array(corrMat)


        # Compute covariance matrix
        covMat = [ [ 0. for i in xrange(nBinsExp) ] for j in xrange(nBinsExp) ]
        for i in xrange(nBinsExp):
            for j in xrange(nBinsExp):
                covMat[i][j] = S_data_errs[i] * S_data_errs[j] * corrMat[i][j]

        print '\nComputed covMat:'
        print numpy.array(covMat)

        covMat_npArray = numpy.array(covMat)
        covMat_inversed_npArray = numpy.linalg.inv( covMat_npArray )


        # ======================================
        # Build simple chi2 function

        VERBOSITY_IN_CHI2 = True
        # VERBOSITY_IN_CHI2 = False

        f = lambda number, width = 6: '{0:+{width}.{decimals}f}'.format( number, width=width, decimals=width-4 )
        def chi2_function( inputTuple ):
            kappac, kappab = inputTuple

            if DO_YUKAWA:
                XSs_perGeV_parametrization = parametrization.evaluate_for_binning(
                    theoryBinBoundaries, expBinning,
                    kappab = kappab, kappac = kappac,
                    returnRatios=False,
                    verbose = True if VERBOSITY_IN_CHI2 else False
                    )
            elif DO_TOP:
                XSs_perGeV_parametrization = parametrization.evaluate_for_binning(
                    theoryBinBoundaries, expBinning,
                    ct=kappac, cg=kappab,
                    returnRatios=False,
                    verbose = True if VERBOSITY_IN_CHI2 else False
                    )

                if DO_HIGH_PT:
                    XSs_perGeV_parametrization_600boundary = parametrization.evaluate_for_binning(
                        theoryBinBoundaries, expBinning[:-1] + [ 600., 10000. ],
                        ct=kappac, cg=kappab,
                        returnRatios=False,
                        verbose = True if VERBOSITY_IN_CHI2 else False
                        )
                    mu_350_600 = XSs_perGeV_parametrization_600boundary[-2] / SMXS_350_600
                    mu_GT600   = XSs_perGeV_parametrization_600boundary[-1] / SMXS_GT600


            XSs_parametrization = [ xs * binWidth for xs, binWidth in zip( XSs_perGeV_parametrization, binWidths ) ]

            # Calculate the shape for the exp binning
            if NORMALIZE_BY_SMXS:
                S_parametrization = [ xs / SMXStot for xs in XSs_parametrization ]
            else:
                S_parametrization = [ xs / sum(XSs_parametrization) for xs in XSs_parametrization ]

            # Calculate the column
            x_column = numpy.array( [ s_parametrization - s_data for s_parametrization, s_data in zip( S_parametrization, S_data ) ] ).T
            chi2 = x_column.T.dot(  covMat_inversed_npArray.dot( x_column )  )


            if PUT_CONSTRAINT_ON_LAST_BIN:
                mu_lastbin = XSs_parametrization[-1] / XSs_data[-1]

                if VARYING_HBB_UNCERTAINTY:
                    deltaq = abs( mu_lastbin - 1.0 ) / hbb_uncertainty
                    # print '    mu_lastbin = {0:.2f}; Adding deltaq = {1:.3f} to chi2 (hbb unc = {2:.2f})'.format( mu_lastbin, deltaq, hbb_uncertainty )
                else:
                    if mu_lastbin > 1.0:
                        deltaq = abs( mu_lastbin - 1.0 ) / 2.5
                    else:
                        deltaq = abs( mu_lastbin - 1.0 ) / 2.3

                chi2 += deltaq


            if VERBOSITY_IN_CHI2:
                print ''
                if DO_YUKAWA:
                    print '    kappac = {0}, kappab = {1}'.format( f(kappac), f(kappab) )
                elif DO_TOP:
                    print '    ct = {0}, cg = {1}'.format( f(kappac), f(kappab) )
                print '    XSs_parametrization: ' + ' | '.join([ f(number) for number in XSs_parametrization ])
                print '    Shape_param:         ' + ' | '.join([ f(number) for number in S_parametrization ])
                print '    Shape_data:          ' + ' | '.join([ f(number) for number in S_data ])
                print '    mu_350_600 =', mu_350_600
                print '    mu_GT600   =', mu_GT600
                print '    Chi2: {0:.8f}'.format(chi2)

            return chi2


        # ======================================
        # Loop - can be only 1 iteration

        if VARYING_HBB_UNCERTAINTY:
            hbb_uncertainties = [ 0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 2.4, 3.0 ]
        else:
            hbb_uncertainties = [ 999 ]

        from scipy.optimize import minimize
        for hbb_uncertainty in hbb_uncertainties:


            # ======================================
            # Do a best fit
            
            print ''
            res = minimize( chi2_function, [ 1., 1. ], method='Nelder-Mead', tol=1e-6 )
            print 'End of minimization'
            print res

            chi2_bestfit = res.fun
            kappac_bestfit = res.x[0]
            kappab_bestfit = res.x[1]


            # ======================================
            # Do a scan

            # kappac <--> ct, kappab <--> cg


            if DO_YUKAWA:
                if NORMALIZE_BY_SMXS:
                    kappacMin = -45.
                    kappacMax = 45.
                    kappabMin = -23.
                    kappabMax = 23.
                else:
                    # kappacMin = -1000.
                    # kappacMax = 1000.
                    # kappabMin = -1000.
                    # kappabMax = 1000.

                    kappacMin = -100000.
                    kappacMax = 100000.
                    kappabMin = -100000.
                    kappabMax = 100000.

            elif DO_TOP:
                kappacMin = -6.8
                kappacMax = 6.8
                kappabMin = -0.6
                kappabMax = 0.6


            kappacNPoints       = 200
            kappabNPoints       = 200
            kappacBinBoundaries = [ kappacMin + i*(kappacMax-kappacMin)/float(kappacNPoints) for i in xrange(kappacNPoints+1) ]
            kappabBinBoundaries = [ kappabMin + i*(kappabMax-kappabMin)/float(kappabNPoints) for i in xrange(kappabNPoints+1) ]
            kappacPoints        = [ 0.5*(kappacBinBoundaries[i]+kappacBinBoundaries[i+1]) for i in xrange(kappacNPoints) ]
            kappabPoints        = [ 0.5*(kappabBinBoundaries[i]+kappabBinBoundaries[i+1]) for i in xrange(kappabNPoints) ]

            H2 = ROOT.TH2F(
                'H2', '',
                kappacNPoints, array( 'f', kappacBinBoundaries ),
                kappabNPoints, array( 'f', kappabBinBoundaries ),
                )

            print '\n\nDoing scan'

            # VERBOSITY_IN_CHI2 = True
            VERBOSITY_IN_CHI2 = False
            iIteration = 0

            for i_kappab, kappabVal in enumerate(kappabPoints):
                for i_kappac, kappacVal in enumerate(kappacPoints):

                    if i_kappab % 50 == 0 and i_kappac % 50 == 0:
                        # print '  ', i_kappab, i_kappac
                        VERBOSITY_IN_CHI2 = True

                    H2.SetBinContent( i_kappac, i_kappab, chi2_function( (kappacVal, kappabVal) ) - chi2_bestfit )

                    VERBOSITY_IN_CHI2 = False


            print ''
            contours_1sigma = TheoryCommands.get_contours_from_TH2( H2, 2.30 )
            contours_2sigma = TheoryCommands.get_contours_from_TH2( H2, 6.18 )


            # ======================================
            # Plotting

            H2.SetTitle('')
            H2.GetXaxis().SetTitle( '#kappa_{c}' )
            H2.GetYaxis().SetTitle( '#kappa_{b}' )
            H2.SetMaximum(7.0)

            c.Clear()
            set_cmargins(
                LeftMargin   = 0.12,
                RightMargin  = 0.10,
                BottomMargin = 0.12,
                TopMargin    = 0.09,
                )
            H2.Draw('COLZ')

            for Tg in contours_1sigma:
                Tg.SetLineWidth(2)
                Tg.Draw('LSAME')
            for Tg in contours_2sigma:
                Tg.SetLineWidth(2)
                Tg.SetLineStyle(2)
                Tg.Draw('LSAME')

            SMpoint = ROOT.TGraph( 1, array( 'f', [1.0] ), array( 'f', [1.0] ) )
            SMpoint.SetMarkerStyle(21)
            SMpoint.SetMarkerSize(2)
            SMpoint.Draw('PSAME')

            bestfitpoint = ROOT.TGraph( 1, array( 'f', [kappac_bestfit] ), array( 'f', [kappab_bestfit] ) )
            bestfitpoint.SetMarkerStyle(34)
            bestfitpoint.SetMarkerSize(1.1)
            bestfitpoint.SetMarkerColor(13)
            bestfitpoint.Draw('PSAME')
            bestfitpoint.SetName('bestfitpoint')

            suffix = ''
            if DO_TOP:
                suffix += '_Top'
            elif DO_YUKAWA:
                suffix += '_Yukawa'
            if PUT_CONSTRAINT_ON_LAST_BIN:
                suffix += '_lastBinConstrained'
            if VARYING_HBB_UNCERTAINTY:
                suffix += '_hbbUnc_{0}'.format( Commands.convert_float_to_str(hbb_uncertainty) )

            save_c(
                'TheoristChi2Fit' + suffix,
                asROOT = True, asPNG=True
                )


    #____________________________________________________________________
    if args.Make2DPlotOfVariableInWS:
        TheoryCommands.set_plot_dir( 'plots_{0}_Yukawa'.format(datestr) )


        # ws = LatestPaths.ws_combined_split_betterYukawa_couplingDependentBR
        ws = LatestPaths.ws_combined_yukawa_couplingDependentBR

        # varToPlot = 'hggBRmodifier'
        # varToPlot = 'hzzBRmodifier'
        varToPlot = 'c7_BRscal_hgg'
        # varToPlot = 'c7_BRscal_hzz'
        # varToPlot = 'Scaling_hgluglu'


        xVar = 'kappac'
        xRange = [ -32., 32. ]
        xNBins = 200

        yVar = 'kappab'
        yRange = [ -15., 15. ]
        yNBins = 200


        # ======================================
        # 

        rootFp = ROOT.TFile.Open( ws )
        w = rootFp.Get('w')

        x = w.var(xVar)
        y = w.var(yVar)

        mH = w.var('MH')
        mH.setVal( 125. )

        for getter in [ 'var', 'function' ]:
            z = getattr( w, getter )(varToPlot)
            try:
                z.getVal()
                break
            except ReferenceError:
                pass
            except TypeError:
                pass

        else:
            print 'Variable \'{0}\' is probably not in the workspace'.format(varToPlot)
            sys.exit()



        def makeBinBoundaries( xMin, xMax, nBins ):
            return [ xMin + i*(xMax-xMin)/nBins for i in xrange(nBins+1) ]

        xBinBoundaries = makeBinBoundaries( xRange[0], xRange[1], xNBins )
        yBinBoundaries = makeBinBoundaries( yRange[0], yRange[1], yNBins )

        xBinCenters = [ 0.5*(xBinBoundaries[i]+xBinBoundaries[i+1]) for i in xrange(xNBins) ]
        yBinCenters = [ 0.5*(yBinBoundaries[i]+yBinBoundaries[i+1]) for i in xrange(yNBins) ]


        H2 = ROOT.TH2F(
            'H2', '',
            xNBins, array( 'f', xBinBoundaries ),
            yNBins, array( 'f', yBinBoundaries ),
            )

        H2.SetTitle(varToPlot)
        H2.GetXaxis().SetTitle( xVar )
        H2.GetYaxis().SetTitle( yVar )

        H2.SetMaximum(1.5)

        for ix in xrange(xNBins):
            for iy in xrange(yNBins):

                xCenter = xBinCenters[ix]
                yCenter = yBinCenters[iy]

                x.setVal(xCenter)
                y.setVal(yCenter)

                zValue = z.getVal()

                H2.SetBinContent( ix, iy, zValue )



        c.Clear()
        set_cmargins(
            LeftMargin   = 0.12,
            RightMargin  = 0.10,
            BottomMargin = 0.12,
            TopMargin    = 0.09,
            )

        H2.Draw('COLZ')

        save_c( '2Dplot_{0}'.format(varToPlot) )


        rootFp.Close()




    #____________________________________________________________________
    if args.PlotBRsInOnePlot:
        TheoryCommands.set_plot_dir( 'plots_{0}_Yukawa'.format(datestr) )

        PLOT_SCALING = False

        LOGSCALE = True
        # LOGSCALE = False


        # ======================================
        # 

        # ws = LatestPaths.ws_combined_split_betterYukawa_couplingDependentBR_profiledTotalXS
        # ws = LatestPaths.ws_combined_split_betterYukawa_couplingDependentBR
        # ws = LatestPaths.ws_combined_yukawa_couplingDependentBR
        ws = LatestPaths.ws_combined_Yukawa_couplingDependentBR

        rootFp = ROOT.TFile.Open( ws )
        w = rootFp.Get('w')

        mH = w.var('MH')
        mH.setVal(125.)


        titleDict = {
            'kappab' : '#kappa_{b}',
            'kappac' : '#kappa_{c}',
            'kappa_V' : '#kappa_{V}',
            }


        SMBR = {
            'hww'     : 2.137E-01,
            'hzz'     : 2.619E-02,
            'htt'     : 6.272E-02,
            'hmm'     : 2.176E-04,
            'hbb'     : 5.824E-01,
            'hcc'     : 2.891E-02,
            'hgg'     : 2.270E-03,
            'hzg'     : 1.533E-03,
            'hgluglu' : 8.187E-02,
            }

        SMcouplings = {
            'kappab' : 1.0,
            'kappac' : 1.0,
            'kappa_V' : 1.0,
            }
        couplings = SMcouplings.keys()

        BRfunctions = {
            'hww'     : 'c7_BRscal_hww',
            'hzz'     : 'c7_BRscal_hzz',
            'htt'     : 'c7_BRscal_htt',
            'hmm'     : 'c7_BRscal_hmm',
            'hbb'     : 'c7_BRscal_hbb',
            'hcc'     : 'c7_BRscal_hcc',
            'hgg'     : 'c7_BRscal_hgg',
            'hzg'     : 'c7_BRscal_hzg',
            'hgluglu' : 'c7_BRscal_hgluglu',
            }

        decaysToPlot = BRfunctions.keys()


        for xVar in couplings:

            x = w.var(xVar)

            # Set other couplings to their SM value
            for otherCoupling in couplings:
                if otherCoupling == xVar: continue
                w.var(otherCoupling).setVal( SMcouplings[otherCoupling] )

            if xVar == 'kappab':
                kappab_min = -10
                kappab_max = 10
            elif xVar == 'kappac':
                kappab_min = -30
                kappab_max = 30
            elif xVar == 'kappa_V':
                kappab_min = -15
                kappab_max = 15

            # ======================================
            # From here on, 'kappab' is the xVar

            kappab = x

            nPoints = 100
            kappab_axis = [ kappab_min + i*(kappab_max-kappab_min)/(nPoints-1.) for i in xrange(nPoints) ]

            BRvalues = { d : [] for d in decaysToPlot }
            for decay in decaysToPlot:
                y = w.function( BRfunctions[decay] )
                for kappab_val in kappab_axis:
                    kappab.setVal(kappab_val)
                    BRvalues[decay].append(y.getVal())


            # ======================================
            # Plotting

            c.Clear()
            set_cmargins()

            yMin = -0.1
            yMax = 5.0 if PLOT_SCALING else 1.05

            if LOGSCALE:
                yMin = 0.0001
                c.SetLogy(True)

            base = get_plot_base(
                xMin = kappab_min,
                xMax = kappab_max,
                yMin = yMin,
                yMax = yMax,
                xTitle = titleDict.get( xVar, xVar ),
                yTitle = 'BR scaling' if PLOT_SCALING else 'BR',
                )
            base.Draw('P')

            leg = ROOT.TLegend(
                1 - c.GetRightMargin() - 0.35,
                1 - c.GetTopMargin() - 0.50,
                1 - c.GetRightMargin(),
                1 - c.GetTopMargin() 
                )
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)


            colorCycle = itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )
            for decay in decaysToPlot:
                color = next(colorCycle)

                if PLOT_SCALING:
                    yValues = BRvalues[decay]
                else:
                    yValues = [ SMBR[decay] * val for val in BRvalues[decay] ]

                Tg = ROOT.TGraph(
                    nPoints,
                    array( 'd', kappab_axis ),
                    array( 'd', yValues ),
                    )
                ROOT.SetOwnership( Tg, False )
                Tg.SetLineWidth(2)
                Tg.SetLineColor(color)
                Tg.Draw('LSAME')
                Tg.SetName(decay)

                leg.AddEntry( Tg.GetName(), Tg.GetName(), 'l' )

            if not PLOT_SCALING:

                sumY = [ 0. for i in xrange(nPoints) ]
                for i in xrange(nPoints):
                    for d in decaysToPlot:
                        sumY[i] +=  SMBR[d] * BRvalues[d][i]

                sumTg = ROOT.TGraph(
                    nPoints,
                    array( 'd', kappab_axis ),
                    array( 'd', sumY ),
                    )
                ROOT.SetOwnership( sumTg, False )
                sumTg.SetLineWidth(4)
                sumTg.SetLineStyle(2)
                sumTg.SetLineColor(1)
                sumTg.Draw('LSAME')


            lineAtOne = ROOT.TLine( 1.0, yMin, 1.0, yMax )
            lineAtOne.Draw('L')

            lineAtKappabOne = ROOT.TLine( kappab_min, 1.0, kappab_max, 1.0 )
            lineAtKappabOne.Draw('L')

            leg.Draw()

            save_c(
                ( 'BRscalings' if PLOT_SCALING else 'BRs' )
                + ( '_logscale' if LOGSCALE else '' )
                + '-' + xVar
                )
            rootFp.Close()





    #____________________________________________________________________
    if args.PlotOfTotalXSInYukawaWS:
        TheoryCommands.set_plot_dir( 'plots_{0}_Yukawa'.format(datestr) )

        # ws = LatestPaths.ws_combined_split_betterYukawa_couplingDependentBR_profiledTotalXS
        # ws = LatestPaths.ws_combined_split_betterYukawa_couplingDependentBR
        ws = LatestPaths.ws_combined_yukawa_couplingDependentBR

        rootFp = ROOT.TFile.Open( ws )
        w = rootFp.Get('w')


        for yVar in [ 'totalXS', 'Scaling_hgluglu' ]:

            y = w.function(yVar)

            kappab = w.var('kappab')
            kappac = w.var('kappac')
            mH     = w.var('MH')

            mH.setVal(125.)

            y_SM = y.getVal()


            minKappab = -5.
            maxKappab = 15.
            nPoints = 100

            kappabAxis = [ minKappab + i*(maxKappab-minKappab)/(nPoints-1.) for i in xrange(nPoints) ]

            yAxis = []
            for kappabVal in kappabAxis:
                kappab.setVal(kappabVal)

                yVal = y.getVal()

                if yVar == 'totalXS':
                    yAxis.append( yVal / y_SM )
                else:
                    yAxis.append( yVal )


            c.Clear()
            set_cmargins()


            # yMin = min(yAxis) - 0.1*(max(yAxis)-min(yAxis)),
            # yMax = max(yAxis) + 0.1*(max(yAxis)-min(yAxis)),

            yMin = 0.7
            yMax = 2.3

            base = get_plot_base(
                xMin = minKappab,
                xMax = maxKappab,
                yMin = yMin,
                yMax = yMax,
                xTitle = '#kappa_{b}',
                yTitle = yVar,
                )
            base.Draw('P')

            base.SetTitle( '{0}_byEvaluatingFunctionInWorkspace'.format(yVar) )


            Tg = ROOT.TGraph(
                nPoints,
                array( 'd', kappabAxis ),
                array( 'd', yAxis ),
                )
            Tg.SetLineWidth(2)
            Tg.SetLineColor(2)
            Tg.Draw('LSAME')

            lineAtOne = ROOT.TLine( minKappab, 1.0, maxKappab, 1.0 )
            lineAtOne.Draw('L')


            save_c( 'kappabDependency_{0}_byEvaluatingFunctionInWorkspace'.format(yVar) )

            rootFp.Close()


            if True:

                minY = min(yAxis)
                minX = kappabAxis[ yAxis.index(minY) ]

                print 'minimum: kappab = {0}, {1} = {2}'.format( minX, yVar, minY )

                kappab.setVal(0.)
                y_0 = y.getVal() / y_SM

                print 'kappab = 0.0, {0} = {1}'.format( yVar, y_0 )



    #____________________________________________________________________
    if args.PlotOfTotalXS_FromParametrization:
        TheoryCommands.set_plot_dir( 'plots_{0}_Yukawa'.format(datestr) )


        # ======================================
        # Gather all the theory shapes

        readList = lambda theoryFiles: [ TheoryFileInterface.read_derived_theory_file( theoryFile ) for theoryFile in theoryFiles ]

        theoryFiles_kappab_kappac = TheoryFileInterface.file_finder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFilesDirectory_YukawaSummed
            )
        containers_kappab_kappac = readList( theoryFiles_kappab_kappac )


        theoryFiles_kappab_kappac_Gluon = TheoryFileInterface.file_finder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFilesDirectory_YukawaGluonInduced
            )
        containers_kappab_kappac_Gluon = readList( theoryFiles_kappab_kappac_Gluon )


        theoryFiles_kappab_kappac_Quark = TheoryFileInterface.file_finder(
            kappab='*', kappac='*', filter='muR',
            directory = LatestPaths.derivedTheoryFilesDirectory_YukawaQuarkInduced
            )
        containers_kappab_kappac_Quark = readList( theoryFiles_kappab_kappac_Quark )


        theoryFiles_kappab_kappac_QuarkScaled = TheoryFileInterface.file_finder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFilesDirectory_YukawaQuarkInducedScaled
            )
        containers_kappab_kappac_QuarkScaled = readList( theoryFiles_kappab_kappac_QuarkScaled )


        # ======================================
        # 

        for key in [ 'Gluon', 'Quark', 'QuarkScaled', 'Summed' ]:

            if key == 'Gluon':
                containers = containers_kappab_kappac_Gluon
                parametrizeByFitting = False
            elif key == 'Quark':
                containers = containers_kappab_kappac_Quark
                parametrizeByFitting = True
            elif key == 'QuarkScaled':
                containers = containers_kappab_kappac_QuarkScaled
                parametrizeByFitting = False
            elif key == 'Summed':
                containers = containers_kappab_kappac
                parametrizeByFitting = False


            SM = [ container for container in containers if container.kappab == 1 and container.kappac == 1 and container.muR == 1 and container.muF ==1 and container.Q == 1 ][0]

            parametrization = Parametrization()
            if parametrizeByFitting:
                parametrization.parametrize_by_fitting( containers, fitWithScipy=True )
            else:
                parametrization.parametrize( containers )
            parametrization.kappac = 1.0

            minKappab = -5.
            maxKappab = 15.
            nPoints = 100
            kappabAxis = [ minKappab + i*(maxKappab-minKappab)/(nPoints-1.) for i in xrange(nPoints) ]
            totalXSaxis = []

            for kappabVal in [ 0., 1. ] + kappabAxis:
                parametrization.kappab = kappabVal
                parametrizationResult = parametrization.get_output_container()

                # Calculate integral
                binWidths = [ 0.5*(parametrizationResult.binBoundaries[i]+parametrizationResult.binBoundaries[i+1]) for i in xrange(len(parametrizationResult.binBoundaries)-1) ]
                integral = sum([
                    width * ratio * SMXS for width, ratio, SMXS in zip( binWidths, parametrizationResult.ratios, SM.crosssection )
                    ])

                totalXSaxis.append( integral )

            totalXS_at_kappab_zero = totalXSaxis.pop(0)
            totalXS_at_kappab_one  = totalXSaxis.pop(0)

            totalXSaxis = [ xs / totalXS_at_kappab_one for xs in totalXSaxis ]


            c.Clear()
            set_cmargins()


            # yMin = min(totalXSaxis) - 0.1*(max(totalXSaxis)-min(totalXSaxis)),
            # yMax = max(totalXSaxis) + 0.1*(max(totalXSaxis)-min(totalXSaxis)),

            yMin = 0.7
            yMax = 2.3

            base = get_plot_base(
                xMin = minKappab,
                xMax = maxKappab,
                yMin = yMin,
                yMax = yMax,
                xTitle = '#kappa_{b}',
                yTitle = 'totalXS',
                )
            base.Draw('P')

            base.SetTitle('totalXS')

            Tg = ROOT.TGraph(
                nPoints,
                array( 'd', kappabAxis ),
                array( 'd', totalXSaxis ),
                )
            Tg.SetLineWidth(2)
            Tg.SetLineColor(2)
            Tg.Draw('LSAME')

            lineAtOne = ROOT.TLine( minKappab, 1.0, maxKappab, 1.0 )
            lineAtOne.Draw('L')


            save_c( 'kappabDependency_totalXS_fromHistograms_{0}'.format(key) )


            if True:

                minY = min(totalXSaxis)
                minX = kappabAxis[ totalXSaxis.index(minY) ]

                print 'minimum: kappab = {0}, {1} = {2}'.format( minX, 'totalXS', minY )

                print 'kappab = 0.0, {0} = {1}'.format( 'totalXS', totalXS_at_kappab_zero / totalXS_at_kappab_one )


    #____________________________________________________________________
    if args.CorrelationMatrixScaleDependence_Yukawa:

        CorrelationMatrices.set_plot_dir( 'plots_CorrelationMatrixCrosscheck_{0}'.format(datestr) )
        TheoryCommands.set_plot_dir( 'plots_CorrelationMatrixCrosscheck_{0}'.format(datestr) )

        expCorrMatrices = []
        theoryCorrMatrices = []

        kappabs = [ -2, -1, 0, 1, 2 ]
        kappacs = [ -10, -5, 0, 1, 5, 10 ]
        # kappabs = [ -1, 0, 1 ]
        # kappacs = [ -5, 1, 5 ]

        for kappab, kappac in itertools.product( kappabs, kappacs ):
            print 'Processing kappab = {0}, kappac = {1}'.format( kappab, kappac )

            variationFiles = TheoryFileInterface.file_finder(
                directory = LatestPaths.derivedTheoryFiles_YukawaGluonInduced,
                # directory = LatestPaths.derivedTheoryFilesDirectory_YukawaSummed,
                kappab = kappab, kappac = kappac
                )

            variations = [
                TheoryFileInterface.read_derived_theory_file( variationFile, returnContainer=True )
                    for variationFile in variationFiles ]

            theoryCorrMatrices.append( CorrelationMatrices.get_correlation_matrix(
                variations,
                makeScatterPlots          = False,
                makeCorrelationMatrixPlot = True,
                outname                   = 'corrMat_theory_kappab_{0}_kappac_{1}'.format( kappab, kappac ),
                verbose                   = True,
                ))

            variations_expbinning = deepcopy(variations)
            for variation in variations_expbinning:
                TheoryCommands.rebin_derived_theory_container( variation, [ 0., 15., 30., 45., 85., 125. ] )

            expCorrMatrices.append( CorrelationMatrices.get_correlation_matrix(
                variations_expbinning,
                makeScatterPlots          = False,
                makeCorrelationMatrixPlot = True,
                outname                   = 'corrMat_exp_kappab_{0}_kappac_{1}'.format( kappab, kappac ),
                verbose                   = True,
                ))


        # ======================================
        # Min-max study

        print '\nCreating min-max matrix'

        toString = lambda number: str(int(number)) if number.is_integer() else '{0:.1f}'.format(number)
        def getBinLabels( variation ):
            binLabels = []
            for iBin in xrange(len(variation.binCenters)):
                if not hasattr( variation, 'binBoundaries' ):
                    binLabel = toString(variation.binCenters[iBin])
                else:
                    binLabel = '{0} - {1}'.format(
                        toString(variation.binBoundaries[iBin] ),
                        toString(variation.binBoundaries[iBin+1] )
                        )
                binLabels.append(binLabel)
            return binLabels

        expBinLabels    = getBinLabels( variations_expbinning[0] )
        theoryBinLabels = getBinLabels( variations[0] )


        for name, corrMatrices, binLabels in [
                ( 'exp', expCorrMatrices, expBinLabels ),
                ( 'theory', theoryCorrMatrices, theoryBinLabels )
                ]:

            nBins = len(corrMatrices[0])

            # Construct minimum and maximum corrMatrix
            minCorrMatrix = [ [ 999.  for j in xrange(nBins) ] for i in xrange(nBins) ]
            maxCorrMatrix = [ [ -999. for j in xrange(nBins) ] for i in xrange(nBins) ]
            degreeOfAsymmetryMatrix = [ [ 0. for j in xrange(nBins) ] for i in xrange(nBins) ]
            for iRow in xrange(nBins):
                for iCol in xrange(nBins):
                    minVal = 999.
                    maxVal = -999.
                    for corrMatrix in corrMatrices:
                        if corrMatrix[iRow][iCol] < minVal:
                            minVal = corrMatrix[iRow][iCol]
                        if corrMatrix[iRow][iCol] > maxVal:
                            maxVal = corrMatrix[iRow][iCol]
                    minCorrMatrix[iRow][iCol] = minVal
                    maxCorrMatrix[iRow][iCol] = maxVal
                    degreeOfAsymmetryMatrix[iRow][iCol] = maxVal - minVal


            # Make a plot

            c.Clear()
            set_cmargins(
                LeftMargin   = 0.18,
                RightMargin  = 0.12,
                TopMargin    = 0.06,
                BottomMargin = 0.16,
                )

            T = ROOT.TH2D(
                'corrMat', '#scale[0.85]{Min and max corr. for p_{T} bins}',
                nBins, 0., nBins,
                nBins, 0., nBins
                )
            ROOT.SetOwnership( T, False )

            for iRow in xrange(nBins):
                for iCol in xrange(nBins):
                    T.SetBinContent( iCol+1, iRow+1, degreeOfAsymmetryMatrix[iRow][iCol] )

            # Bin titles
            for iBin in xrange(nBins):
                if nBins < 20 or iBin % int(0.1*nBins) == 0:
                    T.GetXaxis().SetBinLabel( iBin+1, binLabels[iBin] )
                    T.GetYaxis().SetBinLabel( iBin+1, binLabels[iBin] )

            # Set color range

            # # With negative (to orange)
            # n_stops = 3
            # stops  = [ 0.0, 0.5, 1.0 ]
            # reds   = [ 0.0, 1.0, 1.0 ]
            # blues  = [ 1.0, 1.0, 0.0 ]
            # greens = [ 0.0, 1.0, 0.4 ]
            # zMin   = -1.2
            # zMax   = 1.2

            # With only positive
            n_stops = 2
            stops  = [ 0.0, 1.0 ]
            reds   = [ 1.0, 255./255. ]
            greens = [ 1.0, 191./255. ]
            blues  = [ 1.0, 128./255. ]
            zMin   = 0.
            zMax   = 1.05 * max([ max(row) for row in degreeOfAsymmetryMatrix ])


            ROOT.TColor.CreateGradientColorTable(
                n_stops,
                array('d', stops ),
                array('d', reds ),
                array('d', greens ),
                array('d', blues ),
                255 )
            T.GetZaxis().SetRangeUser( zMin, zMax )

            T.GetXaxis().SetTitle( 'p_{T}' )
            # T.GetXaxis().SetTitleOffset( 1.7 )
            T.GetXaxis().SetTitleOffset( 1.2 )
            T.GetXaxis().SetTitleSize(0.05)

            T.GetYaxis().SetTitle( 'p_{T}' )
            # T.GetYaxis().SetTitleOffset( 2.45 )
            T.GetYaxis().SetTitleOffset( 1.85 )
            T.GetYaxis().SetTitleSize(0.05)

            T.GetXaxis().SetLabelSize(0.045)
            T.GetYaxis().SetLabelSize(0.045)

            T.Draw('COLZ')

            c.cd()
            c.Update()


            # Draw in the min and max numbers

            labelMin = ROOT.TLatex()
            labelMin.SetTextSize(0.035)
            labelMin.SetTextColor(4)
            labelMin.SetTextAlign(23)

            labelMax = ROOT.TLatex()
            labelMax.SetTextSize(0.035)
            labelMax.SetTextColor(2)
            labelMax.SetTextAlign(21)

            xMin = T.GetXaxis().GetBinLowEdge(1)
            xMax = T.GetXaxis().GetBinUpEdge(nBins)
            yMin = T.GetYaxis().GetBinLowEdge(1)
            yMax = T.GetYaxis().GetBinUpEdge(nBins)
            dy = yMax - yMin

            significance = 2
            formatter = lambda number: '{0:+.{significance}f}'.format(
                number, significance=significance )

            if nBins >= 20:
                labelMin.SetTextSize(0.008)
                labelMax.SetTextSize(0.008)
                dy *= 0.0
                significance = 1
                formatter = lambda number: '{0:+.{significance}f}'.format(
                    number, significance=significance ).replace('+0','+').replace('-0','-').replace('1.0','1')


            for iRow in xrange(nBins):
                for iCol in xrange(nBins):

                    if (
                        nBins < 20
                        or ( iRow % int(0.1*nBins) == 0 and iCol % int(0.1*nBins) == 0 )
                        or ( degreeOfAsymmetryMatrix[iRow][iCol] > 0.5 )
                        ):

                        x = T.GetXaxis().GetBinCenter(iCol+1)
                        y = T.GetYaxis().GetBinCenter(iRow+1)

                        labelMin.DrawLatex(
                            x, y - 0.01*dy,
                            formatter( minCorrMatrix[iRow][iCol] )
                            )

                        labelMax.DrawLatex(
                            x, y + 0.01*dy,
                            formatter( maxCorrMatrix[iRow][iCol] )
                            )

            save_c( '_minmax_corrMat_{0}'.format(name), asPNG=True, asROOT=True )


    #____________________________________________________________________
    if args.CreateParametrizationTables:

        titles = {
            'kappab' : '\\kappa_b',
            'kappac' : '\\kappa_c',
            'ct'     : '\\kappa_t',
            'cg'     : '\\kappa_g',
            'cb'     : '\\kappa_b',
            }

        wss = [
            LatestPaths.ws_combined_Yukawa,
            LatestPaths.ws_combined_Top,
            LatestPaths.ws_combined_TopCtCb,
            ]

        for ws in wss:

            with Commands.OpenRootFile( ws ) as rootFp:
                w = rootFp.Get('w')

            theoryBinBoundariesArgList = ROOT.RooArgList( w.set('theoryBinBoundaries') )
            nTheoryBins         = theoryBinBoundariesArgList.getSize()
            theoryBinBoundaries = [ 0. ] + [ theoryBinBoundariesArgList[i].getVal() for i in xrange(nTheoryBins) ] 

            def numToStr( number ):
                if float(number).is_integer():
                    return str(int(number))
                else:
                    return '{0:.2f}'.format(number)

            texText = [ '% ' + Commands.tag_git_commit_and_module() ]
            texText.append( '% Used workspace: {0}'.format(ws) )
            for iParam in xrange(nTheoryBins):
                with Commands.RedirectStdout() as redirected:
                    w.function( 'parametrization{0}'.format(iParam) ).Print()
                    formulaPrinted = redirected.read()

                match = re.search( r'formula=\"(.*?)\"', formulaPrinted )
                if not match:
                    print 'No match for parametrization {0}; continuing'.format(iParam)
                    continue
                formula = match.group(1)

                match = re.search( r'actualVars=\((.*?)\)', formulaPrinted )
                if not match:
                    print 'No match for parametrization {0}; continuing'.format(iParam)
                    continue
                actualVars = match.group(1).split(',')


                formula = formula.replace( '+-', '-' )
                for iVar, actualVar in enumerate(actualVars):

                    formula = formula.replace( '@{0}*@{0}'.format(iVar), '\,' + titles.get( actualVar, actualVar ) + '^2' )
                    formula = formula.replace( '@{0}'.format(iVar),      '\,' + titles.get( actualVar, actualVar ) )

                formula = re.sub( r'(\d\.\d\d\d\d)\d*', r'\1', formula )

                formula = formula.replace( '*', '' )
                if formula.startswith('('): formula = formula[1:]
                if formula.endswith(')'):   formula = formula[:-1]

                formula = re.sub( r'e(\-\d+)', r'\\cdot 10^{\1}', formula )

                texText.append(
                    '{0:3} & $({1}-{2})$ & ${3}$ \\\\'.format(
                        iParam,
                        numToStr(theoryBinBoundaries[iParam]),
                        numToStr(theoryBinBoundaries[iParam+1]),
                        formula
                        )
                    )

            if Commands.is_test_mode():
                print '\n'.join(texText)
            else:

                outdir = 'parametrizationTexs_{0}'.format(datestr)
                if not isdir(outdir): os.makedirs(outdir)
                texFile = join( outdir, 'parametrization_{0}.tex'.format( basename(ws).replace('/','').replace('.root','') ) )

                print 'Writing output to', texFile

                with open(texFile, 'w' ) as texFp:
                    texFp.write( '\n'.join(texText) )



########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'