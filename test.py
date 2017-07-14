#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, itertools, operator, re, argparse, sys
from math import isnan, isinf
from os.path import *
from glob import glob


import Commands
import PhysicsCommands
import OneOfCommands
import TheoryCommands


from time import strftime
datestr = strftime( '%b%d' )


########################################
# Main
########################################

def main():


    import argparse
    parser = argparse.ArgumentParser()
    # parser.add_argument( '--string', type=str, default='default', help='default string' )
    parser.add_argument( '--test', action='store_true' )
    parser.add_argument( '--makeStewartTackmannDatacard', action='store_true' )
    parser.add_argument( '--plot', action='store_true' )
    parser.add_argument( '--plot2D', action='store_true' )
    parser.add_argument( '--parametrize', action='store_true' )
    parser.add_argument( '--t2ws', action='store_true' )
    parser.add_argument( '--couplingT2WS', action='store_true' )
    parser.add_argument( '--combineCards', action='store_true' )
    parser.add_argument( '--couplingBestfit', action='store_true' )
    parser.add_argument( '--couplingImportanceHistogram', action='store_true' )
    parser.add_argument( '--doFastscan', action='store_true' )
    parser.add_argument( '--createDerivedTheoryFiles', action='store_true' )
    parser.add_argument( '--createDerivedTheoryFiles_Yukawa', action='store_true' )
    parser.add_argument( '--rebinnedTheoryPlot', action='store_true' )

    # parser.add_argument( '--list', metavar='N', type=str, nargs='+', help='list of strings' )
    args = parser.parse_args()


    if args.test:
        Commands.TestMode()


    if args.makeStewartTackmannDatacard:

        container = TheoryCommands.ReadDerivedTheoryFile( 'derivedTheoryFiles_Jun22/SM_NNLO.txt', returnContainer=True )
        
        covMat = TheoryCommands.GetStewartTackmannCovarianceMatrix( container )

        TheoryCommands.AddCovarianceMatrixAsNuisanceParameters(
            'suppliedInput/combinedCard_May15.txt',
            covMat
            )





    if args.couplingT2WS:

        # Commands.TestMode()

        Commands.BasicT2WSwithModel(
            # 'suppliedInput/combinedCard_May15.txt',
            'suppliedInput/combinedCard_May15_theoryUncertainties.txt',
            'CouplingModel.py',
            extraOptions = [
                '--PO verbose=1',
                '--PO \'higgsMassRange=123,127\'',
                '--PO SM=[ct=1.0,cg=0.0,file={0}]'.format( abspath('derivedTheoryFiles_May23/SM_NNLO.txt') ),
                '--PO theory=[ct=0.1,cg=0.075,file={0}]'.format( abspath('derivedTheoryFiles_May23/ct_0p1_cg_0p075.txt') ),
                '--PO theory=[ct=0.5,cg=0.042,file={0}]'.format( abspath('derivedTheoryFiles_May23/ct_0p5_cg_0p042.txt') ),
                '--PO theory=[ct=1.5,cg=-0.042,file={0}]'.format( abspath('derivedTheoryFiles_May23/ct_1p5_cg_m0p042.txt') ),
                '--PO theory=[ct=2.0,cg=-0.083,file={0}]'.format( abspath('derivedTheoryFiles_May23/ct_2p0_cg_m0p083.txt') ),
                ]
            )


    if args.couplingImportanceHistogram:

        if args.doFastscan:

            datacard = 'workspaces_May30/combinedCard_May15.root'
            Commands.MultiDimCombineTool(
                datacard,
                nPoints       = 10000,
                nPointsPerJob = 10,
                queue         = 'all.q',
                # notOnBatch    = False,
                notOnBatch    = True,
                jobDirectory  = 'Fastscan_couplings_{0}'.format( datestr ),
                extraOptions  = [
                    '-P ct -P cg',
                    '--setPhysicsModelParameters ct=1.0,cg=0.0',
                    '--setPhysicsModelParameterRanges ct=-0.5,4.2:cg=-0.32,0.1',
                    '--fastScan',
                    '--saveSpecifiedFunc {0}'.format( ','.join(Commands.ListSet( datacard, 'yieldParameters' )) ),
                    ]
                )
            TheoryCommands.WriteTH2DToFile(
                glob( 'Fastscan_couplings_{0}/*.root'.format( datestr ) )
                )

        else:

            TheoryCommands.WriteTH2DToFile(
                # glob( 'Scan_couplings_May31_3/*.root' )
                glob( 'Fastscan_couplings_{0}/*.root'.format( datestr ) )
                )




    if args.couplingBestfit:

        TESTFIT = False

        # datacard = 'workspaces_May30/combinedCard_May15.root'
        datacard = 'workspaces_Jun22/combinedCard_May15_theoryUncertainties.root'

        if TESTFIT:
            Commands.BasicBestfit(
                datacard,
                setPOIs = False,
                extraOptions = [
                    '--setPhysicsModelParameters ct=1.0,cg=0.0'
                    ]
                )

        else:
            
            Commands.MultiDimCombineTool(
                datacard,
                nPoints       = 1600,
                # nPoints       = 100,
                nPointsPerJob = 10,
                queue         = 'all.q',
                notOnBatch    = False,
                # notOnBatch    = True,
                jobDirectory  = 'Scan_couplings_{0}'.format( datestr ),
                extraOptions  = [
                    '--importanceSampling={0}:couplingScan'.format( abspath('scanTH2D_Jun01.root') ),
                    '-P ct -P cg',
                    '--setPhysicsModelParameters ct=1.0,cg=0.0',
                    '--setPhysicsModelParameterRanges ct=-0.5,4.2:cg=-0.32,0.1',
                    # '--fastScan', # TURN THIS OFF FOR REAL RUN!!!!
                    '--saveSpecifiedFunc {0}'.format(','.join(
                        Commands.ListSet( datacard, 'yieldParameters' ) + [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
                    '--squareDistPoiStep',
                    ]
                )

        # datacard = 'workspaces_May15/hzz4l_comb_13TeV_xs.root'
        # Commands.MultiDimCombineTool(
        #     datacard,
        #     # nPoints       = 1600,
        #     nPoints       = 100,
        #     nPointsPerJob = 10,
        #     queue         = 'all.q',
        #     # notOnBatch    = False,
        #     notOnBatch    = True,
        #     jobDirectory  = 'Scan_couplings_{0}_hzz1'.format( datestr ),
        #     extraOptions  = [
        #         '--importanceSampling={0}:couplingScan'.format( abspath('scanTH2D_Jun01.root') ),
        #         # '-P ct -P cg',
        #         '-P r_smH_PTH_0_15 -P r_smH_PTH_15_30',
        #         # '--setPhysicsModelParameters ct=1.0,cg=0.0',
        #         # '--setPhysicsModelParameterRanges ct=-0.5,4.2:cg=-0.32,0.1',
        #         '--fastScan', # TURN THIS OFF FOR REAL RUN!!!!
        #         # '--saveSpecifiedFunc {0}'.format( ','.join(Commands.ListSet( datacard, 'yieldParameters' )) ),
        #         # '--squareDistPoiStep',
        #         ]
        #     )




    if args.createDerivedTheoryFiles:
        TheoryCommands.CreateDerivedTheoryFiles( pattern=r'ct_[mp\d]+_cg_[mp\d]+' )

    if args.createDerivedTheoryFiles_Yukawa:
        TheoryCommands.CreateDerivedTheoryFiles_Yukawa()


    # ======================================
    # combineCards.py

    if args.combineCards:

        Commands.BasicCombineCards(
            'suppliedInput/combinedCard_{0}.txt'.format(datestr),
            'suppliedInput/pT_v5_renamed/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.txt',
            'suppliedInput/PTH/hzz4l_comb_13TeV_xs.txt'
            )


        Commands.BasicCombineCards(
            'suppliedInput/combinedCard_{0}.txt'.format(datestr),
            'suppliedInput/pT_v5_renamed/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.txt',
            'suppliedInput/PTH/hzz4l_comb_13TeV_xs_processesShifted.txt'
            )


    if args.t2ws:

        # ======================================
        # text2workspace

        # Specific for HZZ
        Commands.BasicT2WS(
            'suppliedInput/hzz_ggH_xH_split_Jun26/hzz4l_all_13TeV_xs.txt',
            manualMaps=[
                '--PO \'map=.*/ggH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
                '--PO \'map=.*/xH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
                '--PO \'map=.*/xH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_30_45:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/xH_PTH_30_45:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_45_85:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/xH_PTH_45_85:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_85_125:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/xH_PTH_85_125:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_125_200:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/xH_PTH_125_200:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_200_350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/xH_PTH_200_350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_GT350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/xH_PTH_GT350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
                ],
            )



        # Commands.BasicT2WS( 'suppliedInput/pT_v5_renamed/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.txt' )
        # Commands.BasicT2WS( 'suppliedInput/PTH/hzz4l_comb_13TeV_xs.txt' )
        # Commands.BasicT2WS( 'suppliedInput/fromDavid/PTH_May09/hzz4l_comb_13TeV_xs.txt' )
        # Commands.BasicT2WS( 'suppliedInput/combinedCard_{0}.txt'.format(datestr) )


        # Commands.BasicT2WS( 'suppliedInput/pT_v5_renamed/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.txt' )

        # Commands.BasicT2WS(
        #     'suppliedInput/PTH/hzz4l_comb_13TeV_xs.txt',
        #     # extraOptions = '--X-no-check-norm'
        #     )


        # Commands.hzz_T2WS(
        #     'suppliedInput/PTH/hzz4l_comb_13TeV_xs.txt',
        #     # extraOptions = '--X-no-check-norm'
        #     )

        # Commands.BasicT2WS(
        #     'suppliedInput/combinedCard_{0}.txt'.format(datestr),
        #     # extraOptions = '--X-no-check-norm'
        #     )


    # ======================================
    # Basic best fit - for testing purposes

    # hgg
    # Commands.BasicBestfit(
    #     'workspaces_May15/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.root',
    #     # onBatch = True
    #     )
    
    # hzz
    # Commands.BasicBestfit(
    #     'workspaces_May10/hzz4l_comb_13TeV_xs.root',
    #     # setPOIs = False,
    #     # onBatch = True
    #     )



    # May 15: workspaces with included but faulty _norm

    # extraVars = []
    # nProcesses = len(Commands.ListProcesses( 'suppliedInput/PTH/hzz4l_comb_13TeV_xs.txt' )[0])
    # extraVars.extend( [ 'K1Bin{0}'.format(i) for i in xrange(nProcesses) ] )
    # extraVars.extend( [ 'K2Bin{0}'.format(i) for i in xrange(nProcesses) ] )

    # # hzz
    # Commands.BasicBestfit(
    #     'workspaces_May15/hzz4l_comb_13TeV_xs.root',
    #     setPOIs = False,
    #     # onBatch = True,
    #     extraOptions = (
    #         '--saveSpecifiedFunc {0}'.format( ','.join(extraVars) )
    #         + ' -m 125.09 --floatOtherPOIs=1 --freezeNuisances MH'
    #         )
    #     )



    # Commands.BasicBestfit(
    #     'workspaces_May15/combinedCard_May15.root',
    #     # setPOIs = False,
    #     # onBatch = True
    #     extraOptions = (
    #         ' --floatOtherPOIs=1'
    #         # + ' -m 125.09 --freezeNuisances MH'
    #         )
    #     )



    # ======================================
    # combineTool.py

    # Commands.BasicCombineTool(
    #     'workspaces_Apr27/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.root',
    #     POIpattern    = '30_45',
    #     nPoints       = 10,
    #     nPointsPerJob = 2,
    #     )



    # # Ran on 15 May

    # Commands.BasicCombineTool(
    #     'workspaces_May15/usingDavidsCommand_OAshifted.root',
    #     POIpattern    = '*',
    #     nPoints       = 40,
    #     nPointsPerJob = 8,
    #     # notOnBatch    = True,
    #     jobDirectory  = 'Scan_{0}'.format(datestr),
    #     queue         = 'short.q',
    #     )

    # Commands.BasicCombineTool(
    #     'workspaces_May15/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.root',
    #     POIpattern    = '*',
    #     nPoints       = 40,
    #     nPointsPerJob = 4,
    #     # notOnBatch    = True,
    #     jobDirectory  = 'Scan_{0}'.format(datestr),
    #     queue         = 'short.q',
    #     )

    # Commands.BasicCombineTool(
    #     'workspaces_May15/combinedCard_May15.root',
    #     POIpattern    = '*',
    #     nPoints       = 40,
    #     nPointsPerJob = 4,
    #     # notOnBatch    = True,
    #     jobDirectory  = 'Scan_{0}'.format(datestr),
    #     queue         = 'short.q',
    #     )





    # ======================================
    # Drawing

    if args.plot:

        # res = Commands.ConvertTChainToArray(
        #     'limit',
        #     [
        #         'test001_May03/higgsCombineSCAN_May03_Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.POINTS.0.3.MultiDimFit.mH125.root',
        #         'test001_May03/higgsCombineSCAN_May03_Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.POINTS.12.15.MultiDimFit.mH125.root',
        #         'test001_May03/higgsCombineSCAN_May03_Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.POINTS.16.19.MultiDimFit.mH125.root',
        #         'test001_May03/higgsCombineSCAN_May03_Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.POINTS.4.7.MultiDimFit.mH125.root',
        #         'test001_May03/higgsCombineSCAN_May03_Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.POINTS.8.11.MultiDimFit.mH125.root',
        #         ],
        #     'r_',
        #     )
        # for i in res['r_smH_PTH_30_45']:
        #     print i
        # print
        # print len(res['r_smH_PTH_30_45'])



        hggPOIs = Commands.ListPOIs( 'workspaces_May15/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.root' )
        hggscans = PhysicsCommands.GetScanResults(
            hggPOIs,
            'Scan_May15',
            pattern = 'Datacard_13TeV_differential_pT'
            )
        PhysicsCommands.BasicDrawScanResults( hggPOIs, hggscans, name='hgg' )
        PhysicsCommands.BasicDrawSpectrum( hggPOIs, hggscans, name='hgg' )


        hzzPOIs = Commands.ListPOIs( 'workspaces_May15/usingDavidsCommand_OAshifted.root' )
        hzzscans = PhysicsCommands.GetScanResults(
            hzzPOIs,
            'Scan_May15',
            pattern = 'usingDavidsCommand_OAshifted'
            )
        PhysicsCommands.BasicDrawScanResults( hzzPOIs, hzzscans, name='hzz' )
        PhysicsCommands.BasicDrawSpectrum( hzzPOIs, hzzscans, name='hzz' )


        combinationPOIs = Commands.ListPOIs( 'workspaces_May15/combinedCard_May15.root' )
        combinationscans = PhysicsCommands.GetScanResults(
            combinationPOIs,
            'Scan_May15',
            pattern = 'combinedCard'
            )
        PhysicsCommands.BasicDrawScanResults( combinationPOIs, combinationscans, name='combination' )
        PhysicsCommands.BasicDrawSpectrum( combinationPOIs, combinationscans, name='combination' )


        PhysicsCommands.BasicCombineSpectra(
            ( 'combination', combinationPOIs, combinationscans,
                ( 'SetLineColor', 1 ),
                ( 'SetMarkerStyle', 2 ),
                # ( 'SetFillColorAlpha', 1, 0.2 ),
                ( 'SetFillColor', 13 ),
                # ( 'SetFillStyle', 3544 ),
                ( 'SetFillStyle', 3345 ),
                ),
            ( 'hgg', hggPOIs, hggscans,
                ( 'SetLineColor', 2 ),
                ( 'SetMarkerStyle', 5 ),
                ( 'SetFillColorAlpha', 2, 0.2 ),
                ),
            ( 'hzz', hzzPOIs, hzzscans,
                ( 'SetLineColor', 4 ),
                ( 'SetMarkerStyle', 8 ),
                ( 'SetFillColorAlpha', 4, 0.2 ),
                ),
            printTable=True,
            )


        PhysicsCommands.BasicDrawMultipleParabolas(
            ( 'combination', combinationPOIs, combinationscans,
                ( 'SetLineColor', 1 ),
                # ( 'SetMarkerStyle', 2 ),
                # ( 'SetFillColorAlpha', 1, 0.2 ),
                # ( 'SetFillColor', 13 ),
                # ( 'SetFillStyle', 3544 ),
                # ( 'SetFillStyle', 3345 ),
                ),
            ( 'hgg', hggPOIs, hggscans,
                ( 'SetLineColor', 2 ),
                # ( 'SetMarkerStyle', 5 ),
                # ( 'SetFillColorAlpha', 2, 0.2 ),
                ),
            ( 'hzz', hzzPOIs, hzzscans,
                ( 'SetLineColor', 4 ),
                # ( 'SetMarkerStyle', 8 ),
                # ( 'SetFillColorAlpha', 4, 0.2 ),
                ),
            )


        # Combination only with theory curves
        for pattern in [
            r'^cg_[mp\d]+$',
            r'^ct_[mp\d]+$',
            r'^cb_[mp\d]+$',
            r'ct_[mp\d]+_cb_[mp\d]+_cg_[mp\d]+',
            r'ct_[mp\d]+_cb_([mp\d]+)$', # ct_XX_cb_XX , but not ct_XX_cb_XX_cg_XX
            r'ct_[mp\d]+_cg_[mp\d]+',
            ]:

            theoryTgs = TheoryCommands.LoadTheoryCurves( pattern )

            PhysicsCommands.BasicCombineSpectra(
                ( 'combination', combinationPOIs, combinationscans,
                    ( 'SetLineColor', 1 ),
                    ( 'SetMarkerStyle', 2 ),
                    ( 'SetFillColorAlpha', 1, 0.2 ),
                    # ( 'SetFillColor', 13 ),
                    # ( 'SetFillStyle', 3544 ),
                    # ( 'SetFillStyle', 3345 ),
                    ),
                # ( 'hgg', hggPOIs, hggscans,
                #     ( 'SetLineColor', 2 ),
                #     ( 'SetMarkerStyle', 5 ),
                #     ( 'SetFillColorAlpha', 2, 0.2 ),
                #     ),
                # ( 'hzz', hzzPOIs, hzzscans,
                #     ( 'SetLineColor', 4 ),
                #     ( 'SetMarkerStyle', 8 ),
                #     ( 'SetFillColorAlpha', 4, 0.2 ),
                #     ),
                theoryCurves=theoryTgs
                )


    # ======================================
    # Try parametrization of theory curves

    if args.parametrize:

        # Only to get the pt axis
        pt_variation, ratio_ct_1p1 = TheoryCommands.ReadTheoryFile( 'suppliedInput/numbersFromAgnieszka/HRes_SMEFT_May16/ratio_ctup_new' )

        # Get the following shapes:
        # ct = 0.5, cg = 0.042
        # ct = 0.1, cg = 0.075
        # ct = 2.0, cg = -0.083
        # ct = 1.5, cg = -0.042
        Tgs = TheoryCommands.LoadTheoryCurves( r'ct_[mp\d]+_cg_[mp\d]+' )

        parametrization = TheoryCommands.GetParametrization(
            [ ( Tg.ct**2, Tg.ct*Tg.cg, Tg.cg**2 ) for Tg in Tgs ],
            [ Tg.binValues for Tg in Tgs ],
            # testMode=True
            testMode=False
            )


        # Test it out
        def GetParametrizedTg( ct, cg ):
            yVals = parametrization( ct**2, ct*cg, cg**2 )
            Tg = TheoryCommands.GetTheoryTGraph(
                'ct_{0}_cg_{1}_param'.format( str(ct).replace('.','p').replace('-','m'), str(cg).replace('.','p').replace('-','m') ),
                pt_variation, yVals,
                # prodlist( ratio_SM_NNLO_down, ratio ),
                # prodlist( ratio_SM_NNLO_up, ratio )
                )
            return Tg


        Tgs_copy = Tgs[:]
        Tgs_copy.append( GetParametrizedTg( ct = 1.5, cg = -0.042 ) )
        Tgs_copy.append( GetParametrizedTg( ct = 1.5, cg = -0.043 ) )
        Tgs_copy.append( GetParametrizedTg( ct = 2.0, cg = -0.083 ) )
        TheoryCommands.BasicTheoryPlot( Tgs_copy, drawErrors=False )



        # Scale a theory graph

        experimentalBinning = PhysicsCommands.FigureOutBinning( Commands.ListPOIs( 'workspaces_May15/combinedCard_May15.root' ) )

        ptAxisDummy, SM_NNLO_variationBinned, ratio_SM_NNLO_down_variationBinned, ratio_SM_NNLO_up_variationBinned = \
            TheoryCommands.LoadStandardModelCurves( ptAxis=Tgs[0].binCenters )



        # Make test plot only if other curves are loaded as well
        if args.plot:


            Tgs_expBinning = []
            for Tg in Tgs:

                Tg_expBinning = TheoryCommands.MapPredictionToExperimental(
                    Tg.binBoundaries,
                    map( operator.mul, Tg.binValues, SM_NNLO_variationBinned ),
                    experimentalBinning,
                    makeTGraph = Tg.name + '_rebin',
                    # verbose=True,
                    )

                Tgs_expBinning.append( Tg_expBinning )


            PhysicsCommands.BasicCombineSpectra(
                ( 'combination', combinationPOIs, combinationscans,
                    ( 'SetLineColor', 1 ),
                    ( 'SetMarkerStyle', 2 ),
                    ( 'SetFillColorAlpha', 1, 0.2 ),
                    # ( 'SetFillColor', 13 ),
                    # ( 'SetFillStyle', 3544 ),
                    # ( 'SetFillStyle', 3345 ),
                    ),
                theoryCurves=Tgs_expBinning
                )

                

    if args.plot2D:

        scan = TheoryCommands.PlotCouplingScan2D(
            # 'workspaces_May30/combinedCard_May15.root',
            # glob( 'Scan_couplings_Jun01_last/*.root' ),
            'workspaces_Jun22/combinedCard_May15_theoryUncertainties.root',
            glob( 'Scan_couplings_Jun22/*.root' ),
            )

        scan = [ s for s in scan if not ( isnan(s[2]) or isinf(s[2]) ) ]
        scan.sort( key = lambda i: i[2] )

        # Get first N unique results
        uniquePoints = []
        allNLLs    = [ s[2] for s in scan ]
        uniqueNLLs = list(set( allNLLs ))
        uniqueNLLs.sort()
        for nll in uniqueNLLs:
            iPoint = allNLLs.index(nll)
            uniquePoints.append( scan[iPoint] )

        testcouplings = []
        for cg, ct, deltaNLL in [ uniquePoints[0], uniquePoints[4], uniquePoints[8] ] :
            testcouplings.append( { 'ct' : ct, 'cg' : cg } )


        if args.plot:


            # ======================================
            # Theory curves from files

            SM = TheoryCommands.ReadDerivedTheoryFile( 'derivedTheoryFiles_May23/SM_NNLO.txt' )
            theories = []
            for derivedTheoryFile in glob( 'derivedTheoryFiles_May23/c*.txt' ):
                theory = TheoryCommands.ReadDerivedTheoryFile(derivedTheoryFile)
                theory['name'] = basename(derivedTheoryFile).replace('.txt','')
                theories.append( theory )

            expBinBoundaries = PhysicsCommands.FigureOutBinning( combinationPOIs )
            rebinnedSMXS = TheoryCommands.Rebin(
                SM['binBoundaries'], SM['crosssection'],
                expBinBoundaries,
                )

            Tgs = []
            colorCycle = itertools.cycle( range(2,5) + range(6,10) )
            for theory in theories:
                color = next(colorCycle)

                rebinnedXS = TheoryCommands.Rebin(
                    theory['binBoundaries'], theory['crosssection'],
                    expBinBoundaries,
                    verbose=False,
                    )
                rebinnedRatios = [ xs / smxs for xs, smxs in zip( rebinnedXS, rebinnedSMXS ) ]

                Tg_theoryBinning = TheoryCommands.GetTheoryTGraph(
                    theory['name'] + '_tb',
                    theory['binBoundaries'],
                    theory['ratios'],
                    boundaries = True,
                    )
                Tg_theoryBinning.SetLineColor( color )

                Tg_expBinning = TheoryCommands.GetTheoryTGraph(
                    theory['name'] + '_eb',
                    expBinBoundaries,
                    rebinnedRatios,
                    boundaries = True,
                    )
                Tg_expBinning.SetLineColor( color )

                # Tgs.append( Tg_theoryBinning )
                # Tgs.append( Tg_expBinning )


            # ======================================
            # Theory curves from parametrization in workspace

            # testcouplings = [
            #     { 'ct' : 3.00, 'cg' : -0.18 },
            #     { 'ct' : 2.00, 'cg' : -0.083 },
            #     { 'ct' : 1.40, 'cg' : -0.05 }
            #     # { 'ct' : 1.75, 'cg' : -0.0625 },
            #     # { 'ct' : 0.35, 'cg' : 0.0585 },
            #     ]

            yValuesPerTestCoupling = TheoryCommands.TestParametrizationsInWorkspace(
                # 'workspaces_May30/combinedCard_May15.root',
                'workspaces_Jun22/combinedCard_May15_theoryUncertainties.root',
                testcouplings
                )

            color = next(colorCycle)
            for couplings, ( y, yParametrization ) in zip( testcouplings, yValuesPerTestCoupling ):

                name = ''
                for key, value in couplings.iteritems():
                    name += '_' + key + '_' + '{0:.2f}'.format(value).replace('.','p').replace('-','m')
                name = name[1:]

                color = next(colorCycle)

                Tg = TheoryCommands.GetTheoryTGraph(
                    name,
                    expBinBoundaries,
                    y,
                    boundaries = True,
                    )
                Tg.SetLineColor( color )
                Tgs.append( Tg )

                yParametrization = yParametrization[1:-1]
                # yParametrization = [ xs / smxs for xs, smxs in zip( yParametrization, SM['crosssection'] ) ]

                Tg = TheoryCommands.GetTheoryTGraph(
                    name + '_param',
                    SM['binBoundaries'],
                    yParametrization,
                    boundaries = True,
                    )
                Tg.SetLineColor( color )
                Tgs.append( Tg )


            # ======================================
            # Plot

            PhysicsCommands.BasicCombineSpectra(
                ( 'combination', combinationPOIs, combinationscans,
                    ( 'SetLineColor', 1 ),
                    ( 'SetMarkerStyle', 2 ),
                    ( 'SetFillColorAlpha', 1, 0.2 ),
                    # ( 'SetFillColor', 13 ),
                    # ( 'SetFillStyle', 3544 ),
                    # ( 'SetFillStyle', 3345 ),
                    ),
                theoryCurves=Tgs[:1],
                drawTheoryCurvesAsLines=True,
                autocolor=False,
                filenameSuffix='_withParametrizationsFromWS'
                )


    if args.rebinnedTheoryPlot:

        container = TheoryCommands.ReadDerivedTheoryFile( 'derivedTheoryFiles_Jun22/SM_NNLO.txt', returnContainer=True )
        TheoryCommands.MakeRebinnedTheoryPlot( container )

        
        




########################################
# End of Main
########################################
if __name__ == "__main__":
    main()