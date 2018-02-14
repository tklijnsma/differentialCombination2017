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
from copy import deepcopy
from array import array

import combineCommands
import plotCommands

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

import LatestPaths

from time import strftime
datestr = strftime( '%b%d' )

import ROOT


########################################
# Main
########################################

def main():

    # ======================================
    # Parser

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument( '--test',                            action='store_true' )

    parser.add_argument( '--simplebestfit',                      action='store_true' )
    parser.add_argument( '--toptest',                            action='store_true' )
    parser.add_argument( '--asimovtest',                            action='store_true' )
    parser.add_argument( '--testfit_njets',                            action='store_true' )

    parser.add_argument( '--BRdependency',                            action='store_true' )
    parser.add_argument( '--BRdependency_Yukawa',                         action='store_true' )

    parser.add_argument( '--debugging_Nov08_hzz',                            action='store_true' )
    parser.add_argument( '--debugging_Nov08_TopWS',                            action='store_true' )

    parser.add_argument( '--jscaletest_RenameHggProcesses_alsoSystematics',      action='store_true' )
    parser.add_argument( '--jscaletest_t2ws',      action='store_true' )
    parser.add_argument( '--jscaletest_bestfit',      action='store_true' )

    parser.add_argument( '--kappaVMaxOne_Top',      action='store_true' )

    parser.add_argument( '--testingDoPoints',      action='store_true' )

    parser.add_argument( '--testing2PanelCanvas',      action='store_true' )


    parser.add_argument( '--hbb_bestfit',      action='store_true' )
    parser.add_argument( '--hbb_scan',      action='store_true' )
    parser.add_argument( '--hbb_print',      action='store_true' )
    parser.add_argument( '--hbbOnly_t2ws',      action='store_true' )
    parser.add_argument( '--hbbOnly_scan',      action='store_true' )

    parser.add_argument( '--statonly_test',      action='store_true' )

    parser.add_argument( '--new_model_implementation', action='store_true' )

    combineCommands.AppendParserOptions(parser)
    plotCommands.AppendParserOptions(parser)

    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument( '--latest', dest='latest', action='store_true', default=True )
    group.add_argument( '--older',  dest='latest', action='store_false' )

    args = parser.parse_args()

    print args
    print ''



    if args.test:
        Commands.TestMode()


    ########################################
    # New commands
    ########################################

    Commands.SetTempJobDir( 'plainWStests_{0}'.format(datestr) )

    base = abspath( join( os.environ['CMSSW_BASE'], 'src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/' ))



    if args.new_model_implementation:

        ws = '/mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/v3_Approval/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/out/workspaces_Feb12/combinedCard_smH_Nov07_DifferentialModel_lumiScale.root'

        cmd = [
            'combine',
            ws,
            '-n _debugging_new_model_implementation_{0}'.format(datestr),
            '-M MultiDimFit',
            '-m 125.0',
            '-t -1',
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            '--saveNLL',
            '--saveInactivePOI 1',
            '--floatOtherPOIs=1',
            # '--job-mode psi --task-name _SCAN_ASIMOV_bPOI_r_smH_PTH_GT350_ePOI_combinedCard_smH_Nov07_DifferentialModel_lumiScale --sub-opts='-q short.q' ',
            '--algo=grid',
            '--points=4',
            # '--split-points 5',
            '-P r_smH_PTH_GT350',
            '--setPhysicsModelParameterRanges r_smH_PTH_0_15=0.2,1.8:r_smH_PTH_15_30=0.2,1.8:r_smH_PTH_30_45=0.2,1.8:r_smH_PTH_45_85=0.2,1.8:r_smH_PTH_85_125=0.2,1.8:r_smH_PTH_125_200=0.2,1.8:r_smH_PTH_200_350=0.2,1.8:r_smH_PTH_GT350=0.2,1.8',
            '--setPhysicsModelParameters r_smH_PTH_0_15=1.0,r_smH_PTH_15_30=1.0,r_smH_PTH_30_45=1.0,r_smH_PTH_45_85=1.0,r_smH_PTH_85_125=1.0,r_smH_PTH_125_200=1.0,r_smH_PTH_200_350=1.0,r_smH_PTH_GT350=1.0',
            ]

        Commands.BasicGenericCombineCommand( cmd, onBatch = False, )


    if args.statonly_test:

        ws = '/mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/v3_Approval/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/out/postfitWss_Feb05/higgsCombine_POSTFIT_combinedCard_smH_Nov07_pth_smH.MultiDimFit.mH125.root'

        cmd = [
            'combine',
            ws,
            '-n debugging_statonly_test_01',
            '-M MultiDimFit',
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            '--algo=grid',
            '--floatOtherPOIs=1',
            '-P "r_smH_PTH_45_85"',
            # '--setPhysicsModelParameters r_smH_PTH_200_350=1.0,r_smH_PTH_125_200=1.0,r_smH_PTH_30_45=1.0,r_smH_PTH_0_15=1.0,r_smH_PTH_GT350=1.0,r_smH_PTH_85_125=1.0,r_smH_PTH_15_30=1.0,r_smH_PTH_45_85=1.0',
            '-m 125.00',
            '--squareDistPoi',
            '--saveNLL',
            '--saveInactivePOI 1',
            '--points=10 ',
            # '--setPhysicsModelParameterRanges "r_smH_PTH_45_85"=-1.000,4.000:"r_smH_PTH_15_30"=-1.000,4.000:"r_smH_PTH_85_125"=-1.000,4.000:"r_smH_PTH_GT350"=-1.000,4.000:"r_smH_PTH_0_15"=-1.000,4.000:"r_smH_PTH_30_45"=-1.000,4.000:"r_smH_PTH_125_200"=-1.000,4.000:"r_smH_PTH_200_350"=-1.000,4.000',
            '--setPhysicsModelParameterRanges "r_smH_PTH_45_85"=0.5000,1.500',
            '--snapshotName MultiDimFit',
            '--skipInitialFit',
            '--freezeNuisances rgx{.*},CMS_zz4l_mass,CMS_hgg_mass,CMS_fakeH_p1_1_8,CMS_fakeH_p3_1_8,CMS_fakeH_p1_2_8,CMS_fakeH_p3_2_8,CMS_fakeH_p1_3_8,CMS_fakeH_p3_3_8',
            # ,K1Bin0,K1Bin1,K1Bin2,K1Bin3,K1Bin4,K1Bin5,K1Bin6,K1Bin7,K2Bin0,K2Bin1,K2Bin2,K2Bin3,K2Bin4,K2Bin5,K2Bin6,K2Bin7
            ]

        Commands.BasicGenericCombineCommand( cmd, onBatch = False, )


    #____________________________________________________________________
    if args.hbbOnly_t2ws:

        card = 'suppliedInput/combinedCard_hbb_debuggingTest_Dec19.txt'

        Commands.BasicT2WS(
            card,
            manualMaps=[
                '--PO \'map=.*/ggH_PTH_200_350:r_ggH_PTH_200_350[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_350_600:r_ggH_PTH_350_600[1.0,0.0,10.0]\'',
                '--PO \'map=.*/ggH_PTH_GT600:r_ggH_PTH_GT600[1.0,0.0,10.0]\'',
                ],
            outputWS = basename(card).replace( '.txt', '_xHfixed.root' )
            )

    #____________________________________________________________________
    if args.hbbOnly_scan:

        ws = join( base, 'out/workspaces_Dec19/combinedCard_hbb_debuggingTest_Dec19_xHfixed.root' )

        cmd = [            
            'combine',
            ws,
            '-n _debugging_hbbOnly_scan01',
            '-M MultiDimFit',
            # '--cminDefaultMinimizerType Minuit2',
            # '--cminDefaultMinimizerAlgo migrad',
            '--algo=grid',
            '-P "r_ggH_PTH_350_600"',
            '--setPhysicsModelParameterRanges r_ggH_PTH_350_600=0.9,1.1 ',
            '--setPhysicsModelParameters r_ggH_PTH_200_350=1.0,r_ggH_PTH_GT600=1.0,r_ggH_PTH_350_600=1.0',
            '-m 125.00',
            # '--squareDistPoi',
            '--saveNLL',
            '--saveInactivePOI 1',
            '--points=2 ',
            '-t -1',
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            # '--saveSpecifiedFunc qcdeff',
            '--floatOtherPOIs=1',
            '--saveWorkspace',
            ]

        Commands.BasicGenericCombineCommand( cmd, onBatch = False, )

    #____________________________________________________________________
    if args.hbb_bestfit:

        # ws = join( base, 'out/workspaces_Dec15/combinedCard_hgg_hzz_hbb_ggHxH_Dec15_xHfixed.root' )
        # ws = join( base, 'out/workspaces_Dec21/combinedCard_hgg_hzz_hbb_ggHxH_Dec21_xHfixed.root' )
        # ws = join( base, 'out/workspaces_Jan16/combinedCard_hgg_hzz_hbb_ggHxH_Dec21_xHfixed.root' )
        ws = join( base, 'out/workspaces_Jan19/combinedCard_hgg_hzz_hbb_ggHxH_Jan19_xHfixed.root' )
        
        cmd = [            
            'combine',
            ws,
            '-n _debugging_{0}_hbb_bestfit'.format(datestr),
            '-M MultiDimFit',
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            # '--algo=grid',
            # '-P r_ggH_PTH_350_600',
            # '-P r_ggH_PTH_GT600',
            # '--setPhysicsModelParameterRanges "r_ggH_PTH_350_600"=0.000,6.500 ',
            # 
            # '--setPhysicsModelParameterRanges r1p0=0,5:r2p0=0,5:r3p0=0,5:r0p1=0,5:r1p1=0,5:r2p1=0,5:r3p1=0,5 ',
            '--setPhysicsModelParameterRanges qcdeff=0.0001,10.0:r1p0=0,5:r2p0=0,5:r3p0=0,5:r0p1=0,5:r1p1=0,5:r2p1=0,5:r3p1=0,5 ',
            # 
            '--setPhysicsModelParameters r_ggH_PTH_0_15=1.0,r_ggH_PTH_200_350=1.0,r_ggH_PTH_85_125=1.0,r_ggH_PTH_15_30=1.0,r_ggH_PTH_GT600=1.0,r_ggH_PTH_350_600=1.0,r_ggH_PTH_30_45=1.0,r_ggH_PTH_125_200=1.0,r_ggH_PTH_45_85=1.0',
            '-m 125.00',
            '--squareDistPoi',
            '--saveNLL',
            '--saveInactivePOI 1',
            # '--points=2 ',
            '-t -1',
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            '--saveSpecifiedFunc qcdeff,r1p0,r2p0,r3p0,r0p1,r1p1,r2p1,r3p1',
            '--saveWorkspace',
            '--floatOtherPOIs=1',
            ]

        Commands.BasicGenericCombineCommand( cmd, onBatch = False, )


    #____________________________________________________________________
    if args.hbb_scan:

        # ws = base + 'out/workspaces_Dec15/combinedCard_hgg_hzz_hbb_ggHxH_Dec15_xHfixed.root' 
        # postfit = join( base, 'higgsCombine_debugging_hbb_bestfit_noFloatOtherPOIs.MultiDimFit.mH125.root' )
        # postfit = join( base, 'higgsCombine_debugging_Dec21_hbb_bestfit.MultiDimFit.mH125.root' )
        postfit = join( base, 'higgsCombine_debugging_{0}_hbb_bestfit.MultiDimFit.mH125.root'.format(datestr) )

        cmd = [            
            'combine',
            postfit,
            '-n _debugging_{0}_hbb_scan02'.format(datestr),
            '-M MultiDimFit',
            '--snapshotName MultiDimFit',
            '--skipInitialFit',
            # '--cminDefaultMinimizerType Minuit2',
            # '--cminDefaultMinimizerAlgo migrad',
            '--algo=grid',
            '-P "r_ggH_PTH_350_600"',
            # 
            # '--setPhysicsModelParameterRanges r_ggH_PTH_350_600=0.0,2.0:r1p0=0,5:r2p0=0,5:r3p0=0,5:r0p1=0,5:r1p1=0,5:r2p1=0,5:r3p1=0,5 ',
            '--setPhysicsModelParameterRanges qcdeff=0.0001,10.0:r_ggH_PTH_350_600=0.0,2.0:r1p0=0,5:r2p0=0,5:r3p0=0,5:r0p1=0,5:r1p1=0,5:r2p1=0,5:r3p1=0,5 ',
            # 
            '--setPhysicsModelParameters r_ggH_PTH_0_15=1.0,r_ggH_PTH_200_350=1.0,r_ggH_PTH_85_125=1.0,r_ggH_PTH_15_30=1.0,r_ggH_PTH_GT600=1.0,r_ggH_PTH_350_600=1.0,r_ggH_PTH_30_45=1.0,r_ggH_PTH_125_200=1.0,r_ggH_PTH_45_85=1.0',
            '-m 125.00',
            # '--squareDistPoi',
            '--saveNLL',
            '--saveInactivePOI 1',
            '--points=30 ',
            '-t -1',
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            # 
            # '--saveSpecifiedFunc qcdeff',
            '--saveSpecifiedFunc qcdeff,r1p0,r2p0,r3p0,r0p1,r1p1,r2p1,r3p1',
            # '--freezeNuisances qcdeff,r1p0,r2p0,r3p0,r0p1,r1p1,r2p1,r3p1',
            # '-v 10',
            ]
            # --> Note floatOtherPOIs is missing, which really should be there for a scan
            #     (otherwise all the r's not in -P are not profiled)

        Commands.BasicGenericCombineCommand( cmd, onBatch = False, )

        # Literal command from svn:
        # combine comb_2017_ggHbb.root
        # --cminDefaultMinimizerType Minuit2
        # --cminDefaultMinimizerAlgo migrad
        # --algo=grid
        # --floatOtherPOIs=1
        # -P r_ggH_PTH_350_600
        # --setPhysicsModelParameters r_ggH_PTH_GT600=1.0,r_ggH_PTH_350_600=1.0
        # --squareDistPoi
        # --saveNLL
        # --saveInactivePOI 1
        # -t
        # -1
        # --minimizerStrategy 2
        # --minimizerTolerance 0.001
        # --robustFit 1
        # --minimizerAlgoForMinos Minuit2,Migrad
        # -M MultiDimFit
        # -m 125.00
        # --setPhysicsModelParameterRanges r_ggH_PTH_350_600=-10.000,10.000


    #____________________________________________________________________
    if args.hbb_print:

        with Commands.OpenRootFile( 'higgsCombine_debugging_Jan15_hbb_bestfit.MultiDimFit.mH125.root' ) as rootFp:
            w = rootFp.Get('w')


        print 'r*p* at DEFAULT values'

        w.function( 'qcd_pass_cat1_Bin21_hbb_cat1_pass_cat1' ).Print()

        print '\nSome other vars:'
        w.function( 'Var_RhoPol_Bin_495.0_-1.985_hbb_cat1_pass_cat1' ).Print()
        w.function( 'Var_Pol_Bin_495.0_-1.985_0_hbb_cat1_pass_cat1' ).Print()
        w.function( 'Var_Pol_Bin_495.0_-1.985_1_hbb_cat1_pass_cat1' ).Print()
        w.function( 'Var_Pol_Bin_495.0_-1.985_2_hbb_cat1_pass_cat1' ).Print()
        w.function( 'Var_Pol_Bin_495.0_-1.985_3_hbb_cat1_pass_cat1' ).Print()
        w.function( 'Var_Rho_rescaled_495.0_-1.985' ).Print()


        print '\n\nr*p* at STRANGE values'

        reset = {
            'qcdeff' : 8.85844,
            'r0p1'   : 3.34228,
            'r1p0'   : 4.09749,
            'r1p1'   : 4.97966,
            'r2p0'   : 0.00362817,
            'r2p1'   : 0.0380896,
            'r3p0'   : 0.642744,
            'r3p1'   : 0.0384077,
            # 
            'qcd_fail_cat1_Bin1'  : 65041.5,
            'qcd_fail_cat1_Bin10' : 2563.44,
            'qcd_fail_cat1_Bin11' : 119152,
            'qcd_fail_cat1_Bin12' : 134176,
            'qcd_fail_cat1_Bin13' : 273.105,
            'qcd_fail_cat1_Bin14' : 113898,
            'qcd_fail_cat1_Bin15' : 148850,
            'qcd_fail_cat1_Bin16' : 45569.4,
            'qcd_fail_cat1_Bin17' : 85297.5,
            'qcd_fail_cat1_Bin18' : 49202.3,
            'qcd_fail_cat1_Bin19' : 6533.66,
            'qcd_fail_cat1_Bin2'  : 314905,
            'qcd_fail_cat1_Bin20' : 38023.4,
            'qcd_fail_cat1_Bin21' : 1766.74,
            'qcd_fail_cat1_Bin22' : 8783.92,
            'qcd_fail_cat1_Bin23' : 16825.1,
            'qcd_fail_cat1_Bin3'  : 9330.23,
            'qcd_fail_cat1_Bin4'  : 296098,
            'qcd_fail_cat1_Bin5'  : 7002.19,
            'qcd_fail_cat1_Bin6'  : 143492,
            'qcd_fail_cat1_Bin7'  : 189575,
            'qcd_fail_cat1_Bin8'  : 223431,
            'qcd_fail_cat1_Bin9'  : 75448.3,
            'qcd_fail_cat2_Bin1'  : 0.265255,
            'qcd_fail_cat2_Bin10' : 135.548,
            'qcd_fail_cat2_Bin11' : 38360.7,
            'qcd_fail_cat2_Bin12' : 1445.61,
            'qcd_fail_cat2_Bin13' : 1103.76,
            'qcd_fail_cat2_Bin14' : 26853.1,
            'qcd_fail_cat2_Bin15' : 33561.8,
            'qcd_fail_cat2_Bin16' : 26543.1,
            'qcd_fail_cat2_Bin17' : 37945,
            'qcd_fail_cat2_Bin18' : 24385.4,
            'qcd_fail_cat2_Bin19' : 168.786,
            'qcd_fail_cat2_Bin2'  : 34.5315,
            'qcd_fail_cat2_Bin20' : 25823.3,
            'qcd_fail_cat2_Bin21' : 5739.9,
            'qcd_fail_cat2_Bin22' : 25605.6,
            'qcd_fail_cat2_Bin23' : 13501.4,
            'qcd_fail_cat2_Bin3'  : 60362,
            'qcd_fail_cat2_Bin4'  : 42123.9,
            'qcd_fail_cat2_Bin5'  : 54895.5,
            'qcd_fail_cat2_Bin6'  : 21872.1,
            'qcd_fail_cat2_Bin7'  : 38219.6,
            'qcd_fail_cat2_Bin8'  : 50236.4,
            'qcd_fail_cat2_Bin9'  : 1730.5,
            }

        for varName, value in reset.iteritems():
            w.var(varName).setVal(value)


        for i in xrange(23):
            w.function( 'qcd_pass_cat1_Bin{0}_hbb_cat1_pass_cat1'.format(i+1) ).Print()

        # w.function( 'qcd_pass_cat1_Bin21_hbb_cat1_pass_cat1' ).Print()

        print '\nSome other vars:'
        w.function( 'Var_RhoPol_Bin_495.0_-1.985_hbb_cat1_pass_cat1' ).Print()
        w.function( 'Var_Pol_Bin_495.0_-1.985_0_hbb_cat1_pass_cat1' ).Print()
        w.function( 'Var_Pol_Bin_495.0_-1.985_1_hbb_cat1_pass_cat1' ).Print()
        w.function( 'Var_Pol_Bin_495.0_-1.985_2_hbb_cat1_pass_cat1' ).Print()
        w.function( 'Var_Pol_Bin_495.0_-1.985_3_hbb_cat1_pass_cat1' ).Print()
        w.function( 'Var_Rho_rescaled_495.0_-1.985' ).Print()


    #____________________________________________________________________
    if args.testing2PanelCanvas:

        TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )

        Tg_top = ROOT.TGraph( 2, array( 'f', [ 0.1, 0.9 ] ), array( 'f', [ 0.1, 0.9 ] ) )
        ROOT.SetOwnership( Tg_top, False )

        Tg_bottom = ROOT.TGraph( 2, array( 'f', [ 0.1, 0.9 ] ), array( 'f', [ 0.1, 0.9 ] ) )
        ROOT.SetOwnership( Tg_bottom, False )

        PlotCommands.PlotWithBottomPanel(
            'twopaneltest',
            [ ( Tg_top, 'AL' ) ],
            [ ( Tg_bottom, 'AL' ) ],
            xTitle = 'x',
            yTitleTop = 'y_{top}',
            yTitleBottom = 'y_{bottom}',
            )

    #____________________________________________________________________
    if args.testingDoPoints:

        base = '/mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/v2_NNLOPS/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/'
        ws   = base + 'postfitWSs_Nov22/POSTFIT_hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_Yukawa_withTheoryUncertainties.root'

        cmd = [
            'combine',
            ws,
            '-n _debugging_Nov22_doPoints',
            '--algo=grid',
            '--points=64',
            '-P kappab -P kappac',
            '--setPhysicsModelParameters kappab=1.0,kappac=1.0',
            '--setPhysicsModelParameterRanges kappab=-15.0,15.0:kappac=-35.0,35.0',
            '-M MultiDimFit',
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            '-m 125.0',
            '--saveNLL',
            '--saveInactivePOI 1',
            '--floatOtherPOIs=1',
            '--doPoints 8,14,20,21,45'
            ]

        Commands.BasicGenericCombineCommand(
            cmd,
            onBatch = False,
            )

    #____________________________________________________________________
    if args.kappaVMaxOne_Top:

        base = '/mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/v2_NNLOPS/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/'
        ws = base + 'workspaces_Nov10/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_couplingDependentBR.root'

        cmd = [
            'combine',
            ws,
            '-n debugging_{0}_kappaVMaxOne_Top'.format(datestr),
            '-M MultiDimFit',
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            '--floatOtherPOIs=1',
            '-m 125.00',
            '--saveNLL',
            '--saveInactivePOI 1',
            '--squareDistPoiStep',
            '-P ct -P cg',
            '--setPhysicsModelParameters ct=1.0,cg=0.0,kappa_V=0.99',
            '--setPhysicsModelParameterRanges ct=-1.0,2.0:cg=-0.1,0.2:kappa_V=-100.0,1.0',
            '--floatNuisances kappa_V',
            '--saveSpecifiedFunc r_ggH_PTH_0_15,r_ggH_PTH_15_30,r_ggH_PTH_30_45,r_ggH_PTH_45_85,r_ggH_PTH_85_125,r_ggH_PTH_125_200,r_ggH_PTH_200_350,r_ggH_PTH_GT350,theoryUncertainty_0,theoryUncertainty_1,theoryUncertainty_2,theoryUncertainty_3,theoryUncertainty_4,theoryUncertainty_5,theoryUncertainty_6,r_hgg_ggH_PTH_0_15,r_hgg_ggH_PTH_15_30,r_hgg_ggH_PTH_30_45,r_hgg_ggH_PTH_45_85,r_hgg_ggH_PTH_85_125,r_hgg_ggH_PTH_125_200,r_hgg_ggH_PTH_200_350,r_hgg_ggH_PTH_GT350,r_hgg_ggH_PTH_0_15,r_hgg_ggH_PTH_15_30,r_hgg_ggH_PTH_30_45,r_hgg_ggH_PTH_45_85,r_hgg_ggH_PTH_85_125,r_hgg_ggH_PTH_125_200,r_hgg_ggH_PTH_200_350,r_hgg_ggH_PTH_GT350,ct,kappa_b,kappa_c,kappa_V,kappa_tau,kappa_mu,hggBRmodifier,hzzBRmodifier,xH_modifier,Scaling_hgg,Scaling_hzg,Scaling_hgluglu',
            ]

        Commands.BasicGenericCombineCommand(
            cmd,
            onBatch = False,
            )


    #____________________________________________________________________
    if args.jscaletest_RenameHggProcesses_alsoSystematics:

        # LatestPaths.card_hgg_smH_NJ_unprocessed
        # LatestPaths.card_hzz_smH_NJ

        MergeHGGWDatacards.RenameProcesses_Hgg_nJets(
            LatestPaths.card_hgg_smH_NJ_unprocessed,
            outTag = '_debugging_Nov10',
            globalReplace = [
                ( 'CMS_hgg_JER', 'CMS_scale_j' )
                ]
            )

        hgg_debugging_out = 'suppliedInput/fromVittorio/differential_Njets2p5NNLOPS_Nov10/Datacard_13TeV_differential_Njets2p5NNLOPS_debugging_Nov10.txt'

        Commands.BasicCombineCards(
            'suppliedInput/combinedCard_smH_debugging_Nov10.txt',
            'hgg=' + hgg_debugging_out,
            'hzz=' + LatestPaths.card_hzz_smH_NJ
            )


    #____________________________________________________________________
    if args.jscaletest_t2ws:

        card = 'suppliedInput/combinedCard_smH_debugging_Nov10.txt'
        ws   = card.replace( '.txt', '.root' )

        Commands.BasicT2WS(
            card,
            smartMaps = [
                ( r'.*/smH_NJ_([\d\_GE]+)', r'r_smH_NJ_\1[1.0,-1.0,4.0]' )
                ],
            )

    #____________________________________________________________________
    if args.jscaletest_bestfit:

        ws = abspath( 'workspaces_Nov10/combinedCard_smH_debugging_Nov10.root' )

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
            #     Commands.ListSet( datacard, 'yieldParameters' ) + [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
            '--squareDistPoiStep',
            '-M MultiDimFit',
            '-m 125.00',
            # '--setPhysicsModelParameterRanges kappab=-20.0,20.0:kappac=-50.0,50.0',
            # '--setPhysicsModelParameterRanges kappab=0.5,1.0:kappac=1.0,2.0',
            # '--points 12800',
            # '--firstPoint 0',
            # '--lastPoint 799',
            '-n debugging_Nov10_jscaletest',
            # '-v 3',
            ]

        Commands.BasicGenericCombineCommand(
            cmd,
            onBatch = False,
            )


    #____________________________________________________________________
    if args.debugging_Nov08_TopWS:

        base = '/mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/v2_NNLOPS/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/'

        # ws = base + 'workspaces_Nov08/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties.root'
        # ws = base + 'workspaces_Nov08/combinedCard_Nov03_CouplingModel_Top_noTheoryUncertainties.root'

        ws = base + 'workspaces_Nov08/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties.root'

        cmd = [
            'combine',
            ws,
            '-n _debugging_Nov08_TopWS',
            '-M MultiDimFit',
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            # '--algo=grid',
            '--floatOtherPOIs=1',
            '-m 125.00',
            '--saveNLL',
            '--saveInactivePOI 1',
            # '--points=5 ',
            # '-t -1',
            '-P ct -P cg',
            '--squareDistPoiStep',
            '--setPhysicsModelParameters ct=1.0,cg=0.0',
            # '--setPhysicsModelParameterRanges ct=-1.0,2.0:cg=-0.1,0.2',
            '--saveSpecifiedFunc r_ggH_PTH_0_15,r_ggH_PTH_15_30,r_ggH_PTH_30_45,r_ggH_PTH_45_85,r_ggH_PTH_85_125,r_ggH_PTH_125_200,r_ggH_PTH_200_350,r_ggH_PTH_GT350'
            + ( ',theoryUncertainty_0,theoryUncertainty_1,theoryUncertainty_2,theoryUncertainty_3,theoryUncertainty_4,theoryUncertainty_5,theoryUncertainty_6,theoryUncertainty_7'
                if False else '' ),
            ]

        Commands.BasicGenericCombineCommand(
            cmd,
            onBatch = False,
            )

    #____________________________________________________________________
    if args.debugging_Nov08_hzz:

        base = '/mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/v2_NNLOPS/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/'

        # ws = base + 'workspaces_Nov08/hzz4l_comb_13TeV_xs.root'
        ws = base + 'workspaces_Nov08/hzz4l_comb_13TeV_xs_processesShifted.root'

        cmd = [
            'combine',
            ws,
            '-n _debugging_Nov08_hzz',
            '-M MultiDimFit',
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            '--algo=grid',
            '--floatOtherPOIs=1',
            '-P "r_smH_PTH_0_15"',
            '--setPhysicsModelParameterRanges r_smH_PTH_0_15=0.700,2.000 ',
            '--setPhysicsModelParameters r_smH_PTH_0_15=1.0,r_smH_PTH_85_200=1.0,r_smH_PTH_15_30=1.0,r_smH_PTH_30_85=1.0,r_smH_PTH_GT200=1.0',
            '-m 125.00',
            '--squareDistPoi',
            '--saveNLL',
            '--saveInactivePOI 1',
            '--points=10',
            # 
            # '--freezeNuisances r_smH_PTH_GT200',
            ]


        Commands.BasicGenericCombineCommand(
            cmd,
            onBatch = False,
            )


    if args.BRdependency_Yukawa:

        base = '/mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/'

        ws = 'workspaces_Oct02/combinedCard_Aug21_CouplingModel_Yukawa_withTheoryUncertainties_couplingDependentBR.root'

        cmd = [
            'combine',
            base + ws,
            '-M MultiDimFit',
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            # '--algo=grid',
            '--floatOtherPOIs=1',
            '-m 125.00',
            '--saveNLL',
            '--saveInactivePOI 1',
            # '-t -1',
            '-P kappab -P kappac',
            '--floatNuisances kappa_V',
            '--squareDistPoiStep',
            '--setPhysicsModelParameters kappab=1.0,kappac=1.0,kappa_V=1.0',
            '--setPhysicsModelParameterRanges kappab=-20.0,20.0:kappac=-50.0,50.0:kappa_V=-10.0,10.0',
            '--saveSpecifiedFunc r_ggH_PTH_0_15,r_ggH_PTH_15_30,r_ggH_PTH_30_45,r_ggH_PTH_45_85,r_ggH_PTH_85_125,theoryUncertainty_0,theoryUncertainty_1,theoryUncertainty_2,theoryUncertainty_3,theoryUncertainty_4,r_hgg_ggH_PTH_0_15,r_hgg_ggH_PTH_15_30,r_hgg_ggH_PTH_30_45,r_hgg_ggH_PTH_45_85,r_hgg_ggH_PTH_85_125,r_hgg_ggH_PTH_0_15,r_hgg_ggH_PTH_15_30,r_hgg_ggH_PTH_30_45,r_hgg_ggH_PTH_45_85,r_hgg_ggH_PTH_85_125,kappa_t,kappab,kappac,kappa_V,kappa_tau,kappa_mu,hggBRmodifier,hzzBRmodifier,xHmodifier,Scaling_hgg,Scaling_hzg,Scaling_hgluglu',
            '-n debugging_{0}_BRdependency_Yukawa'.format(datestr),
            ]

        Commands.BasicGenericCombineCommand(
            cmd,
            onBatch = False,
            )


    if args.BRdependency:

        base = '/mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/'

        # ws = 'workspaces_Sep29/combinedCard_Aug21_CouplingModel_Top_withTheoryUncertainties_couplingDependentBR.root'
        ws = 'workspaces_Oct02/combinedCard_Aug21_CouplingModel_Top_withTheoryUncertainties_couplingDependentBR.root'

        cmd = [
            'combine',
            base + ws,
            '-M MultiDimFit',
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            '--floatOtherPOIs=1',
            '-m 125.00',
            '--saveNLL',
            '--saveInactivePOI 1',
            # '-t -1',
            '-P ct -P cg',
            '--floatNuisances kappa_V'
            '--squareDistPoiStep',
            '--setPhysicsModelParameters ct=1.0,cg=0.0,kappa_V=1.0',
            '--setPhysicsModelParameterRanges ct=-1.0,2.0:cg=-0.1,0.2:kappa_V=-10.0,10.0',
            '--saveSpecifiedFunc r_ggH_PTH_0_15,r_ggH_PTH_15_30,r_ggH_PTH_30_45,r_ggH_PTH_45_85,r_ggH_PTH_85_125,r_ggH_PTH_125_200,r_ggH_PTH_200_350,r_ggH_PTH_GT350,theoryUncertainty_0,theoryUncertainty_1,theoryUncertainty_2,theoryUncertainty_3,theoryUncertainty_4,theoryUncertainty_5,theoryUncertainty_6,r_hgg_ggH_PTH_0_15,r_hgg_ggH_PTH_15_30,r_hgg_ggH_PTH_30_45,r_hgg_ggH_PTH_45_85,r_hgg_ggH_PTH_85_125,r_hgg_ggH_PTH_125_200,r_hgg_ggH_PTH_200_350,r_hgg_ggH_PTH_GT350,Scaling_hgg,kappa_V',
            '--saveWorkspace',
            '-n debugging_BRdependency',
            ]

        Commands.BasicGenericCombineCommand(
            cmd,
            onBatch = False,
            )


    if args.testfit_njets:

        ws = abspath( 'workspaces_Sep19/combinedCard_nJets_Sep19.root' )

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
            #     Commands.ListSet( datacard, 'yieldParameters' ) + [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
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

        Commands.BasicGenericCombineCommand(
            cmd,
            onBatch = False,
            )


    if args.simplebestfit:

        Commands.SetTempJobDir( 'plainWStests_{0}'.format(datestr) )

        datacard = ( '/mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/'
              # 'workspaces_Sep14/hzz4l_all_13TeV_xs_processesShifted_CouplingModel_Yukawa_withTheoryUncertainties.root'
              # 'workspaces_Sep15/hzz4l_all_13TeV_xs_processesShifted_CouplingModel_Yukawa_withTheoryUncertainties.root'
              # 'workspaces_Sep15/hzz4l_all_13TeV_xs_processesShifted_CouplingModel_Yukawa.root'
              'workspaces_Sep18/hzz4l_all_13TeV_xs_processesShifted_CouplingModel_Yukawa_withTheoryUncertainties.root'
              )

        cmd = [
            'combine',
            datacard,
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
            '--saveSpecifiedFunc {0}'.format(','.join(
                Commands.ListSet( datacard, 'yieldParameters' ) + [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
            '--squareDistPoiStep',
            '-M MultiDimFit',
            '-m 125.00',
            '--setPhysicsModelParameterRanges kappab=-20.0,20.0:kappac=-50.0,50.0',
            # '--setPhysicsModelParameterRanges kappab=0.5,1.0:kappac=1.0,2.0',
            # '--points 12800',
            # '--firstPoint 0',
            # '--lastPoint 799',
            '-n testjob',
            # '-v 3',
            ]

        Commands.BasicGenericCombineCommand(
            cmd,
            onBatch = False,
            )

    if args.asimovtest:

        Commands.SetTempJobDir( 'plainWStests_{0}'.format(datestr) )

        datacard = LatestPaths.ws_combined_split_yukawa

        cmd = [
            'combine',
            datacard,
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            # '--algo=grid',
            '--floatOtherPOIs=1',
            '--saveNLL',
            '--saveInactivePOI 1',
            '-t -1',
            # '--fastScan',
            # '-P kappab',
            # '-P kappac',
            # '--setPhysicsModelParameters kappab=1.0,kappac=1.0',
            '--saveSpecifiedFunc {0}'.format(','.join(
                Commands.ListSet( datacard, 'yieldParameters' ) + [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
            '--squareDistPoiStep',
            '-M MultiDimFit',
            '-m 125.00',
            '--setPhysicsModelParameterRanges kappab=-20.0,20.0:kappac=-50.0,50.0',
            # '--setPhysicsModelParameterRanges kappab=0.5,1.0:kappac=1.0,2.0',
            # '--points 12800',
            # '--firstPoint 0',
            # '--lastPoint 799',
            '-n testjob',
            # '-v 3',
            ]

        Commands.BasicGenericCombineCommand(
            cmd,
            onBatch = False,
            )


    if args.toptest:

        # datacard = ( '/mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/'
        #     'workspaces_Aug22/hzz4l_all_13TeV_xs_processesShifted_CouplingModel_Top_withTheoryUncertainties.root'
        #     )

        datacard = 'workspaces_Sep27/combinedCard_Aug21_CouplingModel_Top_withTheoryUncertainties_couplingDependentBR.root'

        cmd = [
            'combine',
            datacard,
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            # '--algo=grid',
            '--floatOtherPOIs=1',
            '--saveNLL',
            '--saveInactivePOI 1',
            '-P ct',
            '-P cg',
            '--setPhysicsModelParameters ct=1.0,cg=0.0',
            '--saveSpecifiedFunc {0}'.format( ','.join(
                Commands.ListSet( datacard, 'yieldParameters' )
                + Commands.ListSet( datacard, 'hgg_yieldParameters' )
                + [ 'Scaling_hgg' ]
                + [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]
                ) ),
            '--squareDistPoiStep',
            '-M MultiDimFit',
            '-m 125.00',
            '--setPhysicsModelParameterRanges ct=-1.0,2.0:cg=-0.1,0.2',
            # '--points 6400',
            # '--firstPoint 0',
            # '--lastPoint 79',
            '-n testjob',
            ]

        Commands.BasicGenericCombineCommand(
            cmd,
            onBatch = False,
            )




########################################
# End of Main
########################################
if __name__ == "__main__":
    main()
