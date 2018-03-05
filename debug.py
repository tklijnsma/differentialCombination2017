#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

from OptionHandler import flag_as_option

import LatestPaths
import LatestPathsGetters
import LatestBinning

import differentials

import logging
import copy
import random
random.seed(1002)

import sys
sys.path.append('src')
import TheoryFileInterface

from time import strftime
datestr = strftime('%b%d')

########################################
# Main
########################################


@flag_as_option
def debug_theoryfiles_top(args):
    differentials.theory.ctcg_interpreter.create_all_ctcg()


@flag_as_option
def debug_Scan_hbb(args):
    cmd = [
        'combine',
        # 'out/workspaces_Feb28/ws_pth_ggH_hbb.root',
        'out/postfits_Mar02/higgsCombine_POSTFIT_ws_pth_ggH_hbb.MultiDimFit.mH125.root',
        # 
        '-n _DEBUG_{0}_SCAN_bPOI_r_ggH_PTH_350_600_ePOI_ws_pth_ggH_hbb'.format(datestr),
        '-M MultiDimFit',
        '-m 125.0',
        # '--cminDefaultMinimizerType Minuit2',
        # '--cminDefaultMinimizerAlgo migrad',
        '--minimizerStrategy 2',
        '--minimizerTolerance 0.001',
        '--robustFit 1',
        '--minimizerAlgoForMinos Minuit2,Migrad',
        '--saveNLL',
        '--saveInactivePOI 1',
        '--floatOtherPOIs=1',
        '-P r_ggH_PTH_350_600',
        '--setPhysicsModelParameterRanges r_ggH_PTH_200_350=-10.0,10.0:r_ggH_PTH_350_600=-10.0,10.0:r_ggH_PTH_GT600=-10.0,10.0',
        # :qcdeff=0.001,8.0:r1p0=0.0,8.0:r2p0=0.0,8.0:r3p0=0.0,8.0:r0p1=0.0,8.0:r1p1=0.0,8.0:r2p1=0.0,8.0:r3p1=0.0,8.0',
        # '--setPhysicsModelParameters r_ggH_PTH_200_350=1.0,r_ggH_PTH_350_600=1.0,r_ggH_PTH_GT600=1.0',
        '--algo=grid',
        '--points=3',
        # 
        '--snapshotName MultiDimFit',
        '--skipInitialFit',
        '-v 2',
        ]
    differentials.core.execute(cmd)


@flag_as_option
def debug_combWithHbb_Scan_Mar02(args):
    # ws = 'out/workspaces_Mar01/ws_pth_ggH_combWithHbb.root'

    # ws = 'out/postfits_Mar02_v1/higgsCombine_POSTFIT_ws_pth_ggH_combination.MultiDimFit.mH125.root'
    # ws = 'out/postfits_Mar02_v1/higgsCombine_POSTFIT_ASIMOV_ws_pth_ggH_combination.MultiDimFit.mH125.root'
    # ws = 'out/postfits_Mar02_v1/higgsCombine_POSTFIT_ws_pth_ggH_combWithHbb.MultiDimFit.mH125.root'
    # ws = 'out/postfits_Mar02_v1/higgsCombine_POSTFIT_ASIMOV_ws_pth_ggH_combWithHbb.MultiDimFit.mH125.root'

    ws = 'out/postfits_Mar02/higgsCombine_POSTFIT_ASIMOV_ws_pth_ggH_combWithHbb.MultiDimFit.mH125.root'
    # ws = 'out/postfits_Mar02/higgsCombine_POSTFIT_ASIMOV_ws_pth_ggH_hbb.MultiDimFit.mH125.root'
    # ws = 'out/postfits_Mar02/higgsCombine_POSTFIT_ws_pth_ggH_combWithHbb.MultiDimFit.mH125.root'
    # ws = 'out/postfits_Mar02/higgsCombine_POSTFIT_ws_pth_ggH_hbb.MultiDimFit.mH125.root'

    cmd = [
        'combine',
        ws,
        '-n _DEBUG_{0}_SCAN'.format(datestr),
        '-M MultiDimFit',
        '-m 125.0',
        # 
        '--minimizerStrategy 2',
        '--minimizerTolerance 0.001',
        '--robustFit 1',
        '--minimizerAlgoForMinos Minuit2,Migrad',
        # 
        '--saveNLL',
        '--saveInactivePOI 1',
        '--floatOtherPOIs=1',
        # '--job-mode psi --task-name _SCAN_bPOI_r_ggH_PTH_600_10000_ePOI_ws_pth_ggH_combWithHbb --sub-opts='-q short.q' ',
        # 
        # '-P r_ggH_PTH_GT600',
        '-P r_ggH_PTH_120_200',
        # 
        '--setPhysicsModelParameterRanges r_ggH_PTH_0_15={0},{1}:r_ggH_PTH_15_30={0},{1}:r_ggH_PTH_30_45={0},{1}:r_ggH_PTH_45_80={0},{1}:r_ggH_PTH_80_120={0},{1}:r_ggH_PTH_120_200={0},{1}:r_ggH_PTH_200_350={0},{1}:r_ggH_PTH_350_600={0},{1}:r_ggH_PTH_GT600={0},{1}:qcdeff=0.001,8.0:r1p0=0.0,8.0:r2p0=0.0,8.0:r3p0=0.0,8.0:r0p1=0.0,8.0:r1p1=0.0,8.0:r2p1=0.0,8.0:r3p1=0.0,8.0'.format(
            0.0, 8.0
            ),
        '--setPhysicsModelParameters r_ggH_PTH_0_15=1.0,r_ggH_PTH_15_30=1.0,r_ggH_PTH_30_45=1.0,r_ggH_PTH_45_80=1.0,r_ggH_PTH_80_120=1.0,r_ggH_PTH_120_200=1.0,r_ggH_PTH_200_350=1.0,r_ggH_PTH_350_600=1.0,r_ggH_PTH_GT600=1.0',
        '--algo=grid',
        '--points=3',
        # '--split-points 2',
        '--snapshotName MultiDimFit',
        '--skipInitialFit',
        '-v 2',
        # 
        # '--freezeNuisances qcdeff,r1p0,r2p0,r3p0,r0p1,r1p1,r2p1,r3p1',
        # 
        '--cminDefaultMinimizerType Minuit2',
        '--cminDefaultMinimizerAlgo migrad',
        # 
        '-t -1'
        ]
    differentials.core.execute(cmd)

    # qcdeff r1p0 r2p0 r3p0 r0p1 r1p1 r2p1 r3p1

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
    # -t -1
    # --minimizerStrategy 2
    # --minimizerTolerance 0.001
    # --robustFit 1
    # --minimizerAlgoForMinos Minuit2,Migrad
    # -M MultiDimFit
    # -m 125.00
    # --setPhysicsModelParameterRanges r_ggH_PTH_350_600=-10.000,10.000


@flag_as_option
def debug_hbb_bestfit(args):
    cmd = [
        'combine',
        'out/workspaces_Feb28/ws_pth_ggH_hbb.root',
        '-n _POSTFIT_ASIMOV_ws_pth_ggH_hbb',
        '-M MultiDimFit',
        '-m 125.0',
        '-t -1',
        # 
        '--cminDefaultMinimizerType Minuit2',
        '--cminDefaultMinimizerAlgo migrad',
        # 
        '--minimizerStrategy 2',
        '--minimizerTolerance 0.001',
        '--robustFit 1',
        '--minimizerAlgoForMinos Minuit2,Migrad',
        '--saveNLL',
        '--saveInactivePOI 1',
        '--floatOtherPOIs=1',
        '-P r_ggH_PTH_200_350 -P r_ggH_PTH_350_600 -P r_ggH_PTH_GT600',
        '--setPhysicsModelParameterRanges r_ggH_PTH_200_350=-10.0,10.0:r_ggH_PTH_350_600=-10.0,10.0:r_ggH_PTH_GT600=-10.0,10.0:qcdeff=0.001,8.0:r1p0=0.0,8.0:r2p0=0.0,8.0:r3p0=0.0,8.0:r0p1=0.0,8.0:r1p1=0.0,8.0:r2p1=0.0,8.0:r3p1=0.0,8.0',
        '--setPhysicsModelParameters r_ggH_PTH_200_350=1.0,r_ggH_PTH_350_600=1.0,r_ggH_PTH_GT600=1.0',
        '--saveWorkspace',
        ]
    differentials.core.execute(cmd)


@flag_as_option
def debug_combWithHbb_Mar01(args):
    ws = 'out/workspaces_Mar01/ws_pth_ggH_combWithHbb.root'
    cmd = [
        'combine',
        ws,
        '-n _DEBUG_SCAN_bPOI_r_ggH_PTH_600_10000_ePOI_ws_pth_ggH_combWithHbb',
        '-M MultiDimFit',
        '-m 125.0',
        '--minimizerStrategy 2',
        '--minimizerTolerance 0.001',
        '--robustFit 1',
        '--minimizerAlgoForMinos Minuit2,Migrad',
        '--saveNLL',
        '--saveInactivePOI 1',
        '--floatOtherPOIs=1',
        # '--job-mode psi --task-name _SCAN_bPOI_r_ggH_PTH_600_10000_ePOI_ws_pth_ggH_combWithHbb --sub-opts='-q short.q' ',
        '-P r_ggH_PTH_GT600',
        '--setPhysicsModelParameterRanges r_ggH_PTH_0_15=0.0,8.5:r_ggH_PTH_15_30=0.0,8.5:r_ggH_PTH_30_45=0.0,8.5:r_ggH_PTH_45_80=0.0,8.5:r_ggH_PTH_80_120=0.0,8.5:r_ggH_PTH_120_200=0.0,8.5:r_ggH_PTH_200_350=0.0,8.5:r_ggH_PTH_350_600=0.0,8.5:r_ggH_PTH_GT600=0.0,8.5:r_ggH_PTH_600_10000=0.0,8.5:qcdeff=0.001,8.0:r1p0=0.0,8.0:r2p0=0.0,8.0:r3p0=0.0,8.0:r0p1=0.0,8.0:r1p1=0.0,8.0:r2p1=0.0,8.0:r3p1=0.0,8.0',
        '--setPhysicsModelParameters r_ggH_PTH_0_15=1.0,r_ggH_PTH_15_30=1.0,r_ggH_PTH_30_45=1.0,r_ggH_PTH_45_80=1.0,r_ggH_PTH_80_120=1.0,r_ggH_PTH_120_200=1.0,r_ggH_PTH_200_350=1.0,r_ggH_PTH_350_600=1.0,r_ggH_PTH_GT600=1.0,r_ggH_PTH_600_10000=1.0',
        '--algo=grid',
        '--points=10',
        # '--split-points 2',
        ]
    differentials.core.execute(cmd)

@flag_as_option
def debug_combination_Mar01(args):
    # ws = 'out/workspaces_Mar01/ws_pth_ggH_combination.root'
    cmd = [
        'combine',
        'out/workspaces_Mar01/ws_pth_ggH_combination.root',
        '-n _DEBUG_SCAN_bPOI_r_ggH_PTH_GT600_ePOI_ws_pth_ggH_combination',
        '-M MultiDimFit',
        '-m 125.0',
        '--cminDefaultMinimizerType Minuit2',
        '--cminDefaultMinimizerAlgo migrad',
        '--saveNLL',
        '--saveInactivePOI 1',
        '--floatOtherPOIs=1',
        # '--job-mode psi --task-name _SCAN_bPOI_r_ggH_PTH_GT600_ePOI_ws_pth_ggH_combination --sub-opts='-q short.q' ',
        '-P r_ggH_PTH_GT600',
        '--setPhysicsModelParameterRanges r_ggH_PTH_0_15=0.0,4.0:r_ggH_PTH_15_30=0.0,4.0:r_ggH_PTH_30_45=0.0,4.0:r_ggH_PTH_45_80=0.0,4.0:r_ggH_PTH_80_120=0.0,4.0:r_ggH_PTH_120_200=0.0,4.0:r_ggH_PTH_200_350=0.0,4.0:r_ggH_PTH_350_600=0.0,4.0:r_ggH_PTH_GT600=0.0,4.0',
        '--setPhysicsModelParameters r_ggH_PTH_0_15=1.0,r_ggH_PTH_15_30=1.0,r_ggH_PTH_30_45=1.0,r_ggH_PTH_45_80=1.0,r_ggH_PTH_80_120=1.0,r_ggH_PTH_120_200=1.0,r_ggH_PTH_200_350=1.0,r_ggH_PTH_350_600=1.0,r_ggH_PTH_GT600=1.0',
        '--algo=grid',
        '--points=5',
        # '--split-points 5',
        ]
    differentials.core.execute(cmd)