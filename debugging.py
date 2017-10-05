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

import LatestPaths

from time import strftime
datestr = strftime( '%b%d' )


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
