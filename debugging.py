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
