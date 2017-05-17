#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, itertools
from os.path import *
from glob import glob


import Commands
import PhysicsCommands
import OneOfCommands


from time import strftime
datestr = strftime( '%b%d' )


########################################
# Main
########################################

def main():

    # Commands.TestMode()


    # ======================================
    # Separate tests

    # processes, bins = Commands.ListProcesses( 'suppliedInput/PTH/hzz4l_comb_13TeV_xs.txt' )
    # print 'Processes ({0}):\n    '.format(len(processes)) + '\n    '.join( processes )
    # print 'Bins ({0}):\n    '.format(len(bins)) + '\n    '.join( bins )


    # ======================================
    # Changing OutsideAcceptance from background to a signal
    # OneOfCommands.ChangeOutsideAcceptanceToSignalProcess()



    # ======================================
    # combineCards.py

    # Commands.BasicCombineCards(
    #     'suppliedInput/combinedCard_{0}.txt'.format(datestr),
    #     'suppliedInput/pT_v5_renamed/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.txt',
    #     'suppliedInput/PTH/hzz4l_comb_13TeV_xs.txt'
    #     )


    # Commands.BasicCombineCards(
    #     'suppliedInput/combinedCard_{0}.txt'.format(datestr),
    #     'suppliedInput/pT_v5_renamed/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.txt',
    #     'suppliedInput/PTH/hzz4l_comb_13TeV_xs_processesShifted.txt'
    #     )


    # ======================================
    # text2workspace

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

            # Tg.SetLineColor(color)
            # Tg.SetMarkerColor(color)
            # Tg.SetLineWidth(2)
            # Tg.SetMarkerStyle(5)
        )




########################################
# End of Main
########################################
if __name__ == "__main__":
    main()