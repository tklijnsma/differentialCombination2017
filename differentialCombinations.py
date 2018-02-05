#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, sys
from os.path import *
from glob import glob
from copy import deepcopy

from OptionHandler import flag_as_option

sys.path.append('src')
import Commands
import PhysicsCommands
import TheoryCommands
import LatestPaths
import LatestBinning
from Container import Container
import PlotCommands

from time import strftime
datestr = strftime( '%b%d' )


########################################
# Main
########################################

# Test of new bins for hzz and hbb
@flag_as_option
def RenumberHzzProcesses_Jan24(args):
    MergeHGGWDatacards.RenumberProcessesHZZ_Aug21(
        LatestPaths.card_hzz_ggHxH_PTH_newBins_unprocessed,
        )
@flag_as_option
def CombineCards_Jan24_hzz_hbb(args):
    Commands.BasicCombineCards(
        'suppliedInput/combinedCard_newBins_hzz_hbb_ggHxH_{0}.txt'.format(datestr),
        'hzz=' + LatestPaths.card_hzz_ggHxH_PTH_newBins,
        'hbb=' + LatestPaths.card_hbb_ggHxH_PTH
        )


#____________________________________________________________________

@flag_as_option
def postfit_all(real_args):
    args = deepcopy(real_args)
    def set_all_false(args):
        args.hgg = False
        args.hzz = False
        args.hbb = False
        args.combWithHbb = False

    for decay_channel in [
            'hgg',
            'hzz',
            # Workspaces are not working yet; Wait for Vs input
            # 'hbb',
            # 'combWithHbb',
            'combination'
            ]:
        set_all_false(args)
        setattr(args, decay_channel, True)
        print '\nMaking postfits for {0}'.format(decay_channel)

        try:
            pth_smH_postfit(args)
        except NotImplementedError:
            pass

        try:
            pth_ggH_postfit(args)
        except NotImplementedError:
            pass

        try:
            ptjet_postfit(args)
        except NotImplementedError:
            pass

        try:
            rapidity_postfit(args)
        except NotImplementedError:
            pass

        try:
            njets_postfit(args)
        except NotImplementedError:
            pass


#____________________________________________________________________
# pth_smH

# Preprocessing
@flag_as_option
def RenameHggProcesses_smHcard(args):
    MergeHGGWDatacards.RenameProcesses_Aug21(
        LatestPaths.card_hgg_smH_PTH_unprocessed,
        )
@flag_as_option
def RenumberHzzProcesses_smHcard(args):
    MergeHGGWDatacards.RenumberProcessesHZZ_Aug21(
        LatestPaths.card_hzz_smH_PTH_unprocessed,
        )
@flag_as_option
def CombineCards_smHcard(args):
    Commands.BasicCombineCards(
        'suppliedInput/combinedCard_smH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_smH_PTH,
        'hzz=' + LatestPaths.card_hzz_smH_PTH
        )

def pth_smH_get_card(args):
    if args.combination:
        card = LatestPaths.card_combined_smH_PTH
    elif args.hgg:
        card = LatestPaths.card_hgg_smH_PTH
    elif args.hzz:
        card = LatestPaths.card_hzz_smH_PTH
    elif args.hbb:
        raise NotImplementedError('Case --hbb in pth_smH_get_card is not yet implemented')
    elif args.combWithHbb:
        raise NotImplementedError('Case --combWithHbb in pth_smH_get_card is not yet implemented')
    else:
        raise DecayChannelNotFoundError()
    return card

def pth_smH_get_ws(args):
    if args.combination:
        ws = LatestPaths.ws_combined_smH
    elif args.hgg:
        ws = LatestPaths.ws_hgg_smH
    elif args.hzz:
        ws = LatestPaths.ws_hzz_smH
    elif args.hbb:
        raise NotImplementedError('Case --hbb in pth_smH_get_ws is not yet implemented')
    elif args.combWithHbb:
        raise NotImplementedError('Case --combWithHbb in pth_smH_get_ws is not yet implemented')
    else:
        raise DecayChannelNotFoundError()
    return ws

def pth_smH_get_postfit(args):
    decay_channel = get_decay_channel_tag(args)
    postfit_key = 'postfit_{0}_pth_smH'.format(decay_channel)
    try:
        return getattr(LatestPaths, postfit_key)
    except AttributeError:
        raise NotImplementedError('Postfit key \'{0}\' is not implemented in LatestPaths.py'.format(postfit_key))

@flag_as_option
def pth_smH_t2ws(args):
    card = pth_smH_get_card(args)
    if args.hzz:
        Commands.BasicT2WS(
            card,
            manualMaps = [
                '--PO \'map=.*/smH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_30_45:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_45_85:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_85_125:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_125_200:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_200_350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/smH_PTH_GT350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
                ],
            )
    elif args.hbb:
        raise NotImplementedError('Case --hbb in t2ws_smH is not yet implemented')
    elif args.combWithHbb:
        raise NotImplementedError('Case --combWithHbb in t2ws_smH is not yet implemented')
    else:
        Commands.BasicT2WS(
            LatestPaths.card_combined_smH_PTH,
            smartMaps = [
                ( r'.*/smH_PTH_([\d\_GT]+)', r'r_smH_PTH_\1[1.0,-1.0,4.0]' )
                ],
            )

@flag_as_option
def pth_smH_scan(args):
    ws = pth_smH_get_ws(args)
    jobDirectory = 'out/Scan_pth_smH_{0}_{1}'.format(datestr, get_decay_channel_tag(args))

    nPoints = 45
    nPointsPerJob = 3
    extraOptions = None
    
    if args.statonly:
        nPointsPerJob = 6
        jobDirectory += '_statonly'
        ws = pth_smH_get_postfit(args)
        # variables_to_freeze = get_all_vars_except_POIs(args, ws)
        # extraOptions = [
        #     '--snapshotName MultiDimFit',
        #     '--skipInitialFit',
        #     '--freezeNuisances {0}'.format(','.join(variables_to_freeze))
        #     ]
        extraOptions = [
            '--snapshotName MultiDimFit',
            '--skipInitialFit',
            '--freezeNuisancances=rgx{r_.*}'
            ]
    if args.hzz or args.hbb:
        nPointsPerJob = nPoints

    Commands.BasicCombineTool(
        ws,
        POIpattern    = '*',
        nPoints       = nPoints,
        nPointsPerJob = nPointsPerJob,
        jobDirectory  = jobDirectory,
        queue         = 'short.q',
        asimov        = args.asimov,
        extraOptions  = extraOptions
        )

@flag_as_option
def pth_smH_postfit(args):
    make_postfit(args, pth_smH_get_ws(args), 'pth_smH')

#____________________________________________________________________
# pth_ggH

# Preprocessing
@flag_as_option
def RenameHggProcesses_Aug21(args):
    MergeHGGWDatacards.RenameProcesses_Aug21(
        LatestPaths.card_hgg_ggHxH_PTH_unprocessed,
        )
@flag_as_option
def RenumberHzzProcesses_Aug21(args):
    MergeHGGWDatacards.RenumberProcessesHZZ_Aug21(
        LatestPaths.card_hzz_ggHxH_PTH_unprocessed,
        )
@flag_as_option
def CombineCards_Aug21(args):
    Commands.BasicCombineCards(
        'suppliedInput/combinedCard_ggHxH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_ggHxH_PTH,
        'hzz=' + LatestPaths.card_hzz_ggHxH_PTH
        )
@flag_as_option
def CombineCards_Dec15_hbb(args):
    Commands.BasicCombineCards(
        'suppliedInput/combinedCard_hgg_hzz_hbb_ggHxH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_ggHxH_PTH,
        'hzz=' + LatestPaths.card_hzz_ggHxH_PTH,
        'hbb=' + LatestPaths.card_hbb_ggHxH_PTH
        )

def pth_ggH_get_card(args):
    if args.combination:
        datacard = LatestPaths.card_combined_ggHxH_PTH
    elif args.hgg:
        datacard = LatestPaths.card_hgg_ggHxH_PTH
    elif args.hzz:
        datacard = LatestPaths.card_hzz_ggHxH_PTH
    elif args.hbb:
        datacard = LatestPaths.card_hbb_ggHxH_PTH
    elif args.combWithHbb:
        datacard = LatestPaths.card_combinedWithHbb_ggHxH_PTH
    else:
        raise DecayChannelNotFoundError()
    return datacard

def pth_ggH_get_ws(args):
    if args.combination:
        ws = LatestPaths.ws_combined_ggH_xHfixed
    elif args.hgg:
        ws = LatestPaths.ws_hgg_ggH_xHfixed
    elif args.hzz:
        ws = LatestPaths.ws_hzz_ggH_xHfixed
    elif args.hbb:
        ws = LatestPaths.ws_hbb_ggH_xHfixed
    elif args.combWithHbb:
        ws = LatestPaths.ws_combWithHbb_ggH_xHfixed
    else:
        raise DecayChannelNotFoundError()
    return ws

def pth_ggH_get_postfit(args):
    decay_channel = get_decay_channel_tag(args)
    postfit_key = 'postfit_{0}_pth_ggH'.format(decay_channel)
    try:
        return getattr(LatestPaths, postfit_key)
    except AttributeError:
        raise NotImplementedError('Postfit key \'{0}\' is not implemented in LatestPaths.py'.format(postfit_key))

@flag_as_option
def pth_ggH_t2ws(args):
    datacard = pth_ggH_get_card(args)
    outputWS = basename(datacard).replace( '.txt', '_xHfixed.root' )
    if args.hzz:
        Commands.BasicT2WS(
            datacard,
            manualMaps=[
                '--PO \'map=.*/ggH_PTH_0_15:r_ggH_PTH_0_15[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_15_30:r_ggH_PTH_15_30[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_30_45:r_ggH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_45_85:r_ggH_PTH_30_85[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_85_125:r_ggH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_125_200:r_ggH_PTH_85_200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_200_350:r_ggH_PTH_GT200[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_GT350:r_ggH_PTH_GT200[1.0,0.0,3.0]\'',
                ],
            outputWS=outputWS
            )
    elif args.hbb:
        Commands.BasicT2WS(
            datacard,
            manualMaps=[
                # '--PO \'map=.*/ggH_PTH_200_350:r_ggH_PTH_200_350[1.0,0.0,3.0]\'',
                '--PO \'map=.*/ggH_PTH_350_600:r_ggH_PTH_350_600[1.0,0.0,10.0]\'',
                '--PO \'map=.*/ggH_PTH_GT600:r_ggH_PTH_GT600[1.0,0.0,10.0]\'',
                ],
            outputWS=outputWS
            )
    elif args.combWithHbb:
        Commands.BasicT2WS(
            datacard,
            manualMaps=[
                '--PO \'map=.*/ggH_PTH_GT350:r_ggH_PTH_350_600[1.0,0.0,10.0]\'',
                '--PO \'map=.*/ggH_PTH_350_600:r_ggH_PTH_350_600[1.0,0.0,10.0]\'',
                '--PO \'map=.*/ggH_PTH_GT600:r_ggH_PTH_GT600[1.0,0.0,10.0]\'',
                ],
            smartMaps = [
                ( r'.*/ggH_PTH_([\d\_GT]+)', r'r_ggH_PTH_\1[1.0,-1.0,4.0]' )
                ],
            outputWS=outputWS
            )
    else:
        Commands.BasicT2WS(
            datacard,
            smartMaps = [
                ( r'.*/ggH_PTH_([\d\_GT]+)', r'r_ggH_PTH_\1[1.0,-1.0,4.0]' )
                ],
            outputWS=outputWS
            )

@flag_as_option
def pth_ggH_scan(args):
    ws = pth_ggH_get_ws(args)
    jobDirectory  = 'out/Scan_pth_ggH_{0}_{1}'.format(datestr, get_decay_channel_tag(args))

    nPoints = 45
    nPointsPerJob = 3

    extraOptions = []
    physicsModelParameterRanges = []
    POIRange = None

    if args.statonly:
        nPointsPerJob = 6
        jobDirectory += '_statonly'
        ws = pth_ggH_get_postfit(args)
        variables_to_freeze = get_all_vars_except_POIs(args, ws)
        extraOptions.extend([
            '--snapshotName MultiDimFit',
            '--skipInitialFit',
            '--freezeNuisances {0}'.format(','.join(variables_to_freeze))
            ])

    if args.hzz:
        nPointsPerJob = nPoints
        POIRange = [ 0.0, 4.0 ]
    elif args.hbb:
        nPointsPerJob = nPoints
        POIRange = [ -10.0, 10.0 ]
        extraOptions.extend([
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            ])
    elif args.combWithHbb:
        POIRange = [ 0.0, 8.5 ]
        nPoints = 150
        nPointsPerJob = 2
        extraOptions.extend([
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            ])
        if not args.statonly:
            physicsModelParameterRanges = [
                [ 'qcdeff', 0.001, 8.0 ],
                [ 'r1p0', 0.0, 8.0 ],
                [ 'r2p0', 0.0, 8.0 ],
                [ 'r3p0', 0.0, 8.0 ],
                [ 'r0p1', 0.0, 8.0 ],
                [ 'r1p1', 0.0, 8.0 ],
                [ 'r2p1', 0.0, 8.0 ],
                [ 'r3p1', 0.0, 8.0 ],
                ]

    Commands.BasicCombineTool(
        ws,
        POIpattern    = '*',
        POIRange      = POIRange,
        nPoints       = nPoints,
        nPointsPerJob = nPointsPerJob,
        jobDirectory  = jobDirectory,
        queue         = 'short.q',
        asimov        = args.asimov,
        extraOptions  = extraOptions,
        physicsModelParameterRanges = physicsModelParameterRanges
        )

@flag_as_option
def pth_ggH_postfit(args):
    make_postfit(args, pth_ggH_get_ws(args), 'pth_ggH')

#____________________________________________________________________
# njets

@flag_as_option
def rename_hgg_njets(args):
    MergeHGGWDatacards.RenameProcesses_Hgg_nJets(
        LatestPaths.card_hgg_smH_NJ_unprocessed,
        globalReplace = [
            ( 'CMS_hgg_JEC', 'CMS_scale_j' )
            ]
        )

@flag_as_option
def njets_combineCards(args):
    Commands.BasicCombineCards(
        'suppliedInput/combinedCard_NJ_smH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_smH_NJ,
        'hzz=' + LatestPaths.card_hzz_smH_NJ
        )

def njets_get_ws(args):
    if args.combination:
        ws = LatestPaths.ws_combined_smH_NJ
    elif args.hgg:
        ws = LatestPaths.ws_hgg_smH_NJ
    elif args.hzz:
        ws = LatestPaths.ws_hzz_smH_NJ
    elif args.hbb or args.combWithHbb:
        raise NotImplementedError()
    else:
        raise DecayChannelNotFoundError()
    return ws    

def njets_get_card(args):
    if args.combination:
        datacard = LatestPaths.card_combined_smH_NJ
    if args.hzz:
        datacard = LatestPaths.card_hzz_smH_NJ
    elif args.hgg:
        datacard = LatestPaths.card_hgg_smH_NJ
    elif args.hbb or args.combWithHbb:
        raise NotImplementedError()
    else:
        raise DecayChannelNotFoundError()
    return datacard

@flag_as_option
def njets_t2ws(args):
    datacard = njets_get_card(args)
    if args.hzz:
        Commands.BasicT2WS(
            datacard,
            manualMaps=[
                '--PO \'map=.*/smH_NJ_0:r_smH_NJ_0[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_NJ_1:r_smH_NJ_1[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_NJ_2:r_smH_NJ_2[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_NJ_3:r_smH_NJ_GE3[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_NJ_GE4:r_smH_NJ_GE3[1.0,-1.0,4.0]\'',
                ],
            )
    else:
        Commands.BasicT2WS(
            datacard,
            smartMaps = [
                ( r'.*/smH_NJ_([\d\_GE]+)', r'r_smH_NJ_\1[1.0,-1.0,4.0]' )
                ],
            )

@flag_as_option
def njets_scan(args):
    nPoints       = 39
    nPointsPerJob = 3
    ws = njets_get_ws(args)
    if args.hzz:
        nPoints       = 39
        nPointsPerJob = 39
    Commands.BasicCombineTool(
        ws,
        POIpattern    = '*',
        nPoints       = nPoints,
        nPointsPerJob = nPointsPerJob,
        jobDirectory  = 'Scan_njets_{0}'.format(datestr),
        queue         = 'short.q',
        asimov        = args.asimov,
        )

@flag_as_option
def njets_postfit(args):
    make_postfit(args, njets_get_ws(args), 'njets')

#____________________________________________________________________
# ptjet

@flag_as_option
def ptjet_rename_hgg(args):
    MergeHGGWDatacards.RenameProcesses_Hgg_differentials(
        LatestPaths.card_hgg_smH_PTJ_unprocessed
        )

@flag_as_option
def ptjet_combineCards(args):
    Commands.BasicCombineCards(
        'suppliedInput/combinedCard_PTJ_smH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_smH_PTJ,
        'hzz=' + LatestPaths.card_hzz_smH_PTJ
        )

def ptjet_get_card(args):
    if args.combination:
        datacard = LatestPaths.card_combined_smH_PTJ
    elif args.hgg:
        datacard = LatestPaths.card_hgg_smH_PTJ
    elif args.hzz:
        datacard = LatestPaths.card_hzz_smH_PTJ
    elif args.hbb or args.combWithHbb:
        raise NotImplementedError()
    else:
        raise DecayChannelNotFoundError()
    return card

def ptjet_get_ws(args):
    if args.combination:
        ws = LatestPaths.ws_combined_smH_PTJ
    elif args.hgg:
        ws = LatestPaths.ws_hgg_smH_PTJ
    elif args.hzz:
        ws = LatestPaths.ws_hzz_smH_PTJ
    elif args.hbb or args.combWithHbb:
        raise NotImplementedError()
    else:
        raise DecayChannelNotFoundError()
    return ws

@flag_as_option
def ptjet_t2ws(args):
    datacard = ptjet_get_card(args)
    if args.hzz:
        Commands.BasicT2WS(
            datacard,
            manualMaps=[
                '--PO \'map=.*/smH_PTJ_LT30:r_smH_PTJ_LT30[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTJ_30_55:r_smH_PTJ_30_55[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTJ_55_95:r_smH_PTJ_55_95[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTJ_95_120:r_smH_PTJ_GT95[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTJ_120_200:r_smH_PTJ_GT95[1.0,-1.0,4.0]\'',
                '--PO \'map=.*/smH_PTJ_GT200:r_smH_PTJ_GT95[1.0,-1.0,4.0]\'',
                ],
            suffix = '_ptjet'
            )
    else:
        Commands.BasicT2WS(
            datacard,
            smartMaps = [
                ( r'.*/smH_PTJ_([pm\d\_GELT]+)', r'r_smH_PTJ_\1[1.0,-1.0,4.0]' )
                ],
            # suffix = '_ptjet'
            )

@flag_as_option
def ptjet_scan(args):
    nPoints       = 55
    nPointsPerJob = 5
    if args.hzz:
        nPointsPerJob = nPoints
    Commands.BasicCombineTool(
        ptjet_get_ws(args),
        POIpattern    = '*',
        nPoints       = nPoints,
        nPointsPerJob = nPointsPerJob,
        jobDirectory  = 'Scan_PTJ_{0}'.format(datestr),
        queue         = 'short.q',
        asimov        = args.asimov,
        )

@flag_as_option
def ptjet_postfit(args):
    make_postfit(args, ptjet_get_ws(args), 'ptjet')

#____________________________________________________________________
# rapidity

@flag_as_option
def rename_hgg_rapidity(args):
    MergeHGGWDatacards.RenameProcesses_Hgg_differentials(
        LatestPaths.card_hgg_smH_YH_unprocessed,
        )

@flag_as_option
def rapidity_combineCards(args):
    Commands.BasicCombineCards(
        'suppliedInput/combinedCard_YH_smH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_smH_YH,
        'hzz=' + LatestPaths.card_hzz_smH_YH
        )

def rapidity_get_card(args):
    if args.combination:
        datacard = LatestPaths.card_combined_smH_YH
    elif args.hgg:
        datacard = LatestPaths.card_hgg_smH_YH
    elif args.hzz:
        datacard = LatestPaths.card_hzz_smH_YH
    elif args.hbb or args.combWithHbb:
        raise NotImplementedError()
    else:
        raise DecayChannelNotFoundError()
    return datacard

def rapidity_get_ws(args):
    if args.combination:
        ws = LatestPaths.ws_combined_smH_YH
    elif args.hgg:
        ws = LatestPaths.ws_hgg_smH_YH
    elif args.hzz:
        ws = LatestPaths.ws_hzz_smH_YH
    elif args.hbb or args.combWithHbb:
        raise NotImplementedError()
    else:
        raise DecayChannelNotFoundError()
    return ws    

@flag_as_option
def rapidity_t2ws(args):
    Commands.BasicT2WS(
        rapidity_get_card(args),
        smartMaps = [
            ( r'.*/smH_YH_([pm\d\_GE]+)', r'r_smH_YH_\1[1.0,-1.0,4.0]' )
            ],
        )

@flag_as_option
def rapidity_scan(args):
    nPoints       = 55
    nPointsPerJob = 5
    if args.hzz:
        nPointsPerJob = nPoints
    Commands.BasicCombineTool(
        rapidity_get_ws(args),
        POIpattern    = '*',
        nPoints       = nPoints,
        nPointsPerJob = nPointsPerJob,
        jobDirectory  = 'Scan_YH_{0}'.format(datestr),
        queue         = 'short.q',
        asimov        = args.asimov,
        )

@flag_as_option
def rapidity_postfit(args):
    make_postfit(args, rapidity_get_ws(args), 'rapidity')


########################################
# Plotting
########################################

@flag_as_option
def plot_all_differentials(args):
    pth_plot(args)
    pth_ggH_plot(args)
    pth_ggH_hbb_plot(args)
    njets_plot(args)
    ptjet_plot(args)
    rapidity_plot(args)

#____________________________________________________________________
@flag_as_option
def pth_plot(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
    containers = []

    # hgg = prepare_container('hgg', LatestPaths.ws_hgg_smH, LatestPaths.scan_hgg_PTH)
    # hgg.SMcrosssections = LatestBinning.obs_pth.crosssection_over_binwidth()
    # containers.append(hgg)

    # hzz = prepare_container('hzz', LatestPaths.ws_hzz_smH, LatestPaths.scan_hzz_PTH)
    # hzz.SMcrosssections = LatestBinning.obs_pth_hzzBinning.crosssection_over_binwidth()
    # containers.append(hzz)

    combination = prepare_container('combination', LatestPaths.ws_combined_smH, LatestPaths.scan_combined_PTH)
    combination.SMcrosssections = LatestBinning.obs_pth.crosssection_over_binwidth()
    containers.append(combination)

    combination_statonly = prepare_container('combination_statonly', LatestPaths.postfit_combination_pth_smH, 'out/Scan_pth_smH_Feb05_combination_statonly')
    combination_statonly.SMcrosssections = LatestBinning.obs_pth.crosssection_over_binwidth()
    combination_statonly.color = 2
    containers.append(combination_statonly)

    for container in containers:
        PlotCommands.WriteScansToTable(
            container,
            'pth',
            xTitle = 'p_{T}^{H} (GeV)',
            yTitle = '#Delta#sigma/#Delta p_{T}^{H} (pb/GeV)',
            lastBinIsOverflow = True,
            )

    SM = prepare_SM_container(
        LatestBinning.obs_pth.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True),
        LatestBinning.obs_pth.binning
        )
    containers.append(SM)

    PlotCommands.PlotSpectraOnTwoPanel(
        'twoPanel_pthSpectrum',
        containers,
        xTitle = 'p_{T}^{H} (GeV)',
        yTitleTop = '#Delta#sigma/#Deltap_{T}^{H} (pb/GeV)',
        # 
        # yMinExternalTop = 0.0005,
        # yMaxExternalTop = 110.,
        )

#____________________________________________________________________
@flag_as_option
def pth_ggH_plot(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
    containers = []

    hgg = prepare_container( 'hgg', LatestPaths.ws_hgg_ggH_xHfixed, LatestPaths.scan_hgg_PTH_ggH )
    hgg.SMcrosssections = LatestBinning.obs_pth_ggH.crosssection_over_binwidth()
    containers.append(hgg)

    hzz         = prepare_container( 'hzz', LatestPaths.ws_hzz_ggH_xHfixed, LatestPaths.scan_hzz_PTH_ggH )
    hzz.SMcrosssections = LatestBinning.obs_pth_ggH_hzzBinning.crosssection_over_binwidth()
    containers.append(hzz)

    combination = prepare_container( 'combination', LatestPaths.ws_combined_ggH_xHfixed, LatestPaths.scan_combined_PTH_ggH )
    combination.SMcrosssections = LatestBinning.obs_pth_ggH.crosssection_over_binwidth()
    containers.append(combination)

    for container in containers:
        PlotCommands.WriteScansToTable(
            container,
            'pth_ggh',
            xTitle = 'p_{T}^{H} (GeV)',
            yTitle = '#Delta#sigma/#Delta p_{T}^{H} (pb/GeV)',
            lastBinIsOverflow = True,
            )

    SM = prepare_SM_container(
        LatestBinning.obs_pth_ggH.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True),
        LatestBinning.obs_pth_ggH.binning
        )
    containers.append(SM)

    l = PlotCommands.TLatexMultiPanel(
        lambda c: 1.0 - c.GetRightMargin() - 0.01,
        lambda c: 1.0 - c.GetTopMargin() - 0.14,
        '(non-ggH fixed to SM)'
        )
    l.SetTextSize(0.05)
    l.SetTextAlign(33)

    PlotCommands.PlotSpectraOnTwoPanel(
        'twoPanel_pth_ggH_Spectrum',
        containers,
        xTitle = 'p_{T}^{H} (GeV)',
        yTitleTop = '#Delta#sigma^{ggH}/#Deltap_{T}^{H} (pb/GeV)',
        topPanelObjects = [ ( l, '' ) ],
        )

#____________________________________________________________________
@flag_as_option
def pth_ggH_hbb_plot(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
    containers = []

    if args.asimov:
        combined_scanDir    = LatestPaths.scan_combined_PTH_ggH_asimov
        hgg_scanDir         = LatestPaths.scan_hgg_PTH_ggH_asimov
        hzz_scanDir         = LatestPaths.scan_hzz_PTH_ggH_asimov
        hbb_scanDir         = LatestPaths.scan_hbb_PTH_ggH_asimov
        combWithHbb_scanDir = LatestPaths.scan_combWithHbb_PTH_ggH_asimov
    else:
        combined_scanDir    = LatestPaths.scan_combined_PTH_ggH
        hgg_scanDir         = LatestPaths.scan_hgg_PTH_ggH
        hzz_scanDir         = LatestPaths.scan_hzz_PTH_ggH
        hbb_scanDir         = LatestPaths.scan_hbb_PTH_ggH
        combWithHbb_scanDir = LatestPaths.scan_combWithHbb_PTH_ggH


    hgg = prepare_container( 'hgg', LatestPaths.ws_hgg_ggH_xHfixed, hgg_scanDir )
    hgg.SMcrosssections = LatestBinning.obs_pth_ggH.crosssection_over_binwidth()
    containers.append(hgg)

    hzz = prepare_container( 'hzz', LatestPaths.ws_hzz_ggH_xHfixed, hzz_scanDir )
    hzz.SMcrosssections = LatestBinning.obs_pth_ggH_hzzBinning.crosssection_over_binwidth()
    containers.append(hzz)

    combination = prepare_container( 'combination', LatestPaths.ws_combined_ggH_xHfixed, combined_scanDir )
    combination.SMcrosssections = LatestBinning.obs_pth_ggH.crosssection_over_binwidth()
    containers.append(combination)

    hbb = prepare_container( 'hbb', LatestPaths.ws_hbb_ggH_xHfixed, hbb_scanDir )
    hbb.SMcrosssections = LatestBinning.obs_pth_ggH_hbbBinning.crosssection_over_binwidth()
    containers.append(hbb)

    combWithHbb = prepare_container( 'combWithHbb', LatestPaths.ws_combWithHbb_ggH_xHfixed, combWithHbb_scanDir )
    combWithHbb.SMcrosssections = LatestBinning.obs_pth_ggH_combWithHbbBinning.crosssection_over_binwidth()
    containers.append(combWithHbb)

    check_containers(containers)
    for container in containers:
        PlotCommands.WriteScansToTable(
            container,
            'pth_ggh_whbb',
            xTitle = 'p_{T}^{H} (GeV)',
            yTitle = '#Delta#sigma/#Delta p_{T}^{H} (pb/GeV)',
            lastBinIsOverflow = True,
            )

    SM = prepare_SM_container(
        LatestBinning.obs_pth_ggH_combWithHbbBinning.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True),
        LatestBinning.obs_pth_ggH_combWithHbbBinning.binning
        )
    containers.append(SM)

    l = PlotCommands.TLatexMultiPanel(
        lambda c: 1.0 - c.GetRightMargin() - 0.01,
        lambda c: 1.0 - c.GetTopMargin() - 0.14,
        '(non-ggH fixed to SM)'
        )
    l.SetTextSize(0.05)
    l.SetTextAlign(33)

    PlotCommands.PlotSpectraOnTwoPanel(
        'twoPanel_pth_ggH_hbb_Spectrum' + ( '_asimov' if args.asimov else '' ),
        containers,
        xTitle = 'p_{T}^{H} (GeV)',
        yTitleTop = '#Delta#sigma^{ggH}/#Deltap_{T}^{H} (pb/GeV)',
        # topPanelObjects = [ ( l, '' ) ],
        )

#____________________________________________________________________
@flag_as_option
def njets_plot(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
    containers = []

    hgg = prepare_container( 'hgg', LatestPaths.ws_hgg_smH_NJ,LatestPaths.scan_hgg_NJ )
    hgg.SMcrosssections = LatestBinning.obs_njets.crosssection_over_binwidth()
    containers.append(hgg)

    hzz = prepare_container( 'hzz', LatestPaths.ws_hzz_smH_NJ,LatestPaths.scan_hzz_NJ )
    hzz.SMcrosssections = LatestBinning.obs_njets_hzzBinning.crosssection_over_binwidth()
    containers.append(hzz)

    combination = prepare_container( 'combination', LatestPaths.ws_combined_smH_NJ,LatestPaths.scan_combined_NJ )
    combination.SMcrosssections = LatestBinning.obs_njets.crosssection_over_binwidth()
    containers.append(combination)

    check_containers(containers)
    for container in containers:
        PlotCommands.WriteScansToTable(
            container,
            'njets',
            xTitle = 'N_{jets}',
            yTitle = '#Delta#sigma/#Delta N_{jets} (pb)',
            lastBinIsOverflow = False,
            )

    SM = prepare_SM_container(
        LatestBinning.obs_njets.crosssection_over_binwidth(),
        LatestBinning.obs_njets.binning
        )
    containers.append(SM)

    PlotCommands.PlotSpectraOnTwoPanel(
        'twoPanel_nJetsSpectrum',
        containers,
        xTitle = 'N_{jets}',
        # yMinLimit = 0.07,
        # yMaxExternalTop = 500,
        yTitleTop = '#Delta#sigma/#DeltaN_{jets} (pb)',
        )

#____________________________________________________________________
@flag_as_option
def ptjet_plot(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
    containers = []

    hgg = prepare_container('hgg', LatestPaths.ws_hgg_smH_PTJ, LatestPaths.scan_hgg_PTJ_asimov)
    hgg.SMcrosssections = LatestBinning.obs_ptjet.crosssection_over_binwidth()
    containers.append(hgg)

    hzz = prepare_container('hzz', LatestPaths.ws_hzz_smH_PTJ, LatestPaths.scan_hzz_PTJ_asimov)
    hzz.SMcrosssections = LatestBinning.obs_ptjet_hzzBinning.crosssection_over_binwidth()
    containers.append(hzz)

    combination = prepare_container('combination', LatestPaths.ws_combined_smH_PTJ, LatestPaths.scan_combined_PTJ_asimov)
    combination.SMcrosssections = LatestBinning.obs_ptjet.crosssection_over_binwidth()
    containers.append(combination)

    Commands.Warning( 'Skipping first bin for ptjet (should be the underflow)' )
    for container in containers:
        container.POIs = container.POIs[1:]
        container.Scans = container.Scans[1:]
        container.SMcrosssections = container.SMcrosssections[1:]

    check_containers(containers)
    for container in containers:
        PlotCommands.WriteScansToTable(
            container,
            'ptjet',
            xTitle = 'p_{T}^{jet} (GeV)',
            yTitle = '#Delta#sigma/#Delta p_{T}^{jet} (pb/GeV)',
            lastBinIsOverflow = True,
            )

    SM = prepare_SM_container(
        LatestBinning.obs_ptjet.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True),
        LatestBinning.obs_ptjet.binning
        )
    containers.append(SM)

    PlotCommands.PlotSpectraOnTwoPanel(
        'twoPanel_ptjetSpectrum',
        containers,
        xTitle = 'p_{T}^{jet} (GeV)',
        yTitleTop = '#Delta#sigma/#Deltap_{T}^{jet} (pb/GeV)',
        # yMinLimit = 0.07,
        yMaxExternalTop = 10,
        xMinExternal = 30.0,
        # yMinLimit    = 0.1
        )

#____________________________________________________________________
@flag_as_option
def rapidity_plot(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
    containers = []

    hgg = prepare_container('hgg', LatestPaths.ws_hgg_smH_YH, LatestPaths.scan_hgg_YH_asimov)
    hgg.SMcrosssections = LatestBinning.obs_yh.crosssection_over_binwidth()
    containers.append(hgg)

    hzz = prepare_container('hzz', LatestPaths.ws_hzz_smH_YH, LatestPaths.scan_hzz_YH_asimov)
    hzz.SMcrosssections = LatestBinning.obs_yh.crosssection_over_binwidth()
    containers.append(hzz)

    combination = prepare_container('combination', LatestPaths.ws_combined_smH_YH, LatestPaths.scan_combined_YH_asimov)
    combination.SMcrosssections = LatestBinning.obs_yh.crosssection_over_binwidth()
    containers.append(combination)

    check_containers(containers)
    for container in containers:
        PlotCommands.WriteScansToTable(
            container,
            'rapidity',
            xTitle = '|y_{H}|',
            yTitle = '#Delta#sigma/#Delta|y_{H}| (pb)',
            lastBinIsOverflow = False,
            )

    SM = prepare_SM_container(
        LatestBinning.obs_yh.crosssection_over_binwidth(),
        LatestBinning.obs_yh.binning
        )
    containers.append(SM)

    PlotCommands.PlotSpectraOnTwoPanel(
        'twoPanel_rapiditySpectrum',
        containers,
        xTitle = '|y_{H}|',
        yTitleTop = '#Delta#sigma/#Delta|y_{H}| (pb)',
        # yMinLimit = 0.07,
        # yMaxExternalTop = 500
        lastBinIsNotOverflow=True,
        )

#____________________________________________________________________
# Helpers

def get_decay_channel_tag(args):
    if args.combination:
        tag = 'combination'
    elif args.hgg:
        tag = 'hgg'
    elif args.hzz:
        tag = 'hzz'
    elif args.hbb:
        tag = 'hbb'
    elif args.combWithHbb:
        tag = 'combWithHbb'
    else:
        raise DecayChannelNotFoundError()
    return tag

def get_all_vars_except_POIs(args, postfit):
    with Commands.OpenRootFile(postfit) as postfit_fp:
        w = postfit_fp.Get('w')

    variables = []
    var_set = w.allVars()
    var_iterable = var_set.createIterator()
    for i_var in xrange(var_set.getSize()):
        variables.append(var_iterable.Next().GetName())
    
    POIs = Commands.ListSet(w)

    variables = list(set(variables) - set(POIs))
    variables.sort()
    return variables

def make_postfit(args, ws, observable_tag=None):
    tag = 'POSTFIT_' + basename(ws).replace('.root','')
    if not(observable_tag is None):
        tag += '_' + observable_tag
    if args.asimov:
        tag += '_asimov'
    cmd = [
        'combineTool.py',
        abspath(ws),
        '-n _{0}'.format(tag),
        '-M MultiDimFit',
        '--cminDefaultMinimizerType Minuit2',
        '--cminDefaultMinimizerAlgo migrad',
        '--floatOtherPOIs=1',
        # '-P "{0}"'.format( POI ),
        # '--setPhysicsModelParameterRanges "{0}"={1:.3f},{2:.3f} '.format( POI, POIRange[0], POIRange[1] ),
        # '--setPhysicsModelParameters {0}'.format( ','.join([ iterPOI + '=1.0' for iterPOI in allPOIs ]) ),
        '-m 125.00',
        '--saveNLL',
        '--saveInactivePOI 1',
        '--saveWorkspace',
        '--job-mode psi --task-name {0} --sub-opts=\'-q short.q\' '.format(tag),
        ]
    if args.asimov:
        cmd.append( '-t -1' )

    if not 't3' in os.environ['HOSTNAME']:
        raise NotImplementedError('\'t3\' not in $HOSTNAME; Is this executed from the T3?')

    with Commands.EnterDirectory('out/postfitWss_{0}'.format(datestr)):
        Commands.executeCommand(cmd)


def check_containers(containers):
    for container in containers:
        if not len(container.POIs) == len(container.SMcrosssections):
            Commands.ThrowError(
                'For container {2}, found {0} POIs, but {1} SM cross sections; something is misaligned.'.format(
                    len(container.POIs), len(container.SMcrosssections), container.name )
                + '\n  POIs:  {0}'.format(container.POIs)
                + '\n  SM xs: {0}'.format(container.SMcrosssections)
                )

def prepare_container(
        name,
        ws, scandir,
        title=None,
        verbose=False,
        draw_parabolas=True,
        ):

    POIs = Commands.ListPOIs( ws )
    POIs.sort( key=Commands.POIsorter )
    if verbose: print 'Sorted POIs:', POIs

    scans = PhysicsCommands.GetScanResults(
        POIs,
        scandir,
        # pattern = pattern,
        filterNegatives = True
        )

    if draw_parabolas:
        PhysicsCommands.BasicDrawScanResults( POIs, scans, name=name )

    container = Container(name=name)
    if title is None:
        title = name
    container.title = title

    container.POIs = POIs
    container.Scans = scans

    if name == 'hgg':
        container.color = 2
        container.title = 'H#rightarrow#gamma#gamma'
    if name == 'hzz':
        container.color = 4
        container.title = 'H#rightarrowZZ'
    if name == 'combination':
        container.color = 1
        container.title = 'Combination'
    if name == 'hbb':
        container.color = 8
        container.title = 'H#rightarrowbb'
    if name == 'combWithHbb':
        container.color = 14
        container.title = 'Comb. with H#rightarrowbb'

    return container

def prepare_SM_container(crosssection, binBoundaries):
    SM = Container()
    SM.name = 'SM'
    SM.title = 'SM'
    SM.color = 14
    SM.crosssection = crosssection
    SM.binBoundaries = binBoundaries
    SM.ratios = [ 1.0 for i in xrange(len(crosssection)) ]
    return SM

class DecayChannelNotFoundError(Exception):
    def __init__(self):
        Exception.__init__(self, 'Could not determine the decay channel') 
