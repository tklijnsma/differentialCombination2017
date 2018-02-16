#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, sys, re
from os.path import *
from glob import glob
from copy import deepcopy

from OptionHandler import flag_as_option

import differentials

# sys.path.append('src')
# import Commands
# import PhysicsCommands
# import TheoryCommands
# import LatestPaths
# import LatestPathsGetters
# import LatestBinning
# from Container import Container
# import PlotCommands
# from differentialTools import *
# import CombineToolWrapper

from time import strftime
datestr = strftime( '%b%d' )


########################################
# Main
########################################
        
@flag_as_option
def pth_smH_scan(args):
    ws = LatestPathsGetters.get_ws('pth_smH', args)
    config = differential_config(args, ws, 'pth_smH')
    scan_directly(args, config)

@flag_as_option
def pth_ggH_scan(args):
    ws = LatestPathsGetters.get_ws('pth_ggH', args)
    config = differential_config(args, ws, 'pth_ggH')
    scan_directly(args, config)

@flag_as_option
def njets_scan(args):
    ws = LatestPathsGetters.get_ws('njets', args)
    config = differential_config(args, ws, 'njets')
    scan_directly(args, config)

@flag_as_option
def ptjet_scan(args):
    ws = LatestPathsGetters.get_ws('ptjet', args)
    config = differential_config(args, ws, 'ptjet')
    scan_directly(args, config)

@flag_as_option
def rapidity_scan(args):
    ws = LatestPathsGetters.get_ws('rapidity', args)
    config = differential_config(args, ws, 'rapidity')
    scan_directly(args, config)

@flag_as_option
def scan_all(real_args):
    args = deepcopy(real_args)
    decay_channels = [
            'hgg',
            'hzz',
            # Workspaces are not working yet; Wait for Vs input
            # 'hbb',
            # 'combWithHbb',
            'combination'
            ]
    try:
        decay_channel = get_decay_channel_tag(args)
        Commands.warning('Decay channel {0} was specified, so only scans relating to it are submitted'.format(decay_channel))
        decay_channels = [decay_channel]
    except DecayChannelNotFoundError:
        pass

    scan_functions = [
        pth_smH_scan,
        pth_ggH_scan,
        ptjet_scan,
        rapidity_scan,
        njets_scan,
        ]

    for decay_channel in decay_channels:
        set_all_false(args)
        setattr(args, decay_channel, True)
        print '\nSubmitting scans for {0}'.format(decay_channel)
        for function in scan_functions:
            try:
                function(args)
            except NotImplementedError:
                pass


#____________________________________________________________________
# Helpers

def set_all_false(args):
    args.hgg = False
    args.hzz = False
    args.hbb = False
    args.combWithHbb = False

lumiMultiplier300 = 8.356546
lumiMultiplier3000 = 83.56546
lumiMultiplier = lumiMultiplier300

def differential_config(args, ws, obs_name):
    base_config = CombineToolWrapper.CombineConfig(args)
    base_config.onBatch       = True
    base_config.nPoints       = 55
    base_config.nPointsPerJob = 5
    base_config.queue         = 'short.q'

    if args.asimov:
        base_config.asimov = True
    else:
        base_config.asimov = False

    base_config.decay_channel = get_decay_channel_tag(args)    
    if args.hzz or args.hbb:
        base_config.nPointsPerJob = base_config.nPoints
    if args.combWithHbb:
        base_config.nPoints = 150
        base_config.nPointsPerJob = 2
    if args.hbb or args.combWithHbb:
        base_config.minimizer_settings = [
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            ]

    base_config.datacard = ws
    base_config.subDirectory = 'out/Scan_{0}_{1}_{2}'.format(obs_name, datestr, base_config.decay_channel)
    if args.lumiScale:
        base_config.hardPhysicsModelParameters.append('lumiScale={0}'.format(lumiMultiplier))
        base_config.freezeNuisances.append('lumiScale')
        base_config.subDirectory += '_lumiScale{0}'.format(str(lumiMultiplier).replace('.','p'))
    if args.statonly:
        base_config.subDirectory += '_statonly'
    if args.asimov:
        base_config.subDirectory += '_asimov'

    POIs = Commands.list_pois(ws)
    POIrange = [ -1.0, 4.0 ]
    if args.hzz:
        POIrange[0] = 0.0
    if args.hbb:
        POIrange = [ -10.0, 10.0 ]
    if args.combWithHbb:
        POIrange = [ 0.0, 8.5 ]
    if args.lumiScale:
        POIrange = [ 0.7, 1.3 ]
        base_config.nPoints = 35
        if lumiMultiplier > 10.:
            POIrange = [ 0.9, 1.1 ]

    base_config.POIs = POIs
    base_config.PhysicsModelParameters = [ '{0}=1.0'.format(p) for p in POIs ]
    base_config.PhysicsModelParameterRanges = [ '{0}={1},{2}'.format(p, *POIrange) for p in POIs ]
    if args.hbb or args.combWithHbb:
        Commands.warning(
            'For hbb: qcdeff and r*p* parameters are given ranges! '
            'Should decide on (freezing) behaviour during postfit/scan'
            )
        base_config.PhysicsModelParameterRanges.extend([
                [ 'qcdeff', 0.001, 8.0 ],
                [ 'r1p0', 0.0, 8.0 ],
                [ 'r2p0', 0.0, 8.0 ],
                [ 'r3p0', 0.0, 8.0 ],
                [ 'r0p1', 0.0, 8.0 ],
                [ 'r1p1', 0.0, 8.0 ],
                [ 'r2p1', 0.0, 8.0 ],
                [ 'r3p1', 0.0, 8.0 ],
                ])
    return base_config

def differential_configs_for_scanning(args, original_base_config):
    base_config = deepcopy(original_base_config)

    if args.statonly:
        base_config.freezeNuisances.append('rgx{.*}')

    # Overwrite get_name function to include the POI name
    def new_get_name(self):
        return 'bPOI_{0}_ePOI_{1}'.format(self.POIs[0], basename(self.datacard).replace('.root',''))
    base_config.get_name = new_get_name.__get__(base_config, type(base_config))

    configs_per_POI = []
    for POI in base_config.POIs:
        config = deepcopy(base_config)
        config.POIs = [POI]
        configs_per_POI.append(config)
    return configs_per_POI

def postfit_and_scan(args, postfit_config):
    # Make sure no previous run directory is overwritten
    postfit_config.make_unique_directory()

    postfit = CombineToolWrapper.CombinePostfit(postfit_config)
    postfit.run()
    postfit_ws = postfit.get_output()

    scan_configs = differential_configs_for_scanning(args, postfit_config)
    for config in scan_configs:
        scan = CombineToolWrapper.CombineScanFromPostFit(config)
        scan.run(postfit_ws)

def scan_directly(args, base_config):
    base_config.make_unique_directory()
    scan_configs = differential_configs_for_scanning(args, base_config)
    for config in scan_configs:
        scan = CombineToolWrapper.CombineScan(config)
        scan.run()


########################################
# t2ws
########################################

def t2ws(args, card, extraOptions=None):
    modelName = None
    suffix = None
    
    if args.lumiScale:
        modelName = 'lumiScaleDifferentialModel'
        suffix = 'lumiScale'
    
    Commands.basic_t2ws_with_model(
        card,
        pathToModel = 'physicsModels/DifferentialModel.py',
        modelName = modelName,
        suffix = suffix,
        extraOptions = extraOptions
        )

@flag_as_option
def t2ws_all(args):
    pth_smH_t2ws(args)
    pth_ggH_t2ws(args)
    njets_t2ws(args)
    rapidity_t2ws(args)
    ptjet_t2ws(args)    

@flag_as_option
def pth_smH_t2ws(args):
    card = LatestPathsGetters.get_card('pth_smH', args)
    extraOptions = []
    if args.hzz:
        extraOptions.append('--PO \'binning=0,15,30,85,200\'')
    t2ws(args, card, extraOptions)

@flag_as_option
def pth_ggH_t2ws(args):
    card = LatestPathsGetters.get_card('pth_ggH', args)
    extraOptions = []
    if args.hzz:
        extraOptions.append('--PO \'binning=0,15,30,85,200\'')
    t2ws(args, card, extraOptions)

@flag_as_option
def njets_t2ws(args):
    card = LatestPathsGetters.get_card('njets', args)
    extraOptions = []
    if args.hzz:
        extraOptions.append('--PO \'binning=0,1,2,3\'')
    t2ws(args, card, extraOptions)

@flag_as_option
def rapidity_t2ws(args):
    card = LatestPathsGetters.get_card('rapidity', args)
    t2ws(args, card, extraOptions=None)

@flag_as_option
def ptjet_t2ws(args):
    card = LatestPathsGetters.get_card('ptjet', args)
    extraOptions = []
    if args.hzz:
        raise NotImplementedError(
            '--hzz for ptjet does not work yet!!'
            )
        extraOptions.extend([
            '--PO \'binning=30,55,95\'',
            '--PO \'make_underflow\''
            ])
    t2ws(args, card, extraOptions)


########################################
# Preprocessing
########################################

#____________________________________________________________________
# Test of new bins for hzz and hbb

@flag_as_option
def RenumberHzzProcesses_Jan24(args):
    MergeHGGWDatacards.renumber_processes_hzz_aug21(
        LatestPaths.card_hzz_ggHxH_PTH_newBins_unprocessed,
        )
@flag_as_option
def CombineCards_Jan24_hzz_hbb(args):
    Commands.basic_combine_cards(
        'suppliedInput/combinedCard_newBins_hzz_hbb_ggHxH_{0}.txt'.format(datestr),
        'hzz=' + LatestPaths.card_hzz_ggHxH_PTH_newBins,
        'hbb=' + LatestPaths.card_hbb_ggHxH_PTH
        )

#____________________________________________________________________
# pth_smH

@flag_as_option
def RenameHggProcesses_smHcard(args):
    MergeHGGWDatacards.rename_processes_aug21(
        LatestPaths.card_hgg_smH_PTH_unprocessed,
        )
@flag_as_option
def RenumberHzzProcesses_smHcard(args):
    MergeHGGWDatacards.renumber_processes_hzz_aug21(
        LatestPaths.card_hzz_smH_PTH_unprocessed,
        )
@flag_as_option
def CombineCards_smHcard(args):
    Commands.basic_combine_cards(
        'suppliedInput/combinedCard_smH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_smH_PTH,
        'hzz=' + LatestPaths.card_hzz_smH_PTH
        )

#____________________________________________________________________
# pth_ggH

@flag_as_option
def RenameHggProcesses_Aug21(args):
    MergeHGGWDatacards.rename_processes_aug21(
        LatestPaths.card_hgg_ggHxH_PTH_unprocessed,
        )
@flag_as_option
def RenumberHzzProcesses_Aug21(args):
    MergeHGGWDatacards.renumber_processes_hzz_aug21(
        LatestPaths.card_hzz_ggHxH_PTH_unprocessed,
        )
@flag_as_option
def CombineCards_Aug21(args):
    Commands.basic_combine_cards(
        'suppliedInput/combinedCard_ggHxH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_ggHxH_PTH,
        'hzz=' + LatestPaths.card_hzz_ggHxH_PTH
        )
@flag_as_option
def CombineCards_Dec15_hbb(args):
    Commands.basic_combine_cards(
        'suppliedInput/combinedCard_hgg_hzz_hbb_ggHxH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_ggHxH_PTH,
        'hzz=' + LatestPaths.card_hzz_ggHxH_PTH,
        'hbb=' + LatestPaths.card_hbb_ggHxH_PTH
        )

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
    Commands.basic_combine_cards(
        'suppliedInput/combinedCard_NJ_smH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_smH_NJ,
        'hzz=' + LatestPaths.card_hzz_smH_NJ
        )

#____________________________________________________________________
# ptjet

@flag_as_option
def ptjet_rename_hgg(args):
    MergeHGGWDatacards.rename_processes_hgg_differentials(
        LatestPaths.card_hgg_smH_PTJ_unprocessed
        )

@flag_as_option
def ptjet_combineCards(args):
    Commands.basic_combine_cards(
        'suppliedInput/combinedCard_PTJ_smH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_smH_PTJ,
        'hzz=' + LatestPaths.card_hzz_smH_PTJ
        )

#____________________________________________________________________
# rapidity

@flag_as_option
def rename_hgg_rapidity(args):
    MergeHGGWDatacards.rename_processes_hgg_differentials(
        LatestPaths.card_hgg_smH_YH_unprocessed,
        )

@flag_as_option
def rapidity_combineCards(args):
    Commands.basic_combine_cards(
        'suppliedInput/combinedCard_YH_smH_{0}.txt'.format(datestr),
        'hgg=' + LatestPaths.card_hgg_smH_YH,
        'hzz=' + LatestPaths.card_hzz_smH_YH
        )


########################################
# Correlation matrices
########################################

@flag_as_option
def all_corrMats(args):
    corrMat_pth_smH(args)
    corrMat_pth_ggH(args)
    corrMat_njets(args)
    corrMat_ptjet(args)
    corrMat_rapidity(args)

@flag_as_option
def corrMat_pth_smH(args):
    ws = LatestPathsGetters.get_ws('pth_smH', args)
    Commands.compute_corr_matrix( ws, asimov=args.asimov, onBatch = False )

@flag_as_option
def corrMat_pth_ggH(args):
    ws = LatestPathsGetters.get_ws('pth_ggH', args)
    Commands.compute_corr_matrix( ws, asimov=args.asimov, onBatch = False )

@flag_as_option
def corrMat_njets(args):
    ws = LatestPathsGetters.get_ws('njets', args)
    Commands.compute_corr_matrix( ws, asimov=args.asimov, onBatch = False )

@flag_as_option
def corrMat_ptjet(args):
    ws = LatestPathsGetters.get_ws('ptjet', args)
    Commands.compute_corr_matrix( ws, asimov=args.asimov, onBatch = False )

@flag_as_option
def corrMat_rapidity(args):
    ws = LatestPathsGetters.get_ws('rapidity', args)
    Commands.compute_corr_matrix( ws, asimov=args.asimov, onBatch = False )

@flag_as_option
def plotCorrelationMatrices(args):

    pth_smH = Container()
    # pth_smH.ws           = LatestPaths.ws_combined_smH
    pth_smH.ws           = LatestPathsGetters.get_ws('pth_smH', args)
    pth_smH.corrRootFile = LatestPaths.correlationMatrix_PTH
    pth_smH.xTitle       = 'p_{T}^{H} (GeV)'
    PlotCommands.plot_correlation_matrix( pth_smH )

    pth_ggH = Container()
    # pth_ggH.ws           = LatestPaths.ws_combined_ggH_xHfixed
    pth_ggH.ws           = LatestPathsGetters.get_ws('pth_ggH', args)
    pth_ggH.corrRootFile = LatestPaths.correlationMatrix_PTH_ggH
    pth_ggH.xTitle       = 'p_{T}^{H} (GeV) (non-ggH fixed to SM)'
    PlotCommands.plot_correlation_matrix( pth_ggH )

    njets = Container()
    # njets.ws           = LatestPaths.ws_combined_smH_NJ
    njets.ws           = LatestPathsGetters.get_ws('njets', args)
    njets.corrRootFile = LatestPaths.correlationMatrix_NJ
    njets.xTitle       = 'N_{jets}'
    PlotCommands.plot_correlation_matrix( njets )

    rapidity = Container()
    # rapidity.ws           = LatestPaths.ws_combined_smH_YH
    rapidity.ws           = LatestPathsGetters.get_ws('rapidity', args)
    rapidity.corrRootFile = LatestPaths.correlationMatrix_YH
    rapidity.xTitle       = '|y_{H}|'
    PlotCommands.plot_correlation_matrix( rapidity )

    ptjet = Container()
    # ptjet.ws           = LatestPaths.ws_combined_smH_PTJ
    ptjet.ws           = LatestPathsGetters.get_ws('ptjet', args)
    ptjet.corrRootFile = LatestPaths.correlationMatrix_PTJ
    ptjet.xTitle       = 'p_{T}^{j1}'
    PlotCommands.plot_correlation_matrix( ptjet )

