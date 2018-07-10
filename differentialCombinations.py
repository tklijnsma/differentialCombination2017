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

import logging
from OptionHandler import flag_as_option

import differentials
import differentialutils
import LatestPaths
# import LatestPathsGetters
import LatestBinning

from differentials.core import AttrDict
import differentials.combine.combine as combine

# sys.path.append('src')
# import CombineToolWrapper

from time import strftime
datestr = strftime( '%b%d' )


########################################
# Main
########################################

rescans = AttrDict()
rescans.pth_ggH = AttrDict()
rescans.pth_ggH.observed = [
    # AttrDict(dc='combWithHbb', POI='r_ggH_PTH_350_600', x_min=-10., x_max=10.),
    # AttrDict(dc='combWithHbb', POI='r_ggH_PTH_GT600', x_min=-5., x_max=10.),
    # AttrDict(dc='hgg', POI='r_ggH_PTH_350_600', x_min=-10., x_max=1.),
    # AttrDict(dc='hgg', POI='r_ggH_PTH_GT600', x_min=-110., x_max=-70.),
    # AttrDict(dc='hzz', POI='r_ggH_PTH_GT200', x_min=-10., x_max=1.),
    ]
rescans.pth_ggH.asimov = [
    # AttrDict(dc='combWithHbb', POI='r_ggH_PTH_GT600', x_min=-10., x_max=10.),
    # AttrDict(dc='hgg', POI='r_ggH_PTH_GT600', x_min=-10., x_max=20.),
    # AttrDict(dc='hzz', POI='r_ggH_PTH_GT200', x_min=-1., x_max=4.),
    ]
rescans.pth_smH = AttrDict()
rescans.pth_smH.observed = [
    # AttrDict(dc='combWithHbb', POI='r_smH_PTH_350_600', x_min=-10., x_max=10.),
    # AttrDict(dc='combWithHbb', POI='r_smH_PTH_GT600', x_min=-10., x_max=10.),
    # AttrDict(dc='hgg', POI='r_smH_PTH_GT600', x_min=-10., x_max=4.),
    # AttrDict(dc='hzz', POI='r_smH_PTH_GT200', x_min=-10., x_max=1.),
    ]
rescans.pth_smH.asimov = [
    # AttrDict(dc='combWithHbb', POI='r_smH_PTH_GT600', x_min=-5., x_max=5.),
    # AttrDict(dc='hgg', POI='r_smH_PTH_GT600', x_min=-3., x_max=8.),
    ]

@flag_as_option
def rescan(real_args):
    for asimov in [False, True]:
        for obsname in ['pth_ggH', 'pth_smH']:
            for rescan in rescans[obsname]['asimov' if asimov else 'observed']:
                args = differentialutils.set_one_decay_channel(real_args,
                    rescan.dc,
                    asimov=asimov
                    )
                ws = LatestPaths.ws[obsname][rescan.dc]
                config = differential_config(args, ws, obsname)

                config.POIs = [rescan.POI]
                config.nPoints = 40

                for i in xrange(len(config.PhysicsModelParameterRanges)):
                    range_str = config.PhysicsModelParameterRanges[i]
                    if range_str.startswith(rescan.POI):
                        config.PhysicsModelParameterRanges[i] = '{0}={1},{2}'.format(
                            rescan.POI, rescan.x_min, rescan.x_max
                            )
            
                config.tags.append('rescan')
                scan_directly(args, config)
                # postfit_and_scan(args, config)


@flag_as_option
def pth_smH_scan(args):
    ws = LatestPaths.ws.pth_smH[differentialutils.get_decay_channel_tag(args)]
    config = differential_config(args, ws, 'pth_smH')
    postfit_and_scan(args, config)

@flag_as_option
def pth_ggH_scan(args):
    ws = LatestPaths.ws.pth_ggH[differentialutils.get_decay_channel_tag(args)]
    config = differential_config(args, ws, 'pth_ggH')
    postfit_and_scan(args, config)

@flag_as_option
def njets_scan(args):
    ws = LatestPaths.ws.njets[differentialutils.get_decay_channel_tag(args)]
    config = differential_config(args, ws, 'njets')
    postfit_and_scan(args, config)

@flag_as_option
def ptjet_scan(args):
    ws = LatestPaths.ws.ptjet[differentialutils.get_decay_channel_tag(args)]
    config = differential_config(args, ws, 'ptjet')
    postfit_and_scan(args, config)

@flag_as_option
def rapidity_scan(args):
    ws = LatestPaths.ws.rapidity[differentialutils.get_decay_channel_tag(args)]
    config = differential_config(args, ws, 'rapidity')
    postfit_and_scan(args, config)


@flag_as_option
def pth_ggH_scan_hbb_floatingOOA(args):
    args = differentialutils.set_one_decay_channel(args, 'hbb', asimov=True)
    ws = LatestPaths.ws.pth_ggH['hbb']
    config = differential_config(args, ws, 'pth_ggH')
    config.tags.append('floatingOOA')
    postfit_and_scan(args, config)

@flag_as_option
def pth_ggH_scan_hbb_fixedOOA(args):
    args = differentialutils.set_one_decay_channel(args, 'hbb', asimov=True)
    ws = LatestPaths.ws.pth_ggH['hbb']
    config = differential_config(args, ws, 'pth_ggH')
    config.POIs.pop(config.POIs.index('r_ggH_PTH_200_350'))
    config.FLOAT_OTHER_POIS = False
    config.del_parameter_range('r_ggH_PTH_200_350')
    config.freezeNuisances.append('r_ggH_PTH_200_350')
    config.tags.append('fixedOOA')
    postfit_and_scan(args, config)


#____________________________________________________________________
# Helpers

lumiMultiplier300 = 8.356546
lumiMultiplier3000 = 83.56546
lumiMultiplier = lumiMultiplier300

def differential_config(args, ws, obs_name):
    base_config = combine.CombineConfig(args)
    base_config.onBatch       = True
    base_config.nPoints       = 55
    base_config.nPointsPerJob = 5
    base_config.queue         = 'short.q'

    if args.asimov:
        base_config.asimov = True
    else:
        base_config.asimov = False

    base_config.decay_channel = differentialutils.get_decay_channel_tag(args)    
    if args.hzz or args.hbb:
        base_config.nPointsPerJob = base_config.nPoints
    if args.combWithHbb:
        base_config.nPoints = 100
        base_config.nPointsPerJob = 4
    if args.hbb or args.combWithHbb:
        base_config.minimizer_settings = [
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            ]

    base_config.datacard = ws
    base_config.subDirectory = 'out/Scan_{0}_{1}_{2}'.format(datestr, obs_name, base_config.decay_channel)
    if args.lumiScale:
        base_config.hardPhysicsModelParameters.append('lumiScale={0}'.format(lumiMultiplier))
        base_config.freezeNuisances.append('lumiScale')
        base_config.subDirectory += '_lumiScale{0}'.format(str(lumiMultiplier).replace('.','p'))
    if args.statonly:
        base_config.subDirectory += '_statonly'

    POIs = differentials.core.list_POIs(ws)
    POIrange = [ -1.0, 4.0 ]
    if args.hzz:
        POIrange[0] = 0.0
    if args.hbb:
        POIrange = [ -10.0, 10.0 ]
    if args.combWithHbb:
        POIrange = [ -1.0, 6.0 ]
    if args.lumiScale:
        POIrange = [ 0.7, 1.3 ]
        base_config.nPoints = 35
        if lumiMultiplier > 10.:
            POIrange = [ 0.9, 1.1 ]

    base_config.POIs = POIs
    base_config.PhysicsModelParameters = [ '{0}=1.0'.format(p) for p in POIs ]

    for POI in POIs:
        left, right = POIrange
        if '350_600' in POI or 'GT600' in POI or 'GT200' in POI:
            left = -10.
            right = 10.
            if base_config.decay_channel == 'hgg':
                left = -15.
        base_config.PhysicsModelParameterRanges.append('{0}={1},{2}'.format(POI, left, right))

    if args.hbb or args.combWithHbb:
        logging.warning(
            'For hbb: qcdeff and r*p* parameters are given ranges! '
            'Should decide on (freezing) behaviour during postfit/scan'
            )
        base_config.PhysicsModelParameterRanges.extend([
                'qcdeff=0.001,8.0',
                'r1p0=0.0,8.0',
                'r2p0=0.0,8.0',
                'r3p0=0.0,8.0',
                'r0p1=0.0,8.0',
                'r1p1=0.0,8.0',
                'r2p1=0.0,8.0',
                'r3p1=0.0,8.0',
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
        logging.warning('Disabling output to limit logging size')
        config.extraOptions.append('-v -1')

        configs_per_POI.append(config)
    return configs_per_POI

def postfit_and_scan(args, postfit_config):
    # Make sure no previous run directory is overwritten
    postfit_config.make_unique_directory()

    postfit = combine.CombinePostfit(postfit_config)
    postfit.run()
    postfit_ws = postfit.get_output()

    scan_configs = differential_configs_for_scanning(args, postfit_config)
    for config in scan_configs:
        scan = combine.CombineScanFromPostFit(config)
        scan.run(postfit_ws)

def scan_directly(args, base_config):
    base_config.make_unique_directory()
    scan_configs = differential_configs_for_scanning(args, base_config)
    for config in scan_configs:
        scan = combine.CombineScan(config)
        scan.run()

def only_postfit(args, config):
    postfit = combine.CombineOnlyPostfit(config)
    postfit.run()


########################################
# t2ws
########################################

def basic_t2ws(obsname, decay_channel):
    t2ws = differentials.combine.t2ws.T2WS(
        LatestPaths.card[obsname][decay_channel],
        # model_file='physicsModels/DifferentialModel.py', model_name='differentialModel',
        name = 'ws_' + obsname + '_' + decay_channel
        )
    return t2ws

# @flag_as_option
# def pth_ggH_t2ws(args):
#     t2ws = basic_t2ws('pth_ggH', differentialutils.get_decay_channel_tag(args))
#     if args.hzz:
#         t2ws.extra_options.append('--PO \'binning=0,15,30,80,200\'')
#     if args.hbb:
#         t2ws.extra_options.append('--PO \'binning=200,350,600\'')
#         # t2ws.extra_options.append('--PO \'binning=350,600\'')
#     t2ws.run()

@flag_as_option
def pth_ggH_t2ws(args):
    t2ws = basic_t2ws('pth_ggH', differentialutils.get_decay_channel_tag(args))
    if args.hzz:
        t2ws.make_maps_from_processes(add_overflow=True)
    else:
        t2ws.make_maps_from_processes()
    t2ws.run()


@flag_as_option
def pth_smH_t2ws(args):
    raise NotImplementedError('Do the same thing as for ggH')
    t2ws = basic_t2ws('pth_smH', differentialutils.get_decay_channel_tag(args))
    if args.hzz:
        t2ws.extra_options.append('--PO \'binning=0,15,30,80,200\'')
    if args.hbb:
        t2ws.extra_options.append('--PO \'binning=200,350,600\'')
    if args.hgg or args.combination or args.combWithHbb or args.hbb:
        t2ws.extra_options.append('--PO \'scale_xH_ggH_as_smH=True\'')
    t2ws.run()

@flag_as_option
def njets_t2ws(args):
    t2ws = basic_t2ws('njets', differentialutils.get_decay_channel_tag(args))
    if args.hzz:
        t2ws.make_maps_from_processes(binning=[0,1,2], add_overflow=True)
    else:
        t2ws.make_maps_from_processes()
    t2ws.run()

@flag_as_option
def rapidity_t2ws(args):
    t2ws = basic_t2ws('rapidity', differentialutils.get_decay_channel_tag(args))
    t2ws.make_maps_from_processes()
    t2ws.run()

@flag_as_option
def ptjet_t2ws(args):
    t2ws = basic_t2ws('ptjet', differentialutils.get_decay_channel_tag(args))
    if args.hzz:
        t2ws.make_maps_from_processes(binning=[30,55,95], add_underflow=True, add_overflow=True)
    else:
        t2ws.make_maps_from_processes()
    t2ws.run()



# TODO: Reform t2ws for pth_smH, njets, ptjet, rapidity

# @flag_as_option
# def t2ws_all(args):
#     pth_smH_t2ws(args)
#     pth_ggH_t2ws(args)
#     njets_t2ws(args)
#     rapidity_t2ws(args)
#     ptjet_t2ws(args)    

# @flag_as_option
# def pth_smH_t2ws(args):
#     card = LatestPathsGetters.get_card('pth_smH', args)
#     extraOptions = []
#     if args.hzz:
#         extraOptions.append('--PO \'binning=0,15,30,85,200\'')
#     t2ws(args, card, extraOptions)

# @flag_as_option
# def pth_ggH_t2ws(args):
#     card = LatestPaths.card.pth_ggH[get_decay_channel_tag(args)]
#     extraOptions = []
#     if args.hzz:
#         extraOptions.append('--PO \'binning=0,15,30,85,200\'')
#     t2ws(args, card, extraOptions)

# @flag_as_option
# def njets_t2ws(args):
#     card = LatestPathsGetters.get_card('njets', args)
#     extraOptions = []
#     if args.hzz:
#         extraOptions.append('--PO \'binning=0,1,2,3\'')
#     t2ws(args, card, extraOptions)

# @flag_as_option
# def rapidity_t2ws(args):
#     card = LatestPathsGetters.get_card('rapidity', args)
#     t2ws(args, card, extraOptions=None)

# @flag_as_option
# def ptjet_t2ws(args):
#     card = LatestPathsGetters.get_card('ptjet', args)
#     extraOptions = []
#     if args.hzz:
#         raise NotImplementedError(
#             '--hzz for ptjet does not work yet!!'
#             )
#         extraOptions.extend([
#             '--PO \'binning=30,55,95\'',
#             '--PO \'make_underflow\''
#             ])
#     t2ws(args, card, extraOptions)


########################################
# Correlation matrices
########################################

sys.path.append('src')
import Commands
import PlotCommands

@flag_as_option
def all_corrMats(args):
    corrMat_pth_smH(args)
    corrMat_pth_ggH(args)
    corrMat_njets(args)
    corrMat_ptjet(args)
    corrMat_rapidity(args)

@flag_as_option
def corrMat_pth_smH(args):
    # ws = LatestPathsGetters.get_ws('pth_smH', args)
    ws = LatestPaths.ws.pth_smH.combWithHbb
    Commands.ComputeCorrMatrix( ws, asimov=args.asimov, onBatch = False )

@flag_as_option
def corrMat_pth_ggH(args):
    # ws = LatestPathsGetters.get_ws('pth_ggH', args)
    ws = LatestPaths.ws.pth_ggH.combWithHbb
    Commands.ComputeCorrMatrix( ws, asimov=args.asimov, onBatch = False )

@flag_as_option
def corrMat_njets(args):
    ws = LatestPathsGetters.get_ws('njets', args)
    Commands.ComputeCorrMatrix( ws, asimov=args.asimov, onBatch = False )

@flag_as_option
def corrMat_ptjet(args):
    ws = LatestPathsGetters.get_ws('ptjet', args)
    Commands.ComputeCorrMatrix( ws, asimov=args.asimov, onBatch = False )

@flag_as_option
def corrMat_rapidity(args):
    ws = LatestPathsGetters.get_ws('rapidity', args)
    Commands.ComputeCorrMatrix( ws, asimov=args.asimov, onBatch = False )

@flag_as_option
def plotCorrelationMatrices(args):
    # pth_smH = Container()
    pth_smH = differentials.core.AttrDict()
    # pth_smH.ws           = LatestPathsGetters.get_ws('pth_smH', args)
    pth_smH.ws = 'corrMat_Jun04_ws_pth_smH_combWithHbb/higgsCombine_POSTFIT_ws_pth_smH_combWithHbb.MultiDimFit.mH125.root'
    # pth_smH.corrRootFile = LatestPaths.correlationMatrix_PTH
    pth_smH.corrRootFile = 'corrMat_Jun04_ws_pth_smH_combWithHbb/higgsCombine_CORRMAT_ws_pth_smH_combWithHbb.MultiDimFit.mH125.root'
    pth_smH.xTitle       = 'p_{T}^{H} (GeV)'
    PlotCommands.PlotCorrelationMatrix( pth_smH )

    # pth_ggH = Container()
    pth_ggH = differentials.core.AttrDict()
    # pth_ggH.ws           = LatestPathsGetters.get_ws('pth_ggH', args)
    pth_ggH.ws = 'corrMat_Jun04_ws_pth_ggH_combWithHbb/higgsCombine_POSTFIT_ws_pth_ggH_combWithHbb.MultiDimFit.mH125.root'
    # pth_ggH.corrRootFile = LatestPaths.correlationMatrix_PTH_ggH
    pth_ggH.corrRootFile = 'corrMat_Jun04_ws_pth_ggH_combWithHbb/higgsCombine_CORRMAT_ws_pth_ggH_combWithHbb.MultiDimFit.mH125.root'
    pth_ggH.xTitle       = 'p_{T}^{H} (GeV) (ggH)'
    PlotCommands.PlotCorrelationMatrix( pth_ggH )

    # njets = Container()
    # # njets.ws           = LatestPaths.ws_combined_smH_NJ
    # njets.ws           = LatestPathsGetters.get_ws('njets', args)
    # njets.corrRootFile = LatestPaths.correlationMatrix_NJ
    # njets.xTitle       = 'N_{jets}'
    # PlotCommands.PlotCorrelationMatrix( njets )

    # rapidity = Container()
    # # rapidity.ws           = LatestPaths.ws_combined_smH_YH
    # rapidity.ws           = LatestPathsGetters.get_ws('rapidity', args)
    # rapidity.corrRootFile = LatestPaths.correlationMatrix_YH
    # rapidity.xTitle       = '|y_{H}|'
    # PlotCommands.PlotCorrelationMatrix( rapidity )

    # ptjet = Container()
    # # ptjet.ws           = LatestPaths.ws_combined_smH_PTJ
    # ptjet.ws           = LatestPathsGetters.get_ws('ptjet', args)
    # ptjet.corrRootFile = LatestPaths.correlationMatrix_PTJ
    # ptjet.xTitle       = 'p_{T}^{j1}'
    # PlotCommands.PlotCorrelationMatrix( ptjet )

