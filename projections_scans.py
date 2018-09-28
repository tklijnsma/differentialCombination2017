#!/usr/bin/env python
"""
Thomas Klijnsma
"""

#____________________________________________________________________
# Imports

import logging
import os, sys, re, copy
from OptionHandler import flag_as_option
import differentials
import differentialutils
import LatestBinning
import LatestPaths

lumiscale300 = 8.356546
lumiscale3000 = 83.56546


#____________________________________________________________________
scenario1 = differentials.core.AttrDict()
scenario1.hgg = 'projections/workspaces_Jun28/ws_pth_smH_hgg.root'
scenario1.hzz = 'projections/workspaces_Jun28/ws_pth_smH_hzz.root'
scenario1.hbb = 'projections/workspaces_Jun28/ws_pth_smH_hbb.root'
scenario1.combWithHbb = 'projections/workspaces_Jun28/ws_pth_smH_combWithHbb.root'
scenario1.hbb_postfit_300ifb = 'out/postfits_Jul03/higgsCombine_POSTFIT_ASIMOV_ws_pth_smH_hbb.MultiDimFit.mH125.root'

scenario2 = differentials.core.AttrDict()
scenario2.hgg = 'projections/workspaces_Jul17/ws_pth_smH_hgg_s2.root'
scenario2.hzz = 'projections/workspaces_Jul17/ws_pth_smH_hzz_s2.root'
scenario2.hbb = 'projections/workspaces_Jul17/ws_pth_smH_hbb_s2.root'
scenario2.combWithHbb = 'projections/workspaces_Jul17/ws_pth_smH_combWithHbb_s2.root'


@flag_as_option
def projection_pth_smH_scan(args):
    args = differentialutils.force_asimov(args)
    decay_channel = differentialutils.get_decay_channel_tag(args)
    wsdict = scenario2 if args.scenario2 else scenario1
    ws = wsdict[decay_channel]

    config = differential_config(args, ws, 'pth_smH', lumiscale=lumiscale3000)

    def tight_ranges(config):
        """Very manual tight r_ ranges based on Hgg scan results"""
        for POI in config.POIs:
            if 'GT600' in POI:
                config.set_parameter_range(POI, 0.64, 1.4)
            else:
                config.set_parameter_range(POI, 0.875, 1.14)
    
    verbosity = 1
    if args.hbb:
        set_hbb_parameter_ranges(args, config)
        # config.saveNuisances.extend([
        #     'CMS_PU', 'CMS_res_j', 'CMS_scale_j', 'bbeff', 'eleveto', 'hqq125GenpT2failcat1mcstat', 'hqq125GenpT2passcat1mcstat', 'hqq125GenpT3failcat1mcstat', 'hqq125GenpT3failcat2mcstat', 'hqq125GenpT3passcat1mcstat', 'hqq125GenpT3passcat2mcstat', 'hqq125GenpT4failcat1mcstat', 'hqq125GenpT4failcat2mcstat', 'hqq125GenpT4passcat1mcstat', 'hqq125GenpT4passcat2mcstat', 'hqq125failmuonCRmcstat', 'hqq125pt', 'hqq125ptShape', 'lumi_13TeV_2016', 'muid', 'muiso', 'mutrigger', 'muveto', 'qcdfailmuonCRmcstat', 'qcdpassmuonCRmcstat', 'scale', 'scalept', 'smear', 'stqqfailmuonCRmcstat', 'stqqpassmuonCRmcstat', 'tqqfailcat1mcstat', 'tqqfailcat2mcstat', 'tqqfailmuonCRmcstat', 'tqqpasscat1mcstat', 'tqqpasscat2mcstat', 'tqqpassmuonCRmcstat', 'trigger', 'veff', 'vvqqfailmuonCRmcstat', 'wlnufailmuonCRmcstat', 'wlnupassmuonCRmcstat', 'wqqfailcat1mcstat', 'wqqfailcat2mcstat', 'wqqpasscat1mcstat', 'wqqpasscat2mcstat', 'wznormEW', 'xhqq125GenpT1failcat1mcstat', 'xhqq125GenpT1failcat2mcstat', 'xhqq125GenpT1passcat1mcstat', 'xhqq125GenpT1passcat2mcstat', 'xhqq125GenpT2failcat1mcstat', 'xhqq125GenpT2failcat2mcstat', 'xhqq125GenpT2passcat1mcstat', 'xhqq125GenpT2passcat2mcstat', 'xhqq125GenpT3failcat1mcstat', 'xhqq125GenpT3failcat2mcstat', 'xhqq125GenpT3passcat1mcstat', 'xhqq125GenpT3passcat2mcstat', 'xhqq125GenpT4failcat1mcstat', 'xhqq125GenpT4failcat2mcstat', 'xhqq125GenpT4passcat1mcstat', 'xhqq125GenpT4passcat2mcstat', 'xhqq125failmuonCRmcstat', 'xhqq125passmuonCRmcstat', 'zllfailmuonCRmcstat', 'znormEW', 'znormQ', 'zqqfailcat1mcstat', 'zqqfailcat2mcstat', 'zqqfailmuonCRmcstat', 'zqqpasscat1mcstat', 'zqqpasscat2mcstat', 'tqqfailcat1norm', 'tqqfailcat2norm', 'tqqfailmuonCRnorm', 'tqqpasscat1norm', 'tqqpasscat2norm', 'tqqpassmuonCRnorm', 'tqqnormSF', 'tqqeffSF', 'lumiscale'
        #     ])
        # verbosity = 3
    if args.hgg:
        tight_ranges(config)
    if args.combWithHbb:
        set_hbb_parameter_ranges(args, config)
        tight_ranges(config)

    if args.hbb or args.hzz:
        scan_directly(args, config, verbosity)
    else:
        postfit_and_scan(args, config)


@flag_as_option
def projection_pth_smH_scan_one_bin_locally(args):
    decay_channel = 'combWithHbb'
    args = differentialutils.set_one_decay_channel(args, decay_channel, asimov=True)
    ws = workspaces[decay_channel]
    config = differential_config(args, ws, 'pth_smH')

    POI = 'r_smH_PTH_0_15'
    left = 0.9
    right = 1.1
    n_points = 8

    config.onBatch = False
    config.POIs = [POI]
    config.set_parameter_range(POI, left, right)
    scan_directly(args, config, verbosity=1)




@flag_as_option
def projection_postfit_hbb_300ifb(args):
    decay_channel = 'hbb'
    args = differentialutils.set_one_decay_channel(args, decay_channel, asimov=True)
    ws = scenario1[decay_channel]
    config = differential_config(args, ws, 'pth_smH', lumiscale=lumiscale300)
    only_postfit(args, config)

def get_ranges_from_hbb_postfit_300ifb(args):
    w = differentials.core.get_ws(scenario1.hbb_postfit_300ifb)
    w.loadSnapshot('MultiDimFit')
    keys = [
        'qcdeff',
        'r0p1',
        'r1p0',
        'r1p1',
        'r2p0',
        'r2p1',
        'r3p0',
        'r3p1',
        ]
    out = {}
    for key in keys:
        out[key] = w.var(key).getVal()
    return out

def set_hbb_parameter_ranges(args, config):
    bestfit = get_ranges_from_hbb_postfit_300ifb(args)
    for key, value in bestfit.iteritems():
        config.set_parameter(key, value)
        config.set_parameter_range(key, 0.9*value, 1.1*value)



def differential_config(args, ws, obs_name, lumiscale=1.):
    # base_config = combine.CombineConfig(args)
    base_config = differentials.combine.combine.CombineConfig(args)
    base_config.onBatch       = True
    base_config.queue         = 'short.q'
    base_config.asimov = True if args.asimov else False

    base_config.decay_channel = differentialutils.get_decay_channel_tag(args)    
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
    base_config.subDirectory = 'out/ScanProjection_{0}_{1}_{2}'.format(
        differentials.core.datestr(),
        obs_name,
        base_config.decay_channel
        )

    # Scale lumi
    base_config.freezeNuisances.append('lumiscale')
    base_config.hardPhysicsModelParameters.append('lumiscale={0}'.format(lumiscale))

    if args.hbb or args.combWithHbb:
        # Freeze mcstat nuisances
        logging.info('Freezing mcstat systematics')
        base_config.extraOptions.append('--freezeNuisanceGroups mcstat')

    if lumiscale == 1.:
        base_config.subDirectory += '_36ifb'
    elif lumiscale == lumiscale300:
        base_config.subDirectory += '_300ifb'
    elif lumiscale == lumiscale3000:
        base_config.subDirectory += '_3000ifb'

    if args.statonly:
        base_config.subDirectory += '_statonly'
    if args.scenario2:
        base_config.subDirectory += '_scenario2'

    POIs = differentials.core.list_POIs(ws)
    if args.hbb: POIs.pop(0)
    base_config.POIs = POIs
    base_config.PhysicsModelParameters = [ '{0}=1.0'.format(p) for p in POIs ]

    if args.hbb or args.combWithHbb:
        logging.warning(
            'For hbb: qcdeff and r*p* parameters are given ranges! '
            'Should decide on (freezing) behaviour during postfit/scan'
            )
        base_config.PhysicsModelParameterRanges.extend([
                'qcdeff=0.001,8.0',
                'r1p0=0.1,8.0',
                'r2p0=0.1,8.0',
                'r3p0=0.1,8.0',
                'r0p1=0.1,8.0',
                'r1p1=0.1,8.0',
                'r2p1=0.1,8.0',
                'r3p1=0.1,8.0',
                ])

    base_config.nPoints       = 55
    # base_config.nPoints       = 20
    base_config.nPointsPerJob = 5

    if args.hzz or args.hbb:
        base_config.nPoints       = 90
        base_config.nPointsPerJob = base_config.nPoints

    if args.combWithHbb:
        base_config.nPointsPerJob = 1

    # Ranges
    default_p3000 = differentials.core.AttrDict()
    default_p3000.hzz = [ 0.7, 1.3 ]
    default_p3000.hgg = [ 0.7, 1.3 ]
    default_p3000.hbb = [ 0.3, 1.7 ]
    default_p3000.combWithHbb = [ 0.5, 1.5 ]

    default_ranges = default_p3000
    for POI in POIs:
        ranges = default_ranges[base_config.decay_channel][:]
        if args.hzz and 'GT200' in POI:
            ranges = [ 0.3, 1.7 ]
        if (args.hgg) and ('GT600' in POI or '350_600' in POI):
            ranges = [ 0.4, 1.6 ]
        # if args.combWithHbb and ('GT600' in POI or '350_600' in POI):
        #     ranges = [ 0.6, 1.4 ]
        base_config.set_parameter_range(POI, ranges[0], ranges[1])

    return base_config

def differential_configs_for_scanning(args, original_base_config):
    base_config = copy.deepcopy(original_base_config)

    if args.statonly:
        base_config.freezeNuisances.append('rgx{.*}')

    # Overwrite get_name function to include the POI name
    def new_get_name(self):
        return 'bPOI_{0}_ePOI_{1}'.format(self.POIs[0], os.path.basename(self.datacard).replace('.root',''))
    base_config.get_name = new_get_name.__get__(base_config, type(base_config))

    configs_per_POI = []
    for POI in base_config.POIs:
        config = copy.deepcopy(base_config)
        config.POIs = [POI]
        config.set_verbosity(-1)

        configs_per_POI.append(config)
    return configs_per_POI

def postfit_and_scan(args, postfit_config):
    # Make sure no previous run directory is overwritten
    postfit_config.make_unique_directory()

    postfit = differentials.combine.combine.CombinePostfit(postfit_config)
    postfit.run()
    postfit_ws = postfit.get_output()

    scan_configs = differential_configs_for_scanning(args, postfit_config)
    for config in scan_configs:
        scan = differentials.combine.combine.CombineScanFromPostFit(config)
        scan.run(postfit_ws)

def scan_directly(args, base_config, verbosity=1):
    base_config.make_unique_directory()
    scan_configs = differential_configs_for_scanning(args, base_config)
    for config in scan_configs:
        config.set_verbosity(verbosity)
        scan = differentials.combine.combine.CombineScan(config)
        scan.run()

def only_postfit(args, config):
    postfit = differentials.combine.combine.CombineOnlyPostfit(config)
    postfit.run()

