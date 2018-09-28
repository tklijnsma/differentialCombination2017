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
workspaces = differentials.core.AttrDict()
workspaces.couplingdependentBRs = differentials.core.AttrDict()
workspaces.couplingdependentBRs.scenario2 = differentials.core.AttrDict()
workspaces.floatingBRs = differentials.core.AttrDict()
workspaces.floatingBRs.scenario2 = differentials.core.AttrDict()

def basic_config(args, hurry=False):
    # assert_highpt(args)
    config = combine.CombineConfig(args)
    config.onBatch       = True
    config.queue         = 'all.q'
    config.asimov        = True if args.asimov else False
    config.decay_channel = differentialutils.get_decay_channel_tag(args)

    if args.combWithHbb or args.hbb:
        config.minimizer_settings.extend([
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            ])
    config.deltaNLLCutOff = 70.
    config.POIs = [ 'ct', 'cg' ]
    config.PhysicsModelParameters = [ 'ct=1.0', 'cg=0.0' ]
    config.subDirectory = 'out/Scan_{0}_Top_{1}'.format(datestr, config.decay_channel)
    return config


class CombineConfigKTCGKB(differentials.combine.combine.CombineConfig):
    """docstring for CombineConfigKTCGKB"""
    def __init__(self, args):
        super(CombineConfigKTCGKB, self).__init__(args)
        self.onBatch       = True
        self.queue         = 'all.q'
        self.asimov        = True if args.asimov else False
        self.decay_channel = differentialutils.get_decay_channel_tag(args)
        if args.combWithHbb or args.hbb:
            self.minimizer_settings.extend([
                '--minimizerStrategy 2',
                '--minimizerTolerance 0.001',
                '--robustFit 1',
                '--minimizerAlgoForMinos Minuit2,Migrad',
                ])
        self.lumiscale = 83.56546
        self.hardPhysicsModelParameters.append('lumiscale={0}'.format(self.lumiscale))
        self.deltaNLLCutOff = 150.

    def configure_ktcg(self):
        self.POIs = [ 'ct', 'cg' ]
        self.PhysicsModelParameters = [ 'ct=1.0', 'cg=0.0' ]
        self.subDirectory = 'out/Scan_projection_ktcg_{0}_{1}'.format(differentials.core.datestr(), self.decay_channel)

    def configure_ktkb(self):
        self.POIs = [ 'ct', 'cb' ]
        self.PhysicsModelParameters = [ 'ct=1.0', 'cb=1.0' ]
        self.subDirectory = 'out/Scan_projection_ktkb_{0}_{1}'.format(differentials.core.datestr(), self.decay_channel)

    def range_couplingdependentBRs(self, args):
        if args.hzz:
            self.range_hzz_couplingdependentBRs(args)
        elif args.hgg:
            self.range_hgg_couplingdependentBRs(args)
        elif args.combWithHbb:
            self.range_combWithHbb_couplingdependentBRs(args)

    def range_hzz_couplingdependentBRs(self, args):
        self.nPoints = 100*100
        self.nPointsPerJob = 200
        self.set_parameter_range('ct', -0.3, 3.8)
        self.set_parameter_range('cg', -0.15, 0.08)

    def range_hgg_couplingdependentBRs(self, args):
        self.nPoints = 60*60
        self.nPointsPerJob = 10
        if args.scenario2:
            self.set_parameter_range('ct', 0.7, 1.4)
            self.set_parameter_range('cg', -0.02, 0.015)
        else:
            self.set_parameter_range('ct', 0.6, 1.5)
            self.set_parameter_range('cg', -0.035, 0.035)

    def range_combWithHbb_couplingdependentBRs(self, args):
        self.nPoints = 60*60
        self.nPointsPerJob = 10
        if args.scenario2:
            self.set_parameter_range('ct', 0.7, 1.4)
            self.set_parameter_range('cg', -0.02, 0.015)
        else:
            self.set_parameter_range('ct', 0.7, 1.4)
            self.set_parameter_range('cg', -0.02, 0.015)

        # config.nPoints = 30*30
        # # On combWithHbb, observed: around 4 hours for 10 points, but lots of variance
        # # Good maximum per job is 10 (tiny bit of loss)
        # config.nPointsPerJob = 10
        # config.set_parameter_range('ct', -0.4, 4.4)
        # config.set_parameter_range('cg', -0.28, 0.1)

        # if args.hgg:
        #     config.set_parameter_range('ct', -0.3, 2.6)
        #     config.set_parameter_range('cg', -0.08, 0.12)
        # elif args.hzz:
        #     config.set_parameter_range('ct', -0.8, 3.2)
        #     config.set_parameter_range('cg', -0.15, 0.16)
        #     config.nPointsPerJob = 300
        #     config.queue = 'short.q'

    def range_floatingBRs(self, args):
        if args.hzz:
            self.range_hzz_floatingBRs(args)
        elif args.hgg:
            self.range_hgg_floatingBRs(args)
        elif args.combWithHbb:
            self.range_combWithHbb_floatingBRs(args)

    def range_hzz_floatingBRs(self, args):
        self.nPoints = 40*40
        self.nPointsPerJob = 200
        self.set_parameter_range('ct', -6.0, 6.0) 
        self.set_parameter_range('cg', -0.45, 0.45)

    def range_hgg_floatingBRs(self, args):
        self.nPoints = 40*40
        self.nPointsPerJob = 12
        self.set_parameter_range('ct', -6.0, 6.0) 
        self.set_parameter_range('cg', -0.45, 0.45)

    def range_combWithHbb_floatingBRs(self, args):
        self.range_hgg_floatingBRs()

        # config.nPoints = 30*30
        # Around 4.5 hours (good max of 6 hours) for 8 points (combWithHbb, observed)
        # Reasonably constant, can go up to 12 points
        # Keep same ranges for hgg/hzz for now
        # config.nPointsPerJob = 12
        # config.set_parameter_range('ct', -4.2, 4.2) 
        # config.set_parameter_range('cg', -0.30, 0.30)

        # config.nPoints = 100*100
        # config.nPointsPerJob = 12
        # config.set_parameter_range('ct', -6.0, 6.0) 
        # config.set_parameter_range('cg', -0.45, 0.45)

        # if args.hzz:
        #     config.nPointsPerJob = 300
        #     config.queue = 'short.q'


workspaces.couplingdependentBRs.hgg         = 'out/workspaces_Jul10/projection_ktcgkb_hgg_Jul10_couplingModel_reweighted_couplingdependentBRs.root'
workspaces.couplingdependentBRs.hzz         = 'out/workspaces_Jul10/projection_ktcgkb_hzz_Jul10_couplingModel_reweighted_couplingdependentBRs.root'
workspaces.couplingdependentBRs.combWithHbb = 'out/workspaces_Jul10/projection_ktcgkb_combWithHbb_Jul10_couplingModel_reweighted_couplingdependentBRs.root'

workspaces.couplingdependentBRs.scenario2.hgg         = 'out/workspaces_Jul18/projection_ktcgkb_hgg_s2groups_Jul18_couplingModel_reweighted_scenario2_couplingdependentBRs.root'
workspaces.couplingdependentBRs.scenario2.hzz         = 'out/workspaces_Jul18/projection_ktcgkb_hzz_s2groups_Jul18_couplingModel_reweighted_scenario2_couplingdependentBRs.root'
workspaces.couplingdependentBRs.scenario2.combWithHbb = 'out/workspaces_Jul18/projection_ktcgkb_combWithHbb_s2groups_Jul18_couplingModel_reweighted_scenario2_couplingdependentBRs.root'

@flag_as_option
def projection_ktcg_scan_couplingdependentBRs(args):
    args = differentialutils.force_asimov(args)
    decay_channel = differentialutils.get_decay_channel_tag(args)

    config = CombineConfigKTCGKB(args)
    config.configure_ktcg()
    config.datacard = workspaces.couplingdependentBRs[decay_channel]
    config.tags.append('couplingdependentBRs')
    config.range_couplingdependentBRs(args)

    if args.scenario2:
        config.datacard = workspaces.couplingdependentBRs.scenario2[decay_channel]
        config.tags.append('scenario2')

    differentialutils.run_postfit_fastscan_scan(config)
    # differentialutils.run_postfit_scan(config)


workspaces.floatingBRs.hgg                  = 'out/workspaces_Jul10/projection_ktcgkb_hgg_Jul10_couplingModel_reweighted_floatingBRs.root'
workspaces.floatingBRs.hzz                  = 'out/workspaces_Jul10/projection_ktcgkb_hzz_Jul10_couplingModel_reweighted_floatingBRs.root'
workspaces.floatingBRs.combWithHbb          = 'out/workspaces_Jul10/projection_ktcgkb_combWithHbb_Jul10_couplingModel_reweighted_constrainedbbZZ_floatingBRs.root'

workspaces.floatingBRs.scenario2.hgg         = 'out/workspaces_Jul18/projection_ktcgkb_hgg_s2groups_Jul18_couplingModel_reweighted_scenario2_floatingBRs.root'
workspaces.floatingBRs.scenario2.hzz         = 'out/workspaces_Jul18/projection_ktcgkb_hzz_s2groups_Jul18_couplingModel_reweighted_scenario2_floatingBRs.root'
workspaces.floatingBRs.scenario2.combWithHbb = 'out/workspaces_Jul18/projection_ktcgkb_combWithHbb_s2groups_Jul18_couplingModel_reweighted_constrainedbbZZ_scenario2_floatingBRs.root'

@flag_as_option
def projection_ktcg_scan_floatingBRs(args):
    args = differentialutils.force_asimov(args)
    decay_channel = differentialutils.get_decay_channel_tag(args)

    config = CombineConfigKTCGKB(args)
    config.configure_ktcg()
    config.datacard = workspaces.floatingBRs[decay_channel]
    config.tags.append('floatingBRs')
    config.range_floatingBRs(args)

    if args.scenario2:
        config.datacard = workspaces.floatingBRs.scenario2[decay_channel]
        config.tags.append('scenario2')

    differentialutils.run_postfit_scan(config)

