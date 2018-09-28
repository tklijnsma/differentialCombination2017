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

class CombineConfigKBKC(differentials.combine.combine.CombineConfig):
    """docstring for CombineConfigKBKC"""
    def __init__(self, args):
        super(CombineConfigKBKC, self).__init__(args)
        self.onBatch       = True
        self.queue         = 'all.q'
        self.asimov = True if args.asimov else False
        self.decay_channel = differentialutils.get_decay_channel_tag(args)

        self.POIs = [ 'kappab', 'kappac' ]
        self.PhysicsModelParameters = [ 'kappab=1.0', 'kappac=1.0' ]
        self.subDirectory = 'out/Scan_projection_kbkc_{0}_{1}'.format(differentials.core.datestr(), self.decay_channel)
        self.deltaNLLCutOff = 1000.

        self.lumiscale = 83.56546
        self.hardPhysicsModelParameters.append('lumiscale={0}'.format(self.lumiscale))


    def range_couplingdependentBRs(self, args):
        self.nPoints       = 100*100
        if args.scenario2:
            self.set_parameter_range('kappab', -1.15, 1.15)
            self.set_parameter_range('kappac', -5.25, 5.25)
        else:
            self.set_parameter_range('kappab', -1.3, 1.3)
            self.set_parameter_range('kappac', -6., 6.)

        if args.hzz:
            self.range_hzz_couplingdependentBRs()
        elif args.hgg:
            self.range_hgg_couplingdependentBRs()
        elif args.combination:
            self.range_combination_couplingdependentBRs()

    def range_hzz_couplingdependentBRs(self):
        self.nPointsPerJob = 200
        self.queue = 'short.q'

    def range_hgg_couplingdependentBRs(self):
        self.nPointsPerJob = 20
        self.queue = 'all.q'

    def range_combination_couplingdependentBRs(self):
        self.range_hgg_couplingdependentBRs()


    def range_floatingBRs(self, args):
        if args.hzz:
            self.range_hzz_floatingBRs()
        elif args.hgg:
            self.range_hgg_floatingBRs()
        elif args.combination:
            self.range_combination_floatingBRs()

    def range_hzz_floatingBRs(self):
        self.nPoints = 35*35
        self.nPointsPerJob = 100
        self.set_parameter_range('kappab', -5., 8.)
        self.set_parameter_range('kappac', -16., 18.)

    def range_hgg_floatingBRs(self):
        self.nPoints = 35*35
        self.nPointsPerJob = 5
        self.set_parameter_range('kappab', -2., 3.5)
        self.set_parameter_range('kappac', -8., 10.)

    def range_combination_floatingBRs(self):
        # For asimov, around 15 min per point for the combination
        # with only little variance iirc
        self.range_hgg_floatingBRs()


workspaces.couplingdependentBRs.hgg         = 'out/workspaces_Jul10/projection_yukawa_hgg_Jul09_couplingModel_reweighted_couplingdependentBRs.root'
workspaces.couplingdependentBRs.hzz         = 'out/workspaces_Jul09/projection_yukawa_hzz_Jul09_couplingModel_reweighted_couplingdependentBRs.root'
workspaces.couplingdependentBRs.combination = 'out/workspaces_Jul10/projection_yukawa_combination_Jul09_couplingModel_reweighted_couplingdependentBRs.root'

workspaces.couplingdependentBRs.scenario2.hgg         = 'out/workspaces_Jul17/projection_yukawa_hgg_s2groups_Jul17_couplingModel_reweighted_scenario2_couplingdependentBRs.root'
workspaces.couplingdependentBRs.scenario2.hzz         = 'out/workspaces_Jul17/projection_yukawa_hzz_s2groups_Jul17_couplingModel_reweighted_scenario2_couplingdependentBRs.root'
workspaces.couplingdependentBRs.scenario2.combination = 'out/workspaces_Jul17/projection_yukawa_combination_s2groups_Jul17_couplingModel_reweighted_scenario2_couplingdependentBRs.root'

@flag_as_option
def projection_kbkc_scan_couplingdependentBRs(args):
    args = differentialutils.force_asimov(args)
    decay_channel = differentialutils.get_decay_channel_tag(args)

    config = CombineConfigKBKC(args)
    config.datacard = workspaces.couplingdependentBRs[decay_channel]
    config.tags.append('couplingdependentBRs')
    config.range_couplingdependentBRs(args)

    if args.scenario2:
        config.datacard = workspaces.couplingdependentBRs.scenario2[decay_channel]
        config.tags.append('scenario2')

    # config.tags.append('testscan')
    # config.nPoints = 5*5
    # config.nPointsPerJob = 25
    differentialutils.run_postfit_fastscan_scan(config)
    # differentialutils.run_postfit_scan(config)

@flag_as_option
def test_dopoints_option_implementation(args):
    args = differentialutils.set_one_decay_channel(args, 'hzz', asimov=True)
    decay_channel = differentialutils.get_decay_channel_tag(args)
    config = CombineConfigKBKC(args)
    config.datacard = workspaces.couplingdependentBRs[decay_channel]
    config.tags.append('couplingdependentBRs')
    config.range_couplingdependentBRs(args)
    config.nPoints = 10*10
    config.nPointsPerJob = 50
    differentialutils.run_postfit_fastscan_scan(config)



workspaces.floatingBRs.hgg                  = 'out/workspaces_Jul10/projection_yukawa_hgg_Jul09_couplingModel_reweighted_floatingBRs.root'
workspaces.floatingBRs.hzz                  = 'out/workspaces_Jul10/projection_yukawa_hzz_Jul09_couplingModel_reweighted_floatingBRs.root'
workspaces.floatingBRs.combination          = 'out/workspaces_Jul10/projection_yukawa_combination_Jul09_couplingModel_reweighted_floatingBRs.root'

workspaces.floatingBRs.scenario2.hgg         = 'out/workspaces_Jul17/projection_yukawa_hgg_s2groups_Jul17_couplingModel_reweighted_scenario2_floatingBRs.root'
workspaces.floatingBRs.scenario2.hzz         = 'out/workspaces_Jul17/projection_yukawa_hzz_s2groups_Jul17_couplingModel_reweighted_scenario2_floatingBRs.root'
workspaces.floatingBRs.scenario2.combination = 'out/workspaces_Jul17/projection_yukawa_combination_s2groups_Jul17_couplingModel_reweighted_scenario2_floatingBRs.root'

@flag_as_option
def projection_kbkc_scan_floatingBRs(args):
    args = differentialutils.force_asimov(args)
    decay_channel = differentialutils.get_decay_channel_tag(args)
    config = CombineConfigKBKC(args)
    config.datacard = workspaces.floatingBRs[decay_channel]
    config.tags.append('floatingBRs')
    config.range_floatingBRs(args)
    if args.scenario2:
        config.datacard = workspaces.floatingBRs.scenario2[decay_channel]
        config.tags.append('scenario2')
    differentialutils.run_postfit_scan(config)

