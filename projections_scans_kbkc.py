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
workspaces.floatingBRs = differentials.core.AttrDict()

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
        self.deltaNLLCutOff = 70.

        self.lumiscale = 83.56546
        self.hardPhysicsModelParameters.append('lumiscale={0}'.format(self.lumiscale))


    def range_couplingdependentBRs(self, args):
        if args.hzz:
            self.range_hzz_couplingdependentBRs()
        elif args.hgg:
            self.range_hgg_couplingdependentBRs()
        elif args.combination:
            self.range_combination_couplingdependentBRs()

    def range_hzz_couplingdependentBRs(self):
        self.nPoints = 55*55
        self.nPointsPerJob = 200
        self.queue = 'short.q'
        self.set_parameter_range('kappab', -1.3, 1.3)
        self.set_parameter_range('kappac', -6., 6.)

    def range_hgg_couplingdependentBRs(self):
        self.nPoints       = 50*50
        self.nPointsPerJob = 25
        self.queue = 'all.q'
        self.set_parameter_range('kappab', -1.3, 1.3)
        self.set_parameter_range('kappac', -6., 6.)

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
        self.nPoints = 30*30
        self.nPointsPerJob = 100
        self.set_parameter_range('kappab', -5., 8.)
        self.set_parameter_range('kappac', -16., 18.)

    def range_hgg_floatingBRs(self):
        self.range_hzz_floatingBRs()
        self.nPointsPerJob = 5

    def range_combination_floatingBRs(self):
        # For asimov, around 15 min per point for the combination
        # with only little variance iirc
        self.range_hzz_floatingBRs()
        self.nPointsPerJob = 5


workspaces.couplingdependentBRs.hgg         = 'out/workspaces_Jul10/projection_yukawa_hgg_Jul09_couplingModel_reweighted_couplingdependentBRs.root'
workspaces.couplingdependentBRs.hzz         = 'out/workspaces_Jul09/projection_yukawa_hzz_Jul09_couplingModel_reweighted_couplingdependentBRs.root'
workspaces.couplingdependentBRs.combination = 'out/workspaces_Jul10/projection_yukawa_combination_Jul09_couplingModel_reweighted_couplingdependentBRs.root'

@flag_as_option
def projection_kbkc_scan_couplingdependentBRs(args):
    args = differentialutils.force_asimov(args)
    decay_channel = differentialutils.get_decay_channel_tag(args)

    config = CombineConfigKBKC(args)
    config.datacard = workspaces.couplingdependentBRs[decay_channel]
    config.tags.append('couplingdependentBRs')
    config.range_couplingdependentBRs(args)

    differentialutils.run_postfit_fastscan_scan(config)
    # differentialutils.run_postfit_scan(config)


workspaces.floatingBRs.hgg                  = 'out/workspaces_Jul10/projection_yukawa_hgg_Jul09_couplingModel_reweighted_floatingBRs.root'
workspaces.floatingBRs.hzz                  = 'out/workspaces_Jul10/projection_yukawa_hzz_Jul09_couplingModel_reweighted_floatingBRs.root'
workspaces.floatingBRs.combination          = 'out/workspaces_Jul10/projection_yukawa_combination_Jul09_couplingModel_reweighted_floatingBRs.root'

@flag_as_option
def projection_kbkc_scan_floatingBRs(args):
    args = differentialutils.force_asimov(args)
    decay_channel = differentialutils.get_decay_channel_tag(args)
    config = CombineConfigKBKC(args)
    config.datacard = workspaces.floatingBRs[decay_channel]
    config.tags.append('floatingBRs')
    config.range_floatingBRs(args)
    differentialutils.run_postfit_scan(config)

