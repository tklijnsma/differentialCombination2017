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

workspaces.couplingdependentBRs.hgg         = 'out/workspaces_Jul10/projection_yukawa_hgg_Jul09_couplingModel_reweighted_couplingdependentBRs.root'
workspaces.couplingdependentBRs.hzz         = 'out/workspaces_Jul09/projection_yukawa_hzz_Jul09_couplingModel_reweighted_couplingdependentBRs.root'
workspaces.couplingdependentBRs.combination = 'out/workspaces_Jul10/projection_yukawa_combination_Jul09_couplingModel_reweighted_couplingdependentBRs.root'

workspaces.couplingdependentBRs.scenario2.hgg         = 'out/workspaces_Jul17/projection_yukawa_hgg_s2groups_Jul17_couplingModel_reweighted_scenario2_couplingdependentBRs.root'
workspaces.couplingdependentBRs.scenario2.hzz         = 'out/workspaces_Jul17/projection_yukawa_hzz_s2groups_Jul17_couplingModel_reweighted_scenario2_couplingdependentBRs.root'
workspaces.couplingdependentBRs.scenario2.combination = 'out/workspaces_Jul17/projection_yukawa_combination_s2groups_Jul17_couplingModel_reweighted_scenario2_couplingdependentBRs.root'

workspaces.floatingBRs.hgg                  = 'out/workspaces_Jul10/projection_yukawa_hgg_Jul09_couplingModel_reweighted_floatingBRs.root'
workspaces.floatingBRs.hzz                  = 'out/workspaces_Jul10/projection_yukawa_hzz_Jul09_couplingModel_reweighted_floatingBRs.root'
workspaces.floatingBRs.combination          = 'out/workspaces_Jul10/projection_yukawa_combination_Jul09_couplingModel_reweighted_floatingBRs.root'

workspaces.floatingBRs.scenario2.hgg         = 'out/workspaces_Jul17/projection_yukawa_hgg_s2groups_Jul17_couplingModel_reweighted_scenario2_floatingBRs.root'
workspaces.floatingBRs.scenario2.hzz         = 'out/workspaces_Jul17/projection_yukawa_hzz_s2groups_Jul17_couplingModel_reweighted_scenario2_floatingBRs.root'
workspaces.floatingBRs.scenario2.combination = 'out/workspaces_Jul17/projection_yukawa_combination_s2groups_Jul17_couplingModel_reweighted_scenario2_floatingBRs.root'


class CombineConfigKBKC(differentials.combine.combine.CombineConfig):
    """docstring for CombineConfigKBKC"""

    lumiscale = 83.56546

    def __init__(self, args):
        super(CombineConfigKBKC, self).__init__(args)
        args = differentialutils.force_asimov(args)
        self.onBatch       = True
        self.queue         = 'all.q'
        self.asimov = True if args.asimov else False
        self.decay_channel = differentialutils.get_decay_channel_tag(args)
        self.POIs = [ 'kappab', 'kappac' ]
        self.PhysicsModelParameters = [ 'kappab=1.0', 'kappac=1.0' ]
        self.subDirectory = 'out/Scan_projection_kbkc_{0}_{1}'.format(differentials.core.datestr(), self.decay_channel)
        self.hardPhysicsModelParameters.append('lumiscale={0}'.format(self.lumiscale))
        self.freezeNuisances.append('lumiscale')
        if args.scenario2: self.tags.append('scenario2')
        self.set_ranges()
        self.datacard = self.get_workspacedict()[self.decay_channel]

    def is_combination(self):
        return not(self.decay_channel in ['hgg', 'hzz', 'hbb'])

    def set_ranges(self):
        if self.is_combination():
            self.set_ranges_combination()
        else:
            getattr(self, 'set_ranges_' + self.decay_channel)()

    def set_ranges_combination(self):
        pass

    def set_ranges_hgg(self):
        pass

    def set_ranges_hzz(self):
        pass


class CombineConfigKBKC_couplingdependentBRs(CombineConfigKBKC):
    """docstring for CombineConfigKBKC_couplingdependentBRs"""
    def __init__(self, args):
        super(CombineConfigKBKC_couplingdependentBRs, self).__init__(args)
        self.deltaNLLCutOff = 1000.
        self.tags.append('couplingdependentBRs')

    def get_workspacedict(self):
        return workspaces.couplingdependentBRs.scenario2 if self.args.scenario2 else workspaces.couplingdependentBRs

    def set_ranges_combination(self):
        self.nPointsPerJob = 22
        self.queue = 'all.q'
        self.nPoints = 150*150
        self.set_parameter_range('kappab', -1.3, 1.3)
        self.set_parameter_range('kappac', -6., 6.)

    def set_ranges_hgg(self):
        self.nPointsPerJob = 48
        self.queue = 'all.q'
        self.nPoints = 100*100
        self.set_parameter_range('kappab', -3.50, 3.50)
        self.set_parameter_range('kappac', -8.50, 8.50)

    def set_ranges_hzz(self):
        self.nPointsPerJob = 500
        self.queue = 'short.q'
        self.nPoints = 100*100
        self.set_parameter_range('kappab', -1.3, 1.3)
        self.set_parameter_range('kappac', -6., 6.)


class CombineConfigKBKC_couplingdependentBRs_scenario2(CombineConfigKBKC_couplingdependentBRs):
    """docstring for CombineConfigKBKC_couplingdependentBRs_scenario2"""
    def __init__(self, args):
        super(CombineConfigKBKC_couplingdependentBRs_scenario2, self).__init__(args)

    def set_ranges_combination(self):
        self.nPointsPerJob = 79
        self.queue = 'all.q'
        self.nPoints = 150*150
        self.set_parameter_range('kappab', -1.15, 1.15)
        self.set_parameter_range('kappac', -5.25, 5.25)
        
    def set_ranges_hgg(self):
        self.nPointsPerJob = 147
        self.queue = 'all.q'
        self.nPoints = 100*100
        self.set_parameter_range('kappab', -1.3, 1.3)
        self.set_parameter_range('kappac', -5.25, 5.25)

    def set_ranges_hzz(self):
        self.nPointsPerJob = 500
        self.queue = 'short.q'
        self.nPoints = 100*100
        self.set_parameter_range('kappab', -1.15, 1.15)
        self.set_parameter_range('kappac', -5.25, 5.25)


class CombineConfigKBKC_floatingBRs(CombineConfigKBKC):
    """docstring for CombineConfigKBKC_couplingdependentBRs"""
    def __init__(self, args):
        super(CombineConfigKBKC_floatingBRs, self).__init__(args)
        self.tags.append('floatingBRs')

    def get_workspacedict(self):
        return workspaces.floatingBRs.scenario2 if self.args.scenario2 else workspaces.floatingBRs

    def set_ranges_hzz(self):
        self.queue = 'short.q'
        self.nPoints = 80*80
        self.nPointsPerJob = 500
        self.set_parameter_range('kappab', -3.5, 6.)
        self.set_parameter_range('kappac', -13., 15.)

    def set_ranges_hgg(self):
        self.queue = 'all.q'
        self.nPoints = 45*45
        self.nPointsPerJob = 70
        self.set_parameter_range('kappab', -2.25, 4.75)
        self.set_parameter_range('kappac', -10., 12.5)

    def set_ranges_combination(self):
        self.queue = 'all.q'
        self.nPoints = 45*45
        self.nPointsPerJob = 40
        self.set_parameter_range('kappab', -2.25, 4.75)
        self.set_parameter_range('kappac', -10., 12.5)


class CombineConfigKBKC_floatingBRs_scenario2(CombineConfigKBKC_floatingBRs):
    """docstring for CombineConfigKBKC_floatingBRs_scenario2"""
    def __init__(self, args):
        super(CombineConfigKBKC_floatingBRs_scenario2, self).__init__(args)
        self.tags.append('scenario2')

    def set_ranges_hzz(self):
        super(CombineConfigKBKC_floatingBRs_scenario2, self).set_ranges_hzz()

    def set_ranges_hgg(self):
        super(CombineConfigKBKC_floatingBRs_scenario2, self).set_ranges_hgg()
        self.nPointsPerJob = 81

    def set_ranges_combination(self):
        super(CombineConfigKBKC_floatingBRs_scenario2, self).set_ranges_combination()
        self.nPointsPerJob = 53


@flag_as_option
def projection_kbkc_scan_couplingdependentBRs(args):
    args = differentialutils.force_asimov(args)
    Config = CombineConfigKBKC_couplingdependentBRs_scenario2 if args.scenario2 else CombineConfigKBKC_couplingdependentBRs
    config = Config(args)
    differentialutils.run_postfit_fastscan_scan(config)

@flag_as_option
def projection_kbkc_scan_floatingBRs(args):
    args = differentialutils.force_asimov(args)
    Config = CombineConfigKBKC_floatingBRs_scenario2 if args.scenario2 else CombineConfigKBKC_floatingBRs
    config = Config(args)
    differentialutils.run_postfit_scan(config)




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
