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

from projections_scans_kbkc import CombineConfigKBKC

lumiscale300 = 8.356546
lumiscale3000 = 83.56546


#____________________________________________________________________
workspaces = differentials.core.AttrDict()
workspaces.couplingdependentBRs = differentials.core.AttrDict()
workspaces.couplingdependentBRs.scenario2 = differentials.core.AttrDict()
workspaces.floatingBRs = differentials.core.AttrDict()
workspaces.floatingBRs.scenario2 = differentials.core.AttrDict()


workspaces.couplingdependentBRs.hgg         = 'out/workspaces_Jul10/projection_ktcgkb_hgg_Jul10_couplingModel_reweighted_couplingdependentBRs.root'
workspaces.couplingdependentBRs.hzz         = 'out/workspaces_Jul10/projection_ktcgkb_hzz_Jul10_couplingModel_reweighted_couplingdependentBRs.root'
workspaces.couplingdependentBRs.combWithHbb = 'out/workspaces_Jul10/projection_ktcgkb_combWithHbb_Jul10_couplingModel_reweighted_couplingdependentBRs.root'

workspaces.couplingdependentBRs.scenario2.hgg         = 'out/workspaces_Jul18/projection_ktcgkb_hgg_s2groups_Jul18_couplingModel_reweighted_scenario2_couplingdependentBRs.root'
workspaces.couplingdependentBRs.scenario2.hzz         = 'out/workspaces_Jul18/projection_ktcgkb_hzz_s2groups_Jul18_couplingModel_reweighted_scenario2_couplingdependentBRs.root'
workspaces.couplingdependentBRs.scenario2.combWithHbb = 'out/workspaces_Jul18/projection_ktcgkb_combWithHbb_s2groups_Jul18_couplingModel_reweighted_scenario2_couplingdependentBRs.root'


class CombineConfigKTCGKB(CombineConfigKBKC):
    """docstring for CombineConfigKTCGKB"""

    def __init__(self, args):
        super(CombineConfigKTCGKB, self).__init__(args)

    def minimizer_settings_andrew(self):
        logging.warning('Applying minimizer_settings suggested by Andrew')
        self.minimizer_settings = [
            # '--cminDefaultMinimizerStrategy 0',
            # '--cminDefaultMinimizerTolerance 0.1',
            '--minimizerStrategy 0',
            '--minimizerTolerance 0.1',
            '--cminFallbackAlgo "Minuit2,migrad,0:0.11"',
            '--cminFallbackAlgo "Minuit2,migrad,0:0.12"',
            '--cminApproxPreFitTolerance=100',
            '--X-rtd MINIMIZER_MaxCalls=9999999',
            ]

    def minimizer_settings_overkill(self):
        logging.warning('Applying minimizer_settings that are overkill')
        self.minimizer_settings =[
            '--minimizerStrategy 2',
            '--minimizerTolerance 0.001',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            ]

    def minimizer_settings_reasonable(self):
        logging.warning('Applying minimizer_settings that are reasonable')
        self.minimizer_settings =[
            '--minimizerStrategy 0',
            '--minimizerTolerance 0.01',
            '--robustFit 1',
            '--minimizerAlgoForMinos Minuit2,Migrad',
            '--X-rtd MINIMIZER_MaxCalls=9999999',
            ]

    def freeze_mcstat_uncertainties(self):
        self.freeze_parameters(
            [ 'hqq125GenpT2failcat1mcstat', 'hqq125GenpT2passcat1mcstat', 'hqq125GenpT3failcat1mcstat', 'hqq125GenpT3failcat2mcstat', 'hqq125GenpT3passcat1mcstat', 'hqq125GenpT3passcat2mcstat', 'hqq125GenpT4failcat1mcstat', 'hqq125GenpT4failcat2mcstat', 'hqq125GenpT4passcat1mcstat', 'hqq125GenpT4passcat2mcstat', 'hqq125failmuonCRmcstat', 'qcdfailmuonCRmcstat', 'qcdpassmuonCRmcstat', 'stqqfailmuonCRmcstat', 'stqqpassmuonCRmcstat', 'tqqfailcat1mcstat', 'tqqfailcat2mcstat', 'tqqfailmuonCRmcstat', 'tqqpasscat1mcstat', 'tqqpasscat2mcstat', 'tqqpassmuonCRmcstat', 'vvqqfailmuonCRmcstat', 'wlnufailmuonCRmcstat', 'wlnupassmuonCRmcstat', 'wqqfailcat1mcstat', 'wqqfailcat2mcstat', 'wqqpasscat1mcstat', 'wqqpasscat2mcstat', 'xhqq125GenpT1failcat1mcstat', 'xhqq125GenpT1failcat2mcstat', 'xhqq125GenpT1passcat1mcstat', 'xhqq125GenpT1passcat2mcstat', 'xhqq125GenpT2failcat1mcstat', 'xhqq125GenpT2failcat2mcstat', 'xhqq125GenpT2passcat1mcstat', 'xhqq125GenpT2passcat2mcstat', 'xhqq125GenpT3failcat1mcstat', 'xhqq125GenpT3failcat2mcstat', 'xhqq125GenpT3passcat1mcstat', 'xhqq125GenpT3passcat2mcstat', 'xhqq125GenpT4failcat1mcstat', 'xhqq125GenpT4failcat2mcstat', 'xhqq125GenpT4passcat1mcstat', 'xhqq125GenpT4passcat2mcstat', 'xhqq125failmuonCRmcstat', 'xhqq125passmuonCRmcstat', 'zllfailmuonCRmcstat', 'zqqfailcat1mcstat', 'zqqfailcat2mcstat', 'zqqfailmuonCRmcstat', 'zqqpasscat1mcstat', 'zqqpasscat2mcstat' ]
            )

    def configure_ktcg(self):
        self.POIs = [ 'ct', 'cg' ]
        self.PhysicsModelParameters = [ 'ct=1.0', 'cg=0.0' ]
        self.subDirectory = 'out/Scan_projection_ktcg_{0}_{1}'.format(differentials.core.datestr(), self.decay_channel)

    def configure_ktkb(self):
        self.POIs = [ 'ct', 'cb' ]
        self.PhysicsModelParameters = [ 'ct=1.0', 'cb=1.0' ]
        self.subDirectory = 'out/Scan_projection_ktkb_{0}_{1}'.format(differentials.core.datestr(), self.decay_channel)


class CombineConfigKTCG_couplingdependentBRs(CombineConfigKTCGKB):
    """docstring for CombineConfigKTCG_couplingdependentBRs"""
    def __init__(self, args):
        super(CombineConfigKTCG_couplingdependentBRs, self).__init__(args)
        self.deltaNLLCutOff = 1000.
        self.tags.append('couplingdependentBRs')
        self.configure_ktcg()

    def get_workspacedict(self):
        return workspaces.couplingdependentBRs.scenario2 if self.args.scenario2 else workspaces.couplingdependentBRs

    def set_ranges_hzz(self):
        self.nPointsPerJob = 250
        self.queue = 'short.q'
        self.nPoints = 100*100
        self.set_parameter_range('ct', -0.1, 3.0)
        self.set_parameter_range('cg', -0.2, 0.1)

    def set_ranges_hgg(self):
        self.nPointsPerJob = 32
        self.queue = 'all.q'
        self.nPoints = 50*50
        self.set_parameter_range('ct', 0.7, 1.45)
        self.set_parameter_range('cg', -0.035, 0.03)

    def set_ranges_combination(self):
        self.set_parameter_range('ct', 0.7, 1.35)
        self.set_parameter_range('cg', -0.035, 0.03)

        # self.nPointsPerJob = 2
        # self.queue = 'all.q'
        # self.nPoints = 50*50

        self.nPointsPerJob = 8
        self.queue = 'all.q'
        self.nPoints = 70*70
        # self.minimizer_settings_andrew()


class CombineConfigKTCG_couplingdependentBRs_scenario2(CombineConfigKTCG_couplingdependentBRs):
    """docstring for CombineConfigKTCG_couplingdependentBRs_scenario2"""
    def __init__(self, args):
        super(CombineConfigKTCG_couplingdependentBRs_scenario2, self).__init__(args)

    def set_ranges_hzz(self):
        self.nPointsPerJob = 200
        self.queue = 'short.q'
        self.nPoints = 100*100
        self.set_parameter_range('ct', -0.1, 3.2)
        self.set_parameter_range('cg', -0.2, 0.08)

    def set_ranges_hgg(self):
        self.queue = 'all.q'
        self.nPoints = 60*60
        self.nPointsPerJob = 10
        self.set_parameter_range('ct', 0.7, 1.4)
        self.set_parameter_range('cg', -0.02, 0.015)
        print 'THIS SCAN IS FINE, don\'t redo it'
        sys.exit()

    def set_ranges_combination(self):
        # self.queue = 'all.q'
        # self.nPoints = 40*40
        # self.nPointsPerJob = 2

        self.nPointsPerJob = 8
        self.queue = 'all.q'
        self.nPoints = 70*70

        self.set_parameter_range('ct', 0.8, 1.18)
        self.set_parameter_range('cg', -0.012, 0.012)
        # self.minimizer_settings_andrew()


@flag_as_option
def projection_ktcg_scan_couplingdependentBRs(args):
    Config = CombineConfigKTCG_couplingdependentBRs_scenario2 if args.scenario2 else CombineConfigKTCG_couplingdependentBRs
    config = Config(args)
    if args.combWithHbb:
        config.minimizer_settings_reasonable()
        config.freeze_mcstat_uncertainties()

        # if args.scenario2:
        #     reusable_postfit = 'out/Scan_projection_ktcg_Oct04_combWithHbb_scenario2_couplingdependentBRs_asimov/postfit_and_fastscan/higgsCombine_POSTFIT_ASIMOV_projection_ktcgkb_combWithHbb_s2groups_Jul18_couplingModel_reweighted_scenario2_couplingdependentBRs.MultiDimFit.mH125.root'
        # else:
        #     reusable_postfit = 'out/Scan_projection_ktcg_Oct04_combWithHbb_couplingdependentBRs_asimov/postfit_and_fastscan/higgsCombine_POSTFIT_ASIMOV_projection_ktcgkb_combWithHbb_Jul10_couplingModel_reweighted_couplingdependentBRs.MultiDimFit.mH125.root'
        # differentialutils.run_fastscan_scan_reused_postfit(
        #     config,
        #     reusable_postfit
        #     )

        differentialutils.run_postfit_fastscan_scan(config)

    else:
        differentialutils.run_postfit_fastscan_scan(config)


#____________________________________________________________________

workspaces.floatingBRs.hgg                  = 'out/workspaces_Jul10/projection_ktcgkb_hgg_Jul10_couplingModel_reweighted_floatingBRs.root'
workspaces.floatingBRs.hzz                  = 'out/workspaces_Jul10/projection_ktcgkb_hzz_Jul10_couplingModel_reweighted_floatingBRs.root'
workspaces.floatingBRs.combWithHbb          = 'out/workspaces_Jul10/projection_ktcgkb_combWithHbb_Jul10_couplingModel_reweighted_constrainedbbZZ_floatingBRs.root'

workspaces.floatingBRs.scenario2.hgg         = 'out/workspaces_Jul18/projection_ktcgkb_hgg_s2groups_Jul18_couplingModel_reweighted_scenario2_floatingBRs.root'
workspaces.floatingBRs.scenario2.hzz         = 'out/workspaces_Jul18/projection_ktcgkb_hzz_s2groups_Jul18_couplingModel_reweighted_scenario2_floatingBRs.root'
workspaces.floatingBRs.scenario2.combWithHbb = 'out/workspaces_Jul18/projection_ktcgkb_combWithHbb_s2groups_Jul18_couplingModel_reweighted_constrainedbbZZ_scenario2_floatingBRs.root'


class CombineConfigKTCG_floatingBRs(CombineConfigKTCGKB):
    """docstring for CombineConfigKTCG_floatingBRs"""
    def __init__(self, args):
        super(CombineConfigKTCG_floatingBRs, self).__init__(args)
        self.tags.append('floatingBRs')
        self.configure_ktcg()

    def get_workspacedict(self):
        return workspaces.floatingBRs.scenario2 if self.args.scenario2 else workspaces.floatingBRs

    def set_ranges_hzz(self):
        self.nPointsPerJob = 150
        self.queue = 'short.q'
        self.nPoints = 100*100
        self.set_parameter_range('ct', -6.0, 6.0)
        self.set_parameter_range('cg', -0.45, 0.45)
        self.minimizer_settings_andrew()

    # These ranges not very well tested; no clue for combination at all
    def set_ranges_hgg_and_combination(self):
        self.set_parameter_range('ct', -2.0, 2.0)
        self.set_parameter_range('cg', -0.08, 0.08)
        self.minimizer_settings_andrew()

    def set_ranges_hgg(self):
        self.set_ranges_hgg_and_combination()
        self.nPoints = 40*40
        self.queue = 'all.q'
        self.nPointsPerJob = 8

    def set_ranges_combination(self):
        self.set_ranges_hgg_and_combination()
        self.nPoints = 70*70
        self.queue = 'all.q'
        self.nPointsPerJob = 1

class CombineConfigKTCG_floatingBRs_scenario2(CombineConfigKTCG_floatingBRs):
    """docstring for CombineConfigKTCG_floatingBRs_scenario2"""
    def __init__(self, args):
        super(CombineConfigKTCG_floatingBRs_scenario2, self).__init__(args)

    def set_ranges_hzz(self):
        super(CombineConfigKTCG_floatingBRs_scenario2, self).set_ranges_hzz()
        self.nPointsPerJob = 200


@flag_as_option
def projection_ktcg_scan_floatingBRs(args):
    Config = CombineConfigKTCG_floatingBRs_scenario2 if args.scenario2 else CombineConfigKTCG_floatingBRs
    config = Config(args)
    differentialutils.run_postfit_scan(config)

