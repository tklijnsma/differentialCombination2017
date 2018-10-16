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
import ROOT

ROOT.gStyle.SetEndErrorSize(3)
ROOT.gStyle.SetHatchesLineWidth(2)
style = differentials.plotting.pywrappers.StyleSheet()
style.line_width = 2
style.error_bar_line_width = 1

from projections_plots_kbkc import PlotPatcher

#____________________________________________________________________


couplingdependentBRs = differentials.core.AttrDict()
couplingdependentBRs.hzz = 'out/Scan_projection_ktcg_Oct04_hzz_couplingdependentBRs_asimov'
couplingdependentBRs.hgg = [
    'out/Scan_projection_ktcg_Oct05_hgg_couplingdependentBRs_asimov',
    'out/Scan_projection_ktcg_Oct05_hgg_couplingdependentBRs_asimov_rescan_Oct08'
    ]
couplingdependentBRs.combWithHbb = [
    # 'out/Scan_projection_ktcg_Oct04_combWithHbb_couplingdependentBRs_asimov',
    # Need to rerun... but probably with more clever settings
    # 'out/Scan_projection_ktcg_Oct12_combWithHbb_couplingdependentBRs_asimov'
    # 'out/Scan_projection_ktcg_Oct13_combWithHbb_couplingdependentBRs_asimov'
    # 'out/Scan_projection_ktcg_Oct13_combWithHbb_couplingdependentBRs_asimov_0'
    'out/Scan_projection_ktcg_Oct14_combWithHbb_couplingdependentBRs_asimov'
    ]

couplingdependentBRs.scenario2 = differentials.core.AttrDict()
couplingdependentBRs.scenario2.hzz = 'out/Scan_projection_ktcg_Oct04_hzz_scenario2_couplingdependentBRs_asimov_0'
couplingdependentBRs.scenario2.hgg = 'out/Scan_projection_ktcg_Aug09_hgg_couplingdependentBRs_scenario2_asimov'
couplingdependentBRs.scenario2.combWithHbb = [
    'out/Scan_projection_ktcg_Oct04_combWithHbb_scenario2_couplingdependentBRs_asimov',
    'out/Scan_projection_ktcg_Oct04_combWithHbb_scenario2_couplingdependentBRs_asimov_rescan_Oct06',
    # Rescan not enough
    # Need to rerun... but probably with more clever settings
    ]

# Waiting to be finished:
# out/Scan_projection_ktcg_Oct04_combWithHbb_couplingdependentBRs_asimov
# out/Scan_projection_ktcg_Oct05_hgg_couplingdependentBRs_asimov
# out/Scan_projection_ktcg_Oct04_combWithHbb_scenario2_couplingdependentBRs_asimov_rescan_Oct06

# couplingdependentBRs.combWithHbb probably will need a rescan (scen 2 needed it)
# out/Scan_projection_ktcg_Oct04_combWithHbb_couplingdependentBRs_asimov

# Regarding combination:
# From bkglogs/bkg_2018-10-04-145429.log (scenario 1), bestfit:
# Done in 214.55 min (cpu), 214.80 min (real)
# From bkglogs/bkg_2018-10-04-145434.log (scenario 2), bestfit:
# Done in 171.39 min (cpu), 171.57 min (real)
# So, pretty slow... 2 pnts/job on all.q is not a bad choice
# Too high failure, need 1 pnt/job


class PlotPatcherKTCG(PlotPatcher):
    x_sm=1.0
    y_sm=0.0
    x_var = 'ct'
    y_var = 'cg'

    def patch_hzz(self):
        self.shift_scan()
        return self.scan.to_hist()

    def patch_hgg(self):
        self.shift_scan()
        return self.scan.to_hist()

    def patch_combination(self):
        self.shift_scan()
        return self.scan.to_hist()


class PlotPatcherKTCG_couplingdependentBRs(PlotPatcherKTCG):
    """docstring for PlotPatcherKTCG_couplingdependentBRs"""
    def __init__(self, args, decay_channel):
        super(PlotPatcherKTCG_couplingdependentBRs, self).__init__(args, decay_channel)

    def get_scandict(self):
        return couplingdependentBRs.scenario2 if self.is_scenario2 else couplingdependentBRs

    def patch_combination(self):
        self.shift_scan()
        H = self.scan.to_hist()
        H.add_padding(
            700.,
            x_min = 0.72,
            x_max = 1.44,
            y_min = -0.035,
            y_max = 0.035
            )
        return H


class PlotPatcherKTCG_couplingdependentBRs_scenario2(PlotPatcherKTCG_couplingdependentBRs):
    """docstring for PlotPatcherKTCG_couplingdependentBRs_scenario2"""

    def patch_combination(self):
        self.shift_scan()
        H = self.scan.to_hist()
        H.add_padding(
            700.,
            x_min = 0.7,
            x_max = 1.4,
            y_min = -0.02,
            y_max = 0.025
            )
        return H


@flag_as_option
def projection_ktcg_plot_comparison(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    Patcher = PlotPatcherKTCG_couplingdependentBRs_scenario2 if args.scenario2 else PlotPatcherKTCG_couplingdependentBRs
    Patcher(args, decay_channel).quick_patched_vs_unpatched()


@flag_as_option
def projection_ktcg_plot_couplingdependentBRs(args):
    Patcher = PlotPatcherKTCG_couplingdependentBRs_scenario2 if args.scenario2 else PlotPatcherKTCG_couplingdependentBRs
    scans = []

    combWithHbb = Patcher(args, 'combWithHbb').patch()
    scans.append(combWithHbb)

    hgg = Patcher(args, 'hgg').patch()
    scans.append(hgg)

    hzz = Patcher(args, 'hzz').patch()
    scans.append(hzz)

    x_min, x_max, y_min, y_max = (None, None, None, None)
    if args.scenario2:
        x_min = 0.7
        x_max = 1.4
        y_min = -0.02
        y_max = 0.015

    plot = differentials.plotting.plots.MultiContourPlot(
        ('projection_ktcg_plot_couplingdependentBRs'
            + ('_asimov' if args.asimov else '')
            + ('_scenario2' if args.scenario2 else '')
            ),
        scans,
        x_title = differentials.core.standard_titles['ct'], y_title = differentials.core.standard_titles['cg'],
        x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max
        )
    plot.draw()


#____________________________________________________________________
floatingBRs = differentials.core.AttrDict()
# floatingBRs.hzz = 'out/Scan_projection_ktcg_Jul10_hzz_floatingBRs_asimov'

floatingBRs.scenario2 = differentials.core.AttrDict()

floatingBRs.hzz = 'out/Scan_projection_ktcg_Oct05_hzz_floatingBRs_asimov_0'
floatingBRs.hgg = 'out/Scan_projection_ktcg_Oct04_hgg_floatingBRs_asimov'
floatingBRs.combWithHbb = ''

# Slightly improved stability w.r.t. to the old min settings
# floatingBRs.hzz = 'out/Scan_projection_ktcg_Oct04_hzz_floatingBRs_asimov'

floatingBRs.scenario2.hzz = 'out/Scan_projection_ktcg_Oct04_hzz_scenario2_floatingBRs_asimov'
floatingBRs.scenario2.hgg = 'out/Scan_projection_ktcg_Oct04_hgg_scenario2_floatingBRs_asimov'
floatingBRs.scenario2.combWithHbb = ''

# Waiting to be finished:
# out/Scan_projection_ktcg_Oct06_combWithHbb_floatingBRs_asimov
# out/Scan_projection_ktcg_Oct06_combWithHbb_scenario2_floatingBRs_asimov

# Running:
# out/Scan_projection_ktcg_Oct06_combWithHbb_floatingBRs_asimov
# out/Scan_projection_ktcg_Oct06_combWithHbb_scenario2_floatingBRs_asimov

class PlotPatcherKTCG_floatingBRs(PlotPatcherKTCG):
    """docstring for PlotPatcherKTCG_floatingBRs"""
    def __init__(self, args, decay_channel):
        super(PlotPatcherKTCG_floatingBRs, self).__init__(args, decay_channel)

    def get_scandict(self):
        return floatingBRs.scenario2 if self.args.scenario2 else floatingBRs


class PlotPatcherKTCG_floatingBRs_scenario2(PlotPatcherKTCG_floatingBRs):
    """docstring for PlotPatcherKTCG_floatingBRs_scenario2"""
    pass        



@flag_as_option
def projection_ktcg_plot_floatingBRs(args):
    Patcher = PlotPatcherKTCG_floatingBRs_scenario2 if args.scenario2 else PlotPatcherKTCG_floatingBRs
    scans = []

    # combWithHbb = Patcher(args, 'combWithHbb').patch()
    # scans.append(combWithHbb)

    hzz = Patcher(args, 'hzz').patch()
    scans.append(hzz)

    hgg = Patcher(args, 'hgg').patch()
    scans.append(hgg)

    x_min, x_max, y_min, y_max = (None, None, None, None)
    if args.scenario2:
        x_min = 0.7
        x_max = 1.4
        y_min = -0.02
        y_max = 0.015

    plot = differentials.plotting.plots.MultiContourPlot(
        ('projection_ktcg_plot_floatingBRs'
            + ('_asimov' if args.asimov else '')
            + ('_scenario2' if args.scenario2 else '')
            ),
        scans,
        x_title = differentials.core.standard_titles['ct'], y_title = differentials.core.standard_titles['cg'],
        x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max
        )
    plot.draw()


