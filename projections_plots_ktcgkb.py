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

#____________________________________________________________________


couplingdependentBRs = differentials.core.AttrDict()
couplingdependentBRs.hzz = 'out/Scan_projection_ktcg_Jul10_hzz_couplingdependentBRs_asimov_0'
# couplingdependentBRs.hgg = 
# couplingdependentBRs.combWithHbb = 
couplingdependentBRs.hgg = 'out/Scan_projection_ktcg_Aug09_hgg_couplingdependentBRs_asimov'

couplingdependentBRs.scenario2 = differentials.core.AttrDict()
couplingdependentBRs.scenario2.hzz = 'out/Scan_projection_ktcg_Jul18_hzz_couplingdependentBRs_scenario2_asimov_0'
# couplingdependentBRs.scenario2.hgg = 'out/Scan_projection_ktcg_Jul18_hgg_couplingdependentBRs_scenario2_asimov'
# couplingdependentBRs.scenario2.combWithHbb = 'out/Scan_projection_ktcg_Jul18_combWithHbb_couplingdependentBRs_scenario2_asimov'
couplingdependentBRs.scenario2.hgg = 'out/Scan_projection_ktcg_Aug09_hgg_couplingdependentBRs_scenario2_asimov'


def latest_couplingdependentBRs(args, decay_channel=None):
    args = differentialutils.force_asimov(args)
    if decay_channel is None: decay_channel = differentialutils.get_decay_channel_tag(args)
    scandict = couplingdependentBRs.scenario2 if args.scenario2 else couplingdependentBRs
    scandir = scandict[decay_channel]
    scan = differentials.scans.Scan2D(
        'ktcg_{0}'.format(decay_channel), 'ct', 'cg', scandir = scandir
        )
    scan.title = differentials.core.standard_titles.get(decay_channel, decay_channel)
    scan.color = 1
    scan.deltaNLL_threshold = -10.
    scan.read()
    return scan

@flag_as_option
def projection_ktcg_plot_couplingdependentBRs(args):
    differentials.plotting.pywrappers.CMS_Latex_lumi.CMS_lumi = 3000.
    scans = []

    hzz = latest_couplingdependentBRs(args, decay_channel='hzz')
    hzz.color = differentials.core.safe_colors.blue
    scans.append(hzz)

    hgg = latest_couplingdependentBRs(args, decay_channel='hgg')
    hgg.color = differentials.core.safe_colors.red
    scans.append(hgg)

    # if args.scenario2:
    # combWithHbb = latest_couplingdependentBRs(args, decay_channel='combWithHbb')
    # combWithHbb.color = 1
    # scans.append(combWithHbb)

    plot = differentials.plotting.plots.MultiContourPlot(
        ('projection_ktcg_plot_couplingdependentBRs'
            + ('_asimov' if args.asimov else '')
            + ('_scenario2' if args.scenario2 else '')
            ),
        scans,
        x_title = differentials.core.standard_titles['ct'], y_title = differentials.core.standard_titles['cg'],
        # x_min=couplingdependentBRs_x_min, x_max=couplingdependentBRs_x_max, y_min=couplingdependentBRs_y_min, y_max=couplingdependentBRs_y_max
        )
    plot.draw()


#____________________________________________________________________
floatingBRs = differentials.core.AttrDict()
floatingBRs.hzz = 'out/Scan_projection_ktcg_Jul10_hzz_floatingBRs_asimov'

floatingBRs.scenario2 = differentials.core.AttrDict()
floatingBRs.scenario2.hzz = 'out/Scan_projection_ktcg_Jul18_hzz_floatingBRs_scenario2_asimov_0'
floatingBRs.scenario2.hgg = 'out/Scan_projection_ktcg_Jul18_hgg_floatingBRs_scenario2_asimov'
floatingBRs.scenario2.combWithHbb = 'out/Scan_projection_ktcg_Jul18_combWithHbb_floatingBRs_scenario2_asimov'


def latest_floatingBRs(args, decay_channel=None, asimov=None, splined=False):
    args = differentialutils.force_asimov(args)
    if decay_channel is None: decay_channel = differentialutils.get_decay_channel_tag(args)
    scandict = floatingBRs.scenario2 if args.scenario2 else floatingBRs
    scandir = scandict[decay_channel]

    scan = differentials.scans.Scan2D(
        'ktcg_{0}'.format(decay_channel), 'ct', 'cg', scandir = scandir
        )
    scan.title = differentials.core.standard_titles.get(decay_channel, decay_channel)
    scan.color = 1
    scan.deltaNLL_threshold = -10.
    scan.read()
    return scan

@flag_as_option
def projection_ktcg_plot_floatingBRs(args):
    differentials.plotting.pywrappers.CMS_Latex_lumi.CMS_lumi = 3000.
    scans = []

    hzz = latest_floatingBRs(args, decay_channel='hzz')
    hzz.color = differentials.core.safe_colors.blue
    scans.append(hzz)

    plot = differentials.plotting.plots.MultiContourPlot(
        ('projection_ktcg_plot_floatingBRs'
            + ('_asimov' if args.asimov else '')
            + ('_scenario2' if args.scenario2 else '')
            ),
        scans,
        x_title = differentials.core.standard_titles['ct'], y_title = differentials.core.standard_titles['cg'],
        )
    plot.draw()




