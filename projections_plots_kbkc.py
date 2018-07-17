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
couplingdependentBRs.hzz         = 'out/Scan_projection_kbkc_Jul10_hzz_couplingdependentBRs_asimov_0'
couplingdependentBRs.hgg         = 'out/Scan_projection_kbkc_Jul10_hgg_couplingdependentBRs_asimov'
couplingdependentBRs.combination = 'out/Scan_projection_kbkc_Jul10_combination_couplingdependentBRs_asimov'

def latest_couplingdependentBRs(args, decay_channel=None):
    args = differentialutils.force_asimov(args)
    if decay_channel is None: decay_channel = differentialutils.get_decay_channel_tag(args)
    scandir = couplingdependentBRs[decay_channel]

    scan = differentials.scans.Scan2D(
        'kbkc_{0}'.format(decay_channel), 'kappac', 'kappab', scandir = scandir
        )
    scan.title = differentials.core.standard_titles.get(decay_channel, decay_channel)
    scan.color = 1
    scan.deltaNLL_threshold = -10.
    scan.read()
    return scan


@flag_as_option
def projection_kbkc_plot_couplingdependentBRs(args):
    differentials.plotting.pywrappers.CMS_Latex_lumi.CMS_lumi = 3000.
    scans = []

    hzz = latest_couplingdependentBRs(args, decay_channel='hzz')
    hzz.color = differentials.core.safe_colors.blue
    scans.append(hzz)

    # combination = latest_couplingdependentBRs(args, decay_channel='combination')
    # combination.color = 1
    # scans.append(combination)

    # hgg = latest_couplingdependentBRs(args, decay_channel='hgg')
    # hgg.color = differentials.core.safe_colors.red
    # scans.append(hgg)

    plot = differentials.plotting.plots.MultiContourPlot(
        'projection_kbkc_plot_couplingdependentBRs' + ('_asimov' if args.asimov else ''),
        scans,
        x_title = differentials.core.standard_titles['kappac'], y_title = differentials.core.standard_titles['kappab'],
        # x_min=couplingdependentBRs_x_min, x_max=couplingdependentBRs_x_max, y_min=couplingdependentBRs_y_min, y_max=couplingdependentBRs_y_max
        )
    plot.draw()


#____________________________________________________________________
floatingBRs = differentials.core.AttrDict()
floatingBRs.hzz         = 'out/Scan_projection_kbkc_Jul10_hzz_floatingBRs_asimov_0'
floatingBRs.hgg         = 'out/Scan_projection_kbkc_Jul10_hgg_floatingBRs_asimov'
floatingBRs.combination = 'out/Scan_projection_kbkc_Jul10_combination_floatingBRs_asimov'

def latest_floatingBRs(args, decay_channel=None, asimov=None, splined=False):
    args = differentialutils.force_asimov(args)
    if decay_channel is None: decay_channel = differentialutils.get_decay_channel_tag(args)
    scandir = floatingBRs[decay_channel]

    scan = differentials.scans.Scan2D(
        'kbkc_{0}'.format(decay_channel), 'kappac', 'kappab', scandir = scandir
        )
    scan.title = differentials.core.standard_titles.get(decay_channel, decay_channel)
    scan.color = 1
    scan.deltaNLL_threshold = -10.
    scan.read()
    return scan

@flag_as_option
def projection_kbkc_plot_floatingBRs(args):
    differentials.plotting.pywrappers.CMS_Latex_lumi.CMS_lumi = 3000.
    scans = []

    hzz = latest_floatingBRs(args, decay_channel='hzz')
    hzz.color = differentials.core.safe_colors.blue
    scans.append(hzz)

    plot = differentials.plotting.plots.MultiContourPlot(
        'projection_kbkc_plot_floatingBRs' + ('_asimov' if args.asimov else ''),
        scans,
        x_title = differentials.core.standard_titles['kappac'], y_title = differentials.core.standard_titles['kappab'],
        )
    plot.draw()




