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


# out/Scan_projection_kbkc_Jul18_combination_floatingBRs_scenario2_asimov
# out/Scan_projection_kbkc_Jul18_hgg_floatingBRs_scenario2_asimov

# out/Scan_projection_kbkc_Jul18_hzz_couplingdependentBRs_scenario2_asimov
# out/Scan_projection_kbkc_Jul18_hzz_floatingBRs_scenario2_asimov_0


#____________________________________________________________________
couplingdependentBRs = differentials.core.AttrDict()
couplingdependentBRs.hzz         = 'out/Scan_projection_kbkc_Jul10_hzz_couplingdependentBRs_asimov_0'
# couplingdependentBRs.hgg         = 'out/Scan_projection_kbkc_Jul10_hgg_couplingdependentBRs_asimov'
# couplingdependentBRs.combination = 'out/Scan_projection_kbkc_Jul10_combination_couplingdependentBRs_asimov'
couplingdependentBRs.hgg         = 'out/Scan_projection_kbkc_Aug09_hgg_couplingdependentBRs_asimov'
couplingdependentBRs.combination = 'out/Scan_projection_kbkc_Aug09_combination_couplingdependentBRs_asimov'

couplingdependentBRs.scenario2 = differentials.core.AttrDict()
couplingdependentBRs.scenario2.hzz         = 'out/Scan_projection_kbkc_Jul18_hzz_couplingdependentBRs_scenario2_asimov'
# couplingdependentBRs.scenario2.hgg         = 'out/Scan_projection_kbkc_Jul17_hgg_couplingdependentBRs_scenario2_asimov_0'
# couplingdependentBRs.scenario2.combination = 'out/Scan_projection_kbkc_Jul17_combination_couplingdependentBRs_scenario2_asimov_0'
# couplingdependentBRs.scenario2.hgg         = 'out/Scan_projection_kbkc_Jul18_combination_couplingdependentBRs_scenario2_asimov'
# couplingdependentBRs.scenario2.combination = 'out/Scan_projection_kbkc_Jul18_hgg_couplingdependentBRs_scenario2_asimov'
couplingdependentBRs.scenario2.hgg         = 'out/Scan_projection_kbkc_Aug09_hgg_couplingdependentBRs_scenario2_asimov'
couplingdependentBRs.scenario2.combination = 'out/Scan_projection_kbkc_Aug09_combination_couplingdependentBRs_scenario2_asimov'


def latest_couplingdependentBRs(args, decay_channel=None):
    args = differentialutils.force_asimov(args)
    if decay_channel is None: decay_channel = differentialutils.get_decay_channel_tag(args)
    scandict = couplingdependentBRs.scenario2 if args.scenario2 else couplingdependentBRs
    scandir = scandict[decay_channel]

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
    hzz.only_1sigma_contours = True
    scans.append(hzz)

    combination = latest_couplingdependentBRs(args, decay_channel='combination')
    combination.color = 1
    scans.append(combination)

    hgg = latest_couplingdependentBRs(args, decay_channel='hgg')
    hgg.color = differentials.core.safe_colors.red
    hgg.only_1sigma_contours = True

    hgg_hist = hgg.to_hist()
    hgg_hist.set_value_for_patch(0., -3., 3., -0.75, 0.75)
    scans.append(hgg_hist)


    plot = differentials.plotting.plots.MultiContourPlot(
        ('projection_kbkc_plot_couplingdependentBRs'
            + ('_asimov' if args.asimov else '')
            + ('_scenario2' if args.scenario2 else '')
            ),
        scans,
        x_title = differentials.core.standard_titles['kappac'], y_title = differentials.core.standard_titles['kappab'],
        # x_min=couplingdependentBRs_x_min, x_max=couplingdependentBRs_x_max, y_min=couplingdependentBRs_y_min, y_max=couplingdependentBRs_y_max
        )
    plot.draw()


#____________________________________________________________________
floatingBRs = differentials.core.AttrDict()
floatingBRs.hzz         = 'out/Scan_projection_kbkc_Jul10_hzz_floatingBRs_asimov_0'
# floatingBRs.hgg         = 'out/Scan_projection_kbkc_Jul10_hgg_floatingBRs_asimov'
# floatingBRs.combination = 'out/Scan_projection_kbkc_Jul10_combination_floatingBRs_asimov'
floatingBRs.hgg         = 'out/Scan_projection_kbkc_Aug09_hgg_floatingBRs_asimov'
floatingBRs.combination = 'out/Scan_projection_kbkc_Aug09_combination_floatingBRs_asimov'

floatingBRs.scenario2 = differentials.core.AttrDict()
floatingBRs.scenario2.hzz         = 'out/Scan_projection_kbkc_Jul18_hzz_floatingBRs_scenario2_asimov_0'
# floatingBRs.scenario2.hgg         = 'out/Scan_projection_kbkc_Jul17_hgg_floatingBRs_scenario2_asimov'
# floatingBRs.scenario2.combination = 'out/Scan_projection_kbkc_Jul17_combination_floatingBRs_scenario2_asimov'
# floatingBRs.scenario2.hgg         = 'out/Scan_projection_kbkc_Jul18_hgg_floatingBRs_scenario2_asimov'
# floatingBRs.scenario2.combination = 'out/Scan_projection_kbkc_Jul18_combination_floatingBRs_scenario2_asimov'
floatingBRs.scenario2.hgg         = 'out/Scan_projection_kbkc_Aug09_hgg_floatingBRs_scenario2_asimov'
floatingBRs.scenario2.combination = 'out/Scan_projection_kbkc_Aug09_combination_floatingBRs_scenario2_asimov'



# out/Scan_projection_kbkc_Jul18_combination_floatingBRs_scenario2_asimov
# out/Scan_projection_kbkc_Jul18_hgg_floatingBRs_scenario2_asimov


def latest_floatingBRs(args, decay_channel=None, asimov=None, splined=False):
    args = differentialutils.force_asimov(args)
    if decay_channel is None: decay_channel = differentialutils.get_decay_channel_tag(args)
    scandict = floatingBRs.scenario2 if args.scenario2 else floatingBRs
    scandir = scandict[decay_channel]

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

    combination = latest_floatingBRs(args, decay_channel='combination')
    combination.color = 1
    combination_hist = combination.to_hist()
    combination_hist.add_offset_to_zero()
    scans.append(combination_hist)

    hgg = latest_floatingBRs(args, decay_channel='hgg')
    hgg.color = differentials.core.safe_colors.red
    hgg_hist = hgg.to_hist()
    hgg_hist.add_offset_to_zero()
    scans.append(hgg_hist)

    hzz = latest_floatingBRs(args, decay_channel='hzz')
    hzz.color = differentials.core.safe_colors.blue
    hzz_hist = hzz.to_hist()
    hzz_hist.add_offset_to_zero()
    scans.append(hzz_hist)

    plot = differentials.plotting.plots.MultiContourPlot(
        ('projection_kbkc_plot_floatingBRs'
            + ('_asimov' if args.asimov else '')
            + ('_scenario2' if args.scenario2 else '')
            ),
        scans,
        x_title = differentials.core.standard_titles['kappac'], y_title = differentials.core.standard_titles['kappab'],
        x_min = -17., x_max = 17., y_min = -2.5, y_max = 4.5,
        )
    plot.draw()




