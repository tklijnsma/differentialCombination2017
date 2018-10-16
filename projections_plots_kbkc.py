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

from math import sqrt

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
couplingdependentBRs.scenario2 = differentials.core.AttrDict()

couplingdependentBRs.hzz = 'out/Scan_projection_kbkc_Oct02_hzz_couplingdependentBRs_asimov'
couplingdependentBRs.hgg = [
    'out/Scan_projection_kbkc_Oct02_hgg_couplingdependentBRs_asimov',
    'out/Scan_projection_kbkc_Oct02_hgg_couplingdependentBRs_asimov_rescan_Oct03'
    ]
couplingdependentBRs.combination = [
    'out/Scan_projection_kbkc_Oct02_combination_couplingdependentBRs_asimov',
    'out/Scan_projection_kbkc_Oct02_combination_couplingdependentBRs_asimov_rescan_Oct04',
    ]

couplingdependentBRs.scenario2.hzz = 'out/Scan_projection_kbkc_Oct02_hzz_scenario2_couplingdependentBRs_asimov'
couplingdependentBRs.scenario2.hgg = 'out/Scan_projection_kbkc_Oct04_hgg_scenario2_couplingdependentBRs_asimov'
couplingdependentBRs.scenario2.combination = [
    'out/Scan_projection_kbkc_Oct02_combination_scenario2_couplingdependentBRs_asimov',
    'out/Scan_projection_kbkc_Oct02_combination_scenario2_couplingdependentBRs_asimov_rescan_Oct03'
    ]


class PlotPatcher(object):
    """docstring for PlotPatcher"""

    colors = {
        'hgg' : differentials.core.safe_colors.red,
        'hzz' : differentials.core.safe_colors.blue,
        'hbb' : differentials.core.safe_colors.green,
        'combination' : 1,
        'combWithHbb' : 1,
        }
    x_sm=1.0
    y_sm=1.0
    x_var = 'kappac'
    y_var = 'kappab'

    def __init__(self, args, decay_channel):
        super(PlotPatcher, self).__init__()
        args = differentialutils.force_asimov(args)
        differentials.plotting.pywrappers.CMS_Latex_lumi.CMS_lumi = 3000.
        self.args = args
        self.is_scenario2 = args.scenario2
        self.decay_channel = decay_channel
        self.quick_draw_list = []

        self.scan = self.get_scan()
        self.scan.name = self.decay_channel
        self.scan.color = self.colors[decay_channel]
        if not(self.is_combination()): self.scan.only_1sigma_contours = True
        self.unpatched = self.scan.to_hist()
        self.unpatched.name = self.decay_channel
        self.add_to_quick_draw_list('rawscan', self.unpatched)

    def get_scan(self):
        scandir = self.get_scandict()[self.decay_channel]
        scan = differentials.scans.Scan2D(
            self.decay_channel, self.x_var, self.y_var, scandir = scandir
            )
        scan.read_one_scandir = False
        scan.title = differentials.core.standard_titles.get(self.decay_channel, self.decay_channel)
        scan.deltaNLL_threshold = -10.
        scan.read()
        return scan

    def add_to_quick_draw_list(self, name, H):
        self.quick_draw_list.append((name, H))

    def is_combination(self):
        return not(self.decay_channel in ['hgg', 'hzz', 'hbb'])

    def patch(self):
        if self.is_combination():
            H = self.patch_combination()
        else:
            H = getattr(self, 'patch_' + self.decay_channel)()
        H.name = self.decay_channel
        return H

    def patch_hzz(self):
        return self.scan.to_hist()

    def patch_hgg(self):
        return self.scan.to_hist()

    def patch_combination(self):
        return self.scan.to_hist()

    def get_better_minimum_est(self):
        dr2 = lambda entry: (entry.x - self.x_sm)**2 + (entry.y - self.y_sm)**2
        bestfit = self.scan.bestfit()
        entries = [ e for e in self.scan.entries if not(e.x == bestfit.x and e.y == bestfit.y) ]
        entries.sort(key=dr2)
        closest = entries[:47]
        close_deltaNLLs = [ e.deltaNLL for e in closest ]
        min_delta_nll = min(close_deltaNLLs)
        logging.info(
            'min deltaNLL = {0}, deltaNLL values close to SM: {1}'
            .format(min_delta_nll, close_deltaNLLs)
            )
        return min_delta_nll

    def shift_scan(self):
        new_min_delta_nll = self.get_better_minimum_est()
        logging.info('Shifting to new minimum at deltaNLL = {0}'.format(new_min_delta_nll))
        for i_entry in xrange(len(self.scan.entries)):
            if self.scan.entries[i_entry] == self.scan.bestfit(): continue
            self.scan.entries[i_entry].deltaNLL -= new_min_delta_nll + 0.00001

    def quick_patched_vs_unpatched(self):
        self.add_to_quick_draw_list('zpatched', self.patch())
        for name, H in self.quick_draw_list:
            plotname = 'patchcomp_{0}_{1}'.format(self.decay_channel, name)
            if self.args.scenario2: plotname += '_scenario2'
            H.quickplot(plotname)


class PlotPatcherKBKC_couplingdependentBRs(PlotPatcher):
    """docstring for PlotPatcherKBKC_couplingdependentBRs"""
    def __init__(self, args, decay_channel):
        super(PlotPatcherKBKC_couplingdependentBRs, self).__init__(args, decay_channel)

    def get_scandict(self):
        return couplingdependentBRs.scenario2 if self.is_scenario2 else couplingdependentBRs

    def patch_hzz(self):
        self.shift_scan()
        return self.scan.to_hist()

    def patch_combination(self):
        self.shift_scan()
        self.add_to_quick_draw_list('aftershifting', self.scan.to_hist())

        self.scan.filter(
            lambda e: (
                ( (e.y**2) / (1.1**2)  +  (e.x**2) / (4.5**2) < 1.2**2 )
                and
                ( (e.y**2) / (1.1**2)  +  (e.x**2) / (4.5**2) > 0.8**2 )
                )
            )
        self.add_to_quick_draw_list('ovalfilter', self.scan.to_hist())

        # spline = self.scan.to_spline(
        #     x_min = -5.0, x_max = 5.4, y_min = -1.20, y_max = 1.20,
        #     deltaNLL_cutoff = 10.,
        #     eps = 0.9
        #     )
        # return spline.to_hist(nx=180, ny=180)
        return self.scan.to_hist()

    def patch_hgg(self):
        self.shift_scan()
        H = self.scan.to_hist()
        H.polyfit_patch(x_min=-4.5, x_max=-2.25, y_min=-1.2, y_max=-0.4, order=6)
        H.polyfit_patch(x_min=-6.1, x_max=-3.8, y_min=-0.45, y_max=0.15, order=6)
        H.polyfit_patch(x_min=-3.4, x_max=0.0, y_min=-1.3, y_max=-0.75, order=6)
        H.set_value_for_patch(0., -2.8, 2.6, -0.65, 0.65)
        return H


class PlotPatcherKBKC_couplingdependentBRs_scenario2(PlotPatcherKBKC_couplingdependentBRs):
    """docstring for PlotPatcherKBKC_couplingdependentBRs_scenario2"""
    def __init__(self, args, decay_channel):
        super(PlotPatcherKBKC_couplingdependentBRs_scenario2, self).__init__(args, decay_channel)
        
    def patch_hgg(self):
        self.shift_scan()
        return self.scan.to_hist()

    def patch_combination(self):
        self.shift_scan()
        self.add_to_quick_draw_list('shifted', self.scan.to_hist())

        # spline = self.scan.to_spline(
        #     x_min = -5.0, x_max = 5.4, y_min = -1.15, y_max = 1.15,
        #     # deltaNLL_cutoff = 30.,
        #     eps = 0.7
        #     )
        # return spline.to_hist(nx=180, ny=180)

        return self.scan.to_hist()


@flag_as_option
def projection_kbkc_plot_comparison(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    Patcher = PlotPatcherKBKC_couplingdependentBRs_scenario2 if args.scenario2 else PlotPatcherKBKC_couplingdependentBRs
    Patcher(args, decay_channel).quick_patched_vs_unpatched()


@flag_as_option
def projection_kbkc_plot_couplingdependentBRs(args):
    Patcher = PlotPatcherKBKC_couplingdependentBRs_scenario2 if args.scenario2 else PlotPatcherKBKC_couplingdependentBRs
    scans = []

    combination = Patcher(args, 'combination').patch()
    scans.append(combination)

    hzz = Patcher(args, 'hzz').patch()
    scans.append(hzz)

    hgg = Patcher(args, 'hgg').patch()
    scans.append(hgg)

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
floatingBRs.scenario2 = differentials.core.AttrDict()

floatingBRs.hzz = 'out/Scan_projection_kbkc_Oct03_hzz_floatingBRs_asimov'
floatingBRs.hgg = 'out/Scan_projection_kbkc_Oct03_hgg_floatingBRs_asimov'
floatingBRs.combination = 'out/Scan_projection_kbkc_Oct03_combination_floatingBRs_asimov'

floatingBRs.scenario2.hzz = 'out/Scan_projection_kbkc_Oct03_hzz_scenario2_floatingBRs_scenario2_asimov'
floatingBRs.scenario2.hgg = 'out/Scan_projection_kbkc_Oct03_hgg_scenario2_floatingBRs_scenario2_asimov'
floatingBRs.scenario2.combination = 'out/Scan_projection_kbkc_Oct03_combination_scenario2_floatingBRs_scenario2_asimov'


class PlotPatcherKBKC_floatingBRs(PlotPatcher):
    """docstring for PlotPatcherKBKC_floatingBRs"""
    def __init__(self, args, decay_channel):
        super(PlotPatcherKBKC_floatingBRs, self).__init__(args, decay_channel)

    def get_scandict(self):
        return floatingBRs.scenario2 if self.is_scenario2 else floatingBRs

    def patch_hgg(self):
        self.shift_scan()
        return self.scan.to_hist()

    def patch_hzz(self):
        self.shift_scan()
        return self.scan.to_hist()

    def patch_combination(self):
        self.shift_scan()
        return self.scan.to_hist()



class PlotPatcherKBKC_floatingBRs_scenario2(PlotPatcherKBKC_floatingBRs):
    """docstring for PlotPatcherKBKC_floatingBRs_scenario2"""
    def __init__(self, args, decay_channel):
        super(PlotPatcherKBKC_floatingBRs_scenario2, self).__init__(args, decay_channel)



@flag_as_option
def projection_kbkc_plot_floatingBRs(args):
    Patcher = PlotPatcherKBKC_floatingBRs_scenario2 if args.scenario2 else PlotPatcherKBKC_floatingBRs
    scans = []

    # combination = Patcher(args, 'combination').unpatched
    combination = Patcher(args, 'combination').patch()
    scans.append(combination)

    # hzz = Patcher(args, 'hzz').unpatched
    hzz = Patcher(args, 'hzz').patch()
    scans.append(hzz)

    # hgg = Patcher(args, 'hgg').unpatched
    hgg = Patcher(args, 'hgg').patch()
    scans.append(hgg)

    plot = differentials.plotting.plots.MultiContourPlot(
        ('projection_kbkc_plot_floatingBRs'
            + ('_asimov' if args.asimov else '')
            + ('_scenario2' if args.scenario2 else '')
            ),
        scans,
        x_title = differentials.core.standard_titles['kappac'], y_title = differentials.core.standard_titles['kappab'],
        # x_min=couplingdependentBRs_x_min, x_max=couplingdependentBRs_x_max, y_min=couplingdependentBRs_y_min, y_max=couplingdependentBRs_y_max
        )
    plot.draw()


# @flag_as_option
# def projection_kbkc_plot_floatingBRs(args):
#     differentials.plotting.pywrappers.CMS_Latex_lumi.CMS_lumi = 3000.
#     scans = []

#     combination = latest_floatingBRs(args, decay_channel='combination')
#     combination.color = 1
#     combination_hist = combination.to_hist()
#     combination_hist.add_offset_to_zero()
#     scans.append(combination_hist)

#     hgg = latest_floatingBRs(args, decay_channel='hgg')
#     hgg.color = differentials.core.safe_colors.red
#     hgg_hist = hgg.to_hist()
#     hgg_hist.add_offset_to_zero()
#     scans.append(hgg_hist)

#     hzz = latest_floatingBRs(args, decay_channel='hzz')
#     hzz.color = differentials.core.safe_colors.blue
#     hzz_hist = hzz.to_hist()
#     hzz_hist.add_offset_to_zero()
#     scans.append(hzz_hist)

#     plot = differentials.plotting.plots.MultiContourPlot(
#         ('projection_kbkc_plot_floatingBRs'
#             + ('_asimov' if args.asimov else '')
#             + ('_scenario2' if args.scenario2 else '')
#             ),
#         scans,
#         x_title = differentials.core.standard_titles['kappac'], y_title = differentials.core.standard_titles['kappab'],
#         x_min = -17., x_max = 17., y_min = -2.5, y_max = 4.5,
#         )
#     plot.draw()




