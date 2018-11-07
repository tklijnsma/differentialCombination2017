#!/usr/bin/env python
"""
Thomas Klijnsma
"""

#____________________________________________________________________
# Imports

import logging
import os, sys, re, copy, math
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
couplingdependentBRs.hzz = [
    # 'out/Scan_projection_ktcg_Oct04_hzz_couplingdependentBRs_asimov'
    # Frozen lumiScale:
    'out/Scan_projection_ktcg_Oct17_hzz_couplingdependentBRs_asimov'
    ]
couplingdependentBRs.hgg = [
    # 'out/Scan_projection_ktcg_Oct05_hgg_couplingdependentBRs_asimov',
    # 'out/Scan_projection_ktcg_Oct05_hgg_couplingdependentBRs_asimov_rescan_Oct08'
    # Frozen lumiScale:
    'out/Scan_projection_ktcg_Oct17_hgg_couplingdependentBRs_asimov'
    ]

couplingdependentBRs.combWithHbb = [
    # 'out/Scan_projection_ktcg_Oct04_combWithHbb_couplingdependentBRs_asimov',
    # Need to rerun... but probably with more clever settings
    # 'out/Scan_projection_ktcg_Oct12_combWithHbb_couplingdependentBRs_asimov'
    # 'out/Scan_projection_ktcg_Oct13_combWithHbb_couplingdependentBRs_asimov'
    # 'out/Scan_projection_ktcg_Oct13_combWithHbb_couplingdependentBRs_asimov_0'
    # 'out/Scan_projection_ktcg_Oct16_combWithHbb_couplingdependentBRs_asimov'
    # 'out/Scan_projection_ktcg_Oct17_combWithHbb_couplingdependentBRs_asimov_0'
    # 'out/Scan_projection_ktcg_Oct18_combWithHbb_couplingdependentBRs_asimov'
    'out/Scan_projection_ktcg_Oct18_combWithHbb_couplingdependentBRs_asimov_0'
    ]

couplingdependentBRs.scenario2 = differentials.core.AttrDict()
couplingdependentBRs.scenario2.hzz = [
    # 'out/Scan_projection_ktcg_Oct04_hzz_scenario2_couplingdependentBRs_asimov_0'
    # Frozen lumiScale:
    'out/Scan_projection_ktcg_Oct17_hzz_scenario2_couplingdependentBRs_asimov'
    ]
couplingdependentBRs.scenario2.hgg = 'out/Scan_projection_ktcg_Aug09_hgg_couplingdependentBRs_scenario2_asimov'
couplingdependentBRs.scenario2.combWithHbb = [
    # 'out/Scan_projection_ktcg_Oct04_combWithHbb_scenario2_couplingdependentBRs_asimov',
    # 'out/Scan_projection_ktcg_Oct04_combWithHbb_scenario2_couplingdependentBRs_asimov_rescan_Oct06',
    # Rescan not enough
    # Need to rerun... but probably with more clever settings
    # Frozen lumiScale
    # 'out/Scan_projection_ktcg_Oct17_combWithHbb_scenario2_couplingdependentBRs_asimov_0'
    # 'out/Scan_projection_ktcg_Oct18_combWithHbb_scenario2_couplingdependentBRs_asimov'
    'out/Scan_projection_ktcg_Oct19_combWithHbb_scenario2_couplingdependentBRs_asimov' # Slightly better ranges
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

    def preprocess_combination(self):
        self.deltaNLL_threshold = -1.

    # def patch_combination(self):
    #     self.shift_scan()
    #     H = self.scan.to_hist()
    #     H.add_padding(
    #         700.,
    #         x_min = 0.72,
    #         x_max = 1.44,
    #         y_min = -0.035,
    #         y_max = 0.035
    #         )
    #     return H

    def patch_combination(self):
        self.shift_scan()
        fitter = self.scan
        fitter.__class__ = KTCGFitter
        fitter.fill_bestfit(1.0, 0.0)
        fitter.make_bin_centers()
        x_bin_centers_unrotated = fitter.x_bin_centers
        y_bin_centers_unrotated = fitter.y_bin_centers

        fitter.filter_points()
        fitter.set_rotation_pars(x1=1.0, y1=0.0, x2=1.144, y2=-0.0100)
        fitter.rotate()
        fitter.set_extrema()
        fitter.get_polyfit(order=5)

        fitter.fill_polyfit_unrotated(
            x_min = 0.75, x_max = 1.45,
            y_min = -0.019, y_max = 0.018,
            n_bins = 200
            )

        H = fitter.to_hist()
        H.add_padding(
            700.,
            x_min = 0.7,
            x_max = 1.46,
            y_min = -0.038,
            y_max = 0.054,
            )
        H.set_value_for_patch(700.0, x_min=0.65, x_max=0.9, y_min=-0.02, y_max=-0.01)
        H.set_value_for_patch(700.0, x_min=1.34, x_max=1.42, y_min=-0.01, y_max=0.01)
        return H


def print_bin(H, i):
    print 'Bin {0}:'.format(i)
    print '  left:   {0}'.format(H.H2.GetXaxis().GetBinLowEdge(i))
    print '  center: {0}'.format(H.H2.GetXaxis().GetBinCenter(i))
    print '  right:  {0}'.format(H.H2.GetXaxis().GetBinUpEdge(i))


class PlotPatcherKTCG_couplingdependentBRs_scenario2(PlotPatcherKTCG_couplingdependentBRs):
    """docstring for PlotPatcherKTCG_couplingdependentBRs_scenario2"""

    def preprocess_combination(self):
        self.deltaNLL_threshold = -10.

    def patch_combination(self):
        self.shift_scan()
        # return self.scan.to_hist()

        # self.scan.filter(lambda e:
        #     e.cg > 0.005 -  (0.012/0.2) * e.ct
        #     )

        fitter = self.scan
        fitter.__class__ = KTCGFitter
        fitter.fill_bestfit(1.0, 0.0)
        fitter.make_bin_centers()
        x_bin_centers_unrotated = fitter.x_bin_centers
        y_bin_centers_unrotated = fitter.y_bin_centers

        fitter.filter_points()
        fitter.set_rotation_pars(x1=0.905, y1=0.0070, x2=1.0, y2=0.0)
        fitter.rotate()
        fitter.set_extrema()
        fitter.get_polyfit(order=4)

        fitter.fill_polyfit_unrotated(
            # x_min = 0.8, x_max = 1.2,
            # y_min = -0.014, y_max = 0.014,
            x_min = 0.773357152939,
            x_max = 1.23664283752,
            y_min = -0.0137999998406,
            y_max = 0.0137999998406,
            n_bins = 200
            )

        H = fitter.to_hist()
        H.add_padding(
            700.,
            x_min = 0.7,
            x_max = 1.46,
            y_min = -0.038,
            y_max = 0.054,
            )
        
        H.set_value_for_patch(700.0,
            x_min=0.7, x_max=0.923, y_min=-0.02, y_max=-0.0015
            )        
        H.set_value_for_patch(700.0,
            x_min=1.14, x_max=1.3, y_min=0.005, y_max=0.02
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
    x_min = combWithHbb.x_bin_boundaries[0]
    x_max = combWithHbb.x_bin_boundaries[-1]
    y_min = combWithHbb.y_bin_boundaries[0]
    y_max = combWithHbb.y_bin_boundaries[-1]

    plot = differentials.plotting.plots.MultiContourPlot(
        ('projection_ktcg_plot_couplingdependentBRs'
            + ('_asimov' if args.asimov else '')
            + ('_scenario2' if args.scenario2 else '')
            ),
        scans,
        x_title = differentials.core.standard_titles['ct'], y_title = differentials.core.standard_titles['cg'],
        x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max
        )
    plot.draw(wait=True)
    plot.add_BR_parametrized_text(
        # x = lambda c: 1. - c.GetRightMargin() - 0.01,
        x = lambda c: c.GetLeftMargin() + 0.03,
        y = lambda c: c.GetBottomMargin() + 0.03 + 0.015 + 0.07
        )
    plot.add_text_on_the_fly(
        text_size = 0.040,
        x = lambda c: 1. - c.GetRightMargin() - 0.01,
        y = lambda c: 1. - c.GetTopMargin() - 0.05,
        text = 'w/ YR18 syst. uncert. (S2)' if args.scenario2 else 'w/ Run 2 syst. uncert. (S1)'
        )

    plot.base.GetYaxis().SetTitleOffset(1.1)
    plot.legend.SetBorderSize(0)
    plot.cdl.SetBorderSize(0)

    plot.wrapup()





@flag_as_option
def projection_ktcg_fittesting(args):
    Patcher = PlotPatcherKTCG_couplingdependentBRs_scenario2 if args.scenario2 else PlotPatcherKTCG_couplingdependentBRs
    scans = []

    combWithHbb = Patcher(args, 'combWithHbb')
    combWithHbb.shift_scan()

    fitter = combWithHbb.scan
    fitter.__class__ = KTCGFitter
    fitter.fill_bestfit(1.0, 0.0)
    fitter.make_bin_centers()

    x_bin_centers_unrotated = fitter.x_bin_centers
    y_bin_centers_unrotated = fitter.y_bin_centers

    fitter.plot('ktcg_fittesting_base')

    fitter.filter_points()
    fitter.set_rotation_pars(x1=1.0, y1=0.0, x2=1.144, y2=-0.0100)
    if args.scenario2: fitter.set_rotation_pars(x1=0.905, y1=0.0070, x2=1.0, y2=0.0)

    fitter.rotate()

    fitter.plot_rotated('ktcg_fittesting_rotated')

    fitter.set_extrema()
    fitter.get_polyfit()
    fitter.apply_polyfit()
    fitter.plot_rotated('ktcg_fittesting_rotpolyfit')

    fitter.fill_polyfit()
    fitter.plot_rotated('ktcg_fittesting_rotpolyfitfilled')

    fitter.unrotate()
    fitter.set_bin_centers(x_bin_centers_unrotated, y_bin_centers_unrotated)
    fitter.rebin_entries()
    fitter.plot('ktcg_fittesting_unrotpolyfit')

    print 'x_min = {0}'.format(min(x_bin_centers_unrotated))
    print 'x_max = {0}'.format(max(x_bin_centers_unrotated))
    print 'y_min = {0}'.format(min(y_bin_centers_unrotated))
    print 'y_max = {0}'.format(max(y_bin_centers_unrotated))

    fitter.fill_polyfit_unrotated(
        x_min = min(x_bin_centers_unrotated), x_max = max(x_bin_centers_unrotated),
        y_min = min(y_bin_centers_unrotated), y_max = max(y_bin_centers_unrotated),
        n_bins = 100
        )
    fitter.plot('ktcg_fittesting_unrotpolyfit_fine')



class KTCGFitter(differentials.scans.Scan2D):
    """docstring for KTCGFitter"""

    x_SM = 1.0
    y_SM = 0.0

    def make_bin_centers(self):
        self.x_bin_centers = list(set([ entry.x for entry in self.entries if entry.x != self.x_SM ]))
        self.x_bin_centers.sort()
        self.y_bin_centers = list(set([ entry.y for entry in self.entries if entry.y != self.y_SM ]))
        self.y_bin_centers.sort()

    def set_bin_centers(self, x_bin_centers, y_bin_centers):
        self.x_bin_centers = x_bin_centers
        self.y_bin_centers = y_bin_centers

    def filter_points(self):
        self.filter(
            lambda e: e.deltaNLL < 30.,
            inplace = True
            )

    def set_rotation_pars(self, x1, y1, x2, y2):
        if x1 > x2: x1, y1, x2, y2 = (x2, y2, x1, y1)
        self.a = (y2-y1)/(x2-x1)
        self.b = y2 - self.a * x2
        self.phi = -math.atan(self.a)

    def rotate(self):
        for entry in self.entries:
            entry.x, entry.y = self.get_rotated_xy(entry.x, entry.y)

    def unrotate(self):
        for entry in self.entries:
            entry.x, entry.y = self.get_unrotated_xy(entry.x, entry.y)

    def get_rotated_xy(self, x, y):
        dx = x - self.x_SM
        dy = y - self.y_SM
        dx_new = math.cos(self.phi) * dx - math.sin(self.phi) * dy
        dy_new = math.sin(self.phi) * dx + math.cos(self.phi) * dy
        return dx_new, dy_new

    def get_unrotated_xy(self, x, y):
        x_new = math.cos(-self.phi) * x - math.sin(-self.phi) * y
        y_new = math.sin(-self.phi) * x + math.cos(-self.phi) * y
        x_new += self.x_SM
        y_new += self.y_SM
        return x_new, y_new

    def rebin_entries(self):
        for entry in self.entries:
            logging.trace('x: Searching for {0}'.format(entry.x))
            entry.x = find_closest(entry.x, self.x_bin_centers)
            logging.trace('    Found {0}'.format(entry.x))

            logging.trace('y: Searching for {0}'.format(entry.y))
            entry.y = find_closest(entry.y, self.y_bin_centers)
            logging.trace('    Found {0}'.format(entry.y))

    def set_extrema(self):
        self.x_min = min([e.x for e in self.entries])
        self.x_max = max([e.x for e in self.entries])
        self.y_min = min([e.y for e in self.entries])
        self.y_max = max([e.y for e in self.entries])

    def plot_rotated(self, plotname):
        clone = copy.deepcopy(self)
        clone.set_extrema()
        clone.set_bin_centers(
            x_bin_centers = differentials.core.get_axis(45, clone.x_min, clone.x_max),
            y_bin_centers = differentials.core.get_axis(45, clone.y_min, clone.y_max),
            )
        clone.rebin_entries()
        clone.plot(plotname)

    def get_polyfit(self, order=4, no_padding=False):
        T2D = ROOT.TGraph2D()
        ROOT.SetOwnership(T2D, False)
        for i_entry, entry in enumerate(self.entries):
            T2D.SetPoint(i_entry, entry.x, entry.y, entry.deltaNLL)

        dx = 0.5 * (self.x_max-self.x_min)
        dy = 0.5 * (self.x_max-self.x_min)
        if no_padding: dx, dy = ( 0., 0. )
        polyfit_factory = differentials.spline2d.Spline2DFactory()
        polyfit_factory.x_min = self.x_min - dx
        polyfit_factory.x_max = self.x_max + dx
        polyfit_factory.y_min = self.y_min - dy
        polyfit_factory.y_max = self.y_max + dy
        polyfit_factory.ord_polynomial = order
        self.polyfit = polyfit_factory.make_polyfit(T2D)
        self.polyfit.multiply_by_two = False # Because included in scan.to_hist()
        self.polyfit.negativity_is_zero = True

    def apply_polyfit(self):
        for entry in self.entries:
            entry.deltaNLL = self.polyfit.eval(entry.x, entry.y)

    def eval_polyfit_unrotated(self, x, y):
        x_rot, y_rot = self.get_rotated_xy(x, y)
        return self.polyfit.eval(x_rot, y_rot)

    def fill_polyfit(self, x_min, x_max, y_min, y_max, n_bins=100):
        entries = []
        x_bin_centers = differentials.core.get_axis(n_bins, x_min, x_max)
        y_bin_centers = differentials.core.get_axis(n_bins, y_min, y_max)
        for x in x_bin_centers:
            for y in y_bin_centers:
                entry = differentials.core.AttrDict(x = x, y = y, deltaNLL = self.polyfit.eval(x,y))
                entries.append(entry)
        self.entries = entries

    def fill_polyfit_mirrored_y(self, x_min, x_max, y_min, y_max, n_bins=100):
        entries = []
        x_bin_centers = differentials.core.get_axis(n_bins, x_min, x_max)
        y_bin_centers = differentials.core.get_axis(n_bins, y_min, y_max)
        for x in x_bin_centers:
            for y in y_bin_centers:
                if x < 0.:
                    deltaNLL = self.polyfit.eval(-x,-y)
                else:
                    deltaNLL = self.polyfit.eval(x,y)
                entry = differentials.core.AttrDict(x = x, y = y, deltaNLL = deltaNLL)
                entries.append(entry)
        self.entries = entries

    def fill_polyfit_unrotated(self, x_min, x_max, y_min, y_max, n_bins=100):
        entries = []
        x_bin_centers = differentials.core.get_axis(n_bins, x_min, x_max)
        y_bin_centers = differentials.core.get_axis(n_bins, y_min, y_max)
        for x in x_bin_centers:
            for y in y_bin_centers:
                deltaNLL = self.eval_polyfit_unrotated(x, y)
                entry = differentials.core.AttrDict(x=x, y=y, deltaNLL=deltaNLL)
                entries.append(entry)
        self.entries = entries


def find_closest(x, x_list):
    return min([((x - i)**2, i) for i in x_list ])[1]






#____________________________________________________________________
floatingBRs = differentials.core.AttrDict()
floatingBRs.scenario2 = differentials.core.AttrDict()

# floatingBRs.hzz = 'out/Scan_projection_ktcg_Oct05_hzz_floatingBRs_asimov_0'
# floatingBRs.hgg = 'out/Scan_projection_ktcg_Oct04_hgg_floatingBRs_asimov'
# floatingBRs.combWithHbb = 'out/Scan_projection_ktcg_Oct06_combWithHbb_floatingBRs_asimov'
floatingBRs.hzz = 'out/Scan_projection_ktcg_Oct19_hzz_floatingBRs_asimov'
floatingBRs.hgg = 'out/Scan_projection_ktcg_Oct19_hgg_floatingBRs_asimov'
# floatingBRs.combWithHbb = 'out/Scan_projection_ktcg_Oct19_combWithHbb_floatingBRs_asimov'
floatingBRs.combWithHbb = 'out/Scan_projection_ktcg_Oct21_combWithHbb_floatingBRs_asimov_0'

# floatingBRs.scenario2.hzz = 'out/Scan_projection_ktcg_Oct04_hzz_scenario2_floatingBRs_asimov'
# floatingBRs.scenario2.hgg = 'out/Scan_projection_ktcg_Oct04_hgg_scenario2_floatingBRs_asimov'
floatingBRs.scenario2.hzz = 'out/Scan_projection_ktcg_Oct19_hzz_scenario2_floatingBRs_asimov'
floatingBRs.scenario2.hgg = 'out/Scan_projection_ktcg_Oct19_hgg_scenario2_floatingBRs_asimov'
# floatingBRs.scenario2.combWithHbb = 'out/Scan_projection_ktcg_Oct19_combWithHbb_scenario2_floatingBRs_asimov'
# floatingBRs.scenario2.combWithHbb = 'out/Scan_projection_ktcg_Oct20_combWithHbb_scenario2_floatingBRs_asimov'
# floatingBRs.scenario2.combWithHbb = 'out/Scan_projection_ktcg_Oct20_combWithHbb_scenario2_floatingBRs_asimov_0'
floatingBRs.scenario2.combWithHbb = 'out/Scan_projection_ktcg_Oct21_combWithHbb_scenario2_floatingBRs_asimov_0'

# Done:
# out/Scan_projection_ktcg_Oct19_combWithHbb_scenario2_couplingdependentBRs_asimov
# out/Scan_projection_ktcg_Oct19_combWithHbb_scenario2_floatingBRs_asimov
# out/Scan_projection_ktcg_Oct19_hzz_floatingBRs_asimov
# out/Scan_projection_ktcg_Oct19_hzz_scenario2_floatingBRs_asimov
# out/Scan_projection_ktcg_Oct19_hgg_floatingBRs_asimov
# out/Scan_projection_kbkc_Oct18_hgg_couplingdependentBRs_asimov_rescan_Oct19
# out/Scan_projection_ktcg_Oct19_hgg_scenario2_floatingBRs_asimov

# Waiting:
# out/Scan_projection_ktcg_Oct19_combWithHbb_floatingBRs_asimov



class PlotPatcherKTCG_floatingBRs(PlotPatcherKTCG):
    """docstring for PlotPatcherKTCG_floatingBRs"""
    def __init__(self, args, decay_channel):
        super(PlotPatcherKTCG_floatingBRs, self).__init__(args, decay_channel)

    def get_scandict(self):
        return floatingBRs.scenario2 if self.args.scenario2 else floatingBRs

    def patch_hgg(self):
        self.shift_scan()
        # return self.scan.to_hist()
        fitter = self.scan
        fitter.__class__ = KTCGFitter
        fitter.fill_bestfit(1.0, 0.0)
        fitter.make_bin_centers()
        fitter.x_min = 0.6
        fitter.x_max = 1.5
        fitter.y_min = -0.04
        fitter.y_max = 0.035
        fitter.filter(
            lambda e: all([
                e.deltaNLL >= 0., e.deltaNLL < 20.,
                e.ct < fitter.x_max, e.ct > fitter.x_min, 
                e.cg < fitter.y_max, e.cg > fitter.y_min, 
                ]),
            inplace = True
            )
        fitter.get_polyfit(order=5, no_padding=True)

        # fitter.fill_polyfit(
        #     x_min = fitter.x_min,
        #     x_max = fitter.x_max,
        #     y_min = fitter.y_min,
        #     y_max = fitter.y_max,
        #     n_bins = 300
        #     )

        fitter.fill_polyfit_mirrored_y(
            x_min = -fitter.x_max,
            x_max = fitter.x_max,
            y_min = -fitter.y_max,
            y_max = fitter.y_max,
            n_bins = 300
            )

        H = fitter.to_hist()

        # H.mirror()
        # H.add_padding(
        #     700.,
        #     x_min = -1.7,
        #     x_max = 1.7,
        #     y_min = -0.045,
        #     y_max = 0.045,
        #     )
        return H

    def patch_combination(self):
        self.shift_scan()
        # return self.scan.to_hist()
        fitter = self.scan
        fitter.__class__ = KTCGFitter
        fitter.fill_bestfit(1.0, 0.0)
        fitter.make_bin_centers()
        fitter.filter_points()
        fitter.set_extrema()
        fitter.get_polyfit(order=5) # 5 is max, 6 fits non-existent features due to missing points

        dx = fitter.x_max - fitter.x_min
        dy = fitter.y_max - fitter.y_min
        fitter.fill_polyfit(
            # x_min = fitter.x_min - 0.2*dx,
            # x_max = fitter.x_max + 0.2*dx,
            x_min = -2.21,
            x_max = 2.21,
            y_min = fitter.y_min - 0.2*dy,
            y_max = fitter.y_max + 0.2*dy,
            n_bins = 300
            )

        H = fitter.to_hist()
        H.mirror()

        H.set_value_for_patch(700.0,
            x_min=-0.1, x_max=1.0, y_min=-0.07, y_max=-0.05
            )        
        H.set_value_for_patch(700.0,
            x_min=-1.0, x_max=0.1, y_min=0.05, y_max=0.07
            )        

        H.add_padding(
            700.,
            x_min = -5.9,
            x_max = 2.6,
            y_min = -0.07,
            y_max = 0.07,
            )
        return H


class PlotPatcherKTCG_floatingBRs_scenario2(PlotPatcherKTCG_floatingBRs):
    """docstring for PlotPatcherKTCG_floatingBRs_scenario2"""

    def preprocess_combination(self):
        self.deltaNLL_threshold = -5.

    def patch_combination(self):
        H = super(PlotPatcherKTCG_floatingBRs_scenario2, self).patch_combination()

        H.set_value_for_patch(700.0,
            x_min=1.4, x_max=2.3, y_min=0.05, y_max=0.07
            )        
        H.set_value_for_patch(700.0,
            x_min=-2.3, x_max=-1.4, y_min=-0.07, y_max=-0.05
            )

        return H

    # def patch_hgg(self):
    #     self.shift_scan()
    #     return self.scan.to_hist()


@flag_as_option
def projection_ktcg_plot_floatingBRs(args):
    Patcher = PlotPatcherKTCG_floatingBRs_scenario2 if args.scenario2 else PlotPatcherKTCG_floatingBRs
    scans = []

    combWithHbb = Patcher(args, 'combWithHbb').patch()
    scans.append(combWithHbb)

    hzz = Patcher(args, 'hzz').patch()
    scans.append(hzz)

    hgg = Patcher(args, 'hgg').patch()
    scans.append(hgg)

    x_min, x_max, y_min, y_max = (None, None, None, None)
    # x_min, x_max, y_min, y_max = (-1.9, 1.9, -0.05, 0.05)
    x_min, x_max, y_min, y_max = (-5.9, 2.5, -0.050, 0.050)

    plot = differentials.plotting.plots.MultiContourPlot(
        ('projection_ktcg_plot_floatingBRs'
            + ('_asimov' if args.asimov else '')
            + ('_scenario2' if args.scenario2 else '')
            ),
        scans,
        x_title = differentials.core.standard_titles['ct'], y_title = differentials.core.standard_titles['cg'],
        x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max
        )

    plot.cdl.set(x2 = lambda c: c.GetLeftMargin() + 0.30, y2=lambda c: c.GetBottomMargin() + 0.17)

    plot.draw(wait=True)
    plot.add_BR_floating_text(
        x = lambda c: c.GetLeftMargin() + 0.03,
        y = lambda c: c.GetBottomMargin() + 0.24
        )
    plot.add_text_on_the_fly(
        text_size = 0.040,
        # x = lambda c: 1. - c.GetRightMargin() - 0.01,
        # y = lambda c: 1. - c.GetTopMargin() - 0.05,
        x = lambda c: c.GetLeftMargin() + 0.03,
        y = lambda c: c.GetBottomMargin() + 0.44,
        text = (
            '#splitline{w/ YR18}{syst. uncert. (S2)}'
            if args.scenario2 else
            '#splitline{w/ Run 2}{syst. uncert. (S1)}'
            )
        )
    plot.see_through_legends()
    plot.cdl.SetNColumns(2)
    plot.wrapup()

