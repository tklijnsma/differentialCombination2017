from OptionHandler import flag_as_option, flag_as_parser_options

import os, logging, copy, sys, re
from array import array
from math import sqrt, pi, cos, sin, tan, atan, atan2

import ROOT

import LatestPaths, LatestBinning
import differentials
import differentialutils
from math import sqrt

from differentials.plotting.canvas import c

import fermilabcode
from fermilabcode.read_canvas import CanvasReaderHzz, CanvasReaderHgg

#____________________________________________________________________
# CMS input processing

@flag_as_option
def read_hzz_from_canvas(args):
    cr = CanvasReaderHzz('fermilabcode/input/pT4l_unfoldwith_SM_125_logscale.root')
    cr.read()
    cr.dump('hzz_data')

@flag_as_option
def read_hgg_from_canvas(args):
    cr = CanvasReaderHgg('fermilabcode/input/diffXsec_pt_test_log.root')
    cr.read()
    cr.dump('hgg_data')

@flag_as_option
def naive_pt_combination(args):
    args = differentialutils.force_asimov(args)
    hgg = get_hgg(args)
    hzz = get_hzz(args)

    combination = fermilabcode.minicombine.PtCombination([hgg, hzz])
    combination.color = 1
    combination.plot_scans()

    hgg_hist = fermilabcode.minicombine.data_to_hist(hgg, title='hgg', color=2)
    hzz_hist = fermilabcode.minicombine.data_to_hist(hzz, title='hzz', color=4)

    combination_hist_hzz = combination.get_histogram_bestfit_hzz_binning()
    combination_hist_hzz.color = 3

    plot = differentials.plotting.plots.QuickPlot(
        'spectrum_ptcomb_withHzz',
        x_min=0., x_max=500., y_min=-1., y_max=3.,
        )
    plot.add(hgg_hist, 'repr_point_with_vertical_bar_and_horizontal_bar')
    plot.add(hzz_hist, 'repr_point_with_vertical_bar_and_horizontal_bar')
    plot.add(combination_hist_hzz, 'repr_point_with_horizontal_bar')
    plot.add(combination.to_hist(), 'repr_point_with_vertical_bar_and_horizontal_bar')
    plot.draw()
    plot.wrapup()

    combination.dump('combination_data' + ('_asimov' if args.asimov else ''))


def get_hgg(args):
    hgg = fermilabcode.minicombine.load_data('fermilabcode/input/hgg_data_Mar21.py')
    hgg.name = 'hgg'
    if args.asimov:
        hgg.mu = [ 1.0 for m in hgg.mu ]
    return hgg

def get_hzz(args):
    hzz = fermilabcode.minicombine.load_data('fermilabcode/input/hzz_data_Mar21.py')
    hzz.name = 'hzz'
    if args.asimov:
        hzz.mu = [ 1.0 for m in hzz.mu ]
    return hzz

def get_combination(args):
    if args.asimov:
        combination = fermilabcode.minicombine.load_data('fermilabcode/combination_data_asimov_Mar28.py')
    else:
        combination = fermilabcode.minicombine.load_data('fermilabcode/combination_data_Mar22.py')
    combination.name = 'combination'
    return combination

def get_combination_at_lumi(args, lumi):
    combination = fermilabcode.minicombine.load_data('fermilabcode/combination_data_Mar22.py')
    combination.name = 'combination'
    if not args.asimov: logging.warning('Forcing asimov, even though --asimov was not given')
    combination.mu = [ 1.0 for m in combination.mu ]
    combination.mu_up = [ d * 1/sqrt(lumi/35.9) for d in combination.mu_up ]
    combination.mu_down = [ d * 1/sqrt(lumi/35.9) for d in combination.mu_down ]
    return combination

#____________________________________________________________________
# ATLAS input processing

def get_hgg_ATLAS_3000fb():
    binning = [ 0., 20., 30., 45., 60., 80., 120., 170., 220., 280., 350., 600. ]
    deltas = [ 0.05, 0.058, 0.052, 0.06, 0.064, 0.05, 0.044, 0.05, 0.058, 0.07, 0.078 ]
    hgg = differentials.core.AttrDict(
        binning=binning,
        mu = [ 1.0 for i in xrange(len(binning)) ],
        mu_up = deltas,
        mu_down = deltas
        )
    hgg.deltas = deltas
    return hgg
    
def get_hgg_ATLAS():
    binning = [0., 20., 30., 45., 60., 80., 120., 170., 220., 350.]
    deltas = [0.31, 0.41, 0.42, 0.49, 0.73, 0.45, 0.41, 0.47, 0.46]
    hgg = differentials.core.AttrDict(
        binning=binning,
        mu = [ 1.0 for i in xrange(len(binning)-1) ],
        mu_up = deltas,
        mu_down = deltas
        )
    hgg.deltas = deltas
    return hgg

def get_combination_ATLAS():
    binning = [ 0., 10., 20., 30., 45., 60., 80., 120., 200., 350. ]
    deltas = [ 0.2755905511811024, 0.3228346456692913, 0.3858267716535433, 0.2992125984251969, 0.4409448818897638, 0.3858267716535433, 0.2204724409448819, 0.1889763779527559, 0.2440944881889764 ]
    combination = differentials.core.AttrDict(
        binning=binning,
        mu = [ 1.0 for i in xrange(len(binning)-1) ],
        mu_up = deltas,
        mu_down = deltas
        )
    combination.deltas = deltas
    return combination

def get_hgg_ATLAS_3000fb_mappings():
    hgg_deltas_36fb = get_hgg_ATLAS().deltas
    hgg_deltas_3000fb = get_hgg_ATLAS_3000fb().deltas[:-2] # pretends 280/350 boundary match... not a big deal
    decrease_factors = [ l3000/l36 for l36, l3000 in zip(hgg_deltas_36fb, hgg_deltas_3000fb) ]
    return decrease_factors

def get_combination_ATLAS_3000fb():
    combination = get_combination_ATLAS()
    decrease_factors_hgg = get_hgg_ATLAS_3000fb_mappings()
    # Since binning isnt the same, safer to take the average decrease and use that
    average_decrease_factor = sum(decrease_factors_hgg)/len(decrease_factors_hgg)
    combination.deltas  = [ i*average_decrease_factor for i in combination.deltas ]
    combination.mu_up   = [ i*average_decrease_factor for i in combination.mu_up ]
    combination.mu_down = [ i*average_decrease_factor for i in combination.mu_down ]
    return combination

#____________________________________________________________________
# Plots

@flag_as_option
def plot_naive_pt_combination(args):
    args = differentialutils.force_asimov(args)
    hgg = get_hgg(args)
    hzz = get_hzz(args)
    combination = get_combination(args)

    hgg_hist = fermilabcode.minicombine.data_to_hist(hgg, title=differentials.core.get_standard_title('hgg'), color=2)
    hzz_hist = fermilabcode.minicombine.data_to_hist(hzz, title=differentials.core.get_standard_title('hzz'), color=4)
    combination_hist = fermilabcode.minicombine.data_to_hist(combination, title=differentials.core.get_standard_title('combination'), color=1)

    hgg_hist_xs = fermilabcode.minicombine.data_to_hist(hgg, title=differentials.core.get_standard_title('hgg'), color=2, do_xs=True)
    hzz_hist_xs = fermilabcode.minicombine.data_to_hist(hzz, title=differentials.core.get_standard_title('hzz'), color=4, do_xs=True)
    combination_hist_xs = fermilabcode.minicombine.data_to_hist(combination, title=differentials.core.get_standard_title('combination'), color=1, do_xs=True)

    plot = differentials.plotting.plots.BottomPanelPlot('spectrum_ptcomb')
    plot.y_title_top = '#Delta#sigma(p_{T}^{H})/#Deltap_{T}^{H} (pb/GeV)'
    plot.y_title_bottom = 'ratio to SM'
    plot.x_title = 'p_{T}^{H} (GeV)'
    plot.disable_CMS_labels = True

    plot.top_x_min = 0.0
    plot.top_x_max = 500.0
    plot.bottom_x_min = 0.0
    plot.bottom_x_max = 500.0

    # plot.top_y_max = 10
    # plot.top_y_min = 0.00000001

    plot.add_bottom(hgg_hist, 'repr_point_with_vertical_bar_and_horizontal_bar')
    plot.add_bottom(hzz_hist, 'repr_point_with_vertical_bar_and_horizontal_bar')
    plot.add_bottom(combination_hist, 'repr_point_with_vertical_bar_and_horizontal_bar')

    plot.make_legend()
    plot.add_top(hgg_hist_xs, 'repr_point_with_vertical_bar_and_horizontal_bar', plot.leg)
    plot.add_top(hzz_hist_xs, 'repr_point_with_vertical_bar_and_horizontal_bar', plot.leg)
    plot.add_top(combination_hist_xs, 'repr_point_with_vertical_bar_and_horizontal_bar', plot.leg)
    plot.add_top(differentials.plotting.pywrappers.CMS_Latex_lumi(text_size=0.07), '')

    plot.draw()
    plot.wrapup()

@flag_as_option
def plot_naive_pt_combination_asimov(args):
    combination = fermilabcode.minicombine.load_data('fermilabcode/combination_data_Mar22.py')    
    combination.mu = [ 1.0 for mu in combination.mu ]
    combination_hist = fermilabcode.minicombine.data_to_hist(combination, title=differentials.core.get_standard_title('combination'), color=1)
    combination_hist_xs = fermilabcode.minicombine.data_to_hist(combination, title=differentials.core.get_standard_title('combination'), color=1, do_xs=True)

    plot = differentials.plotting.plots.BottomPanelPlot('spectrum_ptcomb_asimov')
    plot.y_title_top = '#Delta#sigma(p_{T}^{H})/#Deltap_{T}^{H} (pb/GeV)'
    plot.y_title_bottom = 'ratio to SM'
    plot.x_title = 'p_{T}^{H} (GeV)'
    plot.disable_CMS_labels = True

    plot.top_x_min = 0.0
    plot.top_x_max = 500.0
    plot.bottom_x_min = 0.0
    plot.bottom_x_max = 500.0
    plot.top_y_min = 0.001

    plot.make_legend()
    plot.add_bottom(combination_hist, 'repr_point_with_vertical_bar_and_horizontal_bar')
    plot.add_top(combination_hist_xs, 'repr_point_with_vertical_bar_and_horizontal_bar', plot.leg)

    plot.add_top(differentials.plotting.pywrappers.CMS_Latex_lumi(text_size=0.07), '')
    plot.draw()
    plot.wrapup()



@flag_as_option
def plot_naive_pt_combination_asimov_with3000fb(args):
    args = differentialutils.force_asimov(args)

    combination = get_combination(args)
    combination_3000fb = get_combination_at_lumi(args, 3000.)

    combination_hist = fermilabcode.minicombine.data_to_hist(combination, title=differentials.core.get_standard_title('combination'), color=1)
    combination_hist_xs = fermilabcode.minicombine.data_to_hist(combination, title=differentials.core.get_standard_title('combination'), color=1, do_xs=True)
    combination_hist_xs.title = 'CMS 35.9 fb^{-1}'

    combination_3000fb_hist = fermilabcode.minicombine.data_to_hist(combination_3000fb, title=differentials.core.get_standard_title('combination'), color=2)
    combination_3000fb_hist_xs = fermilabcode.minicombine.data_to_hist(combination_3000fb, title=differentials.core.get_standard_title('combination'), color=2, do_xs=True)
    combination_3000fb_hist_xs.title = 'CMS 3 ab^{-1}'

    plot = differentials.plotting.plots.BottomPanelPlot('spectrum_ptcomb_asimov_with_3000fb')
    plot.y_title_top = '#Delta#sigma(p_{T}^{H})/#Deltap_{T}^{H} (pb/GeV)'
    plot.y_title_bottom = 'ratio to SM'
    plot.x_title = 'p_{T}^{H} (GeV)'
    plot.disable_CMS_labels = True

    plot.top_x_min = 0.0
    plot.top_x_max = 500.0
    plot.bottom_x_min = 0.0
    plot.bottom_x_max = 500.0

    # plot.top_y_max = 10
    # plot.top_y_min = 0.00000001

    # plot.add_bottom(combination_hist, 'repr_point_with_vertical_bar_and_horizontal_bar')
    plot.add_bottom(combination_3000fb_hist, 'repr_point_with_vertical_bar_and_horizontal_bar')

    plot.make_legend()
    plot.add_top(combination_hist_xs, 'repr_point_with_vertical_bar_and_horizontal_bar', plot.leg)
    plot.add_top(combination_3000fb_hist_xs, 'repr_point_with_vertical_bar_and_horizontal_bar', plot.leg)

    plot.draw()
    plot.wrapup()


@flag_as_option
def OLD_plot_naive_pt_combination_asimov_with3000fb(args):
    args = differentialutils.force_asimov(args)

    combination = fermilabcode.minicombine.load_data('fermilabcode/combination_data_Mar22.py')    
    combination.mu = [ 1.0 for mu in combination.mu ]
    combination_hist = fermilabcode.minicombine.data_to_hist(combination, title=differentials.core.get_standard_title('combination'), color=1)
    combination_hist_xs = fermilabcode.minicombine.data_to_hist(combination, title=differentials.core.get_standard_title('combination'), color=1, do_xs=True)
    combination_hist_xs.title = '35.9 fb^{-1}'

    combination_3000fb = copy.deepcopy(combination)
    combination_3000fb.mu_up = [ d * 1/sqrt(3000./35.9) for d in combination_3000fb.mu_up ]
    combination_3000fb.mu_down = [ d * 1/sqrt(3000./35.9) for d in combination_3000fb.mu_down ]
    combination_3000fb_hist = fermilabcode.minicombine.data_to_hist(combination_3000fb, title=differentials.core.get_standard_title('combination'), color=2)
    combination_3000fb_hist_xs = fermilabcode.minicombine.data_to_hist(combination_3000fb, title=differentials.core.get_standard_title('combination'), color=2, do_xs=True)
    combination_3000fb_hist_xs.title = '3 ab^{-1}'

    plot = differentials.plotting.plots.BottomPanelPlot('spectrum_ptcomb_highlumi_asimov')
    plot.y_title_top = '#Delta#sigma(p_{T}^{H})/#Deltap_{T}^{H} (pb/GeV)'
    plot.y_title_bottom = 'ratio to SM'
    plot.x_title = 'p_{T}^{H} (GeV)'
    plot.disable_CMS_labels = True

    plot.top_x_min = 0.0
    plot.top_x_max = 500.0
    plot.bottom_x_min = 0.0
    plot.bottom_x_max = 500.0
    plot.top_y_min = 0.001

    plot.make_legend()
    # plot.add_bottom(combination_hist, 'repr_point_with_vertical_bar_and_horizontal_bar')
    plot.add_top(combination_hist_xs, 'repr_point_with_vertical_bar_and_horizontal_bar', plot.leg)
    plot.add_bottom(combination_3000fb_hist, 'repr_point_with_vertical_bar_and_horizontal_bar')
    plot.add_top(combination_3000fb_hist_xs, 'repr_point_with_vertical_bar_and_horizontal_bar', plot.leg)

    plot.add_top(differentials.plotting.pywrappers.CMS_Latex_lumi(text_size=0.07), '')
    plot.draw()
    plot.wrapup()


#____________________________________________________________________
# 2D plots

def quick_single_hist_plot(
        args, hist, name='kappabkappac',
        x_title = '#kappa_{c}', y_title = '#kappa_{b}',
        set_ranges_by_contour = True
        ):
    plot = differentials.plotting.plots.Single2DHistPlot(
        name + ('_asimov' if args.asimov else ''),
        hist,
        )
    plot.x_title = x_title
    plot.y_title = y_title
    plot.set_ranges_by_contour = set_ranges_by_contour
    plot.disable_CMS_labels = True
    plot.draw()
    return plot

def get_kappabkappac(data, c1_min=-35., c1_max=35., c2_min=-10., c2_max=10.):
    kappabkappac = fermilabcode.minicombine.CouplingFit()
    kappabkappac.c1_min = c1_min
    kappabkappac.c1_max = c1_max
    kappabkappac.c2_min = c2_min
    kappabkappac.c2_max = c2_max
    kappabkappac.title = '#kappa_{b}/#kappa_{c}'
    kappabkappac.color = 2
    kappabkappac.setup_kappabkappac(data)
    kappabkappac.get_scan()
    return kappabkappac


def get_kappatkappag(data, c1_min=-3.0, c1_max=3.0, c2_min=-0.20, c2_max=0.20):
    kappatkappag = fermilabcode.minicombine.CouplingFit()
    kappatkappag.c1_min = c1_min
    kappatkappag.c1_max = c1_max
    kappatkappag.c2_min = c2_min
    kappatkappag.c2_max = c2_max

    kappatkappag.c1_n_points = 300
    kappatkappag.c2_n_points = 300

    kappatkappag.title = '#kappa_{t}/#kappa_{g}'
    kappatkappag.color = 2
    kappatkappag.setup_kappatkappag(data)
    kappatkappag.get_scan()
    return kappatkappag


class KappabKappacPlotter(object):
    """docstring for KappabKappacPlotter"""
    def __init__(self, args):
        super(KappabKappacPlotter, self).__init__()
        args = differentialutils.force_asimov(args)
        self.set_run1(args)
        self.set_3ab(args)
    
    def set_run1(self, args):        
        self.CMS_combination = get_combination(args)
        self.CMS_kappabkappac = get_kappabkappac(self.CMS_combination)
        self.CMS_kappabkappac.title = 'CMS 35.9 fb^{-1}'
        self.CMS_kappabkappac.color = 2
        
        self.ATLAS_combination = get_combination_ATLAS()
        self.ATLAS_kappabkappac = get_kappabkappac(self.ATLAS_combination)
        self.ATLAS_kappabkappac.title = 'ATLAS 36.1 fb^{-1}'
        self.ATLAS_kappabkappac.color = 4

        self.summed_kappabkappac = self.CMS_kappabkappac.sum(self.ATLAS_kappabkappac)
        self.summed_kappabkappac.color = 8
        self.summed_kappabkappac.title = 'CMS+ATLAS'

    def set_3ab(self, args):
        self.CMS_combination_3000fb = get_combination_at_lumi(args, 3000.)
        self.CMS_kappabkappac_3000fb = get_kappabkappac(
            self.CMS_combination_3000fb,
            c1_min = -35. * 2./sqrt(3000./35.9),
            c1_max = 35.  * 2./sqrt(3000./35.9),
            c2_min = -10. * 2./sqrt(3000./35.9),
            c2_max = 10.  * 2./sqrt(3000./35.9)
            )
        self.CMS_kappabkappac_3000fb.title = 'CMS 3000 fb^{-1}'
        self.CMS_kappabkappac_3000fb.color = 46

        self.ATLAS_combination_3000fb = get_combination_ATLAS_3000fb()
        self.ATLAS_kappabkappac_3000fb = get_kappabkappac(
            self.ATLAS_combination_3000fb,
            c1_min = -35. * 2./sqrt(3000./35.9),
            c1_max = 35.  * 2./sqrt(3000./35.9),
            c2_min = -10. * 2./sqrt(3000./35.9),
            c2_max = 10.  * 2./sqrt(3000./35.9)
            )
        self.ATLAS_kappabkappac_3000fb.title = 'ATLAS 3000 fb^{-1}'
        self.ATLAS_kappabkappac_3000fb.color = 38

        self.summed_kappabkappac_3000fb = self.CMS_kappabkappac_3000fb.sum(self.ATLAS_kappabkappac_3000fb)
        self.summed_kappabkappac_3000fb.color = 419
        self.summed_kappabkappac_3000fb.title = 'CMS+ATLAS 3000 fb^{-1}'


    def CMS_only_plot(self, args):
        args = differentialutils.force_asimov(args)

        plot = differentials.plotting.plots.MultiContourPlot(
            'KappabKappacPlotter_CMSonly' + ('_asimov' if args.asimov else ''),
            [self.CMS_kappabkappac]
            )
        plot.x_title = '#kappa_{c}'
        plot.y_title = '#kappa_{b}'
        plot.disable_CMS_labels = True
        plot.draw_individual_contours = False
        plot.set_ranges_by_contour = True
        plot.draw()
        plot.wrapup()


    def CMS_ATLAS_plot(self, args):
        args = differentialutils.force_asimov(args)

        plot = differentials.plotting.plots.MultiContourPlot(
            'KappabKappacPlotter_CMSATLAS' + ('_asimov' if args.asimov else ''),
            [self.CMS_kappabkappac, self.ATLAS_kappabkappac, self.summed_kappabkappac]
            )
        plot.x_title = '#kappa_{c}'
        plot.y_title = '#kappa_{b}'
        plot.disable_CMS_labels = True
        plot.draw_individual_contours = False
        plot.set_ranges_by_contour = True
        plot.draw()
        plot.wrapup()


    def summed_only_plot(self, args):
        args = differentialutils.force_asimov(args)

        plot = differentials.plotting.plots.MultiContourPlot(
            'KappabKappacPlotter_summed_3000fb' + ('_asimov' if args.asimov else ''),
            [self.summed_kappabkappac, self.summed_kappabkappac_3000fb]
            )
        plot.x_title = '#kappa_{c}'
        plot.y_title = '#kappa_{b}'
        plot.disable_CMS_labels = True
        plot.draw_individual_contours = False
        plot.set_ranges_by_contour = True
        plot.draw()
        plot.wrapup()

@flag_as_option
def plots_kappabkappac(args):
    plotter = KappabKappacPlotter(args)
    plotter.CMS_only_plot(args)
    plotter.CMS_ATLAS_plot(args)
    plotter.summed_only_plot(args)


class KappatKappagPlotter(object):
    """docstring for KappatKappagPlotter"""
    def __init__(self, args):
        super(KappatKappagPlotter, self).__init__()
        args = differentialutils.force_asimov(args)
        self.set_run1(args)
        self.set_3ab(args)
    
    def set_run1(self, args):        
        self.CMS_combination = get_combination(args)
        self.CMS = get_kappatkappag(self.CMS_combination)
        self.CMS.title = 'CMS 35.9 fb^{-1}'
        self.CMS.color = 2
        
        self.ATLAS_combination = get_combination_ATLAS()
        self.ATLAS = get_kappatkappag(self.ATLAS_combination)
        self.ATLAS.title = 'ATLAS 36.1 fb^{-1}'
        self.ATLAS.color = 4

        self.summed = self.CMS.sum(self.ATLAS)
        self.summed.color = 8
        self.summed.title = 'CMS+ATLAS'

    def set_3ab(self, args):
        self.CMS_3000fb_combination = get_combination_at_lumi(args, 3000.)
        self.CMS_3000fb = get_kappatkappag(
            self.CMS_3000fb_combination,
            c1_min = 0.5,
            c1_max = 1.6,
            c2_min = -0.05,
            c2_max = 0.04,
            )
        self.CMS_3000fb.title = 'CMS 3000 fb^{-1}'
        self.CMS_3000fb.color = 46

        self.ATLAS_3000fb_combination = get_combination_ATLAS_3000fb()
        self.ATLAS_3000fb = get_kappatkappag(
            self.ATLAS_3000fb_combination,
            c1_min = 0.5,
            c1_max = 1.6,
            c2_min = -0.05,
            c2_max = 0.04,
            )
        self.ATLAS_3000fb.title = 'ATLAS 3000 fb^{-1}'
        self.ATLAS_3000fb.color = 38

        self.summed_3000fb = self.CMS_3000fb.sum(self.ATLAS_3000fb)
        self.summed_3000fb.color = 419
        self.summed_3000fb.title = 'CMS+ATLAS 3000 fb^{-1}'

    def CMS_ATLAS_plot(self, args):
        args = differentialutils.force_asimov(args)

        plot = differentials.plotting.plots.MultiContourPlot(
            'KappatKappagPlotter_CMSATLAS' + ('_asimov' if args.asimov else ''),
            [self.summed]
            )
        plot.x_title = '#kappa_{t}'
        plot.y_title = 'c_{g}'
        plot.disable_CMS_labels = True
        plot.draw_individual_contours = False
        plot.set_ranges_by_contour = False
        plot.draw()
        plot.wrapup()


    def CMS_ATLAS_3000fb_plot(self, args):
        plot = differentials.plotting.plots.MultiContourPlot(
            'KappatKappagPlotter_CMSATLAS_3000fb' + ('_asimov' if args.asimov else ''),
            [self.summed, self.summed_3000fb]
            )

        plot.legend.SetNColumns(1)
        plot.legend.set(
            x1 = c.GetLeftMargin() + 0.01,
            y1 = c.GetBottomMargin() + 0.02,
            x2 = c.GetLeftMargin() + 0.40,
            y2 = c.GetBottomMargin() + 0.26
            )

        plot.x_min = 0.5
        plot.x_max = 1.6
        plot.y_min = -0.05
        plot.y_max = 0.04
        plot.x_title = '#kappa_{t}'
        plot.y_title = 'c_{g}'
        plot.disable_CMS_labels = True
        plot.draw_individual_contours = False
        plot.set_ranges_by_contour = False
        plot.draw()
        plot.wrapup()

@flag_as_option
def plots_kappatkappag(args):
    plotter = KappatKappagPlotter(args)
    plotter.CMS_ATLAS_plot(args)
    plotter.CMS_ATLAS_3000fb_plot(args)



#____________________________________________________________________
# Same but for kappat kappag

@flag_as_option
def naive_kappatkappag_combination(args):
    args = differentialutils.force_asimov(args)

    CMS_combination = get_combination(args)
    CMS_kappatkappag = get_kappatkappag(CMS_combination)
    CMS_kappatkappag.title = 'CMS 35.9 fb^{-1}'
    CMS_kappatkappag.color = 2

    CMS_kappatkappag_hist = CMS_kappatkappag.to_hist()
    plot = quick_single_hist_plot(
        args, CMS_kappatkappag_hist, name='kappatkappag', x_title = '#kappa_{t}', y_title = 'c_{g}', set_ranges_by_contour=False
        ).wrapup()

    ATLAS_combination = get_combination_ATLAS()
    ATLAS_kappatkappag = get_kappatkappag(ATLAS_combination)
    ATLAS_kappatkappag.title = 'ATLAS 36.1 fb^{-1}'
    ATLAS_kappatkappag.color = 4

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_kappatkappag_CMS_ATLAS' + ('_asimov' if args.asimov else ''),
        [
            CMS_kappatkappag, ATLAS_kappatkappag,
            ]
        )
    plot.x_title = '#kappa_{t}'
    plot.y_title = 'c_{g}'
    plot.disable_CMS_labels = True
    plot.draw_individual_contours = False
    plot.set_ranges_by_contour = False
    plot.draw()
    plot.wrapup()


    CMS_3000fb_combination = get_combination_at_lumi(args, 3000.)
    CMS_3000fb_kappatkappag = get_kappatkappag(
        CMS_3000fb_combination,
        c1_min = 0.5,
        c1_max = 1.6,
        c2_min = -0.05,
        c2_max = 0.04,
        )
    CMS_3000fb_kappatkappag.title = 'CMS 3000 fb^{-1}'
    CMS_3000fb_kappatkappag.color = 46

    ATLAS_3000fb_combination = get_combination_ATLAS_3000fb()
    ATLAS_3000fb_kappatkappag = get_kappatkappag(
        ATLAS_3000fb_combination,
        c1_min = 0.5,
        c1_max = 1.6,
        c2_min = -0.05,
        c2_max = 0.04,
        )
    ATLAS_3000fb_kappatkappag.title = 'ATLAS 3000 fb^{-1}'
    ATLAS_3000fb_kappatkappag.color = 38

    # Add also the combination of both
    summed = CMS_3000fb_kappatkappag.sum(ATLAS_3000fb_kappatkappag)
    summed.title = 'CMS+ATLAS 3000 fb^{-1}'
    summed.color = 8

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_kappatkappag_CMS_ATLAS_with3000fb_withCombination' + ('_asimov' if args.asimov else ''),
        [
            CMS_kappatkappag, ATLAS_kappatkappag,
            CMS_3000fb_kappatkappag, ATLAS_3000fb_kappatkappag,
            summed
            ]
        )
    plot.legend.SetNColumns(1)
    plot.legend.set(
        x1 = c.GetLeftMargin() + 0.01,
        y1 = c.GetBottomMargin() + 0.02,
        x2 = c.GetLeftMargin() + 0.40,
        y2 = c.GetBottomMargin() + 0.44
        )

    plot.x_min = 0.5
    plot.x_max = 1.6
    plot.y_min = -0.05
    plot.y_max = 0.04
    plot.x_title = '#kappa_{t}'
    plot.y_title = 'c_{g}'
    plot.disable_CMS_labels = True
    plot.draw_individual_contours = False
    plot.set_ranges_by_contour = False
    plot.draw()
    plot.wrapup()


@flag_as_option
def thetaplot_kappatkappag(args):
    args = differentialutils.force_asimov(args)

    CMS_combination = get_combination(args)
    CMS_kappatkappag = get_kappatkappag(CMS_combination)
    CMS_kappatkappag.title = 'CMS 35.9 fb^{-1}'
    CMS_kappatkappag.color = 2

    CMS_kappatkappag_hist = CMS_kappatkappag.to_hist()
    plot = quick_single_hist_plot(
        args, CMS_kappatkappag_hist, name='kappatkappag', x_title = '#kappa_{t}', y_title = 'c_{g}', set_ranges_by_contour=False
        ).wrapup()

    # Get SM normalization
    binning = CMS_kappatkappag.chi2.binning()
    smxs = CMS_kappatkappag.chi2.evaluate_parametrization_xs(CMS_kappatkappag.bestfit.pois)
    inc_smxs = differentials.integral.Integrator(binning, smxs).integral(-10000., 800.)
    logging.info('Will normalize to SM xs: {0} pb'.format(inc_smxs))

    theta_1sigma_left = atan2(0.03, 0.6)
    theta_1sigma_right = atan2(-0.12, 2.4)
    thetas = [ 0.5*pi, theta_1sigma_left, 0.0, theta_1sigma_right, -0.5*pi ]

    theta_hists = []
    for theta in thetas:
        ct = cos(theta)
        cg = sin(theta)
        pois = [ ct, cg ]

        xss = CMS_kappatkappag.chi2.evaluate_parametrization_xs(pois)
        inc_xs = differentials.integral.Integrator(binning, xss).integral(-10000., 800.)

        # Normalize to the SM normalization
        xss = [ xs * (inc_smxs/inc_xs) for xs in xss ]

        # Now get the mu's
        mus = [ xs / xs_sm for xs, xs_sm in zip(xss, smxs) ]

        histogram = differentials.plotting.pywrappers.Histogram(
            differentials.plotting.plotting_utils.get_unique_rootname(),
            '#theta = {0:.2f} #pi'.format(theta/pi),
            binning,
            mus
            )
        # histogram.color = 2
        theta_hists.append(histogram)

    CMS_combination_hist = fermilabcode.minicombine.data_to_hist(CMS_combination, title='Data', color=1)


    plot = differentials.plotting.plots.QuickPlot(
        'spectra_theta',
        x_min=0., x_max=500., y_min=-1., y_max=5.,
        )
    plot.x_title = 'p_{T}^{H} (GeV)'
    plot.y_title = '#mu'
    for theta_hist in theta_hists:
        plot.add(theta_hist, 'repr_basic_histogram')
    plot.add(CMS_combination_hist, 'repr_point_with_vertical_bar_and_horizontal_bar')

    plot.draw()
    plot.wrapup()


@flag_as_option
def thetaplot_kappatkappag_ATLAS(args):
    args = differentialutils.force_asimov(args)

    ATLAS_combination = get_combination_ATLAS()
    ATLAS_kappatkappag = get_kappatkappag(ATLAS_combination)
    ATLAS_kappatkappag.title = 'ATLAS 36.1 fb^{-1}'
    ATLAS_kappatkappag.color = 2

    ATLAS_kappatkappag_hist = ATLAS_kappatkappag.to_hist()
    plot = quick_single_hist_plot(
        args, ATLAS_kappatkappag_hist, name='ATLAS_kappatkappag', x_title = '#kappa_{t}', y_title = 'c_{g}', set_ranges_by_contour=False
        ).wrapup()

    # Get SM normalization
    binning = ATLAS_kappatkappag.chi2.binning()
    smxs = ATLAS_kappatkappag.chi2.evaluate_parametrization_xs(ATLAS_kappatkappag.bestfit.pois)
    inc_smxs = differentials.integral.Integrator(binning, smxs).integral(-10000., 350.)
    logging.info('Will normalize to SM xs: {0} pb'.format(inc_smxs))

    # thetas = fermilabcode.minicombine.get_axis(-0.5*pi, 0.5*pi, 12)

    theta_1sigma_left = atan2(0.03, 0.6)
    theta_1sigma_right = atan2(-0.12, 2.4)
    thetas = [ 0.5*pi, theta_1sigma_left, 0.0, theta_1sigma_right, -0.5*pi ]

    theta_hists = []
    for theta in thetas:
        ct = cos(theta)
        cg = sin(theta)
        pois = [ ct, cg ]

        xss = ATLAS_kappatkappag.chi2.evaluate_parametrization_xs(pois)
        inc_xs = differentials.integral.Integrator(binning, xss).integral(-10000., 350.)

        # Normalize to the SM normalization
        xss = [ xs * (inc_smxs/inc_xs) for xs in xss ]

        # Now get the mu's
        mus = [ xs / xs_sm for xs, xs_sm in zip(xss, smxs) ]

        histogram = differentials.plotting.pywrappers.Histogram(
            differentials.plotting.plotting_utils.get_unique_rootname(),
            '#theta = {0:.2f} #pi'.format(theta/pi),
            binning,
            mus
            )
        # histogram.color = 2
        theta_hists.append(histogram)

    ATLAS_combination_hist = fermilabcode.minicombine.data_to_hist(ATLAS_combination, title='Data', color=1, x_max=350.)


    plot = differentials.plotting.plots.QuickPlot(
        'ATLAS_spectra_theta',
        x_min=0., x_max=350., y_min=-1., y_max=5.,
        )
    plot.x_title = 'p_{T}^{H} (GeV)'
    plot.y_title = '#mu'
    for theta_hist in theta_hists:
        plot.add(theta_hist, 'repr_basic_histogram')
    plot.add(ATLAS_combination_hist, 'repr_point_with_vertical_bar_and_horizontal_bar')

    plot.draw()
    plot.wrapup()

