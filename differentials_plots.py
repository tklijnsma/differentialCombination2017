#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import ROOT
import logging
import os, sys, re, copy
from os.path import *
from glob import glob
from copy import deepcopy
from math import sqrt
import traceback

from OptionHandler import flag_as_option

import LatestPaths
import LatestBinning

import differentials

from time import strftime
datestr = strftime('%b%d')


########################################
# Plotting
########################################

#____________________________________________________________________
@flag_as_option
def plot_all_differentials(args):
    pth_smH_plot(args)
    pth_ggH_plot(args)
    njets_plot(args)
    ptjet_plot(args)
    rapidity_plot(args)

@flag_as_option
def plot_pth(args):
    pth_smH_plot(args)
    pth_ggH_plot(args)


def get_sm_histograms(observable, normalize_by_second_to_last_bin_width, x_max=None):
    binning = observable.binning
    if not(x_max is None):
        if x_max < binning[-2]: raise ValueError('x_max {0} kills order of binning {1}'.format(x_max, binning))
        binning[-1] = x_max

    style = differentials.plotting.pywrappers.StyleSheet(color=16, plot_priority=-10, bin_center_offset=-0.34, fill_style=3257)

    # Create SM spectra
    sm_xs = differentials.plotting.pywrappers.Histogram(
        'auto', differentials.core.standard_titles['SM_Vittorio'],
        binning,
        observable.crosssection_over_binwidth(normalize_by_second_to_last_bin_width)
        )
    sm_xs.set_err_up(observable.unc_xs_over_binwidth(normalize_by_second_to_last_bin_width))
    sm_xs.set_err_down(observable.unc_xs_over_binwidth(normalize_by_second_to_last_bin_width))
    sm_xs.add_stylesheet(style)

    sm_ratio = differentials.plotting.pywrappers.Histogram(
        'auto', differentials.core.standard_titles['SM_Vittorio'],
        binning,
        [ 1. for i in observable.shape ]
        )
    sm_ratio.set_err_up(observable.unc_fraction)
    sm_ratio.set_err_down(observable.unc_fraction)
    sm_ratio.add_stylesheet(style)

    return sm_xs, sm_ratio


#____________________________________________________________________
ROOT.gStyle.SetEndErrorSize(3)
ROOT.gStyle.SetHatchesLineWidth(2)
style = differentials.plotting.pywrappers.StyleSheet()
style.line_width = 2
style.error_bar_line_width = 1

#____________________________________________________________________
@flag_as_option
def pth_smH_plot(args):
    spectra = []
    obs_name = 'pth_smH'
    obstuple = LatestBinning.obstuple_pth_smH
    scandict = LatestPaths.scan.pth_smH.asimov if args.asimov else LatestPaths.scan.pth_smH.observed

    APPLY_FIXED_BINNING = True

    # Load scans
    hgg = differentials.scans.DifferentialSpectrum('hgg', scandict.hgg)
    hgg.set_sm(obstuple.hgg.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    hgg.add_stylesheet(style.copy(color=differentials.core.safe_colors.red, marker_style=26))
    spectra.append(hgg)

    hzz = differentials.scans.DifferentialSpectrum('hzz', scandict.hzz)
    hzz.set_sm(obstuple.hzz.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    hzz.add_stylesheet(style.copy(color=differentials.core.safe_colors.blue, marker_style=32))
    spectra.append(hzz)

    hbb = differentials.scans.DifferentialSpectrum('hbb', scandict.hbb)
    hbb.drop_first_bin()
    hbb.set_sm(obstuple.hbb.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    hbb.add_stylesheet(style.copy(color=differentials.core.safe_colors.green, marker_style=27, marker_size=1.5))
    spectra.append(hbb)

    combWithHbb = differentials.scans.DifferentialSpectrum('combWithHbb', scandict.combWithHbb)
    combWithHbb.no_overflow_label = True
    combWithHbb.set_sm(obstuple.combWithHbb.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    combWithHbb.add_stylesheet(style.copy(color=1, plot_priority=20))
    spectra.append(combWithHbb)

    # Align the right boundary of all the spectra, but not at 10000
    x_max = max([ 2*s.binning()[-2]-s.binning()[-3] for s in spectra ])
    for s in spectra:
        s.read()
        if not args.table: s.give_x_max(x_max)
        s.draw_method = 'repr_vertical_bar_with_horizontal_lines_dashed_onlymerged'

    # Get syst only shape
    combWithHbb_statonly = differentials.scans.DifferentialSpectrum('combWithHbb_statonly', scandict.combWithHbb_statonly)
    combWithHbb_statonly.set_sm(obstuple.combWithHbb.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    combWithHbb_statonly.read()
    systshapemaker = differentials.systshapemaker.SystShapeMaker()
    systonly_histogram, systonly_histogram_xs = systshapemaker.get_systonly_histogram(combWithHbb, combWithHbb_statonly)

    # Get SM histograms
    sm_xs, sm_ratio = get_sm_histograms(obstuple.combWithHbb, normalize_by_second_to_last_bin_width=True, x_max=x_max)

    if args.table:
        rowproducer = differentials.plotting.tableproducer.SpectrumRowProducer(combWithHbb.binning())
        rowproducer.do_xs = True
        table = differentials.plotting.newtables.BaseTable()
        table.latex_mode(True)
        table.append(rowproducer.produce_binning_row('$\\pth$ (GeV)'))
        table.append(rowproducer.produce(hgg))
        table.append(rowproducer.produce(hzz))
        table.append(rowproducer.produce(hbb))
        table.append(rowproducer.produce(combWithHbb))
        print table.produce_table_string()
        table.produce_to_file('tables_{0}/pth_smH.tex'.format(differentials.core.datestr()))
        return

    # Start compiling plot
    plotname = 'spectra_{0}'.format(obs_name) + ('_asimov' if args.asimov else '')
    plot = differentials.plotting.plots.SpectraPlot(plotname, spectra)
    plot.draw_multiscans = True
    plot.obsname = obs_name
    plot.obsunit = 'GeV'

    # Add the SM and syst-only histograms
    if systshapemaker.success:
        plot.add_top(systonly_histogram_xs, systonly_histogram_xs.draw_method, plot.leg)
        plot.add_bottom(systonly_histogram, systonly_histogram.draw_method)
    plot.add_top(sm_xs, 'repr_basic_with_full_fill', plot.leg)
    plot.add_bottom(sm_ratio, 'repr_basic_with_full_fill')

    # Some ranges
    plot.top_y_min = 0.9*10e-6
    plot.top_y_max = 10.
    plot.bottom_y_min = -1.0
    plot.bottom_y_max = 5.0

    if APPLY_FIXED_BINNING:
        # Apply fixed binning
        reference_binning = combWithHbb.binning()
        plot.make_fixed_widths(reference_binning)
        plot.top_x_max = len(reference_binning)-1
        plot.bottom_x_max = len(reference_binning)-1
        hgg.style().bin_center_offset = -0.17
        hzz.style().bin_center_offset = 0.17
        hbb.style().bin_center_offset = 0.17
        hzz.style().plot_priority = 8
        plot.overflow_label_base_offset = 0.60
    else:
        hgg.style().bin_center_offset = -0.17
        hzz.style().bin_center_offset = 0.17
        hbb.style().bin_center_offset = 0.17
        plot.set_reference_bounds(combWithHbb.binning())
        plot.add_lines_at_bin_boundaries()
        plot.overflow_label_base_offset = 0.29

    lw = 0.42
    lh = 0.54 * 0.9

    if APPLY_FIXED_BINNING:
        xshift = 0.02
        yshift  = 0.09 + 0.01
        plot.leg.set(
            lambda c: c.GetLeftMargin() + xshift,
            lambda c: c.GetBottomMargin() + yshift,
            lambda c: c.GetLeftMargin() + xshift + lw,
            lambda c: c.GetBottomMargin() + yshift + lh,
            )
    else:
        lh = 0.54 * 0.6
        xshift = 0.24
        yshift = 0.46
        plot.leg.set(
            lambda c: c.GetLeftMargin() + xshift,
            lambda c: c.GetBottomMargin() + 0.09 + yshift,
            lambda c: c.GetLeftMargin() + xshift + 0.02 + lw,
            lambda c: c.GetBottomMargin() + 0.09 + lh + yshift,
            )

    plot.leg.SetNColumns(1)
    plot.draw()

    l = differentials.plotting.pywrappers.Latex(
        lambda c: c.GetLeftMargin() + 0.02 + xshift,
        lambda c: c.GetBottomMargin() - 0.04 + yshift - 0.009,
        '#sigma_{SM} from CYRM-2017-002'
        )
    if not APPLY_FIXED_BINNING:
        l.x = lambda c: c.GetLeftMargin() + 0.04 + xshift
        l.y = lambda c: c.GetBottomMargin() + 0.05 + yshift
    l.SetNDC()
    l.SetTextAlign(11)
    l.SetTextFont(42) 
    l.SetTextSize(0.038 + 0.005)
    l.Draw()

    l2 = differentials.plotting.pywrappers.Latex(
        lambda c: 1. - c.GetRightMargin() - 0.01,
        lambda c: plot.overflow_label_base_offset + 0.08*3,
        'Overflow norm.'
        )
    l2.SetNDC()
    l2.SetTextAlign(31)
    l2.SetTextFont(42) 
    l2.SetTextSize(0.04)
    l2.Draw()

    if APPLY_FIXED_BINNING: plot.replace_bin_labels([ '0', '15', '30', '45', '80', '120', '200', '350', '600', '#infty' ])
    plot.wrapup()

#____________________________________________________________________
@flag_as_option
def pth_ggH_plot(args):
    spectra = []
    obs_name = 'pth_ggH'
    obstuple = LatestBinning.obstuple_pth_ggH
    scandict = LatestPaths.scan.pth_ggH.asimov if args.asimov else LatestPaths.scan.pth_ggH.observed

    # Load scans
    combWithHbb = differentials.scans.DifferentialSpectrum('combWithHbb', scandict.combWithHbb)
    # combWithHbb.no_overflow_label = True
    combWithHbb.set_sm(obstuple.combWithHbb.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    combWithHbb.add_stylesheet(style.copy(color=1, plot_priority=20))
    spectra.append(combWithHbb)

    # Align the right boundary of all the spectra, but not at 10000
    x_max = max([ 2*s.binning()[-2]-s.binning()[-3] for s in spectra ])
    for s in spectra:
        s.read()
        s.give_x_max(x_max)
        s.draw_method = 'repr_vertical_bar_with_horizontal_lines_dashed_onlymerged'

    if args.table:
        # table = differentials.plotting.tables.SpectraTable('pth_ggH', [combWithHbb])
        # # table.print_only_symm_unc = True
        # logging.info('Table:\n{0}'.format( table.repr_terminal() ))
        # return
        rowproducer = differentials.plotting.tableproducer.SpectrumRowProducer(combWithHbb.binning(), last_bin_is_overflow=True)
        rowproducer.do_xs = True
        table = differentials.plotting.newtables.BaseTable()
        table.latex_mode(True)
        table.append(rowproducer.produce_binning_row('$\\pth$ (GeV)'))
        # table.append(rowproducer.produce(hgg))
        # table.append(rowproducer.produce(hzz))
        # table.append(rowproducer.produce(hbb))
        table.append(rowproducer.produce(combWithHbb))
        print table.produce_table_string()
        table.produce_to_file('tables_{0}/pth_ggH.tex'.format(differentials.core.datestr()))
        return

    # Get syst only shape
    combWithHbb_statonly = differentials.scans.DifferentialSpectrum('combWithHbb_statonly', scandict.combWithHbb_statonly)
    combWithHbb_statonly.set_sm(obstuple.combWithHbb.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    combWithHbb_statonly.read()
    systshapemaker = differentials.systshapemaker.SystShapeMaker()
    systonly_histogram, systonly_histogram_xs = systshapemaker.get_systonly_histogram(combWithHbb, combWithHbb_statonly)

    # Get SM histograms
    sm_xs, sm_ratio = get_sm_histograms(obstuple.combWithHbb, normalize_by_second_to_last_bin_width=True, x_max=x_max)

    # Start compiling plot
    plotname = 'spectra_{0}'.format(obs_name) + ('_asimov' if args.asimov else '')
    plot = differentials.plotting.plots.SpectraPlot(plotname, spectra)
    plot.draw_multiscans = True
    plot.obsname = obs_name
    plot.obsunit = 'GeV'
    plot.obstitle = differentials.core.get_standard_title('pth_smH')
    plot.y_title_top = '#Delta#sigma^{{ggH}}/#Delta{0} (pb/{1})'.format(plot.obstitle, plot.obsunit)
    plot.overflow_label_base_offset = 0.55

    # Add the SM and syst-only histograms
    if systshapemaker.success:
        plot.add_top(systonly_histogram_xs, systonly_histogram_xs.draw_method, plot.leg)
        plot.add_bottom(systonly_histogram, systonly_histogram.draw_method)
    plot.add_top(sm_xs, 'repr_basic_with_full_fill', plot.leg)
    plot.add_bottom(sm_ratio, 'repr_basic_with_full_fill')

    # Some ranges
    plot.top_y_min = 0.9*10e-6
    plot.top_y_max = 10.
    plot.bottom_y_min = -1.0
    plot.bottom_y_max = 4.0

    # Apply fixed binning
    reference_binning = combWithHbb.binning()
    plot.make_fixed_widths(reference_binning)
    plot.top_x_max = len(reference_binning)-1
    plot.bottom_x_max = len(reference_binning)-1
    plot.add_lines_at_bin_boundaries(range(1,len(reference_binning)-1))

    # lw = 0.38
    # lh = 0.41 * 0.5
    lw = 0.42
    lh = 0.54 * 0.5
    leg_yshift = 0.008
    plot.leg.set(
        lambda c: c.GetLeftMargin() + 0.02,
        lambda c: c.GetBottomMargin() + 0.09 + leg_yshift,
        lambda c: c.GetLeftMargin() + 0.02 + lw,
        lambda c: c.GetBottomMargin() + 0.09 + lh + leg_yshift,
        )
    plot.leg.SetNColumns(1)
    plot.draw()

    l = differentials.plotting.pywrappers.Latex(
        lambda c: c.GetLeftMargin() + 0.04,
        lambda c: c.GetBottomMargin() + 0.045 - 0.006 + leg_yshift,
        '#sigma_{SM} from CYRM-2017-002'
        )
    l.SetNDC()
    l.SetTextAlign(11)
    l.SetTextFont(42) 
    l.SetTextSize(0.038 + 0.005)
    l.Draw()

    lggh = differentials.plotting.pywrappers.Latex(
        lambda c: 1.0 - c.GetRightMargin() - 0.05,
        lambda c: 1.0 - c.GetTopMargin() - 0.05,
        'gg #rightarrow H'
        )
    lggh.SetNDC()
    lggh.SetTextSize(0.06)
    lggh.SetTextAlign(33)
    lggh.SetTextFont(42)
    lggh.Draw()

    l2 = differentials.plotting.pywrappers.Latex(
        lambda c: 1. - c.GetRightMargin() - 0.01,
        lambda c: plot.overflow_label_base_offset + 0.08,
        'Overflow norm.'
        )
    l2.SetNDC()
    l2.SetTextAlign(31)
    l2.SetTextFont(42) 
    l2.SetTextSize(0.04)
    l2.Draw()

    plot.replace_bin_labels([ '0', '15', '30', '45', '80', '120', '200', '350', '600', '#infty' ])
    plot.wrapup()
    return plot

#____________________________________________________________________
@flag_as_option
def njets_plot(args):
    spectra = []
    obs_name = 'njets'
    obstuple = LatestBinning.obstuple_njets
    scandict = LatestPaths.scan.njets.asimov if args.asimov else LatestPaths.scan.njets.observed

    hgg = differentials.scans.DifferentialSpectrum('hgg', scandict.hgg)
    hgg.set_sm(obstuple.hgg.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    hgg.add_stylesheet(style.copy(color=differentials.core.safe_colors.red, marker_style=26))
    spectra.append(hgg)

    hzz = differentials.scans.DifferentialSpectrum('hzz', scandict.hzz)
    hzz.set_sm(obstuple.hzz.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    hzz.add_stylesheet(style.copy(color=differentials.core.safe_colors.blue, marker_style=32))
    spectra.append(hzz)

    combination = differentials.scans.DifferentialSpectrum('combination', scandict.combination)
    combination.no_overflow_label = True
    combination.set_sm(obstuple.combination.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    combination.add_stylesheet(style.copy(color=1, plot_priority=20))
    spectra.append(combination)

    # Align the right boundary of all the spectra, but not at 10000
    x_max = max([ 2*s.binning()[-2]-s.binning()[-3] for s in spectra ])
    for s in spectra:
        s.no_overflow_label = True
        s.read()
        s.give_x_max(x_max)
        s.draw_method = 'repr_vertical_bar_with_horizontal_lines_dashed_onlymerged'

    if args.table:
        # table = differentials.plotting.tables.SpectraTable('njets',
        #     spectra
        #     )
        # table.print_only_symm_unc = True
        # table.add_symm_improvement_row(hgg, combination)
        # # table.add_symm_improvement_row(hgg, combination)
        # # table.add_symm_improvement_row(combination, combWithHbb)
        # logging.info('Table:\n{0}'.format( table.repr_terminal() ))
        # return
        rowproducer = differentials.plotting.tableproducer.SpectrumRowProducer(combination.binning(), last_bin_is_overflow=True)
        rowproducer.do_xs = True
        table = differentials.plotting.newtables.BaseTable()
        table.latex_mode(True)
        table.append(rowproducer.produce_row_given_labels([
            '$\\njets$', '0', '1', '2', '3', '$\\ge$ 4'
            ]))
        table.append(rowproducer.produce(hgg))
        table.append(rowproducer.produce(hzz))
        table.append(rowproducer.produce(combination))
        print table.produce_table_string()
        table.produce_to_file('tables_{0}/njets.tex'.format(differentials.core.datestr()))
        return


    # Get syst only shape
    combination_statonly = differentials.scans.DifferentialSpectrum('combination_statonly', scandict.combination_statonly)
    combination_statonly.set_sm(obstuple.combination.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    combination_statonly.read()
    systshapemaker = differentials.systshapemaker.SystShapeMaker()
    systonly_histogram, systonly_histogram_xs = systshapemaker.get_systonly_histogram(combination, combination_statonly)

    # Get SM histograms
    sm_xs, sm_ratio = get_sm_histograms(obstuple.combination, normalize_by_second_to_last_bin_width=True, x_max=x_max)

    # Start compiling plot
    plotname = 'spectra_{0}'.format(obs_name) + ('_asimov' if args.asimov else '')
    plot = differentials.plotting.plots.SpectraPlot(plotname, spectra)
    plot.draw_multiscans = True
    plot.obsname = obs_name
    # plot.obsunit = 'GeV'
    plot.overflow_label_base_offset = 0.35

    # Add the SM and syst-only histograms
    if systshapemaker.success:
        plot.add_top(systonly_histogram_xs, systonly_histogram_xs.draw_method, plot.leg)
        plot.add_bottom(systonly_histogram, systonly_histogram.draw_method)
    plot.add_top(sm_xs, 'repr_basic_with_full_fill', plot.leg)
    plot.add_bottom(sm_ratio, 'repr_basic_with_full_fill')

    # Some ranges
    plot.top_y_min = 2*10e-2
    # plot.top_y_max = 500.
    plot.top_y_max = 1500.
    plot.bottom_y_min = -0.9
    plot.bottom_y_max = 4.0

    # Apply fixed binning
    reference_binning = combination.binning()
    plot.make_fixed_widths(reference_binning)
    plot.top_x_max = len(reference_binning)-1
    plot.bottom_x_max = len(reference_binning)-1
    hgg.style().bin_center_offset = -0.17
    hzz.style().bin_center_offset = 0.17
    hzz.style().plot_priority = 8
    plot.add_lines_at_bin_boundaries(range(1,len(reference_binning)-1))

    # leg_xshift = -0.115
    leg_xshift = 0.02
    leg_yshift = -0.00
    # legw = 0.38
    # legh = 0.41 * 5./6.
    legw = 0.42
    legh = 0.54 * 4./6.
    plot.leg.set(
        lambda c: 1-c.GetRightMargin() + leg_xshift - legw,
        lambda c: 1-c.GetTopMargin() + leg_yshift - legh,
        lambda c: 1-c.GetRightMargin() + leg_xshift,
        lambda c: 1-c.GetTopMargin() + leg_yshift,
        )
    plot.leg.SetNColumns(1)
    plot.draw()

    l = differentials.plotting.pywrappers.Latex(
        lambda c: 1-c.GetRightMargin() + leg_xshift + 0.02 - legw,        
        lambda c: 1-c.GetTopMargin()   + leg_yshift - legh - 0.041 - 0.009,
        '#sigma_{SM} from CYRM-2017-002'
        )
    l.SetNDC()
    l.SetTextAlign(11)
    l.SetTextFont(42) 
    l.SetTextSize(0.035 + 0.005)
    l.Draw()

    plot.replace_bin_labels([ '0', '1', '2', '3', '#geq4' ], offset=0.1)
    plot.base_bottom.GetXaxis().SetNdivisions(500)
    plot.base_top.GetXaxis().SetNdivisions(500)
    plot.wrapup()


#____________________________________________________________________
@flag_as_option
def ptjet_plot(args):
    spectra = []
    obs_name = 'ptjet'
    obstuple = LatestBinning.obstuple_ptjet
    scandict = LatestPaths.scan.ptjet.asimov if args.asimov else LatestPaths.scan.ptjet.observed

    hgg = differentials.scans.DifferentialSpectrum('hgg', scandict.hgg)
    hgg.set_sm(obstuple.hgg.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    hgg.add_stylesheet(style.copy(color=differentials.core.safe_colors.red, marker_style=26))
    spectra.append(hgg)

    hzz = differentials.scans.DifferentialSpectrum('hzz', scandict.hzz)
    hzz.set_sm(obstuple.hzz.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    hzz.add_stylesheet(style.copy(color=differentials.core.safe_colors.blue, marker_style=32))
    spectra.append(hzz)

    combination = differentials.scans.DifferentialSpectrum('combination', scandict.combination)
    combination.no_overflow_label = True
    combination.set_sm(obstuple.combination.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    combination.add_stylesheet(style.copy(color=1, plot_priority=20))
    spectra.append(combination)

    # Align the right boundary of all the spectra, but not at 10000
    x_max = max([ 2*s.binning()[-2]-s.binning()[-3] for s in spectra ])
    for s in spectra:
        s.read()
        s.drop_first_bin()
        s.give_x_max(x_max)
        s.draw_method = 'repr_vertical_bar_with_horizontal_lines_dashed_onlymerged'

    if args.table:
        rowproducer = differentials.plotting.tableproducer.SpectrumRowProducer(combination.binning(), last_bin_is_overflow=True)
        rowproducer.do_xs = True
        table = differentials.plotting.newtables.BaseTable()
        table.latex_mode(True)
        table.append(rowproducer.produce_binning_row('$\\ptjet$ (GeV)'))
        table.append(rowproducer.produce(hgg))
        table.append(rowproducer.produce(hzz))
        table.append(rowproducer.produce(combination))
        print table.produce_table_string()
        table.produce_to_file('tables_{0}/ptjet.tex'.format(differentials.core.datestr()))
        return

    # Get syst only shape
    combination_statonly = differentials.scans.DifferentialSpectrum('combination_statonly', scandict.combination_statonly)
    combination_statonly.set_sm(obstuple.combination.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    combination_statonly.read()
    combination_statonly.drop_first_bin()
    systshapemaker = differentials.systshapemaker.SystShapeMaker()
    systonly_histogram, systonly_histogram_xs = systshapemaker.get_systonly_histogram(combination, combination_statonly)

    # Get SM histograms
    obs = copy.deepcopy(obstuple.combination)
    obs.drop_first_bin()
    sm_xs, sm_ratio = get_sm_histograms(obs, normalize_by_second_to_last_bin_width=True, x_max=x_max)

    # Start compiling plot
    plotname = 'spectra_{0}'.format(obs_name) + ('_asimov' if args.asimov else '')
    plot = differentials.plotting.plots.SpectraPlot(plotname, spectra)
    plot.draw_multiscans = True
    plot.obsname = obs_name
    plot.obsunit = 'GeV'
    plot.overflow_label_base_offset = 0.33

    # Add the SM and syst-only histograms
    if systshapemaker.success:
        plot.add_top(systonly_histogram_xs, systonly_histogram_xs.draw_method, plot.leg)
        plot.add_bottom(systonly_histogram, systonly_histogram.draw_method)
    plot.add_top(sm_xs, 'repr_basic_with_full_fill', plot.leg)
    plot.add_bottom(sm_ratio, 'repr_basic_with_full_fill')

    # Some ranges
    plot.top_y_min = 0.8*10e-3
    plot.top_y_max = 50.
    plot.bottom_y_min = -0.8
    plot.bottom_y_max = 3.8

    # Apply fixed binning
    reference_binning = combination.binning()
    # reference_binning[0] = -1000. # Underflow left is somehow defined differently
    plot.make_fixed_widths(reference_binning)
    plot.top_x_min = 0.
    plot.bottom_x_min = 0.
    plot.top_x_max = len(reference_binning)-1
    plot.bottom_x_max = len(reference_binning)-1
    hgg.style().bin_center_offset = -0.17
    hzz.style().bin_center_offset = 0.17
    hzz.style().plot_priority = 8
    plot.add_lines_at_bin_boundaries(range(1,len(reference_binning)-1))

    # leg_xshift = -0.115
    leg_xshift = 0.02
    leg_yshift = -0.00
    # legw = 0.38
    # legh = 0.41 * 5./6.
    legw = 0.42
    legh = 0.54 * 4./6.
    # plot.leg.set(
    #     lambda c: 1-c.GetRightMargin() + leg_xshift - legw,
    #     lambda c: 1-c.GetTopMargin() + leg_yshift - legh,
    #     lambda c: 1-c.GetRightMargin() + leg_xshift,
    #     lambda c: 1-c.GetTopMargin() + leg_yshift,
    #     )
    plot.leg.set( # Move to left to accommodate more text for overflow
        lambda c: c.GetLeftMargin() + leg_xshift,
        lambda c: 1-c.GetTopMargin() + leg_yshift - legh,
        lambda c: c.GetLeftMargin() + leg_xshift + legw,
        lambda c: 1-c.GetTopMargin() + leg_yshift,
        )
    plot.leg.SetNColumns(1)
    plot.draw()

    l = differentials.plotting.pywrappers.Latex(
        # lambda c: 1-c.GetRightMargin() + leg_xshift + 0.02 - legw,        
        # lambda c: 1-c.GetTopMargin()   + leg_yshift - legh - 0.041 - 0.009,
        lambda c: c.GetLeftMargin() + leg_xshift + 0.02,
        lambda c: 1-c.GetTopMargin() + leg_yshift - legh - 0.041 - 0.009,
        '#sigma_{SM} from CYRM-2017-002'
        )
    l.SetNDC()
    l.SetTextAlign(11)
    l.SetTextFont(42) 
    l.SetTextSize(0.035 + 0.005)
    l.Draw()

    l2 = differentials.plotting.pywrappers.Latex(
        lambda c: 1. - c.GetRightMargin() - 0.01,
        lambda c: plot.overflow_label_base_offset + 0.08*2,
        'Overflow norm.'
        )
    l2.SetNDC()
    l2.SetTextAlign(31)
    l2.SetTextFont(42) 
    l2.SetTextSize(0.04)
    l2.Draw()

    plot.replace_bin_labels([ '30', '55', '95', '120', '200', '#infty' ])
    plot.wrapup()


#____________________________________________________________________
@flag_as_option
def rapidity_plot(args):
    spectra = []
    obs_name = 'rapidity'
    obstuple = LatestBinning.obstuple_rapidity
    scandict = LatestPaths.scan.rapidity.asimov if args.asimov else LatestPaths.scan.rapidity.observed

    hgg = differentials.scans.DifferentialSpectrum('hgg', scandict.hgg)
    hgg.set_sm(obstuple.hgg.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=False))
    hgg.add_stylesheet(style.copy(color=differentials.core.safe_colors.red, marker_style=26))
    spectra.append(hgg)

    hzz = differentials.scans.DifferentialSpectrum('hzz', scandict.hzz)
    hzz.set_sm(obstuple.hzz.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=False))
    hzz.add_stylesheet(style.copy(color=differentials.core.safe_colors.blue, marker_style=32))
    spectra.append(hzz)

    combination = differentials.scans.DifferentialSpectrum('combination', scandict.combination)
    combination.set_sm(obstuple.combination.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=False))
    combination.add_stylesheet(style.copy(color=1, plot_priority=20))
    spectra.append(combination)

    # Align the right boundary of all the spectra, but not at 10000
    x_max = 2.5
    for s in spectra:
        s.no_overflow_label = True
        s.read()
        s.give_x_max(x_max)
        s.draw_method = 'repr_vertical_bar_with_horizontal_lines_dashed_onlymerged'

    if args.table:
        rowproducer = differentials.plotting.tableproducer.SpectrumRowProducer(combination.binning(), last_bin_is_overflow=False)
        rowproducer.do_xs = True
        table = differentials.plotting.newtables.BaseTable()
        table.latex_mode(True)
        table.append(rowproducer.produce_binning_row('$\\absy$'))
        table.append(rowproducer.produce(hgg))
        table.append(rowproducer.produce(hzz))
        table.append(rowproducer.produce(combination))
        print table.produce_table_string()
        table.produce_to_file('tables_{0}/absy.tex'.format(differentials.core.datestr()))
        return

    # Get syst only shape
    combination_statonly = differentials.scans.DifferentialSpectrum('combination_statonly', scandict.combination_statonly)
    combination_statonly.set_sm(obstuple.combination.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=False))
    combination_statonly.read()
    systshapemaker = differentials.systshapemaker.SystShapeMaker()
    systonly_histogram, systonly_histogram_xs = systshapemaker.get_systonly_histogram(combination, combination_statonly)

    # Get SM histograms
    sm_xs, sm_ratio = get_sm_histograms(obstuple.combination, normalize_by_second_to_last_bin_width=False, x_max=x_max)

    # Start compiling plot
    plotname = 'spectra_{0}'.format(obs_name) + ('_asimov' if args.asimov else '')
    plot = differentials.plotting.plots.SpectraPlot(plotname, spectra)
    plot.draw_multiscans = True
    plot.obsname = obs_name
    # plot.obsunit = '#Delta|y_{H}|'
    plot.overflow_label_base_offset = 0.35

    # Add the SM and syst-only histograms
    if systshapemaker.success:
        plot.add_top(systonly_histogram_xs, systonly_histogram_xs.draw_method, plot.leg)
        plot.add_bottom(systonly_histogram, systonly_histogram.draw_method)
    plot.add_top(sm_xs, 'repr_basic_with_full_fill', plot.leg)
    plot.add_bottom(sm_ratio, 'repr_basic_with_full_fill')

    # Some ranges
    plot.top_y_min = 5.0
    plot.top_y_max = 900.
    plot.bottom_y_min = 0.4
    plot.bottom_y_max = 2.2

    # Apply fixed binning
    reference_binning = combination.binning()
    plot.make_fixed_widths(reference_binning)
    plot.top_x_max = len(reference_binning)-1
    plot.bottom_x_max = len(reference_binning)-1
    hgg.style().bin_center_offset = -0.17
    hzz.style().bin_center_offset = 0.17
    hzz.style().plot_priority = 8
    plot.add_lines_at_bin_boundaries(range(1,len(reference_binning)-1))

    leg_xshift = 0.02
    leg_yshift = -0.00
    # legw = 0.38
    # legh = 0.41 * 5./6.
    legw = 0.42
    legh = 0.54 * 4./6.
    plot.leg.set(
        lambda c: c.GetLeftMargin() + leg_xshift,
        lambda c: 1-c.GetTopMargin() + leg_yshift - legh,
        lambda c: c.GetLeftMargin() + leg_xshift + legw,
        lambda c: 1-c.GetTopMargin() + leg_yshift,
        )
    plot.leg.SetNColumns(1)
    plot.draw()

    plot.base_bottom.GetYaxis().SetNdivisions(505)

    l = differentials.plotting.pywrappers.Latex(
        lambda c: c.GetLeftMargin() + leg_xshift + 0.02,
        lambda c: 1-c.GetTopMargin() + leg_yshift - legh - 0.041 - 0.009,
        '#sigma_{SM} from CYRM-2017-002'
        )
    l.SetNDC()
    l.SetTextAlign(11)
    l.SetTextFont(42) 
    l.SetTextSize(0.035 + 0.005)
    l.Draw()

    plot.replace_bin_labels([ '0.0', '0.15', '0.3', '0.6', '0.9', '1.2', '2.5' ])
    plot.wrapup()


########################################
# Other plots
########################################
#____________________________________________________________________
@flag_as_option
def all_tables(args):
    pth_smH_tables(args)
    pth_ggH_tables(args)
    ptjet_tables(args)
    njets_tables(args)
    rapidity_tables(args)

@flag_as_option
def pth_smH_tables(args):
    differentialTable = DifferentialTable.DifferentialTable(name='pth_smH', last_bin_is_overflow=True)
    for decay_channel in ['hgg', 'hzz', 'combination']:
        statsyst = read_container(args, 'pth_smH', pth_smH_obs(decay_channel).crosssection_over_binwidth(), decay_channel, statonly=False )
        statonly = read_container(args, 'pth_smH', pth_smH_obs(decay_channel).crosssection_over_binwidth(), decay_channel, statonly=True )
        differentialTable.calculate_stat_syst(statsyst.name, statsyst, statonly)
    print differentialTable.repr_twiki_symm()
    print
    # differentialTable.do_xs = True
    # print differentialTable.repr_twiki_symm()

@flag_as_option
def pth_ggH_tables(args):
    differentialTable = DifferentialTable.DifferentialTable(name='pth_ggH', last_bin_is_overflow=True)
    for decay_channel in ['hgg', 'hzz', 'combination']:
        statsyst = read_container(args, 'pth_ggH', pth_ggH_obs(decay_channel).crosssection_over_binwidth(), decay_channel, statonly=False )
        statonly = read_container(args, 'pth_ggH', pth_ggH_obs(decay_channel).crosssection_over_binwidth(), decay_channel, statonly=True )
        differentialTable.calculate_stat_syst(statsyst.name, statsyst, statonly)
    print differentialTable.repr_twiki_symm()
    print
    # differentialTable.do_xs = True
    # print differentialTable.repr_twiki_symm()

@flag_as_option
def ptjet_tables(args):
    differentialTable = DifferentialTable.DifferentialTable(name='ptjet', last_bin_is_overflow=True)
    for decay_channel in ['hgg', 'hzz', 'combination']:
        statsyst = read_container(args, 'ptjet', ptjet_obs(decay_channel).crosssection_over_binwidth(), decay_channel, statonly=False )
        statonly = read_container(args, 'ptjet', ptjet_obs(decay_channel).crosssection_over_binwidth(), decay_channel, statonly=True )
        differentialTable.calculate_stat_syst(statsyst.name, statsyst, statonly)
    print differentialTable.repr_twiki_symm()
    print
    # differentialTable.do_xs = True
    # print differentialTable.repr_twiki_symm()

@flag_as_option
def njets_tables(args):
    differentialTable = DifferentialTable.DifferentialTable(name='njets', last_bin_is_overflow=True)
    for decay_channel in ['hgg', 'hzz', 'combination']:
        statsyst = read_container(args, 'njets', njets_obs(decay_channel).crosssection_over_binwidth(), decay_channel, statonly=False )
        statonly = read_container(args, 'njets', njets_obs(decay_channel).crosssection_over_binwidth(), decay_channel, statonly=True )
        differentialTable.calculate_stat_syst(statsyst.name, statsyst, statonly)
    print differentialTable.repr_twiki_symm()
    print
    # differentialTable.do_xs = True
    # print differentialTable.repr_twiki_symm()

@flag_as_option
def rapidity_tables(args):
    differentialTable = DifferentialTable.DifferentialTable(name='rapidity', last_bin_is_overflow=False)
    for decay_channel in ['hgg', 'hzz', 'combination']:
        statsyst = read_container(args, 'rapidity', rapidity_obs(decay_channel).crosssection_over_binwidth(), decay_channel, statonly=False )
        statonly = read_container(args, 'rapidity', rapidity_obs(decay_channel).crosssection_over_binwidth(), decay_channel, statonly=True )
        differentialTable.calculate_stat_syst(statsyst.name, statsyst, statonly)
    print differentialTable.repr_twiki_symm()
    print
    # differentialTable.do_xs = True
    # print differentialTable.repr_twiki_symm()
