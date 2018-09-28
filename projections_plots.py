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
p3000 = differentials.core.AttrDict()
p3000.hgg = 'out/ScanProjection_Jun28_pth_smH_hgg_asimov_0'
p3000.hzz = 'out/ScanProjection_Jun28_pth_smH_hzz_asimov_1'
# p3000.hbb = 'out/ScanProjection_Jun28_pth_smH_hbb_asimov_0'
# p3000.hbb = 'out/ScanProjection_Jul02_pth_smH_hbb_asimov_0'
# p3000.hbb = 'out/ScanProjection_Jul03_pth_smH_hbb_300ifb_asimov'
p3000.hbb = 'out/ScanProjection_Jul03_pth_smH_hbb_3000ifb_asimov'
# p3000.hbb = 'out/ScanProjection_Jul03_pth_smH_hbb_3000ifb_asimov_0'
# p3000.combWithHbb = 'out/ScanProjection_Jun29_pth_smH_combWithHbb_asimov'
p3000.combWithHbb = 'out/ScanProjection_Jul03_pth_smH_combWithHbb_3000ifb_asimov'

p3000.hgg = 'out/ScanProjection_Jul03_pth_smH_hgg_3000ifb_asimov'
p3000.hzz = 'out/ScanProjection_Jul03_pth_smH_hzz_3000ifb_asimov'
p3000.hbb = 'out/ScanProjection_Jul03_pth_smH_hbb_3000ifb_asimov_1'
p3000.combWithHbb = 'out/ScanProjection_Jul03_pth_smH_combWithHbb_3000ifb_asimov_0'

p36 = differentials.core.AttrDict()
p36.hgg = 'out/ScanProjection_Jul04_pth_smH_hgg_36ifb_asimov'
p36.hzz = 'out/ScanProjection_Jul04_pth_smH_hzz_36ifb_asimov'
p36.hbb = 'out/ScanProjection_Jul04_pth_smH_hbb_36ifb_asimov'
p36.combWithHbb = 'out/ScanProjection_Jul04_pth_smH_combWithHbb_36ifb_asimov'

p3000_statonly = differentials.core.AttrDict()
# p3000_statonly.hgg = 'out/ScanProjection_Jun29_pth_smH_hgg_statonly_asimov'
# p3000_statonly.hzz = 'out/ScanProjection_Jun29_pth_smH_hzz_statonly_asimov'
# p3000_statonly.hbb = 'out/ScanProjection_Jun29_pth_smH_hbb_statonly_asimov'
p3000_statonly.hgg = 'out/ScanProjection_Jul04_pth_smH_hgg_3000ifb_statonly_asimov'
p3000_statonly.hzz = 'out/ScanProjection_Jul04_pth_smH_hzz_3000ifb_statonly_asimov'
p3000_statonly.hbb = 'out/ScanProjection_Jul04_pth_smH_hbb_3000ifb_statonly_asimov'
p3000_statonly.combWithHbb = 'out/ScanProjection_Jul04_pth_smH_combWithHbb_3000ifb_statonly_asimov'

p3000_s2 = differentials.core.AttrDict()
p3000_s2.hgg = 'out/ScanProjection_Jul17_pth_smH_hgg_3000ifb_scenario2_asimov'
p3000_s2.hzz = 'out/ScanProjection_Jul17_pth_smH_hzz_3000ifb_scenario2_asimov'
p3000_s2.hbb = 'out/ScanProjection_Jul17_pth_smH_hbb_3000ifb_scenario2_asimov'
p3000_s2.combWithHbb = 'out/ScanProjection_Jul17_pth_smH_combWithHbb_3000ifb_scenario2_asimov'

#____________________________________________________________________
@flag_as_option
def projection_pth_smH_plot(args):
    differentials.scans.Scan.deltaNLL_threshold = -10.
    differentials.plotting.pywrappers.CMS_Latex_lumi.CMS_lumi = 3000.
    spectra = []
    obs_name = 'pth_smH'
    obstuple = LatestBinning.obstuple_pth_smH
    scandict = p3000 

    # APPLY_FIXED_BINNING = False
    APPLY_FIXED_BINNING = True

    PLOT_SYSTEMATIC_ONLY = False

    DO_STAT_ONLY = True
    # DO_STAT_ONLY = False

    # DO_S2 = True
    DO_S2 = False

    # DO_LUMI_36 = True
    DO_LUMI_36 = False


    if DO_STAT_ONLY:
        scandict = p3000_statonly

    if DO_S2:
        scandict = p3000_s2

    if DO_LUMI_36:
        scandict = p36


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
    # hbb.drop_first_bin()
    hbb.set_sm(obstuple.hbb.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    hbb.add_stylesheet(style.copy(color=differentials.core.safe_colors.green, marker_style=27, marker_size=1.5))
    spectra.append(hbb)

    # hbb.read()
    # filter_hbb(hbb)
    # sys.exit()

    # combWithHbb_given = False
    combWithHbb = differentials.scans.DifferentialSpectrum('combWithHbb', scandict.combWithHbb)
    combWithHbb.no_overflow_label = True
    combWithHbb.set_sm(obstuple.combWithHbb.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    combWithHbb.add_stylesheet(style.copy(color=1, plot_priority=20))
    spectra.append(combWithHbb)
    combWithHbb_given = True

    # Align the right boundary of all the spectra, but not at 10000
    x_max = max([ 2*s.binning()[-2]-s.binning()[-3] for s in spectra ])
    for s in spectra:
        s.read()
        s.give_x_max(x_max)
        s.draw_method = 'repr_vertical_bar_with_horizontal_lines_dashed_onlymerged'

    if combWithHbb_given and PLOT_SYSTEMATIC_ONLY:
        # Get syst only shape
        combWithHbb_statonly = differentials.scans.DifferentialSpectrum('combWithHbb_statonly', scandict.combWithHbb_statonly)
        combWithHbb_statonly.set_sm(obstuple.combWithHbb.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
        combWithHbb_statonly.read()
        systshapemaker = differentials.systshapemaker.SystShapeMaker()
        systonly_histogram, systonly_histogram_xs = systshapemaker.get_systonly_histogram(combWithHbb, combWithHbb_statonly)

    # Get SM histograms
    sm_xs, sm_ratio = get_sm_histograms(obstuple.combWithHbb, normalize_by_second_to_last_bin_width=True, x_max=x_max)

    if args.table:

        table = differentials.plotting.newtables.BaseTable()
        table.end_line_with_tab_sep = True

        rowproducer = differentials.plotting.tableproducer.SpectrumRowProducer(combWithHbb.binning())
        rowproducer.do_xs = True
        rowproducer.normalize = True
        table.append(rowproducer.produce_binning_row('pT (GeV)'))
        table.append(rowproducer.produce(hgg))
        table.append(rowproducer.produce(hzz))
        table.append(rowproducer.produce(hbb))
        table.append(rowproducer.produce(combWithHbb))

        print table.produce_table_string()
        return

        # table = differentials.plotting.tables.SpectraTable('pth_smH', [s for s in spectra])
        # table.print_only_symm_unc = True
        # table.add_symm_improvement_row(hgg, combWithHbb)
        # logging.info('Table:\n{0}'.format( table.repr_terminal() ))
        # return

    # Start compiling plot
    plotname = (
        'projectionspectra_{0}'.format(obs_name)
        + ('_asimov' if args.asimov else '')
        + ( '_nonfixedbinwidth' if not APPLY_FIXED_BINNING else '' )
        + ( '_statonly' if DO_STAT_ONLY else '' )
        + ( '_lumi36' if DO_LUMI_36 else '' )
        + ( '_scenario2' if DO_S2 else '' )
        )
    plot = differentials.plotting.plots.SpectraPlot(plotname, spectra)
    plot.draw_multiscans = True
    plot.obsname = obs_name
    plot.obsunit = 'GeV'

    # Add the SM and syst-only histograms
    if combWithHbb_given and PLOT_SYSTEMATIC_ONLY and systshapemaker.success:
        plot.add_top(systonly_histogram_xs, systonly_histogram_xs.draw_method, plot.leg)
        plot.add_bottom(systonly_histogram, systonly_histogram.draw_method)
    plot.add_top(sm_xs, 'repr_basic_with_full_fill', plot.leg)
    plot.add_bottom(sm_ratio, 'repr_basic_with_full_fill')

    # Some ranges
    plot.top_y_min = 0.9*10e-6
    plot.top_y_max = 10.
    # plot.bottom_y_min = 0.2
    # plot.bottom_y_max = 1.7

    plot.scans_x_min = 0.5
    plot.scans_x_max = 1.5
    # plot.scans_y_min = -1.0
    # plot.scans_y_max = 5.0

    if APPLY_FIXED_BINNING:
        # Apply fixed binning
        reference_binning = hgg.binning()
        plot.make_fixed_widths(reference_binning)
        plot.top_x_max = len(reference_binning)-1
        plot.bottom_x_max = len(reference_binning)-1
        hgg.style().bin_center_offset = -0.17
        hzz.style().bin_center_offset = 0.17
        hbb.style().bin_center_offset = 0.17
        hzz.style().plot_priority = 8
        plot.overflow_label_base_offset = 0.65
        plot.add_lines_at_bin_boundaries(range(1, len(reference_binning)-1))
    else:
        hgg.style().bin_center_offset = -0.17
        hzz.style().bin_center_offset = 0.17
        hbb.style().bin_center_offset = 0.17
        plot.set_reference_bounds(hgg.binning())
        plot.add_lines_at_bin_boundaries()
        plot.overflow_label_base_offset = 0.29

    lw = 0.42
    lh = 0.54 * 0.9

    if APPLY_FIXED_BINNING:
        plot.leg.set(
            lambda c: c.GetLeftMargin() + 0.02,
            lambda c: c.GetBottomMargin() + 0.09,
            lambda c: c.GetLeftMargin() + 0.02 + lw,
            lambda c: c.GetBottomMargin() + 0.09 + lh,
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
    plot.lumi_text_size = 0.06
    plot.draw()

    l = differentials.plotting.pywrappers.Latex(
        lambda c: c.GetLeftMargin() + 0.04,
        lambda c: c.GetBottomMargin() + 0.05,
        '#sigma_{SM} from DOI: 10.23731/CYRM-2017-002'
        )
    if not APPLY_FIXED_BINNING:
        l.x = lambda c: c.GetLeftMargin() + 0.04 + xshift
        l.y = lambda c: c.GetBottomMargin() + 0.05 + yshift
    l.SetNDC()
    l.SetTextAlign(11)
    l.SetTextFont(42) 
    l.SetTextSize(0.038)
    l.Draw()

    plot.base_bottom.GetYaxis().SetNdivisions(505)

    if APPLY_FIXED_BINNING: plot.replace_bin_labels([ '0', '15', '30', '45', '80', '120', '200', '350', '600', '#infty' ])
    plot.wrapup()
    differentials.plotting.pywrappers.CMS_Latex_lumi.CMS_lumi = 35.9




def filter_hbb(hbb):

    graphs = []
    for scan in hbb.scans:
        scanfilter = differentials.onedimscanfilter.OneDimScanFilter(scan.x(), scan.y())
        scanfilter.filter_clear_nonsense()
        graph = scanfilter.to_graph()
        graph.title = scan.x_variable
        graphs.append(graph)

    plot = differentials.plotting.plots.MultiScanPlot('hbb_filter_test')
    plot.scans.extend(hbb.scans)
    plot.manual_graphs.extend(graphs)
    plot.x_min = 0.0
    plot.x_max = 2.0
    plot.draw()
    plot.wrapup()





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
