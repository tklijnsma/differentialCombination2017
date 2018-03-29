#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

from OptionHandler import flag_as_option

import copy
import LatestPaths
import LatestPathsGetters
import LatestBinning

import differentials
import differentialutils

from time import strftime
datestr = strftime( '%b%d' )

########################################
# Main
########################################

x_coupling = 'ct'
y_coupling = 'cg'

top_x_min = -5.0
top_x_max = 5.0
top_y_min = -0.40
top_y_max = 0.40

@flag_as_option
def preapproval_plots_Top(args):
    points_on_contour_Top(args)
    multicont_Top_reweighted(args)
    multicont_Top_lumi300fb(args)
    args = differentialutils.force_asimov(args)
    multicont_Top_reweighted(args)

@flag_as_option
def multicont_Top_reweighted(args):
    scans = []
    scandict = LatestPaths.scan.top.reweighted.asimov if args.asimov else LatestPaths.scan.top.reweighted.observed

    if not args.asimov:
        hgg = differentials.scans.Scan2D('hgg', x_coupling, y_coupling, scandir = scandict.hgg)
        hgg.color = 2
        hgg.read()
        scans.append(hgg)

        hzz = differentials.scans.Scan2D('hzz', x_coupling, y_coupling, scandir = scandict.hzz)
        hzz.color = 4
        hzz.read()
        scans.append(hzz)

    combWithHbb = differentials.scans.Scan2D('combWithHbb', x_coupling, y_coupling, scandir = scandict.combWithHbb)
    combWithHbb.color = 1
    combWithHbb.read()
    scans.append(combWithHbb)
    
    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Top_reweighted' + ('_asimov' if args.asimov else ''),
        scans,
        x_min=top_x_min, x_max=top_x_max, y_min=top_y_min, y_max=top_y_max,
        )
    # plot.legend.SetNColumns(2)
    # plot.only_1sigma_contours = True
    plot.draw()


@flag_as_option
def multicont_Top_profiledTotalXS(args):
    args = differentialutils.force_asimov(args)

    combWithHbb_profiledTotalXS = differentials.scans.Scan2D('combWithHbb_profiledTotalXS', x_coupling, y_coupling,
        scandir = 'out/Scan_Mar13_Top_combWithHbb_profiledTotalXS_asimov'
        )
    combWithHbb_profiledTotalXS.color = 1
    combWithHbb_profiledTotalXS.title = 'Profiled #sigma_{tot}'
    combWithHbb_profiledTotalXS.read()
    
    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Top_profiledTotalXS' + ('_asimov' if args.asimov else ''),
        [combWithHbb_profiledTotalXS],
        x_min=top_x_min, x_max=top_x_max, y_min=top_y_min, y_max=top_y_max,
        )
    plot.only_1sigma_contours = True
    plot.draw()


@flag_as_option
def multicont_Top_lumi300fb(args):
    args = differentialutils.force_asimov(args)

    combWithHbb = differentials.scans.Scan2D('combWithHbb', x_coupling, y_coupling,
        scandir = LatestPaths.scan.top.reweighted.asimov.combWithHbb
        )
    combWithHbb.color = 1
    combWithHbb.title = '35.9 fb^{-1}'
    combWithHbb.read()

    lumi300fb = differentials.scans.Scan2D('lumi300fb', x_coupling, y_coupling,
        scandir = LatestPaths.scan.top.lumi300fb
        )
    lumi300fb.color = 4
    lumi300fb.title = '300 fb^{-1}'
    # lumi300fb.contour_filter_method = 'max_distance_to_com'
    lumi300fb.read()
    
    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Top_lumi300fb' + ('_asimov' if args.asimov else ''),
        [combWithHbb, lumi300fb],
        x_min=top_x_min, x_max=top_x_max, y_min=top_y_min, y_max=top_y_max,
        )
    # plot.only_1sigma_contours = True
    plot.draw()


@flag_as_option
def multicont_Top_hbbinclusion_plot(args):
    args = differentialutils.force_asimov(args)
    scans = []

    # Black
    combination_last1dropped = differentials.scans.Scan2D('combination_last1dropped', x_coupling, y_coupling,
        scandir = 'out/Scan_Mar07_Top_combination_asimov_5'
        )
    combination_last1dropped.color = 1
    # combination_last1dropped.title = 'GT600 dropped'
    # combination_last1dropped.title = 'No Hbb, no bound@600'
    combination_last1dropped.title = 'Low p_{T}'
    combination_last1dropped.read()
    scans.append(combination_last1dropped)

    # Red
    combination_nothingdropped = differentials.scans.Scan2D('combination_nothingdropped', x_coupling, y_coupling,
        scandir = 'out/Scan_Mar13_Top_combination_noBinsDropped_asimov'
        )
    combination_nothingdropped.color = 2
    # combination_nothingdropped.title = 'Nothing dropped'
    # combination_nothingdropped.title = 'No Hbb, with bound@600'
    combination_nothingdropped.title = 'Low p_{T} + high p_{T} H#rightarrow#gamma#gamma'
    combination_nothingdropped.read()
    scans.append(combination_nothingdropped)

    # Blue
    combWithHbb_last1dropped = differentials.scans.Scan2D('combWithHbb_last1dropped', x_coupling, y_coupling,
        # scandir = 'out/Scan_Mar07_Top_combWithHbb_asimov'
        scandir = 'out/Scan_Mar12_Top_combWithHbb_last2BinsDropped_asimov'
        )
    combWithHbb_last1dropped.color = 4
    # combWithHbb_last1dropped.title = 'GT600 dropped'
    # combWithHbb_last1dropped.title = 'With Hbb, no bound@600'
    combWithHbb_last1dropped.title = 'Low p_{T} + high p_{T} H#rightarrowbb'
    combWithHbb_last1dropped.read()
    scans.append(combWithHbb_last1dropped)

    # Green
    combWithHbb_nothingdropped = differentials.scans.Scan2D('combWithHbb_nothingdropped', x_coupling, y_coupling,
        scandir = 'out/Scan_Mar13_Top_combWithHbb_noBinsDropped_asimov'
        )
    combWithHbb_nothingdropped.color = 8
    # combWithHbb_nothingdropped.title = 'Nothing dropped'
    # combWithHbb_nothingdropped.title = 'With Hbb, with bound@600'
    combWithHbb_nothingdropped.title = 'Full combination'
    combWithHbb_nothingdropped.read()
    scans.append(combWithHbb_nothingdropped)

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Top_hbbinclusion_plot' + ('_asimov' if args.asimov else ''),
        scans,
        x_min=top_x_min, x_max=top_x_max, y_min=top_y_min, y_max=top_y_max,
        )
    plot.legend.SetNColumns(2)
    plot.only_1sigma_contours = True
    plot.draw()


@flag_as_option
def multicont_Top_droppingbins_combWithHbb(args):
    args = differentialutils.force_asimov(args)
    scans = []

    combWithHbb_nothingdropped = differentials.scans.Scan2D('combWithHbb_nothingdropped', x_coupling, y_coupling,
        # scandir = 'out/Scan_Mar06_Top_combWithHbb_asimov'
        scandir = 'out/Scan_Mar13_Top_combWithHbb_noBinsDropped_asimov'
        )
    combWithHbb_nothingdropped.color = 1
    combWithHbb_nothingdropped.title = 'Nothing dropped'
    combWithHbb_nothingdropped.read()
    scans.append(combWithHbb_nothingdropped)

    combWithHbb_last1dropped = differentials.scans.Scan2D('combWithHbb_last1dropped', x_coupling, y_coupling,
        scandir = 'out/Scan_Mar07_Top_combWithHbb_asimov'
        )
    combWithHbb_last1dropped.color = 2
    combWithHbb_last1dropped.title = 'GT600 dropped'
    combWithHbb_last1dropped.read()
    scans.append(combWithHbb_last1dropped)

    combWithHbb_last2dropped = differentials.scans.Scan2D('combWithHbb_last2dropped', x_coupling, y_coupling,
        scandir = 'out/Scan_Mar12_Top_combWithHbb_last2BinsDropped_asimov'
        )
    combWithHbb_last2dropped.color = 4
    combWithHbb_last2dropped.title = '350-600 and GT600 dropped'
    combWithHbb_last2dropped.read()
    scans.append(combWithHbb_last2dropped)

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Top_droppingbins_combWithHbb' + ('_asimov' if args.asimov else ''),
        scans,
        x_min=top_x_min, x_max=top_x_max, y_min=top_y_min, y_max=top_y_max,
        )
    plot.only_1sigma_contours = True
    plot.draw()


@flag_as_option
def multicont_Top_droppingbins_combination(real_args):
    args = copy.deepcopy(real_args)
    args.asimov = True
    scans = []

    combination_nothingdropped = differentials.scans.Scan2D('combination_nothingdropped', x_coupling, y_coupling,
        # scandir = 'out/Scan_Mar06_Top_combination_asimov'
        scandir = 'out/Scan_Mar13_Top_combination_noBinsDropped_asimov'
        )
    combination_nothingdropped.color = 1
    combination_nothingdropped.title = 'Nothing dropped'
    combination_nothingdropped.read()
    scans.append(combination_nothingdropped)

    combination_last1dropped = differentials.scans.Scan2D('combination_last1dropped', x_coupling, y_coupling,
        scandir = 'out/Scan_Mar07_Top_combination_asimov_5'
        )
    combination_last1dropped.color = 2
    combination_last1dropped.title = 'GT600 dropped'
    combination_last1dropped.read()
    scans.append(combination_last1dropped)

    combination_last2dropped = differentials.scans.Scan2D('combination_last2dropped', x_coupling, y_coupling,
        scandir = 'out/Scan_Mar12_Top_combination_last2BinsDropped_asimov'
        )
    combination_last2dropped.color = 4
    combination_last2dropped.title = '350-600 and GT600 dropped'
    combination_last2dropped.read()
    scans.append(combination_last2dropped)

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Top_droppingbins_combination' + ('_asimov' if args.asimov else ''),
        scans,
        x_min=top_x_min, x_max=top_x_max, y_min=top_y_min, y_max=top_y_max,
        )
    plot.only_1sigma_contours = True
    plot.draw()


@flag_as_option
def multicont_Top_nominal(args):
    scans = []

    hgg = differentials.scans.Scan2D('hgg', x_coupling, y_coupling,
        scandir = 'out/Scan_Mar06_Top_hgg_asimov'  if args.asimov else 'out/Scan_Mar03_Top_hgg'
        )
    hgg.color = 2
    hgg.read()
    scans.append(hgg)

    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling,
        scandir = 'out/Scan_Mar06_Top_combination_asimov'  if args.asimov else 'out/Scan_Mar03_Top_combination'
        )
    combination.color = 1
    combination.read()
    scans.append(combination)

    combWithHbb = differentials.scans.Scan2D('combWithHbb', x_coupling, y_coupling,
        scandir = 'out/Scan_Mar06_Top_combWithHbb_asimov'  if args.asimov else 'out/Scan_Mar03_Top_combWithHbb'
        )
    combWithHbb.color = 9
    combWithHbb.read()
    scans.append(combWithHbb)
   

    
    # hgg = differentials.scans.Scan2D('hgg', x_coupling, y_coupling,
    #     scandir = get_nominal(args).hgg
    #     )
    # hgg.color = 2
    # hgg.read()

    # hzz = differentials.scans.Scan2D('hzz', x_coupling, y_coupling,
    #     scandir = get_nominal(args).hzz
    #     )
    # hzz.color = 4
    # hzz.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Top_nominal' + ('_asimov' if args.asimov else ''),
        scans,
        x_min=top_x_min, x_max=top_x_max, y_min=top_y_min, y_max=top_y_max,
        )
    plot.draw()


@flag_as_option
def points_on_contour_Top(args):
    obs_name = 'pth_ggH'
    obstuple = LatestBinning.obstuple_pth_ggH

    scandict_2D = LatestPaths.scan.top.reweighted.asimov if args.asimov else LatestPaths.scan.top.reweighted.observed
    Top_combWithHbb = differentials.scans.Scan2D('combWithHbb', x_coupling, y_coupling, scandir = scandict_2D.combWithHbb)
    Top_combWithHbb.color = 1
    Top_combWithHbb.read()

    # obs = LatestBinning.obs_pth_ggH
    # obs.drop_bins_up_to_value(125.)

    scandict = LatestPaths.scan.pth_ggH.asimov if args.asimov else LatestPaths.scan.pth_ggH.observed
    combWithHbb = differentials.scans.DifferentialSpectrum('combWithHbb', scandict.combWithHbb)
    combWithHbb.color = 1
    combWithHbb.no_overflow_label = True
    combWithHbb.draw_method = 'repr_point_with_vertical_bar'
    combWithHbb.set_sm(obstuple.combWithHbb.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    combWithHbb.read()

    # ws = 'out/workspaces_Dec11/combinedCard_Nov03_CouplingModel_TopHighPt_withTheoryUncertainties.root'
    # ws = 'out/workspaces_Mar06/combWithHbb_Top_reweighted_nominal.root'
    ws = LatestPaths.ws.top.nominal.combWithHbb

    # ======================================
    # Load into plot

    plot = differentials.plotting.plots.BottomPanelPlotWithParametrizations('points_on_contour_Top')
    plot.scan2D = Top_combWithHbb
    plot.ws_file = ws
    plot.ptspectrum = combWithHbb
    plot.obs = obstuple.combWithHbb

    plot.default_points_xy_maxima = False

    # plot.top_y_min = top_y_min
    plot.top_y_max = 500.

    plot.y_SM = 0.0

    obs_title = 'p_{T}'
    obs_unit  = 'GeV'
    plot.x_title = obs_title + (' ({0})'.format(obs_unit) if not(obs_unit is None) else '')
    plot.y_title_top = '#Delta#sigma/#Delta{0} (pb{1})'.format(
        obs_title,
        '/' + obs_unit if not(obs_unit is None) else ''
        )
    plot.y_title_bottom = '#mu'

    plot.draw()


#____________________________________________________________________
# kappat kappab

@flag_as_option
def multicont_TopCtCb(args):
    scans = []
    scandict = LatestPaths.scan.topctcb.reweighted.asimov if args.asimov else LatestPaths.scan.topctcb.reweighted.observed
    y_coupling = 'cb'

    hgg = differentials.scans.Scan2D('hgg', x_coupling, y_coupling, scandir = scandict.hgg)
    hgg.color = 2
    hgg.read()
    scans.append(hgg)

    hzz = differentials.scans.Scan2D('hzz', x_coupling, y_coupling, scandir = scandict.hzz)
    hzz.color = 4
    hzz.read()
    scans.append(hzz)

    combWithHbb = differentials.scans.Scan2D('combWithHbb', x_coupling, y_coupling, scandir = scandict.combWithHbb)
    combWithHbb.color = 1
    combWithHbb.read()
    scans.append(combWithHbb)
    
    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_TopCtCb' + ('_asimov' if args.asimov else ''),
        scans,
        # x_min=top_x_min, x_max=top_x_max, y_min=top_y_min, y_max=top_y_max,
        )
    # plot.legend.SetNColumns(2)
    # plot.only_1sigma_contours = True
    plot.draw()
