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

from time import strftime
datestr = strftime( '%b%d' )

########################################
# Main
########################################

x_coupling = 'ct'
y_coupling = 'cg'


@flag_as_option
def multicont_Top_chi2XC(real_args):
    args = copy.deepcopy(real_args)
    args.asimov = True
    scans = []

    combination_lastBinDroppedHgg = differentials.scans.Scan2D('combination_lastBinDroppedHgg', x_coupling, y_coupling,
        scandir = 'out/Scan_Mar07_Top_combination_asimov_5'
        )
    combination_lastBinDroppedHgg.color = 1
    combination_lastBinDroppedHgg.title = 'Old combination'
    combination_lastBinDroppedHgg.read()
    scans.append(combination_lastBinDroppedHgg)

    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling,
        scandir = 'out/Scan_Mar06_Top_combination_asimov'  if args.asimov else 'out/Scan_Mar03_Top_combination'
        )
    combination.color = 2
    combination.title = 'H#gamma#gamma new bins'
    combination.read()
    scans.append(combination)

    combWithHbb_lastBinDroppedHgg = differentials.scans.Scan2D('combWithHbb_lastBinDroppedHgg', x_coupling, y_coupling,
        scandir = 'out/Scan_Mar07_Top_combWithHbb_asimov'
        )
    combWithHbb_lastBinDroppedHgg.title = 'Hbb incl.'
    combWithHbb_lastBinDroppedHgg.color = 4
    combWithHbb_lastBinDroppedHgg.read()
    scans.append(combWithHbb_lastBinDroppedHgg)

    combWithHbb = differentials.scans.Scan2D('combWithHbb', x_coupling, y_coupling,
        scandir = 'out/Scan_Mar06_Top_combWithHbb_asimov'  if args.asimov else 'out/Scan_Mar03_Top_combWithHbb'
        )
    combWithHbb.title = 'Hgg new b. + Hbb'
    combWithHbb.color = 8
    combWithHbb.read()
    scans.append(combWithHbb)

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Top_chi2XC' + ('_asimov' if args.asimov else ''),
        scans,
        x_min=None, x_max=None, y_min=None, y_max=None
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
        x_min=None, x_max=None, y_min=None, y_max=None
        )
    plot.draw()


@flag_as_option
def points_on_contour_Top(args):
    obs_name = 'pth_ggH'
    obstuple = LatestBinning.obstuple_pth_ggH

    Top_combWithHbb = differentials.scans.Scan2D('combWithHbb', x_coupling, y_coupling,
        scandir = 'out/Scan_Mar06_Top_combWithHbb_asimov'  if args.asimov else 'out/Scan_Mar03_Top_combWithHbb'
        )
    Top_combWithHbb.color = 9
    Top_combWithHbb.read()

    # obs = LatestBinning.obs_pth_ggH
    # obs.drop_bins_up_to_value(125.)

    dir_combWithHbb = 'out/Scan_Mar02_pth_ggH_combWithHbb'
    if args.asimov:
        dir_combWithHbb = 'out/Scan_Mar02_pth_ggH_combWithHbb_asimov'
    combWithHbb = differentials.scans.DifferentialSpectrum('combWithHbb', dir_combWithHbb)
    combWithHbb.scandirs.append('out/Scan_Mar07_pth_ggH_combWithHbb_asimov')
    # combWithHbb.color = 47
    combWithHbb.color = 1
    combWithHbb.no_overflow_label = True
    combWithHbb.draw_method = 'repr_point_with_vertical_bar'
    combWithHbb.set_sm(obstuple.combWithHbb.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    combWithHbb.read()


    # ws = 'out/workspaces_Dec11/combinedCard_Nov03_CouplingModel_TopHighPt_withTheoryUncertainties.root'
    ws = 'out/workspaces_Mar06/combWithHbb_Top_reweighted_nominal.root'


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
