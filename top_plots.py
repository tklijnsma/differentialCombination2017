#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

from OptionHandler import flag_as_option

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
def points_on_contour_Top(args):

    combination = differentials.scans.Scan2D(
        'combination', x_coupling, y_coupling,
        scandir=LatestPaths.scan_combined_TopHighPt
        )
    combination.color = 1
    combination.read()

    obs = LatestBinning.obs_pth_ggH
    # obs.drop_bins_up_to_value(125.)

    obs_name = 'pth_ggH'
    pth_smH_combination = differentials.scans.DifferentialSpectrum('combination',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='combination', statonly=False),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='combination')
        )
    pth_smH_combination.color = 1
    pth_smH_combination.no_overflow_label = True
    pth_smH_combination.draw_method = 'repr_point_with_vertical_bar'
    pth_smH_combination.set_sm(obs.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    # pth_smH_combination.drop_bins_up_to_value(125.)
    pth_smH_combination.read()

    ws = 'out/workspaces_Dec11/combinedCard_Nov03_CouplingModel_TopHighPt_withTheoryUncertainties.root'


    # ======================================
    # Load into plot

    plot = differentials.plotting.multipanel.BottomPanelPlotWithParametrizations('points_on_contour_Top')
    plot.scan2D = combination
    plot.ws_file = ws
    plot.ptspectrum = pth_smH_combination
    plot.obs = obs

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
