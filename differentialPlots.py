#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import logging
import os, sys, re
from os.path import *
from glob import glob
from copy import deepcopy
from math import sqrt

from OptionHandler import flag_as_option

sys.path.append('src')
import Commands
import PhysicsCommands
import TheoryCommands
import LatestPaths
import LatestPathsGetters
import LatestBinning
from Container import Container
import PlotCommands
import DifferentialTable
from differentialTools import *

import differentials
# import differentials.plotting

from time import strftime
datestr = strftime('%b%d')


########################################
# Plotting
########################################


@flag_as_option
def plot_pt_combination(args):
    spectra = []
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
    obs_name = 'pth_ggH'
    obstuple = LatestBinning.obstuple_pth_ggH

    if args.asimov:
        dir_combWithHbb = 'out/Scan_Mar02_pth_ggH_combWithHbb_asimov'
        dir_combination = 'out/Scan_Mar02_pth_ggH_combination_asimov'
        dir_hbb = 'out/Scan_Mar02_pth_ggH_hbb_asimov'
        dir_hgg = 'out/Scan_Mar02_pth_ggH_hgg_asimov'
        dir_hzz = 'out/Scan_Mar02_pth_ggH_hzz_asimov'
    else:
        dir_combWithHbb = 'out/Scan_Mar02_pth_ggH_combWithHbb'
        dir_combination = 'out/Scan_Mar02_pth_ggH_combination'
        dir_hbb = 'out/Scan_Mar02_pth_ggH_hbb'
        dir_hgg = 'out/Scan_Mar02_pth_ggH_hgg'
        dir_hzz = 'out/Scan_Mar02_pth_ggH_hzz'
    
    hgg = differentials.scans.DifferentialSpectrum('hgg', dir_hgg)
    hgg.color = 2
    hgg.no_overflow_label = True
    hgg.set_sm(obstuple.hgg.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    hgg.read()
    spectra.append(hgg)

    hzz = differentials.scans.DifferentialSpectrum('hzz', dir_hzz)
    hzz.color = 4
    hzz.no_overflow_label = True
    hzz.set_sm(obstuple.hzz.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    hzz.read()
    spectra.append(hzz)

    hbb = differentials.scans.DifferentialSpectrum('hbb', dir_hbb)
    hbb.drop_first_bin()
    hbb.color = 30
    hbb.no_overflow_label = True
    hbb.set_sm(obstuple.hbb.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    hbb.read()
    spectra.append(hbb)

    combination = differentials.scans.DifferentialSpectrum('combination', dir_combination)
    combination.color = 1
    combination.no_overflow_label = True
    combination.draw_method = 'repr_point_with_vertical_bar'
    combination.set_sm(obstuple.combination.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    combination.read()
    # spectra.append(combination)

    combWithHbb = differentials.scans.DifferentialSpectrum('combWithHbb', dir_combWithHbb)
    combWithHbb.scandirs.append('out/Scan_Mar07_pth_ggH_combWithHbb_asimov')
    # combWithHbb.color = 47
    combWithHbb.color = 1
    combWithHbb.no_overflow_label = True
    combWithHbb.draw_method = 'repr_point_with_vertical_bar'
    combWithHbb.set_sm(obstuple.combWithHbb.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    combWithHbb.read()
    spectra.append(combWithHbb)

    if args.table:
        table = differentials.plotting.tables.SpectraTable('pth_ggH', [s for s in spectra if not s is hbb])
        table.add_symm_improvement_row(hgg, combWithHbb)
        table.add_symm_improvement_row(hgg, combination)
        table.add_symm_improvement_row(combination, combWithHbb)
        logging.info('Table:\n{0}'.format( table.repr_terminal() ))
        return

    l = differentials.plotting.pywrappers.Latex(
        lambda c: 1.0 - c.GetRightMargin() - 0.01,
        lambda c: 1.0 - c.GetTopMargin() - 0.14,
        'gluon fusion cross section'
        )
    l.SetNDC()
    l.SetTextSize(0.05)
    l.SetTextAlign(33)

    plot = differentials.plotting.plots.SpectraPlot(
        'debugspectra_{0}'.format(obs_name) + ('_asimov' if args.asimov else ''),
        spectra
        )
    plot.draw_multiscans = True
    plot.obsname = obs_name
    plot.obsunit = 'GeV'
    plot.add_top(l, '')
    plot.leg.SetNColumns(3)
    plot.draw()
    plot.wrapup()


def get_obs(obs_name, decay_channel):
    return {
        'pth_smH'  : pth_smH_obs,
        'pth_ggH'  : pth_ggH_obs,
        'njets'    : njets_obs,
        'ptjet'    : ptjet_obs,
        'rapidity' : rapidity_obs,
        }[obs_name](decay_channel)

def pth_smH_obs(decay_channel):
    return {
        'hgg' : LatestBinning.obs_pth,
        'hzz' : LatestBinning.obs_pth_hzzBinning,
        'combination' : LatestBinning.obs_pth,
        }[decay_channel]

def pth_ggH_obs(decay_channel):
    return {
        'hgg' : LatestBinning.obs_pth_ggH,
        'hzz' : LatestBinning.obs_pth_ggH_hzzBinning,
        'combination' : LatestBinning.obs_pth_ggH,
        }[decay_channel]

def njets_obs(decay_channel):
    return {
        'hgg' : LatestBinning.obs_njets,
        'hzz' : LatestBinning.obs_njets_hzzBinning,
        'combination' : LatestBinning.obs_njets,
        }[decay_channel]

def ptjet_obs(decay_channel):
    return {
        'hgg' : LatestBinning.obs_ptjet,
        'hzz' : LatestBinning.obs_ptjet_hzzBinning,
        'combination' : LatestBinning.obs_ptjet,
        }[decay_channel]

def rapidity_obs(decay_channel):
    return {
        'hgg' : LatestBinning.obs_yh,
        'hzz' : LatestBinning.obs_yh,
        'combination' : LatestBinning.obs_yh,
        }[decay_channel]


#____________________________________________________________________
@flag_as_option
def plot_all_differentials(args):
    pth_smH_plot(args)
    pth_ggH_plot(args)
    # pth_ggH_hbb_plot(args)
    njets_plot(args)
    ptjet_plot(args)
    rapidity_plot(args)

#____________________________________________________________________
@flag_as_option
def pth_smH_plot(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
    obs_name = 'pth_smH'
    obstuple = LatestBinning.obstuple_pth_smH
    
    hgg = differentials.scans.DifferentialSpectrum('hgg',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='hgg', statonly=False),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='hgg')
        )
    hgg.color = 2
    hgg.set_sm(obstuple.hgg.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))

    hzz = differentials.scans.DifferentialSpectrum('hzz',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='hzz', statonly=False),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='hzz')
        )
    hzz.color = 4
    hzz.set_sm(obstuple.hzz.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))

    combination = differentials.scans.DifferentialSpectrum('combination',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='combination', statonly=False),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='combination')
        )
    combination.color = 1
    combination.no_overflow_label = True
    combination.draw_method = 'repr_point_with_vertical_bar'
    combination.set_sm(obstuple.combination.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))

    plot = differentials.plotting.plots.SpectraPlot('spectra_{0}'.format(obs_name),
        [ hgg, hzz, combination ]
        )
    plot.obsname = obs_name
    plot.obsunit = 'GeV'
    plot.draw()
    plot.wrapup()


#____________________________________________________________________
@flag_as_option
def pth_ggH_plot(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
    obs_name = 'pth_ggH'

    hgg = differentials.scans.DifferentialSpectrum('hgg',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='hgg', statonly=False),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='hgg')
        )
    hgg.color = 2
    hgg.set_sm( get_obs(obs_name, 'hgg').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )

    hzz = differentials.scans.DifferentialSpectrum('hzz',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='hzz', statonly=False),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='hzz')
        )
    hzz.color = 4
    hzz.set_sm( get_obs(obs_name, 'hzz').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )

    combination = differentials.scans.DifferentialSpectrum('combination',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='combination', statonly=False),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='combination')
        )
    combination.color = 1
    combination.no_overflow_label = True
    combination.draw_method = 'repr_point_with_vertical_bar'
    combination.set_sm( get_obs(obs_name, 'combination').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )

    l = differentials.plotting.pywrappers.Latex(
        lambda c: 1.0 - c.GetRightMargin() - 0.01,
        lambda c: 1.0 - c.GetTopMargin() - 0.14,
        'gluon fusion cross section'
        )
    l.SetNDC()
    l.SetTextSize(0.05)
    l.SetTextAlign(33)

    plot = differentials.plotting.plots.SpectraPlot('spectra_{0}'.format(obs_name),
        [ hgg, hzz, combination ]
        )
    plot.obsname = obs_name
    plot.obsunit = 'GeV'
    plot.add_top(l, '')
    plot.draw()
    plot.wrapup()



# out/Scan_njets_Feb06_combination_asimov
# out/Scan_njets_Feb06_hgg_asimov
# out/Scan_njets_Feb06_hzz_asimov
# out/Scan_pth_ggH_Feb06_combination_asimov
# out/Scan_pth_ggH_Feb06_hgg_asimov
# out/Scan_pth_ggH_Feb06_hzz_asimov
# out/Scan_pth_smH_Feb06_combination_asimov
# out/Scan_pth_smH_Feb06_hgg_asimov
# out/Scan_pth_smH_Feb06_hzz_asimov
# out/Scan_ptjet_Feb06_combination_asimov
# out/Scan_ptjet_Feb06_hgg_asimov
# out/Scan_ptjet_Feb06_hzz_asimov
# out/Scan_rapidity_Feb06_combination_asimov
# out/Scan_rapidity_Feb06_hgg_asimov
# out/Scan_rapidity_Feb06_hzz_asimov

def get_POIs_oldstyle_scandir(scandir):
    root_files = glob(join(scandir, '*.root'))
    POIs = []
    for root_file in root_files:
        POI = basename(root_file).split('_')[1:6]
        try:
            differentials.core.str_to_float(POI[-1])
        except ValueError:
            POI = POI[:-1]
        POI = '_'.join(POI)
        POIs.append(POI)
    POIs = list(set(POIs))
    POIs.sort(key=differentials.core.range_sorter)
    logging.info('Retrieved following POIs from oldstyle {0}:\n{1}'.format(scandir, POIs))
    return POIs

#____________________________________________________________________
@flag_as_option
def njets_plot(args):
    obs_name = 'njets'

    if args.asimov:
        combination_dir = 'out/Scan_njets_Feb06_combination_asimov'
        hgg_dir = 'out/Scan_njets_Feb06_hgg_asimov'
        hzz_dir = 'out/Scan_njets_Feb06_hzz_asimov'
    else:
        raise NotImplementedError

    hgg = differentials.scans.DifferentialSpectrum('hgg', scandir=hgg_dir)
    hgg.POIs = get_POIs_oldstyle_scandir(hgg_dir)
    hgg.color = 2
    hgg.set_sm( get_obs(obs_name, 'hgg').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )
    hgg.read()

    hzz = differentials.scans.DifferentialSpectrum('hzz', scandir=hzz_dir)
    hzz.POIs = get_POIs_oldstyle_scandir(hzz_dir)
    hzz.color = 4
    hzz.set_sm( get_obs(obs_name, 'hzz').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )
    hzz.read()

    combination = differentials.scans.DifferentialSpectrum('combination', scandir=combination_dir)
    combination.POIs = get_POIs_oldstyle_scandir(combination_dir)
    combination.color = 1
    combination.no_overflow_label = True
    combination.draw_method = 'repr_point_with_vertical_bar'
    combination.set_sm( get_obs(obs_name, 'combination').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )
    combination.read()

    if args.table:
        table = differentials.plotting.tables.SpectraTable('njets', [ hgg, hzz, combination ])
        table.add_symm_improvement_row(hgg, combination)
        logging.info('Table:\n{0}'.format(table.repr_terminal()))
        return

    plot = differentials.plotting.plots.SpectraPlot('spectra_{0}'.format(obs_name),
        [ hgg, hzz, combination ]
        )
    plot.obsname = obs_name
    plot.draw()
    plot.wrapup()


#____________________________________________________________________
@flag_as_option
def ptjet_plot(args):
    obs_name = 'ptjet'

    if args.asimov:
        combination_dir = 'out/Scan_ptjet_Feb06_combination_asimov'
        hgg_dir = 'out/Scan_ptjet_Feb06_hgg_asimov'
        hzz_dir = 'out/Scan_ptjet_Feb06_hzz_asimov'
    else:
        raise NotImplementedError

    hgg = differentials.scans.DifferentialSpectrum('hgg', scandir=hgg_dir)
    hgg.POIs = get_POIs_oldstyle_scandir(hgg_dir)
    hgg.color = 2
    hgg.set_sm( get_obs(obs_name, 'hgg').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )

    hzz = differentials.scans.DifferentialSpectrum('hzz', scandir=hzz_dir)
    hzz.POIs = get_POIs_oldstyle_scandir(hzz_dir)
    hzz.color = 4
    hzz.set_sm( get_obs(obs_name, 'hzz').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )

    combination = differentials.scans.DifferentialSpectrum('combination', scandir=combination_dir)
    combination.POIs = get_POIs_oldstyle_scandir(combination_dir)
    combination.color = 1
    combination.no_overflow_label = True
    combination.draw_method = 'repr_point_with_vertical_bar'
    combination.set_sm( get_obs(obs_name, 'combination').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )

    plot = differentials.plotting.plots.SpectraPlot('spectra_{0}'.format(obs_name),
        [ hgg, hzz, combination ]
        )
    plot.top_y_max = 10.
    plot.obsname = obs_name
    plot.obsunit = 'GeV'
    plot.draw()
    plot.wrapup()


#____________________________________________________________________
@flag_as_option
def rapidity_plot(args):
    obs_name = 'rapidity'

    if args.asimov:
        combination_dir = 'out/Scan_rapidity_Feb06_combination_asimov'
        hgg_dir = 'out/Scan_rapidity_Feb06_hgg_asimov'
        hzz_dir = 'out/Scan_rapidity_Feb06_hzz_asimov'
    else:
        raise NotImplementedError

    hgg = differentials.scans.DifferentialSpectrum('hgg', scandir=hgg_dir)
    hgg.POIs = get_POIs_oldstyle_scandir(hgg_dir)
    hgg.color = 2
    hgg.set_sm( get_obs(obs_name, 'hgg').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=False) )

    hzz = differentials.scans.DifferentialSpectrum('hzz', scandir=hzz_dir)
    hzz.POIs = get_POIs_oldstyle_scandir(hzz_dir)
    hzz.color = 4
    hzz.set_sm( get_obs(obs_name, 'hzz').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=False) )

    combination = differentials.scans.DifferentialSpectrum('combination', scandir=combination_dir)
    combination.POIs = get_POIs_oldstyle_scandir(combination_dir)
    combination.color = 1
    combination.no_overflow_label = True
    combination.draw_method = 'repr_point_with_vertical_bar'
    combination.set_sm( get_obs(obs_name, 'combination').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=False) )

    plot = differentials.plotting.plots.SpectraPlot('spectra_{0}'.format(obs_name),
        [ hgg, hzz, combination ]
        )
    plot.top_y_max = 300
    plot.obsname = obs_name
    plot.draw()
    plot.wrapup()


#____________________________________________________________________
@flag_as_option
def ptjet_plot_old(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
    obs_name = 'ptjet'
    
    hgg = differentials.scans.DifferentialSpectrum('hgg',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='hgg', statonly=False),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='hgg')
        )
    hgg.color = 2
    hgg.set_sm( get_obs(obs_name, 'hgg').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )

    hzz = differentials.scans.DifferentialSpectrum('hzz',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='hzz', statonly=False),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='hzz')
        )
    hzz.color = 4
    hzz.set_sm( get_obs(obs_name, 'hzz').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )

    combination = differentials.scans.DifferentialSpectrum('combination',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='combination', statonly=False),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='combination')
        )
    combination.color = 1
    combination.no_overflow_label = True
    combination.draw_method = 'repr_point_with_vertical_bar'
    combination.set_sm( get_obs(obs_name, 'combination').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )

    plot = differentials.plotting.plots.SpectraPlot('spectra_{0}'.format(obs_name),
        [ hgg, hzz, combination ]
        )
    plot.obsname = obs_name
    plot.obsunit = 'GeV'
    plot.draw()
    plot.wrapup()


#____________________________________________________________________
@flag_as_option
def rapidity_plot_old(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
    obs_name = 'rapidity'
    
    hgg = differentials.scans.DifferentialSpectrum('hgg',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='hgg', statonly=False),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='hgg')
        )
    hgg.color = 2
    hgg.set_sm( get_obs(obs_name, 'hgg').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )

    hzz = differentials.scans.DifferentialSpectrum('hzz',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='hzz', statonly=False),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='hzz')
        )
    hzz.color = 4
    hzz.set_sm( get_obs(obs_name, 'hzz').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )

    combination = differentials.scans.DifferentialSpectrum('combination',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='combination', statonly=False),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='combination')
        )
    combination.color = 1
    combination.no_overflow_label = True
    combination.draw_method = 'repr_point_with_vertical_bar'
    combination.set_sm( get_obs(obs_name, 'combination').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )

    plot = differentials.plotting.plots.SpectraPlot('spectra_{0}'.format(obs_name),
        [ hgg, hzz, combination ]
        )
    plot.obsname = obs_name
    plot.obsunit = 'GeV'
    plot.draw()
    plot.wrapup()


########################################
# Other plots
########################################

#____________________________________________________________________
@flag_as_option
def pth_smH_plot_statsyst(args):
    obs_name = 'pth_smH'

    sm_crosssections = get_obs(obs_name, 'combination').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True)

    combination = differentials.scans.DifferentialSpectrum('combination',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='combination', statonly=False),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='combination')
        )
    combination.color = 9
    combination.no_overflow_label = True
    combination.draw_method = 'repr_point_with_vertical_bar_and_horizontal_bar'
    combination.set_sm(sm_crosssections)
    combination.title = 'Total unc.'
    combination.read()

    combination_statonly = differentials.scans.DifferentialSpectrum('combination_statonly',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='combination', statonly=True),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='combination')
        )
    combination_statonly.set_sm(sm_crosssections)
    combination_statonly.read()

    # Calculate the syst-only errors
    combination_histogram = combination.to_hist()
    statonly_histogram = combination_statonly.to_hist()

    systonly_histogram = deepcopy(combination_histogram)
    systonly_histogram.title = 'Syst. unc.'
    systonly_histogram.draw_method = 'repr_uncertainties_fully_filled_area'
    systonly_histogram_xs = combination.to_hist_xs()
    systonly_histogram_xs.title = 'Syst. unc.'
    systonly_histogram_xs.draw_method = 'repr_uncertainties_fully_filled_area'

    syst_err_up = []
    syst_err_down = []
    syst_err_up_xs = []
    syst_err_down_xs = []
    for i in xrange(combination_histogram.n_bins):
        up_tot = abs(combination_histogram.errs_up[i])
        down_tot = abs(combination_histogram.errs_down[i])
        symm_tot = 0.5*(up_tot+down_tot)
        symm_stat = 0.5*(abs(statonly_histogram.errs_up[i]) + abs(statonly_histogram.errs_down[i]))
        sm_crosssection = sm_crosssections[i]

        if symm_tot>symm_stat:
            symm_syst = sqrt(symm_tot**2-symm_stat**2)
        else:
            logging.warning(
                'For bin {0}, symm_tot={1} > symm_stat={2}. '
                'Taking sqrt(symm_stat**2-symm_tot**2) instead'
                .format(i, symm_tot, symm_stat)
                )
            symm_syst = sqrt(symm_stat**2-symm_tot**2)
        syst_err_up.append(up_tot * symm_syst/symm_tot)
        syst_err_down.append(down_tot * symm_syst/symm_tot)
        syst_err_up_xs.append(up_tot * symm_syst/symm_tot * sm_crosssection)
        syst_err_down_xs.append(down_tot * symm_syst/symm_tot * sm_crosssection)
    systonly_histogram.set_err_up(syst_err_up)
    systonly_histogram.set_err_down(syst_err_down)
    systonly_histogram_xs.set_err_up(syst_err_up_xs)
    systonly_histogram_xs.set_err_down(syst_err_down_xs)

    # combination_histogram_xs = combination.to_hist_xs()
    # statonly_histogram_xs = combination_statonly.to_hist_xs()
    # systonly_histogram_xs = deepcopy(combination_histogram_xs)
    # systonly_histogram_xs.title = 'Syst. only'
    # systonly_histogram_xs.set_err_up([ 
    #     sqrt(tot**2-stat**2) for tot, stat in zip(combination_histogram_xs.errs_up, statonly_histogram_xs.errs_up)
    #     ])
    # systonly_histogram_xs.set_err_down([ 
    #     sqrt(tot**2-stat**2) for tot, stat in zip(combination_histogram_xs.errs_down, statonly_histogram_xs.errs_down)
    #     ])
    # systonly_histogram_xs.draw_method = 'repr_horizontal_bar_and_narrow_fill'

    plot = differentials.plotting.plots.SpectraPlot('spectra_{0}_statsyst'.format(obs_name), [combination])
    plot.obsname = obs_name
    plot.obsunit = 'GeV'
    plot.add_top(systonly_histogram_xs, systonly_histogram_xs.draw_method, plot.leg)
    plot.add_bottom(systonly_histogram, systonly_histogram.draw_method)
    plot.draw()
    plot.wrapup()


#____________________________________________________________________
@flag_as_option
def pth_smH_plot_lumiscale(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
    containers = []

    lumi35 = prepare_container('lumi35', LatestPaths.ws_combined_smH, LatestPaths.scan_combination_pth_smH_asimov)
    lumi35.SMcrosssections = LatestBinning.obs_pth.crosssection_over_binwidth()
    lumi35.color = 1
    containers.append(lumi35)

    # lumi35_new = prepare_container('lumi35_new', LatestPaths.ws_combined_smH, 'out/Scan_pth_smH_Feb12_combination_asimov')
    # lumi35_new.SMcrosssections = LatestBinning.obs_pth.crosssection_over_binwidth()
    # lumi35_new.color = 4
    # containers.append(lumi35_new)

    lumi300 = prepare_container('lumi300', LatestPaths.ws_hgg_smH, 'out/Scan_pth_smH_Feb12_combination_lumiScale_asimov_1')
    lumi300.SMcrosssections = LatestBinning.obs_pth.crosssection_over_binwidth()
    lumi300.color = 2
    containers.append(lumi300)

    lumi3000 = prepare_container(
        'lumi3000', LatestPaths.ws_hgg_smH, 'out/Scan_pth_smH_Feb12_combination_lumiScale_asimov_1',
        scale_scans = 10.
        )
    lumi3000.SMcrosssections = LatestBinning.obs_pth.crosssection_over_binwidth()
    lumi3000.color = 4
    containers.append(lumi3000)

    # lumi35 = prepare_container('lumi35', LatestPaths.ws_hgg_smH, LatestPaths.scan_hgg_PTH)
    # lumi35.SMcrosssections = LatestBinning.obs_pth.crosssection_over_binwidth()
    # lumi35.color = 1
    # containers.append(lumi35)

    for container in containers:
        draw_parabolas(container)

    SM = prepare_SM_container(
        LatestBinning.obs_pth.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True),
        LatestBinning.obs_pth.binning
        )
    containers.append(SM)

    PlotCommands.PlotSpectraOnTwoPanel(
        'twoPanel_pthSpectrum' + ('_statsyst' if args.statsyst else ''),
        containers,
        xTitle = 'p_{T}^{H} (GeV)',
        yTitleTop = '#Delta#sigma/#Deltap_{T}^{H} (pb/GeV)',
        # 
        # yMinExternalTop = 0.0005,
        # yMaxExternalTop = 110.,
        )


#____________________________________________________________________
@flag_as_option
def pth_ggH_hbb_plot(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
    containers = []

    if args.asimov:
        combined_scanDir    = LatestPaths.scan_combined_PTH_ggH_asimov
        hgg_scanDir         = LatestPaths.scan_hgg_PTH_ggH_asimov
        hzz_scanDir         = LatestPaths.scan_hzz_PTH_ggH_asimov
        hbb_scanDir         = LatestPaths.scan_hbb_PTH_ggH_asimov
        combWithHbb_scanDir = LatestPaths.scan_combWithHbb_PTH_ggH_asimov
    else:
        combined_scanDir    = LatestPaths.scan_combined_PTH_ggH
        hgg_scanDir         = LatestPaths.scan_hgg_PTH_ggH
        hzz_scanDir         = LatestPaths.scan_hzz_PTH_ggH
        hbb_scanDir         = LatestPaths.scan_hbb_PTH_ggH
        combWithHbb_scanDir = LatestPaths.scan_combWithHbb_PTH_ggH


    hgg = prepare_container( 'hgg', LatestPaths.ws_hgg_ggH_xHfixed, hgg_scanDir )
    hgg.SMcrosssections = LatestBinning.obs_pth_ggH.crosssection_over_binwidth()
    containers.append(hgg)

    hzz = prepare_container( 'hzz', LatestPaths.ws_hzz_ggH_xHfixed, hzz_scanDir )
    hzz.SMcrosssections = LatestBinning.obs_pth_ggH_hzzBinning.crosssection_over_binwidth()
    containers.append(hzz)

    combination = prepare_container( 'combination', LatestPaths.ws_combined_ggH_xHfixed, combined_scanDir )
    combination.SMcrosssections = LatestBinning.obs_pth_ggH.crosssection_over_binwidth()
    containers.append(combination)

    hbb = prepare_container( 'hbb', LatestPaths.ws_hbb_ggH_xHfixed, hbb_scanDir )
    hbb.SMcrosssections = LatestBinning.obs_pth_ggH_hbbBinning.crosssection_over_binwidth()
    containers.append(hbb)

    combWithHbb = prepare_container( 'combWithHbb', LatestPaths.ws_combWithHbb_ggH_xHfixed, combWithHbb_scanDir )
    combWithHbb.SMcrosssections = LatestBinning.obs_pth_ggH_combWithHbbBinning.crosssection_over_binwidth()
    containers.append(combWithHbb)

    check_containers(containers)
    for container in containers:
        PlotCommands.WriteScansToTable(
            container,
            'pth_ggh_whbb',
            xTitle = 'p_{T}^{H} (GeV)',
            yTitle = '#Delta#sigma/#Delta p_{T}^{H} (pb/GeV)',
            lastBinIsOverflow = True,
            )

    SM = prepare_SM_container(
        LatestBinning.obs_pth_ggH_combWithHbbBinning.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True),
        LatestBinning.obs_pth_ggH_combWithHbbBinning.binning
        )
    containers.append(SM)

    l = PlotCommands.TLatexMultiPanel(
        lambda c: 1.0 - c.GetRightMargin() - 0.01,
        lambda c: 1.0 - c.GetTopMargin() - 0.14,
        '(non-ggH fixed to SM)'
        )
    l.SetTextSize(0.05)
    l.SetTextAlign(33)

    PlotCommands.PlotSpectraOnTwoPanel(
        'twoPanel_pth_ggH_hbb_Spectrum' + ( '_asimov' if args.asimov else '' ) + ('_statsyst' if args.statsyst else ''),
        containers,
        xTitle = 'p_{T}^{H} (GeV)',
        yTitleTop = '#Delta#sigma^{ggH}/#Deltap_{T}^{H} (pb/GeV)',
        # topPanelObjects = [ ( l, '' ) ],
        )


#____________________________________________________________________
@flag_as_option
def all_tables(args):
    Commands.DisableWarnings()
    pth_smH_tables(args)
    pth_ggH_tables(args)
    ptjet_tables(args)
    njets_tables(args)
    rapidity_tables(args)
    Commands.DisableWarnings(False)

@flag_as_option
def pth_smH_tables(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
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
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
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
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
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
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
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
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
    differentialTable = DifferentialTable.DifferentialTable(name='rapidity', last_bin_is_overflow=False)
    for decay_channel in ['hgg', 'hzz', 'combination']:
        statsyst = read_container(args, 'rapidity', rapidity_obs(decay_channel).crosssection_over_binwidth(), decay_channel, statonly=False )
        statonly = read_container(args, 'rapidity', rapidity_obs(decay_channel).crosssection_over_binwidth(), decay_channel, statonly=True )
        differentialTable.calculate_stat_syst(statsyst.name, statsyst, statonly)
    print differentialTable.repr_twiki_symm()
    print
    # differentialTable.do_xs = True
    # print differentialTable.repr_twiki_symm()
