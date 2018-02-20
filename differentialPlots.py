#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, sys, re
from os.path import *
from glob import glob
from copy import deepcopy

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
datestr = strftime( '%b%d' )


########################################
# Helpers
########################################

def check_containers(containers):
    for container in containers:
        if not len(container.POIs) == len(container.SMcrosssections):
            Commands.ThrowError(
                'For container {2}, found {0} POIs, but {1} SM cross sections; something is misaligned.'.format(
                    len(container.POIs), len(container.SMcrosssections), container.name )
                + '\n  POIs:  {0}'.format(container.POIs)
                + '\n  SM xs: {0}'.format(container.SMcrosssections)
                )

def draw_parabolas(container):
    PhysicsCommands.BasicDrawScanResults(container.POIs, container.Scans, name=container.name)

def prepare_container(
        name,
        ws, scandir,
        title=None,
        verbose=False,
        draw_parabolas=False,
        scale_scans=None
        ):
    POIs = Commands.ListPOIs( ws )
    POIs.sort( key=Commands.POIsorter )
    if verbose: print 'Sorted POIs:', POIs

    scans = PhysicsCommands.GetScanResults(
        POIs,
        scandir,
        # pattern = pattern,
        filterNegatives = True
        )
    if scale_scans:
        for scan in scans:
            for i in xrange(len(scan)):
                scan[i] = ( scan[i][0], scale_scans*scan[i][1] )

    if draw_parabolas:
        draw_parabolas(container)

    container = Container(name=name)
    if title is None:
        title = name
    container.title = title

    container.POIs = POIs
    container.Scans = scans

    if name == 'hgg':
        container.color = 2
        container.title = 'H#rightarrow#gamma#gamma'
    if name == 'hzz':
        container.color = 4
        container.title = 'H#rightarrowZZ'
    if name == 'combination':
        container.color = 1
        container.title = 'Combination'
    if name == 'hbb':
        container.color = 8
        container.title = 'H#rightarrowbb'
    if name == 'combWithHbb':
        container.color = 14
        container.title = 'Comb. with H#rightarrowbb'

    return container

def prepare_SM_container(crosssection, binBoundaries):
    SM = Container()
    SM.name = 'SM'
    SM.title = 'SM'
    SM.color = 16
    SM.crosssection = crosssection
    SM.binBoundaries = binBoundaries
    SM.ratios = [ 1.0 for i in xrange(len(crosssection)) ]
    return SM

def get_scan_dir(args, obs_name):
    decay_channel = get_decay_channel_tag(args)
    key = 'scan_{0}_{1}'.format(decay_channel, obs_name)
    if args.statonly:
        key += '_statonly'
    if args.asimov:
        key += '_asimov'
    return getattr(LatestPaths, key)

def read_container_base(args, obs_name, SMcrosssections, decay_channel, statonly=False):
    ws = LatestPathsGetters.get_ws(obs_name, args, decay_channel=decay_channel)
    scan = LatestPathsGetters.get_scan(obs_name, args, decay_channel=decay_channel, statonly=statonly)
    container = prepare_container(decay_channel, ws, scan)
    container.name += '_statonly' if statonly else ''
    container.SMcrosssections = SMcrosssections
    container.statonly = statonly
    container.suppress_text = statonly
    container.error_line = not(statonly)
    return container

def read_container(args, obs_name, decay_channel=['combination']):
    SMcrosssections = get_obs(obs_name, decay_channel).crosssection_over_binwidth()
    container = read_container_base(args, obs_name, SMcrosssections, decay_channel)
    return container

def read_containers_statsyst(args, obs_name, decay_channels=['combination']):
    containers = []
    for statonly in [False, True]:
        for decay_channel in decay_channels:
            container = read_container_base(
                args,
                obs_name,
                get_obs(obs_name, decay_channel).crosssection_over_binwidth(),
                decay_channel,
                statonly
                )
            container.suppress_text = False
            container.title = 'stat.' + ( ' #oplus syst.' if not statonly else '' )
            containers.append(container)
    # for container in containers:
    #     draw_parabolas(container)
    return containers


########################################
# Plotting
########################################

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

    differentials.scans.plot_spectra(
        'spectra_{0}'.format(obs_name),
        [ hgg, hzz, combination ],
        obs_name, obs_title='p_{T}', obs_unit='GeV'
        )

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

    plot = differentials.scans.plot_spectra(
        'spectra_{0}'.format(obs_name),
        [ hgg, hzz, combination ],
        obs_name, obs_title='p_{T}^{ggH}', obs_unit='GeV', inplace=False
        )

    l = differentials.plotting.pywrappers.Latex(
        lambda c: 1.0 - c.GetRightMargin() - 0.01,
        lambda c: 1.0 - c.GetTopMargin() - 0.14,
        'gluon fusion cross section'
        )
    l.SetNDC()
    l.SetTextSize(0.05)
    l.SetTextAlign(33)
    plot.add_top(l, '')

    plot.draw()


#____________________________________________________________________
@flag_as_option
def njets_plot(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )
    obs_name = 'njets'
    
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

    differentials.scans.plot_spectra(
        'spectra_{0}'.format(obs_name),
        [ hgg, hzz, combination ],
        obs_name, obs_title='N_{jets}'
        )

#____________________________________________________________________
@flag_as_option
def ptjet_plot(args):
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

    differentials.scans.plot_spectra(
        'spectra_{0}'.format(obs_name),
        [ hgg, hzz, combination ],
        obs_name, obs_title='p_{T}^{jet}', obs_unit='GeV'
        )


#____________________________________________________________________
@flag_as_option
def rapidity_plot(args):
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

    differentials.scans.plot_spectra(
        'spectra_{0}'.format(obs_name),
        [ hgg, hzz, combination ],
        obs_name, obs_title='|y_{H}|'
        )


########################################
# Other plots
########################################

#____________________________________________________________________
@flag_as_option
def pth_smH_plot_statsyst(args):
    obs_name = 'pth_smH'

    combination = differentials.scans.DifferentialSpectrum('combination',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='combination', statonly=False),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='combination')
        )
    combination.color = 9
    combination.no_overflow_label = True
    combination.draw_method = 'repr_point_with_vertical_bar'
    combination.set_sm( get_obs(obs_name, 'combination').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )

    combination_statonly = differentials.scans.DifferentialSpectrum('combination_statonly',
        scandir = LatestPathsGetters.get_scan(obs_name, args, decay_channel='combination', statonly=True),
        datacard = LatestPathsGetters.get_ws(obs_name, args, decay_channel='combination')
        )
    combination_statonly.color = 9
    combination_statonly.no_overflow_label = True
    combination_statonly.draw_method = 'repr_horizontal_bar_and_narrow_fill'
    combination_statonly.set_sm( get_obs(obs_name, 'combination').crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True) )

    differentials.scans.plot_spectra(
        'spectra_{0}_statsyst'.format(obs_name),
        [ combination, combination_statonly ],
        obs_name, obs_title='p_{T}', obs_unit='GeV'
        )


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

