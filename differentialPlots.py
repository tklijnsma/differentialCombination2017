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
from differentialTools import *

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
    SM.color = 14
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


def read_containers(args, obs_name, SMcrosssections, decay_channels=['combination', 'hgg', 'hzz']):
    containers = {}
    for statonly in [True, False]:
        containers[statonly] = []
        for decay_channel in decay_channels:
            ws = LatestPathsGetters.get_ws(obs_name, args, decay_channel=decay_channel)
            scan = LatestPathsGetters.get_scan(obs_name, args, decay_channel=decay_channel, statonly=statonly)
            container = prepare_container(decay_channel, ws, scan)
            container.name += '_statonly' if statonly else ''
            container.SMcrosssections = SMcrosssections[decay_channel]
            container.statonly = statonly
            container.suppress_text = statonly
            container.error_line = not(statonly)
            containers[statonly].append(container)
    # Sort to full first, statonly second
    containers = containers[False] + containers[True]

    for container in containers:
        draw_parabolas(container)
    return containers

########################################
# Plotting
########################################

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

    if args.statsyst:
        SMcrosssections = {
            'hgg' : LatestBinning.obs_pth.crosssection_over_binwidth(),
            'hzz' : LatestBinning.obs_pth_hzzBinning.crosssection_over_binwidth(),
            'combination' : LatestBinning.obs_pth.crosssection_over_binwidth(),
            }
        containers = read_containers(args, 'pth_smH', SMcrosssections)

    else:
        containers = []

        hgg = prepare_container('hgg', LatestPaths.ws_hgg_smH, LatestPaths.scan_hgg_PTH)
        hgg.SMcrosssections = LatestBinning.obs_pth.crosssection_over_binwidth()
        containers.append(hgg)

        hzz = prepare_container('hzz', LatestPaths.ws_hzz_smH, LatestPaths.scan_hzz_PTH)
        hzz.SMcrosssections = LatestBinning.obs_pth_hzzBinning.crosssection_over_binwidth()
        containers.append(hzz)

        combination = prepare_container('combination', LatestPaths.ws_combined_smH, LatestPaths.scan_combined_PTH)
        combination.SMcrosssections = LatestBinning.obs_pth.crosssection_over_binwidth()
        containers.append(combination)

        for container in containers:
            PlotCommands.WriteScansToTable(
                container,
                'pth',
                xTitle = 'p_{T}^{H} (GeV)',
                yTitle = '#Delta#sigma/#Delta p_{T}^{H} (pb/GeV)',
                lastBinIsOverflow = True,
                )

    SM = prepare_SM_container(
        LatestBinning.obs_pth.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True),
        LatestBinning.obs_pth.binning
        )
    containers.append(SM)

    PlotCommands.PlotSpectraOnTwoPanel(
        'twoPanel_pthSpectrum',
        containers,
        xTitle = 'p_{T}^{H} (GeV)',
        yTitleTop = '#Delta#sigma/#Deltap_{T}^{H} (pb/GeV)',
        # 
        # yMinExternalTop = 0.0005,
        # yMaxExternalTop = 110.,
        )

#____________________________________________________________________
@flag_as_option
def pth_ggH_plot(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )

    if args.statsyst:
        SMcrosssections = {
            'hgg' : LatestBinning.obs_pth_ggH.crosssection_over_binwidth(),
            'hzz' : LatestBinning.obs_pth_ggH_hzzBinning.crosssection_over_binwidth(),
            'combination' : LatestBinning.obs_pth_ggH.crosssection_over_binwidth(),
            }
        containers = read_containers(args, 'pth_ggH', SMcrosssections)

    else:
        containers = []

        hgg = prepare_container( 'hgg', LatestPaths.ws_hgg_ggH_xHfixed, LatestPaths.scan_hgg_PTH_ggH )
        hgg.SMcrosssections = LatestBinning.obs_pth_ggH.crosssection_over_binwidth()
        containers.append(hgg)

        hzz         = prepare_container( 'hzz', LatestPaths.ws_hzz_ggH_xHfixed, LatestPaths.scan_hzz_PTH_ggH )
        hzz.SMcrosssections = LatestBinning.obs_pth_ggH_hzzBinning.crosssection_over_binwidth()
        containers.append(hzz)

        combination = prepare_container( 'combination', LatestPaths.ws_combined_ggH_xHfixed, LatestPaths.scan_combined_PTH_ggH )
        combination.SMcrosssections = LatestBinning.obs_pth_ggH.crosssection_over_binwidth()
        containers.append(combination)

        for container in containers:
            PlotCommands.WriteScansToTable(
                container,
                'pth_ggh',
                xTitle = 'p_{T}^{H} (GeV)',
                yTitle = '#Delta#sigma/#Delta p_{T}^{H} (pb/GeV)',
                lastBinIsOverflow = True,
                )

    SM = prepare_SM_container(
        LatestBinning.obs_pth_ggH.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True),
        LatestBinning.obs_pth_ggH.binning
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
        'twoPanel_pth_ggH_Spectrum',
        containers,
        xTitle = 'p_{T}^{H} (GeV)',
        yTitleTop = '#Delta#sigma^{ggH}/#Deltap_{T}^{H} (pb/GeV)',
        topPanelObjects = [ ( l, '' ) ],
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
        'twoPanel_pth_ggH_hbb_Spectrum' + ( '_asimov' if args.asimov else '' ),
        containers,
        xTitle = 'p_{T}^{H} (GeV)',
        yTitleTop = '#Delta#sigma^{ggH}/#Deltap_{T}^{H} (pb/GeV)',
        # topPanelObjects = [ ( l, '' ) ],
        )

#____________________________________________________________________
@flag_as_option
def njets_plot(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )

    if args.statsyst:
        SMcrosssections = {
            'hgg' : LatestBinning.obs_njets.crosssection_over_binwidth(),
            'hzz' : LatestBinning.obs_njets_hzzBinning.crosssection_over_binwidth(),
            'combination' : LatestBinning.obs_njets.crosssection_over_binwidth(),
            }
        containers = read_containers(args, 'njets', SMcrosssections)

    else:
        containers = []
        hgg = prepare_container( 'hgg', LatestPaths.ws_hgg_smH_NJ,LatestPaths.scan_hgg_NJ )
        hgg.SMcrosssections = LatestBinning.obs_njets.crosssection_over_binwidth()
        containers.append(hgg)

        hzz = prepare_container( 'hzz', LatestPaths.ws_hzz_smH_NJ,LatestPaths.scan_hzz_NJ )
        hzz.SMcrosssections = LatestBinning.obs_njets_hzzBinning.crosssection_over_binwidth()
        containers.append(hzz)

        combination = prepare_container( 'combination', LatestPaths.ws_combined_smH_NJ,LatestPaths.scan_combined_NJ )
        combination.SMcrosssections = LatestBinning.obs_njets.crosssection_over_binwidth()
        containers.append(combination)

        check_containers(containers)
        for container in containers:
            PlotCommands.WriteScansToTable(
                container,
                'njets',
                xTitle = 'N_{jets}',
                yTitle = '#Delta#sigma/#Delta N_{jets} (pb)',
                lastBinIsOverflow = False,
                )

    SM = prepare_SM_container(
        LatestBinning.obs_njets.crosssection_over_binwidth(),
        LatestBinning.obs_njets.binning
        )
    containers.append(SM)

    PlotCommands.PlotSpectraOnTwoPanel(
        'twoPanel_nJetsSpectrum',
        containers,
        xTitle = 'N_{jets}',
        # yMinLimit = 0.07,
        # yMaxExternalTop = 500,
        yTitleTop = '#Delta#sigma/#DeltaN_{jets} (pb)',
        )

#____________________________________________________________________
@flag_as_option
def ptjet_plot(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )

    if args.statsyst:
        SMcrosssections = {
            'hgg' : LatestBinning.obs_ptjet.crosssection_over_binwidth(),
            'hzz' : LatestBinning.obs_ptjet_hzzBinning.crosssection_over_binwidth(),
            'combination' : LatestBinning.obs_ptjet.crosssection_over_binwidth(),
            }
        containers = read_containers(args, 'ptjet', SMcrosssections)

    else:
        containers = []
        hgg = prepare_container('hgg', LatestPaths.ws_hgg_smH_PTJ, LatestPaths.scan_hgg_PTJ_asimov)
        hgg.SMcrosssections = LatestBinning.obs_ptjet.crosssection_over_binwidth()
        containers.append(hgg)

        hzz = prepare_container('hzz', LatestPaths.ws_hzz_smH_PTJ, LatestPaths.scan_hzz_PTJ_asimov)
        hzz.SMcrosssections = LatestBinning.obs_ptjet_hzzBinning.crosssection_over_binwidth()
        containers.append(hzz)

        combination = prepare_container('combination', LatestPaths.ws_combined_smH_PTJ, LatestPaths.scan_combined_PTJ_asimov)
        combination.SMcrosssections = LatestBinning.obs_ptjet.crosssection_over_binwidth()
        containers.append(combination)

    Commands.Warning( 'Skipping first bin for ptjet (should be the underflow)' )
    for container in containers:
        container.POIs = container.POIs[1:]
        container.Scans = container.Scans[1:]
        container.SMcrosssections = container.SMcrosssections[1:]

    if not args.statsyst:
        check_containers(containers)
        for container in containers:
            PlotCommands.WriteScansToTable(
                container,
                'ptjet',
                xTitle = 'p_{T}^{jet} (GeV)',
                yTitle = '#Delta#sigma/#Delta p_{T}^{jet} (pb/GeV)',
                lastBinIsOverflow = True,
                )

    SM = prepare_SM_container(
        LatestBinning.obs_ptjet.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True),
        LatestBinning.obs_ptjet.binning
        )
    containers.append(SM)

    PlotCommands.PlotSpectraOnTwoPanel(
        'twoPanel_ptjetSpectrum',
        containers,
        xTitle = 'p_{T}^{jet} (GeV)',
        yTitleTop = '#Delta#sigma/#Deltap_{T}^{jet} (pb/GeV)',
        # yMinLimit = 0.07,
        yMaxExternalTop = 10,
        xMinExternal = 30.0,
        # yMinLimit    = 0.1
        )

#____________________________________________________________________
@flag_as_option
def rapidity_plot(args):
    TheoryCommands.SetPlotDir( 'plots_{0}'.format(datestr) )

    if args.statsyst:
        SMcrosssections = {
            'hgg' : LatestBinning.obs_yh.crosssection_over_binwidth(),
            'hzz' : LatestBinning.obs_yh.crosssection_over_binwidth(),
            'combination' : LatestBinning.obs_yh.crosssection_over_binwidth(),
            }
        containers = read_containers(args, 'rapidity', SMcrosssections)

    else:
        containers = []
        hgg = prepare_container('hgg', LatestPaths.ws_hgg_smH_YH, LatestPaths.scan_hgg_YH_asimov)
        hgg.SMcrosssections = LatestBinning.obs_yh.crosssection_over_binwidth()
        containers.append(hgg)

        hzz = prepare_container('hzz', LatestPaths.ws_hzz_smH_YH, LatestPaths.scan_hzz_YH_asimov)
        hzz.SMcrosssections = LatestBinning.obs_yh.crosssection_over_binwidth()
        containers.append(hzz)

        combination = prepare_container('combination', LatestPaths.ws_combined_smH_YH, LatestPaths.scan_combined_YH_asimov)
        combination.SMcrosssections = LatestBinning.obs_yh.crosssection_over_binwidth()
        containers.append(combination)

        check_containers(containers)
        for container in containers:
            PlotCommands.WriteScansToTable(
                container,
                'rapidity',
                xTitle = '|y_{H}|',
                yTitle = '#Delta#sigma/#Delta|y_{H}| (pb)',
                lastBinIsOverflow = False,
                )

    SM = prepare_SM_container(
        LatestBinning.obs_yh.crosssection_over_binwidth(),
        LatestBinning.obs_yh.binning
        )
    containers.append(SM)

    PlotCommands.PlotSpectraOnTwoPanel(
        'twoPanel_rapiditySpectrum',
        containers,
        xTitle = '|y_{H}|',
        yTitleTop = '#Delta#sigma/#Delta|y_{H}| (pb)',
        # yMinLimit = 0.07,
        # yMaxExternalTop = 500
        lastBinIsNotOverflow=True,
        )

