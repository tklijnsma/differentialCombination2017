#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

from OptionHandler import flag_as_option

import LatestPaths
import LatestBinning

import differentials
import differentialutils
import logging
from collections import namedtuple

from time import strftime
datestr = strftime( '%b%d' )

import ROOT
from array import array
import sys
import copy

########################################
# Main
########################################

x_coupling = 'kappac'
y_coupling = 'kappab'
titles = { 'kappac': '#kappa_{c}', 'kappab' : '#kappa_{b}' }

yukawa_x_min = -35
yukawa_x_max = 35
yukawa_y_min = -11
yukawa_y_max = 12

couplingdependentBRs_x_min = -8.
couplingdependentBRs_x_max = 8.
couplingdependentBRs_y_min = -2.
couplingdependentBRs_y_max = 2.

floatingBRs_x_min = -75.
floatingBRs_x_max = 75.
floatingBRs_y_min = -30.
floatingBRs_y_max = 40.


@flag_as_option
def all_plots_Yukawa(args):
    multicont_Yukawa(args)
    points_on_contour_Yukawa(args)
    args = differentialutils.force_asimov(args)
    multicont_Yukawa(args)
    multicont_Yukawa_highLumi(args)
    multicont_Yukawa_profiledTotalXS(args)

scandict_G = differentials.core.AttrDict()
# scandict_G.G0A = LatestPaths.scan.yukawa.reweighted.asimov.combination
scandict_G.G0A = 'out/Scan_Yukawa_May11_combination_G0A_asimov'
# scandict_G.G0B = 'out/Scan_Yukawa_May11_combination_G0B_asimov_0'
scandict_G.G1A = 'out/Scan_Yukawa_May11_combination_G1A_asimov'
# scandict_G.G1B = 'out/Scan_Yukawa_May11_combination_G1B_asimov_0'
# scandict_G.G2A = 'out/Scan_Yukawa_May11_combination_G2A_asimov'
scandict_G.G1BKV = 'out/Scan_Yukawa_May11_combination_G1BKV_asimov_0'

scandict_G.G0B = 'out/Scan_Yukawa_May14_combination_G0B_asimov'
scandict_G.G1B = 'out/Scan_Yukawa_May14_combination_G1B_asimov'
scandict_G.G2A = 'out/Scan_Yukawa_May14_combination_G2A_asimov'

scandict_G.G1A_unreweighted = 'out/Scan_Yukawa_May14_combination_G1A_unreweighted_asimov/'
scandict_G.G1B_reweighted = 'out/Scan_Yukawa_May14_combination_G1B_reweighted_asimov'

# scandict_G.G0B_reweighted = 'out/Scan_Yukawa_May16_combination_G0B_reweighted_asimov'
scandict_G.G0B_reweighted = 'out/Scan_Yukawa_May17_combination_G0B_reweighted_asimov'

# approval = differentials.core.AttrDict()
# approval.scalingbbH_couplingdependentBRs = differentials.core.AttrDict()
# approval.scalingbbH_couplingdependentBRs.data = 'out/Scan_Yukawa_May24_combination_scalingbbH_couplingdependentBRs'
# approval.scalingbbH_couplingdependentBRs.asimov = 'out/Scan_Yukawa_May22_combination_scalingbbH_couplingdependentBRs_asimov'
# approval.scalingbbH_fixedSMBRs = 'out/Scan_Yukawa_May18_combination_scalingbbH_asimov'


approval = differentials.core.AttrDict.create_tree(['fixedBRs', 'couplingdependentBRs', 'floatingBRs'], ['asimov', 'observed'])

approval.couplingdependentBRs.asimov.combination   = 'out/Scan_Yukawa_Jun05_combination_NONscalingbbH_couplingdependentBRs_asimov'
approval.couplingdependentBRs.observed.combination = 'out/Scan_Yukawa_May30_combination_NONscalingbbH_couplingdependentBRs'
approval.couplingdependentBRs.observed.hgg         = 'out/Scan_Yukawa_Jun07_hgg_NONscalingbbH_couplingdependentBRs'
approval.couplingdependentBRs.observed.hzz         = 'out/Scan_Yukawa_Jun11_hzz_NONscalingbbH_couplingdependentBRs'

# approval.floatingBRs.asimov.combination   = 'out/Scan_Yukawa_May24_combination_NONscalingbbH_floatingBRs_asimov'
approval.floatingBRs.asimov.combination   = 'out/Scan_Yukawa_Jun09_combination_NONscalingbbH_floatingBRs_asimov'
approval.floatingBRs.observed.combination = 'out/Scan_Yukawa_May30_combination_NONscalingbbH_floatingBRs'
approval.floatingBRs.observed.hgg         = 'out/Scan_Yukawa_Jun07_hgg_NONscalingbbH_floatingBRs'
approval.floatingBRs.observed.hzz         = 'out/Scan_Yukawa_Jun11_hzz_NONscalingbbH_floatingBRs_0'


#____________________________________________________________________
def latest_couplingdependentBRs(args, decay_channel=None, asimov=None, splined=False):
    if not(asimov is None): args = differentialutils.force_asimov(args, asimov)
    if decay_channel is None: decay_channel = differentialutils.get_decay_channel_tag(args)
    scandir = approval.couplingdependentBRs['asimov' if args.asimov else 'observed'][decay_channel]

    scan = differentials.scans.Scan2D(
        'kbkc_{0}'.format(decay_channel), x_coupling, y_coupling,
        scandir = scandir
        )
    scan.title = differentials.core.standard_titles.get(decay_channel, decay_channel)
    scan.color = 1
    scan.read()
    if splined: return spline_couplingdependentBRs(scan)
    return scan

def spline_couplingdependentBRs(scan):
    # Careful with ranges; too far outside causes noise
    spline = scan.to_spline(
        x_min = -8.,
        x_max = 8.,
        y_min = -1.8,
        y_max = 2.1,
        )
    spline.disallow_negativity = False
    hist = spline.to_hist(nx=300, ny=300)
    # hist = spline.to_hist(nx=600, ny=600)
    hist.color = scan.color
    hist.name  = scan.name + '_splined'
    hist.title = scan.title
    return hist

def latest_floatingBRs(args, decay_channel=None, asimov=None, splined=False):
    if not(asimov is None): args = differentialutils.force_asimov(args, asimov)
    if not(decay_channel is None): args = differentialutils.set_one_decay_channel(args, decay_channel)
    decay_channel = differentialutils.get_decay_channel_tag(args)
    scandir = approval.floatingBRs['asimov' if args.asimov else 'observed'][decay_channel]

    scan = differentials.scans.Scan2D(
        'kbkc_{0}'.format(decay_channel), x_coupling, y_coupling,
        scandir = scandir
        )
    scan.title = differentials.core.standard_titles.get(decay_channel, decay_channel)
    scan.color = 1
    scan.read()
    if splined: return spline_floatingBRs(args, scan)
    return scan

def spline_floatingBRs(args, scan):
    # Careful with ranges; too far outside causes noise
    x_min = -75.
    x_max = 75.
    y_min = -30.
    y_max = 40.
    eps = 2.2
    deltaNLL_cutoff = 30.
    if args.hzz:
        x_min = -35.
        x_max = 40.
        y_min = -9.
        y_max = 19.
        eps   = 2.5
        deltaNLL_cutoff = 10.
    spline = scan.to_spline(
        x_min = x_min,
        x_max = x_max,
        y_min = y_min,
        y_max = y_max,
        eps = eps,
        deltaNLL_cutoff = deltaNLL_cutoff
        )
    spline.disallow_negativity = False
    hist = spline.to_hist(nx=300, ny=300)
    hist.color = scan.color
    hist.name  = scan.name + '_splined'
    hist.title = scan.title
    return hist

def latest_scenario3(args, splined=False):
    if args.asimov:
        scandir = 'out/Scan_Yukawa_May11_combination_G0A_asimov'
    else:
        scandir = LatestPaths.scan.yukawa.reweighted.observed.combination

    scan = differentials.scans.Scan2D(
        'scenario3', x_coupling, y_coupling,
        scandir = scandir
        )
    scan.title = differentials.core.standard_titles['hgg'] + ' + ' + differentials.core.standard_titles['hzz']
    scan.color = 1
    scan.read()

    if not splined: return scan

    x_min = -18.
    x_max = 18.
    y_min = -5.
    y_max = 6.5
    eps = 2.2
    deltaNLL_cutoff = 30.
    spline = scan.to_spline(
        x_min = x_min,
        x_max = x_max,
        y_min = y_min,
        y_max = y_max,
        # eps = eps,
        # deltaNLL_cutoff = deltaNLL_cutoff
        )
    spline.disallow_negativity = False
    hist = spline.to_hist(nx=180, ny=180, x_min=1.1*x_min, x_max=1.1*x_max, y_min=1.2*y_min, y_max=1.3*y_max)
    hist.color = scan.color
    hist.name  = scan.name + '_splined'
    hist.title = scan.title

    return hist

#____________________________________________________________________
@flag_as_option
def multicont_Yukawa_NONscalingbbH_floatingBRs(args):
    scans = []
    scans.append(latest_floatingBRs(args, decay_channel='combination', splined=True))
    if not args.asimov:
        # hzz = latest_floatingBRs(args, decay_channel='hzz', splined=True)
        hzz = latest_floatingBRs(args, decay_channel='hzz', splined=False)
        hzz.color = differentials.core.safe_colors.blue
        scans.append(hzz)
        hgg = latest_floatingBRs(args, decay_channel='hgg', splined=True)
        hgg.color = differentials.core.safe_colors.red
        scans.append(hgg)

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_floatingBRs' + ('_asimov' if args.asimov else ''),
        scans,
        x_title = differentials.core.standard_titles['kappac'], y_title = differentials.core.standard_titles['kappab'],
        x_min=floatingBRs_x_min, x_max=floatingBRs_x_max, y_min=floatingBRs_y_min, y_max=floatingBRs_y_max
        )
    plot.draw(wait=True)
    plot.add_BR_floating_text()
    plot.wrapup()


@flag_as_option
def multicont_Yukawa_NONscalingbbH_couplingdependentBRs(args):
    do_spline = True
    scans = []
    scans.append(latest_couplingdependentBRs(args, decay_channel='combination', splined=do_spline))
    if not args.asimov:
        hzz = latest_couplingdependentBRs(args, decay_channel='hzz', splined=do_spline)
        hzz.color = differentials.core.safe_colors.blue
        scans.append(hzz)
        hgg = latest_couplingdependentBRs(args, decay_channel='hgg', splined=do_spline)
        hgg.color = differentials.core.safe_colors.red
        scans.append(hgg)

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_couplingdependentBRs' + ('_asimov' if args.asimov else ''),
        scans,
        x_title = differentials.core.standard_titles['kappac'], y_title = differentials.core.standard_titles['kappab'],
        x_min=couplingdependentBRs_x_min, x_max=couplingdependentBRs_x_max, y_min=couplingdependentBRs_y_min, y_max=couplingdependentBRs_y_max
        )
    plot.draw(wait=True)
    plot.add_BR_parametrized_text()
    plot.wrapup()

@flag_as_option
def multicont_Yukawa_NONscalingbbH(args):
    args = differentialutils.force_asimov(args)
    hist = latest_scenario3(args, splined=True)

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_scenario3' + ('_asimov' if args.asimov else ''),
        [ hist ],
        x_title = differentials.core.standard_titles['kappac'], y_title = differentials.core.standard_titles['kappab'],
        # x_min=1.1*x_min, x_max=1.1*x_max, y_min=1.1*y_min, y_max=1.1*y_max
        )
    plot.legend.set(
        x1 = lambda c: c.GetLeftMargin() + 0.02,
        x2 = lambda c: c.GetLeftMargin() + 0.30,
        y1 = lambda c: 1. - c.GetTopMargin() - 0.11,
        y2 = lambda c: 1. - c.GetTopMargin() - 0.01,
        )
    plot.draw()




#____________________________________________________________________
@flag_as_option
def onedimscans_Yukawa_scalingbbH_couplingdependentBRs(args):
    args = differentialutils.set_one_decay_channel(args, 'combination')
    obs2D = latest_couplingdependentBRs(args, splined=True)
    exp2D = latest_couplingdependentBRs(args, splined=True, asimov=True)
    onedimscans_Yukawa('kappab', exp2D, obs2D, x_min=-4.0, x_max=4.0, apply_smoothing=True, tag='couplingdependentBRs')
    onedimscans_Yukawa('kappac', exp2D, obs2D, x_min=-15., x_max=15., apply_smoothing=True, tag='couplingdependentBRs')

@flag_as_option
def onedimscans_Yukawa_scalingbbH_floatingBRs(args):
    args = differentialutils.set_one_decay_channel(args, 'combination')
    obs2D = latest_floatingBRs(args, splined=True)
    exp2D = latest_floatingBRs(args, splined=True, asimov=True)
    onedimscans_Yukawa('kappab', exp2D, obs2D, x_min=-30.0, x_max=30.0, tag='floatingBRs')
    onedimscans_Yukawa('kappac', exp2D, obs2D, x_min=-75., x_max=75., tag='floatingBRs')

@flag_as_option
def onedimscans_Yukawa_scalingbbH_scenario3(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    # obs2D = latest_floatingBRs(args, splined=True)
    exp2D = latest_scenario3(args, splined=True)
    onedimscans_Yukawa('kappab', exp2D, x_min=-10.0, x_max=10.0, tag='scenario3')
    onedimscans_Yukawa('kappac', exp2D, x_min=-20., x_max=20., tag='scenario3')

def f2digits(number):
    if abs(number) >= 10.:
        n_decimals = 0
    elif abs(number) >= 1.:
        n_decimals = 1
    else:
        n_decimals = 2
    return '{0:.{n_decimals}f}'.format(number, n_decimals=n_decimals)

def onedimscans_Yukawa(kappa, exp2D, obs2D=None, x_min=-10., x_max=10., apply_smoothing=False, tag=None):
    do_95percent_CL = True
    if not obs2D is None:
        onedimscanner_obs = differentials.onedimscanner.OneDimScanner(
            obs2D,
            'kappac', 'kappab',
            do_95percent_CL = do_95percent_CL
            )
        obs1D = onedimscanner_obs.get_1d(kappa)
        obs1D.draw_style = 'repr_smooth_line'
        obs1D.title = '{0} observed; ({1} - {2})  ({3}% CL)'.format(
            differentials.core.standard_titles[kappa],
            f2digits(obs1D.unc.left_bound), f2digits(obs1D.unc.right_bound),
            '95' if do_95percent_CL else '68'
            )
        obs1D.color = 1
        if apply_smoothing: obs1D.smooth_y(10)

    onedimscanner_exp = differentials.onedimscanner.OneDimScanner(
        exp2D,
        'kappac', 'kappab',
        do_95percent_CL = do_95percent_CL
        )
    exp1D = onedimscanner_exp.get_1d(kappa)
    exp1D.draw_style = 'repr_dashed_line'
    exp1D.title = '{0} expected; ({1} - {2})  ({3}% CL)'.format(
        differentials.core.standard_titles[kappa],
        f2digits(exp1D.unc.left_bound), f2digits(exp1D.unc.right_bound),
        '95' if do_95percent_CL else '68'
        )
    exp1D.color = 1
    if apply_smoothing: exp1D.smooth_y(10)

    plotname = 'onekappascan_kbkc'
    if not(tag is None): plotname += '_' + tag
    plotname += '_' + kappa

    differentials.plotting.canvas.c.resize_temporarily(850, 800)
    plot = differentials.plotting.plots.MultiScanPlot(plotname)
    plot.do_95percent_CL = do_95percent_CL

    plot.y_max    = 7.0
    plot.y_cutoff = 5.0

    if not obs2D is None: plot.manual_graphs.append(obs1D)
    plot.manual_graphs.append(exp1D)
    plot.x_title = differentials.core.standard_titles[kappa]
    plot.x_min = x_min
    plot.x_max = x_max
    plot.leg.SetNColumns(1)
    plot.leg._y1 = lambda c: 1. - c.GetTopMargin() - 0.20
    plot.draw()
    plot.wrapup()

    texstr_var = lambda cmdname, value: '\\newcommand{{\\{0}}}{{{1}}}'.format(cmdname, value)
    print texstr_var('{0}LeftAsimov'.format(kappa), f2digits(exp1D.unc.left_bound))
    print texstr_var('{0}RightAsimov'.format(kappa), f2digits(exp1D.unc.right_bound))
    if not obs2D is None:
        print texstr_var('{0}LeftObserved'.format(kappa), f2digits(obs1D.unc.left_bound))
        print texstr_var('{0}RightObserved'.format(kappa), f2digits(obs1D.unc.right_bound))



#____________________________________________________________________

@flag_as_option
def multicont_Yukawa_scalingbbH_couplingdependentBRs(args):
    if args.asimov:
        scandir = approval.scalingbbH_couplingdependentBRs.asimov
    else:
        scandir = approval.scalingbbH_couplingdependentBRs.data

    scalingbbH = differentials.scans.Scan2D(
        'scalingbbH', x_coupling, y_coupling,
        scandir = scandir
        )
    scalingbbH.title = 'Combination (incl. bbH / BR(#vec{#kappa}))'
    scalingbbH.color = 1
    scalingbbH.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_scalingbbH_couplingdependentBRs' + ('_asimov' if args.asimov else ''),
        [ scalingbbH ],
        x_min=yukawa_x_min, x_max=yukawa_x_max, y_min=yukawa_y_min, y_max=yukawa_y_max
        )
    plot.draw()


@flag_as_option
def multicont_Yukawa_scalingbbH_floatingBRs(args):
    if args.asimov:
        scandir = 'out/Scan_Yukawa_May23_combination_scalingbbH_floatingBRs_asimov_1'
    else:
        scandir = 'out/Scan_Yukawa_May24_combination_scalingbbH_floatingBRs'

    scalingbbH = differentials.scans.Scan2D(
        'scalingbbH', x_coupling, y_coupling,
        scandir = scandir
        )
    scalingbbH.title = 'Combination (incl. bbH / BRs free)'
    scalingbbH.color = 1
    scalingbbH.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_scalingbbH_floatingBRs' + ('_asimov' if args.asimov else ''),
        [ scalingbbH ],
        x_min=yukawa_x_min, x_max=yukawa_x_max, y_min=yukawa_y_min, y_max=yukawa_y_max
        )
    plot.draw()


@flag_as_option
def multicont_Yukawa_scalingbbH(args):
    scans = []
    compare_bbH = False

    if compare_bbH:
        G0A = differentials.scans.Scan2D('nominal', x_coupling, y_coupling, scandir = scandict_G.G0A)
        G0A.title = 'Nominal'
        G0A.color = 2
        G0A.read()
        scans.append(G0A)

    scalingbbH = differentials.scans.Scan2D('scalingbbH', x_coupling, y_coupling, scandir = 'out/Scan_Yukawa_May18_combination_scalingbbH_asimov')
    scalingbbH.title = 'incl. scaling bbH'
    scalingbbH.color = 1
    scalingbbH.read()
    scans.append(scalingbbH)

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_scalingbbH' + ( '_comparebbH' if compare_bbH else '' ),
        scans,
        x_min=yukawa_x_min, x_max=yukawa_x_max, y_min=yukawa_y_min, y_max=yukawa_y_max
        )
    plot.draw()


#____________________________________________________________________
@flag_as_option
def multicont_Yukawa_G2A(args):
    G2A = differentials.scans.Scan2D('G2A', x_coupling, y_coupling, scandir = scandict_G.G2A)
    G2A.color = 2
    G2A.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_G2A',
        [ 
            G2A,
            ],
        x_min=yukawa_x_min, x_max=yukawa_x_max, y_min=yukawa_y_min, y_max=yukawa_y_max
        )
    plot.draw()


@flag_as_option
def multicont_Yukawa_G0(args):
    G0A = differentials.scans.Scan2D('G0A', x_coupling, y_coupling, scandir = scandict_G.G0A)
    G0A.color = 4
    G0A.read()

    G0B_reweighted = differentials.scans.Scan2D('G0B_reweighted', x_coupling, y_coupling, scandir = scandict_G.G0B_reweighted)
    G0B_reweighted.title = 'G0B'
    G0B_reweighted.color = 2
    G0B_reweighted.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_G0',
        [ 
            G0A,
            G0B_reweighted,
            ],
        x_min=yukawa_x_min, x_max=yukawa_x_max, y_min=yukawa_y_min, y_max=yukawa_y_max
        )
    plot.draw()

@flag_as_option
def multicont_Yukawa_G1(args):
    G1A = differentials.scans.Scan2D('G1A', x_coupling, y_coupling, scandir = scandict_G.G1A)
    G1A.color = 38
    G1A.read()

    G1B_reweighted = differentials.scans.Scan2D('G1B_reweighted', x_coupling, y_coupling, scandir = scandict_G.G1B_reweighted)
    G1B_reweighted.color = 46
    G1B_reweighted.read()

    G1A_unreweighted = differentials.scans.Scan2D('G1A_unreweighted', x_coupling, y_coupling, scandir = scandict_G.G1A_unreweighted)
    G1A_unreweighted.color = 4
    G1A_unreweighted.read()

    G1B = differentials.scans.Scan2D('G1B', x_coupling, y_coupling, scandir = scandict_G.G1B)
    G1B.color = 2
    G1B.read()

    # G1BKV = differentials.scans.Scan2D('G1BKV', x_coupling, y_coupling, scandir = scandict_G.G1BKV)
    # G1BKV.color = 9
    # G1BKV.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_G1',
        [ 
            G1A_unreweighted,
            G1B,
            G1A,
            G1B_reweighted,
            ],
        x_min=0.35*yukawa_x_min, x_max=0.35*yukawa_x_max, y_min=0.35*yukawa_y_min, y_max=0.35*yukawa_y_max
        )
    plot.only_1sigma_contours = True
    plot.draw()


scandict_G.G1A_noTheoryUnc = 'out/Scan_Yukawa_May15_combination_G1A_reweighted_noTheoryUnc_asimov'
# scandict_G.G1B_noTheoryUnc = 'out/Scan_Yukawa_May15_combination_G1B_reweighted_noTheoryUnc_asimov'
# scandict_G.G1A_noTheoryUnc_scaledByMuTotalXS = 'out/Scan_Yukawa_May16_combination_G1A_reweighted_noTheoryUnc_scaledByMuTotalXS_asimov'
scandict_G.G1B_noTheoryUnc = 'out/Scan_Yukawa_May16_combination_G1B_reweighted_noTheoryUnc_asimov'
scandict_G.G1A_noTheoryUnc_scaledByMuTotalXS = 'out/Scan_Yukawa_May16_combination_G1A_reweighted_noTheoryUnc_scaledByMuTotalXS_asimov_0'

@flag_as_option
def multicont_Yukawa_G1_noTheoryUnc(args):
    G1A = differentials.scans.Scan2D('G1A', x_coupling, y_coupling, scandir = scandict_G.G1A_noTheoryUnc)
    G1A.color = 2
    G1A.read()

    G1B = differentials.scans.Scan2D('G1B', x_coupling, y_coupling, scandir = scandict_G.G1B_noTheoryUnc)
    G1B.color = 4
    G1B.read()

    # G1A_scaledByMuTotalXS = differentials.scans.Scan2D('G1A_scaledByMuTotalXS', x_coupling, y_coupling, scandir = scandict_G.G1A_noTheoryUnc_scaledByMuTotalXS)
    # G1A_scaledByMuTotalXS.color = 3
    # G1A_scaledByMuTotalXS.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_G1_noTheoryUnc',
        [ 
            G1A,
            G1B,
            # G1A_scaledByMuTotalXS,
            ],
        x_min=0.35*yukawa_x_min, x_max=0.35*yukawa_x_max, y_min=0.35*yukawa_y_min, y_max=0.35*yukawa_y_max
        )
    plot.only_1sigma_contours = True
    plot.draw()




@flag_as_option
def multicont_Yukawa_compareBRuncertainties(args):
    args = differentialutils.set_one_decay_channel(args, 'combination', asimov=True)
    scandict = LatestPaths.scan.yukawa.reweighted.asimov

    # combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling, scandir = scandict.combination)
    # combination.color = 1
    # combination.read()
    
    # hgg = differentials.scans.Scan2D('hgg', x_coupling, y_coupling, scandir = scandict.hgg)
    # hgg.color = 2
    # hgg.read()

    hzz = differentials.scans.Scan2D('hzz', x_coupling, y_coupling, scandir = scandict.hzz)
    hzz.color = 4
    hzz.read()

    hzz_BRU = differentials.scans.Scan2D('hzz_BRU', x_coupling, y_coupling, scandir = 'out/Scan_Yukawa_May08_hzz_withBRuncertainties_asimov')
    hzz_BRU.color = 2
    hzz_BRU.title = differentials.core.standard_titles['hzz'] + ' with BR unc.'
    hzz_BRU.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_compareBRuncertainties' + ('_asimov' if args.asimov else ''),
        [hzz, hzz_BRU],
        x_min=yukawa_x_min, x_max=yukawa_x_max, y_min=yukawa_y_min, y_max=yukawa_y_max
        )
    plot.draw()


@flag_as_option
def multicont_Yukawa(args):
    multicont_Yukawa_reweighted(args)
    multicont_Yukawa_unreweighted(args)

@flag_as_option
def multicont_Yukawa_reweighted(args):
    scandict = LatestPaths.scan.yukawa.reweighted.asimov if args.asimov else LatestPaths.scan.yukawa.reweighted.observed

    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling, scandir = scandict.combination)
    combination.color = 1
    combination.read()
    
    hgg = differentials.scans.Scan2D('hgg', x_coupling, y_coupling, scandir = scandict.hgg)
    hgg.color = 2
    hgg.read()

    hzz = differentials.scans.Scan2D('hzz', x_coupling, y_coupling, scandir = scandict.hzz)
    hzz.color = 4
    hzz.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_reweighted' + ('_asimov' if args.asimov else ''),
        [combination, hgg, hzz],
        x_min=yukawa_x_min, x_max=yukawa_x_max, y_min=yukawa_y_min, y_max=yukawa_y_max
        )
    plot.draw()


@flag_as_option
def multicont_Yukawa_reweighted_paper(args):

    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling, scandir = scandict.combination)
    combination.color = 1
    combination.read()
    
    hgg = differentials.scans.Scan2D('hgg', x_coupling, y_coupling, scandir = scandict.hgg)
    hgg.color = 2
    hgg.read()

    hzz = differentials.scans.Scan2D('hzz', x_coupling, y_coupling, scandir = scandict.hzz)
    hzz.color = 4
    hzz.read()

    plot = differentials.plotting.plots.Single2DHistPlot(
        'multicont_Yukawa_paper' + ('_asimov' if args.asimov else ''),
        combination.to_hist()
        )

    # plot = differentials.plotting.plots.MultiContourPlot(
    #     'multicont_Yukawa_reweighted' + ('_asimov' if args.asimov else ''),
    #     [combination, hgg, hzz],
    #     x_min=yukawa_x_min, x_max=yukawa_x_max, y_min=yukawa_y_min, y_max=yukawa_y_max
    #     )
    plot.draw()
    plot.wrapup()



@flag_as_option
def multicont_Yukawa_unreweighted(args):
    scandict = LatestPaths.scan.yukawa.unreweighted.asimov if args.asimov else LatestPaths.scan.yukawa.unreweighted.observed

    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling, scandir = scandict.combination)
    combination.color = 1
    combination.read()
    
    hgg = differentials.scans.Scan2D('hgg', x_coupling, y_coupling, scandir = scandict.hgg)
    hgg.color = 2
    hgg.read()

    hzz = differentials.scans.Scan2D('hzz', x_coupling, y_coupling, scandir = scandict.hzz)
    hzz.color = 4
    hzz.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_unreweighted' + ('_asimov' if args.asimov else ''),
        [combination, hgg, hzz],
        x_min=yukawa_x_min, x_max=yukawa_x_max, y_min=yukawa_y_min, y_max=yukawa_y_max
        )
    plot.draw()


#____________________________________________________________________
# Extra studies

@flag_as_option
def multicont_Yukawa_highLumi(args):

    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling, scandir=LatestPaths.scan.yukawa.reweighted.asimov.combination)
    combination.color = 1
    combination.title = '35.9 fb^{-1}'
    combination.read()

    lumi300fb = differentials.scans.Scan2D('lumi300fb', x_coupling, y_coupling, scandir=LatestPaths.scan.yukawa.lumi300fb)
    lumi300fb.color = 4
    lumi300fb.title = '300 fb^{-1}'
    lumi300fb.read()
    lumi300fb.contour_filter_method = 'max_distance_to_com'

    # raise NotImplementedError('Need to multiply values from 300 fb-1')
    # lumi_3000 = differentials.scans.Scan2D('lumi_3000', x_coupling, y_coupling,
    #     scandir=LatestPaths.scan_combined_Yukawa_lumiStudy_asimov
    #     )
    # lumi_3000.color = 4
    # lumi_3000.title = '3000 fb^{-1}'
    # lumi_3000.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_highLumi',
        [combination, lumi300fb],
        x_min=yukawa_x_min, x_max=yukawa_x_max, y_min=yukawa_y_min, y_max=yukawa_y_max
        )
    plot.draw()


@flag_as_option
def multicont_Yukawa_profiledTotalXS(args):
    # Range too small... range should be higher to observe starlike pattern.
    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling, scandir=LatestPaths.scan.yukawa.reweighted.asimov.combination)
    combination.color = 1
    combination.title = 'Nominal'
    combination.read()

    profiled_total_XS = differentials.scans.Scan2D('profiled_total_XS', x_coupling, y_coupling, scandir=LatestPaths.scan.yukawa.profiledTotalXS)
    profiled_total_XS.color = 4
    profiled_total_XS.title = '#sigma_{tot} profiled'
    profiled_total_XS.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_profiledTotalXS',
        [combination, profiled_total_XS],
        x_min=yukawa_x_min, x_max=yukawa_x_max, y_min=yukawa_y_min, y_max=yukawa_y_max
        )
    plot.draw()


#____________________________________________________________________
# To redo

@flag_as_option
def multicont_Yukawa_ratio_of_BRs(args):

    containers = []

    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling,
        scandir=LatestPaths.scan_combined_Yukawa_asimov
        )
    combination.color = 1
    combination.title = 'SM BR'
    combination.read()

    ratio_of_BRs = differentials.scans.Scan2D('ratio_of_BRs', x_coupling, y_coupling,
        scandir=LatestPaths.scan_combined_Yukawa_ratioOfBRs_asimov
        )
    ratio_of_BRs.color = 4
    ratio_of_BRs.title = 'floating BR_{H#gamma#gamma}/BR_{HZZ}'

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_ratio_of_BRs',
        [combination, ratio_of_BRs],
        x_min=yukawa_x_min, x_max=yukawa_x_max, y_min=yukawa_y_min, y_max=yukawa_y_max
        )
    plot.draw()


@flag_as_option
def multicont_Yukawa_onlyNormalization(args):

    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling,
        scandir=LatestPaths.scan_combined_Yukawa_asimov
        )
    combination.color = 1
    combination.title = 'Nominal'
    combination.read()

    only_normalization = differentials.scans.Scan2D('only_normalization', x_coupling, y_coupling,
        scandir=LatestPaths.scan_combined_Yukawa_onlyNormalization_asimov
        )
    only_normalization.color = 4
    only_normalization.title = 'Only normalization'
    only_normalization.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_onlyNormalization',
        [combination, only_normalization],
        x_min=yukawa_x_min, x_max=yukawa_x_max, y_min=yukawa_y_min, y_max=yukawa_y_max
        )
    plot.draw()


@flag_as_option
def multicont_Yukawa_theoryCrossCheck(args):

    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling,
        scandir=LatestPaths.scan_combined_Yukawa_asimov
        )
    combination.color = 1
    combination.title = 'Nominal'
    combination.read()

    no_theory_unc = differentials.scans.Scan2D('no_theory_unc', x_coupling, y_coupling,
        scandir=LatestPaths.scan_combined_Yukawa_noTheoryUncertainties_asimov 
        )
    no_theory_unc.color = 4
    no_theory_unc.title = 'No th. uncertainty'
    no_theory_unc.read()

    uncorr_theory_unc = differentials.scans.Scan2D('uncorr_theory_unc', x_coupling, y_coupling,
        scandir=LatestPaths.scan_combined_Yukawa_uncorrelatedTheoryUncertainties_asimov 
        )
    uncorr_theory_unc.color = 2
    uncorr_theory_unc.title = 'Uncorr. uncertainty'
    uncorr_theory_unc.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_theoryCrossCheck',
        [combination, no_theory_unc, uncorr_theory_unc],
        x_min=yukawa_x_min, x_max=yukawa_x_max, y_min=yukawa_y_min, y_max=yukawa_y_max
        )
    plot.draw()


@flag_as_option
def multicont_Yukawa_BRdependencyComparison(args):

    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling,
        scandir=LatestPaths.scan_combined_Yukawa_asimov
        )
    combination.color = 1
    combination.title = 'Nominal'
    combination.read()

    brscal_kappaV_fixed = differentials.scans.Scan2D('brscal_kappaV_fixed', x_coupling, y_coupling,
        scandir=LatestPaths.scan_combined_Yukawa_couplingDependentBR_fixedKappaV_asimov
        )
    brscal_kappaV_fixed.color = 4
    brscal_kappaV_fixed.title = 'BR(#kappa_{b}, #kappa_{c})' # + ' (#kappa_{V} fixed)'
    brscal_kappaV_fixed.read()

    brscal_profiledTotalXS = differentials.scans.Scan2D('brscal_profiledTotalXS', x_coupling, y_coupling,
        scandir='out/Scan_Yukawa_Feb07_combination_couplingDependentBR_profiledTotalXS'
        )
    brscal_profiledTotalXS.color = 2
    brscal_profiledTotalXS.title = 'BR(#kappa_{b}, #kappa_{c}, #sigma_{tot})' # + ' (#kappa_{V} fixed)'
    brscal_profiledTotalXS.read()

    brscal_kappaV_max1 = differentials.scans.Scan2D('brscal_kappaV_max1', x_coupling, y_coupling,
        scandir=LatestPaths.scan_combined_Yukawa_couplingDependentBR_kappaVMaxOne_asimov
        )
    brscal_kappaV_max1.color = 8
    brscal_kappaV_max1.title = 'BR(#kappa_{t}, #kappa_{V}#leq1)'
    brscal_kappaV_max1.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_BRdependencyComparison',
        [combination, brscal_kappaV_fixed, brscal_profiledTotalXS, brscal_kappaV_max1],
        x_min=yukawa_x_min, x_max=yukawa_x_max, y_min=yukawa_y_min, y_max=yukawa_y_max
        )
    plot.draw()


@flag_as_option
def multicont_Yukawa_atfchi2(args):

    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling,
        scandir=LatestPaths.scan_combined_Yukawa_noTheoryUncertainties_asimov
        )
    combination.color = 1
    combination.title = 'Full likelihood (no theor. unc.)'
    combination.read()

    raise NotImplementedError('Have to fix the \'as theorist\' scan')

    # Get relevant info from a canvas that has the histogram and the bestfit point
    atfchi2_rootFp = ROOT.TFile.Open( 'plots_Nov14_Yukawa/AfterTheFactChi2Fit.root' )
    canvas = atfchi2_rootFp.Get('ctc')

    # Get TH2F
    atfchi2_H = canvas.GetPrimitive('H2')
    ROOT.SetOwnership( atfchi2_H, False )

    # Get values of x and y for bestfit
    bestfitpoint_Tg = canvas.GetPrimitive('bestfitpoint')
    x_Double = ROOT.Double(0)
    y_Double = ROOT.Double(0)
    bestfitpoint_Tg.GetPoint( 0, x_Double, y_Double )
    xBestfit = float(x_Double)
    yBestfit = float(y_Double)

    atfchi2_rootFp.Close()

    # Construct container
    atfchi2 = Container()
    atfchi2.H2       = atfchi2_H
    atfchi2.xBestfit = xBestfit
    atfchi2.yBestfit = yBestfit
    atfchi2.color    = 2
    atfchi2.name     = 'atfchi2'
    atfchi2.title    = 'Simple #chi^{2} fit'
    containers.append( atfchi2 )


    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_atfchi2',
        [combination, atfchi2],
        x_min=yukawa_x_min, x_max=yukawa_x_max, y_min=yukawa_y_min, y_max=yukawa_y_max
        )
    plot.draw()


# ======================================
# Other plots

@flag_as_option
def points_on_contour_Yukawa(args):
    if args.asimov:
        scandir = LatestPaths.scan.yukawa.reweighted.asimov.combination
    else:
        # scandir = LatestPaths.scan.yukawa.reweighted.observed.combination
        scandir = approval.scalingbbH_couplingdependentBRs.data

    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling, scandir=scandir)
    combination.color = 1
    combination.read()

    obs = LatestBinning.obs_pth_ggH
    obs.drop_bins_up_to_value(125.)

    obs_name = 'pth_ggH'

    scandict = LatestPaths.scan.pth_ggH.asimov if args.asimov else LatestPaths.scan.pth_ggH.observed
    pth_ggH_combination = differentials.scans.DifferentialSpectrum('combination', scandict.combination)
    pth_ggH_combination.color = 1
    pth_ggH_combination.no_overflow_label = True
    pth_ggH_combination.draw_method = 'repr_point_with_vertical_bar'
    pth_ggH_combination.set_sm(obs.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=False))
    pth_ggH_combination.drop_bins_up_to_value(125.)
    pth_ggH_combination.read()

    # ws = LatestPaths.ws_combined_TopHighPt
    # ws = LatestPaths.ws_combined_Yukawa
    # ws = 'out/workspaces_Feb01/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties.root'
    # ws = LatestPaths.ws.yukawa.nominal.combination
    ws = 'out/workspaces_May30/combination_Yukawa_reweighted_couplingdependentBRs.root'

    # ======================================
    # Load into plot

    plot = differentials.plotting.plots.BottomPanelPlotWithParametrizations('points_on_contour_Yukawa')
    plot.scan2D = combination
    plot.ws_file = ws
    plot.ptspectrum = pth_ggH_combination
    plot.obs = obs

    # plot.top_y_min = top_y_min
    plot.top_y_max = 100.

    obs_title = 'p_{T}'
    obs_unit  = 'GeV'
    plot.x_title = obs_title + (' ({0})'.format(obs_unit) if not(obs_unit is None) else '')
    plot.y_title_top = '#Delta#sigma/#Delta{0} (pb{1})'.format(
        obs_title,
        '/' + obs_unit if not(obs_unit is None) else ''
        )
    plot.y_title_bottom = '#mu'

    plot.draw()



@flag_as_option
def one_kappa_scan(args):
    one_kappa_scan_plot_manualProfile(args, 'kappab')
    one_kappa_scan_plot(args, 'kappac')

def one_kappa_scan_plot(args, kappa):
    differentials.plotting.canvas.c.resize_temporarily(850, 800)
    scandict = LatestPaths.scan.yukawa.onekappa
    scans = []

    observed = differentials.scans.Scan(kappa, scandir=scandict.observed[kappa])
    observed.color = 1
    observed.title = '{0} observed'.format(differentials.core.standard_titles[kappa])
    observed.read()
    observed.create_uncertainties()
    observed.title += '; ({0:.2f} - {1:.2f}) @ 68% CL'.format(observed.unc.left_bound, observed.unc.right_bound)
    scans.append(observed)

    expected = differentials.scans.Scan(kappa, scandir=scandict.asimov[kappa])
    expected.color = 1
    expected.title = '{0} expected'.format(differentials.core.standard_titles[kappa])
    expected.read()
    expected.create_uncertainties()
    expected.title += '; ({0:.2f} - {1:.2f}) @ 68% CL'.format(expected.unc.left_bound, expected.unc.right_bound)
    expected.draw_style = 'repr_dashed_line'
    scans.append(expected)

    plot = differentials.plotting.plots.MultiScanPlot('onekappascan_{0}'.format(kappa))
    plot.scans.extend(scans)

    plot.x_title = differentials.core.standard_titles[kappa]
    if kappa == 'kappab':
        plot.x_min = -7.
        plot.x_max = 9.
    else:
        plot.x_min = -20.
        plot.x_max = 20.

    plot.leg.SetNColumns(1)
    plot.leg._y1 = lambda c: 1. - c.GetTopMargin() - 0.20

    plot.draw()
    plot.wrapup()


def one_kappa_scan_plot_manualProfile(args, kappa):
    differentials.plotting.canvas.c.resize_temporarily(850, 800)
    scandict = LatestPaths.scan.yukawa.onekappa

    expected = differentials.scans.Scan(kappa, scandir=scandict.asimov[kappa])
    expected.color = 1
    expected.title = '{0} expected'.format(differentials.core.standard_titles[kappa])
    expected.read()
    expected.create_uncertainties()
    expected.title += '; ({0:.2f} - {1:.2f}) @ 68% CL'.format(expected.unc.left_bound, expected.unc.right_bound)
    expected.draw_style = 'repr_dashed_line'

    # Get the 2D scan
    observed2D = differentials.scans.Scan2D('observed2D', x_coupling, y_coupling,
        scandir = LatestPaths.scan.yukawa.reweighted.observed.combination
        )
    observed2D.color = 1
    observed2D.read()

    # Read the observed from the 2D scan; take a minimum per slice
    observed1D = observed2D.get_1d(kappa)
    observed1D.draw_style = 'repr_smooth_line'
    observed1D.title = '{0} observed; ({1:.2f} - {2:.2f}) @ 68% CL'.format(
        differentials.core.standard_titles[kappa],
        observed1D.unc.left_bound, observed1D.unc.right_bound
        )

    plot = differentials.plotting.plots.MultiScanPlot('onekappascan_{0}'.format(kappa))
    plot.scans.append(expected)
    plot.manual_graphs.append(observed1D)
    plot.x_title = differentials.core.standard_titles[kappa]
    if kappa == 'kappab':
        plot.x_min = -7.
        plot.x_max = 9.
    else:
        plot.x_min = -20.
        plot.x_max = 20.
    plot.leg.SetNColumns(1)
    plot.leg._y1 = lambda c: 1. - c.GetTopMargin() - 0.20
    plot.draw()
    plot.wrapup()


