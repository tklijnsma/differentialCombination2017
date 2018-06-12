#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

from OptionHandler import flag_as_option

import copy, math, sys, ROOT
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

ktcg_couplingdependentBRs_x_min = -0.1
ktcg_couplingdependentBRs_x_max = 2.0
ktcg_couplingdependentBRs_y_min = -0.05
ktcg_couplingdependentBRs_y_max = 0.085

ktcg_floatingBRs_x_min = -2.4
ktcg_floatingBRs_x_max = 2.4
ktcg_floatingBRs_y_min = -0.15
ktcg_floatingBRs_y_max = 0.15


#____________________________________________________________________
# Results

approval = differentials.core.AttrDict.create_tree(['ktcg', 'ktkb'], ['observed', 'asimov'], ['fixedBRs', 'couplingdependentBRs', 'floatingBRs'])

approval.ktcg.asimov.fixedBRs.combWithHbb = 'out/Scan_May18_Top_combWithHbb_scalingttH_asimov'
approval.ktcg.observed.fixedBRs.combWithHbb = 'out/Scan_Jun03_Top_combWithHbb_scalingttH'

approval.ktcg.asimov.couplingdependentBRs.combWithHbb = 'out/Scan_May22_Top_combWithHbb_scalingttH_couplingdependentBRs_asimov'
approval.ktcg.observed.couplingdependentBRs.combWithHbb = 'out/Scan_May31_Top_combWithHbb_scalingttH_couplingdependentBRs'
approval.ktcg.observed.couplingdependentBRs.hgg = 'out/Scan_Jun09_Top_hgg_scalingttH_couplingdependentBRs'
approval.ktcg.observed.couplingdependentBRs.hzz = 'out/Scan_Jun09_Top_hzz_scalingttH_couplingdependentBRs'

approval.ktcg.asimov.floatingBRs.combWithHbb    = 'out/Scan_Jun11_Top_combWithHbb_scalingttH_floatingBRs_constrainedbbZZ_asimov'
approval.ktcg.asimov.floatingBRs.hgg            = 'out/Scan_Jun11_Top_hgg_scalingttH_floatingBRs_constrainedbbZZ_asimov'
approval.ktcg.observed.floatingBRs.combWithHbb  = 'out/Scan_May31_Top_combWithHbb_scalingttH_floatingBRs_constrainedbbZZ_0'
approval.ktcg.observed.floatingBRs.hgg          = 'out/Scan_Jun09_Top_hgg_scalingttH_floatingBRs_constrainedbbZZ'
approval.ktcg.observed.floatingBRs.hzz          = 'out/Scan_Jun09_Top_hzz_scalingttH_floatingBRs_constrainedbbZZ'

approval.ktkb.observed.floatingBRs.combWithHbb  = 'out/Scan_May31_TopCtCb_combWithHbb_scalingbbHttH_floatingBRs_constrainedbbZZ'
# approval.ktkb.observed.floatingBRs.combWithHbb = 'out/Scan_Jun10_TopCtCb_combWithHbb_scalingbbHttH_floatingBRs_constrainedbbZZ'
approval.ktkb.observed.floatingBRs.hgg          = 'out/Scan_Jun09_TopCtCb_hgg_scalingbbHttH_floatingBRs_constrainedbbZZ'
approval.ktkb.observed.floatingBRs.hzz          = 'out/Scan_Jun09_TopCtCb_hzz_scalingbbHttH_floatingBRs_constrainedbbZZ'

approval.ktkb.observed.couplingdependentBRs.combWithHbb = 'out/Scan_May29_TopCtCb_combWithHbb_scalingbbHttH_couplingdependentBRs_0'
# approval.ktkb.observed.couplingdependentBRs.hgg = 'out/Scan_Jun09_TopCtCb_hgg_scalingbbHttH_couplingdependentBRs'
approval.ktkb.observed.couplingdependentBRs.hgg = 'out/Scan_Jun10_TopCtCb_hgg_scalingbbHttH_couplingdependentBRs'
approval.ktkb.observed.couplingdependentBRs.hzz = 'out/Scan_Jun09_TopCtCb_hzz_scalingbbHttH_couplingdependentBRs'


def latest_ktkg_couplingdependentBRs(args, decay_channel=None, asimov=None, splined=False):
    if not(asimov is None): args = differentialutils.force_asimov(args, asimov)
    if not(decay_channel is None): args = differentialutils.set_one_decay_channel(args, decay_channel)
    decay_channel = differentialutils.get_decay_channel_tag(args)
    scandir = approval.ktcg['asimov' if args.asimov else 'observed'].couplingdependentBRs[decay_channel]

    scan = differentials.scans.Scan2D(
        'ktcg_{0}'.format(decay_channel), x_coupling, y_coupling,
        scandir = scandir
        )
    scan.title = differentials.core.standard_titles.get(decay_channel, decay_channel)
    scan.color = 1
    scan.read()

    if splined:
        return spline_ktkg_couplingdependentBRs(args, scan)
    return scan

def spline_ktkg_couplingdependentBRs(args, scan):
    x_min = ktcg_couplingdependentBRs_x_min
    x_max = ktcg_couplingdependentBRs_x_max
    y_min = ktcg_couplingdependentBRs_y_min
    y_max = ktcg_couplingdependentBRs_y_max
    deltaNLL_cutoff = 30.

    if args.combWithHbb:
        if args.asimov:
            deltaNLL_cutoff = 10.
            x_max = 2.5
            y_min = -0.11

    spline = scan.to_spline(
        x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max,
        deltaNLL_cutoff = deltaNLL_cutoff
        )

    if args.combWithHbb:
        if args.asimov:
            spline.add_noise_selector(lambda ct, cg: (cg  <  (1./12.)-0.02 - (1./12.)*ct))
            spline.add_noise_selector(lambda ct, cg: (cg  >  (1./12.)+0.03 - (1./13.)*ct))
        else:
            spline.add_noise_selector(lambda ct, cg: (cg  <  (1./12.)-0.02 - (1./12.)*ct))
            spline.add_noise_selector(lambda ct, cg: (cg  >  (1./12.)+0.04 - (1./13.)*ct))
    if args.hgg:
        spline.negativity_is_zero = True
    
    hist = spline.to_hist(nx=180, ny=180)
    hist.color = scan.color
    hist.name  = scan.name + '_splined'
    hist.title = scan.title
    return hist

@flag_as_option
def multicont_Top_scalingttH_couplingdependentBRs(args):
    scans = []
    scans.append(latest_ktkg_couplingdependentBRs(args, decay_channel='combWithHbb', splined=True))
    if not(args.asimov):
        hgg = latest_ktkg_couplingdependentBRs(args, decay_channel='hgg', splined=True)
        hgg.color = differentials.core.safe_colors.red
        scans.append(hgg)
        hzz = latest_ktkg_couplingdependentBRs(args, decay_channel='hzz')
        hzz.color = differentials.core.safe_colors.blue
        scans.append(hzz)
    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_ktcg_couplingdependentBRs' + ('_asimov' if args.asimov else ''),
        scans,
        x_title = differentials.core.standard_titles['ct'], y_title = differentials.core.standard_titles['cg'],
        x_min = ktcg_couplingdependentBRs_x_min, x_max = ktcg_couplingdependentBRs_x_max, y_min = ktcg_couplingdependentBRs_y_min, y_max = ktcg_couplingdependentBRs_y_max
        )
    plot.y_SM = 0.0
    x_shift = 0.48
    plot.legend.set(
        lambda c: c.GetLeftMargin() + 0.02 + x_shift,
        lambda c: 1. - c.GetTopMargin() - 0.20,
        lambda c: c.GetLeftMargin() + 0.24 + x_shift,
        lambda c: 1. - c.GetTopMargin() - 0.01,
        )
    plot.draw()

@flag_as_option
def quicktest_splining_ktcg_couplingdependentBRs(args):
    unsplined = latest_ktkg_couplingdependentBRs(args, decay_channel='combWithHbb')
    splined = latest_ktkg_couplingdependentBRs(args, decay_channel='combWithHbb', splined=True)
    splined.quickplot('quicktest_splining_ktcg_couplingdependentBRs_splined' + ('_asimov' if args.asimov else ''))
    unsplined.to_hist().quickplot('quicktest_splining_ktcg_couplingdependentBRs_unsplined' + ('_asimov' if args.asimov else ''))

@flag_as_option
def onedimscans_ktcg_couplingdependentBRs(args):
    obs2D = latest_ktkg_couplingdependentBRs(args, decay_channel='combWithHbb', splined=True)
    exp2D = latest_ktkg_couplingdependentBRs(args, decay_channel='combWithHbb', splined=True, asimov=True)
    onedimscans_ktcg('ct', obs2D, exp2D, x_min=0.0, x_max=2.0, tag='ktcg_couplingdependentBRs')
    onedimscans_ktcg('cg', obs2D, exp2D, x_min=-0.05, x_max=0.07, tag='ktcg_couplingdependentBRs')

def onedimscans_ktcg(kappa, obs2D, exp2D, x_min, x_max, apply_smoothing=False, tag=None):
    onedimscanner_obs = differentials.onedimscanner.OneDimScanner(
        obs2D,
        'ct', 'cg'
        )
    obs1D = onedimscanner_obs.get_1d(kappa)
    obs1D.draw_style = 'repr_smooth_line'
    obs1D.title = '{0} observed; ({1:.2f} - {2:.2f}) @ 68% CL'.format(
        differentials.core.standard_titles[kappa],
        obs1D.unc.left_bound, obs1D.unc.right_bound
        )
    obs1D.color = 1
    if apply_smoothing: obs1D.smooth_y(10)

    onedimscanner_exp = differentials.onedimscanner.OneDimScanner(
        exp2D,
        'ct', 'cg'
        )
    exp1D = onedimscanner_exp.get_1d(kappa)
    exp1D.draw_style = 'repr_dashed_line'
    exp1D.title = '{0} expected; ({1:.2f} - {2:.2f}) @ 68% CL'.format(
        differentials.core.standard_titles[kappa],
        exp1D.unc.left_bound, exp1D.unc.right_bound
        )
    exp1D.color = 1
    if apply_smoothing: exp1D.smooth_y(10)

    tag = '' if tag is None else '_'+tag
    plotname = 'onekappascan{0}_{1}'.format(tag, kappa)

    differentials.plotting.canvas.c.resize_temporarily(850, 800)
    plot = differentials.plotting.plots.MultiScanPlot(plotname)
    plot.manual_graphs.append(obs1D)
    plot.manual_graphs.append(exp1D)
    plot.x_title = differentials.core.standard_titles[kappa]
    plot.x_min = x_min
    plot.x_max = x_max
    plot.leg.SetNColumns(1)
    plot.leg._y1 = lambda c: 1. - c.GetTopMargin() - 0.20
    plot.draw()
    plot.wrapup()


#____________________________________________________________________
def latest_ktkg_floatingBRs(args, decay_channel=None, asimov=None, splined=False):
    if not(asimov is None): args = differentialutils.force_asimov(args, asimov)
    if not(decay_channel is None): args = differentialutils.set_one_decay_channel(args, decay_channel)
    decay_channel = differentialutils.get_decay_channel_tag(args)
    scandir = approval.ktcg['asimov' if args.asimov else 'observed'].floatingBRs[decay_channel]

    scan = differentials.scans.Scan2D(
        'ktcg_{0}'.format(decay_channel), x_coupling, y_coupling,
        scandir = scandir
        )
    scan.title = differentials.core.standard_titles.get(decay_channel, decay_channel)
    scan.color = 1
    scan.read()

    if splined:
        return spline_ktkg_floatingBRs(args, scan)
    return scan

def spline_ktkg_floatingBRs(args, scan):
    x_min = ktcg_floatingBRs_x_min
    x_max = ktcg_floatingBRs_x_max
    y_min = ktcg_floatingBRs_y_min
    y_max = ktcg_floatingBRs_y_max
    deltaNLL_cutoff = 15.
    eps = 2.2
    if args.hgg:
        x_min = -2.9
        x_max = 2.9
        y_min = -0.15
        y_max = 0.12
        deltaNLL_cutoff = 50.
        eps = 1.2

    spline = scan.to_spline(
        x_min = x_min,
        x_max = x_max,
        y_min = y_min,
        y_max = y_max,
        deltaNLL_cutoff = deltaNLL_cutoff,
        eps = eps
        )
    spline.disallow_negativity = False
    # spline.add_noise_selector(
    #     lambda ct, cg: (cg  <  (1./12.)-0.02 - (1./12.)*ct)
    #     )
    # spline.add_noise_selector(
    #     lambda ct, cg: (cg  >  (1./12.)+0.04 - (1./13.)*ct)
    #     )
    hist = spline.to_hist(nx=180, ny=180)
    hist.color = 1
    hist.name  = scan.name + '_splined'
    hist.title = scan.title
    return hist

@flag_as_option
def multicont_Top_scalingttH_floatingBRs_constrainedbbZZ(args):
    scans = []

    if args.asimov:
        scans.append(latest_ktkg_floatingBRs(args, decay_channel='combWithHbb', splined=False))
        hgg = latest_ktkg_floatingBRs(args, decay_channel='hgg', splined=False)
        hgg.color = differentials.core.safe_colors.red
        scans.append(hgg)
    else:
        scans.append(latest_ktkg_floatingBRs(args, decay_channel='combWithHbb', splined=True))
        hgg = latest_ktkg_floatingBRs(args, decay_channel='hgg', splined=True)
        hgg.color = differentials.core.safe_colors.red
        scans.append(hgg)
        hzz = latest_ktkg_floatingBRs(args, decay_channel='hzz')
        hzz.color = differentials.core.safe_colors.blue
        scans.append(hzz)
    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_ktcg_floatingBRs' + ('_asimov' if args.asimov else ''),
        scans,
        x_title = differentials.core.standard_titles['ct'], y_title = differentials.core.standard_titles['cg'],
        x_min = ktcg_floatingBRs_x_min, x_max = ktcg_floatingBRs_x_max, y_min = ktcg_floatingBRs_y_min, y_max = ktcg_floatingBRs_y_max
        )
    plot.y_SM = 0.0
    plot.draw()

@flag_as_option
def quicktest_splining_ktcg_floatingBRs(args):
    unsplined = latest_ktkg_floatingBRs(args, decay_channel='combWithHbb')
    splined = latest_ktkg_floatingBRs(args, decay_channel='combWithHbb', splined=True)
    splined.quickplot('quicktest_splining_ktcg_floatingBRs_splined' + ('_asimov' if args.asimov else ''))
    unsplined.to_hist().quickplot('quicktest_splining_ktcg_floatingBRs_unsplined' + ('_asimov' if args.asimov else ''))



#____________________________________________________________________
@flag_as_option
def multicont_Top_scalingttH(args):
    args = differentialutils.set_one_decay_channel(args, 'combWithHbb', asimov=True)
    # compare = True
    compare = False
    scans = []
    
    if compare:
        combWithHbb = differentials.scans.Scan2D('combWithHbb', x_coupling, y_coupling, scandir=LatestPaths.scan.top.reweighted.asimov.combWithHbb)
        combWithHbb.color = 2
        combWithHbb.read()
        scans.append(combWithHbb)

    combWithHbb_scalingttH = differentials.scans.Scan2D('combWithHbb_scalingttH', x_coupling, y_coupling, scandir='out/Scan_May18_Top_combWithHbb_scalingttH_asimov')
    combWithHbb_scalingttH.title = 'incl. scaling ttH'
    combWithHbb_scalingttH.color = 1
    combWithHbb_scalingttH.read()
    scans.append(combWithHbb_scalingttH)
    
    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Top_scalingttH' + ('_asimov' if args.asimov else '') + ('_asimov' if args.asimov else ''),
        scans,
        x_min=top_x_min, x_max=top_x_max, y_min=top_y_min, y_max=top_y_max,
        )
    # plot.legend.SetNColumns(2)
    # plot.only_1sigma_contours = True
    plot.draw()


@flag_as_option
def singlehist_Top_scalingttH(args):
    args = differentialutils.set_one_decay_channel(args, 'combWithHbb')
    if args.asimov:
        scandir = 'out/Scan_May18_Top_combWithHbb_scalingttH_asimov'
    else:
        scandir = 'out/Scan_Jun03_Top_combWithHbb_scalingttH'

    scalingttH = differentials.scans.Scan2D(
        'scalingttH', x_coupling, y_coupling,
        scandir = scandir
        )
    scalingttH.title = 'Combination (incl. bbH / BR(#vec{#kappa}))'
    scalingttH.color = 1
    scalingttH.read()

    x_min = -0.2
    x_max = 3.8
    y_min = -0.23
    y_max = 0.12
    spline = scalingttH.to_spline(
        x_min = x_min,
        x_max = x_max,
        y_min = y_min,
        y_max = y_max,
        )

    if args.asimov:
        spline.add_noise_selector(lambda ct, cg: (cg  <  (1./12.)-0.03 - (1./12.)*ct) )
        spline.add_noise_selector(lambda ct, cg: cg > 0.1 )
        spline.add_noise_selector(lambda ct, cg: (cg  >  (1./12.)+0.04 - (1./13.)*ct) )
    else:
        spline.add_noise_selector(lambda ct, cg: (cg  <  (1./12.)-0.02 - (1./12.)*ct) )
        spline.add_noise_selector(lambda ct, cg: cg > 0.095 )
        spline.add_noise_selector(lambda ct, cg: (cg  >  (1./12.)+0.04 - (1./13.)*ct) )
    hist = spline.to_hist(nx=180, ny=180)
    hist.color = 1

    plot = differentials.plotting.plots.Single2DHistPlot(
        'singlehist_Top_scalingttH' + ('_asimov' if args.asimov else ''),
        hist,
        x_min = x_min,
        x_max = x_max,
        y_min = y_min,
        y_max = y_max,
        )
    plot.set_ranges_by_contour = False
    plot.x_SM = 1.0
    plot.y_SM = 0.0
    plot.x_title = '#kappa_{t}'
    plot.y_title = 'c_{g}'
    plot.draw()
    plot.wrapup()

#____________________________________________________________________
def latest_ktkb_couplingdependentBRs(args, decay_channel=None, asimov=None, splined=False):
    if not(asimov is None): args = differentialutils.force_asimov(args, asimov)
    if not(decay_channel is None): args = differentialutils.set_one_decay_channel(args, decay_channel)
    decay_channel = differentialutils.get_decay_channel_tag(args)
    scandir = approval.ktkb['asimov' if args.asimov else 'observed'].couplingdependentBRs[decay_channel]

    scan = differentials.scans.Scan2D('ktkb_{0}'.format(decay_channel), x_coupling, 'cb', scandir = scandir)
    scan.title = differentials.core.standard_titles.get(decay_channel, decay_channel)
    scan.color = 1
    scan.read()

    if splined:
        return spline_ktkb_couplingdependentBRs(args, scan)
    return scan

def spline_ktkb_couplingdependentBRs(args, scan):
    if args.combWithHbb:
        x_min = 0.2
        x_max = 2.0
        y_min = -2.0
        y_max = 2.0
    elif args.hzz:
        x_min = -3.0
        x_max = 3.0
        y_min = -3.0
        y_max = 3.0
    elif args.hgg:
        x_min = 0.2
        x_max = 3.2
        y_min = -2.0
        y_max = 2.0
    spline = scan.to_spline(
        x_min = x_min,
        x_max = x_max,
        y_min = y_min,
        y_max = y_max,
        )
    if args.combWithHbb:
        spline.add_noise_selector(lambda kt, kb: (kt > 1.4) and (kt < 2.0) and (kb > -0.5) and (kb < 0.75))
        spline.add_noise_selector(lambda kt, kb: (kt < 0.6) and (kb < -1.0))
        spline.add_noise_selector(lambda kt, kb: (kt < 0.75) and (kb > 1.05))
        spline.add_noise_selector(lambda kt, kb: (kb > 1.6))
    if args.hzz:
        spline.disallow_negativity = False
    if args.hgg:
        x1 = 0.62
        x2 = 1.05
        y1 = -0.5
        y2 = -0.86
        a = (y2-y1)/(x2-x1)
        b = y2 - a*x2
        spline.add_signal_selector(lambda kt, kb: (kb < -a*kt-b + 0.1) and (kb > -a*kt-b - 0.1) and (kt>0.95) and (kt<1.05) )
        spline.add_signal_selector(lambda kt, kb: (kb < a*kt+b + 0.1) and (kb > a*kt+b - 0.1) and (kt>x1) and (kt<x2) )

    hist = spline.to_hist(nx=180, ny=180, x_min=-3.0, x_max=3.0, y_min=-3.0, y_max=3.0)
    hist.color = 1
    hist.name  = scan.name + '_splined'
    hist.title = scan.title
    return hist

@flag_as_option
def multicont_ktkb_scalingbbHttH_couplingdependentBRs(args):
    y_coupling = 'cb'
    scans = []
    scans.append(latest_ktkb_couplingdependentBRs(args, decay_channel='combWithHbb', splined=True))
    if not(args.asimov):
        hgg = latest_ktkb_couplingdependentBRs(args, decay_channel='hgg', splined=True)
        hgg.color = differentials.core.safe_colors.red
        scans.append(hgg)
        hzz = latest_ktkb_couplingdependentBRs(args, decay_channel='hzz', splined=True)
        hzz.color = differentials.core.safe_colors.blue
        scans.append(hzz)
    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_ktkb_couplingdependentBRs' + ('_asimov' if args.asimov else ''),
        scans,
        x_title = differentials.core.standard_titles['ct'], y_title = differentials.core.standard_titles['cb'],
        )
    plot.draw()


#____________________________________________________________________

def latest_ktkb_floatingBRs(args, decay_channel=None, asimov=None, splined=False):
    if not(asimov is None): args = differentialutils.force_asimov(args, asimov)
    if not(decay_channel is None): args = differentialutils.set_one_decay_channel(args, decay_channel)
    decay_channel = differentialutils.get_decay_channel_tag(args)
    scandir = approval.ktkb['asimov' if args.asimov else 'observed'].floatingBRs[decay_channel]

    scan = differentials.scans.Scan2D('ktkb_{0}'.format(decay_channel), x_coupling, 'cb', scandir = scandir)
    scan.title = differentials.core.standard_titles.get(decay_channel, decay_channel)
    scan.color = 1
    scan.read()

    if splined:
        return spline_ktkb_floatingBRs(args, scan)
    return scan

def spline_ktkb_floatingBRs(args, scan):
    x_min = -3.0
    x_max = 3.0
    y_min = -15.0
    y_max = 15.0

    if args.combWithHbb:
        deltaNLL_cutoff = 7.
        eps = 1.9

    if args.hgg:
        deltaNLL_cutoff = 50.
        eps = 1.0

    if args.hzz:
        deltaNLL_cutoff = 15.
        eps = 1.9

    spline = scan.to_spline(
        x_min = x_min,
        x_max = x_max,
        y_min = y_min,
        y_max = y_max,
        deltaNLL_cutoff = deltaNLL_cutoff,
        eps = eps,
        )

    if args.hzz:
        spline.add_signal_selector(lambda kt, kb: kt > 0.8 and kb > -2.)
        spline.add_signal_selector(lambda kt, kb: kt < -0.8 and kb < 2.)

    if args.combWithHbb:
        spline.disallow_negativity = False
    else:
        spline.negativity_is_zero = True

    hist = spline.to_hist(nx=180, ny=180)
    hist.color = 1
    hist.name  = scan.name + '_splined'
    hist.title = scan.title
    return hist

@flag_as_option
def multicont_ktkb_scalingbbHttH_floatingBRs_constrainedbbZZ(args):
    y_coupling = 'cb'
    scans = []

    scans = []
    scans.append(latest_ktkb_floatingBRs(args, decay_channel='combWithHbb', splined=True))
    if not(args.asimov):
        hgg = latest_ktkb_floatingBRs(args, decay_channel='hgg', splined=True)
        hgg.color = differentials.core.safe_colors.red
        scans.append(hgg)
        hzz = latest_ktkb_floatingBRs(args, decay_channel='hzz', splined=True)
        hzz.color = differentials.core.safe_colors.blue
        scans.append(hzz)
    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_ktkb_floatingBRs' + ('_asimov' if args.asimov else ''),
        scans, x_title = differentials.core.standard_titles['ct'], y_title = differentials.core.standard_titles['cb'],
        )
    plot.draw()


@flag_as_option
def quicktest_splining_ktkb_floatingBRs(args):
    unsplined = latest_ktkb_floatingBRs(args, decay_channel='combWithHbb')
    splined = latest_ktkb_floatingBRs(args, decay_channel='combWithHbb', splined=True)
    splined.quickplot('quicktest_splining_ktkb_floatingBRs_splined' + ('_asimov' if args.asimov else ''))
    unsplined.to_hist().quickplot('quicktest_splining_ktcb_floatingBRs_unsplined' + ('_asimov' if args.asimov else ''))


# @flag_as_option
# def multicont_TopCtCb_scalingbbHttH_floatingBRs_constrainedbbZZ(args):
#     args = differentialutils.set_one_decay_channel(args, 'combWithHbb')
#     y_coupling = 'cb'

#     # scalingbbHttH = differentials.scans.Scan2D(
#     #     'scalingbbHttH', x_coupling, y_coupling,
#     #     scandir = 'out/Scan_May29_TopCtCb_combWithHbb_scalingbbHttH_floatingBRs_0'
#     #     )
#     # scalingbbHttH.title = 'Combination (incl. bbH / BRs free)'
#     # scalingbbHttH.color = 1
#     # scalingbbHttH.read()

#     scalingbbHttH_bbZZconstraint = differentials.scans.Scan2D(
#         'scalingbbHttH_bbZZconstraint', x_coupling, y_coupling,
#         # scandir = 'out/Scan_May29_TopCtCb_combWithHbb_scalingbbHttH_floatingBRs_constrainedbbZZ'
#         scandir = 'out/Scan_May31_TopCtCb_combWithHbb_scalingbbHttH_floatingBRs_constrainedbbZZ'
#         )
#     scalingbbHttH_bbZZconstraint.title = 'bb/ZZ constrained'
#     scalingbbHttH_bbZZconstraint.color = 1
#     scalingbbHttH_bbZZconstraint.read()

#     plot = differentials.plotting.plots.MultiContourPlot(
#         'multicont_TopCtCb_scalingbbHttH_floatingBRs_constrainedbbZZ' + ('_asimov' if args.asimov else ''),
#         [
#             # scalingbbHttH,
#             scalingbbHttH_bbZZconstraint
#             ],
#         # x_min=top_x_min, x_max=top_x_max, y_min=top_y_min, y_max=top_y_max,
#         )
#     plot.draw()


@flag_as_option
def multicont_TopCtCb_scalingbbHttH(args):
    args = differentialutils.set_one_decay_channel(args, 'combWithHbb')
    y_coupling = 'cb'

    if args.asimov:
        scandir = 'out/Scan_May26_TopCtCb_combWithHbb_scalingbbHttH_asimov'
    else:
        # scandir = 'out/Scan_May26_TopCtCb_combWithHbb_scalingbbHttH'
        # scandir = 'out/Scan_May28_TopCtCb_combWithHbb_scalingbbHttH'
        scandir = 'out/Scan_May29_TopCtCb_combWithHbb_scalingbbHttH_0'

    scalingbbHttH = differentials.scans.Scan2D(
        'scalingbbHttH', x_coupling, y_coupling,
        scandir = scandir
        )
    scalingbbHttH.title = 'Combination (incl. bbH / BRs free)'
    scalingbbHttH.color = 1
    scalingbbHttH.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_TopCtCb_scalingbbHttH' + ('_asimov' if args.asimov else ''),
        [scalingbbHttH],
        # x_min=top_x_min, x_max=top_x_max, y_min=top_y_min, y_max=top_y_max,
        )
    plot.draw()


#____________________________________________________________________
@flag_as_option
def multicont_Top_scalingttH_floatingBRs(args):
    args = differentialutils.set_one_decay_channel(args, 'combWithHbb')
    scans = []
    
    if args.asimov:
        scandir = 'out/Scan_May23_Top_combWithHbb_scalingttH_floatingBRs_asimov'
    else:
        scandir = 'out/Scan_May24_Top_combWithHbb_scalingttH_floatingBRs'
        # scandir = 'out/Scan_May25_Top_combWithHbb_scalingttH_floatingBRs'

    scalingttH = differentials.scans.Scan2D(
        'scalingttH', x_coupling, y_coupling,
        scandir = scandir
        )
    scalingttH.title = 'Combination (incl. bbH / BRs free)'
    scalingttH.color = 1
    scalingttH.read()
    scans.append(scalingttH)
    
    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Top_scalingttH_floatingBRs' + ('_asimov' if args.asimov else ''),
        scans,
        x_min=top_x_min, x_max=top_x_max, y_min=top_y_min, y_max=top_y_max,
        )
    plot.y_SM = 0.0
    plot.draw()

@flag_as_option
def thetascan_Top_scalingttH_floatingBRs(args):
    args = differentialutils.set_one_decay_channel(args, 'combWithHbb')

    exp2D = differentials.scans.Scan2D(
        'exp2D', x_coupling, y_coupling,
        scandir = 'out/Scan_May23_Top_combWithHbb_scalingttH_floatingBRs_asimov'
        )
    exp2D.title = 'Combination (incl. bbH / BRs free)'
    exp2D.color = 1
    exp2D.read()

    thetas1, deltaNLLs1 = exp2D.get_thetas_along_x(0.0, 6.0, y=-0.2)
    thetas2, deltaNLLs2 = exp2D.get_thetas_along_y(-0.2, 0.1, x=6.0)
    thetas3, deltaNLLs3 = exp2D.get_thetas_along_x(6.0, 0.0, y=0.1)
    thetas = thetas1 + thetas2 + thetas3
    deltaNLLs = deltaNLLs1 + deltaNLLs2 + deltaNLLs3
    theta_graph_exp = differentials.plotting.pywrappers.Graph(
        differentials.plotting.plotting_utils.get_unique_rootname(),
        'theta_exp',
        thetas,
        deltaNLLs,
        color=1
        )
    theta_graph_exp.title = 'Expected'

    obs2D = differentials.scans.Scan2D(
        'obs2D', x_coupling, y_coupling,
        scandir = 'out/Scan_May24_Top_combWithHbb_scalingttH_floatingBRs'
        )
    obs2D.title = 'Combination (incl. bbH / BRs free)'
    obs2D.color = 1
    obs2D.read()

    thetas1, deltaNLLs1 = obs2D.get_thetas_along_y(-0.08, 0.1, x=1.2)
    thetas2, deltaNLLs2 = obs2D.get_thetas_along_x(1.2, 0.0, y=0.1)
    thetas = thetas1 + thetas2
    deltaNLLs = deltaNLLs1 + deltaNLLs2
    theta_graph_obs = differentials.plotting.pywrappers.Graph(
        differentials.plotting.plotting_utils.get_unique_rootname(),
        'theta_obs',
        thetas,
        deltaNLLs,
        color=1
        )
    theta_graph_obs.title = 'Observed'

    x_min = -0.2
    x_max = 0.75*math.pi
    y_min = 0.0
    y_max = 9.
    plot = differentials.plotting.plots.QuickPlot(
        'thetaplot_test',
        x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max
        )
    # plot.do_legend = False
    plot.leg._y1 = lambda c: 1. - c.GetTopMargin() - 0.07
    plot.x_title = '#theta = tan(c_{g}/#kappa_{t})'
    plot.y_title = '2#DeltaNLL'
    plot.add(theta_graph_exp, 'repr_dashed_line')
    plot.add(theta_graph_obs, 'repr_smooth_line')
    plot.draw()

    l = ROOT.TLine(x_min, 1.0, x_max, 1.0)
    ROOT.SetOwnership(l, False)
    l.SetLineColor(14)
    l.Draw()

    plot.wrapup()


@flag_as_option
def multicont_TopCtCb_scalingbbHttH_floatingBRs(args):
    args = differentialutils.set_one_decay_channel(args, 'combWithHbb')
    y_coupling = 'cb'

    # scandir = 'out/Scan_May26_TopCtCb_combWithHbb_scalingbbHttH_floatingBRs'
    # scandir = 'out/Scan_May28_TopCtCb_combWithHbb_scalingbbHttH_floatingBRs'
    # scandir = 'out/Scan_May28_TopCtCb_combWithHbb_scalingbbHttH_floatingBRs_0'
    scandir = 'out/Scan_May29_TopCtCb_combWithHbb_scalingbbHttH_floatingBRs_0'

    scalingbbHttH = differentials.scans.Scan2D(
        'scalingbbHttH', x_coupling, y_coupling,
        scandir = scandir
        )
    scalingbbHttH.title = 'Combination (incl. bbH / BRs free)'
    scalingbbHttH.color = 1
    scalingbbHttH.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_TopCtCb_scalingbbHttH_floatingBRs' + ('_asimov' if args.asimov else ''),
        [scalingbbHttH],
        # x_min=top_x_min, x_max=top_x_max, y_min=top_y_min, y_max=top_y_max,
        )
    plot.draw()


@flag_as_option
def ratioctcg_Top_floatingBRs(args):

    scan = differentials.scans.Scan(
        x_variable = 'cg', y_variable='nll',
        scandir = 'out/Scan_May26_TopPoints_combWithHbb_asimov',
        )
    scan.save_all_variables = True
    scan.read()

    for entry in scan.entries:
        print entry.cg, entry.ct, entry.nll


#____________________________________________________________________

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

    if args.asimov:
        scandir = LatestPaths.scan.top.reweighted.asimov.combWithHbb
    else:
        # scandir = LatestPaths.scan.top.reweighted.observed.combWithHbb
        scandir = 'out/Scan_May31_Top_combWithHbb_scalingttH_couplingdependentBRs'
    Top_combWithHbb = differentials.scans.Scan2D('combWithHbb', x_coupling, y_coupling, scandir=scandir)
    Top_combWithHbb.color = 1
    Top_combWithHbb.read()

    # obs = LatestBinning.obs_pth_ggH
    # obs.drop_bins_up_to_value(125.)

    if args.asimov:
        scandir = LatestPaths.scan.pth_ggH.asimov.combWithHbb
    else:
        scandir = LatestPaths.scan.pth_ggH.observed.combWithHbb
    combWithHbb = differentials.scans.DifferentialSpectrum('combWithHbb', scandir)
    combWithHbb.color = 1
    combWithHbb.no_overflow_label = True
    combWithHbb.draw_method = 'repr_point_with_vertical_bar'
    combWithHbb.set_sm(obstuple.combWithHbb.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    combWithHbb.read()

    # ws = 'out/workspaces_Dec11/combinedCard_Nov03_CouplingModel_TopHighPt_withTheoryUncertainties.root'
    # ws = 'out/workspaces_Mar06/combWithHbb_Top_reweighted_nominal.root'
    # ws = LatestPaths.ws.top.nominal.combWithHbb
    ws = 'out/workspaces_May31/combWithHbb_Top_reweighted_scalingttH_couplingdependentBRs.root'

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

@flag_as_option
def multicont_TopCtCb_lumi300fb(args):
    args = differentialutils.force_asimov(args)
    y_coupling = 'cb'

    combWithHbb = differentials.scans.Scan2D('combWithHbb', x_coupling, y_coupling,
        scandir = LatestPaths.scan.topctcb.reweighted.asimov.combWithHbb
        )
    combWithHbb.color = 1
    combWithHbb.title = '35.9 fb^{-1}'
    combWithHbb.read()

    lumi300fb = differentials.scans.Scan2D('lumi300fb', x_coupling, y_coupling,
        scandir = LatestPaths.scan.topctcb.lumi300fb
        )
    lumi300fb.color = 4
    lumi300fb.title = '300 fb^{-1}'
    # lumi300fb.contour_filter_method = 'max_distance_to_com'
    lumi300fb.read()
    
    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_TopCtCb_lumi300fb' + ('_asimov' if args.asimov else ''),
        [combWithHbb, lumi300fb],
        x_min = -1.8,
        x_max = 1.8,
        y_min = -20.,
        y_max = 20.,
        )
    # plot.only_1sigma_contours = True
    plot.draw()


@flag_as_option
def points_on_contour_TopCtCb(args):
    obs_name = 'pth_ggH'
    obstuple = LatestBinning.obstuple_pth_ggH
    y_coupling = 'cb'

    # Load 2D scan
    if args.asimov:
        scandir = LatestPaths.scan.topctcb.reweighted.asimov.combWithHbb
    else:
        # scandir = LatestPaths.scan.topctcb.reweighted.observed.combWithHbb
        scandir = 'out/Scan_May29_TopCtCb_combWithHbb_scalingbbHttH_couplingdependentBRs_0'

    TopCtCb_combWithHbb = differentials.scans.Scan2D('combWithHbb', x_coupling, y_coupling, scandir = scandir)
    TopCtCb_combWithHbb.color = 1
    TopCtCb_combWithHbb.read()

    # Load pt combination
    if args.asimov:
        scandir = LatestPaths.scan.pth_ggH.asimov.combWithHbb
    else:
        scandir = LatestPaths.scan.pth_ggH.observed.combWithHbb

    combWithHbb = differentials.scans.DifferentialSpectrum('combWithHbb', scandir)
    combWithHbb.color = 1
    combWithHbb.no_overflow_label = True
    combWithHbb.draw_method = 'repr_point_with_vertical_bar'
    combWithHbb.set_sm(obstuple.combWithHbb.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    combWithHbb.read()

    # Load ws to get parametrization from
    # ws = LatestPaths.ws.topctcb.nominal.combWithHbb
    ws = 'out/workspaces_May29/combWithHbb_TopCtCb_reweighted_scalingbbHttH_couplingdependentBRs.root'

    # ======================================
    # Load into plot

    plot = differentials.plotting.plots.BottomPanelPlotWithParametrizations('points_on_contour_TopCtCb')
    plot.scan2D = TopCtCb_combWithHbb
    plot.ws_file = ws
    plot.ptspectrum = combWithHbb
    plot.obs = obstuple.combWithHbb

    plot.get_points_method = 'extrema_y_and_x_max'

    # plot.topctcb_y_min = topctcb_y_min
    plot.top_y_max = 150.
    plot.legend.set(
        x1 = lambda c: c.GetLeftMargin()+0.05,
        y1 = lambda c: 1.-c.GetTopMargin()-0.35,
        x2 = lambda c: c.GetLeftMargin()+0.38,
        y2 = lambda c: 1.-c.GetTopMargin()-0.02,
        )

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



@flag_as_option
def thetaplot_test(args):
    differentials.plotting.canvas.c.resize_temporarily(850, 800)

    # Get the 2D scan
    scalingttH = differentials.scans.Scan2D('scalingttH', x_coupling, y_coupling,
        scandir = 'out/Scan_May22_Top_combWithHbb_scalingttH_couplingdependentBRs_asimov'
        )
    scalingttH.color = 1
    scalingttH.read()

    # Read the observed from the 2D scan; take a minimum per slice
    theta_graph = scalingttH.get_thetas()
    theta_graph.draw_style = 'repr_smooth_line'
    # theta_graph.title = '{0} observed; ({1:.2f} - {2:.2f}) @ 68% CL'.format(
    #     differentials.core.standard_titles[kappa],
    #     theta_graph.unc.left_bound, theta_graph.unc.right_bound
    #     )

    plot = differentials.plotting.plots.QuickPlot(
        'thetaplot_test',
        x_min = -0.26*math.pi, x_max = 0.76*math.pi, y_min = 0.0, y_max = 20.
        )
    plot.do_legend = False
    plot.x_title = '#theta = tan(c_{g}/#kappa_{t})'
    plot.y_title = '2#DeltaNLL'
    plot.add(theta_graph, 'repr_smooth_line')
    plot.draw()
    plot.wrapup()

    # plot.scans.append(expected)
    # plot.manual_graphs.append(observed1D)
    # plot.x_title = differentials.core.standard_titles[kappa]
    # if kappa == 'kappab':
    #     plot.x_min = -7.
    #     plot.x_max = 9.
    # else:
    #     plot.x_min = -20.
    #     plot.x_max = 20.
    # plot.leg.SetNColumns(1)
    # plot.leg._y1 = lambda c: 1. - c.GetTopMargin() - 0.20
    # plot.draw()
    # plot.wrapup()

