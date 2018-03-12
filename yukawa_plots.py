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
import logging
from collections import namedtuple

from time import strftime
datestr = strftime( '%b%d' )

import ROOT
from array import array

########################################
# Main
########################################

x_coupling = 'kappac'
y_coupling = 'kappab'
titles = { 'kappac': '#kappa_{c}', 'kappab' : '#kappa_{b}' }


Nominal = namedtuple('Nominal', ['hgg', 'hzz', 'combination'])
nominal_scans = Nominal(
    hgg = 'out/Scan_Yukawa_Feb24_hgg',
    hzz = 'out/Scan_Yukawa_Feb24_hzz',
    combination = 'out/Scan_Yukawa_Feb24_combination',
    )
nominal_scans_asimov = Nominal(
    # hgg = 'out/Scan_Yukawa_Feb24_hgg_0',
    # # hzz = 'out/Scan_Yukawa_Feb24_hzz_0',
    # hzz = 'out/Scan_Yukawa_Feb24_hzz_asimov',
    # combination = 'out/Scan_Yukawa_Feb24_combination_0',
    hgg = 'out/Scan_Yukawa_Nov03_hgg_asimov',
    hzz = 'out/Scan_Yukawa_Nov03_hzz_asimov',
    combination = 'out/Scan_Yukawa_Nov03_asimov',
    # out/Scan_Yukawa_Nov03_asimov
    # out/Scan_Yukawa_Nov03_hgg_asimov
    # out/Scan_Yukawa_Nov03_hzz_asimov
    )
def get_nominal(args):
    if args.asimov:
        return nominal_scans_asimov
    else:
        return nominal_scans


@flag_as_option
def multicont_Yukawa_Mar09(args):
    multicont_Yukawa_Mar09_reweighted(args)
    multicont_Yukawa_Mar09_unreweighted(args)

@flag_as_option
def multicont_Yukawa_Mar09_reweighted(args):
    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling,
        scandir = 'out/Scan_Yukawa_Mar09_combination_asimov' if args.asimov else 'out/Scan_Yukawa_Mar09_combination'
        )
    combination.color = 1
    combination.read()
    
    hgg = differentials.scans.Scan2D('hgg', x_coupling, y_coupling,
        scandir = 'out/Scan_Yukawa_Mar09_hgg_asimov' if args.asimov else 'out/Scan_Yukawa_Mar09_hgg'
        )
    hgg.color = 2
    hgg.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_Mar09_reweighted' + ('_asimov' if args.asimov else ''),
        [combination, hgg],
        x_min=-25, x_max=25, y_min=-9, y_max=12
        )
    plot.draw()

@flag_as_option
def multicont_Yukawa_Mar09_unreweighted(args):
    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling,
        scandir = 'out/Scan_Yukawa_Mar09_combination_asimov_0' if args.asimov else 'out/Scan_Yukawa_Mar09_combination_0'
        )
    combination.color = 1
    combination.read()
    
    hgg = differentials.scans.Scan2D('hgg', x_coupling, y_coupling,
        scandir = 'out/Scan_Yukawa_Mar09_hgg_asimov_0' if args.asimov else 'out/Scan_Yukawa_Mar09_hgg_0'
        )
    hgg.color = 2
    hgg.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_Mar09_unreweighted' + ('_asimov' if args.asimov else ''),
        [combination, hgg],
        x_min=-25, x_max=25, y_min=-9, y_max=12
        )
    plot.draw()



@flag_as_option
def multicont_Yukawa_nominal(args):

    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling,
        # scandir=LatestPaths.scan_combined_Yukawa_old
        scandir = get_nominal(args).combination
        )
    combination.color = 1
    combination.read()
    
    hgg = differentials.scans.Scan2D('hgg', x_coupling, y_coupling,
        # scandir=LatestPaths.scan_hgg_Yukawa_old
        scandir = get_nominal(args).hgg
        )
    hgg.color = 2
    hgg.read()

    hzz = differentials.scans.Scan2D('hzz', x_coupling, y_coupling,
        # scandir=LatestPaths.scan_hzz_Yukawa_old
        scandir = get_nominal(args).hzz
        )
    hzz.color = 4
    hzz.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_nominal' + ('_asimov' if args.asimov else ''),
        [combination, hgg, hzz],
        x_min=None, x_max=None, y_min=None, y_max=None
        )
    plot.x_min = -25
    plot.x_max = 25
    plot.y_min = -9
    plot.y_max = 12
    plot.draw()


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
        x_min=None, x_max=None, y_min=None, y_max=None
        )
    plot.draw()


@flag_as_option
def multicont_Yukawa_profiledTotalXS(args):

    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling,
        # scandir=LatestPaths.scan_combined_Yukawa_asimov
        scandir=LatestPaths.scan_combined_Yukawa_reweighted_asimov
        )
    combination.color = 1
    combination.title = 'Nominal'
    combination.read()

    profiled_total_XS = differentials.scans.Scan2D('profiled_total_XS', x_coupling, y_coupling,
        scandir=LatestPaths.scan_combined_Yukawa_profiledTotalXS_asimov
        )
    profiled_total_XS.color = 4
    profiled_total_XS.title = '#sigma_{tot} profiled'
    profiled_total_XS.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_profiledTotalXS',
        [combination, profiled_total_XS],
        x_min=None, x_max=None, y_min=None, y_max=None
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
        x_min=None, x_max=None, y_min=None, y_max=None
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
        x_min=None, x_max=None, y_min=None, y_max=None
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
        x_min=None, x_max=None, y_min=None, y_max=None
        )
    plot.draw()



@flag_as_option
def multicont_Yukawa_highLumi(args):

    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling,
        scandir=LatestPaths.scan_combined_Yukawa_asimov
        )
    combination.color = 1
    combination.title = '35.9 fb^{-1}'
    combination.read()

    lumi_300 = differentials.scans.Scan2D('lumi_300', x_coupling, y_coupling,
        scandir=LatestPaths.scan_combined_Yukawa_lumiStudy_asimov
        )
    lumi_300.color = 4
    lumi_300.title = '300 fb^{-1}'
    lumi_300.read()

    raise NotImplementedError('Need to multiply values from 300 fb-1')
    lumi_3000 = differentials.scans.Scan2D('lumi_3000', x_coupling, y_coupling,
        scandir=LatestPaths.scan_combined_Yukawa_lumiStudy_asimov
        )
    lumi_3000.color = 4
    lumi_3000.title = '3000 fb^{-1}'
    lumi_3000.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_highLumi',
        [combination, lumi_300, lumi_3000],
        x_min=None, x_max=None, y_min=None, y_max=None
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
        x_min=None, x_max=None, y_min=None, y_max=None
        )
    plot.draw()


# ======================================
# Other plots

@flag_as_option
def points_on_contour_Yukawa(args):

    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling,
        scandir=get_nominal(args).combination
        )
    combination.color = 1
    combination.read()

    obs = LatestBinning.obs_pth_ggH
    obs.drop_bins_up_to_value(125.)

    obs_name = 'pth_ggH'
    scandir = 'out/Scan_pth_ggH_Feb06_combination_asimov' if args.asimov else 'out/Scan_pth_ggH_Feb06_combination'
    pth_ggH_combination = differentials.scans.DifferentialSpectrum('combination', scandir=scandir)
    pth_ggH_combination.POIs = differentials.core.get_POIs_oldstyle_scandir(scandir)
    pth_ggH_combination.color = 1
    pth_ggH_combination.no_overflow_label = True
    pth_ggH_combination.draw_method = 'repr_point_with_vertical_bar'
    pth_ggH_combination.set_sm(obs.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=False))
    pth_ggH_combination.drop_bins_up_to_value(125.)
    pth_ggH_combination.read()

    # ws = LatestPaths.ws_combined_TopHighPt
    # ws = LatestPaths.ws_combined_Yukawa
    ws = 'out/workspaces_Feb01/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties.root'


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
    one_kappa_scan_kappab(args)
    one_kappa_scan_kappac(args)


@flag_as_option
def one_kappa_scan_kappab(args):
    differentials.plotting.canvas.c.resize_temporarily(850, 800)
    scans = []

    # observed = differentials.scans.Scan('kappab', scandir=LatestPaths.scan_combined_Yukawa_oneKappa_kappab)
    # observed.color = 38
    # observed.title = '#kappa_{b} observed'
    # observed.read()
    # observed.create_uncertainties()
    # observed.title += '; ({0:.2f} - {1:.2f}) @ 68% CL'.format(observed.unc.left_bound, observed.unc.right_bound)
    # scans.append(observed)

    expected = differentials.scans.Scan('kappab', scandir='out/Scan_Yukawa_Feb07_combined_oneKappa_kappab_asimov')
    expected.color = 4
    expected.title = '#kappa_{b} expected'
    expected.read()
    expected.create_uncertainties()
    expected.title += '; ({0:.2f} - {1:.2f}) @ 68% CL'.format(expected.unc.left_bound, expected.unc.right_bound)
    scans.append(expected)

    plot = differentials.plotting.plots.MultiScanPlot('onekappascan_kappab')
    plot.scans.extend(scans)
    plot.x_min = -6.
    plot.x_max = 6.
    plot.x_title = '#kappa_{b}'

    plot.leg.SetNColumns(1)
    plot.leg._y1 = lambda c: 1. - c.GetTopMargin() - 0.20

    plot.draw()
    plot.wrapup()


@flag_as_option
def one_kappa_scan_kappac(args):
    differentials.plotting.canvas.c.resize_temporarily(850, 800)
    scans = []

    # observed = differentials.scans.Scan('kappac', scandir=LatestPaths.scan_combined_Yukawa_oneKappa_kappac)
    # observed.color = 38
    # observed.title = '#kappa_{c} observed'
    # observed.read()
    # observed.create_uncertainties()
    # observed.title += '; ({0:.2f} - {1:.2f}) @ 68% CL'.format(observed.unc.left_bound, observed.unc.right_bound)
    # observed.draw_style = 'repr_smooth_line'
    # scans.append(observed)

    expected = differentials.scans.Scan('kappac', scandir='out/Scan_Yukawa_Feb07_combined_oneKappa_kappac_asimov')
    expected.color = 4
    expected.title = '#kappa_{c} expected'
    expected.read()
    expected.create_uncertainties()
    expected.title += '; ({0:.2f}, {1:.2f}) @ 68% CL'.format(expected.unc.left_bound, expected.unc.right_bound)
    expected.draw_style = 'repr_smooth_line'
    scans.append(expected)

    plot = differentials.plotting.plots.MultiScanPlot('onekappascan_kappac')
    plot.scans.extend(scans)
    plot.x_min = -15.
    plot.x_max = 15.
    plot.x_title = '#kappa_{c}'

    plot.leg.SetNColumns(1)
    plot.leg._y1 = lambda c: 1. - c.GetTopMargin() - 0.20

    plot.draw()
    plot.wrapup()


