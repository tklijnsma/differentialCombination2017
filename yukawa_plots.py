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
import differentialutils
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

yukawa_x_min = -35
yukawa_x_max = 35
yukawa_y_min = -11
yukawa_y_max = 12

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
scandict_G.G0B = 'out/Scan_Yukawa_May11_combination_G0B_asimov_0'
scandict_G.G1A = 'out/Scan_Yukawa_May11_combination_G1A_asimov'
scandict_G.G1B = 'out/Scan_Yukawa_May11_combination_G1B_asimov'
scandict_G.G2A = 'out/Scan_Yukawa_May11_combination_G2A_asimov'
scandict_G.G1BKV = 'out/Scan_Yukawa_May11_combination_G1BKV_asimov'

@flag_as_option
def multicont_Yukawa_G0(args):
    G0A = differentials.scans.Scan2D('G0A', x_coupling, y_coupling, scandir = scandict_G.G0A)
    G0A.color = 1
    G0A.read()

    G0B = differentials.scans.Scan2D('G0B', x_coupling, y_coupling, scandir = scandict_G.G0B)
    G0B.color = 14
    G0B.read()

    G1A = differentials.scans.Scan2D('G1A', x_coupling, y_coupling, scandir = scandict_G.G1A)
    G1A.color = 8
    G1A.read()

    G1B = differentials.scans.Scan2D('G1B', x_coupling, y_coupling, scandir = scandict_G.G1B)
    G1B.color = 3
    G1B.read()

    # G1BKV = differentials.scans.Scan2D('G1BKV', x_coupling, y_coupling, scandir = scandict_G.G1BKV)
    # G1BKV.color = 9
    # G1BKV.read()

    G2A = differentials.scans.Scan2D('G2A', x_coupling, y_coupling, scandir = scandict_G.G2A)
    G2A.color = 2
    G2A.read()

    plot = differentials.plotting.plots.MultiContourPlot(
        'multicont_Yukawa_G0',
        [G0A, G0B, G1A, G1B, G2A],
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
        'multicont_Yukaw_paper' + ('_asimov' if args.asimov else ''),
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
    scandict_2D = LatestPaths.scan.yukawa.reweighted.asimov if args.asimov else LatestPaths.scan.yukawa.reweighted.observed
    combination = differentials.scans.Scan2D('combination', x_coupling, y_coupling, scandir=scandict_2D.combination)
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
    ws = LatestPaths.ws.yukawa.nominal.combination

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


