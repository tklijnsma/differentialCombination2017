import differentials.core as core
import differentials.logger as logger

from canvas import c
import plotting_utils as utils
import multipanel
import pywrappers

from array import array

import ROOT


def plot_spectra(plotname, spectra, obs_name, obs_title, obs_unit=None, inplace=True):
    plot = multipanel.BottomPanelPlot(plotname)
    plot.x_title = obs_title + (' ({0})'.format(obs_unit) if not(obs_unit is None) else '')
    plot.y_title_top = '#Delta#sigma/#Delta{0} (pb{1})'.format(
        obs_title,
        '/' + obs_unit if not(obs_unit is None) else ''
        )
    plot.y_title_bottom = '#mu'

    leg = pywrappers.Legend(
        lambda c: c.GetLeftMargin() + 0.01,
        lambda c: 1 - c.GetTopMargin() - 0.10,
        lambda c: 1 - c.GetRightMargin() - 0.01,
        lambda c: 1 - c.GetTopMargin()
        )

    for spectrum in spectra:
        if not spectrum.smxs_set:
            raise RuntimeError('SM cross sections were not provided for spectrum')
        spectrum.read() 

    if spectra[0].last_bin_is_overflow():
        # Make sure overflow bins are aligned
        x_max = max([ 2*s.binning()[-2]-s.binning()[-3] for s in spectra ])
        for spectrum in spectra:
            spectrum.hard_x_max = x_max

    plot.make_SM_line(spectra, leg=leg)
    for spectrum in spectra:
        plot.add_bottom(spectrum.to_hist(), spectrum.draw_method)
        plot.add_top(spectrum.to_hist_xs(), spectrum.draw_method, leg=leg)

    plot.make_labels_for_overflow_spectra(spectra, obs_title)

    plot.add_top(leg, '')

    if inplace:
        plot.draw()
    else:
        return plot


def plot_multi_contour(
        plotname, scans,
        x_min=None, x_max=None, y_min=None, y_max=None,
        x_SM=1.0, y_SM=1.0,
        draw_individual_contours=True,
        inplace=True
        ):

    c.Clear()
    c.set_margins_2D()

    x_title = scans[0].x_title
    y_title = scans[0].y_title
    base = pywrappers.Base(x_title=x_title, y_title=y_title)
    base.Draw()

    leg = pywrappers.Legend(
        lambda c: c.GetLeftMargin() + 0.01,
        lambda c: c.GetBottomMargin() + 0.02,
        lambda c: 1 - c.GetRightMargin() - 0.01,
        lambda c: c.GetBottomMargin() + 0.09
        )

    histograms = []
    for scan in scans:
        histogram2D = scan.to_hist()
        histogram2D._legend = leg
        histogram2D.Draw('repr_contours')
        histograms.append(histogram2D)

    leg.Draw()

    if x_min is None: x_min = min([H.x_min() for H in histograms])
    if y_min is None: y_min = min([H.y_min() for H in histograms])
    if x_max is None: x_max = max([H.x_max() for H in histograms])
    if y_max is None: y_max = max([H.y_max() for H in histograms])
    base.set(x_min=x_min, y_min=y_min, x_max=x_max, y_max=y_max)

    base.GetXaxis().SetTitle(x_title)
    base.GetYaxis().SetTitle(y_title)
    base.GetXaxis().SetTitleSize(0.06)
    base.GetXaxis().SetLabelSize(0.05)
    base.GetYaxis().SetTitleSize(0.06)
    base.GetYaxis().SetLabelSize(0.05)

    SM_point = pywrappers.Point(x_SM, y_SM, color=12, marker_style=21)
    SM_point.Draw()

    pywrappers.CMS_Latex_type().Draw()
    pywrappers.CMS_Latex_lumi().Draw()

    pywrappers.ContourDummyLegend(
        c.GetLeftMargin() + 0.01,
        1. - c.GetTopMargin() - 0.1,
        1. - c.GetRightMargin() - 0.01,
        1. - c.GetTopMargin() - 0.01,
        ).Draw()

    c.Update()
    c.RedrawAxis()
    if inplace:
        c.save(plotname)
        if draw_individual_contours:
            for scan in scans:
                scan.plot(plotname + '_' + scan.name)
