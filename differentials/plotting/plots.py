import ROOT
import plotting_utils as utils
import pywrappers
from canvas import c

import differentials.parametrization

import logging
from math import isnan, isinf, log10, sqrt
from array import array
from collections import namedtuple


class PlotBase(object):
    """docstring for PlotBase"""
    def __init__(self, plotname):
        self.plotname = plotname
        c.Clear()
        c.set_margins()
        
    def pre_draw(self):
        pass

    def post_draw(self):
        pass

    def draw(self):
        self.pre_draw()

    def wrapup(self):
        self.post_draw()
        self.save()

    def save(self):
        c.save(self.plotname)


class MultiContourPlot(PlotBase):
    """docstring for MultiContourPlot"""
    def __init__(self, plotname, scans, x_min=None, x_max=None, y_min=None, y_max=None):
        super(MultiContourPlot, self).__init__(plotname)

        self.scans = scans

        self.x_min = x_min
        self.y_min = y_min
        self.x_max = x_max
        self.y_max = y_max

        self.draw_individual_contours=True
        self.x_SM = 1.0
        self.y_SM = 1.0

        self.x_title = scans[0].x_title
        self.y_title = scans[0].y_title
        self.base = pywrappers.Base(x_title=self.x_title, y_title=self.y_title)

        self.legend = pywrappers.Legend(
            lambda c: c.GetLeftMargin() + 0.01,
            lambda c: c.GetBottomMargin() + 0.02,
            lambda c: 1 - c.GetRightMargin() - 0.01,
            lambda c: c.GetBottomMargin() + 0.09
            )

    def pre_draw(self):
        super(MultiContourPlot, self).pre_draw()
        self.base.Draw()

    def draw(self):
        super(MultiContourPlot, self).draw()

        self.histograms = []
        for scan in self.scans:
            histogram2D = scan.to_hist()
            histogram2D._legend = self.legend
            histogram2D.Draw('repr_contours')
            self.histograms.append(histogram2D)

        self.legend.Draw()

        if self.x_min is None: self.x_min = min([H.x_min() for H in self.histograms])
        if self.y_min is None: self.y_min = min([H.y_min() for H in self.histograms])
        if self.x_max is None: self.x_max = max([H.x_max() for H in self.histograms])
        if self.y_max is None: self.y_max = max([H.y_max() for H in self.histograms])
        self.base.set(x_min=self.x_min, y_min=self.y_min, x_max=self.x_max, y_max=self.y_max)

        self.base.GetXaxis().SetTitle(self.x_title)
        self.base.GetYaxis().SetTitle(self.y_title)
        self.base.GetXaxis().SetTitleSize(0.06)
        self.base.GetXaxis().SetLabelSize(0.05)
        self.base.GetYaxis().SetTitleSize(0.06)
        self.base.GetYaxis().SetLabelSize(0.05)

        SM_point = pywrappers.Point(self.x_SM, self.y_SM)
        SM_point.Draw('repr_SM_point')

        pywrappers.CMS_Latex_type().Draw()
        pywrappers.CMS_Latex_lumi().Draw()

        pywrappers.ContourDummyLegend(
            c.GetLeftMargin() + 0.01,
            1. - c.GetTopMargin() - 0.1,
            1. - c.GetRightMargin() - 0.01,
            1. - c.GetTopMargin() - 0.01,
            ).Draw()

        self.wrapup()

    def wrapup(self):
        c.Update()
        c.RedrawAxis()
        self.save()        
        if self.draw_individual_contours:
            for scan in self.scans:
                scan.plot(self.plotname + '_' + scan.name)


