import ROOT
import plotting_utils as utils
import pywrappers
from canvas import c

import differentials
import differentials.core as core
import plots

import logging
from math import isnan, isinf, log10, sqrt
from array import array
from collections import namedtuple


class CorrelationMatrixPlot(plots.PlotBase):

    def __init__(self, plotname):
        super(CorrelationMatrixPlot, self).__init__(plotname)
        self.x_title = 'x'
        self.y_title = None
        self.corrmat = None
        self.binlabels = []
        self.n_bins = None

        self.manybin_cutoff = 20

        self.add_zaxis_label = True
        self.zaxis_label_text = 'Correlation'

        self.marker_size = 1.4

    def get_bin_labels(self):
        if len(self.binlabels) == 0:
            ['bin{0}'.format(i) for i in xrange(self.n_bins)]
        else:
            return self.binlabels

    def draw(self):
        if self.n_bins is None:
            self.n_bins = len(self.corrmat)

        c.Clear()
        c.set_margins(
            LeftMargin  = 0.21,
            TopMargin   = 0.08,
            RightMargin = 0.14,
            BottomMargin = 0.19,
            )
        utils.set_color_palette('correlation_matrix')

        H = ROOT.TH2D(
            utils.get_unique_rootname(),
            # '#scale[0.85]{{Bin-to-bin correlation matrix for {0}}}'.format(observableName),
            '',
            self.n_bins, 0., self.n_bins,
            self.n_bins, 0., self.n_bins
            )
        ROOT.SetOwnership(H, False)
        H.SetContour(100)

        for i_row in xrange(self.n_bins):
            for i_col in xrange(self.n_bins):
                H.SetBinContent(i_col+1, i_row+1, self.corrmat[i_row][i_col])
        H.GetZaxis().SetRangeUser(-1.0,1.0)

        binning_labels = self.get_bin_labels()
        for i in xrange(self.n_bins):
            if self.n_bins < self.manybin_cutoff or i % int(0.1*self.n_bins) == 0:
                H.GetXaxis().SetBinLabel(i+1, binning_labels[i])
                H.GetYaxis().SetBinLabel(i+1, binning_labels[i])

        H.GetXaxis().SetTitle(self.x_title)

        H.GetXaxis().SetTitleSize(0.05)
        H.GetXaxis().SetTitleOffset(1.8)
        # H.GetXaxis().SetLabelSize(0.045)
        # H.GetYaxis().SetLabelSize(0.045)

        # H.GetXaxis().SetTitleSize(0.06)
        # H.GetXaxis().SetTitleOffset(1.6)
        H.GetXaxis().SetLabelSize(0.05)
        H.GetYaxis().SetLabelSize(0.05)

        if self.n_bins < self.manybin_cutoff:
            H.Draw('COLZ TEXT')
        else:
            H.Draw('COLZ')

        ROOT.gStyle.SetHistMinimumZero() # To draw the "0", otherwise ROOT leaves it empty
        ROOT.gStyle.SetPaintTextFormat('1.2g')

        H.SetMarkerSize(self.marker_size)

        pywrappers.CMS_Latex_type('Supplementary').Draw()
        pywrappers.CMS_Latex_lumi().Draw()

        if self.add_zaxis_label:
            l = differentials.plotting.pywrappers.Latex(
                lambda c: 1.0 - 0.01,
                lambda c: 1.0 - c.GetTopMargin(),
                self.zaxis_label_text
                )
            l.SetTextFont(42)
            l.SetNDC()
            l.SetTextAngle(-90)
            l.SetTextSize(0.04)
            l.SetTextAlign(13)
            l.Draw()

        self.H = H


class CorrelationMatrixFromCombinePlot(CorrelationMatrixPlot):
    """docstring for CorrelationMatrixFromCombinePlot"""
    def __init__(self, plotname, ws):
        super(CorrelationMatrixFromCombinePlot, self).__init__(plotname)
        self.ws = ws
        
    def get_pois_from_ws(self):
        self.pois = [ p for p in core.read_set(self.ws, 'POI') if not p=='MH' ]
        self.pois.sort()

    def get_correlation_matrix_from_ws(self):
        self.n_pois = len(self.pois)
        with core.openroot(self.ws) as ws_fp:
            self.fit = ws_fp.Get('fit')

        corrmat = [[ 999. for j in xrange(self.n_pois)] for i in xrange(self.n_pois) ]
        for i, poi1 in enumerate(self.pois):
            for j, poi2 in enumerate(self.pois):
                if i == j:
                    corrmat[i][j] = 1.0
                else:
                    corrmat[i][j] = self.fit.correlation(poi1, poi2)
        self.corrmat = corrmat

    def get_bin_labels(self):
        if len(self.binlabels) == 0:
            return self.pois
        else:
            return self.binlabels

    def draw(self):
        self.n_bins = self.n_pois
        if self.corrmat is None:
            self.get_correlation_matrix_from_ws()
        super(CorrelationMatrixFromCombinePlot, self).draw()

        self.H.GetYaxis().SetTitle(self.H.GetXaxis().GetTitle())
        self.H.GetYaxis().SetTitleSize(self.H.GetXaxis().GetTitleSize())
        self.H.GetYaxis().SetTitleOffset(2.15)

