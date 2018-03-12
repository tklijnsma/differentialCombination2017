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

    def __init__(self, plotname, ws=None):
        super(CorrelationMatrixPlot, self).__init__(plotname)
        self.ws = ws

        self.x_title = 'x'
        self.y_title = 'y'

    def get_fit_from_ws(self):
        with core.openroot(self.ws) as ws_fp:
            self.fit = ws_fp.Get('fit')
        self.pois = [ p for p in core.read_set(self.ws, 'POI') if not p=='MH' ]
        self.pois.sort()
        self.n_pois = len(self.pois)

    def get_correlation_matrix_from_ws(self):
        self.get_fit_from_ws()
        corrmat = [[ 999. for j in xrange(self.n_pois)] for i in xrange(self.n_pois) ]
        for i, poi1 in enumerate(self.pois):
            for j, poi2 in enumerate(self.pois):
                if i == j:
                    corrmat[i][j] = 1.0
                else:
                    corrmat[i][j] = self.fit.correlation(poi1, poi2)
        return corrmat

    def get_bin_labels(self):
        return self.pois

    def draw(self):
        corrmat = self.get_correlation_matrix_from_ws()

        c.Clear()
        c.set_margins(
            TopMargin   = 0.08,
            RightMargin = 0.14,
            BottomMargin = 0.17,
            )
        utils.set_color_palette('correlation_matrix')

        H = ROOT.TH2D(
            utils.get_unique_rootname(),
            # '#scale[0.85]{{Bin-to-bin correlation matrix for {0}}}'.format(observableName),
            '',
            self.n_pois, 0., self.n_pois,
            self.n_pois, 0., self.n_pois
            )
        ROOT.SetOwnership(H, False)
        H.SetContour(100)

        for i_row in xrange(self.n_pois):
            for i_col in xrange(self.n_pois):
                H.SetBinContent(i_col+1, i_row+1, corrmat[i_row][i_col])
        H.GetZaxis().SetRangeUser(-1.0,1.0)

        binning_labels = self.get_bin_labels()
        for i in xrange(self.n_pois):
            H.GetXaxis().SetBinLabel(i+1, binning_labels[i])
            H.GetYaxis().SetBinLabel(i+1, binning_labels[i])
        H.GetXaxis().SetTitle(self.x_title)
        H.GetXaxis().SetTitleSize(0.05)
        H.GetXaxis().SetTitleOffset(1.6)
        H.GetXaxis().SetLabelSize(0.045)
        H.GetYaxis().SetLabelSize(0.045)

        H.Draw('COLZ TEXT')

        ROOT.gStyle.SetHistMinimumZero() # To draw the "0", otherwise ROOT leaves it empty
        ROOT.gStyle.SetPaintTextFormat('1.2g')

        pywrappers.CMS_Latex_type().Draw()
        pywrappers.CMS_Latex_lumi().Draw()

        self.wrapup()

