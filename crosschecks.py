from OptionHandler import flag_as_option, flag_as_parser_options

import os, logging, copy, sys
from array import array

import ROOT

import LatestPaths, LatestBinning
import differentials
import differentialutils

from differentials.plotting.canvas import c


class Spectrum(object):
    """docstring for Spectrum"""
    def __init__(self, name, values, isSM, title=None, color=4, line_style=1):
        self.values = values
        self.name = name
        self.isSM = isSM
        self.color = color
        self.line_style = line_style
        if title is None:
            self.title = name
        else:
            self.title = title

class NormalizationCrossCheck(object):
    """docstring for NormalizationCrossCheck"""
    color_cycle = differentials.plotting.canvas.global_color_cycle
 
    def __init__(self):
        self.logscale = True
        self.last_bin_is_overflow = False
        self.spectra = []

    def add_spectrum(self, name, values, isSM=False, color=None, line_style=1):
        if len(values) < self.n_bins:
            raise ValueError(
                'The passed spectrum \'{0}\' has {1} bins, but the number of passed bins is {2}'.format(
                    name, len(values), self.n_bins) +
                '\n    values: {0}'.format(values) +
                '\n    bin_boundaries: {0}'.format(self.bin_boundaries)
                )
        elif len(values) > self.n_bins:
            Commands.Warning(
                'The passed spectrum \'{0}\' has {1} bins, but the number of passed bins is {2}'.format(
                    name, len(values), self.n_bins) +
                '\n    Will keep only the first {0} bins'.format(self.n_bins)
                )
            values = values[:self.n_bins]

        if color is None:
            color = next(self.color_cycle)
        self.spectra.append(Spectrum(name, values, isSM, color=color, line_style=line_style))

    def set_bin_boundaries(self, bin_boundaries, add_overflow=False):
        self.bin_boundaries = bin_boundaries
        if add_overflow:
            self.n_bins = len(self.bin_boundaries)
            self.last_bin_is_overflow = True
            self.bin_boundaries.append(10000.)
        else:
            self.n_bins = len(self.bin_boundaries)-1
    
    def get_SM(self):
        return [s for s in self.spectra if s.isSM][0]

    def plot_spectra(self, tag=None):
        SM = self.get_SM()

        x_min = self.bin_boundaries[0]
        if self.last_bin_is_overflow:
            x_max = self.bin_boundaries[-2] + (self.bin_boundaries[-2]-self.bin_boundaries[-3])
        else:
            x_max = self.bin_boundaries[-1]

        if self.logscale:
            y_minAbs = min([ y for y in SM.values if y >= 0.0 ])
            y_maxAbs = max(SM.values)
            y_min = 0.01 * y_minAbs
            y_max = 2.0 * y_maxAbs
        else:
            y_minAbs = min(SM.values)
            y_maxAbs = max(SM.values)
            y_min = y_minAbs - 0.1*(y_maxAbs-y_minAbs)
            y_max = y_maxAbs + 0.1*(y_maxAbs-y_minAbs)

        c.Clear()
        c.set_margins()

        base = differentials.plotting.pywrappers.Base(
            x_min = x_min,
            x_max = x_max,
            y_min = y_min,
            y_max = y_max,
            x_title = 'p_{T} (GeV)',
            y_title = '#sigma (pb/GeV)',
            )
        base.Draw()

        leg = ROOT.TLegend(
            1 - c.GetRightMargin() - 0.3,
            1 - c.GetTopMargin()   - 0.25,
            1 - c.GetRightMargin(),
            1 - c.GetTopMargin() 
            )
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)

        for spectrum in self.spectra:
            H = self.spectrum_to_hist(spectrum)
            H.Draw('HISTSAME')
            leg.AddEntry( H.GetName(), spectrum.title, 'l' )

        leg.Draw()

        if self.logscale: c.SetLogy()

        outname = 'normalizationCrosscheck_reimpl'
        if not(tag is None):
            outname += '_' + tag
        c.save(outname)


    def spectrum_to_hist(self, s):
        H = ROOT.TH1F(
            s.name+'_'+ differentials.plotting.plotting_utils.get_unique_rootname(), s.title,
            self.n_bins, array('f', self.bin_boundaries)
            )
        ROOT.SetOwnership( H, False )
        for i_bin in xrange(self.n_bins):
            H.SetBinContent(i_bin+1, s.values[i_bin])
        H.SetLineWidth(2)
        H.SetLineColor(s.color)
        H.SetLineStyle(s.line_style)
        return H



@flag_as_option
def normalization_xcheck_yukawa(args):

    normalizationCrossCheck = NormalizationCrossCheck()

    from LatestBinning import obs_pth_ggH as obs_pth
    obs_pth.drop_bins_up_to_value(120.)

    normalizationCrossCheck.set_bin_boundaries(obs_pth.binning, add_overflow=False)
    bin_boundaries = normalizationCrossCheck.bin_boundaries # Convenience

    # Add SM
    SMXSs = obs_pth.crosssection_over_binwidth()
    normalizationCrossCheck.add_spectrum('YR4, shape Vitt.', SMXSs, isSM=True, color=4)

    # # Add raw from Pier
    # SM_Yukawa = TheoryFileInterface.FileFinder(
    #     kappab=1, kappac=1, muR=1, muF=1, Q=1,
    #     expectOneFile=True, loadImmediately=True,
    #     directory=LatestPaths.derivedTheoryFiles_YukawaSummed
    #     )
    # xs_Yukawa = TheoryCommands.Rebin(SM_Yukawa.binBoundaries, SM_Yukawa.crosssection, obs_pth.binning, lastBinIsOverflow=False)
    # normalizationCrossCheck.add_spectrum('Pier', xs_Yukawa, color=2)

    SM_Pier = differentials.theory.theory_utils.FileFinder(
        kappab=1.0, kappac=1.0, muR=1.0, muF=1.0, Q=1.0,
        directory=LatestPaths.theory.yukawa.filedir
        ).get_one()
    SM_Pier_rebinned = differentials.theory.theory_utils.rebin_theory(SM_Pier, [0., 15., 30., 45., 80., 120.])
    normalizationCrossCheck.add_spectrum('Pier', SM_Pier_rebinned.crosssection, color=2)

    ws = LatestPaths.ws.yukawa.nominal.combination
    with differentials.core.openroot(ws) as ws_fp:
        w = ws_fp.Get('w')
        reweightors = []
        for left, right in zip(obs_pth.binning[:-1], obs_pth.binning[1:]):
            name = 'reweightor_ggH_PTH_{0:d}_{1:d}'.format(int(left), int(right))
            reweightor = w.function(name)
            if reweightor == None:
                raise ValueError('Workspace {0} does not have a function called {1}'.format(ws, name))
            reweightors.append(reweightor)
        reweighting_factors = [ rew.getVal() for rew in reweightors ]
        # ws_smxs = [ xspar.getVal() for xspar in differentials.core.read_set(w, 'SMXS', return_names=False) ]
        # widths = [ r-l for l,r in zip(bin_boundaries[:-1], bin_boundaries[1:]) ]
        # ws_xs = [ rew*xs/width for rew, xs, width in zip(reweighting_factors, ws_smxs, widths) ]
        # normalizationCrossCheck.add_spectrum('Pier reweighted', ws_xs, line_style=2, color=2)
        Vitt_reweighted = [ rew*xs for rew, xs in zip(reweighting_factors, SMXSs)]
        normalizationCrossCheck.add_spectrum('Vitt. reweighted', Vitt_reweighted, line_style=2, color=4)

    # # Add SM from the ws (what is actually used)
    # ws = 'out/workspaces_Jan29/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_couplingDependentBR_reweighted.root'
    # with Commands.OpenRootFile(ws) as wsFp:
    #     w = wsFp.Get('w')
    #     reweightors = []
    #     for left, right in zip(obs_pth.binning[:-1], obs_pth.binning[1:]):
    #         name = 'reweightor_ggH_PTH_{0:d}_{1:d}'.format(int(left), int(right))
    #         try:
    #             reweightors.append( w.function(name).getVal() )
    #         except TypeError:
    #             print 'Error: left = {0}, right = {1}'.format(left, right)
    #             print name
    #             raise
    #     SMXS_ws = [ w.var(y).getVal() for y in Commands.ListSet(w, 'SMXS') ]
    # binWidths = [ right-left for left, right in zip(normalizationCrossCheck.bin_boundaries[:-1], normalizationCrossCheck.bin_boundaries[1:]) ]
    # ws_xs = [ reweightor * SMXS / binWidth for reweightor, SMXS, binWidth in zip(reweightors, SMXS_ws, binWidths) ]
    # normalizationCrossCheck.add_spectrum('Pier reweighted', ws_xs, line_style=2)

    normalizationCrossCheck.plot_spectra('Yukawa')

