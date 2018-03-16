from OptionHandler import flag_as_option, flag_as_parser_options

import os, logging, copy, sys, re
from array import array

import ROOT

import LatestPaths, LatestBinning
import differentials
import differentialutils

from differentials.plotting.canvas import c


#____________________________________________________________________

class PDFDrawer(object):
    """docstring for PDFDrawer"""
    def __init__(self, ws):
        super(PDFDrawer, self).__init__()
        self.ws = ws
        self.w = differentials.core.get_ws(self.ws)

        self.MH = self.w.var('MH')
        self.MH.setVal(125.)

        self.mH_hgg   = self.w.var('CMS_hgg_mass')
        ROOT.SetOwnership( self.mH_hgg, False )
        self.default_mH_hgg_range = [ self.mH_hgg.getMin(), self.mH_hgg.getMax() ]
        self.default_mH_hgg_nBins = self.mH_hgg.getBins()

        self.all_pdfs = []
        all_pdfs_arglist = ROOT.RooArgList(self.w.allPdfs())
        for i_pdf in xrange(all_pdfs_arglist.getSize()):
            self.all_pdfs.append( all_pdfs_arglist[i_pdf].GetName() )
        self.all_pdfs.sort()

        self.bin_width = 0.25

        self.hgg_cat = self.w.cat('CMS_channel')
        ROOT.SetOwnership(self.hgg_cat, False )


    def set_testable_yield_parameter(self, name):
        self.mu = self.w.var(name)
        if self.mu == None:
            raise ValueError
        self.mu.setRange(-100., 100.)

    def get_histogram(self, pdf_name):
        mH = self.mH_hgg
        default_mH_range = self.default_mH_hgg_range

        # Get bkg and sig pdfs
        mH.setRange( *default_mH_range )
        mH.setBins( int((default_mH_range[1]-default_mH_range[0])/self.bin_width) )

        pdf = self.w.pdf(pdf_name)
        if pdf == None:
            raise ValueError('Pdf {0} does not exist in {1}'.format(pdf_name, self.ws))

        H = pdf.createHistogram('H_'+pdf_name, mH)
        # Hbkg.Scale(self.bin_width) # default is /GeV
        ROOT.SetOwnership(H, False)

        return H

    def set_category(self, some_string):
        i_cat = int(re.search(r'SigmaMpTTag_(\d)', some_string).group(1))
        # cat_str = 'SigmaMpTTag_{0}_recoPt_600p0to10000p0_13TeV'.format(i_cat)
        cat_str = 'recoPt_600p0_10000p0_SigmaMpTTag_0_13TeV'

        # self.hgg_cat.setIndex(i_cat)
        self.hgg_cat.setLabel(cat_str)

    def get_data_histogram(self, pdf_name):
        self.set_category(pdf_name)
        dataset = self.w.data('data_obs')
        ROOT.SetOwnership( dataset, False )
        Hdata = self.roodataset_to_hist(dataset, self.mH_hgg, self.hgg_cat)
        return Hdata

    def roodataset_to_hist(self, dataset, xVar, category):
        ctemp = ROOT.TCanvas( 'ctemp', 'ctemp', 1000, 800 )
        ctemp.cd()
        ctemp.Clear()

        frame = xVar.frame()

        reduceStr = '{0}=={1}'.format( category.GetName(), category.getIndex() )

        dataset_reduced = dataset.reduce( reduceStr )
        dataset_reduced.plotOn( frame )
        frame.Draw()

        l = ctemp.GetListOfPrimitives()
        for i in xrange(l.GetEntries()):
            if isinstance( l.At(i), ROOT.RooHist ):
                H = l.At(i)
                break
        else:
            raise ValueError( 'ERROR: did not find a histogram', throwException=True )

        Hcopy = ROOT.RooHist( H )
        ROOT.SetOwnership( Hcopy, False )

        Hcopy.SetName( differentials.plotting.plotting_utils.get_unique_rootname() )

        # ctemp.SaveAs( 'plots_{0}_onetimeplots/roodatasetplottest.pdf'.format(datestr) )
        del ctemp
        del frame

        c.cd()

        return Hcopy

        

@flag_as_option
def hgg_pdf_xcheck(args):

    ws_smH = 'out/workspaces_Mar14/ws_pth_smH_hgg.root'
    ws_ggH = 'out/workspaces_Mar02/ws_pth_ggH_hgg.root'

    yp = {
        ws_smH : 'r_smH_PTH_GT600',
        ws_ggH : 'r_ggH_PTH_GT600'
        }

    for ws in [ ws_smH, ws_ggH ]:

        pdfdrawer = PDFDrawer(ws)
        pdfdrawer.set_testable_yield_parameter(yp[ws])

        c.Clear()
        c.set_margins()

        base = differentials.plotting.pywrappers.Base(
                x_min = 100.,
                x_max = 180.,
                y_min = 0.,
                # y_max = 0.20,
                y_max = 1.20,
                x_title = 'mH',
                y_title = '',
                )
        base.Draw()

        pdf_name = 'pdf_binrecoPt_600p0_10000p0_SigmaMpTTag_0_13TeV'

        for mu_val in [ 0.5, 1.0, 2.0, -1.0, -2.0, -3.0, -5.0, -10.0 ]:
            pdfdrawer.mu.setVal(mu_val)
            # H = pdfdrawer.get_histogram('pdf_binrecoPt_600p0_10000p0_SigmaMpTTag_0_13TeV_obsOnly')
            H = pdfdrawer.get_histogram(pdf_name)
            H.Draw('HISTSAME')

        H_data = pdfdrawer.get_data_histogram(pdf_name)
        H_data.Draw('HISTSAME')

        outname = 'hgg_pdf_xcheck_' + os.path.basename(ws).replace('.root','')
        c.save(outname)



#____________________________________________________________________
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

    def get_extrema(self):
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

        return x_min, y_min, x_max, y_max


    def plot_spectra(self, tag=None):
        outname = 'normalizationCrosscheck_reimpl'
        if not(tag is None):
            outname += '_' + tag
        plot = differentials.plotting.plots.BottomPanelPlot(outname)

        x_min, y_min, x_max, y_max = self.get_extrema()
        plot.top_x_min = x_min
        plot.top_y_min = y_min
        plot.top_x_max = x_max
        plot.top_y_max = y_max
        plot.bottom_x_max = x_max
        plot.bottom_y_max = y_max

        for spectrum in self.spectra:
            H = self.spectrum_to_hist(spectrum)
            # H.Draw('HISTSAME')
            # leg.AddEntry( H.GetName(), spectrum.title, 'l' )
            plot.add_top(H, 'HISTSAME')

        plot.draw()
        plot.wrapup()



    def plot_spectra_old(self, tag=None):
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

    # Add raw from Pier, only rebinned
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

    normalizationCrossCheck.plot_spectra('Yukawa')

