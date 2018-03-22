import os
import importlib
import logging
import differentials

import scipy
scipyversion = scipy.__version__.split('.')[1]
if int(scipyversion) < 15:
    logging.error('Failed to load scipy.optimize; need to pass the python environment')
    _optimize_loaded = False
else:
    from scipy.optimize import minimize
    _optimize_loaded = True

#____________________________________________________________________

# SM cross section in the public hgg binning scheme, taken from differentialCombination V3
def normalize(l, normalization=1.0):
    s = float(sum(l))
    return [ e/s*normalization for e in l ]
smxs = normalize(
    [ 27.46607434,  29.44961617,  19.6934684,   26.59903008,  10.46543842, 6.31602595, 2.0582112, 0.28218741 ],
    normalization = 55.620335
    )

def load_data(module):
    module = module.replace('/','.').replace('.py','')
    mod = importlib.import_module(module)
    return mod.data

def data_to_hist(data, do_xs=False, title='', x_max=500., color=None):
    if do_xs:
        getstr = 'xs'
    else:
        getstr = 'mu'
    histogram = differentials.plotting.pywrappers.Histogram(
        differentials.plotting.plotting_utils.get_unique_rootname(),
        title,
        data.binning,
        data[getstr]
        )
    histogram.set_err_up(data[getstr + '_up'])
    histogram.set_err_down(data[getstr + '_down'])
    histogram.set_last_bin_is_overflow(method='HARDVALUE', hard_value=500.)
    if not(color is None):
        histogram.color = color
    return histogram

def get_axis(x_min, x_max, n_points):
    dx = (x_max-x_min)/(n_points-1)
    return [ x_min + i*dx for i in xrange(n_points) ]

def get_kappab_kappac_parametrization():
    import LatestPaths
    coupling_variations = differentials.theory.theory_utils.FileFinder(
        muR=1.0, muF=1.0, Q=1.0, directory=LatestPaths.theory.yukawa.filedir
        ).get()
    sm = [ v for v in coupling_variations if v.kappab==1.0 and v.kappac==1.0 ][0]
    coupling_variations.pop(coupling_variations.index(sm))

    parametrization = differentials.parametrization.Parametrization()
    parametrization.parametrize_by_matrix_inversion = True
    parametrization.c1_name = 'kappab'
    parametrization.c2_name = 'kappac'
    parametrization.from_theory_dicts(coupling_variations)
    parametrization.bin_boundaries = sm.binBoundaries
    return parametrization

#____________________________________________________________________


class CouplingFit(object):
    """docstring for CouplingFit"""
    def __init__(self):
        super(CouplingFit, self).__init__()
        self.c1_min = -10.
        self.c1_max = 10.
        self.c1_n_points = 70
        self.c2_min = -35.
        self.c2_max = 35.
        self.c2_n_points = 70
        self.title = ''
        self.color = None
        self.bestfit = None
        
    def setup_kappabkappac(self, data):
        self.chi2 = Chi2CouplingFitter()
        self.chi2.c1_name = 'kappab'
        self.chi2.c2_name = 'kappac'
        self.chi2.set_parametrization_kappabkappac()
        self.chi2.set_data(data)
        self.chi2.build()

    def get_bestfit(self):
        self.bestfit = self.chi2.minimize()

    def get_scan(self):
        if self.bestfit is None: self.get_bestfit()
        self.c1_bin_boundaries = get_axis(self.c1_min, self.c1_max, self.c1_n_points+1)
        self.c2_bin_boundaries = get_axis(self.c2_min, self.c2_max, self.c2_n_points+1)
        self.c1_bin_centers = [ 0.5*(l+r) for l, r in zip(self.c1_bin_boundaries[:-1], self.c1_bin_boundaries[1:]) ]
        self.c2_bin_centers = [ 0.5*(l+r) for l, r in zip(self.c2_bin_boundaries[:-1], self.c2_bin_boundaries[1:]) ]

        scan = [[ 999. for j in xrange(self.c2_n_points) ] for i in xrange(self.c1_n_points)]
        for i_c1, c1 in enumerate(self.c1_bin_centers):
            for i_c2, c2 in enumerate(self.c2_bin_centers):
                scan[i_c1][i_c2] = self.chi2.evaluate([c1, c2]) - self.bestfit.chi2
        self.scan = scan

    def to_hist(self):
        histogram2D = differentials.plotting.pywrappers.Histogram2D(
            differentials.plotting.plotting_utils.get_unique_rootname(), self.title, self.color
            )
        histogram2D.fill_with_matrix(self.scan, self.c1_bin_boundaries, self.c2_bin_boundaries)
        histogram2D.fill_bestfit(self.bestfit.pois[0], self.bestfit.pois[1])
        histogram2D.x_title = self.chi2.c1_name
        histogram2D.y_title = self.chi2.c2_name
        histogram2D.z_title = '#chi^2'
        return histogram2D

    # def to_hist(self):
    #     histogram2D = plotting.pywrappers.Histogram2D(
    #         utils.get_unique_rootname(), getattr(self, 'title', ''), self.color
    #         )
    #     histogram2D.x_title = self.x_title
    #     histogram2D.y_title = self.y_title
    #     histogram2D.z_title = self.z_title
    #     histogram2D.fill_from_entries(self.entries)
    #     if hasattr(self, 'contour_filter_method'):
    #         histogram2D.contour_filter_method = self.contour_filter_method
    #     return histogram2D


class Chi2(object):
    """docstring for Chi2"""
    def __init__(self):
        super(Chi2, self).__init__()
        self.freeze_pois = []
        self.freeze_poi_values = []
        
    def freeze_poi(self, i_poi, value):
        if not i_poi in self.freeze_pois:
            self.freeze_pois.append(i_poi)
            self.freeze_poi_values.append(value)
        else:
            # The poi is already frozen, only the value should be updated
            i_frozen_poi = self.freeze_pois.index(i_poi)
            self.freeze_poi_values[i_frozen_poi] = value

    def unfreeze_poi(self, i_poi):
        if not i_poi in self.freeze_pois:
            logging.warning('poi {0} was not frozen in the first place'.format(i_poi))
            return
        i_frozen_poi = self.freeze_pois.index(i_poi)
        self.freeze_pois.pop(i_frozen_poi)
        self.freeze_poi_values.pop(i_frozen_poi)

    def unfreeze_all(self):
        self.freeze_pois = []
        self.freeze_poi_values = []

    def get_all_pois(self, floating_pois):
        if len(self.freeze_pois) == 0:
            return floating_pois
        all_pois = []
        floating_pois_iter = iter(floating_pois)
        frozen_pois_iter = iter(self.freeze_poi_values)
        for i in xrange(len(self.pois)):
            if i in self.freeze_pois:
                all_pois.append(frozen_pois_iter.next())
            else:
                all_pois.append(floating_pois_iter.next())
        return all_pois

    def minimize(self):
        pois = [ p for i, p in enumerate(self.pois) if not i in self.freeze_pois ] # pick only floating pois
        fit = minimize(self.evaluate, pois, method='Nelder-Mead', tol=1e-6)
        chi2 = fit.fun
        fitted_pois = self.get_all_pois(list(fit.x)) # get fit, and plug back in the frozen pois as well
        return differentials.core.AttrDict(chi2=chi2, pois=fitted_pois)


class Chi2CouplingFitter(Chi2):
    """docstring for Chi2CouplingFitter; hard-coded for maximum 2 couplings"""
    def __init__(self):
        super(Chi2CouplingFitter, self).__init__()
        self.parametrization = None
        self.pois = [ 1., 1. ]

    def build(self):
        self.get_rebinner()
        self.get_smxs_parametrization()

    def get_smxs_parametrization(self):
        self.smxs_parametrization = self.evaluate_parametrization_xs(self.pois) # Should give SM values!

    def get_rebinner(self):
        self.rebinner = differentials.integral.Rebinner(
            bin_boundaries_old = self.parametrization.bin_boundaries,
            bin_boundaries_new = self.binning()
            )

    def set_data(self, data):
        self.data = data
        self.data.delta = [ 0.5*(abs(up)+abs(down)) for up, down in zip(self.data.mu_up, self.data.mu_down) ]

    def binning(self):
        return self.data.binning

    def set_parametrization_kappabkappac(self):
        self.parametrization = get_kappab_kappac_parametrization()
        self.c1_name = self.parametrization.c1_name
        self.c2_name = self.parametrization.c2_name

    def evaluate_parametrization_xs(self, pois):
        xs_per_GeV_theory_binning = self.parametrization.evaluate(*pois)
        xs_per_GeV = self.rebinner.rebin_values(xs_per_GeV_theory_binning)
        return xs_per_GeV

    def evaluate_parametrization(self, pois):
        xs_per_GeV = self.evaluate_parametrization_xs(pois)
        mus = [ xs/xs_at_sm if xs_at_sm != 0. else 0. for xs, xs_at_sm in zip(xs_per_GeV, self.smxs_parametrization) ]
        return mus

    def evaluate(self, pois):
        mus_parametrization = self.evaluate_parametrization(pois)
        chi2 = 0.0
        for mu_data, mu_param, delta_data in zip(self.data.mu, mus_parametrization, self.data.delta):
            chi2 += ((mu_data-mu_param)**2) / (delta_data**2)
        return chi2

    # def minimize(self):
    #     fit = super(Chi2CouplingFitter, self).minimize()
    #     fit['mu'] = self.evaluate_parametrization(fit.pois)
    #     fit['xs'] = self.evaluate_parametrization_xs(fit.pois)
    #     return fit


class Scan(object):
    """docstring for Scan"""
    def __init__(self, xs, ys):
        super(Scan, self).__init__()
        self.xs = xs
        self.ys = ys
        self.unc = None
        self.title = ''
    
    def to_graph(self):
        name = differentials.plotting.plotting_utils.get_unique_rootname()
        graph = differentials.plotting.pywrappers.Graph(name, self.title, self.xs, self.ys)
        if hasattr(self, 'unc'): graph.unc = self.unc
        if hasattr(self, 'color'): graph.color = self.color
        if hasattr(self, 'draw_style'): graph.draw_style = self.draw_style
        return graph

class PtCombination(object):
    """docstring for PtCombination"""
    def __init__(self, spectra=None):
        super(PtCombination, self).__init__()
        self.chi2 = Chi2PtCombination()
        self.uncertaintycalculator = differentials.uncertaintycalculator.UncertaintyCalculator()
        self.uncertaintycalculator.cutoff = 1.0 # Default is 0.5, which is tailored to combine output

        self.scans = []
        self.title = 'ptcombination'
        self.hard_x_max = 500. # For plotting only
        self.outdir = 'fermilabcode'

        if not(spectra is None):
            self.spectra = spectra
            for spectrum in self.spectra:
                self.chi2.add_spectrum(spectrum)
            self.chi2.build()
        self.bestfit_done = False

    def get_bestfit(self):
        if self.bestfit_done:
            # Do nothing
            pass
        else:
            self.bestfit = self.chi2.minimize()
            self.bestfit_done = True

    def get_scan(self, i):
        self.get_bestfit() # Make sure the bestfit is there
        poi_axis = get_axis(-1.5, 3.0, 70)
        chi2_axis = []
        for poi in poi_axis:
            self.chi2.freeze_poi(i, poi)
            fit = self.chi2.minimize()
            chi2_axis.append(fit.chi2 - self.bestfit.chi2)
        self.chi2.unfreeze_all()

        # Insert the bestfit values into the scan as well (at the right place)
        for i_poi_axis in xrange(len(poi_axis)-1):
            if poi_axis[i_poi_axis] <= self.bestfit.pois[i] and self.bestfit.pois[i] <= poi_axis[i_poi_axis+1]:
                poi_axis.insert(i_poi_axis+1, self.bestfit.pois[i])
                chi2_axis.insert(i_poi_axis+1, 0.0)
                break
        else:
            logging.error('Could not plug in the bestfit in the chi2-scan!')

        # for poi, chi2 in zip(poi_axis, chi2_axis):
        #     print 'poi = {0:+9.3f}; chi2 = {1:+9.3f}'.format(poi, chi2)

        unc = self.uncertaintycalculator.create_uncertainties(poi_axis, chi2_axis)
        if unc.is_hopeless:
            logging.error(
                'Uncertainty determination for poi {0} failed; increase the range'.format(i)
                )
            raise RuntimeError

        scan = Scan(poi_axis, chi2_axis)
        scan.unc = unc
        scan.title = 'poi{0}'.format(i)
        return scan

    def get_scans(self):
        self.scans = []
        for i in xrange(self.chi2.n_bins):
            scan = self.get_scan(i)
            self.scans.append(scan)

    def plot_scans(self, plotname=None):
        if len(self.scans) == 0:
            logging.warning('Scans were not yet calculated; scanning now...')
            self.get_scans()
        if plotname is None:
            plotname = 'scans_naivept'
        plot = differentials.plotting.plots.MultiScanPlot(plotname)
        plot.manual_graphs = [ scan.to_graph() for scan in self.scans ]
        plot.x_min = -2.
        plot.x_max = 5.
        plot.draw()
        plot.wrapup()

    def binning(self):
        return self.chi2.binning

    def to_hist(self):
        histogram = differentials.plotting.pywrappers.Histogram(
            differentials.plotting.plotting_utils.get_unique_rootname(),
            self.title,
            self.binning(),
            [ s.unc.x_min for s in self.scans ]
            )
        histogram.set_err_up([ s.unc.right_error for s in self.scans ])
        histogram.set_err_down([ s.unc.left_error for s in self.scans ])
        histogram.set_last_bin_is_overflow(method='HARDVALUE', hard_value=self.hard_x_max)
        if hasattr(self, 'color'): histogram.color = self.color
        return histogram

    def get_histogram_bestfit_hzz_binning(self):
        hzz_chi2 = self.chi2.chi2_spectra[1]
        hzz_binning = hzz_chi2.binning
        mu_hzz_binning = hzz_chi2.binmerger.evaluate(self.bestfit.pois)

        histogram = differentials.plotting.pywrappers.Histogram(
            differentials.plotting.plotting_utils.get_unique_rootname(),
            self.title + '_hzzBins',
            hzz_binning,
            mu_hzz_binning,
            )
        histogram.set_last_bin_is_overflow(method='HARDVALUE', hard_value=self.hard_x_max)
        if hasattr(self, 'color'): histogram.color = self.color
        return histogram

    def dump(self, name):
        out_file = os.path.join(self.outdir, name + '_{0}.py'.format(differentials.core.datestr()))
        contents = [
            '# ' + differentials.core.gittag(),
            'from differentials.core import AttrDict',
            'data = AttrDict(',
            '    binning = {0},'.format(self.binning()),
            '    mu      = {0},'.format(self.bestfit.pois),
            '    mu_up   = {0},'.format([ s.unc.right_error for s in self.scans ]),
            '    mu_down = {0},'.format([ s.unc.left_error for s in self.scans ]),
            '    )'
            ]
        contents = '\n'.join(contents)
        logging.info(
            'Writing following contents to {0}:\n{1}'
            .format(out_file, contents)
            )
        if not differentials.core.is_testmode():
            with open(out_file, 'w') as out_fp:
                out_fp.write(contents)


class Chi2PtCombination(Chi2):
    """docstring for Chi2PtCombination"""
    def __init__(self):
        super(Chi2PtCombination, self).__init__()
        self.spectra = []
        self.chi2_spectra = []
        self.pois = [ 1. for i in xrange(self.n_bins) ]
       
    def add_spectrum(self, spectrum):
        self.spectra.append(spectrum)

    def build(self):
        # Pick finest binning (i.e. most bins)
        finest = sorted(self.spectra, key=lambda s: -len(s.binning))[0]

        self.binning = finest.binning
        self.n_bins = len(self.binning)-1

        self.chi2_spectra = [] # Make sure it's empty, so .build() can be re-called
        for spectrum in self.spectra:
            chi2_spectrum = Chi2Spectrum(spectrum)
            chi2_spectrum.set_finest(finest)
            self.chi2_spectra.append(chi2_spectrum)

    def evaluate(self, pois):
        # logging.trace('Evaluating combination; pois = {0}'.format(pois))
        pois = self.get_all_pois(pois) # Add frozen pois as well; does nothing if there is nothing frozen
        chi2 = 0.0
        for chi2_spectrum in self.chi2_spectra:
            chi2_this_spectrum = chi2_spectrum.evaluate(pois)
            # logging.trace('Adding evaluation of {0}: {1}'.format(chi2_spectrum.spectrum.name, chi2_this_spectrum))
            chi2 += chi2_this_spectrum
        return chi2


class Chi2Spectrum(object):
    """docstring for Chi2Spectrum"""
    def __init__(self, spectrum):
        super(Chi2Spectrum, self).__init__()
        self.spectrum = spectrum

        self.binning = self.spectrum.binning
        self.n_bins = len(self.binning)-1

        self.mus = self.spectrum.mu
        self.deltas = [ 0.5*(abs(l)+abs(r)) for l, r in zip(self.spectrum.mu_down, self.spectrum.mu_up) ]

        self.needs_rebinning = False

    def set_finest(self, finest):
        if finest.binning == self.spectrum.binning:
            self.needs_rebinning = False
            return
        self.needs_rebinning = True
        self.finest = finest
        self.binmerger = BinMerger(self.finest, self.spectrum)
        self.binmerger.build()

    def evaluate(self, pois):
        if self.needs_rebinning:
            pois = self.binmerger.evaluate(pois)

        chi2 = 0.0
        for poi, mu, delta in zip(pois, self.mus, self.deltas):
            chi2 += ((poi-mu)**2) / (delta**2)
        return chi2


class BinMerger(object):
    """docstring for BinMerger"""
    def __init__(self, fine, coarse):
        super(BinMerger, self).__init__()
        self.fine = fine
        self.n_fine = len(self.fine.binning)-1
        self.coarse = coarse
        self.n_coarse = len(self.coarse.binning)-1
        self.merge_map = [ [] for i in xrange(self.n_coarse) ]
        self.bin_expressions = []

    def build(self):
        self.get_merge_map()
        self.fill_expressions()

    def get_merge_map(self):
        for i_fine in xrange(self.n_fine):
            left_fine = self.fine.binning[i_fine]
            right_fine = self.fine.binning[i_fine+1]
            i_coarse = self.get_corresponding_coarse_bin(left_fine, right_fine)
            self.merge_map[i_coarse].append(i_fine)

    def get_corresponding_coarse_bin(self, left_fine, right_fine):
        for i_coarse in xrange(self.n_coarse):
            left_coarse = self.coarse.binning[i_coarse]
            right_coarse = self.coarse.binning[i_coarse+1]
            if left_coarse <= left_fine and right_fine <= right_coarse:
                break
        else:
            raise RuntimeError('Could not find a coarse bin')
        return i_coarse

    def fill_expressions(self):

        # fine_widths = [ self.fine.binning[i+1]-self.fine.binning[i] for i in xrange(len(self.fine.binning)) ]
        # # Last bin is overflow; give it the same weight as the second-to-last bin
        # fine_widths[-1] = fine_widths[-2]

        for i_coarse in xrange(self.n_coarse):
            fine_bins = self.merge_map[i_coarse]
            # If unmerged, dont do anything
            if len(fine_bins) == 1:
                i_fine = fine_bins[0]
                expr = lambda pois, i_fine=i_fine: pois[i_fine]
            else:
                # WRONG: Not multiplied by bin width
                # xs_fine = []
                # for i_fine in fine_bins:
                #     xs_fine.append(self.fine.xs[i_fine])
                # xs_total = sum(xs_fine)
                # expr = lambda pois, xs_fine=xs_fine, xs_total=xs_total, fine_bins=fine_bins: sum(
                #     [ xs_fine[i]/xs_total * pois[i_fine] for i, i_fine in enumerate(fine_bins) ]
                #     )
                
                # WRONG: Using the fine xs instead of SM
                # xs_fine = [
                #     self.fine.xs[i_fine] * (self.fine.binning[i_fine+1]-self.fine.binning[i_fine])
                #     for i_fine in fine_bins
                #     ]
                # xs_fine_total = sum(xs_fine)
                # coefficients = [ xs/xs_fine_total for xs in xs_fine ]
                # expr = lambda pois, coefficients=coefficients, fine_bins=fine_bins: sum(
                #     [ coeff * pois[i_fine] for coeff, i_fine in zip(coefficients, fine_bins) ]
                #     )

                # Using already integrated smxs; no need for multiplying by bin width
                xs_fine = [ smxs[i_fine] for i_fine in fine_bins ]
                xs_fine_total = sum(xs_fine)
                coefficients = [ xs/xs_fine_total for xs in xs_fine ]
                expr = lambda pois, coefficients=coefficients, fine_bins=fine_bins: sum(
                    [ coeff * pois[i_fine] for coeff, i_fine in zip(coefficients, fine_bins) ]
                    )

                print 'Merged bins {0}'.format(fine_bins)
                print '  xs: {0}'.format(xs_fine)
                print '  coefficients: {0}'.format(coefficients)
                print '  sum(coefficients): {0}'.format(sum(coefficients))
                mock_pois = [ 1.0 for i in xrange(len(self.fine.xs)) ]
                mock_pois[fine_bins[1]] = 2.0
                print '  testing mock pois: {0}'.format(mock_pois)
                print '    resulting poi: {0}'.format(expr(mock_pois))

            self.bin_expressions.append(expr)

    def evaluate(self, pois):
        merged_pois = [ expr(pois) for expr in self.bin_expressions ]
        return merged_pois



