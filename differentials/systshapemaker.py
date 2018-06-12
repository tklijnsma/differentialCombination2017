import os, logging, traceback, copy
import core, scans

from math import sqrt

class SystShapeMaker(object):
    """docstring for SystShapeMaker"""
    def __init__(self):
        super(SystShapeMaker, self).__init__()
        self.success = False
        self.set_smxs = False

    def set_sm(self, smxs):
        self.smxs = smxs
        self.set_smxs = True

    def scandir_to_spectrum(self, scandir):
        if not self.set_smxs:
            raise RuntimeError(
                'Required to set the SM first, use set_sm(list of cross sections)'
                )
        if not os.path.isdir(scandir):
            raise ValueError(
                'Path \'{0}\' does not exist'.format(scandir)
                )
        spectrum = scans.DifferentialSpectrum(
            'name_{0}'.format(core.__uniqueid__().next()), scandir
            )
        spectrum.set_sm(self.smxs)
        spectrum.read()
        return spectrum
       
    def get_systonly_histogram(self, combination, combination_statonly):
        # Convert to Spectrum if necessary
        if isinstance(combination, basestring):
            combination = self.scandir_to_spectrum(combination)
        if isinstance(combination_statonly, basestring):
            combination_statonly = self.scandir_to_spectrum(combination_statonly)

        try:
            logging.info('Getting syst-only shape for combination')
            systonly_histogram, systonly_histogram_xs = self.get_systonly_histogram_notry(combination, combination_statonly)

            sheet = combination.style().copy(
                plot_priority = combination.style().plot_priority - 1,
                marker_size=0,
                line_width = 0,
                line_color = 0,
                fill_style = 1001,
                fill_color = core.safe_colors['lightblue']
                )
            systonly_histogram.add_stylesheet(sheet)
            systonly_histogram_xs.add_stylesheet(sheet)

            # systonly_histogram.draw_method = 'repr_uncertainties_narrow_filled_area'
            # systonly_histogram_xs.draw_method = 'repr_uncertainties_narrow_filled_area'
            systonly_histogram.draw_method = 'repr_narrow_bar_onlyfill_legend'
            systonly_histogram_xs.draw_method = 'repr_narrow_bar_onlyfill_legend'

            # systonly_histogram_xs.move_to_bottom_of_legend = True
            self.success = True
            return systonly_histogram, systonly_histogram_xs
        except Exception as exc:
            self.success = False
            logging.error('Getting syst-only shape FAILED:')
            print traceback.format_exc()
            print exc
            return False, False

    def get_systonly_histogram_notry(self, combination, combination_statonly):
        # Calculate the syst-only errors
        combination_histogram = combination.to_hist()
        statonly_histogram    = combination_statonly.to_hist()

        systonly_histogram = copy.deepcopy(combination_histogram)
        systonly_histogram.title = 'Syst. unc.'
        systonly_histogram.draw_method = 'repr_uncertainties_fully_filled_area'
        systonly_histogram_xs = combination.to_hist_xs()
        systonly_histogram_xs.title = 'Syst. unc.'
        systonly_histogram_xs.draw_method = 'repr_uncertainties_fully_filled_area'

        syst_err_up = []
        syst_err_down = []
        syst_err_up_xs = []
        syst_err_down_xs = []
        for i in xrange(combination_histogram.n_bins):
            up_tot = abs(combination_histogram.errs_up[i])
            down_tot = abs(combination_histogram.errs_down[i])
            symm_tot = 0.5*(up_tot+down_tot)
            if symm_tot == 0.0:
                logging.error(
                    'Found symm_tot == 0.0 (probably because of faulty uncertainty determination). '
                    'Will now use 999 instead.'
                    )
                symm_tot = 999.

            symm_stat = 0.5*(abs(statonly_histogram.errs_up[i]) + abs(statonly_histogram.errs_down[i]))
            sm_crosssection = combination.smxs[i]

            if symm_tot>symm_stat:
                symm_syst = sqrt(symm_tot**2-symm_stat**2)
            else:
                logging.warning(
                    'For bin {0}, symm_tot={1} > symm_stat={2}. '
                    'Taking sqrt(symm_stat**2-symm_tot**2) instead'
                    .format(i, symm_tot, symm_stat)
                    )
                symm_syst = sqrt(symm_stat**2-symm_tot**2)
            syst_err_up.append(up_tot * symm_syst/symm_tot)
            syst_err_down.append(down_tot * symm_syst/symm_tot)
            syst_err_up_xs.append(up_tot * symm_syst/symm_tot * sm_crosssection)
            syst_err_down_xs.append(down_tot * symm_syst/symm_tot * sm_crosssection)
        systonly_histogram.set_err_up(syst_err_up)
        systonly_histogram.set_err_down(syst_err_down)
        systonly_histogram_xs.set_err_up(syst_err_up_xs)
        systonly_histogram_xs.set_err_down(syst_err_down_xs)
        return systonly_histogram, systonly_histogram_xs
