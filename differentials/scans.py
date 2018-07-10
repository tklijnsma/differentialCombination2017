import os.path
import glob, re, copy, math, sys

import logging
import core
import plotting
from plotting.canvas import c
import plotting.plotting_utils as utils
from uncertaintycalculator import UncertaintyCalculator
from spline2d import Spline2DFactory
from onedimscanfilter import OneDimScanFilter

from collections import namedtuple
from array import array

import ROOT


def glob_rootfiles(d):
    if not d.endswith('/'): d += '/'
    return glob.glob(d + '*.root')


class DifferentialSpectrum(object):
    """Essentially a collection of Scan instances, 1 per POI that was scanned"""
    standard_titles = {
        'hgg' : 'H#rightarrow#gamma#gamma',
        'hzz' : 'H#rightarrowZZ',
        'combination' : 'Combination',
        'hbb' : 'H#rightarrowbb',
        'combWithHbb' : 'Combination',
        }

    default_style = plotting.pywrappers.StyleSheet()

    def __init__(self, name, scandirs=None, POIs=None):
        self.name = name

        self.scandirs = []
        if not(scandirs is None):
            if isinstance(scandirs, basestring): scandirs = [scandirs]
            self.scandirs = scandirs

        self.scans = []

        if POIs is None:
            self.auto_determine_POIs = True
            self.POIs = []
        else:
            self.auto_determine_POIs = False
            self.POIs = POIs

        self.title = self.standard_titles.get(name, name)
        self.smxs = []
        self.smxs_set = False
        self.color = None
        self.draw_method = 'repr_horizontal_bar_and_narrow_fill'
        self.hard_x_max = None

        self.get_POIs()
        self._is_read = False

        self.scans_x_min = -10.0
        self.scans_x_max = 10.0

        self.has_x_max = False
        self.stylesheets = [ DifferentialSpectrum.default_style.copy() ]

        self.root_files_for_POI = {}

    def add_POI(self, POI, root_files):
        self.POIs.append(POI)
        if isinstance(root_files, basestring): root_files = [root_files]
        self.root_files_for_POI[POI] = root_files

    def style(self):
        return self.stylesheets[-1]

    def add_stylesheet(self, sheet):
        self.stylesheets.append(sheet)

    def set_sm(self, smxs):
        self.smxs = smxs
        self.smxs_set = True

    def add_scandir(self, scandir):
        self.scandirs.append(scandir)
        self.get_POIs()

    def get_POIs_from_scandir(self):
        root_files = []
        for scandir in self.scandirs:
            root_files.extend(glob_rootfiles(scandir))
        POIs = []
        for root_file in root_files:
            match = re.search(r'bPOI_(\w+)_ePOI', os.path.basename(root_file))
            if not match:
                continue
            POIs.append(match.group(1))
        POIs = list(set(POIs))
        POIs.sort(key=core.range_sorter)
        logging.debug(
            'Found the following set of POIs from {0}: {1}'
            .format(', '.join(self.scandirs), ', '.join(POIs))
            )
        return POIs

    def get_POIs_from_datacard(self, datacard=None):
        if datacard is None:
            datacard = self.datacard
        POIs = core.list_POIs(datacard)
        POIs.sort(key=core.range_sorter)
        return POIs

    def get_POIs(self):
        if self.auto_determine_POIs:
            self.POIs = self.get_POIs_from_scandir()
        # else:
        #     self.POIs = self.get_POIs_from_datacard()

    def read(self):
        if self._is_read:
            raise RuntimeError('Scan {0} is already read'.format(self.name))
        else:
            if len(self.POIs)==0:
                self.get_POIs()
            for POI in self.POIs:
                scan = Scan(
                    x_variable=POI,
                    y_variable='deltaNLL',
                    globpat=POI
                    )
                scan.scandirs.extend(self.scandirs)
                scan.root_files.extend(self.root_files_for_POI.get(POI, []))
                scan.read()
                scan.create_uncertainties()
                self.scans.append(scan)
            self._is_read = True

    def plot_scans(self, plotname=None):
        if plotname is None:
            plotname = 'scans_{0}'.format(self.name)
        plot = plotting.plots.MultiScanPlot(plotname)
        plot.scans = self.scans
        plot.x_min = self.scans_x_min
        plot.x_max = self.scans_x_max
        plot.y_min = self.scans_y_min
        plot.y_max = self.scans_y_max
        plot.draw()
        plot.wrapup()

    def binning(self):
        binning = core.binning_from_POIs(self.POIs)
        if self.has_x_max:
            if self.x_max < binning[-2]:
                raise ValueError(
                    'The histogram is given an x_max of {0}, but overwriting the last'
                    'bin boundary would break the order of the boundaries: {1}'
                    .format(self.x_max, binning)
                    )
            binning[-1] = self.x_max
        return binning

    def give_x_max(self, x_max):
        self.has_x_max = True
        self.x_max = x_max

    def drop_bins_up_to_value(self, value):
        _POIs_before = self.POIs[:]
        while self.binning()[-1] > value:
            logging.debug('Bin bound {0} > {1}; dropping POI {2}'.format(self.binning()[-1], value, self.POIs[-1]))
            self.POIs = self.POIs[:-1]
        logging.info(
            'Dropped POIs after {0}'
            '\nOld POIs: {1}'
            '\nNew POIs: {2}'
            .format(value, _POIs_before, self.POIs)
            )

    def drop_first_bin(self):
        dropped_POI = self.POIs.pop(0)
        logging.info('Dropped POI {0}'.format(dropped_POI))
        if self._is_read:
            self.scans.pop(0)
        if self.smxs_set:
            self.smxs.pop(0)

    def last_bin_is_overflow(self):
        if self.binning()[-1] == 10000.:
            return True
        return core.last_bin_is_overflow(self.POIs)

    def first_bin_is_underflow(self):
        return core.first_bin_is_underflow(self.POIs)

    def process_overflow_and_underflow(self, histogram):
        # if self.last_bin_is_overflow():
        #     if not(self.hard_x_max is None):
        #         histogram.set_last_bin_is_overflow(method='HARDVALUE', hard_value=self.hard_x_max)
        #     else:
        #         histogram.set_last_bin_is_overflow()
        if self.first_bin_is_underflow():
            histogram.drop_first_bin()

    def to_hist(self):
        histogram = plotting.pywrappers.Histogram(
            utils.get_unique_rootname(),
            self.title,
            self.binning(),
            [ s.unc.x_min for s in self.scans ]
            )
        histogram.set_err_up([ s.unc.right_error for s in self.scans ])
        histogram.set_err_down([ s.unc.left_error for s in self.scans ])
        self.process_overflow_and_underflow(histogram)
        # if not(self.color is None): histogram.color = self.color
        histogram.add_stylesheet(self.style())
        return histogram

    def to_hist_xs(self):
        histogram = plotting.pywrappers.Histogram(
            utils.get_unique_rootname(),
            self.title,
            self.binning(),
            [ s.unc.x_min * xs for s, xs in zip(self.scans, self.smxs) ]
            )
        histogram.set_err_up([ s.unc.right_error * xs for s, xs in zip(self.scans, self.smxs) ])
        histogram.set_err_down([ s.unc.left_error * xs for s, xs in zip(self.scans, self.smxs) ])
        self.process_overflow_and_underflow(histogram)
        if not(self.color is None): histogram.color = self.color
        histogram.add_stylesheet(self.style())
        return histogram

    def to_hist_sm(self):
        histogram = plotting.pywrappers.Histogram(
            utils.get_unique_rootname(),
            'SM',
            self.binning(),
            [1.0 for i in self.smxs]
            )
        # histogram.set_err_up([ s.unc.right_error * xs for s, xs in zip(self.scans, self.smxs) ])
        # histogram.set_err_down([ s.unc.left_error * xs for s, xs in zip(self.scans, self.smxs) ])
        self.process_overflow_and_underflow(histogram)
        histogram.add_stylesheet(self.style().copy(
            color=16,
            plot_priority = self.style().plot_priority - 5
            ))
        return histogram

    def to_hist_smxs(self):
        histogram = plotting.pywrappers.Histogram(
            utils.get_unique_rootname(),
            core.standard_titles['SM_Vittorio'],
            self.binning(),
            self.smxs
            )
        # histogram.set_err_up([ s.unc.right_error * xs for s, xs in zip(self.scans, self.smxs) ])
        # histogram.set_err_down([ s.unc.left_error * xs for s, xs in zip(self.scans, self.smxs) ])
        self.process_overflow_and_underflow(histogram)
        histogram.add_stylesheet(self.style().copy(
            color=16,
            plot_priority = self.style().plot_priority - 5
            ))
        return histogram


    def plot_single_spectrum(self):
        histogram = self.to_hist()
        histogram.color = 2

        c.Clear()
        c.set_margins()
        base = utils.get_plot_base(
            x_min=histogram.x_min(), x_max=histogram.x_max(),
            y_min=histogram.y_min(), y_max=histogram.y_max(),
            x_title='Observable', y_title='#mu'
            )
        base.Draw('P')

        for obj, draw_str in histogram.repr_uncertainties_narrow_filled_area():
            obj.Draw(draw_str)
        for obj, draw_str in histogram.repr_horizontal_bars():
            obj.Draw(draw_str)

        outname = 'spectrum_' + os.path.basename(self.scandir).replace('/','')
        c.save(outname)



class ScanPrimitive(object):

    standard_titles = copy.copy(core.standard_titles)
    legacy_standard_titles = {
        'hgg' : 'H#rightarrow#gamma#gamma',
        'hzz' : 'H#rightarrowZZ',
        'combination' : 'Combination',
        'hbb' : 'H#rightarrowbb',
        # 'combWithHbb' : 'Comb. with H#rightarrowbb',
        # 
        'kappac' : '#kappa_{c}',
        'kappab' : '#kappa_{b}',
        'kappat' : '#kappa_{t}',
        'ct' : '#kappa_{t}',
        'cg' : 'c_{g}',
        }
    for key in legacy_standard_titles:
        if not key in standard_titles:
            standard_titles[key] = legacy_standard_titles[key]

    tree_name = 'limit'
    filter_negatives = True
    deltaNLL_threshold = -0.01

    def __init__(self):
        self.scandirs = []
        self.root_files = []
        self.save_all_variables = False
        self.entries = []
        self.globpat = '*'

        self.read_one_scandir = True


    def collect_root_files(self):
        root_files = copy.copy(self.root_files)

        if len(self.scandirs)>0:
            for scandir in self.scandirs[::-1]:
                if not scandir.endswith('/'): scandir += '/'
                root_files_this_scandir = glob.glob(scandir + self.globpat + '.root')
                n_found = len(root_files_this_scandir)
                logging.debug('Found {0} root files in {1}'.format(n_found, scandir))
                root_files.extend(root_files_this_scandir)
                if self.read_one_scandir and len(root_files_this_scandir) > 0:
                    logging.warning('Using only root files found in {0} (ignoring others for globpat {1})'.format(scandir, self.globpat))
                    self.scandirs = [scandir]
                    break
            logging.debug('Found {0} root files in {1}'.format(len(root_files), ', '.join(self.scandirs)))

        if len(root_files) == 0:
            raise RuntimeError(
                'Attemped to retrieve scan for x:{0} y:{1}, '
                'but no .root files were found. Passed list of dirs to look in:\n'
                .format(self.x_variable, self.y_variable)
                + '\n'.join(self.scandirs)
                )

        logging.trace('List of root files:\n' + '\n'.join(root_files))
        return root_files

    def get_list_of_variables_in_tree(self, root_files, accept_pat='*'):
        for root_file in root_files:
            with core.openroot(root_file) as root_fp:
                if not root_fp.GetListOfKeys().Contains(self.tree_name):
                    found_tree = False
                else:
                    variables = []
                    root_array = root_fp.Get(self.tree_name).GetListOfBranches()
                    for i_var in xrange(root_array.GetEntries()):
                        var_name = root_array[i_var].GetTitle()
                        if accept_pat == '*':
                            variables.append( var_name.split('/',1)[0] )
                        else:
                            if re.search( variablePattern, var_name ):
                                variables.append( var_name.split('/',1)[0] )
                    found_tree = True
            if found_tree: break
        else:
            if len(root_files) > 10:
                print_root_files = root_files[:3] + ['...'] + root_files[-4:]
            else:
                print_root_files = root_files
            raise RuntimeError(
                'Not a single root file had a tree called {0}. '
                'List of root files that were checked:\n'
                .format(self.tree_name)
                + '\n'.join(root_files)
                )
        return variables                            

    def read_chain(self, root_files, variables, filter_x=False, return_chain=False):
        if len(root_files) == 0:
            raise RuntimeError(
                'No root files were passed to read_chain; scandirs = {0}'.format(self.scandirs)
                )

        chain = ROOT.TChain(self.tree_name)
        for root_file in root_files:
            chain.Add(root_file)

        # Entry = namedtuple('Entry', variables + ['x', 'y'] + (['z'] if hasattr(self, 'z_variable') else []))

        entries = []
        for event in chain:
            if filter_x and getattr(event, self.x_variable) in [e.x for e in entries]: continue
            # entry_dict = {}
            entry_dict = core.AttrDict()
            for var_name in variables:
                entry_dict[var_name] = getattr(event, var_name)
                if var_name == self.x_variable:
                    entry_dict['x'] = entry_dict[var_name]
                if var_name == self.y_variable:
                    entry_dict['y'] = entry_dict[var_name]
                if hasattr(self, 'z_variable') and var_name == self.z_variable:
                    entry_dict['z'] = entry_dict[var_name]

            if self.y_variable == 'deltaNLL' and entry_dict['y'] == 9990.0 or entry_dict['y']>1e9:
                logging.trace(
                    'Found deltaNLL=={3}, which indicates a misfit ({0}={1}, {2}={3}). '
                    'Skipping point.'
                    .format(self.x_variable, entry_dict['x'], self.y_variable, entry_dict['y'])
                    )
                continue
            # entries.append(Entry(**entry_dict))
            entries.append(entry_dict)

        if len(entries) == 0:
            raise RuntimeError(
                'No entries were found; looked in root files such as {0}'.format(root_files[0])
                )

        entries.sort(key=lambda entry: (entry.x, entry.y))

        if return_chain:
            return entries, chain
        else:
            return entries


    def filter_entries(self, inplace=True):
        passed_entries = []
        for entry in self.entries:
            if entry.deltaNLL < self.deltaNLL_threshold:
                # POIs = [ k for k in entry._fields if k.startswith('r_') ]
                POIs = [ k for k in entry.keys() if k.startswith('r_') ]
                if len(POIs) == 0:
                    POI = 'x'
                else:
                    POI = POIs[0]
                if self.filter_negatives:
                    logging.warning(
                        'deltaNLL<{0}; Dropping entry (deltaNLL={1:+10.4f}, {2}={3:+10.4f}) (scandirs: {4})'
                        .format(
                            self.deltaNLL_threshold,
                            entry.deltaNLL,
                            POI,
                            getattr(entry, POI),
                            self.scandirs)
                        )
                    continue
                else:
                    raise RuntimeError('Not allowed to filter negatives, but found:',entry)
            passed_entries.append(entry)

        if inplace:
            self.entries = passed_entries
        else:
            return passed_entries

    def x(self):
        return [ entry.x for entry in self.entries ]

    def y(self):
        return [ entry.y for entry in self.entries ]

    def z(self):
        return [ entry.z for entry in self.entries ]

    def deltaNLL(self):
        return [ entry.deltaNLL for entry in self.entries ]

    def two_times_deltaNLL(self):
        return [ 2.*entry.deltaNLL for entry in self.entries ]

    def bestfit(self):
        deltaNLLs = self.deltaNLL()
        if 0.0 in deltaNLLs:
            i_bestfit = deltaNLLs.index(0.0)
        else:
            abs_deltaNLLs = map(abs, deltaNLLs)
            i_bestfit = abs_deltaNLLs.index(min(abs_deltaNLLs))
            logging.error(
                'Could not find deltaNLL == 0.; taking value closest to 0.0: {0} (scandirs: {1})'
                .format(deltaNLLs[i_bestfit], self.scandirs)
                )
        return self.entries[i_bestfit]

    def filter(self, fn):
        bestfit = self.bestfit()
        new_entries = []
        for entry in self.entries:
            if entry is bestfit or fn(entry): new_entries.append(entry)
        self.entries = new_entries


class Scan(ScanPrimitive):
    """docstring for Scan"""

    uncertainty_calculator = UncertaintyCalculator()

    def __init__(self, x_variable, y_variable='deltaNLL', scandir=None, globpat=None, read_immediately=False):
        super(Scan, self).__init__()
        self.x_variable = x_variable
        self.y_variable = y_variable
        if not(scandir is None): self.scandirs.append(scandir)

        if not(globpat is None):
            self.globpat = '*' + globpat + '*'

        # self.x_title = self.standard_titles.get(x_variable, x_variable)
        self.x_title = core.standard_titles.get(x_variable, x_variable)
        self.y_title = '2#DeltaNLL'
        self.title   = self.x_title

        logging.debug(
            'New scan instance:'
            '\nx_variable = {0}'
            '\ny_variable = {1}'
            '\nscandirs   = {2}'
            '\nglobpat    = {3}'
            .format(self.x_variable, self.y_variable, ', '.join(self.scandirs), self.globpat)
            )

        if read_immediately: self.read()


    def read(self, keep_chain=False):
        root_files = self.collect_root_files()
        if self.save_all_variables:
            variables = self.get_list_of_variables_in_tree(root_files)
            if not(self.x_variable in variables) or not(self.y_variable in variables):
                raise RuntimeError(
                    'Variables x:{0} and/or y:{1} are not in the found list of variables:\n'
                    .format(self.x_variable, self.y_variable)
                    + '\n'.join(variables)
                    )
        else:
            variables = [ self.x_variable, self.y_variable ]
        if keep_chain:
            self.entries, self.chain = self.read_chain(root_files, variables, return_chain=True)
        else:
            self.entries = self.read_chain(root_files, variables)
        self.filter_entries()


    def filter_entries(self):
        super(Scan, self).filter_entries()

        scanfilter = OneDimScanFilter(self.x(), self.y())
        scanfilter.filter_clear_nonsense()
        scanfilter.raise_by_minimum()

        # Rebuild from the filtered x
        new_entries = []
        for ix, x in enumerate(scanfilter.xs):
            for entry in self.entries:
                if entry.x == x:
                    entry.y = scanfilter.ys[ix]
                    entry.deltaNLL = scanfilter.ys[ix]
                    new_entries.append(entry)
        self.entries = new_entries


    def create_uncertainties(self, inplace=True):
        logging.debug('Determining uncertainties for scan (x={0}, y={1})'.format(self.x_variable, self.y_variable))
        unc = self.uncertainty_calculator.create_uncertainties(xs = self.x(), deltaNLLs = self.deltaNLL())
        if unc.is_hopeless:
            logging.error(
                'Hopeless interpolation case: Unable to determine uncertainties for '
                'x = {0}, y = {1}, scandirs = {2}'
                .format(self.x_variable, self.y_variable, self.scandirs)
                )
        if inplace:
            self.unc = unc
        else:
            return unc

    def multiply_x_by_constant(self, constant):
        # Very ugly way of doing things... should not have used a namedtuple in the first place.
        for i in xrange(len(self.entries)):
            self.entries[i] = self.entries[i]._replace(x = self.entries[i].x * constant)

    def to_graph(self):
        name = utils.get_unique_rootname()
        graph = plotting.pywrappers.Graph(name, self.title, self.x(), self.two_times_deltaNLL())
        if hasattr(self, 'color'): graph.color = self.color
        if hasattr(self, 'draw_style'): graph.draw_style = self.draw_style
        return graph

    def get_spline_factory(self, x_min, x_max, cutstring_addition=''):
        factory = Spline2DFactory()
        factory.x_var = self.x_variable
        factory.z_var = self.y_variable
        factory.x_min = x_min
        factory.x_max = x_max
        factory.cutstring_addition = cutstring_addition
        factory.tree = self.chain
        bestfit = self.bestfit()
        factory.fill_bestfit(bestfit.x)
        return factory

    def to_spline(self, x_min, x_max, eps=2.2, deltaNLL_cutoff=30., cutstring_addition=''):
        factory = self.get_spline_factory(x_min, x_max, cutstring_addition)
        factory.eps = eps
        factory.deltaNLL_cutoff = deltaNLL_cutoff
        spline = factory.make_spline_1D()
        return spline



class Scan2D(ScanPrimitive):
    """docstring for Scan2D"""

    x_sm = 1.0
    y_sm = 1.0

    def __init__(self, name, x_variable, y_variable, z_variable='deltaNLL', scandir=None, globpat='*', color=None):
        super(Scan2D, self).__init__()
        self.x_variable = x_variable
        self.y_variable = y_variable
        self.z_variable = z_variable
        if not(scandir is None): self.scandirs.append(scandir)

        self.name = name
        self.title = self.standard_titles.get(self.name, self.name)

        self.x_title = self.standard_titles.get(x_variable, x_variable)
        self.y_title = self.standard_titles.get(y_variable, y_variable)
        self.z_title = '2#DeltaNLL'

        if globpat == '*':
            self.globpat = '*'
        else:
            self.globpat = '*' + globpat + '*'
        self.color = color

        self.draw_bestfit_point = True

        logging.debug(
            'New scan instance:'
            '\nx_variable = {0}'
            '\ny_variable = {1}'
            '\nz_variable = {2}'
            '\nscandirs   = {3}'
            '\nglobpat    = {4}'
            .format(self.x_variable, self.y_variable, self.z_variable, ', '.join(self.scandirs), self.globpat)
            )

    def read(self):
        root_files = self.collect_root_files()

        if self.save_all_variables:
            variables = self.get_list_of_variables_in_tree(self, root_files)
            if not(self.x_variable in variables) or not(self.y_variable in variables) or not(self.z_variable in variables):
                raise RuntimeError(
                    'Variables x:{0} and/or y:{1} and/or z:{2} are not in the found list of variables:\n'
                    .format(self.x_variable, self.y_variable, self.z_variable)
                    + '\n'.join(variables)
                    )
        else:
            variables = [ self.x_variable, self.y_variable, self.z_variable ]

        self.entries, self.chain = self.read_chain(root_files, variables, return_chain=True)
        self.filter_entries()

    # def bestfit(self):
    #     return self.entries[self.deltaNLL().index(0.0)]


    def get_spline_factory(self, x_min, x_max, y_min, y_max, cutstring_addition=''):
        factory = Spline2DFactory()
        factory.x_var = self.x_variable
        factory.y_var = self.y_variable
        factory.z_var = self.z_variable
        factory.x_min = x_min
        factory.x_max = x_max
        factory.y_min = y_min
        factory.y_max = y_max
        factory.cutstring_addition = cutstring_addition
        factory.tree = self.chain
        bestfit = self.bestfit()
        factory.fill_bestfit(bestfit.x, bestfit.y)
        return factory

    def to_spline(self, x_min, x_max, y_min, y_max, eps=2.2, deltaNLL_cutoff=30., cutstring_addition=''):
        factory = self.get_spline_factory(x_min, x_max, y_min, y_max, cutstring_addition)
        factory.eps = eps
        factory.deltaNLL_cutoff = deltaNLL_cutoff
        spline = factory.make_spline()
        return spline

    def to_polyfit(self, x_min, x_max, y_min, y_max, cutstring_addition=''):
        factory = self.get_spline_factory(x_min, x_max, y_min, y_max, cutstring_addition)
        polyfit = factory.make_polyfit()
        return polyfit


    def to_hist(self):
        histogram2D = plotting.pywrappers.Histogram2D(
            utils.get_unique_rootname(), getattr(self, 'title', ''), self.color
            )
        histogram2D.x_title = self.x_title
        histogram2D.y_title = self.y_title
        histogram2D.z_title = self.z_title
        histogram2D.fill_from_entries(self.entries)
        if hasattr(self, 'contour_filter_method'):
            histogram2D.contour_filter_method = self.contour_filter_method
        return histogram2D


    def plot(self, plotname, draw_style='repr_2D_with_contours'):
        c.Clear()
        c.set_margins_2D()

        # leg = plotting.pywrappers.Legend(
        #     c.GetLeftMargin() + 0.01,
        #     c.GetBottomMargin() + 0.02,
        #     1 - c.GetRightMargin() - 0.01,
        #     c.GetBottomMargin() + 0.09
        #     )
        histogram2D = self.to_hist()
        # histogram2D._legend = leg
        histogram2D.Draw(draw_style)

        base = histogram2D.H2
        base.GetXaxis().SetTitle(self.x_title)
        base.GetYaxis().SetTitle(self.y_title)
        base.GetXaxis().SetTitleSize(0.06)
        base.GetXaxis().SetLabelSize(0.05)
        base.GetYaxis().SetTitleSize(0.06)
        base.GetYaxis().SetLabelSize(0.05)

        plotting.pywrappers.CMS_Latex_type().Draw()
        plotting.pywrappers.CMS_Latex_lumi().Draw()
        plotting.pywrappers.Point(self.x_sm, self.y_sm).Draw('repr_SM_point')

        cdl = plotting.pywrappers.ContourDummyLegend(
            c.GetLeftMargin() + 0.01,
            1. - c.GetTopMargin() - 0.1,
            1. - c.GetRightMargin() - 0.01,
            1. - c.GetTopMargin() - 0.01,
            )
        if draw_style == 'repr_2D_with_contours_no_bestfit':
            cdl.disable_bestfit = True
        cdl.Draw()

        c.Update()
        c.RedrawAxis()
        c.save(plotname)


    def get_1d(self, variable):
        logging.info('Creating 1D profile from 2D histogram for variable {0}'.format(variable))
        histogram2D = self.to_hist()
        logging.trace('Will use 2D histogram with the following array:\n{0}'.format(histogram2D.H2_array))

        if variable == self.x_variable:
            logging.debug('Making 1D graph for x variable {0}'.format(variable))
            graph = self.get_1d_x(histogram2D)
        elif variable == self.y_variable:
            logging.debug('Making 1D graph for y variable {0}'.format(variable))
            graph = self.get_1d_y(histogram2D)
        else:
            raise ValueError('Variable {0} is not x {1} or y {2}'.format(variable, self.x_variable, self.y_variable))

        logging.trace('Found x values for 1D graph: {0}'.format(graph.xs))
        logging.trace('Found deltaNLL values for 1D graph: {0}'.format(graph.ys))

        logging.debug('Creating uncertainties for graph')
        graph.unc = Scan.uncertainty_calculator.create_uncertainties(
            graph.xs,
            [ 0.5*y for y in graph.ys ] # Histogram has 2dNLL, but unc calculator expects dNLL
            )
        return graph

    def get_1d_x(self, histogram2D):
        xs = []
        deltaNLLs = []
        for i_center, center in enumerate(histogram2D.x_bin_centers):
            deltaNLL = min(histogram2D.H2_array[i_center][:])
            xs.append(center)
            deltaNLLs.append(deltaNLL)

        graph = plotting.pywrappers.Graph(
            utils.get_unique_rootname(),
            self.x_variable,
            xs,
            deltaNLLs,
            color=self.color
            )
        return graph

    def get_1d_y(self, histogram2D):
        xs = []
        deltaNLLs = []
        for i_center, center in enumerate(histogram2D.y_bin_centers):
            deltaNLL = min([ row[i_center] for row in histogram2D.H2_array ])
            xs.append(center)
            deltaNLLs.append(deltaNLL)

        graph = plotting.pywrappers.Graph(
            utils.get_unique_rootname(),
            self.y_variable,
            xs,
            deltaNLLs,
            color=self.color
            )
        return graph

    #____________________________________________________________________
    # Not nice code here

    def get_thetas(self):
        histogram2D = self.to_hist()
        xs = []
        ys = []
        twodeltaNLLs = []
        for i_center, x in enumerate(histogram2D.x_bin_centers):
            all_deltaNNLs = histogram2D.H2_array[i_center][:]
            all_ys = histogram2D.y_bin_centers

            min_deltaNLL = 9999.
            for y, deltaNLL in zip(all_ys, all_deltaNNLs):
                y_on_line = -1./12. * x
                if y < 0.8 * y_on_line: continue
                if deltaNLL < min_deltaNLL:
                    min_deltaNLL = deltaNLL
                    min_y = y

            ys.append(min_y)
            xs.append(x)
            twodeltaNLLs.append(min_deltaNLL)
        thetas = [ math.atan2(y, x) for x, y in zip(xs, ys) ]
        thetas, twodeltaNLLs = [ list(l) for l in zip(*sorted(zip(thetas, twodeltaNLLs))) ]

        graph = plotting.pywrappers.Graph(
            utils.get_unique_rootname(),
            'theta',
            thetas,
            twodeltaNLLs,
            color=self.color
            )
        return graph

    def get_thetas_along_x(self, x_min, x_max, y):
        scanner = ThetaScanner(self.to_hist())
        return scanner.get_thetas_along_x(x_min, x_max, y)

    def get_thetas_along_y(self, y_min, y_max, x):
        scanner = ThetaScanner(self.to_hist())
        return scanner.get_thetas_along_y(y_min, y_max, x)


class ThetaScanner(object):
    """docstring for ThetaScanner"""
    def __init__(self, H):
        super(ThetaScanner, self).__init__()
        self.H = H

    def closest_match(self, x, list_of_bin_centers):
        dx_min = 9999.
        i_min = -1
        x_best = 9999
        for i, x_bin_center in enumerate(list_of_bin_centers):
            dx = abs(x-x_bin_center)
            if dx < dx_min:
                dx_min = dx
                i_min = i
                x_best = x_bin_center
        return i_min, x_best

    def closest_match_x(self, x):
        return self.closest_match(x, self.H.x_bin_centers)

    def closest_match_y(self, y):
        return self.closest_match(y, self.H.y_bin_centers)

    def theta(self, x, y):
        return math.atan2(y, x)

    def sort_thetas(self, thetas, deltaNLLs):
        thetas, deltaNLLs = [ list(l) for l in zip(*sorted(zip(thetas, deltaNLLs))) ]
        return thetas, deltaNLLs

    def get_thetas_along_x(self, x_min, x_max, y):
        ix_min, x_min = self.closest_match_x(x_min)
        ix_max, x_max = self.closest_match_x(x_max)
        iy, y = self.closest_match_y(y)
        logging.info(
            'Getting theta from x = {0} to x = {1} at y = {2}'
            .format(x_min, x_max, y)
            )
        if ix_min > ix_max:
            (ix_min, ix_max) = (ix_max, ix_min)

        thetas = []
        deltaNLLs = []
        for ix in xrange(ix_min, ix_max+1):
            x = self.H.x_bin_centers[ix]
            deltaNLL = self.H.H2_array[ix][iy]

            thetas.append(self.theta(x, y))
            deltaNLLs.append(deltaNLL)

        thetas, deltaNLLs = self.sort_thetas(thetas, deltaNLLs)
        return thetas, deltaNLLs


    def get_thetas_along_y(self, y_min, y_max, x):
        iy_min, y_min = self.closest_match_y(y_min)
        iy_max, y_max = self.closest_match_y(y_max)
        ix, x = self.closest_match_x(x)
        logging.info(
            'Getting theta from y = {0} to y = {1} at x = {2}'
            .format(y_min, y_max, x)
            )
        if iy_min > iy_max:
            (iy_min, iy_max) = (iy_max, iy_min)

        thetas = []
        deltaNLLs = []
        for iy in xrange(iy_min, iy_max+1):
            y = self.H.y_bin_centers[iy]
            deltaNLL = self.H.H2_array[ix][iy]

            thetas.append(self.theta(x, y))
            deltaNLLs.append(deltaNLL)

        thetas, deltaNLLs = self.sort_thetas(thetas, deltaNLLs)
        return thetas, deltaNLLs




