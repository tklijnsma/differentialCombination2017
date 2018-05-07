import os.path
import glob, re, copy

import logging
import core
import plotting
from plotting.canvas import c
import plotting.plotting_utils as utils
from uncertaintycalculator import UncertaintyCalculator

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
        'combWithHbb' : 'Comb. with H#rightarrowbb',
        }

    def __init__(self, name, scandirs, auto_determine_POIs=True, datacard=None):
        self.name = name
        if isinstance(scandirs, basestring): scandirs = [scandirs]
        self.scandirs = scandirs
        self.scans = []

        self.auto_determine_POIs = auto_determine_POIs
        self.datacard = datacard
        if not(datacard is None):
            self.auto_determine_POIs = False
        elif not(auto_determine_POIs) and datacard is None:
            raise ValueError('Either auto_determine_POIs must be True, or a workspace must be passed')

        self.POIs = []
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
        else:
            self.POIs = self.get_POIs_from_datacard()

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
        plot.draw()
        plot.wrapup()

    def binning(self):
        binning = core.binning_from_POIs(self.POIs)
        return binning

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

    def last_bin_is_overflow(self):
        if self.binning()[-1] == 10000.:
            return True
        return core.last_bin_is_overflow(self.POIs)

    def first_bin_is_underflow(self):
        return core.first_bin_is_underflow(self.POIs)

    def process_overflow_and_underflow(self, histogram):
        if self.last_bin_is_overflow():
            if not(self.hard_x_max is None):
                histogram.set_last_bin_is_overflow(method='HARDVALUE', hard_value=self.hard_x_max)
            else:
                histogram.set_last_bin_is_overflow()
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
        if not(self.color is None): histogram.color = self.color
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
        histogram.color = 16
        return histogram

    def to_hist_smxs(self):
        histogram = plotting.pywrappers.Histogram(
            utils.get_unique_rootname(),
            'SM',
            self.binning(),
            self.smxs
            )
        # histogram.set_err_up([ s.unc.right_error * xs for s, xs in zip(self.scans, self.smxs) ])
        # histogram.set_err_down([ s.unc.left_error * xs for s, xs in zip(self.scans, self.smxs) ])
        self.process_overflow_and_underflow(histogram)
        histogram.color = 16
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
        'combWithHbb' : 'Comb. with H#rightarrowbb',
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

    def read_chain(self, root_files, variables, filter_x=False):
        chain = ROOT.TChain(self.tree_name)
        for root_file in root_files:
            chain.Add(root_file)

        Entry = namedtuple('Entry', variables + ['x', 'y'] + (['z'] if hasattr(self, 'z_variable') else []))
        entries = []
        for event in chain:
            if filter_x and getattr(event, self.x_variable) in [e.x for e in entries]: continue
            entry_dict = {}
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
            entries.append(Entry(**entry_dict))

        if len(entries) == 0:
            raise RuntimeError(
                'No entries were found'
                )

        entries.sort(key=lambda entry: (entry.x, entry.y))
        return entries

    def filter_entries(self, inplace=True):
        passed_entries = []
        for entry in self.entries:
            if entry.deltaNLL < self.deltaNLL_threshold:
                POIs = [ k for k in entry._fields if k.startswith('r_') ]
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


    def read(self):
        root_files = self.collect_root_files()

        if self.save_all_variables:
            variables = self.get_list_of_variables_in_tree(self, root_files)
            if not(self.x_variable in variables) or not(self.y_variable in variables):
                raise RuntimeError(
                    'Variables x:{0} and/or y:{1} are not in the found list of variables:\n'
                    .format(self.x_variable, self.y_variable)
                    + '\n'.join(variables)
                    )
        else:
            variables = [ self.x_variable, self.y_variable ]

        self.entries = self.read_chain(root_files, variables)
        self.filter_entries()

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

        self.entries = self.read_chain(root_files, variables)
        self.filter_entries()

    # def bestfit(self):
    #     return self.entries[self.deltaNLL().index(0.0)]

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

        plotting.pywrappers.ContourDummyLegend(
            c.GetLeftMargin() + 0.01,
            1. - c.GetTopMargin() - 0.1,
            1. - c.GetRightMargin() - 0.01,
            1. - c.GetTopMargin() - 0.01,
            ).Draw()

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

