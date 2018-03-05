import os.path
import glob, re, copy

import logging
import core
import plotting
from plotting.canvas import c
import plotting.plotting_utils as utils

from collections import namedtuple
from array import array

import ROOT


def glob_rootfiles(d):
    if not d.endswith('/'): d += '/'
    return glob.glob(d + '*.root')

def rindex( someList, val ):
    # Regular list.index() finds first instance in list, this function finds the last
    return len(someList) - someList[::-1].index(val) - 1


class DifferentialSpectrum(object):
    """Essentially a collection of Scan instances, 1 per POI that was scanned"""
    standard_titles = {
        'hgg' : 'H#rightarrow#gamma#gamma',
        'hzz' : 'H#rightarrowZZ',
        'combination' : 'Combination',
        'hbb' : 'H#rightarrowbb',
        'combWithHbb' : 'Comb. with H#rightarrowbb',
        }

    def __init__(self, name, scandir, auto_determine_POIs=True, datacard=None):
        self.name = name
        self.scandir = scandir
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

    def set_sm(self, smxs):
        self.smxs = smxs
        self.smxs_set = True

    def get_POIs_from_scandir(self):
        root_files = glob_rootfiles(self.scandir)
        POIs = []
        for root_file in root_files:
            match = re.search(r'bPOI_(\w+)_ePOI', os.path.basename(root_file))
            if not match:
                continue
            POIs.append(match.group(1))
        POIs = list(set(POIs))
        POIs.sort(key=core.range_sorter)
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
                    scandir=self.scandir,
                    globpat=POI
                    )
                scan.read()
                scan.create_uncertainties()
                self.scans.append(scan)
            self._is_read = True

    def plot_scans(self, plotname=None):
        if plotname is None:
            plotname = 'scans_{0}'.format(self.name)
        plot = plotting.plots.MultiScanPlot(plotname)
        plot.scans = self.scans
        if 'hbb' in self.name:
            plot.x_min = -10.0
            plot.x_max = 10.0
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

    standard_titles = {
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

    tree_name = 'limit'
    filter_negatives = True
    deltaNLL_threshold = -0.01

    def __init__(self):
        self.scandirs = []
        self.root_files = []
        self.save_all_variables = False
        self.entries = []
        self.globpat = '*'


    def collect_root_files(self):
        root_files = copy.copy(self.root_files)
        for scandir in self.scandirs:
            if not scandir.endswith('/'): scandir += '/'
            root_files.extend(glob.glob(scandir + self.globpat + '.root'))

        if len(root_files) == 0:
            raise RuntimeError(
                'Attemped to retrieve scan for x:{0} y:{1}, '
                'but no .root files were found. Passed list of dirs to look in:\n'
                .format(self.x_variable, self.y_variable)
                + '\n'.join(self.scandirs)
                )

        logging.debug('Found {0} root files in {1}'.format(len(root_files), ', '.join(self.scandirs)))
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
                if self.filter_negatives:
                    logging.warning(
                        'Dropping entry (deltaNLL<{0}:'.format(self.deltaNLL_threshold)
                        + entry.__repr__()
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
        return self.entries[self.deltaNLL().index(0.0)]


class Scan(ScanPrimitive):
    """docstring for Scan"""

    def __init__(self, x_variable, y_variable='deltaNLL', scandir=None, globpat=None):
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
        xs = self.x()
        deltaNLLs = self.deltaNLL()

        min_deltaNLL   = min(deltaNLLs)
        i_min_deltaNLL = rindex(deltaNLLs, min_deltaNLL)
        x_min          = xs[i_min_deltaNLL]

        logging.debug('Determining uncertainties for scan (x={0}, y={1})'.format(self.x_variable, self.y_variable))
        logging.debug('Found minimum at index {0}: x={1}, deltaNLL={2}'.format(i_min_deltaNLL, x_min, min_deltaNLL))

        unc_dict = {
            'min_deltaNLL' : min_deltaNLL,
            'i_min' : i_min_deltaNLL,
            'x_min' : x_min,
            'left_bound' : -0,
            'left_error' : -0,
            'right_bound' : 0,
            'right_error' : 0,
            'well_defined_left_bound' : False,
            'well_defined_right_bound' : False,
            }
        Unc = namedtuple('Unc', unc_dict.keys())

        # Process left uncertainty
        if i_min_deltaNLL < 3:
            logging.debug('Not enough points on the left side for a well defined left bound')
            well_defined_left_bound = False
        else:
            xs_left = xs[:i_min_deltaNLL+1]
            deltaNLLs_left = deltaNLLs[:i_min_deltaNLL+1]
            if min(deltaNLLs_left) > 0.5 or max(deltaNLLs_left) < 0.5:
                logging.debug('Requested dNLL interpolation point is outside the range: min dNLL={0}, max dNLL={1}'.format(min(deltaNLLs_left), max(deltaNLLs_left)))
                well_defined_left_bound = False
            else:
                left_bound = self.interpolate(xs_left, deltaNLLs_left, 0.5)
                if left_bound is False:
                    well_defined_left_bound = False
                else:
                    well_defined_left_bound = True

        # Process right uncertainty
        if i_min_deltaNLL > len(xs)-3:
            logging.debug('Not enough points on the right side for a well defined right bound')
            well_defined_right_bound = False
        else:
            xs_right = xs[i_min_deltaNLL:]
            deltaNLLs_right = deltaNLLs[i_min_deltaNLL:]
            if min(deltaNLLs_right) > 0.5 or max(deltaNLLs_right) < 0.5:
                logging.debug('Requested dNLL interpolation point is outside the range: min dNLL={0}, max dNLL={1}'.format(min(deltaNLLs_right), max(deltaNLLs_right)))
                well_defined_right_bound = False
            else:
                right_bound = self.interpolate(xs_right, deltaNLLs_right, 0.5)
                if right_bound is False:
                    well_defined_right_bound = False
                else:
                    well_defined_right_bound = True

        is_hopeless = False
        if well_defined_left_bound and well_defined_right_bound:
            pass
        elif well_defined_left_bound and not well_defined_right_bound:
            right_bound = x_min + (x_min - left_bound)
        elif well_defined_right_bound and not well_defined_left_bound:
            left_bound  = x_min - (right_bound - x_min)
        else:
            logging.error(
                'Hopeless interpolation case; unable to determine uncertainties for '
                'x = {0}, y = {1}'
                .format(self.x_variable, self.y_variable)
                )
            is_hopeless = True

        if not is_hopeless:
            unc_dict['well_defined_left_bound'] = well_defined_left_bound
            unc_dict['well_defined_right_bound'] = well_defined_right_bound
            left_error = abs(x_min - left_bound)
            right_error = abs(x_min - right_bound)
            unc_dict['left_bound']  = left_bound
            unc_dict['left_error']  = left_error
            unc_dict['right_bound'] = right_bound
            unc_dict['right_error'] = right_error

        unc = Unc(**unc_dict)
        if inplace:
            self.unc = unc
        else:
            return unc

    def interpolate(self, ys, xs, x_value):
        logging.debug('Interpolating for x_value={0}'.format(x_value))
        logging.trace('  x  /  y:')
        for x, y in zip(xs, ys): logging.trace('    {0:+7.2f}  /  {1:+7.2f}'.format(x, y))

        if min(xs) > x_value or max(xs) < x_value:
            logging.debug('  Requested interpolation {0} is outside the range: {1} to {2}'.format(x_value, min(xs), max(xs)))
            return False
        Tg = ROOT.TGraph(len(xs), array('f', xs), array('f', ys))
        y_value = Tg.Eval(x_value)
        # if y_value < min(ys) or y_value > max(ys):
        #     return False
        if y_value is False:
            logging.debug('  Could not interpolate properly')
        else:
            logging.debug('  Interpolated y_value {0} for x_value {1}'.format(y_value, x_value))
        return y_value


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
        return histogram2D


    def plot(self, plotname):
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
        histogram2D.Draw('repr_2D_with_contours')

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

