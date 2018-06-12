import logging
import core
import scans
import plotting.pywrappers
import uncertaintycalculator

class OneDimScanner(object):
    """docstring for OneDimScanner"""

    def __init__(self, histogram2D=None, x_variable=None, y_variable=None):
        super(OneDimScanner, self).__init__()
        if isinstance(histogram2D, scans.Scan2D):
            self.histogram2D = histogram2D.to_hist()
            self.x_variable = histogram2D.x_variable
            self.y_variable = histogram2D.y_variable
        else:
            self.histogram2D = histogram2D    
            self.x_variable = x_variable
            self.y_variable = y_variable

        self.uncertaintycalculator = uncertaintycalculator.UncertaintyCalculator()


    def get_minima_along_x(self):
        xs = []
        ys = []
        deltaNLLs = []
        for i_center, center in enumerate(self.histogram2D.x_bin_centers):
            deltaNLL_thiscolumn = self.histogram2D.H2_array[i_center][:]
            deltaNLL = min(deltaNLL_thiscolumn)
            y =  self.histogram2D.y_bin_centers[ deltaNLL_thiscolumn.index(deltaNLL) ]

            xs.append(center)
            ys.append(y)
            deltaNLLs.append(deltaNLL)
        return xs, ys, deltaNLLs

    def get_minima_along_y(self):
        xs = []
        ys = []
        deltaNLLs = []
        for i_center, center in enumerate(self.histogram2D.y_bin_centers):
            deltaNLL_thisrow = [ row[i_center] for row in self.histogram2D.H2_array ]
            deltaNLL = min(deltaNLL_thisrow)
            x =  self.histogram2D.x_bin_centers[ deltaNLL_thisrow.index(deltaNLL) ]

            xs.append(x)
            ys.append(center)
            deltaNLLs.append(deltaNLL)
        return xs, ys, deltaNLLs


    def get_1d(self, variable):
        logging.info('Creating 1D profile from 2D histogram for variable {0}'.format(variable))
        logging.trace('Will use 2D histogram with the following array:\n{0}'.format(self.histogram2D.H2_array))

        if variable == self.x_variable:
            logging.debug('Making 1D graph for x variable {0}'.format(variable))
            xs, ys, deltaNLLs = self.get_minima_along_x()
            graph = plotting.pywrappers.Graph(
                'auto',
                core.standard_titles.get(self.x_variable, self.x_variable),
                xs,
                deltaNLLs
                )

        elif variable == self.y_variable:
            logging.debug('Making 1D graph for y variable {0}'.format(variable))
            xs, ys, deltaNLLs = self.get_minima_along_y()
            graph = plotting.pywrappers.Graph(
                'auto',
                core.standard_titles.get(self.y_variable, self.y_variable),
                ys,
                deltaNLLs
                )

        else:
            raise ValueError('Variable {0} is not x {1} or y {2}'.format(variable, self.x_variable, self.y_variable))

        logging.trace('Found x values for 1D graph: {0}'.format(graph.xs))
        logging.trace('Found deltaNLL values for 1D graph: {0}'.format(graph.ys))

        logging.debug('Creating uncertainties for graph')
        self.uncertaintycalculator.cutoff = 1.0 # uncertaintycalculator by default expects deltaNLL, hist has 2deltaNLL
        graph.unc = self.uncertaintycalculator.create_uncertainties(graph.xs, graph.ys)
        return graph

