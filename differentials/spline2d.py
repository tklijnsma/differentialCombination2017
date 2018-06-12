import sys
import core
import ROOT

from plotting.plotting_utils import get_unique_rootname
from plotting.pywrappers import Histogram2D


class Spline2DFactory(object):
    """docstring for Spline2DFactory"""
    def __init__(self):
        super(Spline2DFactory, self).__init__()

        self.x_var = 'ct'
        self.y_var = 'cg'
        self.z_var = 'deltaNLL'

        self.x_min = -2.0
        self.x_max = 2.0
        self.y_min = -2.0
        self.y_max = 2.0

        self.tree = None

        # Probably some splining parameter
        # In combine code it's 1.7
        self.eps = 2.2
        self.deltaNLL_cutoff = 30.
        self.filled_bestfit = False


    def fill_bestfit(self, x, y):
        self.filled_bestfit = True
        self.x_bestfit = x
        self.y_bestfit = y

    def make_var_unique_name(self, name, x_min, x_max):
        x = ROOT.RooRealVar(
            # name + get_unique_rootname(),
            name,
            name,
            x_min,
            x_max
            )
        x.name = name
        ROOT.SetOwnership(x, False)
        return x

    def make_spline(self):
        x = self.make_var_unique_name(self.x_var, self.x_min, self.x_max)
        y = self.make_var_unique_name(self.y_var, self.y_min, self.y_max)

        s = ROOT.RooArgList(x, y)
        ROOT.SetOwnership(s, False)

        spline = ROOT.RooSplineND(
            'spline_' + get_unique_rootname(), 'spline2D',
            s,
            self.tree,
            self.z_var,
            self.eps,
            True, # boolean rescale
            'deltaNLL>0 && deltaNLL<{0}'.format(self.deltaNLL_cutoff),
            )
        splinewrapper = Spline2DWrapper(spline, x, y,
            x_min = self.x_min,
            x_max = self.x_max,
            y_min = self.y_min,
            y_max = self.y_max,
            )

        if self.filled_bestfit:
            splinewrapper.fill_bestfit(self.x_bestfit, self.y_bestfit)

        return splinewrapper


class Spline2DWrapper(object):
    """docstring for Spline2DWrapper"""
    def __init__(self, spline, x, y, x_min=None, x_max=None, y_min=None, y_max=None):
        super(Spline2DWrapper, self).__init__()
        self.spline = spline
        self.x = x
        self.y = y
        self.multiply_by_two = True
        self.disallow_negativity = True
        self.negativity_is_zero = False
        self.noise_selectors = []
        self.signal_selectors = []
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.filled_bestfit = False

    def fill_bestfit(self, x, y):
        self.filled_bestfit = True
        self.x_bestfit = x
        self.y_bestfit = y

    def eval(self, x, y):
        if not(self.x_min is None) and x < self.x_min: return 999.
        if not(self.x_max is None) and x > self.x_max: return 999.
        if not(self.y_min is None) and y < self.y_min: return 999.
        if not(self.y_max is None) and y > self.y_max: return 999.

        if len(self.signal_selectors) > 0:
            for selector in self.signal_selectors:
                if selector(x, y):
                    return 0.0

        if len(self.noise_selectors) > 0:
            for selector in self.noise_selectors:
                if selector(x, y):
                    return 999.0

        self.x.setVal(x)
        self.y.setVal(y)
        r = self.spline.getVal()
        if self.multiply_by_two:
            r *= 2

        if r < 0.:
            if self.negativity_is_zero:
                r = 0.0
            elif self.disallow_negativity:
                r = 999.
        return r

    def add_noise_selector(self, selector):
        self.noise_selectors.append(selector)

    def add_signal_selector(self, selector):
        self.signal_selectors.append(selector)

    def to_graph(self, nx=100, ny=100):
        graph = TGraph2D();
        graph.SetName( self.spline.GetName().replace('spline_', 'graph_') )

        x_axis = core.get_axis(nx, self.x.getMin(), self.x.getMax())
        y_axis = core.get_axis(ny, self.y.getMin(), self.y.getMax())

        i_point = 0
        for ix, x in enumerate(x_axis):
            for iy, y in enumerate(y_axis):
                graph.SetPoint(i_point, x, y, self.eval(x, y))
                i_point += 1

        return graph

    def to_hist(self, nx=100, ny=100, x_min=None, x_max=None, y_min=None, y_max=None):
        """Take standard spline ranges by default, but allow smaller or bigger rangers"""
        name = self.spline.GetName().replace('spline_', 'hist_')
        H = Histogram2D(name, name)

        if x_min is None: x_min = self.x_min
        if x_max is None: x_max = self.x_max
        if y_min is None: y_min = self.y_min
        if y_max is None: y_max = self.y_max

        x_boundaries = core.get_axis(nx+1, x_min, x_max) # nx would be # of bin centers, not bounds
        y_boundaries = core.get_axis(ny+1, y_min, y_max)
        x_centers = [ 0.5*(l+r) for l, r in zip(x_boundaries[:-1], x_boundaries[1:]) ]
        y_centers = [ 0.5*(l+r) for l, r in zip(y_boundaries[:-1], y_boundaries[1:]) ]

        # Initialize with 0
        mat = [ [0.0 for iy in xrange(ny)] for ix in xrange(nx) ]

        for ix, x in enumerate(x_centers):
            for iy, y in enumerate(y_centers):
                mat[ix][iy] = self.eval(x,y)

        H.fill_with_matrix(mat, x_boundaries, y_boundaries)

        if self.filled_bestfit:
            H.fill_bestfit(self.x_bestfit, self.y_bestfit)
        return H











