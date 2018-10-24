import sys
import core
import ROOT
import logging

from plotting.plotting_utils import get_unique_rootname
from plotting.pywrappers import Histogram2D, Graph


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

        self.cutstring_addition = ''

        # Options for 2d fit
        self.ord_polynomial = 7


    def fill_bestfit(self, x, y=None):
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
            'deltaNLL>0 && deltaNLL<{0}'.format(self.deltaNLL_cutoff) + self.cutstring_addition,
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

    def make_spline_1D(self):
        x = self.make_var_unique_name(self.x_var, self.x_min, self.x_max)
        s = ROOT.RooArgList(x)
        ROOT.SetOwnership(s, False)

        spline = ROOT.RooSplineND(
            'spline_' + get_unique_rootname(), 'spline1D',
            s,
            self.tree,
            self.z_var,
            self.eps,
            True, # boolean rescale
            'deltaNLL>0 && deltaNLL<{0}'.format(self.deltaNLL_cutoff) + self.cutstring_addition,
            )
        splinewrapper = Spline1DWrapper(spline, x,
            x_min = self.x_min,
            x_max = self.x_max,
            )

        if self.filled_bestfit:
            splinewrapper.fill_bestfit(self.x_bestfit)

        return splinewrapper


    #____________________________________________________________________
    def compile_2d_fit_string(self):
        fstring = '('
        param = 0
        for o1 in range(self.ord_polynomial+1):
            if (o1==0):
                fstring += "["+str(param)+"]"
                param+=1
                continue
            else:
                fstring += "+["+str(param)+"]*pow(x,"+str(o1)+")"
                param+=1
                fstring += "+["+str(param)+"]*pow(y,"+str(o1)+")"
                param+=1

                o2=1
                while o2<(o1+1):
                    fstring += "+["+str(param)+"]*pow(x,"+str(o1)+")*pow(y,"+str(o2)+")"
                    param+=1
                    if (not o1==o2):
                        fstring += "+["+str(param)+"]*pow(x,"+str(o2)+")*pow(y,"+str(o1)+")"
                        param+=1
                    o2+=1
        fstring += ")"
        return fstring

    def get_graph(self):
        # graph = plot.TGraph2DFromTree(
        #     limit,
        #     param_x,
        #     param_y,
        #     '2*deltaNLL',
        #     'quantileExpected > -0.5 && deltaNLL > 0 && 2.0*deltaNLL<10.0'
        #     )
        cutstring = 'deltaNLL > 0. && 2.0*deltaNLL<{0}'.format(self.deltaNLL_cutoff) + self.cutstring_addition
        logging.info('Creating TGraph2D; cutstring is "{0}"'.format(cutstring))
        graph = TGraph2DFromTree(
            self.tree,
            self.x_var,
            self.y_var,
            '2.0*deltaNLL',
            cutstring
            )
        return graph

    def make_polyfit(self, graph=None):
        logging.info('Creating polyfit')
        fstring = self.compile_2d_fit_string()
        if graph is None:
            graph = self.get_graph()

        logging.info('fit string is: ' + fstring)

        f2D = ROOT.TF2(
            'f2D' + get_unique_rootname(), fstring,
            self.x_min, self.x_max,
            self.y_min, self.y_max
            )
        logging.info(
            'Fitting TGraph2D ({0} x {1} = {2} points) with TF2...'
            .format(graph.GetNpx(), graph.GetNpy(), graph.GetN())
            )
        graph.Fit(f2D.GetName())

        polyfit = PolyFit2DWrapper(f2D,
            x_min = self.x_min,
            x_max = self.x_max,
            y_min = self.y_min,
            y_max = self.y_max,
            )
        polyfit.multiply_by_two = True

        if self.filled_bestfit:
            polyfit.fill_bestfit(self.x_bestfit, self.y_bestfit)

        return polyfit

    def make_polyfit_for_patch(x_min, x_max, y_min, y_max):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.cutstring_addition = (
            '&& {x}>{x_min} && {x}<{x_max} && {y}>{y_min} && {y}<{y_max}'
            .format(
                x = self.x_var,
                y = self.y_var,
                x_min = self.x_min,
                x_max = self.x_max,
                y_min = self.y_min,
                y_max = self.y_max,
                )
            )
        return self.make_polyfit()


def TGraph2DFromTree(tree, xvar,  yvar, zvar, selection):
    tree.Draw(xvar + ':' + yvar + ':' + zvar, selection, 'goff')
    gr = ROOT.TGraph2D(
        tree.GetSelectedRows(), tree.GetV1(), tree.GetV2(), tree.GetV3())
    ROOT.SetOwnership(gr, False)
    return gr


#____________________________________________________________________
class BaseWrapper(object):
    """docstring for BaseWrapper"""

    multiply_by_two = True
    disallow_negativity = True
    negativity_is_zero = False

    def __init__(self):
        super(BaseWrapper, self).__init__()
        self.noise_selectors = []
        self.signal_selectors = []
        self.filled_bestfit = False

    def add_noise_selector(self, selector):
        self.noise_selectors.append(selector)

    def add_signal_selector(self, selector):
        self.signal_selectors.append(selector)

    def name(self):
        return get_unique_rootname()

class Base1DWrapper(BaseWrapper):
    """docstring for Base1DWrapper"""
    def __init__(self, x_min, x_max):
        super(Base1DWrapper, self).__init__()
        self.x_min = x_min
        self.x_max = x_max
        
    def fill_bestfit(self, x):
        self.filled_bestfit = True
        self.x_bestfit = x

    def eval(self, x):
        if not(self.x_min is None) and x < self.x_min: return 999.
        if not(self.x_max is None) and x > self.x_max: return 999.

        if len(self.signal_selectors) > 0:
            for selector in self.signal_selectors:
                if selector(x): return 0.0

        if len(self.noise_selectors) > 0:
            for selector in self.noise_selectors:
                if selector(x): return 999.0

        r = self.eval_interp(x)

        if self.multiply_by_two: r *= 2

        if r < 0.:
            if self.negativity_is_zero:
                r = 0.00001
            elif self.disallow_negativity:
                r = 999.
        return r


class Spline1DWrapper(Base1DWrapper):
    """docstring for Spline1DWrapper"""
    def __init__(self, spline, x, x_min, x_max):
        super(Spline1DWrapper, self).__init__(x_min, x_max)
        self.spline = spline
        self.x = x

    def name(self):
        return self.spline.GetName().replace('spline_', 'hist_')

    def eval_interp(self, x):
        self.x.setVal(x)
        r = self.spline.getVal()
        return r

    def eval(self, x):
        return super(Spline1DWrapper, self).eval(x)

    def to_graph(self, nx=100):
        x_axis = core.get_axis(nx, self.x_min, self.x_max)
        y_axis = []
        for x in x_axis:
            y_axis.append(self.eval(x))

        graph = Graph('auto', 'new_spline', x_axis, y_axis)
        return graph
        
        

class Base2DWrapper(BaseWrapper):
    def __init__(self, x_min=None, x_max=None, y_min=None, y_max=None):
        super(Base2DWrapper, self).__init__()
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.offset = None

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

        # Not implemented in Base!
        r = self.eval_interp(x, y)

        if not(self.offset is None):
            r += self.offset

        if self.multiply_by_two:
            r *= 2

        if r < 0.:
            if self.negativity_is_zero:
                r = 0.00001
            elif self.disallow_negativity:
                r = 999.
        return r

    def to_hist(self, nx=100, ny=100, x_min=None, x_max=None, y_min=None, y_max=None):
        """Take standard spline ranges by default, but allow smaller or bigger rangers"""
        name = self.name()
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

        if hasattr(self, '_scan'):
            H.color = self._scan.color
            H.name  = self._scan.name + '_splined'
            H.title = self._scan.title
        return H


class PolyFit2DWrapper(Base2DWrapper):
    """docstring for PolyFit2DWrapper"""
    def __init__(self, f2D, x_min=None, x_max=None, y_min=None, y_max=None):
        super(PolyFit2DWrapper, self).__init__(x_min, x_max, y_min, y_max)
        self.f2D = f2D

    def name(self):
        return 'f2D_' + get_unique_rootname()

    def eval_interp(self, x, y):
        return self.f2D.Eval(x, y)


class Spline2DWrapper(Base2DWrapper):
    """docstring for Spline2DWrapper"""
    def __init__(self, spline, x, y, x_min=None, x_max=None, y_min=None, y_max=None):
        super(Spline2DWrapper, self).__init__(x_min, x_max, y_min, y_max)
        self.spline = spline
        self.x = x
        self.y = y

    def name(self):
        return self.spline.GetName().replace('spline_', 'hist_')

    def eval_interp(self, x, y):
        self.x.setVal(x)
        self.y.setVal(y)
        r = self.spline.getVal()
        return r

    def eval(self, x, y):
        return super(Spline2DWrapper, self).eval(x, y)

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






