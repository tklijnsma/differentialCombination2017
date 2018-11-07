#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

from OptionHandler import flag_as_option

import LatestPaths
import LatestBinning

import differentials
import differentialutils

import logging

from differentials.theory.theory_utils import FileFinder

from time import strftime
datestr = strftime('%b%d')

########################################
# Main
########################################


@flag_as_option
def test_param1(args):
    drawer = ParametrizationDrawer('ct', 'cg')
    drawer.ktcg_variations()
    drawer.make_parametrization()
    drawer.plot('parametrizations_nolinearterms')

    drawer_linterms = ParametrizationDrawer('ct', 'cg')
    drawer_linterms.do_linear_terms = True
    drawer_linterms.ktcg_variations()
    drawer_linterms.make_parametrization()
    drawer_linterms.plot('parametrizations_linearterms')

    # print 'w/o  interference; kt = 1.0, cg = 0: incl xs = {0}'.format(drawer.inclusive_xs_for_parametrization(1.0, 0.0))
    # print 'with interference; kt = 1.0, cg = 0: incl xs = {0}'.format(drawer_linterms.inclusive_xs_for_parametrization(1.0, 0.0))
    # print 'w/o  interference; kt = 2.0, cg = 0: incl xs = {0}'.format(drawer.inclusive_xs_for_parametrization(2.0, 0.0))
    # print 'with interference; kt = 2.0, cg = 0: incl xs = {0}'.format(drawer_linterms.inclusive_xs_for_parametrization(2.0, 0.0))

    print '\nCoefficients:'
    print 'w/o  interference: sigma = A*{0}^2 + B*{1}^2 + C*{0}{1}'.format('ct', 'cg')
    print 'with interference: sigma = A*{0}^2 + B*{1}^2 + C*{0}{1} + D*{0} + E*{1} + F'.format('ct', 'cg')
    for i in xrange(len(drawer.parametrization.parametrizations)):
        coeff_wo_intf = drawer.parametrization.parametrizations[i].coefficients
        coeff_with_intf = drawer_linterms.parametrization.parametrizations[i].coefficients
        print (
            '{0:{numformat}} {c1}^2 + {1:{numformat}} {c2}^2 + {2:{numformat}} {c1}{c2}  |  '
            '{3:{numformat}} {c1}^2 + {4:{numformat}} {c2}^2 + {5:{numformat}} {c1}{c2} + {6:{numformat}} {c1} + {7:{numformat}} {c2} + {8:{numformat}}'
            .format(
                *(coeff_wo_intf + coeff_with_intf),
                numformat = '+8.3f',
                c1 = 'ct', c2 = 'cg'
                )
            )


@flag_as_option
def test_param2(args):
    drawer = ParametrizationDrawer('ct', 'cg')
    drawer.ktcg_variations()
    drawer.make_parametrization()

    drawer_linterms = ParametrizationDrawer('ct', 'cg')
    drawer_linterms.do_linear_terms = True
    drawer_linterms.ktcg_variations()
    drawer_linterms.make_parametrization()


    plot = differentials.plotting.plots.BottomPanelPlot('param_crosscheck')
    plot.make_legend()
    plot.leg.SetNColumns(2)
    plot.top_x_min = drawer.sm['binCenters'][0]
    plot.top_x_max = drawer.sm['binCenters'][-1]
    plot.bottom_x_min = drawer.sm['binCenters'][0]
    plot.bottom_x_max = drawer.sm['binCenters'][-1]
    plot.top_y_min = 0.00001
    plot.top_y_max = 1.1
    plot.bottom_y_min = -0.5
    plot.bottom_y_max = 2.5

    points = [
        # ct   cg      cb
        ( 0.0, 1./12., 0.0 ),
        ( 1.0, 0.0,    1.0 ),
        ( 1.0, 0.0,    0.0 ),
        ]

    differentials.parametrization.Parabola.c3 = 0.0
    inf_mtop_xs, _ = drawer_linterms.get_eval_param_graph(0.0, 1./12.)

    for i_point, point in enumerate(points):
        ct, cg, cb = point
        color = drawer.colors[i_point]

        differentials.parametrization.Parabola.c3 = cb
        graph_xs, graph_ratio = drawer_linterms.get_eval_param_graph(ct, cg, color=color)

        # Don't take ratio to SM, but ratio to inf mtop limit
        graph_ratio.ys = [
            y / y_inf_mtop for y, y_inf_mtop in
            zip(graph_xs.ys, inf_mtop_xs.ys)
            ]

        graph_xs.title = 'ct = {0:+.1f}, cg = {1:+.3f}, cb = {2:+.1f}'.format(ct, cg, cb)
        plot.add_top(graph_xs, 'repr_basic_line', plot.leg)
        plot.add_bottom(graph_ratio, 'repr_basic_line')

    plot.draw()
    plot.wrapup()





@flag_as_option
def test_param3(args):
    drawer = ParametrizationDrawer('ct', 'cg', 'cb')
    drawer.do_linear_terms = False
    drawer.ktcgkb_variations()
    drawer.make_parametrization()

    plot = differentials.plotting.plots.BottomPanelPlot('param_crosscheck')
    plot.make_legend()
    plot.leg.SetNColumns(2)
    plot.top_x_min = drawer.sm['binCenters'][0]
    # plot.top_x_max = drawer.sm['binCenters'][-1]
    plot.bottom_x_min = drawer.sm['binCenters'][0]
    # plot.bottom_x_max = drawer.sm['binCenters'][-1]

    plot.top_x_max =300.
    plot.bottom_x_max = 300.

    plot.top_y_min = 0.00001
    plot.top_y_max = 1.1
    plot.bottom_y_min = 0.85
    plot.bottom_y_max = 1.35

    points = [
        # ct   cg      cb
        ( 0.0, 1./12., 0.0 ),
        ( 1.0, 0.0,    1.0 ),
        ( 1.0, 0.0,    0.0 ),
        ]

    inf_mtop_xs, _ = drawer.get_eval_param_graph(ct=0.0, cg=1./12., cb=0.0)

    for i_point, point in enumerate(points):
        ct, cg, cb = point
        color = drawer.colors[i_point]

        graph_xs, graph_ratio = drawer.get_eval_param_graph(ct=ct, cg=cg, cb=cb, color=color)

        # Don't take ratio to SM, but ratio to inf mtop limit
        graph_ratio.ys = [
            y / y_inf_mtop for y, y_inf_mtop in
            zip(graph_xs.ys, inf_mtop_xs.ys)
            ]

        graph_xs.title = 'ct = {0:+.1f}, cg = {1:+.3f}, cb = {2:+.1f}'.format(ct, cg, cb)
        plot.add_top(graph_xs, 'repr_basic_line', plot.leg)
        plot.add_bottom(graph_ratio, 'repr_basic_line')

    plot.draw()
    plot.wrapup()



@flag_as_option
def test_param4(args):
    """Test the interference-term-neglected vs. full for a bunch of points """

    drawer = ParametrizationDrawer('ct', 'cg', 'cb')
    drawer.do_linear_terms = False
    drawer.ktcgkb_variations(include_sm=True)
    drawer.make_parametrization()

    drawer_nointerference = ParametrizationDrawer('ct', 'cg')
    drawer_nointerference.do_linear_terms = False
    drawer_nointerference.ktcg_variations(include_sm=True)
    drawer_nointerference.make_parametrization()

    # drawer_linterminterference = ParametrizationDrawer('ct', 'cg')
    # drawer_linterminterference.do_linear_terms = True
    # drawer_linterminterference.ktcg_variations()
    # drawer_linterminterference.make_parametrization()

    plot = differentials.plotting.plots.BottomPanelPlot('param_crosscheck_for_points')
    plot.make_legend()
    plot.leg.SetNColumns(2)
    plot.top_x_min    = drawer.sm['binCenters'][0]
    plot.bottom_x_min = drawer.sm['binCenters'][0]
    plot.top_x_max    = drawer.sm['binCenters'][-1]
    plot.bottom_x_max = drawer.sm['binCenters'][-1]

    plot.top_y_min = 0.00001
    plot.top_y_max = 1.1
    plot.bottom_y_min = 0.0
    plot.bottom_y_max = 2.0

    points = [
        (1.0, 0.0), # SM resolved loops
        (0.0, 1./12.), # SM point coupling
        (0.55,  0.044), # NW (long)
        (1.52, -0.024), # SE (long)
        (1.14,  0.017), # NE (short)
        (0.94,  0.004), # SW (short)
        ]
    for i_point, point in enumerate(points):
        ct, cg = point
        cb = 1.
        color = drawer.colors[i_point]

        # With interference
        graph_xs, graph_ratio = drawer.get_eval_param_graph(ct=ct, cg=cg, cb=cb, color=color)
        graph_xs.title = 'ct = {0:+.1f}, cg = {1:+.3f}, cb = {2:+.1f}'.format(ct, cg, cb)
        plot.add_top(graph_xs, 'repr_basic_line', plot.leg)
        plot.add_bottom(graph_ratio, 'repr_basic_line')

        # # With interference but only by adding linear terms
        # graph_xs, graph_ratio = drawer_linterminterference.get_eval_param_graph(ct=ct, cg=cg, color=color)
        # graph_xs.title = 'ct = {0:+.1f}, cg = {1:+.3f}, cb = {2:+.1f}'.format(ct, cg, cb)
        # plot.add_top(graph_xs, 'repr_markers')
        # plot.add_bottom(graph_ratio, 'repr_markers')

        # Without interference
        graph_xs, graph_ratio = drawer_nointerference.get_eval_param_graph(ct=ct, cg=cg, color=color)
        graph_xs.title = 'ct = {0:+.1f}, cg = {1:+.3f}, cb = {2:+.1f}'.format(ct, cg, cb)
        plot.add_top(graph_xs, 'repr_dashed_line')
        plot.add_bottom(graph_ratio, 'repr_dashed_line')

    plot.draw()
    plot.wrapup()



@flag_as_option
def test_get_real_cg_sm_value(args):
    """Vary cg until it matches the total XS from ct=cb=1"""

    drawer = ParametrizationDrawer('ct', 'cg', 'cb')
    drawer.do_linear_terms = False
    drawer.ktcgkb_variations(include_sm=True)
    drawer.make_parametrization()

    incl_xs = drawer.parametrization.incl_xs(ct=1.0, cg=0.0, cb=1.0)


    # Function to minimize
    def chi2(cg):
        return (incl_xs - drawer.parametrization.incl_xs(ct=0.0, cg=cg, cb=0.0))**2

    dx_derivative = 0.000001
    def dchi2_dx(cg):
        return (chi2(cg+dx_derivative) - chi2(cg-dx_derivative)) / (2.*dx_derivative)

    x = 0.1
    r = -0.00000001
    for i in xrange(1000):
        print 'Iter {0:4}: cg = {1:+7.4f}, chi2 = {2:+7.4f}'.format(
            i, x, chi2(x)
            )
        x += r * dchi2_dx(x)


    print 'Found: cg =', x





class ParametrizationDrawer(object):
    """docstring for ParametrizationDrawer"""
    def __init__(self, *coupling_names):
        super(ParametrizationDrawer, self).__init__()
        self.coupling_names = coupling_names
        self.colors = (
            [ differentials.core.safe_colors[kw] for kw in ['red','blue','green','lightblue'] ]
            + [ 5, 6, 9, 42, 46 ]
            )
        self.central_style = differentials.plotting.pywrappers.StyleSheet()
        self.do_linear_terms = False

    def print_sm(self):
        if not hasattr(self, 'sm_variation'):
            coupling_vals = { name : self.sm[name] for name in self.coupling_names }
            self.sm_variation = differentials.parametrization.Variation(
                self.sm.crosssection, **coupling_vals
                )
        print 'SM variation:', self.sm_variation
        
    def ktcg_variations(self, include_sm=False):
        coupling_variations = FileFinder(
            cb=1.0, muR=1.0, muF=1.0, Q=1.0, directory=LatestPaths.theory.top.filedir
            ).get()
        sm = [ v for v in coupling_variations if v.ct==1.0 and v.cg==0.0 ][0]
        coupling_variations.pop(coupling_variations.index(sm))
        self.sm = sm
        self.coupling_variations = coupling_variations
        self.bin_boundaries = self.sm['binBoundaries']
        self.bin_widths = [ right - left for left, right in zip(self.bin_boundaries[:-1], self.bin_boundaries[1:]) ]
        if include_sm:
            self.coupling_variations.insert(0, sm)
        self.print_sm()

    def ktcgkb_variations(self, include_sm=False):
        coupling_variations = FileFinder(
            muR=1.0, muF=1.0, Q=1.0, directory=LatestPaths.theory.top.filedir
            ).get()
        sm = [ v for v in coupling_variations if v.ct==1.0 and v.cg==0.0 and v.cb==1.0 ][0]
        coupling_variations.pop(coupling_variations.index(sm))
        coupling_variations = [ v for v in coupling_variations if v.ct!=1.0 and v.cg!=0.0 and v.cb!=1.0 ]
        self.sm = sm
        self.coupling_variations = coupling_variations
        self.bin_boundaries = self.sm['binBoundaries']
        self.bin_widths = [ right - left for left, right in zip(self.bin_boundaries[:-1], self.bin_boundaries[1:]) ]
        if include_sm:
            self.coupling_variations.insert(0, sm)
        self.print_sm()

    def ktkb_variations(self, include_sm=False):
        coupling_variations = FileFinder(
            cg=0.0, muR=1.0, muF=1.0, Q=1.0, directory=LatestPaths.theory.top.filedir
            ).get()
        sm = [ v for v in coupling_variations if v.ct==1.0 and v.cb==1.0 ][0]
        coupling_variations.pop(coupling_variations.index(sm))
        self.sm = sm
        self.coupling_variations = coupling_variations
        if include_sm:
            self.coupling_variations.insert(0, sm)
        self.print_sm()

    def kbkc_variations(self, include_sm=False):
        coupling_variations = FileFinder(
            muR=1.0, muF=1.0, Q=1.0, directory=LatestPaths.theory.yukawa.filedir
            ).get()

        sm = [ v for v in coupling_variations if v.kappab==1.0 and v.kappac==1.0 ][0]
        coupling_variations.pop(coupling_variations.index(sm))
        self.sm = sm
        self.coupling_variations = coupling_variations
        if include_sm:
            self.coupling_variations.insert(0, sm)
        self.print_sm()

    def make_parametrization_2D(self):
        self.parametrization = differentials.parametrization.Parametrization2Dim()
        self.parametrization.do_linear_terms = self.do_linear_terms
        self.parametrization.parametrize_by_matrix_inversion = True
        self.parametrization.c1_name = self.c1
        self.parametrization.c2_name = self.c2
        self.parametrization.c1_SM = 0.0 if self.c1 == 'cg' else 1.0
        self.parametrization.c2_SM = 0.0 if self.c2 == 'cg' else 1.0
        self.parametrization.from_theory_dicts(self.coupling_variations)

    def make_parametrization(self):
        self.parametrization = differentials.parametrization.ParametrizationMultiDim(len(self.coupling_names))
        self.parametrization.do_linear_terms = self.do_linear_terms
        self.parametrization.parametrize_by_matrix_inversion = True
        self.parametrization.set_coupling_names(*self.coupling_names)
        self.parametrization.from_theory_dicts(self.coupling_variations)
        self.parametrization.set_binning(self.sm.binBoundaries)

    def inclusive_xs_for_parametrization(self, c1, c2):
        crosssection_per_bin = self.parametrization.evaluate(c1, c2)
        incl = sum([ xs * width for xs, width in zip(crosssection_per_bin, self.bin_widths) ])
        return incl

    def variation_to_graph(self, variation, xs=True, color=None):
        title = []
        for c in [ 'ct', 'cg', 'cb', 'kappab', 'kappac', 'kappat' ]:
            if c in variation:
                title.append(
                    '{0}={1:.{ndec}f}'
                    .format(differentials.core.standard_titles[c], variation[c], ndec=3 if c=='cg' else 1 )
                    )
        title = ' '.join(title)

        y_attr = 'crosssection' if xs else 'ratios'
        graph = differentials.plotting.pywrappers.Graph(
            'auto', title, variation['binCenters'], variation[y_attr]
            )

        style = self.central_style.copy()
        if not(color is None): style.color = color
        graph.stylesheets.append(style)
        return graph

    def add_variation(self, variation, color):
        xs = self.variation_to_graph(variation, color=color)
        ratio = self.variation_to_graph(variation, xs=False, color=color)
        self.plot.add_top(xs, 'repr_markers', self.plot.leg)
        self.plot.add_bottom(ratio, 'repr_markers')

    # def get_eval_param_graph_for_variation(self, variation, color=None):
    #     c1 = variation[self.c1]
    #     c2 = variation[self.c2]
    #     return self.get_eval_param_graph(c1, c2, color=color)

    def get_eval_param_graph_for_variation(self, variation, color=None):
        coupling_vals = { name : variation[name] for name in self.coupling_names }
        return self.get_eval_param_graph(color=color, **coupling_vals)

    def get_eval_param_graph(self, color=None, **coupling_vals):
        style = self.central_style.copy()
        if not(color is None): style.color = color

        crosssection = self.parametrization.evaluate(**coupling_vals)
        ratios = [ xs_param / xs_sm for xs_param, xs_sm in zip(crosssection, self.sm.crosssection) ]

        bin_centers = self.sm['binCenters']

        graph_xs = differentials.plotting.pywrappers.Graph(
            'auto', 'notimportant',
            bin_centers,
            crosssection,
            )
        graph_xs.stylesheets.append(style)

        graph_ratio = differentials.plotting.pywrappers.Graph(
            'auto', 'notimportant',
            bin_centers,
            ratios,
            )
        graph_ratio.stylesheets.append(style)

        return graph_xs, graph_ratio


    def plot(self,
            plotname,
            top_y_min = 0.00001,
            top_y_max = 1.1,
            bottom_y_min = -0.5,
            bottom_y_max = 2.5,
            ):

        self.plot = differentials.plotting.plots.BottomPanelPlot(plotname)
        self.plot.make_legend()
        self.plot.top_x_min = self.sm['binCenters'][0]
        self.plot.top_x_max = self.sm['binCenters'][-1]
        self.plot.bottom_x_min = self.sm['binCenters'][0]
        self.plot.bottom_x_max = self.sm['binCenters'][-1]
        self.plot.top_y_min = top_y_min
        self.plot.top_y_max = top_y_max
        self.plot.bottom_y_min = bottom_y_min
        self.plot.bottom_y_max = bottom_y_max

        self.add_variation(self.sm, color=1)
        for i_variation, variation in enumerate(self.coupling_variations):
            color = self.colors[i_variation]
            self.add_variation(variation, color=color)
            
            graph_xs, graph_ratio = self.get_eval_param_graph_for_variation(variation, color=color)
            self.plot.add_top(graph_xs, 'repr_basic_line')
            self.plot.add_bottom(graph_ratio, 'repr_basic_line')


        self.plot.draw()
        self.plot.wrapup()




