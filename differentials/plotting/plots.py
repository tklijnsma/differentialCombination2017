import ROOT
import plotting_utils as utils
import pywrappers
from canvas import c

import differentials
import differentials.core as core

import logging, sys
from math import isnan, isinf, log10, sqrt
from array import array
from collections import namedtuple


class PlotBase(object):
    """docstring for PlotBase"""
    def __init__(self, plotname):
        self.plotname = plotname
        self.disable_CMS_labels = False
        c.Clear()
        c.set_margins()
        
    def pre_draw(self):
        pass

    def post_draw(self):
        pass

    def draw(self):
        self.pre_draw()

    def wrapup(self):
        self.post_draw()
        self.save()

    def save(self):
        c.save(self.plotname)



class QuickPlot(PlotBase):
    """docstring for QuickPlot"""
    def __init__(self, plotname, x_min=-1.0, x_max=4.0, y_min=0.0, y_max=5.0):
        super(QuickPlot, self).__init__(plotname)
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.leg = pywrappers.Legend()
        self.objs = []

        self.x_title = 'x'
        self.y_title = 'y'
        differentials.plotting.canvas.reset_global_color_cyle()

        self.do_legend = True

    def add(self, obj, draw_str):
        self.objs.append([obj, draw_str])

    def clear(self):
        self.objs = []
        self.leg.clear()

    def draw(self):
        c.Clear()
        c.set_margins()

        self.base = pywrappers.Base(
            x_min=self.x_min, x_max=self.x_max, y_min=self.y_min, y_max=self.y_max,
            x_title=self.x_title, y_title=self.y_title
            )
        self.base.Draw()

        for obj, draw_str in self.objs:
            obj.legend = self.leg
            obj.Draw(draw_str)

    def wrapup(self):
        if self.do_legend: self.leg.Draw()
        super(QuickPlot, self).wrapup()


class Single2DHistPlot(PlotBase):
    """docstring for Single2DHistPlot"""
    def __init__(
            self, plotname, H2, draw_str='repr_2D_with_contours',
            x_min=None, x_max=None, y_min=None, y_max=None
            ):
        super(Single2DHistPlot, self).__init__(plotname)
        self.H2 = H2
        self.draw_str = draw_str

        self.leg = pywrappers.Legend()
        self.x_title = 'x'
        self.y_title = 'y'
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        differentials.plotting.canvas.reset_global_color_cyle()
        self.set_ranges_by_contour = True

        self.x_SM = None
        self.y_SM = None

    def draw(self):
        c.Clear()
        c.set_margins_2D()

        leg = pywrappers.Legend(
            c.GetLeftMargin() + 0.01,
            c.GetBottomMargin() + 0.02,
            1 - c.GetRightMargin() - 0.01,
            c.GetBottomMargin() + 0.09
            )

        self.H2.legend = leg
        self.H2.Draw(self.draw_str)

        H = self.H2.H2
        H.GetXaxis().SetTitle(self.x_title)
        H.GetYaxis().SetTitle(self.y_title)
        H.GetXaxis().SetTitleSize(0.06)
        H.GetXaxis().SetLabelSize(0.05)
        H.GetYaxis().SetTitleSize(0.06)
        H.GetYaxis().SetLabelSize(0.05)

        if self.set_ranges_by_contour:
            extrema = self.H2.get_extrema_from_contour('2sigma')
            x_min_abs = extrema.x_min
            x_max_abs = extrema.x_max
            dx_abs = x_max_abs-x_min_abs
            y_min_abs = extrema.y_min
            y_max_abs = extrema.y_max
            dy_abs = y_max_abs-y_min_abs
            # Add some margin
            self.x_min = x_min_abs - 0.2*dx_abs
            self.x_max = x_max_abs + 0.2*dx_abs
            self.y_min = y_min_abs - 0.2*dy_abs
            self.y_max = y_max_abs + 0.3*dy_abs

        if not(self.x_min is None) and not(self.x_max is None):
            H.GetXaxis().SetRangeUser(self.x_min, self.x_max)
        elif not(self.x_min is None):
            H.GetXaxis().SetRangeUser(self.x_min, H.GetXaxis().GetXmax())
        elif not(self.x_max is None):
            H.GetXaxis().SetRangeUser(H.GetXaxis().GetXmin(), self.x_max)
        if not(self.y_min is None) and not(self.y_max is None):
            H.GetYaxis().SetRangeUser(self.y_min, self.y_max)
        elif not(self.y_min is None):
            H.GetYaxis().SetRangeUser(self.y_min, H.GetYaxis().GetXmax())
        elif not(self.y_max is None):
            H.GetYaxis().SetRangeUser(H.GetYaxis().GetXmin(), self.y_min)

        if not self.disable_CMS_labels:
            pywrappers.CMS_Latex_type().Draw()
            pywrappers.CMS_Latex_lumi().Draw()

        if self.x_SM is None: self.x_SM = self.H2.bestfit().x
        if self.y_SM is None: self.y_SM = self.H2.bestfit().y
        pywrappers.Point(self.x_SM, self.y_SM).Draw('repr_SM_point')

        self.H2.Draw('repr_bestfitpoint') # Redraw since the SM frequently overlaps the BF for Asimov

        pywrappers.ContourDummyLegend(
            c.GetLeftMargin() + 0.11,
            1. - c.GetTopMargin() - 0.1,
            1. - c.GetRightMargin() - 0.01,
            1. - c.GetTopMargin() - 0.01,
            ).Draw()

        c.Update()
        c.RedrawAxis()

        if self.set_ranges_by_contour:
            # To make sure that a second call to .draw() doesn't break stuff
            self.x_min = None
            self.x_max = None
            self.y_min = None
            self.y_max = None

    def wrapup(self):
        self.leg.Draw()
        super(Single2DHistPlot, self).wrapup()


class MultiScanPlot(PlotBase):
    """docstring for MultiScanPlot"""
    def __init__(self, plotname):
        super(MultiScanPlot, self).__init__(plotname)
        self.x_min = -1.0
        self.x_max = 4.0
        self.y_min = 0.0
        self.y_max = 5.0
        self.leg = pywrappers.Legend()
        self.x_title = 'POI'
        differentials.plotting.canvas.reset_global_color_cyle()

        self.y_cutoff = 3.15
        self.do_95percent_CL = False

        self.scans = []
        self.manual_graphs = []


    def add_scan(self, scan):
        self.scans.append(scan)

    def draw(self):
        c.Clear()
        c.set_margins()

        self.base = utils.get_plot_base(
            x_min=self.x_min, x_max=self.x_max, y_min=self.y_min, y_max=self.y_max,
            x_title=self.x_title, y_title=core.get_standard_title('dnll')
            )
        self.base.Draw('P')

        if self.do_95percent_CL:
            sigma1_line = ROOT.TLine(self.x_min, 3.841, self.x_max, 3.841)
        else:
            sigma1_line = ROOT.TLine(self.x_min, 1.0, self.x_max, 1.0)
        ROOT.SetOwnership(sigma1_line, False)
        sigma1_line.SetLineColor(14)
        sigma1_line.SetLineStyle(2)
        sigma1_line.Draw()

        # sigma2_line = ROOT.TLine(self.x_min, 2.0, self.x_max, 2.0)
        # ROOT.SetOwnership(sigma2_line, False)
        # sigma2_line.SetLineColor(14)
        # sigma2_line.SetLineStyle(2)
        # sigma2_line.Draw()

        for scan in self.scans:
            graph = scan.to_graph()
            graph.filter(y_max=self.y_cutoff)
            graph.legend = self.leg
            graph.Draw(getattr(graph, 'draw_style', 'repr_basic_line'))

            left_point, right_point = self.get_unc_points(scan, color=graph.color)

            if not hasattr(scan, 'no_bestfit_line'):
                line_bestfit = ROOT.TGraph(2,
                    array('f', [scan.bestfit().x, scan.bestfit().x]),
                    array('f', [0.0, 1.0]),
                    )
                ROOT.SetOwnership(line_bestfit, False)
                line_bestfit.SetLineColor(graph.color)
                if 'dashed' in getattr(graph, 'draw_style', 'repr_basic_line'):
                    line_bestfit.SetLineStyle(2)
                line_bestfit.Draw('SAMEL')


        for graph in self.manual_graphs:
            graph.flat_cut(y_top=self.y_cutoff)
            graph.legend = self.leg
            graph.Draw(getattr(graph, 'draw_style', 'repr_basic_line'))

            left_point, right_point = self.get_unc_points(graph, color=graph.style().color)

            if getattr(graph, 'draw_bestfit', False):
                x_bestfit = graph.unc.x_min
                line_bestfit = ROOT.TGraph(2,
                    array('f', [x_bestfit, x_bestfit]),
                    array('f', [0.0, self.y_cutoff]),
                    )
                ROOT.SetOwnership(line_bestfit, False)
                line_bestfit.SetLineColor(graph.style().color)
                if 'dashed' in getattr(graph, 'draw_style', 'repr_basic_line'):
                    line_bestfit.SetLineStyle(2)
                line_bestfit.SetLineWidth(2)
                line_bestfit.Draw('SAMEL')


        pywrappers.CMS_Latex_type().Draw()
        pywrappers.CMS_Latex_lumi().Draw()

        utils.format_plot_base(self.base)


    def get_unc_points(self, scan, color=1):
        """Should work in Graph as well"""

        point_height = 1.0
        if self.do_95percent_CL:
            point_height = 3.841

        if not hasattr(scan, 'unc'):
            scan.create_uncertainties(do_95percent_CL=self.do_95percent_CL)

        left_point = ROOT.TGraph(1, array('f', [scan.unc.left_bound]), array('f', [point_height]))
        ROOT.SetOwnership(left_point, False)
        left_point.SetMarkerColor(color)
        left_point.SetMarkerSize(1.1)
        left_point.SetMarkerStyle(8)
        if not scan.unc.well_defined_left_bound:
            left_point.SetMarkerStyle(5)
        left_point.Draw('PSAME')

        right_point = ROOT.TGraph(1, array('f', [scan.unc.right_bound]), array('f', [point_height]))
        ROOT.SetOwnership(right_point, False)
        right_point.SetMarkerColor(color)
        right_point.SetMarkerSize(1.1)
        right_point.SetMarkerStyle(8)
        if not scan.unc.well_defined_right_bound:
            right_point.SetMarkerStyle(5)
        right_point.Draw('PSAME')

        return left_point, right_point


    def wrapup(self):
        self.leg.Draw()
        super(MultiScanPlot, self).wrapup()



class MultiContourPlot(PlotBase):
    """docstring for MultiContourPlot"""

    z_axis_title = '-2#Delta ln L'

    def __init__(self, plotname, scans, x_min=None, x_max=None, y_min=None, y_max=None, x_title=None, y_title=None):
        super(MultiContourPlot, self).__init__(plotname)

        self.scans = scans

        self.x_min = x_min
        self.y_min = y_min
        self.x_max = x_max
        self.y_max = y_max

        self.draw_individual_contours = True
        self.only_1sigma_contours = False
        self.set_ranges_by_contour = False
        self.draw_bestfit_point = True

        self.use_first_scan_as_base = True
        self.base_is_set_manually = False

        self.only_1sigma_contours_for_secondary = True

        self.x_SM = 1.0
        self.y_SM = 1.0

        if x_title is None:
            self.x_title = getattr(scans[0], 'x_title', 'x')
        else:
            self.x_title = x_title

        if y_title is None:    
            self.y_title = getattr(scans[0], 'y_title', 'y')
        else:
            self.y_title = y_title

        self.base = pywrappers.Base(x_title=self.x_title, y_title=self.y_title)

        self.legend = pywrappers.Legend(
            lambda c: c.GetLeftMargin() + 0.07,
            lambda c: c.GetBottomMargin() + 0.02,
            lambda c: 1 - c.GetRightMargin() - 0.01,
            lambda c: c.GetBottomMargin() + 0.09
            )
        self.legend.set(
            lambda c: c.GetLeftMargin() + 0.02,
            lambda c: 1. - c.GetTopMargin() - 0.20,
            lambda c: c.GetLeftMargin() + 0.24,
            lambda c: 1. - c.GetTopMargin() - 0.01,
            )
        self.legend.SetBorderSize(1)
        self.legend.SetFillStyle(1001)
        self.legend.SetNColumns(1)

        self.cdl = pywrappers.ContourDummyLegend(
            lambda c: c.GetLeftMargin() + 0.07,
            lambda c: 1. - c.GetTopMargin() - 0.1,
            lambda c: 1. - c.GetRightMargin() - 0.01,
            lambda c: 1. - c.GetTopMargin() - 0.01,
            )
        # self.cdl.set(
        #     lambda c: c.GetLeftMargin() + 0.02,
        #     lambda c: 1. - c.GetTopMargin() - 0.22,
        #     lambda c: c.GetLeftMargin() + 0.16,
        #     lambda c: 1. - c.GetTopMargin() - 0.01,
        #     )
        self.cdl.set(
            lambda c: c.GetLeftMargin() + 0.03,
            lambda c: c.GetBottomMargin() + 0.03,
            lambda c: c.GetLeftMargin() + 0.49,
            lambda c: c.GetBottomMargin() + 0.08,
            )
        self.cdl.SetBorderSize(1)
        self.cdl.SetFillStyle(1001)
        # self.cdl.SetNColumns(1)

        self.texts = []


    # def pre_draw(self):
    #     super(MultiContourPlot, self).pre_draw()
    #     self.base.Draw()

    def set_base(self, hist, draw_str='COLZ'):
        if isinstance(hist, differentials.plotting.pywrappers.Histogram2D):
            self.base, draw_str = hist.repr_2D()[0]
            self.base_draw_str = draw_str
        else:
            self.base = hist
            self.base_draw_str = draw_str
        self.use_first_scan_as_base = False
        self.base_is_set_manually = True

    def extrema_from_contours(self, H2s):
        x_min = 10e9
        x_max = -10e9
        y_min = 10e9
        y_max = -10e9
        for H2 in H2s:
            extrema = H2.get_extrema_from_contour('2sigma')
            if extrema.x_min < x_min: x_min = extrema.x_min
            if extrema.x_max > x_max: x_max = extrema.x_max
            if extrema.y_min < y_min: y_min = extrema.y_min
            if extrema.y_max > y_max: y_max = extrema.y_max
        self.x_min = x_min - 0.4*(x_max-x_min)
        self.x_max = x_max + 0.4*(x_max-x_min)
        self.y_min = y_min - 0.3*(y_max-y_min)
        self.y_max = y_max + 0.4*(y_max-y_min)


    def add_BR_parametrized_text(self, x=None, y=None, coupling_x='kappac', coupling_y='kappab'):
        text = '{0}({1}, {2})'.format(
            core.standard_titles['BR'], core.standard_titles[coupling_x], core.standard_titles[coupling_y]
            )
        self.add_BR_text(x, y, text)

    def add_BR_floating_text(self, x=None, y=None, coupling_x='kappac', coupling_y='kappab'):
        text = '{0} unconstr.'.format(
            core.standard_titles['BR'], core.standard_titles[coupling_x], core.standard_titles[coupling_y]
            )
        self.add_BR_text(x, y, text)

    def add_BR_text(self, x, y, text):
        if x is None: x = lambda c: 1. - c.GetRightMargin() - 0.01
        if y is None: y = lambda c: c.GetBottomMargin() + 0.03 + 0.015
        self.br_text = pywrappers.Latex(
            x(c), y(c), text
            )
        self.br_text.SetTextSize(0.045)
        if x(c) > 0.5:
            self.br_text.SetTextAlign(31)
        else:
            self.br_text.SetTextAlign(11)
        self.br_text.SetNDC()
        self.br_text.Draw()

    def add_text_on_the_fly(self,
            x, y, text,
            text_size = 0.045,
            text_align = 'auto',
            ndc = True,
            ):
        if callable(x): x = x(c)
        if callable(y): y = y(c)

        latex = pywrappers.Latex(x, y, text)
        self.texts.append(latex)

        if ndc: latex.SetNDC()
        latex.SetTextSize(text_size)
        if text_align == 'auto':
            if x > 0.5:
                latex.SetTextAlign(31)
            else:
                latex.SetTextAlign(11)
        else:
            latex.SetTextAlign(text_align)

        latex.Draw()
        return latex

    def set_ranges(self):
        if self.set_ranges_by_contour:
            self.extrema_from_contours(self.histograms)
        else:
            if self.x_min is None: self.x_min = min([H.x_min() for H in self.histograms])
            if self.y_min is None: self.y_min = min([H.y_min() for H in self.histograms])
            if self.x_max is None: self.x_max = max([H.x_max() for H in self.histograms])
            if self.y_max is None: self.y_max = max([H.y_max() for H in self.histograms])

    def draw(self, wait=False):
        super(MultiContourPlot, self).draw()
        c.resize_temporarily(870, 800)
        c.SetRightMargin(0.13)

        if self.use_first_scan_as_base:
            if isinstance(self.scans[0], differentials.plotting.pywrappers.Histogram2D):
                self.base_hist = self.scans[0]
            else:
                self.base_hist = self.scans[0].to_hist()
            self.base, draw_str = self.base_hist.repr_2D()[0]
            self.base.SetName('base_' + draw_str)
            self.base.Draw(draw_str)
        elif self.base_is_set_manually:
            self.base.Draw(self.base_draw_str)
        else:
            self.base.Draw()

        self.base.GetZaxis().SetTitle(self.z_axis_title)
        self.base.GetZaxis().SetTitleSize(0.04)
        self.base.GetZaxis().SetTitleOffset(0.9)

        self.histograms = []
        for scan in self.scans[::-1]:
            if isinstance(scan, differentials.plotting.pywrappers.Histogram2D):
                histogram2D = scan
            else:
                histogram2D = scan.to_hist()
            histogram2D.legend = self.legend
            if self.only_1sigma_contours or getattr(scan, 'only_1sigma_contours', False):
                histogram2D.Draw('repr_1sigma_contours')
            elif self.only_1sigma_contours_for_secondary and scan != self.scans[0]:
                histogram2D.Draw('repr_1sigma_contours')
            else:
                histogram2D.Draw('repr_contours_no_bestfit')
            self.histograms.append(histogram2D)

        self.legend.Draw()

        self.set_ranges()
        
        if self.use_first_scan_as_base or self.base_is_set_manually:
            self.base.GetXaxis().SetRangeUser(self.x_min, self.x_max)
            self.base.GetYaxis().SetRangeUser(self.y_min, self.y_max)
        else:
            self.base.set(x_min=self.x_min, y_min=self.y_min, x_max=self.x_max, y_max=self.y_max)

        self.base.GetXaxis().SetTitle(self.x_title)
        self.base.GetYaxis().SetTitle(self.y_title)
        self.base.GetXaxis().SetTitleSize(0.065)
        self.base.GetXaxis().SetLabelSize(0.05)
        self.base.GetYaxis().SetTitleSize(0.065)
        self.base.GetYaxis().SetLabelSize(0.05)

        SM_point = pywrappers.Point(self.x_SM, self.y_SM)
        SM_point.Draw('repr_SM_point')

        if self.draw_bestfit_point:
            for histogram in self.histograms:
                histogram.Draw('repr_bestfitpoint') # Redraw since the SM frequently overlaps the BF for Asimov

        if not self.disable_CMS_labels:
            pywrappers.CMS_Latex_type().Draw()
            pywrappers.CMS_Latex_lumi().Draw()

        if not(self.draw_bestfit_point): self.cdl.disable_bestfit = True
        self.cdl.Draw()

        if not(wait): self.wrapup()

    def see_through_legends(self):
        self.legend.SetFillStyle(0)
        self.legend.SetBorderSize(0)
        self.cdl.SetFillStyle(0)
        self.cdl.SetBorderSize(0)

    def wrapup(self):
        c.Update()
        c.RedrawAxis()

        l = c.GetLeftMargin()
        r = 1. - c.GetRightMargin()
        b = c.GetBottomMargin()
        t = 1. - c.GetTopMargin()

        top_line = ROOT.TLine(l, t, r, t)
        top_line.SetNDC()
        ROOT.SetOwnership(top_line, False)
        top_line.Draw()
        
        right_line = ROOT.TLine(r, b, r, t)
        right_line.SetNDC()
        ROOT.SetOwnership(right_line, False)
        right_line.Draw()

        self.save()        
        if self.draw_individual_contours:
            if self.draw_bestfit_point:
                draw_style = 'repr_2D_with_contours'
            else:
                draw_style = 'repr_2D_with_contours_no_bestfit'
            for scan in self.scans:
                plotname = self.plotname + '_' + scan.name
                if isinstance(scan, differentials.scans.Scan2D):
                    scan.x_sm = self.x_SM
                    scan.y_sm = self.y_SM
                    scan.plot(plotname, draw_style=draw_style)
                elif isinstance(scan, differentials.plotting.pywrappers.Histogram2D):
                    scan.quickplot(
                        plotname, draw_method=draw_style,
                        x_sm = self.x_SM, y_sm = self.y_SM,
                        x_title = self.x_title, y_title = self.y_title
                        )
                else:
                    logging.error('Object {0} could not be plotted individually'.format(scan))


class BottomPanelPlot(PlotBase):
    """docstring for BottomPanelPlot"""
    def __init__(self, plotname):
        super(BottomPanelPlot, self).__init__(plotname)

        self.top_objects = []
        self.bottom_objects = []

        self.x_title = ''
        self.y_title_top = ''
        self.y_title_bottom = ''

        self.pad_split_point = 0.33

        self.debug_lines = False

        self.top_bottom_margin = 0.02
        self.top_top_margin    = 0.10
        self.top_left_margin   = 0.19
        self.top_right_margin  = 0.03

        # self.bottom_bottom_margin = 0.30
        self.bottom_bottom_margin = 0.35
        self.bottom_top_margin    = 0.00
        # self.bottom_left_margin   = 0.14
        self.bottom_left_margin   = 0.19
        self.bottom_right_margin  = 0.03

        self.top_log_scale = True
        self.CMS_labels = True

        self.top_x_min = None
        self.top_x_max = None
        self.top_y_min = None
        self.top_y_max = None
        self.bottom_x_min = None
        self.bottom_x_max = None
        self.bottom_y_min = None
        self.bottom_y_max = None
        self.has_legend = False

        self.fixed_widths = False
        self.add_raster_bottom = False
        self.lumi_text_size = 0.07


    def make_fixed_widths(self, reference_bounds):
        self.fixed_widths = True
        self.reference_bounds = reference_bounds

    def make_legend(self):
        self.has_legend = True
        self.leg = pywrappers.Legend(
            lambda c: c.GetLeftMargin() + 0.01,
            lambda c: 1 - c.GetTopMargin() - 0.10,
            lambda c: 1 - c.GetRightMargin() - 0.01,
            lambda c: 1 - c.GetTopMargin()
            )

    def add_top(self, obj, draw_str=None, leg=None):
        if not(leg is None):
            obj.legend = leg
        self.top_objects.append((obj, draw_str))

    def add_bottom(self, obj, draw_str=None, leg=None):
        if not(leg is None):
            obj.legend = leg
        self.bottom_objects.append((obj, draw_str))

    def get_obj_extremum(self, obj, extremum_fn, only_positive=False):
        logging.trace('Gettting extremum \'{0}\' for {2}; only_positive={1}'.format(extremum_fn, only_positive, obj))
        if hasattr(obj, extremum_fn):
            logging.trace('Method \'{0}()\' found'.format(extremum_fn))
            if only_positive:
                logging.trace('Trying to pass \'only_positive=True\'...')
                try:
                    extremum = getattr(obj, extremum_fn)(only_positive=True)
                    logging.trace('Trying to pass \'only_positive=True\' succeeded')
                except TypeError:
                    logging.trace('Trying to pass \'only_positive=True\' failed; calling without')
                    extremum = getattr(obj, extremum_fn)()
                    if extremum < 0.0:
                        logging.trace(
                            'Found {0}={1}; Overwriting with None as positive is required'
                            .format(extremum_fn, extremum)
                            )
                        extremum = None
            else:
                extremum = getattr(obj, extremum_fn)()
        else:
            logging.trace('No method \'{0}()\' found'.format(extremum_fn))
            extremum = None
        logging.debug('Found {0}={1}'.format(extremum_fn, extremum))
        return extremum

    def get_extremum(self, extremum_fn, default_value, objs, possible_overwrite_key='non_existing999', minimum=True, only_positive=False):
        logging.debug('Getting {0}'.format(extremum_fn))
        if getattr(self, possible_overwrite_key, None) is None:
            logging.debug('No pre-supplied value detected; looping over objects')
            extrema = []
            for obj in objs:
                extremum = self.get_obj_extremum(obj, extremum_fn, only_positive)
                if extremum is None: continue
                extrema.append(extremum)
            if len(extrema) == 0:
                logging.debug('Found no valid object extrema; passing default {0}'.format(default_value))
                extremum = default_value
            else:
                extremum = (min if minimum else max)(extrema)
        else:
            extremum = getattr(self, possible_overwrite_key)
            logging.debug('Using self.{2}: {0}={1}'.format(extremum_fn, extremum, possible_overwrite_key))
        logging.debug('Found overal {0}={1}'.format(extremum_fn, extremum))
        return extremum

    def get_top_extrema(self):
        logging.debug('Getting extrema for top objects')
        objs = [obj for obj, _ in self.top_objects]
        x_min = self.get_extremum('x_min', 0.001, objs, possible_overwrite_key='top_x_min', minimum=True)
        x_max = self.get_extremum('x_max', 1.000, objs, possible_overwrite_key='top_x_max', minimum=False)

        add_margin = False
        if self.top_y_min is None and self.top_y_max is None:
            add_margin = True

        y_min_abs = self.get_extremum('y_min', 0.001, objs, possible_overwrite_key='top_y_min', minimum=True, only_positive=True)
        y_max_abs = self.get_extremum('y_max', 1.000, objs, possible_overwrite_key='top_y_max', minimum=False, only_positive=True)

        if add_margin:
            # Add some margins
            if self.top_log_scale:
                y_max = y_max_abs + (y_max_abs/y_min_abs)**0.2 * y_max_abs
                y_min = 0.5*y_min_abs
            else:
                dy = y_max_abs - y_min_abs
                y_min -= 0.1*dy
                y_max += 0.1*dy
        else:
            y_min = y_min_abs
            y_max = y_max_abs

        return x_min, y_min, x_max, y_max

    def get_bottom_extrema(self):
        logging.debug('Getting extrema for bottom objects')
        objs = [obj for obj, _ in self.bottom_objects]

        add_margin = False
        if self.bottom_y_min is None and self.bottom_y_max is None:
            add_margin = True

        x_min = self.get_extremum('x_min', 0.001, objs,  possible_overwrite_key='bottom_x_min', minimum=True)
        x_max = self.get_extremum('x_max', 1.000, objs,  possible_overwrite_key='bottom_x_max', minimum=False)
        y_min = self.get_extremum('y_min', 0.001, objs,  possible_overwrite_key='bottom_y_min', minimum=True, only_positive=False)
        y_max = self.get_extremum('y_max', 1.000, objs,  possible_overwrite_key='bottom_y_max', minimum=False, only_positive=False)

        if add_margin:
            dy = y_max - y_min
            y_min -= 0.1*dy
            y_max += 0.1*dy

        return x_min, y_min, x_max, y_max


    def sort_objs_priority(self, objs):
        def get_priority(obj_str_tuple):
            obj = obj_str_tuple[0]
            if not hasattr(obj, 'style'):
                return 10
            else:
                return obj.style().plot_priority
        objs.sort(key=get_priority) # Highest priority should be last in list (drawn list = on top)

    def add_lines_at_bin_boundaries(self, bounds=None, stylesheet=None):
        self.add_raster_bottom = True
        if bounds is None:
            bounds = self.spectra[0].binning()[1:-1]
        self.add_raster_bottom_bounds = bounds
        if stylesheet is None:
            stylesheet = differentials.plotting.pywrappers.StyleSheet(line_width=1, color=17)
        self.add_raster_bottom_style = stylesheet

    def draw_lines_at_bin_boundaries(self, bounds, y_min, y_max, stylesheet):
        for bound in bounds:
            line = ROOT.TLine(bound, y_min, bound, y_max)
            ROOT.SetOwnership(line, False)
            stylesheet.apply(line)
            line.Draw()

    def draw(self):
        super(BottomPanelPlot, self).draw()
        c.cd()
        c.Clear()

        self._tmp_width = c.GetWindowWidth()
        self._tmp_height = c.GetWindowHeight()
        c.SetCanvasSize( 800, 900 )

        toppad_bottom    = self.pad_split_point
        bottompad_top    = self.pad_split_point
        height_ratio = float( 1.0 - toppad_bottom ) / float( bottompad_top )

        toppad = ROOT.TPad(
            utils.get_unique_rootname(), '',
            0.0, toppad_bottom, 1.0, 1.0
            )
        ROOT.SetOwnership(toppad, False)
        toppad.SetBottomMargin(self.top_bottom_margin) # Distance to the bottom panel
        toppad.SetTopMargin(self.top_top_margin)     # Space for labels
        toppad.SetLeftMargin(self.top_left_margin)
        toppad.SetRightMargin(self.top_right_margin)
        toppad.Draw()

        bottompad = ROOT.TPad(
            utils.get_unique_rootname(), '',
            0.0, 0.0, 1.0, bottompad_top
            )
        ROOT.SetOwnership(bottompad, False)
        bottompad.SetBottomMargin(self.bottom_bottom_margin)
        bottompad.SetTopMargin(self.bottom_top_margin)
        bottompad.SetLeftMargin(self.bottom_left_margin)
        bottompad.SetRightMargin(self.bottom_right_margin)
        bottompad.Draw()

        # Draw some lines helpful for debugging
        if self.debug_lines:
            for pad in [ toppad, bottompad ]:
                pad.cd()
                for x1, y1, x2, y2 in [
                        ( 0.0, 0.0, 1.0, 0.0 ),
                        ( 1.0, 0.0, 1.0, 1.0 ),
                        ( 0.0, 0.0, 0.0, 1.0 ),
                        ( 0.0, 1.0, 1.0, 1.0 ),
                        ]:
                    line = ROOT.TLine( x1, y1, x2, y2 )
                    ROOT.SetOwnership( line, False )
                    line.Draw()

            c.cd()
            line = ROOT.TLine( 0.0, bottompad_top, 1.0, toppad_bottom )
            ROOT.SetOwnership( line, False )
            line.Draw()

        # Create bases for top and bottom
        toppad.cd()
        top_x_min, top_y_min, top_x_max, top_y_max = self.get_top_extrema()

        base_top = utils.get_plot_base(
            x_min=top_x_min, x_max=top_x_max, y_min=top_y_min, y_max=top_y_max,
            x_title='Observable', y_title='#sigma (unit)'
            )
        base_top.Draw('P')

        bottompad.cd()
        bottom_x_min, bottom_y_min, bottom_x_max, bottom_y_max = self.get_bottom_extrema()

        # Add some margins
        dy = bottom_y_max - bottom_y_min
        bottom_y_min -= 0.1*dy
        bottom_y_max += 0.1*dy

        base_bottom = utils.get_plot_base(
            x_min=bottom_x_min, x_max=bottom_x_max, y_min=bottom_y_min, y_max=bottom_y_max,
            x_title='Observable', y_title='#mu'
            )
        base_bottom.Draw('P')


        #____________________________________________________________________
        # Draw the actual objects

        toppad.cd()
        if self.top_log_scale:
            toppad.SetLogy()

        self.sort_objs_priority(self.top_objects)
        self.sort_objs_priority(self.bottom_objects)

        logging.debug('top_objects: {0}'.format(self.top_objects))
        for obj, drawStr in self.top_objects:
            logging.debug('Attempting to draw obj: {0}, drawStr: {1}'.format(obj, drawStr))

            if self.fixed_widths and hasattr(obj, 'map_binning_to_fixed_binwidth'):
                obj.map_binning_to_fixed_binwidth(self.reference_bounds)

            if drawStr is None:
                obj.Draw()
            else:
                obj.Draw(drawStr)

        bottompad.cd()

        # Draw raster lines first if necessary:
        if self.add_raster_bottom:
            self.draw_lines_at_bin_boundaries(
                self.add_raster_bottom_bounds, bottom_y_min, bottom_y_max, self.add_raster_bottom_style
                )

        logging.debug('bottom_objects: {0}'.format(self.bottom_objects))
        for obj, drawStr in self.bottom_objects:
            logging.debug('Attempting to draw obj: {0}, drawStr: {1}'.format(obj, drawStr))

            if self.fixed_widths and hasattr(obj, 'map_binning_to_fixed_binwidth'):
                obj.map_binning_to_fixed_binwidth(self.reference_bounds)

            if drawStr is None:
                obj.Draw()
            else:
                obj.Draw(drawStr)

        #____________________________________________________________________
        # Process titles, ticks, labels
        toppad.cd()
        # base_top = self.top_objects[0][0]
        base_top.GetXaxis().SetLabelOffset(999.)
        base_top.GetXaxis().SetLabelSize(0.05)
        base_top.GetYaxis().SetLabelSize(0.05)

        base_top.GetYaxis().SetTitle(self.y_title_top)
        base_top.GetYaxis().SetTitleSize(0.08)
        base_top.GetYaxis().SetTitleOffset(0.9)

        bottompad.cd()
        # base_bottom = self.bottom_objects[0][0]
        # Set sizes of labels/ticks equal to top panel (undo automatic scaling by ROOT)
        base_bottom.GetXaxis().SetLabelSize(base_top.GetXaxis().GetLabelSize() * height_ratio)
        base_bottom.GetYaxis().SetLabelSize(base_top.GetYaxis().GetLabelSize() * height_ratio)
        base_bottom.GetXaxis().SetTickLength(base_top.GetXaxis().GetTickLength() * height_ratio)

        base_bottom.GetYaxis().SetTitle(self.y_title_bottom)
        base_bottom.GetXaxis().SetTitle(self.x_title)
        # base_bottom.GetXaxis().SetTitleSize(base_top.GetXaxis().GetTitleSize() * height_ratio)
        base_bottom.GetXaxis().SetTitleSize(0.075 * height_ratio)
        base_bottom.GetYaxis().SetTitleSize(base_top.GetYaxis().GetTitleSize() * height_ratio)

        # FIX for the longer "ratio to prediction" text
        # Title slightly smaller and fix offset a bit
        base_bottom.GetYaxis().SetTitleSize(0.05 * height_ratio)
        base_bottom.GetYaxis().SetTitleOffset(1.2 * 1./height_ratio)


        if not self.disable_CMS_labels:
            toppad.cd()
            pywrappers.CMS_Latex_type(text_size=0.08).Draw()
            pywrappers.CMS_Latex_lumi(text_size=self.lumi_text_size).Draw()


        self.height_ratio = height_ratio
        self.base_top = base_top
        self.base_bottom = base_bottom
        self.bottompad = bottompad
        self.toppad = toppad
        toppad.cd()


    def replace_bin_labels(self, new_labels, offset=0.0):
        self.bottompad.cd()

        # 'Delete' original labels
        self.base_bottom.GetXaxis().SetLabelOffset(999.)

        # text_size = 0.035*self.height_ratio
        text_size = self.base_top.GetYaxis().GetLabelSize() * self.height_ratio

        y = lambda c: c.GetBottomMargin() - 0.01*self.height_ratio

        x_min = self.base_bottom.GetXaxis().GetXmin()
        x_max = self.base_bottom.GetXaxis().GetXmax()
        dx = x_max - x_min
        left_margin = self.bottompad.GetLeftMargin()
        right_margin = self.bottompad.GetRightMargin()
        x_to_NDC = lambda x: left_margin + ((x/dx+offset) * (1.-left_margin-right_margin))

        for i, lbl in enumerate(new_labels):
            x = x_to_NDC(i)
            l = differentials.plotting.pywrappers.Latex(
                x, y,
                lbl
                )
            l.SetNDC()
            l.SetTextSize(text_size)
            l.SetTextAlign(23)
            l.SetTextFont(42)
            l.Draw()


    def wrapup(self):
        if self.has_legend:
            self.toppad.cd()
            self.leg.Draw()
        super(BottomPanelPlot, self).save()
        # c.SetCanvasSize(self._tmp_width, self._tmp_height)



class SpectraPlot(BottomPanelPlot):

    def __init__(self, plotname, spectra, x_min=None, x_max=None, y_min=None, y_max=None):
        super(SpectraPlot, self).__init__(plotname)

        self.spectra = spectra
        self.plotname = plotname
        # self.plot = differentials.plotting.multipanel.BottomPanelPlot(plotname)
        # self.y_title_bottom = '#mu'
        self.y_title_bottom = 'Data / prediction'

        self.obsname = 'obs_name'
        self.obsunit = None
        self.obstitle = None

        self.draw_multiscans = False
        self.draw_multiscans_all_in_one_plot = False

        self.scans_x_min = -10.0
        self.scans_x_max = 10.0
        self.scans_y_min = 0.
        self.scans_y_max = 5.0

        self.fix_x_max_for_spectra = True

        self.overflow_label_base_offset = 0.5

        self.leg = pywrappers.Legend(
            lambda c: c.GetLeftMargin() + 0.01,
            lambda c: 1 - c.GetTopMargin() - 0.10,
            lambda c: 1 - c.GetRightMargin() - 0.01,
            lambda c: 1 - c.GetTopMargin()
            )

        self.has_reference_bounds = False
        self.reference_bounds = None

    def set_reference_bounds(self, bounds):
        self.has_reference_bounds = True
        self.reference_bounds = bounds

    def make_SM_line(self, spectra_original, leg=None):
        spectra = spectra_original[:]
        spectra.sort(key=lambda s: -len(s.scans))
        smxs_histogram = spectra[0].to_hist_smxs()
        self.add_top(smxs_histogram, 'repr_basic_histogram', leg)
        sm_histogram = spectra[0].to_hist_sm()
        self.add_bottom(sm_histogram, 'repr_basic_histogram')

    def make_labels_for_overflow_spectra(self, spectra, obs_name):
        x_min, y_min, x_max, y_max = self.get_top_extrema()
        text_size = 0.04
        
        # Also take second to last bin width as the text overlaps
        y_overflow = max(
            max([ s.scans[-1].unc.right_bound * s.smxs[-1] for s in spectra ]),
            max([ s.scans[-2].unc.right_bound * s.smxs[-2] for s in spectra ])
            )
        
        i_label = 0
        for spectrum in spectra:
            if getattr(spectrum, 'no_overflow_label', False): continue

            y_overflow_NDC = lambda c, i_lambda=i_label: (
                self.overflow_label_base_offset + i_lambda * 2.0*text_size
                )
            # y_overflow_NDC = lambda c, i_lambda=i_label: (
            #     c.GetBottomMargin()
            #     + log10(y_overflow/y_min)/log10(y_max/y_min) * ( 1.0 - c.GetTopMargin() - c.GetBottomMargin() )
            #     + i_lambda * 2.0*text_size
            #     + 1.5*text_size
            #     )

            binning = spectrum.binning()
            scale = binning[-2] - binning[-3]
            scale_text = str(int(scale)) if float(scale).is_integer() else '{0:.2f}'.format(scale)
            left_bin = binning[-2]
            left_bin_text = str(int(left_bin)) if float(left_bin).is_integer() else '{0:.2f}'.format(left_bin)

            l = pywrappers.Latex(
                    lambda c: 1. - c.GetRightMargin() - 0.01,
                    y_overflow_NDC,
                    '#Delta#sigma({1} > {0}) / {2}'.format(
                        left_bin_text,
                        obs_name,
                        scale_text
                        )
                    )

            l.SetNDC()
            l.SetTextSize(text_size)
            l.SetTextColor(spectrum.style().color)
            l.SetTextAlign(31)
            # l.SetTextAlign(33)
            l.SetTextFont(42)
            self.add_top(l, '')
            i_label += 1


    def draw_multiscans_fn(self):
        if self.draw_multiscans_all_in_one_plot or self.draw_multiscans:
            differentials.plotting.canvas.reset_global_color_cyle()
            differentials.plotting.pywrappers.Graph.color_cycle = differentials.plotting.canvas.global_color_cycle

        if self.draw_multiscans_all_in_one_plot:
            plot = MultiScanPlot(self.plotname + '_scans_all')
            for spectrum in self.spectra:
                plot.scans.extend(spectrum.scans)
            plot.x_min = self.scans_x_min
            plot.x_max = self.scans_x_max
            plot.y_min = self.scans_y_min
            plot.y_max = self.scans_y_max
            plot.draw()
            plot.wrapup()
        elif self.draw_multiscans:
            for spectrum in self.spectra:
                spectrum.scans_x_min = self.scans_x_min
                spectrum.scans_x_max = self.scans_x_max
                spectrum.scans_y_min = self.scans_y_min
                spectrum.scans_y_max = self.scans_y_max
                spectrum.plot_scans(self.plotname + '_scans_' + spectrum.name)


    def draw(self):
        self.draw_multiscans_fn()

        if self.obstitle is None:
            self.obstitle = core.standard_titles.get(self.obsname, self.obsname)

        self.x_title = self.obstitle 
        if self.obsunit:
            self.x_title += ' ({0})'.format(self.obsunit)

        if self.y_title_top == '':
            self.y_title_top = '#Delta#sigma/#Delta{0} (pb{1})'.format(
                self.obstitle,
                '/' + self.obsunit if not(self.obsunit is None) else ''
                )

        for spectrum in self.spectra:
            if not spectrum.smxs_set:
                raise RuntimeError('SM cross sections were not provided for spectrum')
            if not spectrum._is_read:
                spectrum.read()

        if self.spectra[0].last_bin_is_overflow() and self.fix_x_max_for_spectra:
            # Make sure overflow bins are aligned
            x_max = max([ 2*s.binning()[-2]-s.binning()[-3] for s in self.spectra ])
            for spectrum in self.spectra:
                spectrum.hard_x_max = x_max

        for spectrum in self.spectra:
            hist_ratio = spectrum.to_hist()
            hist_xs = spectrum.to_hist_xs()
            if self.has_reference_bounds:
                hist_ratio.get_merged_bins_from_reference(self.reference_bounds)
                hist_xs.get_merged_bins_from_reference(self.reference_bounds)
            self.add_bottom(hist_ratio, spectrum.draw_method)
            self.add_top(hist_xs, spectrum.draw_method, leg=self.leg)

        self.make_labels_for_overflow_spectra(self.spectra, self.obstitle)
        self.add_top(self.leg, '')
        super(SpectraPlot, self).draw()



class BottomPanelPlotWithParametrizations(BottomPanelPlot):
    """docstring for BottomPanelPlotWithParametrizations"""
    def __init__(self, *args, **kwargs):
        super(BottomPanelPlotWithParametrizations, self).__init__(*args, **kwargs)
        self.scan2D = None
        self.ptspectrum = None
        self.ws_file = None
        self.obs = None
        self.points = []
        self.parametrized_graphs = []
        self.color_cycle = utils.new_color_cycle()
        self.last_bin_is_overflow = False

        self.x_SM = 1.0
        self.y_SM = 1.0

        self.default_points_xy_maxima = True

        self.get_points_method = 'extrema'

        self.legend = pywrappers.Legend(
            lambda c: c.GetLeftMargin()+0.01,
            lambda c: 1.-c.GetTopMargin()-0.40,
            lambda c: c.GetLeftMargin()+0.42,
            lambda c: 1.-c.GetTopMargin()-0.01,
            )
        self.legend.SetNColumns(1)
        # self.legend.SetBorderSize(1)

    def get_points(self):
        contour = self.scan2D.to_hist().get_most_probable_1sigma_contour()
        if self.get_points_method == 'extrema':
            bestfit_point = pywrappers.Point(self.scan2D.bestfit().x, self.scan2D.bestfit().y, color=self.color_cycle.next())
            bestfit_point._is_bestfit = True
            self.points.append(bestfit_point)
            self.points.append(pywrappers.Point(contour.x_min, contour.y[contour.x.index(contour.x_min)], color=self.color_cycle.next()))
            self.points.append(pywrappers.Point(contour.x_max, contour.y[contour.x.index(contour.x_max)], color=self.color_cycle.next()))
            self.points.append(pywrappers.Point(contour.x[contour.y.index(contour.y_min)], contour.y_min, color=self.color_cycle.next()))
            self.points.append(pywrappers.Point(contour.x[contour.y.index(contour.y_max)], contour.y_max, color=self.color_cycle.next()))

        elif self.get_points_method == 'max_distance_from_bestfit':
            Point = namedtuple('Point', ['x', 'y'])

            bestfit = Point(self.scan2D.bestfit().x, self.scan2D.bestfit().y)
            points = [ Point(x, y) for x, y in zip(contour.x, contour.y) ]

            # Function that returns distance between two points in x-y plane
            d = lambda p1, p2: sqrt( (p1.x-p2.x)**2 + (p1.y-p2.y)**2 )

            # Sort according to decreasing distance to bestfit
            points.sort(key=lambda p: -d(p, bestfit))
            # First selected point is point furthest away from bestfit
            p1 = points.pop(0)

            # Sort according to decreasing distance to first selected point
            points.sort(key=lambda p: -d(p, p1))
            p2 = points.pop(0)

            p3 = Point(
                x = bestfit.x+0.5*(p1.x-bestfit.x),
                y = bestfit.y+0.5*(p1.y-bestfit.y),
                )

            self.points.append(pywrappers.Point(bestfit.x, bestfit.y, color=self.color_cycle.next()))
            self.points.append(pywrappers.Point(p1.x, p1.y, color=self.color_cycle.next()))
            self.points.append(pywrappers.Point(p2.x, p2.y, color=self.color_cycle.next()))
            self.points.append(pywrappers.Point(p3.x, p3.y, color=self.color_cycle.next()))

        elif self.get_points_method == 'extrema_y_and_x_max':
            bestfit_point = pywrappers.Point(self.scan2D.bestfit().x, self.scan2D.bestfit().y, color=self.color_cycle.next())
            bestfit_point._is_bestfit = True
            self.points.append(bestfit_point)
            self.points.append(pywrappers.Point(contour.x_max, contour.y[contour.x.index(contour.x_max)], color=self.color_cycle.next()))
            self.points.append(pywrappers.Point(contour.x[contour.y.index(contour.y_min)], contour.y_min, color=self.color_cycle.next()))
            self.points.append(pywrappers.Point(contour.x[contour.y.index(contour.y_max)], contour.y_max, color=self.color_cycle.next()))


    def draw(self):
        self.last_bin_is_overflow = self.ptspectrum.last_bin_is_overflow()
        logging.info('Determined that last_bin_is_overflow={0} based on the ptspectrum {1}'.format(self.last_bin_is_overflow, self.ptspectrum))

        self.top_x_min    = self.obs.binning[0]
        self.bottom_x_min = self.obs.binning[0]

        if self.last_bin_is_overflow:
            self.top_x_max    = 2.*self.obs.binning[-2] - self.obs.binning[-3]
            self.bottom_x_max = 2.*self.obs.binning[-2] - self.obs.binning[-3]
        else:
            self.top_x_max    = self.obs.binning[-1]
            self.bottom_x_max = self.obs.binning[-1]

        self.get_points()
        self.draw_parametrizations()
        self.draw_ptscan()
        self.add_top(self.legend, '')

        super(BottomPanelPlotWithParametrizations, self).draw()

        self.draw_small_pad()

        c.Update()
        self.wrapup()


    def draw_parametrizations(self):

        SM_mu = pywrappers.Histogram(
            utils.get_unique_rootname(),
            'SM',
            self.obs.binning,
            [1.0 for i in xrange(len(self.obs.binning)-1)],
            color=16
            )
        self.add_bottom(SM_mu, 'repr_basic_histogram', self.legend)
        SM_xs = pywrappers.Histogram(
            utils.get_unique_rootname(),
            'SM',
            self.obs.binning,
            self.obs.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=self.last_bin_is_overflow),
            color=16
            )
        self.add_top(SM_xs, 'repr_basic_histogram', self.legend)

        parametrization = differentials.parametrization.WSParametrization(self.ws_file)
        parametrization.try_to_include_brmodifiers = True

        for point in self.points:
            parametrization.set(self.scan2D.x_variable, point.x)
            parametrization.set(self.scan2D.y_variable, point.y)
            mus = parametrization.get_mus_exp()

            if hasattr(point, '_is_bestfit'):
                title = 'Best fit'
            else:
                title = '{0}={1:.1f}, {2}={3:.1f}'.format(
                    self.scan2D.x_title, point.x, self.scan2D.y_title, point.y
                    )

            mu_histogram = pywrappers.Histogram(
                utils.get_unique_rootname(),
                title,
                self.obs.binning,
                mus,
                color=point.color
                )
            self.add_bottom(mu_histogram, 'repr_basic_histogram')

            xs_histogram = pywrappers.Histogram(
                utils.get_unique_rootname(),
                title,
                self.obs.binning,
                [ mu * xs for mu, xs in zip(mus, self.obs.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=self.last_bin_is_overflow)) ],
                color=point.color
                )
            logging.trace('Setting {0}.legend to {1}'.format(xs_histogram, self.legend))
            xs_histogram.legend = self.legend
            self.add_top(xs_histogram, 'repr_basic_histogram')

    def draw_ptscan(self):
        if self.ptspectrum.last_bin_is_overflow():
            binning = self.ptspectrum.binning()
            self.ptspectrum.hard_x_max = 2.*binning[-2] - binning[-3]

        ptscan_histogram = self.ptspectrum.to_hist()
        self.add_bottom(ptscan_histogram, 'repr_point_with_vertical_bar')

        ptscan_histogram_xs = self.ptspectrum.to_hist_xs()
        ptscan_histogram_xs.title = 'p_{T} combination'
        self.add_top(ptscan_histogram_xs, 'repr_point_with_vertical_bar', self.legend)


    def draw_small_pad(self):
        # Contour finding switches to c, slightly confusing... gather objects before
        histogram2D = self.scan2D.to_hist()
        small_pad_objs = histogram2D.repr_2D_with_contours_no_bestfit()

        # Draw the subplot
        self.toppad.cd()
        # c.cd()
        cw = 1.0 - self.toppad.GetLeftMargin() - self.toppad.GetRightMargin()
        ch = 1.0 - self.toppad.GetBottomMargin() - self.toppad.GetTopMargin()
        small_pad = ROOT.TPad(
            utils.get_unique_rootname(), '',
            self.toppad.GetLeftMargin() + 0.50*cw, self.toppad.GetBottomMargin() + 0.50*ch,
            self.toppad.GetLeftMargin() + 0.99*cw, self.toppad.GetBottomMargin() + 0.99*ch,
            )
        ROOT.SetOwnership(small_pad, False)
        small_pad.SetBottomMargin( 0.14 )
        small_pad.SetTopMargin(    0.03 )
        small_pad.SetLeftMargin(   0.14 )
        # small_pad.SetRightMargin(  0.10 )
        small_pad.SetRightMargin(  0.10 )
        small_pad.Draw()
        small_pad.cd()

        for obj, draw_str in small_pad_objs:
            obj.Draw(draw_str)

        axis_holder = small_pad_objs[0][0]
        axis_holder.GetXaxis().SetTitle(self.scan2D.x_title)
        axis_holder.GetYaxis().SetTitle(self.scan2D.y_title)

        axis_holder.GetXaxis().SetLabelSize(0.06)
        axis_holder.GetYaxis().SetLabelSize(0.06)
        axis_holder.GetXaxis().SetTitleSize(0.07)
        axis_holder.GetYaxis().SetTitleSize(0.07)


        for point in self.points:
            point.SetMarkerSize(2.0)
            point.Draw('repr_diamond_with_border')

        sm_point = pywrappers.Point(self.x_SM, self.y_SM)
        sm_point.SetMarkerSize(1.0)
        sm_point.Draw('repr_SM_point')

        dummy_legend = pywrappers.ContourDummyLegend(
            small_pad.GetLeftMargin() + 0.11,
            1. - small_pad.GetTopMargin() - 0.1,
            1. - small_pad.GetRightMargin() - 0.01,
            1. - small_pad.GetTopMargin() - 0.01,
            )
        dummy_legend.disable_bestfit = True
        # dummy_legend.disable_SM = True
        dummy_legend.Draw()

