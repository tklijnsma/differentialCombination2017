import ROOT
import plotting_utils as utils
import pywrappers
from canvas import c

import differentials
import differentials.core as core

import logging
from math import isnan, isinf, log10, sqrt
from array import array
from collections import namedtuple


class PlotBase(object):
    """docstring for PlotBase"""
    def __init__(self, plotname):
        self.plotname = plotname
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

        self.scans = []
        self.manual_graphs = []


    def add_scan(self, scan):
        self.scans.append(scan)

    def draw(self):
        c.Clear()
        c.set_margins()

        base = utils.get_plot_base(
            x_min=self.x_min, x_max=self.x_max, y_min=self.y_min, y_max=self.y_max,
            x_title=self.x_title, y_title='2#DeltaNLL'
            )
        base.Draw('P')

        sigma1_line = ROOT.TLine(self.x_min, 1.0, self.x_max, 1.0)
        ROOT.SetOwnership(sigma1_line, False)
        sigma1_line.SetLineColor(14)
        sigma1_line.SetLineStyle(2)
        sigma1_line.Draw()

        sigma2_line = ROOT.TLine(self.x_min, 2.0, self.x_max, 2.0)
        ROOT.SetOwnership(sigma2_line, False)
        sigma2_line.SetLineColor(14)
        sigma2_line.SetLineStyle(2)
        sigma2_line.Draw()

        for scan in self.scans:
            graph = scan.to_graph()
            graph.filter(y_max=3.5)
            graph._legend = self.leg
            graph.Draw(getattr(graph, 'draw_style', 'repr_basic_line'))

            left_point, right_point = self.get_unc_points(scan)

            line_bestfit = ROOT.TGraph(2,
                array('f', [scan.bestfit().x, scan.bestfit().x]),
                array('f', [0.0, 3.0]),
                )
            ROOT.SetOwnership(line_bestfit, False)
            line_bestfit.SetLineColor(graph.color)
            if 'dashed' in getattr(graph, 'draw_style', 'repr_basic_line'):
                line_bestfit.SetLineStyle(2)
            line_bestfit.Draw('SAMEL')


        for graph in self.manual_graphs:
            graph.filter(y_max=3.5)
            graph._legend = self.leg
            graph.Draw(getattr(graph, 'draw_style', 'repr_basic_line'))

            left_point, right_point = self.get_unc_points(graph)

            x_bestfit = graph.unc.x_min
            line_bestfit = ROOT.TGraph(2,
                array('f', [x_bestfit, x_bestfit]),
                array('f', [0.0, 3.0]),
                )
            ROOT.SetOwnership(line_bestfit, False)
            line_bestfit.SetLineColor(graph.color)
            if 'dashed' in getattr(graph, 'draw_style', 'repr_basic_line'):
                line_bestfit.SetLineStyle(2)
            line_bestfit.Draw('SAMEL')


        pywrappers.CMS_Latex_type().Draw()
        pywrappers.CMS_Latex_lumi().Draw()


    def get_unc_points(self, scan):
        """Should work in Graph as well"""

        left_point = ROOT.TGraph(1, array('f', [scan.unc.left_bound]), array('f', [1.0]))
        ROOT.SetOwnership(left_point, False)
        left_point.SetMarkerColor(scan.color)
        left_point.SetMarkerSize(1.1)
        left_point.SetMarkerStyle(8)
        if not scan.unc.well_defined_left_bound:
            left_point.SetMarkerStyle(5)
        left_point.Draw('PSAME')

        right_point = ROOT.TGraph(1, array('f', [scan.unc.right_bound]), array('f', [1.0]))
        ROOT.SetOwnership(right_point, False)
        right_point.SetMarkerColor(scan.color)
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
    def __init__(self, plotname, scans, x_min=None, x_max=None, y_min=None, y_max=None):
        super(MultiContourPlot, self).__init__(plotname)

        self.scans = scans

        self.x_min = x_min
        self.y_min = y_min
        self.x_max = x_max
        self.y_max = y_max

        self.draw_individual_contours=True
        self.only_1sigma_contours=False

        self.x_SM = 1.0
        self.y_SM = 1.0

        self.x_title = scans[0].x_title
        self.y_title = scans[0].y_title
        self.base = pywrappers.Base(x_title=self.x_title, y_title=self.y_title)

        self.legend = pywrappers.Legend(
            lambda c: c.GetLeftMargin() + 0.01,
            lambda c: c.GetBottomMargin() + 0.02,
            lambda c: 1 - c.GetRightMargin() - 0.01,
            lambda c: c.GetBottomMargin() + 0.09
            )

    def pre_draw(self):
        super(MultiContourPlot, self).pre_draw()
        self.base.Draw()

    def draw(self):
        super(MultiContourPlot, self).draw()

        self.histograms = []
        for scan in self.scans:
            histogram2D = scan.to_hist()
            histogram2D._legend = self.legend
            if self.only_1sigma_contours:
                histogram2D.Draw('repr_1sigma_contours_with_bestfit')
            else:
                histogram2D.Draw('repr_contours')
            self.histograms.append(histogram2D)

        self.legend.Draw()

        if self.x_min is None: self.x_min = min([H.x_min() for H in self.histograms])
        if self.y_min is None: self.y_min = min([H.y_min() for H in self.histograms])
        if self.x_max is None: self.x_max = max([H.x_max() for H in self.histograms])
        if self.y_max is None: self.y_max = max([H.y_max() for H in self.histograms])
        self.base.set(x_min=self.x_min, y_min=self.y_min, x_max=self.x_max, y_max=self.y_max)

        self.base.GetXaxis().SetTitle(self.x_title)
        self.base.GetYaxis().SetTitle(self.y_title)
        self.base.GetXaxis().SetTitleSize(0.06)
        self.base.GetXaxis().SetLabelSize(0.05)
        self.base.GetYaxis().SetTitleSize(0.06)
        self.base.GetYaxis().SetLabelSize(0.05)

        SM_point = pywrappers.Point(self.x_SM, self.y_SM)
        SM_point.Draw('repr_SM_point')

        pywrappers.CMS_Latex_type().Draw()
        pywrappers.CMS_Latex_lumi().Draw()

        pywrappers.ContourDummyLegend(
            c.GetLeftMargin() + 0.01,
            1. - c.GetTopMargin() - 0.1,
            1. - c.GetRightMargin() - 0.01,
            1. - c.GetTopMargin() - 0.01,
            ).Draw()

        self.wrapup()

    def wrapup(self):
        c.Update()
        c.RedrawAxis()
        self.save()        
        if self.draw_individual_contours:
            for scan in self.scans:
                scan.plot(self.plotname + '_' + scan.name)


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
        self.top_left_margin   = 0.14
        self.top_right_margin  = 0.02

        self.bottom_bottom_margin = 0.30
        self.bottom_top_margin    = 0.00
        self.bottom_left_margin   = 0.14
        self.bottom_right_margin  = 0.02

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


    def add_top(self, obj, draw_str=None, leg=None):
        if not(leg is None):
            obj._legend = leg
        self.top_objects.append((obj, draw_str))

    def add_bottom(self, obj, draw_str=None, leg=None):
        if not(leg is None):
            obj._legend = leg
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
        y_min = self.get_extremum('y_min', 0.001, objs, possible_overwrite_key='top_y_min', minimum=True, only_positive=True)
        y_max = self.get_extremum('y_max', 1.000, objs, possible_overwrite_key='top_y_max', minimum=False, only_positive=True)
        return x_min, y_min, x_max, y_max

    def get_bottom_extrema(self):
        logging.debug('Getting extrema for bottom objects')
        objs = [obj for obj, _ in self.bottom_objects]
        x_min = self.get_extremum('x_min', 0.001, objs,  possible_overwrite_key='bottom_x_min', minimum=True)
        x_max = self.get_extremum('x_max', 1.000, objs,  possible_overwrite_key='bottom_x_max', minimum=False)
        y_min = self.get_extremum('y_min', 0.001, objs,  possible_overwrite_key='bottom_y_min', minimum=True, only_positive=False)
        y_max = self.get_extremum('y_max', 1.000, objs,  possible_overwrite_key='bottom_y_max', minimum=False, only_positive=False)
        return x_min, y_min, x_max, y_max

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

        # Add some margins
        if self.top_log_scale:
            y_max_abs = top_y_max
            y_min_abs = top_y_min
            top_y_max = y_max_abs + (y_max_abs/y_min_abs)**0.2 * y_max_abs
            top_y_min = 0.5*y_min_abs
        else:
            dy = top_y_max - top_y_min
            top_y_min -= 0.1*dy
            top_y_max += 0.1*dy

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

        # Draw the actual objects
        toppad.cd()
        if self.top_log_scale:
            toppad.SetLogy()

        logging.debug('top_objects: {0}'.format(self.top_objects))
        for obj, drawStr in self.top_objects:
            logging.debug('Attempting to draw obj: {0}, drawStr: {1}'.format(obj, drawStr))
            if drawStr is None:
                obj.Draw()
            else:
                obj.Draw(drawStr)

        bottompad.cd()
        logging.debug('bottom_objects: {0}'.format(self.bottom_objects))
        for obj, drawStr in self.bottom_objects:
            logging.debug('Attempting to draw obj: {0}, drawStr: {1}'.format(obj, drawStr))
            if drawStr is None:
                obj.Draw()
            else:
                obj.Draw(drawStr)


        # Process titles, ticks, labels
        toppad.cd()
        # base_top = self.top_objects[0][0]
        base_top.GetXaxis().SetLabelOffset(999.)
        base_top.GetYaxis().SetTitle(self.y_title_top)

        bottompad.cd()
        # base_bottom = self.bottom_objects[0][0]
        # Set sizes of labels/ticks equal to top panel (undo automatic scaling by ROOT)
        base_bottom.GetXaxis().SetLabelSize(base_top.GetXaxis().GetLabelSize() * height_ratio)
        base_bottom.GetYaxis().SetLabelSize(base_top.GetYaxis().GetLabelSize() * height_ratio)
        base_bottom.GetXaxis().SetTickLength(base_top.GetXaxis().GetTickLength() * height_ratio)

        base_bottom.GetYaxis().SetTitle(self.y_title_bottom)
        base_bottom.GetXaxis().SetTitle(self.x_title)
        base_bottom.GetXaxis().SetTitleSize(base_top.GetXaxis().GetTitleSize() * height_ratio)
        base_bottom.GetYaxis().SetTitleSize(base_top.GetYaxis().GetTitleSize() * height_ratio)
        base_bottom.GetYaxis().SetTitleOffset(1./height_ratio)

        if self.CMS_labels:
            toppad.cd()
            pywrappers.CMS_Latex_type(text_size=0.08).Draw()
            pywrappers.CMS_Latex_lumi(text_size=0.07).Draw()

        self.bottompad = bottompad
        self.toppad = toppad


    def wrapup(self):
        super(BottomPanelPlot, self).save()
        c.SetCanvasSize(self._tmp_width, self._tmp_height)



class SpectraPlot(BottomPanelPlot):

    def __init__(self, plotname, spectra, x_min=None, x_max=None, y_min=None, y_max=None):
        super(SpectraPlot, self).__init__(plotname)

        self.spectra = spectra
        self.plotname = plotname
        # self.plot = differentials.plotting.multipanel.BottomPanelPlot(plotname)
        self.y_title_bottom = '#mu'

        self.obsname = 'obs_name'
        self.obsunit = None

        self.draw_multiscans = False
        self.scans_x_min = -10.0
        self.scans_x_max = 10.0

        self.leg = pywrappers.Legend(
            lambda c: c.GetLeftMargin() + 0.01,
            lambda c: 1 - c.GetTopMargin() - 0.10,
            lambda c: 1 - c.GetRightMargin() - 0.01,
            lambda c: 1 - c.GetTopMargin()
            )

    def make_SM_line(self, spectra_original, leg=None):
        spectra = spectra_original[:]
        spectra.sort(key=lambda s: -len(s.scans))
        smxs_histogram = spectra[0].to_hist_smxs()
        self.add_top(smxs_histogram, 'repr_basic_histogram', leg)
        sm_histogram = spectra[0].to_hist_sm()
        self.add_bottom(sm_histogram, 'repr_basic_histogram')

    def make_labels_for_overflow_spectra(self, spectra, obs_name):
        x_min, y_min, x_max, y_max = self.get_top_extrema()
        y_overflow = max([ s.scans[-1].unc.right_bound * s.smxs[-1] for s in spectra ])
        text_size = 0.03

        i_label = 0
        for spectrum in spectra:
            if getattr(spectrum, 'no_overflow_label', False): continue
            y_overflow_NDC = lambda c, i_lambda=i_label: (
                c.GetBottomMargin()
                + log10(y_overflow/y_min)/log10(y_max/y_min) * ( 1.0 - c.GetTopMargin() - c.GetBottomMargin() )
                + i_lambda * 4.5*text_size
                + 2.5*text_size
                )

            binning = spectrum.binning()
            scale = binning[-2] - binning[-3]
            scale_text = str(int(scale)) if float(scale).is_integer() else '{0:.2f}'.format(scale)
            l = pywrappers.Latex(
                    lambda c: 1. - c.GetRightMargin() - 0.01,
                    y_overflow_NDC,
                    '#int_{{#lower[0.3]{{{0}}}}}^{{#lower[0.0]{{#infty}}}}#sigma({1}) d{1} / {2}'.format(
                        binning[-2],
                        obs_name,
                        scale_text
                        )
                    )
            l.SetNDC()
            l.SetTextSize(text_size)
            l.SetTextColor(spectrum.color)
            l.SetTextAlign(31)
            self.add_top(l, '')
            i_label += 1

    def draw(self):
        if self.draw_multiscans:
            for spectrum in self.spectra:
                differentials.plotting.canvas.reset_global_color_cyle()
                differentials.plotting.pywrappers.Graph.color_cycle = differentials.plotting.canvas.global_color_cycle
                spectrum.scans_x_min = self.scans_x_min
                spectrum.scans_x_max = self.scans_x_max
                spectrum.plot_scans(self.plotname + '_scans_' + spectrum.name)

        self.obstitle = core.standard_titles.get(self.obsname, self.obsname)

        self.x_title = self.obstitle 
        if self.obsunit:
            self.x_title += ' ({0})'.format(self.obsunit)

        self.y_title_top = '#Delta#sigma/#Delta{0} (pb{1})'.format(
            self.obstitle,
            '/' + self.obsunit if not(self.obsunit is None) else ''
            )

        for spectrum in self.spectra:
            if not spectrum.smxs_set:
                raise RuntimeError('SM cross sections were not provided for spectrum')
            if not spectrum._is_read:
                spectrum.read()

        if self.spectra[0].last_bin_is_overflow():
            # Make sure overflow bins are aligned
            x_max = max([ 2*s.binning()[-2]-s.binning()[-3] for s in self.spectra ])
            for spectrum in self.spectra:
                spectrum.hard_x_max = x_max

        self.make_SM_line(self.spectra, leg=self.leg)
        for spectrum in self.spectra:
            self.add_bottom(spectrum.to_hist(), spectrum.draw_method)
            self.add_top(spectrum.to_hist_xs(), spectrum.draw_method, leg=self.leg)

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
        if self.default_points_xy_maxima:
            bestfit_point = pywrappers.Point(self.scan2D.bestfit().x, self.scan2D.bestfit().y, color=self.color_cycle.next())
            bestfit_point._is_bestfit = True
            self.points.append(bestfit_point)
            self.points.append(pywrappers.Point(contour.x_min, contour.y[contour.x.index(contour.x_min)], color=self.color_cycle.next()))
            self.points.append(pywrappers.Point(contour.x_max, contour.y[contour.x.index(contour.x_max)], color=self.color_cycle.next()))
            self.points.append(pywrappers.Point(contour.x[contour.y.index(contour.y_min)], contour.y_min, color=self.color_cycle.next()))
            self.points.append(pywrappers.Point(contour.x[contour.y.index(contour.y_max)], contour.y_max, color=self.color_cycle.next()))
        else:
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
            logging.trace('Setting {0}._legend to {1}'.format(xs_histogram, self.legend))
            xs_histogram._legend = self.legend
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

