import ROOT
import plotting_utils as utils
import pywrappers
from canvas import c

import differentials.logger as logger
import differentials.parametrization

from math import isnan, isinf, log10
from array import array


class BottomPanelPlot(object):
    """docstring for BottomPanelPlot"""
    def __init__(self, name):
        self.name = name

        self.top_objects = []
        self.bottom_objects = []

        self.x_title = ''
        self.y_title_top = ''
        self.y_title_bottom = ''

        self.pad_split_point = 0.33

        self.debug_lines = False

        self.top_bottom_margin = 0.02
        self.top_top_margin    = 0.10
        self.top_left_margin   = 0.12
        self.top_right_margin  = 0.02

        self.bottom_bottom_margin = 0.30
        self.bottom_top_margin    = 0.00
        self.bottom_left_margin   = 0.12
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
        obj._legend = leg
        self.top_objects.append((obj, draw_str))

    def add_bottom(self, obj, draw_str=None, leg=None):
        obj._legend = leg
        self.bottom_objects.append((obj, draw_str))

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


    def get_obj_extrema(self, obj, only_positive=False):
        x_min = None
        x_max = None
        y_min = None
        y_max = None
        if hasattr(obj, 'x_min'):
            x_min = obj.x_min()
        if hasattr(obj, 'x_max'):
            x_max = obj.x_max()
        if hasattr(obj, 'y_min'):
            if only_positive:
                try:
                    y_min = obj.y_min(only_positive=True)
                except TypeError:
                    y_min = obj.y_min()
                    if y_min < 0.0:
                        y_min = None
            else:
                y_min = obj.y_min()
        if hasattr(obj, 'y_max'):
            if only_positive:
                try:
                    y_max = obj.y_max(only_positive=True)
                except TypeError:
                    y_max = obj.y_max()
                    if y_max < 0.0:
                        y_max = None
            else:
                y_max = obj.y_max()
        return x_min, x_max, y_min, y_max


    def get_top_extrema(self):
        xs_min = []
        xs_max = []
        ys_min = []
        ys_max = []
        for obj in [obj for obj, _ in self.top_objects]:
            x_min, x_max, y_min, y_max = self.get_obj_extrema(obj, only_positive=self.top_log_scale)
            if not(x_min is None): xs_min.append(x_min)
            if not(x_max is None): xs_max.append(x_max)
            if not(y_min is None): ys_min.append(y_min)
            if not(y_max is None): ys_max.append(y_max)

        y_min = min(ys_min) if len(ys_min) else 0.001
        x_min = min(xs_min) if len(xs_min) else 0.001
        y_max = max(ys_max) if len(ys_max) else 1.0
        x_max = max(xs_max) if len(xs_max) else 1.0

        if not(self.top_x_min is None): x_min = self.top_x_min
        if not(self.top_y_min is None): y_min = self.top_y_min
        if not(self.top_x_max is None): x_max = self.top_x_max
        if not(self.top_y_max is None): y_max = self.top_y_max

        if self.top_log_scale:
            y_max_abs = y_max
            y_min_abs = y_min
            y_max = y_max_abs + (y_max_abs/y_min_abs)**0.2 * y_max_abs
            y_min = 0.5*y_min_abs
        else:
            dy = y_max - y_min
            y_min -= 0.1*dy
            y_max += 0.1*dy

        return x_min, y_min, x_max, y_max


    def get_bottom_extrema(self):
        xs_min = []
        xs_max = []
        ys_min = []
        ys_max = []
        for obj in [obj for obj, _ in self.bottom_objects]:
            x_min, x_max, y_min, y_max = self.get_obj_extrema(obj)
            if not(x_min is None): xs_min.append(x_min)
            if not(x_max is None): xs_max.append(x_max)
            if not(y_min is None): ys_min.append(y_min)
            if not(y_max is None): ys_max.append(y_max)

        y_min = min(ys_min) if len(ys_min) else 0.001
        x_min = min(xs_min) if len(xs_min) else 0.001
        y_max = max(ys_max) if len(ys_max) else 1.0
        x_max = max(xs_max) if len(xs_max) else 1.0

        if not(self.bottom_x_min is None): x_min = self.bottom_x_min
        if not(self.bottom_y_min is None): y_min = self.bottom_y_min
        if not(self.bottom_x_max is None): x_max = self.bottom_x_max
        if not(self.bottom_y_max is None): y_max = self.bottom_y_max

        dy = y_max - y_min
        y_min -= 0.1*dy
        y_max += 0.1*dy

        return x_min, y_min, x_max, y_max


    def draw(self, inplace=True):
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
        base_bottom = utils.get_plot_base(
            x_min=bottom_x_min, x_max=bottom_x_max, y_min=bottom_y_min, y_max=bottom_y_max,
            x_title='Observable', y_title='#mu'
            )
        base_bottom.Draw('P')

        # Draw the actual objects
        toppad.cd()
        if self.top_log_scale:
            toppad.SetLogy()

        logger.debug('top_objects: {0}'.format(self.top_objects))
        for obj, drawStr in self.top_objects:
            logger.debug('Attempting to draw obj: {0}, drawStr: {1}'.format(obj, drawStr))
            if drawStr is None:
                obj.Draw()
            else:
                obj.Draw(drawStr)

        bottompad.cd()
        logger.debug('bottom_objects: {0}'.format(self.bottom_objects))
        for obj, drawStr in self.bottom_objects:
            logger.debug('Attempting to draw obj: {0}, drawStr: {1}'.format(obj, drawStr))
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

        if inplace:
            self.wrapup()
        else:
            self.bottompad = bottompad
            self.toppad = toppad


    def wrapup(self):
        c.save(self.name)
        c.SetCanvasSize(self._tmp_width, self._tmp_height)



class BottomPanelPlotWithParametrizations(BottomPanelPlot):
    """docstring for BottomPanelPlotWithParametrizations"""
    def __init__(self, *args, **kwargs):
        super(BottomPanelPlotWithParametrizations, self).__init__(*args, **kwargs)
        self.scan2D = None
        self.pt_scan = None
        self.ws_file = None
        self.obs = None
        self.points = []
        self.parametrized_graphs = []
        self.color_cycle = utils.new_color_cycle()


    def get_points(self):
        contour = self.scan2D.to_hist().get_most_probable_1sigma_contour()
        self.points.append(pywrappers.Point(self.scan2D.bestfit().x, self.scan2D.bestfit().y, color=self.color_cycle.next()))
        self.points.append(pywrappers.Point(contour.x_min, contour.y[contour.x.index(contour.x_min)], color=self.color_cycle.next()))
        self.points.append(pywrappers.Point(contour.x_max, contour.y[contour.x.index(contour.x_max)], color=self.color_cycle.next()))
        self.points.append(pywrappers.Point(contour.x[contour.y.index(contour.y_min)], contour.y_min, color=self.color_cycle.next()))
        self.points.append(pywrappers.Point(contour.x[contour.y.index(contour.y_max)], contour.y_max, color=self.color_cycle.next()))


    def draw(self):
        self.top_x_min    = self.obs.binning[0]
        self.bottom_x_min = self.obs.binning[0]
        self.top_x_max    = self.obs.binning[-1]
        self.bottom_x_max = self.obs.binning[-1]

        self.bottom_y_min = 0.0
        self.bottom_y_max = 2.0

        self.get_points()
        self.draw_parametrizations()

        super(BottomPanelPlotWithParametrizations, self).draw(inplace=False)

        self.draw_small_pad()

        c.Update()
        self.wrapup()


    def draw_parametrizations(self):
        parametrization = differentials.parametrization.WSParametrization(self.ws_file)
        for point in self.points:
            parametrization.set(self.scan2D.x_variable, point.x)
            parametrization.set(self.scan2D.y_variable, point.y)
            mus = parametrization.get_mus_exp()

            mu_histogram = pywrappers.Histogram(
                utils.get_unique_rootname(),
                'title',
                self.obs.binning,
                mus
                )

            for obj, draw_str in mu_histogram.repr_basic_histogram():
                self.add_bottom(obj, draw_str)


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
        small_pad.SetLeftMargin(   0.12 )
        # small_pad.SetRightMargin(  0.10 )
        small_pad.SetRightMargin(  0.10 )
        small_pad.Draw()
        small_pad.cd()

        for obj, draw_str in small_pad_objs:
            obj.Draw(draw_str)

        for point in self.points:
            point.SetMarkerSize(2.0)
            point.Draw('repr_diamond_with_border')

            # point_edges = point.Clone() # This will be a TGraph then, not a wrapper
            # point_edges.SetMarkerSize(4.5)
            # point_edges.SetMarkerStyle(27)
            # point_edges.SetMarkerColor(1)
            # point_edges.Draw('PSAME')

        dummy_legend = pywrappers.ContourDummyLegend(
            small_pad.GetLeftMargin() + 0.01,
            1. - small_pad.GetTopMargin() - 0.1,
            1. - small_pad.GetRightMargin() - 0.01,
            1. - small_pad.GetTopMargin() - 0.01,
            )
        dummy_legend.disable_bestfit = True
        dummy_legend.disable_SM = True
        dummy_legend.Draw()





