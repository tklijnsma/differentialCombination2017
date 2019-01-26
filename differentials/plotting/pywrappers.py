import itertools, copy, sys

import ROOT
import plotting_utils as utils
from canvas import c, global_color_cycle

import logging
import differentials
import differentials.logger

import numpy

import scipy.ndimage

from array import array

class Legend(object):
    """Wrapper for TLegend class that allows flexible drawing of a legend on a multi panel plot"""

    def __init__(
            self,
            x1=None, y1=None, x2=None, y2=None,
            ):

        self.clear()

        # first set defaults
        self._x1 = lambda c: c.GetLeftMargin()
        self._x2 = lambda c: 1. - c.GetRightMargin()
        self._y1 = lambda c: 1. - c.GetTopMargin() - 0.15
        self._y2 = lambda c: 1. - c.GetTopMargin()
        self.set(x1, y1, x2, y2)

        self.auto_n_columns = True

        self.stylesheet = StyleSheet()
        self.stylesheet.plot_priority = 9000


    def __getattr__(self, name):
        """
        Reroutes calls TLegendMultiPanel.xxxx to TLegendMultiPanel.legend.xxxx
        This method should only be called if the attribute could not be found in TLegendMultiPanel
        """
        return getattr(self.legend, name)

    def style(self):
        return self.stylesheet

    def clear(self):
        self._entries = []
        self.set_new_legend()

    def set_new_legend(self):
        self.legend = ROOT.TLegend(0., 1., 0., 1.)
        self.legend.SetBorderSize(0)
        self.legend.SetFillStyle(0)
        ROOT.SetOwnership(self.legend, False)

    def set(self, x1=None, y1=None, x2=None, y2=None):
        if not(x1 is None):
            self._x1 = x1
        if not(x2 is None):
            self._x2 = x2
        if not(y1 is None):
            self._y1 = y1
        if not(y2 is None):
            self._y2 = y2

    def AddEntry(self, *args):
        """Save entries in a python list, but add them to the actual legend later"""
        self._entries.insert(0, args)

    def SetNColumns(self, value):
        self.auto_n_columns = False
        self.legend.SetNColumns(value)

    def n_columns_heuristic(self):
        n_entries = len(self._entries)
        return min(n_entries, 4)

    def Draw(self, drawStr=''):
        logging.debug('Drawing Legend {0} on gPad {1} with {2} entries'.format(self, ROOT.gPad.GetName(), len(self._entries)))
        x1 = self._x1(ROOT.gPad) if callable(self._x1) else self._x1
        y1 = self._y1(ROOT.gPad) if callable(self._y1) else self._y1
        x2 = self._x2(ROOT.gPad) if callable(self._x2) else self._x2
        y2 = self._y2(ROOT.gPad) if callable(self._y2) else self._y2
        logging.debug('Coordinates: x1 = {0}, y1 = {1}, x2 = {2}, y2 = {3}'.format(x1, y1, x2, y2))

        self.legend.SetX1(x1)
        self.legend.SetY1(y1)
        self.legend.SetX2(x2)
        self.legend.SetY2(y2)

        if self.auto_n_columns:
            self.legend.SetNColumns(self.n_columns_heuristic())

        for args in self._entries:
            if 'MOVE_TO_BOTTOM' in args: continue
            logging.debug('Adding entry with args: {0}'.format(args))
            self.legend.AddEntry(*args)
        for args in self._entries:
            if not('MOVE_TO_BOTTOM' in args): continue
            logging.debug('Adding entry with args: {0}'.format(args))
            args = [ arg for arg in args if arg != 'MOVE_TO_BOTTOM' ]
            self.legend.AddEntry(*args)

        self.legend.Draw(drawStr)


class Box(object):
    def __init__(self, x1=None, y1=None, x2=None, y2=None, color=None):
        super(Box, self).__init__()

        self.fill_alpha = 0.2
        if color is None:
            self.color = 14
        else:
            self.color = color

        self.set(x1, y1, x2, y2)

    def set(self, x1=None, y1=None, x2=None, y2=None):
        if x1 is None:
            self.x1 = 0.0
        else:
            self.x1 = x1
        if x2 is None:
            self.x2 = 1.0
        else:
            self.x2 = x2
        if y1 is None:
            self.y1 = 0.0
        else:
            self.y1 = y1
        if y2 is None:
            self.y2 = 1.0
        else:
            self.y2 = y2

    def Draw(self, draw_str=''):
        self.box = ROOT.TBox(self.x1, self.y1, self.x2, self.y2)
        ROOT.SetOwnership(self.box, False)
        self.box.SetLineWidth(0)
        self.box.SetFillColorAlpha(self.color, self.fill_alpha)
        self.box.Draw()

class ContourDummyLegend(Legend):
    """Special instance of Legend that creates some default objects and stores them in the legend"""

    dummy1sigma = ROOT.TGraph( 1, array( 'f' , [-999.] ), array( 'f' , [-999.] )  )
    dummy1sigma.SetLineWidth(2)
    dummy1sigma.SetName('dummy1sigma')
    ROOT.SetOwnership( dummy1sigma, False )

    dummy2sigma = ROOT.TGraph( 1, array( 'f' , [-999.] ), array( 'f' , [-999.] )  )
    dummy2sigma.SetLineWidth(2)
    dummy2sigma.SetLineStyle(2)
    dummy2sigma.SetName('dummy2sigma')
    ROOT.SetOwnership( dummy2sigma, False )

    dummySM = ROOT.TGraph( 1, array( 'f' , [-999.] ), array( 'f' , [-999.] )  )
    dummySM.SetMarkerSize(2.5)
    # dummySM.SetMarkerStyle(21)
    # dummySM.SetMarkerColor(16)
    dummySM.SetMarkerStyle(29)
    dummySM.SetMarkerColor(1)
    dummySM.SetName('dummySM')
    ROOT.SetOwnership( dummySM, False )

    dummybestfit = ROOT.TGraph( 1, array( 'f' , [-999.] ), array( 'f' , [-999.] )  )
    dummybestfit.SetMarkerSize(2)
    dummybestfit.SetMarkerStyle(34)
    dummybestfit.SetName('dummybestfit')
    ROOT.SetOwnership( dummybestfit, False )

    def __init__(self, *args, **kwargs):
        super(ContourDummyLegend, self).__init__(*args, **kwargs)
        self.disable_1sigma = False
        self.disable_2sigma = False
        self.disable_SM = False
        self.disable_bestfit = False

    def Draw(self, draw_str=''):
        if not self.disable_1sigma:
            self.AddEntry(self.dummy1sigma.GetName(), '1 #sigma', 'l')
            self.dummy1sigma.Draw('P')
        if not self.disable_2sigma:
            self.AddEntry(self.dummy2sigma.GetName(), '2 #sigma', 'l')
            self.dummy2sigma.Draw('P')
        if not self.disable_SM:
            self.AddEntry(self.dummySM.GetName(), 'SM', 'p')
            self.dummySM.Draw('P')
        if not self.disable_bestfit:
            self.AddEntry(self.dummybestfit.GetName(), 'Best fit', 'p')
            self.dummybestfit.Draw('P')        
        super(ContourDummyLegend, self).Draw(draw_str)


class Latex(object):
    """docstring"""
    def __init__(self, x, y, text):
        self.x = x
        self.y = y
        self.text = text
        self.tlatex = ROOT.TLatex()
        ROOT.SetOwnership(self.tlatex, False)

    def __getattr__(self, name):
        """
        Reroutes calls TLatexMultiPanel.xxxx to TLatexMultiPanel.tlatex.xxxx
        This method should only be called if the attribute could not be found in TLatexMultiPanel
        """
        return getattr(self.tlatex, name)

    def Draw(self, draw_str=None):
        x = self.x(ROOT.gPad) if callable(self.x) else self.x
        y = self.y(ROOT.gPad) if callable(self.y) else self.y
        logging.debug(
            'Drawing {0}; x={1}, y={2}, text={3}'
            .format(self, x, y, self.text)
            )
        self.tlatex.DrawLatex(x, y, self.text)


class CMS_Latex_type(Latex):
    """
    Specific implementation of Latex, that prints "CMS Preliminary" or
    "CMS Supplementary" at the default positions w.r.t. a TPad
    """
    CMS_type_str = 'Preliminary'
    apply_text_offset = True
    text_size    = 0.06

    def __init__(self, type_str=None, text_size=None):
        if not(text_size is None):
            self.text_size = text_size

        if type_str is None: type_str = self.CMS_type_str
        text = '#bf{{CMS}} #it{{#scale[0.75]{{{0}}}}}'.format(type_str)
        x = lambda c: c.GetLeftMargin()
        y = self._y
        super(CMS_Latex_type, self).__init__(x, y, text)

    def _y(self, c):
        """Internal use only"""
        return 1.-c.GetTopMargin() + self.get_text_offset()

    def get_text_offset(self):
        if self.apply_text_offset:
            return 0.25 * self.text_size

    def Draw(self, *args, **kwargs):
        self.tlatex.SetNDC()
        self.tlatex.SetTextAlign(11)
        self.tlatex.SetTextFont(42)
        self.tlatex.SetTextSize(self.text_size)
        super(CMS_Latex_type, self).Draw(*args, **kwargs)


class CMS_Latex_lumi(Latex):
    """
    Specific implementation of Latex, that prints "CMS Preliminary" or
    "CMS Supplementary" at the default positions w.r.t. a TPad
    """
    CMS_lumi = 35.9
    apply_text_offset = True
    text_size    = 0.05

    def __init__(self, lumi=None, text_size=None):
        if lumi is None: lumi = self.CMS_lumi
        if not(text_size is None):
            self.text_size = text_size

        if isinstance(lumi, int):
            text = '{0} fb^{{-1}} (13 TeV)'.format(lumi)
        else:
            text = '{0:.1f} fb^{{-1}} (13 TeV)'.format(lumi)
        x = lambda c: 1.-c.GetRightMargin()
        y = self._y
        super(CMS_Latex_lumi, self).__init__(x, y, text)

    def _y(self, c):
        """Internal use only"""
        return 1.-c.GetTopMargin() + self.get_text_offset()

    def get_text_offset(self):
        if self.apply_text_offset:
            return 0.25 * self.text_size

    def Draw(self, *args, **kwargs):
        self.tlatex.SetNDC()
        self.tlatex.SetTextAlign(31)
        self.tlatex.SetTextFont(42)
        self.tlatex.SetTextSize(self.text_size)
        super(CMS_Latex_lumi, self).Draw(*args, **kwargs)


class Base(object):
    """Alternative for get_plot_base"""

    def __init__(self, x_min=0., y_min=0., x_max=1., y_max=1., x_title='', y_title=''):
        self.x_min   = x_min
        self.x_max   = x_max
        self.y_min   = y_min
        self.y_max   = y_max
        self.x_title = x_title
        self.y_title = y_title
        self.H = utils.get_plot_base(
            x_min   = self.x_min,
            x_max   = self.x_max,
            y_min   = self.y_min,
            y_max   = self.y_max,
            x_title = self.x_title,
            y_title = self.y_title,
            )

    def __getattr__(self, name):
        return getattr(self.H, name)

    def set_x_min(self, x_min):
        self.x_min = x_min
        self.H.GetXaxis().SetLimits(self.x_min, self.x_max)

    def set_x_max(self, x_max):
        self.x_max = x_max
        self.H.GetXaxis().SetLimits(self.x_min, self.x_max)

    def set_y_min(self, y_min):
        self.y_min = y_min
        self.H.SetMinimum(self.y_min)

    def set_y_max(self, y_max):
        self.y_max = y_max
        self.H.SetMaximum(self.y_max)

    def set(self, x_min=None, x_max=None, y_min=None, y_max=None):
        if not(x_min is None): self.set_x_min(x_min)
        if not(x_max is None): self.set_x_max(x_max)
        if not(y_min is None): self.set_y_min(y_min)
        if not(y_max is None): self.set_y_max(y_max)

    def Draw(self, draw_str=''):
        self.H.Draw('P')




class StyleSheet(object):
    """doc"""
    def __init__(self, **kwargs):
        super(StyleSheet, self).__init__()

        self.color = 1
        self.line_color = None
        self.marker_color = None
        self.fill_color = None

        self.fill_color_alpha = 1.0

        self.line_style = 1
        self.marker_style = 8
        self.fill_style = None

        self.line_width = 2
        self.error_bar_line_width = 1

        self.marker_size = 1

        self.plot_priority = 10

        self.bin_center_offset = 0.0

        self.x_width = 0.0

        for key, value in kwargs.iteritems():
            setattr(self, key, value)


    def get_line_color(self):
        if self.line_color is None:
            return self.color
        return self.line_color

    def get_marker_color(self):
        if self.marker_color is None:
            return self.color
        return self.marker_color

    def get_fill_color(self):
        if self.fill_color is None:
            return self.color
        return self.fill_color


    def apply(self, O):
        # Colors, falls back to self.color if None
        if hasattr(O, 'SetLineColor'): O.SetLineColor(self.get_line_color())
        if hasattr(O, 'SetMarkerColor'): O.SetMarkerColor(self.get_marker_color())

        if self.fill_color_alpha != 1.0:
            if hasattr(O, 'SetFillColorAlpha'):
                O.SetFillColorAlpha(self.get_fill_color(), self.fill_color_alpha)
            elif hasattr(O, 'SetFillColor'):
                O.SetFillColor(self.get_fill_color())
        elif not(self.fill_style is None):
            if hasattr(O, 'SetFillColor'): O.SetFillColor(self.get_fill_color())

        if hasattr(O, 'SetLineStyle'): O.SetLineStyle(self.line_style)
        if hasattr(O, 'SetMarkerStyle'): O.SetMarkerStyle(self.marker_style)
        if not(self.fill_style is None) and hasattr(O, 'SetFillStyle'): O.SetFillStyle(self.fill_style)

        if hasattr(O, 'SetMarkerSize'): O.SetMarkerSize(self.marker_size)

        if hasattr(O, 'SetLineWidth'): O.SetLineWidth(self.line_width)

    def copy(self, **kwargs):
        new = copy.deepcopy(self)
        for key, value in kwargs.iteritems():
            setattr(new, key, value)
        return new


class BasicDrawable(object):
    """docstring for BasicDrawable"""

    default_stylesheet = StyleSheet()

    def __init__(self):
        super(BasicDrawable, self).__init__()
        self.legend = None
        self.move_to_bottom_of_legend = False
        self.stylesheets = [ BasicDrawable.default_stylesheet.copy() ]

    def has_legend(self):
        return not(self.legend is None)

    def add_to_legend(self, leg, *args):
        if leg is None: return
        if self.move_to_bottom_of_legend:
            args = list(args) + ['MOVE_TO_BOTTOM']
        leg.AddEntry(*args)

    def apply_style(self, O):
        for sheet in self.stylesheets:
            sheet.apply(O)

    def style(self):
        # For easy access to overwriting simple attributes
        return self.stylesheets[-1]

    def add_stylesheet(self, sheet):
        self.stylesheets.append(sheet)


class Histogram(BasicDrawable):
    """docstring for Histogram"""

    color_cycle = global_color_cycle
    fill_style_cycle = itertools.cycle([ 3245, 3254, 3205 ])

    def __init__(self, name, title, bin_boundaries, bin_values, color=None):
        super(Histogram, self).__init__()
        if name == 'auto':
            self.name = utils.get_unique_rootname()
        else:
            self.name = name
        self.title = title

        if len(bin_boundaries) != len(bin_values)+1:
            raise ValueError(
                'Histogram {4}: Inconsistent lengths; len(bin_boundaries)={0} != len(bin_values)+1={1}'
                '\nbin_boundaries: {2}'
                '\nbin_values: {3}'
                .format(
                    len(bin_boundaries), len(bin_values)+1,
                    bin_boundaries, bin_values,
                    self.title
                    )
                )

        self.bin_values = bin_values[:]
        self.bin_boundaries = bin_boundaries[:]
        self.n_bins = len(bin_boundaries)-1

        self.has_uncertainties = False

        self.style().color = color if not(color is None) else self.color_cycle.next()
        # self.style().fill_style = self.fill_style_cycle.next()

        # if color is None:
        #     self.color = self.color_cycle.next()
        # else:
        #     self.color = color
        # self.fill_style = self.fill_style_cycle.next()
        # self.fill_color = self.color
        # self.marker_style = 8
        # self.marker_size = 0
        # self.marker_color = self.color
        # self.line_color = self.color
        # self.line_width = 2

        self.last_bin_is_overflow = False
        self.mergemap = None


    def map_binning_to_fixed_binwidth(self, reference_bounds=None):
        if reference_bounds is None:
            reference_bounds = self.bin_boundaries[:]

        new_bounds = []
        for bound in self.bin_boundaries:
            if not bound in reference_bounds:
                raise ValueError(
                    'Bound {0} is not in the reference_bounds {1}.'
                    ' name = {2}, bin_boundaries = {3}'
                    .format(bound, reference_bounds, self.name, self.bin_boundaries)
                    )
            i = reference_bounds.index(bound)
            new_bounds.append(float(i))
        # print 'Mapped {0} to {1} for {2}'.format(
        #     self.bin_boundaries, new_bounds, self.name
        #     )
        self.reference_bounds = reference_bounds
        self.bin_boundaries = new_bounds
        self.get_merged_bins_from_reference(range(len(reference_bounds)))


    def get_merged_bins_from_reference(self, reference_bounds):
        self.mergemap = {}
        self.mergemap_indexed = {}
        for left, right in zip(self.bin_boundaries[:-1], self.bin_boundaries[1:]):
            i_left = reference_bounds.index(left)
            i_right = reference_bounds.index(right)

            self.mergemap[left] = reference_bounds[i_left:i_right+1]
            self.mergemap_indexed[self.bin_boundaries.index(left)] = range(i_left, i_right+1)


    def set_last_bin_is_overflow(self, flag=True, method='SECONDTOLASTBINWIDTH', hard_value=None):
        logging.debug('Last bin is specified to be overflow, so the last bin boundary will be modified')
        self.last_bin_is_overflow = True
        if method == 'SECONDTOLASTBINWIDTH':
            new_last_bin_boundary = self.bin_boundaries[-2] + (self.bin_boundaries[-2] - self.bin_boundaries[-3])
        elif method == 'HARDVALUE':
            if hard_value is None:
                raise TypeError('method \'HARDVALUE\' expects argument hard_value to be passed')
            new_last_bin_boundary = hard_value
        else:
            raise ValueError('Method \'{0}\' is not implemented'.format(method))
        logging.debug('Will replace last bin boundary {0} by {1}'.format(self.bin_boundaries[-1], new_last_bin_boundary))
        self.bin_boundaries[-1] = new_last_bin_boundary

    def x_min(self):
        return self.bin_boundaries[0]

    def x_max(self):
        return self.bin_boundaries[-1]

    def y_min(self, only_positive=False):
        if self.has_uncertainties:
            if only_positive:
                try:
                    return min([b for b in self.bounds_down if b>0.])
                except ValueError:
                    return 99999.
            else:
                return min(self.bounds_down)
        else:
            if only_positive:
                return min([b for b in self.bin_values if b>0.])
            else:
                return min(self.bin_values)            

    def y_max(self, only_positive=False):
        if self.has_uncertainties:
            if only_positive:
                return max([b for b in self.bounds_up if b>0.])
            else:
                return max(self.bounds_up)
        else:
            if only_positive:
                return max([b for b in self.bin_values if b>0.])
            else:
                return max(self.bin_values)            

    def drop_first_bin(self):
        logging.debug('Dropping first bin from histogram (name={0})'.format(self.name))
        self.bin_boundaries.pop(0)
        self.bin_values.pop(0)
        self.n_bins -= 1
        if self.has_uncertainties:
            self.errs_up.pop(0)
            self.bounds_up.pop(0)
            self.errs_down.pop(0)
            self.bounds_down.pop(0)

    def set_err_up(self, errs_up):
        self.errs_up = [ abs(i) for i in errs_up ]
        self.bounds_up = [ c+e for c, e in zip(self.bin_values, self.errs_up) ]
        self.has_uncertainties = True

    def set_err_down(self, errs_down):
        self.errs_down = [ abs(i) for i in errs_down ]
        self.bounds_down = [ c-abs(e) for c, e in zip(self.bin_values, self.errs_down) ]
        self.has_uncertainties = True

    def get_bin_centers(self):
        return [ 0.5*(left+right) for left, right in zip(self.bin_boundaries[:-1], self.bin_boundaries[1:]) ]

    def get_offset_bin_centers(self):
        if self.style().bin_center_offset == 0.0:
            bin_centers = self.get_bin_centers()
        else:
            bin_centers = []
            for center, width in zip(self.get_bin_centers(), self.get_bin_widths()):
                bin_centers.append( center + self.style().bin_center_offset * width )
        return bin_centers

    def get_offset_bin_centers_onlynonmerged(self):
        if self.mergemap is None:
            return self.get_offset_bin_centers()
        bin_centers = []
        for left, right in zip(self.bin_boundaries[:-1], self.bin_boundaries[1:]):
            width = right - left
            center = 0.5*(left+right)
            containing_bounds = self.mergemap[left]
            if len(containing_bounds) <= 2:
                # Only offset if bin is not merged
                center += self.style().bin_center_offset * width
            bin_centers.append(center)
        return bin_centers

    def get_bin_widths(self):
        return [ right-left for left, right in zip(self.bin_boundaries[:-1], self.bin_boundaries[1:]) ]

    def get_half_bin_widths(self):
        return [ 0.5*i for i in self.get_bin_widths() ]

    def get_zeroes(self):
        return [0.0 for i in xrange(self.n_bins)]

    def get_half_bin_widths_offsetcorrected(self, centers):
        width_left = []
        width_right = []
        for left, right, center in zip(self.bin_boundaries[:-1], self.bin_boundaries[1:], centers):
            width_left.append(center-left)
            width_right.append(right-center)
        return width_left, width_right

    #____________________________________________________________________

    # Standard TGraph with settable x uncertainties
    def y_error_TGraph(self, centers, x_widths_left, x_widths_right):
        Tg = ROOT.TGraphAsymmErrors(
            self.n_bins,
            array( 'f', centers ),
            array( 'f', self.bin_values ),
            array( 'f', x_widths_left ),
            array( 'f', x_widths_right ),
            array( 'f', self.errs_down ),
            array( 'f', self.errs_up ),
            )
        ROOT.SetOwnership(Tg, False)
        Tg.SetName(utils.get_unique_rootname())
        return Tg

    def y_error_TGraph_no_offset(self, x_widths):
        centers = self.get_bin_centers()
        return self.y_error_TGraph(centers, x_widths, x_widths)

    def y_error_TGraph_offset(self, x_widths):
        centers = self.get_offset_bin_centers()
        return self.y_error_TGraph(centers, x_widths, x_widths)

    def y_error_TGraph_offset_onlynonmerged(self, x_widths):
        centers = self.get_offset_bin_centers_onlynonmerged()
        return self.y_error_TGraph(centers, x_widths, x_widths)

    # Legacy
    def repr_basic_histogram(self, leg=None):
        return self.repr_basic(leg)

    #____________________________________________________________________
    # Draw methods
    def Draw(self, draw_style):
        logging.debug('Drawing Histogram {0} with draw_style {1}; legend: {2}'.format(self, draw_style, self.legend))
        for obj, draw_str in getattr(self, draw_style)(self.legend):
            obj.Draw(draw_str)

    # Basic line
    def repr_basic(self, leg=None):
        H = ROOT.TH1F(
            utils.get_unique_rootname(), '',
            len(self.bin_boundaries)-1, array( 'f', self.bin_boundaries)
            )
        ROOT.SetOwnership( H, False )
        for i_bin in xrange(self.n_bins):
            H.SetBinContent( i_bin+1, self.bin_values[i_bin] )

        _tmp_fill_style = self.style().fill_style
        self.style().fill_style = None
        self.apply_style(H)
        self.style().fill_style = _tmp_fill_style

        self.add_to_legend(leg, H.GetName(), self.title, 'l')
        return [ (H, 'HISTSAME') ]

    def repr_basic_dashed(self, leg=None):
        H, _ = self.repr_basic(leg)[0]
        H.SetLineStyle(2)
        return [(H, _)]


    #____________________________________________________________________
    # Basic vertical bar representations
    def repr_vertical_bar_nolegend(self, leg=None):
        Tg = self.y_error_TGraph_offset_onlynonmerged(x_widths = self.get_zeroes())
        self.apply_style(Tg)
        Tg.SetLineWidth(self.style().error_bar_line_width)
        return [(Tg, 'PSAME0')]

    def repr_narrow_bar(self, leg=None):
        if self.style().fill_style is None and self.style().fill_color_alpha == 1.0:
            self.style().fill_style = next(self.fill_style_cycle)
        Tg = self.y_error_TGraph_offset_onlynonmerged(x_widths = [ 0.1*w for w in self.get_bin_widths() ])
        self.apply_style(Tg)
        return [ (Tg, 'E2PSAME') ]

    def repr_full_bar(self, leg=None):
        if self.style().fill_style is None and self.style().fill_color_alpha == 1.0:
            self.style().fill_style = next(self.fill_style_cycle)
        Tg = self.y_error_TGraph_no_offset(x_widths = [ 0.5*w for w in self.get_bin_widths() ])
        self.apply_style(Tg)
        Tg.SetMarkerSize(0)
        self.add_to_legend(leg,  Tg.GetName(), self.title, 'LF' )
        return [ (Tg, 'E2PSAME') ]

    def repr_full_bar_with_connecting_line(self, leg=None):
        return [ (self.repr_full_bar(leg)[0][0], 'LE2SAME') ]

    # Add legend representations
    def repr_vertical_bar(self, leg=None):
        Tg, _ = self.repr_vertical_bar_nolegend()[0]
        self.add_to_legend(leg, Tg.GetName(), self.title, 'PE' )
        return [(Tg, _)]

    def repr_narrow_bar_onlyfill_legend(self, leg=None):
        Tg, _ = self.repr_narrow_bar()[0]
        self.add_to_legend(leg,  Tg.GetName(), self.title, 'F' )
        return [(Tg, _)]

    def repr_narrow_bar_linefill_legend(self, leg=None):
        Tg, _ = self.repr_narrow_bar()[0]
        self.add_to_legend(leg,  Tg.GetName(), self.title, 'LF' )
        return [(Tg, _)]

    def repr_full_bar_onlyfill_legend(self, leg=None):
        Tg, _ = self.repr_full_bar()[0]
        self.add_to_legend(leg,  Tg.GetName(), self.title, 'F' )
        return [(Tg, _)]

    def repr_full_bar_linefill_legend(self, leg=None):
        Tg, _ = self.repr_full_bar()[0]
        self.add_to_legend(leg,  Tg.GetName(), self.title, 'LF' )
        return [(Tg, _)]

    # Combinations with basic
    def repr_basic_with_narrow_fill(self, leg=None):
        return self.repr_basic_histogram() + self.repr_narrow_bar_linefill_legend(leg=leg)
    
    def repr_basic_with_full_fill(self, leg=None):
        return self.repr_basic_histogram() + self.repr_full_bar_linefill_legend(leg=leg)


    # Horizontal bars
    def repr_horizontal_bars(self, leg=None):
        Tg = ROOT.TGraphErrors(
            self.n_bins,
            array( 'f', self.get_bin_centers() ),
            array( 'f', self.bin_values ),
            array( 'f', self.get_half_bin_widths() ),
            array( 'f', self.get_zeroes() ),
            )
        ROOT.SetOwnership( Tg, False )
        Tg.SetName( utils.get_unique_rootname() )
        self.apply_style(Tg)
        Tg.SetMarkerSize(0)
        self.add_to_legend(leg, Tg.GetName(), self.title, 'PE' )
        return [ (Tg, 'EPSAME') ]

    def repr_horizontal_bars_dashed(self, leg=None):
        Tg, _ = self.repr_horizontal_bars(leg)[0]
        Tg.SetLineStyle(2)
        return [ (Tg, _) ]

    def repr_horizontal_bars_for_merged_bins(self, leg=None):
        if self.mergemap is None:
            logging.info('No mergemap for {0}'.format(self.title))
            return []
        r = []
        for bound, value in zip(self.bin_boundaries[:-1], self.bin_values):
            containing_bounds = self.mergemap[bound]
            if len(containing_bounds) > 2:
                line = ROOT.TLine(
                    containing_bounds[0], value, containing_bounds[-1], value
                    )
                ROOT.SetOwnership(line, False)
                self.apply_style(line)
                line.SetLineStyle(2)
                r.append((line, ''))
        return r


    # Horizontal + vertical representations
    def repr_vertical_bar_with_horizontal_lines_dashed(self, leg=None):
        return self.repr_vertical_bar(leg) + self.repr_horizontal_bars_dashed()

    def repr_vertical_bar_with_horizontal_lines_dashed_onlymerged(self, leg=None):
        return self.repr_vertical_bar(leg) + self.repr_horizontal_bars_for_merged_bins()



    #____________________________________________________________________

    def repr_point_with_vertical_bar_with_dashed_horizontal_bars(self, leg=None):
        return self.repr_point_with_vertical_bar(leg) + self.repr_horizontal_bars_dashed()

    def repr_point_with_vertical_bar_with_horizontal_bars_for_merged_bins(self, leg=None):
        return self.repr_point_with_vertical_bar(leg) + self.repr_horizontal_bars_for_merged_bins()

    def repr_point_with_horizontal_bar(self, leg=None):
        Tg, draw_str = self.repr_horizontal_bars()[0]
        Tg.SetMarkerSize(self.style().marker_size)
        return [(Tg, draw_str)]

    def repr_horizontal_bar_and_narrow_fill(self, leg=None):
        return self.repr_horizontal_bars() + self.repr_uncertainties_narrow_filled_area(leg)

    def repr_uncertainties_filled_area(self, leg=None):
        Tg = ROOT.TGraphAsymmErrors(
            self.n_bins,
            array( 'f', self.get_bin_centers() ),
            array( 'f', self.bin_values ),
            array( 'f', [ 0.45*w for w in self.get_bin_widths() ] ),
            array( 'f', [ 0.45*w for w in self.get_bin_widths() ] ),
            array( 'f', self.errs_down ),
            array( 'f', self.errs_up ),
            )
        ROOT.SetOwnership( Tg, False )
        Tg.SetName( utils.get_unique_rootname() )
        self.apply_style(Tg)
        self.add_to_legend(leg,  Tg.GetName(), self.title, 'LF' )
        return [ (Tg, 'E2PSAME') ]

    def repr_point_with_vertical_bar(self, leg=None):
        Tg = ROOT.TGraphAsymmErrors(
            self.n_bins,
            array( 'f', self.get_offset_bin_centers() ),
            array( 'f', self.bin_values ),
            array( 'f', [ 0.0 for i in xrange(self.n_bins) ] ),
            array( 'f', [ 0.0 for i in xrange(self.n_bins) ] ),
            array( 'f', self.errs_down ),
            array( 'f', self.errs_up ),
            )
        ROOT.SetOwnership( Tg, False )
        Tg.SetName( utils.get_unique_rootname() )
        self.apply_style(Tg)
        Tg.SetLineWidth(self.style().error_bar_line_width)
        self.add_to_legend(leg,  Tg.GetName(), self.title, 'PE' )
        return [(Tg, 'PSAME0')]

    def repr_point_with_vertical_bar_and_horizontal_bar(self, leg=None):
        return self.repr_point_with_vertical_bar(leg) + self.repr_horizontal_bars()

    def repr_uncertainties_narrow_filled_area(self, leg=None):
        Tg = ROOT.TGraphAsymmErrors(
            self.n_bins,
            array( 'f', self.get_bin_centers() ),
            array( 'f', self.bin_values ),
            array( 'f', [ 0.1*w for w in self.get_bin_widths() ] ),
            array( 'f', [ 0.1*w for w in self.get_bin_widths() ] ),
            array( 'f', self.errs_down ),
            array( 'f', self.errs_up ),
            )
        ROOT.SetOwnership( Tg, False )
        Tg.SetName( utils.get_unique_rootname() )

        self.style().fill_color_alpha = 0.3
        self.apply_style(Tg)

        self.add_to_legend(leg,  Tg.GetName(), self.title, 'F' )
        return [ (Tg, 'E2PSAME') ]


    def repr_uncertainties_fully_filled_area(self, leg=None):
        Tg = ROOT.TGraphAsymmErrors(
            self.n_bins,
            array( 'f', self.get_bin_centers() ),
            array( 'f', self.bin_values ),
            array( 'f', self.get_half_bin_widths() ),
            array( 'f', self.get_half_bin_widths() ),
            array( 'f', self.errs_down ),
            array( 'f', self.errs_up ),
            )
        ROOT.SetOwnership( Tg, False )
        Tg.SetName( utils.get_unique_rootname() )
    
        self.style().fill_color_alpha = 0.3
        self.apply_style(Tg)

        self.add_to_legend(leg,  Tg.GetName(), self.title, 'LF' )
        return [ (Tg, 'E2PSAME') ]



class Graph(BasicDrawable):
    """docstring for Graph"""

    color_cycle = global_color_cycle
    fill_style_cycle = itertools.cycle([ 3245, 3254, 3205 ])

    def __init__(self, name, title, xs, ys, color=None):
        super(Graph, self).__init__()        
        if name == 'auto':
            self.name = utils.get_unique_rootname()
        else:
            self.name = name
        self.title = title

        self.xs = xs
        self.ys = ys

        self.has_uncertainties = False
        if color is None:
            self.color = self.color_cycle.next()
        else:
            self.color = color

        self.fill_style = self.fill_style_cycle.next()
        self.line_width = 2
        self.line_style = 1

        self._filled_bestfit = False

    def fill_bestfit(self, x):
        self.x_bestfit = x
        self._filled_bestfit = True

    def multiply_x_by_constant(self, c):
        self.xs = [ c*x for x in self.xs ]
        if self._filled_bestfit:
            self.x_bestfit *= c

    def smooth_y(self, window_size=3):
        ys = numpy.array(self.ys)
        window = numpy.ones(window_size) / window_size
        ys_smooth = numpy.convolve(ys, window, mode='same')
        self.ys = list(ys_smooth)

    def SetLineWidth(self, width):
        self.line_width = width

    def SetLineStyle(self, style):
        self.line_style = style
        
    def set_err_up(self, errs_up):
        self.errs_up = [ abs(i) for i in errs_up ]
        self.bounds_up = [ c+e for c, e in zip(self.bin_values, self.errs_up) ]
        self.has_uncertainties = True

    def set_err_down(self, errs_down):
        self.errs_down = [ abs(i) for i in errs_down ]
        self.bounds_down = [ c+e for c, e in zip(self.bin_values, self.errs_down) ]
        self.has_uncertainties = True

    def flat_cut(self, y_top):
        # Get the two cuts at y_top, using the uncertainty calculator
        uncertaintycalculator = differentials.uncertaintycalculator.UncertaintyCalculator()
        uncertaintycalculator.cutoff = y_top
        unc = uncertaintycalculator.create_uncertainties(self.xs, self.ys)
        self.insert(unc.left_bound, y_top)
        self.insert(unc.right_bound, y_top)
        self.filter(y_max=y_top + 0.000001)

    def insert(self, x_new, y_new):
        if x_new < self.xs[0]:
            self.xs.insert(0, x_new)
            self.ys.insert(0, y_new)
            return
        elif x_new > self.xs[-1]:
            self.xs.append(x_new)
            self.ys.append(y_new)
            return
        for i, x in enumerate(self.xs):
            if x >= x_new:
                self.xs.insert(i, x_new)
                self.ys.insert(i, y_new)
                return

    def filter(self, x_min=-10e9, x_max=10e9, y_min=-10e9, y_max=10e9, inplace=True):
        passed_x = []
        passed_y = []
        for x, y in zip(self.xs, self.ys):
            if (x < x_min or x > x_max) or (y < y_min or y > y_max):
                continue
            else:
                passed_x.append(x)
                passed_y.append(y)
        if inplace:
            self.xs = passed_x
            self.ys = passed_y
        else:
            return passed_x, passed_y

    def check_input_sanity(self):
        if len(self.xs) == 0:
            raise ValueError(
                'Graph {0} ({1}) has zero entries in self.xs'
                .format(self.name, self.title)
                )
        if len(self.xs) != len(self.ys):
            raise ValueError(
                'Graph {0} ({1}) has unidentical lengths xs and ys: xs={2}, ys={3}'
                .format(self.name, self.title, len(self.xs), len(self.ys))
                )

    def repr_basic_line(self, leg=None):
        self.check_input_sanity()

        Tg = ROOT.TGraph(len(self.xs), array('f', self.xs), array('f', self.ys))
        ROOT.SetOwnership(Tg, False)
        Tg.SetName(utils.get_unique_rootname())

        Tg.SetLineColor(self.color)
        Tg.SetLineWidth(self.line_width)
        Tg.SetLineStyle(self.line_style)

        self.apply_style(Tg)

        if not(leg is None):
            self.add_to_legend(leg, Tg.GetName(), self.title, 'L' )

        return [(Tg, 'SAMEL')]


    def repr_markers(self, leg=None):
        Tg, _ = self.repr_basic_line()[0]
        if not(leg is None):
            self.add_to_legend(leg, Tg.GetName(), self.title, 'P' )
        return [(Tg, 'SAMEP')]

    def repr_dashed_line(self, leg=None):
        ret = self.repr_basic_line(leg)
        ret[0][0].SetLineStyle(2)
        return ret

    def repr_smooth_line(self, leg=None):
        Tg, _ = self.repr_basic_line(leg)[0]
        return [(Tg, 'SAMEC')]

    def repr_vertical_line_at_minimum(self, leg=None):
        x_at_minimum = self.xs[self.ys.index(min(self.ys))]
        logging.info('Vertical line: x minimum = {0}'.format(x_at_minimum))
        Tg = ROOT.TGraph(1,
            array('f', [x_at_minimum, x_at_minimum]),
            array('f', [0.0, 3.0])
            )
        ROOT.SetOwnership(Tg, False)
        Tg.SetName(utils.get_unique_rootname())
        Tg.SetLineColor(self.color)
        Tg.SetLineWidth(1)
        return [(Tg, 'LSAME')]

    def repr_smooth_and_vertical_line(self, leg=None):
        return self.repr_smooth_line(leg) + self.repr_vertical_line_at_minimum()

    def Draw(self, draw_style):
        for obj, draw_str in getattr(self, draw_style)(self.legend):
            obj.Draw(draw_str)

    def create_uncertainties(self, inplace=True, do_95percent_CL=False):
        logging.debug('Determining uncertainties for Graph')
        uncertaintycalculator = differentials.uncertaintycalculator.UncertaintyCalculator()
        uncertaintycalculator.cutoff = 1.0 # Graph will usually be filled with 2*deltaNLL
        if do_95percent_CL:
            uncertaintycalculator.cutoff = 3.841
        unc = uncertaintycalculator.create_uncertainties(xs = self.xs, deltaNLLs = self.ys)
        if unc.is_hopeless:
            logging.error(
                'Hopeless interpolation case: Unable to determine uncertainties for '
                'x = {0}, y = {1}, scandirs = {2}'
                .format(
                    self.x_variable if hasattr(self, 'x_variable') else '?',
                    self.y_variable if hasattr(self, 'y_variable') else '?',
                    self.scandirs if hasattr(self, 'scandirs') else '?'
                    )
                )
        if inplace:
            self.unc = unc
        else:
            return unc


class GraphDummy(Graph):
    """docstring for GraphDummy"""
    def __init__(self, title):
        super(GraphDummy, self).__init__('auto', title, [-999.], [-999.])
        


class Point(BasicDrawable):
    """docstring for Point"""

    color_cycle = global_color_cycle

    def __init__(self, x, y, color=None, marker_style=21, size=2):
        super(Point, self).__init__()
        self.x = x
        self.y = y
        self.legend = None

        if color is None:
            self.color = self.color_cycle.next()
        else:
            self.color = color
        self.marker_style = marker_style # 21 for SM, 34 for bestfit

        self.Tg = ROOT.TGraph(1, array('f', [self.x]), array('f', [self.y]))
        ROOT.SetOwnership(self.Tg, False)
        self.Tg.SetName(utils.get_unique_rootname())

        self.Tg.SetMarkerSize(size)


    def __getattr__(self, name):
        return getattr(self.Tg, name)

    def set_x(self, x):
        self.x = x
        self.Tg.SetPoint(1, self.x, self.y)

    def set_y(self, y):
        self.y = y
        self.Tg.SetPoint(1, self.x, self.y)

    def SetMarkerStyle(self, marker_style):
        self.marker_style = marker_style

    def SetMarkerColor(self, color):
        self.color = color

    def repr_basic(self, leg=None):
        self.Tg.SetMarkerStyle(self.marker_style)
        self.Tg.SetMarkerColor(self.color)
        Tg_copy = self.Tg.Clone()
        ROOT.SetOwnership(Tg_copy, False)
        return [(Tg_copy, 'PSAME')]

    def repr_filled_diamond(self, leg=None):
        self.Tg.SetMarkerStyle(33)
        self.Tg.SetMarkerColor(self.color)
        Tg_copy = self.Tg.Clone()
        ROOT.SetOwnership(Tg_copy, False)
        return [(Tg_copy, 'PSAME')]

    def repr_SM_point(self, leg=None):
        # self.Tg.SetMarkerStyle(21)
        # self.Tg.SetMarkerColor(16)
        self.Tg.SetMarkerSize(2.5)
        self.Tg.SetMarkerStyle(29)
        self.Tg.SetMarkerColor(1)
        Tg_copy = self.Tg.Clone()
        ROOT.SetOwnership(Tg_copy, False)
        return [(Tg_copy, 'PSAME')]

    def repr_empty_diamond(self, leg=None):
        self.Tg.SetMarkerStyle(27)
        self.Tg.SetMarkerColor(1)
        Tg_copy = self.Tg.Clone()
        ROOT.SetOwnership(Tg_copy, False)
        return [(Tg_copy, 'PSAME')]

    def repr_diamond_with_border(self, leg=None):
        return self.repr_filled_diamond() + self.repr_empty_diamond()

    def Draw(self, draw_style='repr_basic'):
        for obj, draw_str in getattr(self, draw_style)(self.legend):
            obj.Draw(draw_str)
        


class Histogram2D(BasicDrawable):
    """docstring for Histogram2D"""

    color_cycle = global_color_cycle
    default_value = 999.
    contour_filter = utils.ContourFilter()

    def __init__(
            self,
            name, title,
            color=None
            ):
        super(Histogram2D, self).__init__()
        if name == 'auto':
            self.name = utils.get_unique_rootname()
        else:
            self.name = name
        self.title = title

        self.has_uncertainties = False
        if color is None:
            self.color = self.color_cycle.next()
        else:
            self.color = color

        self.legend = None
        self.H2 = None
        self.H2_array = None
        self.entries = []

        self.contour_filter_method = None
        self._filled_bestfit = False


    def infer_bin_boundaries(self, bin_centers):
        bin_boundaries = []
        for i_bin in xrange(len(bin_centers)-1):
            bin_boundaries.append( 0.5*(bin_centers[i_bin]+bin_centers[i_bin+1]) )
        bin_boundaries = (
            [ bin_centers[0] - (bin_boundaries[0]-bin_centers[0]) ] +
            bin_boundaries +
            [ bin_centers[-1] + (bin_centers[-1]-bin_boundaries[-1]) ]
            )
        return bin_boundaries


    def bestfit(self):
        if self._filled_bestfit:
            return self._bestfit
        return self.entries[self.deltaNLL().index(0.0)]

    def x(self):
        return [e.x for e in self.entries]
    def y(self):
        return [e.y for e in self.entries]
    def z(self):
        return [e.z for e in self.entries]
    def deltaNLL(self):
        return [e.deltaNLL for e in self.entries]
    def two_times_deltaNLL(self):
        return [2.*e.deltaNLL for e in self.entries]

    def x_min(self):
        if len(self.entries) == 0: return self.x_bin_boundaries[0]
        return min(self.x())
    def x_max(self):
        if len(self.entries) == 0: return self.x_bin_boundaries[-1]
        return max(self.x())
    def y_min(self):
        if len(self.entries) == 0: return self.y_bin_boundaries[0]
        return min(self.y())
    def y_max(self):
        if len(self.entries) == 0: return self.y_bin_boundaries[-1]
        return max(self.y())

    def set_binning_from_entries(self):
        bestfit = self.bestfit()
        self.x_bin_centers = [ x for x in list(set(self.x())) if not x == bestfit.x ]
        self.y_bin_centers = [ y for y in list(set(self.y())) if not y == bestfit.y ]
        self.x_bin_centers.sort()
        self.y_bin_centers.sort()

        logging.trace('Found the following x_bin_centers:\n{0}'.format(self.x_bin_centers))
        logging.trace('Found the following y_bin_centers:\n{0}'.format(self.y_bin_centers))

        self.n_bins_x = len(self.x_bin_centers)
        self.n_bins_y = len(self.y_bin_centers)
        self.x_bin_boundaries = self.infer_bin_boundaries(self.x_bin_centers)
        self.y_bin_boundaries = self.infer_bin_boundaries(self.y_bin_centers)

        logging.trace('Found the following x_bin_boundaries:\n{0}'.format(self.x_bin_boundaries))
        logging.trace('Found the following y_bin_boundaries:\n{0}'.format(self.y_bin_boundaries))

    def fill_from_entries(self, entries=None):
        if not(entries is None):
            self.entries = entries
        bestfit = self.bestfit()
        self.set_binning_from_entries()

        self.H2 = ROOT.TH2D(
            utils.get_unique_rootname(), '',
            self.n_bins_x, array('d', self.x_bin_boundaries),
            self.n_bins_y, array('d', self.y_bin_boundaries),
            )
        ROOT.SetOwnership(self.H2, False)
        self.H2_array = [ [self.default_value for j in xrange(self.n_bins_y)] for i in xrange(self.n_bins_x) ]

        for i_x in xrange(self.n_bins_x):
            for i_y in xrange(self.n_bins_y):
                self.H2.SetBinContent(i_x+1, i_y+1, self.default_value)

        logging.debug('Filling {0} entries'.format(len(self.entries)))
        for entry in self.entries:
            if entry.x == bestfit.x and entry.y == bestfit.y: continue
            try:
                i_bin_x = self.x_bin_centers.index(entry.x)
            except ValueError:
                logging.error(
                    '{0} could not be filled - x={1} does not match any bin'
                    .format(entry, entry.x)
                    )
            try:
                i_bin_y = self.y_bin_centers.index(entry.y)
            except ValueError:
                logging.error(
                    '{0} could not be filled - y={1} does not match any bin'
                    .format(entry, entry.y)
                    )

            logging.trace(
                'Filling i_x={0} (x={1}) / i_y={2} (y={3}) with 2*deltaNLL={4}'
                .format(i_bin_x, entry.x, i_bin_y, entry.y, 2.*entry.deltaNLL)
                )

            self.H2.SetBinContent(i_bin_x+1, i_bin_y+1, 2.*entry.deltaNLL)
            self.H2_array[i_bin_x][i_bin_y] = 2.*entry.deltaNLL

    def fill_with_matrix(self, matrix, x_bin_boundaries, y_bin_boundaries):
        self.x_bin_boundaries = x_bin_boundaries
        self.y_bin_boundaries = y_bin_boundaries
        self.n_bins_x = len(self.x_bin_boundaries)-1
        self.n_bins_y = len(self.y_bin_boundaries)-1
        self.x_bin_centers = [ 0.5*(r+l) for l, r in zip(self.x_bin_boundaries[:-1], self.x_bin_boundaries[1:]) ]
        self.y_bin_centers = [ 0.5*(r+l) for l, r in zip(self.y_bin_boundaries[:-1], self.y_bin_boundaries[1:]) ]

        self.H2_array = matrix

        self.H2 = ROOT.TH2D(
            utils.get_unique_rootname(), '',
            self.n_bins_x, array('d', self.x_bin_boundaries),
            self.n_bins_y, array('d', self.y_bin_boundaries),
            )
        ROOT.SetOwnership(self.H2, False)

        min_z = 9999.
        min_x = 9999.
        min_y = 9999.
        for i_x in xrange(self.n_bins_x):
            for i_y in xrange(self.n_bins_y):
                z = self.H2_array[i_x][i_y]
                self.H2.SetBinContent(i_x+1, i_y+1, z)
                if z < min_z: # Simultaneously look for minimum
                    min_z = z
                    min_x = self.x_bin_centers[i_x]
                    min_y = self.y_bin_centers[i_y]
        self.fill_bestfit(min_x, min_y)

    def smooth_2d(self):
        smoothed = scipy.ndimage.gaussian_filter(
            self.H2_array, sigma = 2,
            )

        # kernel = numpy.array([
        #     [ 1., 1., 1. ],
        #     [ 1., 9., 1. ],
        #     [ 1., 1., 1. ],
        #     ])
        # kernel = kernel / numpy.sum(kernel)

        # print kernel
        # sys.exit()

        # smoothed = scipy.signal.convolve2d(
        #     self.H2_array, kernel, boundary='symm'
        #     )

        for i_x in xrange(self.n_bins_x):
            for i_y in xrange(self.n_bins_y):
                self.H2.SetBinContent(i_x+1, i_y+1, smoothed[i_x][i_y])




    def add_padding(self, value, x_min=-1000., x_max=1000., y_min=-1000., y_max=1000.):
        if x_min > self.x_bin_boundaries[0]:
            logging.error(
                'Cannot pad histogram: x_min {0} > self.x_bin_boundaries[0] {1}'
                .format(x_min, self.x_bin_boundaries[0])
                )
            return
        if x_max < self.x_bin_boundaries[-1]:
            logging.error(
                'Cannot pad histogram: x_max {0} < self.x_bin_boundaries[-1] {1}'
                .format(x_max, self.x_bin_boundaries[-1])
                )
            return
        if y_min > self.y_bin_boundaries[0]:
            logging.error(
                'Cannot pad histogram: y_min {0} > self.y_bin_boundaries[0] {1}'
                .format(y_min, self.y_bin_boundaries[0])
                )
            return
        if y_max < self.y_bin_boundaries[-1]:
            logging.error(
                'Cannot pad histogram: y_max {0} < self.y_bin_boundaries[-1] {1}'
                .format(y_max, self.y_bin_boundaries[-1])
                )
            return

        self.x_bin_boundaries = [x_min] + self.x_bin_boundaries + [x_max]
        self.y_bin_boundaries = [y_min] + self.y_bin_boundaries + [y_max]
        self.n_bins_x = len(self.x_bin_boundaries)-1
        self.n_bins_y = len(self.y_bin_boundaries)-1
        self.x_bin_centers = [ 0.5*(r+l) for l, r in zip(self.x_bin_boundaries[:-1], self.x_bin_boundaries[1:]) ]
        self.y_bin_centers = [ 0.5*(r+l) for l, r in zip(self.y_bin_boundaries[:-1], self.y_bin_boundaries[1:]) ]

        # Open up new hist with the padding bins
        H = ROOT.TH2D(
            utils.get_unique_rootname(), '',
            self.n_bins_x, array('d', self.x_bin_boundaries),
            self.n_bins_y, array('d', self.y_bin_boundaries),
            )
        ROOT.SetOwnership(H, False)

        # Fill in the new hist using old hist and the padding value
        for i_x in xrange(self.n_bins_x):
            for i_y in xrange(self.n_bins_y):
                if i_x == 0 or i_x == self.n_bins_x-1 or i_y == 0 or i_y == self.n_bins_y-1:
                    fill_val = value
                else:
                    fill_val = self.H2.GetBinContent(i_x, i_y)
                H.SetBinContent(i_x+1, i_y+1, fill_val)
        self.H2 = H

        # Also correct the H2_array
        self.H2_array = (
            [ value for i in xrange(self.n_bins_y) ]
            + [ [value] + r + [value] for r in self.H2_array ]
            + [ value for i in xrange(self.n_bins_y) ]
            )


    def mirror(self, del_y_middle=False):
        x_new_bounds = [ x for x in self.x_bin_boundaries if x > 0. ]
        x_new_bounds = [ -b for b in x_new_bounds[::-1] ] + x_new_bounds
        x_new_centers = [ 0.5*(r+l) for l, r in zip(x_new_bounds[:-1], x_new_bounds[1:]) ]
        x_new_nbins  = len(x_new_bounds)-1

        y_new_bounds = [ y for y in self.y_bin_boundaries if y > 0. ]
        y_new_bounds = [ -b for b in y_new_bounds[::-1] ] + y_new_bounds

        if del_y_middle:
            N = len(y_new_bounds) / 2
            y_new_bounds = y_new_bounds[:N] + [0.] + y_new_bounds[N+2:]
        y_new_centers = [ 0.5*(r+l) for l, r in zip(y_new_bounds[:-1], y_new_bounds[1:]) ]
        y_new_nbins  = len(y_new_bounds)-1

        H = ROOT.TH2D(
            utils.get_unique_rootname(), '',
            x_new_nbins, array('d', x_new_bounds),
            y_new_nbins, array('d', y_new_bounds),
            )
        ROOT.SetOwnership(H, False)

        # Set default value
        H_array = [ [ 888. for i_y in xrange(y_new_nbins) ] for i_x in xrange(x_new_nbins) ]
        for i_x in xrange(x_new_nbins):
            for i_y in xrange(y_new_nbins):
                H.SetBinContent(i_x+1, i_y+1, 888.)

        for i_x in xrange(self.n_bins_x):
            for i_y in xrange(self.n_bins_y):
                x_center = self.x_bin_centers[i_x]
                y_center = self.y_bin_centers[i_y]
                val = self.H2.GetBinContent(i_x+1, i_y+1)

                # Fill positive quadrant
                i_x_new = differentials.core.get_closest_match(x_center, x_new_centers)[1]
                i_y_new = differentials.core.get_closest_match(y_center, y_new_centers)[1]
                H.SetBinContent(i_x_new+1, i_y_new+1, val)
                H_array[i_x_new][i_y_new] = val

                # Fill negative quadrant
                i_x_new = differentials.core.get_closest_match(-x_center, x_new_centers)[1]
                i_y_new = differentials.core.get_closest_match(-y_center, y_new_centers)[1]
                H.SetBinContent(i_x_new+1, i_y_new+1, val)
                H_array[i_x_new][i_y_new] = val

        # Overwrite
        self.H2 = H
        self.H2_array = H_array
        self.x_bin_boundaries = x_new_bounds
        self.x_bin_centers    = x_new_centers
        self.n_bins_x         = x_new_nbins
        self.y_bin_boundaries = y_new_bounds
        self.y_bin_centers    = y_new_centers
        self.n_bins_y         = y_new_nbins


        # all_x_bin_boundaries = self.x_bin_boundaries
        # all_n_bins_x = self.n_bins_x

        # bounds = [ x for x in self.x_bin_boundaries if x > 0. ]
        # self.x_bin_boundaries = bounds[::-1] + bounds

        # self.n_bins_x = len(self.x_bin_boundaries)-1
        # self.x_bin_centers = [ 0.5*(r+l) for l, r in zip(self.x_bin_boundaries[:-1], self.x_bin_boundaries[1:]) ]

        # # self.H2_array = self.H2_array[-self.n_bins_x:]

        # H = ROOT.TH2D(
        #     utils.get_unique_rootname(), '',
        #     self.n_bins_x, array('d', self.x_bin_boundaries),
        #     self.n_bins_y, array('d', self.y_bin_boundaries),
        #     )
        # ROOT.SetOwnership(H, False)

        # # Fill in the new hist using old hist and the padding value
        # for i_x in xrange(self.n_bins_x):
        #     for i_y in xrange(self.n_bins_y):
        #         c = self.x_bin_boundaries[i_x]
        #         if c < 0.:
        #             pass

        #         # HIER VERDER
        #         # spiegelen afmaken

        #         H.SetBinContent(i_x+1, i_y+1, fill_val)
        # self.H2 = H



    def set_value_for_patch(self, value, x_min, x_max, y_min, y_max):
        x_min, ix_min = differentials.core.get_closest_match(x_min, self.x_bin_centers)
        x_max, ix_max = differentials.core.get_closest_match(x_max, self.x_bin_centers)
        y_min, iy_min = differentials.core.get_closest_match(y_min, self.y_bin_centers)
        y_max, iy_max = differentials.core.get_closest_match(y_max, self.y_bin_centers)
        # Insert patch
        for ix in xrange(ix_min, ix_max+1):
            for iy in xrange(iy_min, iy_max+1):
                self.H2.SetBinContent(ix+1, iy+1, value)

    def add_offset(self, value):
        for ix in xrange(self.n_bins_x):
            for iy in xrange(self.n_bins_y):
                self.H2.SetBinContent(
                    ix+1, iy+1,
                     self.H2.GetBinContent(ix+1, iy+1) + value
                    )

    def add_offset_to_zero(self):
        z_min = 10e9
        ix_min = -1
        iy_min = -1
        for ix in xrange(self.n_bins_x):
            for iy in xrange(self.n_bins_y):
                z = self.H2.GetBinContent(ix+1, iy+1)
                if z < z_min:
                    z_min = z
                    ix_min = ix
                    iy_min = iy
        self.add_offset(-z_min)

    def smooth_patch(self, x_min, x_max, y_min, y_max):
        # Get relevant patch
        x_min, ix_min = differentials.core.get_closest_match(x_min, self.x_bin_centers)
        x_max, ix_max = differentials.core.get_closest_match(x_max, self.x_bin_centers)
        y_min, iy_min = differentials.core.get_closest_match(y_min, self.y_bin_centers)
        y_max, iy_max = differentials.core.get_closest_match(y_max, self.y_bin_centers)

        # Get smoothed scan for whole matrix
        smoothed = scipy.ndimage.gaussian_filter(
            self.H2_array, sigma = 2,
            )

        # Insert patch
        for ix in xrange(ix_min, ix_max+1):
            for iy in xrange(iy_min, iy_max+1):
                self.H2.SetBinContent(ix+1, iy+1, smoothed[ix][iy])


    def polyfit_patch(self, x_min, x_max, y_min, y_max, order=7):
        if x_min > x_max: x_max, x_min = ( x_min, x_max )
        if y_min > y_max: y_max, y_min = ( y_min, y_max )
        x_min, ix_min = differentials.core.get_closest_match(x_min, self.x_bin_centers)
        x_max, ix_max = differentials.core.get_closest_match(x_max, self.x_bin_centers)
        y_min, iy_min = differentials.core.get_closest_match(y_min, self.y_bin_centers)
        y_max, iy_max = differentials.core.get_closest_match(y_max, self.y_bin_centers)

        # N = (ix_max+1-ix_min) * (iy_max+1-iy_min)

        T2D = ROOT.TGraph2D()
        ROOT.SetOwnership(T2D, False)
        i_point = 0
        for ix in xrange(ix_min, ix_max+1):
            for iy in xrange(iy_min, iy_max+1):
                T2D.SetPoint(i_point, self.x_bin_centers[ix], self.y_bin_centers[iy], self.H2_array[ix][iy])
                i_point += 1

        polyfit_factory = differentials.spline2d.Spline2DFactory()
        polyfit_factory.x_min = x_min
        polyfit_factory.x_max = x_max
        polyfit_factory.y_min = y_min
        polyfit_factory.y_max = y_max
        polyfit_factory.ord_polynomial = order
        polyfit = polyfit_factory.make_polyfit(T2D)
        polyfit.multiply_by_two = False

        for ix in xrange(ix_min, ix_max+1):
            for iy in xrange(iy_min, iy_max+1):
                x = self.x_bin_centers[ix]
                y = self.y_bin_centers[iy]
                self.H2.SetBinContent(ix+1, iy+1, polyfit.eval(x, y))


    def fill_bestfit(self, x, y):
        self._filled_bestfit = True
        self._bestfit = differentials.core.AttrDict(x=x, y=y)

    def repr_2D(self, leg=None):
        utils.set_color_palette()
        self.H2.SetMinimum(0.0)
        self.H2.SetMaximum(7.0)
        return [(self.H2, 'COLZ')]

    def repr_2D_rainbow(self, leg=None):
        utils.set_color_palette('rainbow')
        self.H2.SetMaximum(300.0)
        return [(self.H2, 'COLZ')]

    def repr_high_contours(self, leg=None):
        ret = []

        color_cycle = itertools.cycle([1, 2, 4])
        for level in [ 20., 50., 70., 100, 200., 500., 1000. ]:
            Tgs = utils.get_contours_from_H2(self.H2, level)
            labels = []
            color = color_cycle.next()
            for Tg in Tgs:
                Tg.SetLineColor(color)
                Tg.SetLineWidth(2)

                xs, ys = utils.get_x_y_from_TGraph(Tg)
                y = max(ys)
                x = xs[ys.index(y)]
                l = Latex(x, y, '{0:.1f}'.format(level))
                l.SetTextSize(0.02)
                l.SetTextColor(color)
                labels.append(l)
                ret.extend([ (Tg, 'LSAME') for Tg in Tgs ])
            ret.extend([ (l, '') for l in labels ])
        return ret

    def repr_2D_rainbow_high_contours(self, leg=None):
        return self.repr_2D_rainbow(leg) + self.repr_high_contours(leg)

    def repr_bestfitpoint(self, leg=None):
        bestfit = self.bestfit()
        Tg = ROOT.TGraph(1, array('f', [bestfit.x]), array('f', [bestfit.y]))
        ROOT.SetOwnership(Tg, False)
        Tg.SetMarkerSize(2)
        Tg.SetMarkerStyle(34)
        Tg.SetMarkerColor(self.color)
        Tg.SetName(utils.get_unique_rootname())
        return [(Tg, 'PSAME')]

    def repr_1sigma_contours(self, leg=None):
        Tgs = utils.get_contours_from_H2(self.H2, 2.30)
        if not(self.contour_filter_method is None):
            Tgs = getattr(self.contour_filter, self.contour_filter_method)(Tgs)[:1]
        if len(Tgs) == 0:
            logging.error('Could not extract 1 sigma contours from {0} ({1})'.format(self.name, self.title))
            return []
        for Tg in Tgs:
            Tg.SetLineColor(self.color)
            Tg.SetLineWidth(2)
        if not(leg is None):
            Tg = Tgs[0]
            Tg.SetName(utils.get_unique_rootname())
            leg.AddEntry(Tg.GetName(), self.title, 'l')
        return [ (Tg, 'LSAME') for Tg in Tgs ]

    def repr_2sigma_contours(self, leg=None):
        Tgs = utils.get_contours_from_H2(self.H2, 6.18)
        if not(self.contour_filter_method is None):
            Tgs = getattr(self.contour_filter, self.contour_filter_method)(Tgs)[:1]
        if len(Tgs) == 0:
            logging.error('Could not extract 2 sigma contours from {0} ({1})'.format(self.name, self.title))
            return []
        for Tg in Tgs:
            Tg.SetLineColor(self.color)
            Tg.SetLineWidth(2)
            Tg.SetLineStyle(2)
        if not(leg is None):
            Tg = Tgs[0]
            Tg.SetName(utils.get_unique_rootname())
            leg.AddEntry(Tg.GetName(), self.title, 'l')
        return [ (Tg, 'LSAME') for Tg in Tgs ]

    def repr_contours(self, leg=None):
        return self.repr_bestfitpoint() + self.repr_1sigma_contours(leg) + self.repr_2sigma_contours()

    def repr_contours_no_bestfit(self, leg=None):
        return self.repr_1sigma_contours(leg) + self.repr_2sigma_contours()

    def repr_1sigma_contours_with_bestfit(self, leg=None):
        return self.repr_bestfitpoint() + self.repr_1sigma_contours(leg)

    def repr_2D_with_contours(self, leg=None):
        return self.repr_2D() + self.repr_1sigma_contours(leg) + self.repr_2sigma_contours() + self.repr_bestfitpoint()

    def repr_2D_with_contours_no_bestfit(self, leg=None):
        return self.repr_2D() + self.repr_1sigma_contours() + self.repr_2sigma_contours()

    def Draw(self, draw_style):
        for obj, draw_str in getattr(self, draw_style)(self.legend):
            obj.Draw(draw_str)


    def get_most_probable_1sigma_contour(self, cutoff=2.30):
        allcontours = utils.get_contours_from_H2(self.H2, cutoff)
        if len(allcontours) == 0:
            raise RuntimeError('No contours at all found for cutoff {0}'.format(cutoff))
        candidatecontours = []

        bestfit = self.bestfit()

        for Tg in allcontours:
            Tg.x, Tg.y = utils.get_x_y_from_TGraph(Tg)
            Tg.x_min = min(Tg.x)
            Tg.x_max = max(Tg.x)
            Tg.y_min = min(Tg.y)
            Tg.y_max = max(Tg.y)

            # Check if bestfit is at least inside minima and maxima
            if not (
                    Tg.x_min < bestfit.x
                    and Tg.x_max > bestfit.x
                    and Tg.y_min < bestfit.y
                    and Tg.y_max > bestfit.y
                    ):
                continue

            # Compute some numbers that may help in selection: minimum distance to bestfit
            Tg.minDist = min([ ( Tg.x[i] - bestfit.x )**2 + ( Tg.y[i] - bestfit.y )**2 for i in xrange(Tg.GetN()) ])

            # Distance to bestfit should be close to half the difference between xMax-xMin (or yMax-yMin)
            # Compute ratio of this, minus 1 - result *should* be close to zero
            Tg.distRatio = Tg.minDist / ( 0.5 * min( Tg.y_max-Tg.y_min, Tg.x_max-Tg.x_min ) )  -  1.0
            # if abs(Tg.distRatio) > 0.8: continue

            candidatecontours.append( Tg )

        if len(candidatecontours) == 0:
            logging.error('Basic contour selection yielded no contours; returning first of all contours')
            return allcontours[0]
        elif len(candidatecontours) > 1:
            candidatecontours.sort( key = lambda Tg: Tg.minDist )
            candidatecontours = candidatecontours[:2]

        # # Pick the contour with the highest ratio (more likely to be the 'outer shell')
        # candidatecontours.sort( key = lambda Tg: Tg.distRatio, reverse=True )

        # Actually, pick the 'inner' shell, outer shell is too likely to be a misfit
        candidatecontours.sort( key = lambda Tg: Tg.distRatio )    
        contour = candidatecontours[0]
        return contour

    def get_most_probable_2sigma_contour(self):
        return self.get_most_probable_1sigma_contour(6.18)
        
    def get_extrema_from_contour(self, cutoff=2.30):
        if cutoff == '2sigma': cutoff = 6.18
        contour = self.get_most_probable_1sigma_contour(cutoff)
        xs, ys = utils.get_x_y_from_TGraph(contour)
        extrema = differentials.core.AttrDict(x_min=min(xs), x_max=max(xs), y_min=min(ys), y_max=max(ys))
        logging.info(
            'Found extreme at z={0}: x_min={1}, x_max={2}, y_min={3}, y_max={4}'
            .format(
                cutoff,
                extrema.x_min,
                extrema.x_max,
                extrema.y_min,
                extrema.y_max,
                )
            )
        return extrema



    def quickplot(
        self,
        plotname='auto',
        x_min=None, x_max=None, y_min=None, y_max=None,
        draw_method='repr_2D_with_contours',
        x_title='x', y_title='y',
        x_sm=1.0, y_sm=1.0,
        ):
        """Meant for a quick test plot, not finalized plots"""

        from differentials.plotting.plots import Single2DHistPlot
        if plotname == 'auto':
            plotname = 'quickhistplot_{0}'.format(self.name)
        if x_min is None:
            x_min = self.x_bin_boundaries[0]
        if x_max is None:
            x_max = self.x_bin_boundaries[-1]
        if y_min is None:
            y_min = self.y_bin_boundaries[0]
        if y_max is None:
            y_max = self.y_bin_boundaries[-1]

        c.resize_temporarily(870, 800)
        c.SetRightMargin(0.13)

        plot = differentials.plotting.plots.Single2DHistPlot(
            plotname,
            self,
            x_min = x_min,
            x_max = x_max,
            y_min = y_min,
            y_max = y_max,
            draw_str = draw_method
            )
        plot.set_ranges_by_contour = False
        plot.x_title = x_title
        plot.y_title = y_title
        plot.x_SM = x_sm
        plot.y_SM = y_sm
        plot.draw()
        plot.wrapup()



