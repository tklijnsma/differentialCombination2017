import ROOT
import plotting_utils as utils
from canvas import c

class Legend(object):
    """Wrapper for TLegend class that allows flexible drawing of a legend on a multi panel plot"""
    def __init__(
            self,
            x1, y1, x2, y2
            ):
        self._x1 = x1
        self._y1 = y1
        self._x2 = x2
        self._y2 = y2
        self._entries = []
        self.legend = ROOT.TLegend(0., 1., 0., 1.)
        ROOT.SetOwnership(self.legend, False)

    def __getattr__(self, name):
        """
        Reroutes calls TLegendMultiPanel.xxxx to TLegendMultiPanel.legend.xxxx
        This method should only be called if the attribute could not be found in TLegendMultiPanel
        """
        return getattr(self.legend, name)

    def AddEntry(self, *args):
        """Save entries in a python list, but add them to the actual legend later"""
        self._entries.append( args )

    def Draw(self, drawStr=''):
        x1 = self._x1(ROOT.gPad) if callable(self._x1) else self._x1
        y1 = self._y1(ROOT.gPad) if callable(self._y1) else self._y1
        x2 = self._x2(ROOT.gPad) if callable(self._x2) else self._x2
        y2 = self._y2(ROOT.gPad) if callable(self._y2) else self._y2

        self.legend.SetX1NDC(x1)
        self.legend.SetY1NDC(y1)
        self.legend.SetX2NDC(x2)
        self.legend.SetY2NDC(y2)

        for args in self._entries:
            self.legend.AddEntry(*args)
        self.legend.Draw(drawStr)


class Latex(object):
    """docstring"""
    def __init__(self, x, y, text):
        self.x = x
        self.y = y
        self.text = text
        self.tlatex = ROOT.TLatex()
        ROOT.SetOwnership( self.tlatex, False )

    def __getattr__(self, name):
        """
        Reroutes calls TLatexMultiPanel.xxxx to TLatexMultiPanel.tlatex.xxxx
        This method should only be called if the attribute could not be found in TLatexMultiPanel
        """
        return getattr(self.tlatex, name)

    def Draw(self, draw_str=None):
        x = self.x(ROOT.gPad) if callable(self.x) else self.x
        y = self.y(ROOT.gPad) if callable(self.y) else self.y
        print 'In Draw():'
        print '  x    = ', self.x
        print '  y    = ', self.y
        print '  text = ', self.text
        self.tlatex.DrawLatex(x, y, self.text)


class CMS_Latex_type(Latex):
    """
    Specific implementation of Latex, that prints "CMS Preliminary" or
    "CMS Supplementary" at the default positions w.r.t. a TPad
    """
    CMS_type_str = 'Preliminary'
    apply_text_offset = True
    text_size    = 0.06

    def __init__(self, text=None, text_size=None):
        if text is None: text = self.CMS_type_str
        if not(text_size is None):
            self.text_size = text_size

        text = '#bf{{CMS}} #it{{#scale[0.75]{{{0}}}}}'.format(self.CMS_type_str)
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
    text_size    = 0.06

    def __init__(self, lumi=None, text_size=None):
        if lumi is None: lumi = self.CMS_lumi
        if not(text_size is None):
            self.text_size = text_size

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


class Histogram(object):
    """docstring for GenericHistogram"""

    colorCycle = Commands.newColorCycle()
    fillStyleCycle = itertools.cycle([ 3245, 3254, 3205 ])

    def __init__(self, name, title, bin_boundaries, bin_values, color=None):
        super(GenericHistogram, self).__init__()
        self.name = name
        self.title = title
        self.bin_boundaries = bin_boundaries
        self.n_bins = len(bin_boundaries)-1
        self.bin_values = bin_values
        self.bins = []
        self._expecting_asymm_unc = False
        if color is None:
            self.color = self.colorCycle.next()
        else:
            self.color = color

        self.fill_style = fillStyleCycle.next()
        self.fill_color = self.color
        self.marker_style
        self.marker_size
        self.marker_color = self.color
        self.line_color = self.color
        self.line_width

        
    def set_err_up(self, errs_up):
        self.errs_up = [ abs(i) for i in errs_up ]
        self.bounds_up = [ c+e for c, e in zip(self.bin_values, self.errs_up) ]
        self._expecting_asymm_err = True

    def set_err_down(self, errs_down):
        self.errs_down = [ abs(i) for i in errs_down ]
        self.bounds_down = [ c+e for c, e in zip(self.bin_values, self.errs_down) ]
        self._expecting_asymm_err = True

    def get_bin_centers(self):
        return [ 0.5*(left+right) for left, right in zip(self.bin_boundaries[:-1], self.bin_boundaries[1:]) ]

    def get_bin_widths(self):
        return [ right-left for left, right in zip(self.bin_boundaries[:-1], self.bin_boundaries[1:]) ]

    def get_half_bin_widths(self):
        return [ 0.5*i for i in self.get_bin_widths() ]

    def get_zeroes(self):
        return [0.0 for i in xrange(self.n_bins)]


    def repr_basic_histogram(self, leg=None):
        H = ROOT.TH1F(
            GetUniqueRootName(), '',
            len(self.bin_boundaries)-1, array( 'f', self.bin_boundaries)
            )
        ROOT.SetOwnership( H, False )
        H.SetLineColor(self.color)
        H.SetLineWidth(2)
        for i_bin in xrange(self.n_bins):
            H.SetBinContent( i_bin+1, self.bin_values[i_bin] )

        if not(leg is None):
            leg.AddEntry(H.GetName(), self.title, 'l')

        return [ (H, 'HISTSAME') ]


    def repr_horizontal_bars(self):
        Tg = ROOT.TGraphErrors(
            self.n_bins,
            array( 'f', self.get_bin_centers() ),
            array( 'f', self.bin_values ),
            array( 'f', self.get_half_bin_widths() ),
            array( 'f', self.get_zeroes() ),
            )
        ROOT.SetOwnership( Tg, False )
        Tg.SetName( GetUniqueRootName() )

        Tg.SetLineColor(   getattr(self, 'setLineColor',   self.color ) )
        # Tg.SetLineWidth(   getattr(self, 'setLineWidth',   2 ) )

        return [ (Tg, 'EPSAME') ]


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
        Tg.SetName( GetUniqueRootName() )
    
        Tg.SetFillStyle(   getattr(self, 'setFillStyle',   3245 ) )
        Tg.SetMarkerStyle( getattr(self, 'setMarkerStyle', 8 ) )
        Tg.SetMarkerSize(  getattr(self, 'setMarkerSize',  0 ) )
        Tg.SetFillColor(   getattr(self, 'setFillColor',   self.color ) )
        Tg.SetMarkerColor( getattr(self, 'setMarkerColor', self.color ) )
        Tg.SetLineColor(   getattr(self, 'setLineColor',   self.color ) )

        if not(leg is None):
            leg.AddEntry( Tg.GetName(), self.title, 'LF' )

        return [ (Tg, 'E2PSAME') ]


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
        Tg.SetName( GetUniqueRootName() )
    
        Tg.SetMarkerStyle( getattr(self, 'setMarkerStyle', 8 ) )
        Tg.SetMarkerSize(  getattr(self, 'setMarkerSize',  0 ) )
        Tg.SetMarkerColor( getattr(self, 'setMarkerColor', self.color ) )
        Tg.SetLineColor(   getattr(self, 'setLineColor',   self.color ) )

        # Tg.SetFillStyle(   getattr(self, 'setFillStyle',   3245 ) )
        # Tg.SetFillColor(   getattr(self, 'setFillColor',   self.color ) )
        Tg.SetFillColorAlpha(self.color, 0.30)

        if not(leg is None):
            leg.AddEntry( Tg.GetName(), self.title, 'LF' )

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
        Tg.SetName( GetUniqueRootName() )
    
        Tg.SetMarkerStyle( getattr(self, 'setMarkerStyle', 8 ) )
        Tg.SetMarkerSize(  getattr(self, 'setMarkerSize',  0 ) )
        Tg.SetMarkerColor( getattr(self, 'setMarkerColor', self.color ) )
        Tg.SetLineColor(   getattr(self, 'setLineColor',   self.color ) )

        # Tg.SetFillStyle(   getattr(self, 'setFillStyle',   3245 ) )
        # Tg.SetFillColor(   getattr(self, 'setFillColor',   self.color ) )
        Tg.SetFillColorAlpha(self.color, 0.30)

        if not(leg is None):
            leg.AddEntry( Tg.GetName(), self.title, 'LF' )

        return [ (Tg, 'E2PSAME') ]


    def repr_point_with_vertical_bar(self, leg=None):

        Tg = ROOT.TGraphAsymmErrors(
            self.n_bins,
            array( 'f', self.get_bin_centers() ),
            array( 'f', self.bin_values ),
            array( 'f', [ 0.0 for i in xrange(self.n_bins) ] ),
            array( 'f', [ 0.0 for i in xrange(self.n_bins) ] ),
            array( 'f', self.errs_down ),
            array( 'f', self.errs_up ),
            )
        ROOT.SetOwnership( Tg, False )
        Tg.SetName( GetUniqueRootName() )

        Tg.SetMarkerStyle( getattr(self, 'setMarkerStyle', 8 ) )
        Tg.SetFillColor(   getattr(self, 'setFillColor',   self.color ) )
        Tg.SetMarkerColor( getattr(self, 'setMarkerColor', self.color ) )
        Tg.SetLineColor(   getattr(self, 'setLineColor',   self.color ) )

        if not(leg is None):
            leg.AddEntry( Tg.GetName(), self.title, 'PE' )

        return [(Tg, 'PSAME')]


    def Draw(self, draw_style):
        for obj, draw_str in getattr(self, draw_style):
            obj.Draw(draw_str)
