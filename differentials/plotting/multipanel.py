import ROOT
import plotting_utils as utils
import pywrappers
from canvas import c

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

        self.bottom_top_margin    = 0.00
        self.bottom_bottom_margin = 0.30
        self.bottom_left_margin   = 0.12
        self.bottom_right_margin  = 0.02

        self.top_log_scale = True
        self.CMS_labels = True


    def add_top(self, obj, draw_str=None):
        self.top_objects.append((obj, draw_str))

    def add_bottom(self, obj, draw_str=None):
        self.bottom_objects.append((obj, draw_str))

    def draw(self):
        c.Clear()

        _width = c.GetWindowWidth()
        _height = c.GetWindowHeight()
        c.SetCanvasSize( 800, 900 )

        toppad_bottom    = padSplitPoint
        bottompad_top    = padSplitPoint
        height_ratio = float( 1.0 - toppad_bottom ) / float( bottompad_top )

        toppad = ROOT.TPad(
            GetUniqueRootName(), '',
            0.0, toppad_bottom, 1.0, 1.0
            )
        toppad.SetBottomMargin( toppadBottomMargin) # Distance to the bottom panel
        toppad.SetTopMargin(    toppadTopMargin)     # Space for labels
        toppad.SetLeftMargin(   toppadLeftMargin)
        toppad.SetRightMargin(  toppadRightMargin)
        toppad.Draw()

        bottompad = ROOT.TPad(
            GetUniqueRootName(), '',
            0.0, 0.0, 1.0, bottompad_top
            )
        bottompad.SetTopMargin(    bottompadTopMargin)    # Distance to the bottom panel
        bottompad.SetBottomMargin( bottompadBottomMargin) # Space for labels
        bottompad.SetLeftMargin(   bottompadLeftMargin)
        bottompad.SetRightMargin(  bottompadRightMargin)
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

        # Draw the actual objects
        toppad.cd()
        if self.top_log_scale:
            toppad.SetLogy()
        for obj, drawStr in self.top_objects:
            obj.Draw(drawStr)

        bottompad.cd()
        for obj, drawStr in self.bottom_objects:
            obj.Draw(drawStr)


        # Process titles, ticks, labels
        toppad.cd()
        axisHolderTop = self.top_objects[0][0]
        axisHolderTop.GetXaxis().SetLabelOffset(999.)
        axisHolderTop.GetYaxis().SetTitle(self.y_title_top)

        bottompad.cd()
        axisHolderBottom = self.bottom_objects[0][0]
        # Set sizes of labels/ticks equal to top panel (undo automatic scaling by ROOT)
        axisHolderBottom.GetXaxis().SetLabelSize(axisHolderTop.GetXaxis().GetLabelSize() * height_ratio)
        axisHolderBottom.GetYaxis().SetLabelSize(axisHolderTop.GetYaxis().GetLabelSize() * height_ratio)
        axisHolderBottom.GetXaxis().SetTickLength(axisHolderTop.GetXaxis().GetTickLength() * height_ratio)

        axisHolderBottom.GetYaxis().SetTitle(self.y_title_bottom)
        axisHolderBottom.GetXaxis().SetTitle(self.x_title)
        axisHolderBottom.GetXaxis().SetTitleSize(axisHolderTop.GetXaxis().GetTitleSize() * height_ratio)
        axisHolderBottom.GetYaxis().SetTitleSize(axisHolderTop.GetYaxis().GetTitleSize() * height_ratio)
        axisHolderBottom.GetYaxis().SetTitleOffset(1./height_ratio)

        if self.CMS_labels:
            toppad.cd()
            pywrappers.CMS_Latex_type(text_size=0.08).Draw()
            pywrappers.CMS_Latex_lumi(text_size=0.07).Draw()

        c.save(self.name)
        c.SetCanvasSize( _width, _height )

