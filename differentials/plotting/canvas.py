import os.path

import ROOT
import plotting_utils as utils

from time import strftime
datestr = strftime( '%b%d' )


class Canvas(object):
    """Python wrapper for a ROOT TCanvas, that has a few extra functionalities"""

    save_pdf = True
    save_png = False
    save_root = False
    save_png_through_convert = False

    default_width = 1000
    default_height = 800

    def __init__(self):
        self.name = utils.get_unique_rootname()
        self.canvas = ROOT.TCanvas('ctc', 'ctc', self.default_width, self.default_height)
        self.plotdir = 'plots_{0}'.format(datestr)
        self._is_resized_temporarily = False

    def __getattr__(self, name):
        """
        Reroutes calls Canvas.xxxx to Canvas.canvas.xxxx
        This method should only be called if the attribute could not be found in Canvas
        """
        return getattr(self.canvas, name)

    def resize(self, width=None, height=None):
        if width is None:
            width = c.GetWindowWidth()
        if height is None:
            height = c.GetWindowHeight()
        self.canvas.SetCanvasSize(width, height)

    def resize_temporarily(self, width=None, height=None):
        self._tmp_width = c.GetWindowWidth()
        self._tmp_height = c.GetWindowHeight()
        self.resize(width, height)
        self._is_resized_temporarily = True

    def set_margins(
            self,
            LeftMargin   = 0.15,
            RightMargin  = 0.03,
            BottomMargin = 0.15,
            TopMargin    = 0.09,
            ):
        self.canvas.SetLeftMargin( LeftMargin )
        self.canvas.SetRightMargin( RightMargin )
        self.canvas.SetBottomMargin( BottomMargin )
        self.canvas.SetTopMargin( TopMargin )

    def set_margins_2D(
            self,
            LeftMargin   = 0.12,
            RightMargin  = 0.10,
            BottomMargin = 0.12,
            TopMargin    = 0.09, # Maybe 0.08
            ):
        self.canvas.SetLeftMargin( LeftMargin )
        self.canvas.SetRightMargin( RightMargin )
        self.canvas.SetBottomMargin( BottomMargin )
        self.canvas.SetTopMargin( TopMargin )

    def set_plotdir(self, newdir):
        self.plotdir = newdir

    def save(self, outname, pdf=True, png=False, root=False, png_through_convert=False ):
        # Check if a '/' was passed in the outname; if so, create a sub directory
        # Only allow 1 additional slash (otherwise may accidentally create deep tree structures)
        outdir = self.plotdir
        subdir = ''
        if len(outname.split('/')) == 2:
            subdir = outname.split('/')[0]
            outdir = os.path.join(outdir, subdir)

        if not os.path.isdir(outdir): os.makedirs(outdir)
        outname = os.path.join(outdir, os.path.basename(outname).replace('.pdf','').replace('.png',''))

        if (pdf or self.save_pdf) or (png_through_convert or self.save_png_through_convert):
            self.canvas.SaveAs(outname+'.pdf')
        if png or self.save_png:
            self.canvas.SaveAs(outname+'.png')
        if root or self.save_root:
            self.canvas.SaveAs(outname+'.root')
        if png_through_convert or self.save_png_through_convert:
            # See: https://stackoverflow.com/a/6605085/9209944
            cmd = 'convert -density 300 -quality 100 {0}.pdf -trim {0}.png'.format(outname)
            os.system(cmd)

        if self._is_resized_temporarily:
            self.resize(self._tmp_width, self._tmp_height)

# Create one instance that can be called from anywhere
c = Canvas()
global_color_cycle = utils.new_color_cycle()
