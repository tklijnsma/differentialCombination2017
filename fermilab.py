from OptionHandler import flag_as_option, flag_as_parser_options

import os, logging, copy, sys, re
from array import array

import ROOT

import LatestPaths, LatestBinning
import differentials
import differentialutils

from differentials.plotting.canvas import c

import fermilabcode
from fermilabcode.read_canvas import CanvasReaderHzz, CanvasReaderHgg

#____________________________________________________________________
@flag_as_option
def read_hzz_from_canvas(args):
    cr = CanvasReaderHzz('fermilabcode/input/pT4l_unfoldwith_SM_125_logscale.root')
    cr.read()
    cr.dump('hzz_data')

@flag_as_option
def read_hgg_from_canvas(args):
    cr = CanvasReaderHgg('fermilabcode/input/diffXsec_pt_test_log.root')
    cr.read()
    cr.dump('hgg_data')


@flag_as_option
def naive_pt_combination(args):
    hgg = fermilabcode.minicombine.load_data('fermilabcode/input/hgg_data_Mar21.py')
    hgg.name = 'hgg'
    hzz = fermilabcode.minicombine.load_data('fermilabcode/input/hzz_data_Mar21.py')
    hzz.name = 'hzz'

    combination = fermilabcode.minicombine.PtCombination([hgg, hzz])
    combination.color = 1
    combination.plot_scans()

    hgg_hist = fermilabcode.minicombine.data_to_hist(hgg, title='hgg', color=2)
    hzz_hist = fermilabcode.minicombine.data_to_hist(hzz, title='hzz', color=4)

    combination_hist_hzz = combination.get_histogram_bestfit_hzz_binning()
    combination_hist_hzz.color = 3

    plot = differentials.plotting.plots.QuickPlot(
        'spectrum_ptcomb',
        x_min=0., x_max=500., y_min=-1., y_max=3.,
        )
    plot.add(hgg_hist, 'repr_point_with_vertical_bar_and_horizontal_bar')
    plot.add(hzz_hist, 'repr_point_with_vertical_bar_and_horizontal_bar')
    plot.add(combination_hist_hzz, 'repr_point_with_horizontal_bar')
    plot.add(combination.to_hist(), 'repr_point_with_vertical_bar_and_horizontal_bar')
    plot.draw()
    plot.wrapup()

    combination.dump('combination_data')


@flag_as_option
def naive_kappabkappac_combination(args):
    # hgg = fermilabcode.minicombine.load_data('fermilabcode/input/hgg_data_Mar21.py')
    # hgg.name = 'hgg'
    # hzz = fermilabcode.minicombine.load_data('fermilabcode/input/hzz_data_Mar21.py')
    # hzz.name = 'hzz'
    combination = fermilabcode.minicombine.load_data('fermilabcode/combination_data_Mar22.py')
    combination.name = 'combination'

    kappabkappac = fermilabcode.minicombine.CouplingFit()
    kappabkappac.setup_kappabkappac(combination)
    kappabkappac.get_scan()
    kappabkappac_hist = kappabkappac.to_hist()

    # Hier verder
    # Make plot prettier. Probably better to just start from c.Clear() instead of a plot script.
    plot = differentials.plotting.plots.QuickPlot(
        'kappabkappac',
        # x_min=0., x_max=500., y_min=-1., y_max=3.,
        )
    differentials.plotting.canvas.c.set_margins_2D()
    plot.add(kappabkappac_hist, 'repr_2D_with_contours')
    plot.draw()
    plot.wrapup()

    # # def plot(self, plotname, draw_style='repr_2D_with_contours'):

    # c.Clear()
    # c.set_margins_2D()

    # # leg = plotting.pywrappers.Legend(
    # #     c.GetLeftMargin() + 0.01,
    # #     c.GetBottomMargin() + 0.02,
    # #     1 - c.GetRightMargin() - 0.01,
    # #     c.GetBottomMargin() + 0.09
    # #     )
    # histogram2D = self.to_hist()
    # # histogram2D._legend = leg
    # histogram2D.Draw(draw_style)

    # base = histogram2D.H2
    # base.GetXaxis().SetTitle(self.x_title)
    # base.GetYaxis().SetTitle(self.y_title)
    # base.GetXaxis().SetTitleSize(0.06)
    # base.GetXaxis().SetLabelSize(0.05)
    # base.GetYaxis().SetTitleSize(0.06)
    # base.GetYaxis().SetLabelSize(0.05)

    # plotting.pywrappers.CMS_Latex_type().Draw()
    # plotting.pywrappers.CMS_Latex_lumi().Draw()
    # plotting.pywrappers.Point(self.x_sm, self.y_sm).Draw('repr_SM_point')

    # plotting.pywrappers.ContourDummyLegend(
    #     c.GetLeftMargin() + 0.01,
    #     1. - c.GetTopMargin() - 0.1,
    #     1. - c.GetRightMargin() - 0.01,
    #     1. - c.GetTopMargin() - 0.01,
    #     ).Draw()

    # c.Update()
    # c.RedrawAxis()
    # c.save(plotname)







@flag_as_option
def param_test(args):
    parametrization = get_kappab_kappac_parametrization()


def get_kappab_kappac_parametrization():
    coupling_variations = differentials.theory.theory_utils.FileFinder(
        muR=1.0, muF=1.0, Q=1.0, directory=LatestPaths.theory.yukawa.filedir
        ).get()
    sm = [ v for v in coupling_variations if v.kappab==1.0 and v.kappac==1.0 ][0]
    coupling_variations.pop(coupling_variations.index(sm))

    parametrization = differentials.parametrization.Parametrization()
    parametrization.parametrize_by_matrix_inversion = True
    parametrization.c1_name = 'kappab'
    parametrization.c2_name = 'kappac'
    parametrization.from_theory_dicts(coupling_variations)
    return parametrization