from OptionHandler import flag_as_option, flag_as_parser_options

import os, logging, copy, sys, re, numpy
from array import array
from math import sqrt

import ROOT

import LatestPaths, LatestBinning
import differentials
import differentialutils

from differentials.plotting.canvas import c








#____________________________________________________________________
class WSParametrizationCtCb(differentials.parametrization.WSParametrization):

    def __init__(self, ws):
        super(WSParametrizationCtCb, self).__init__(ws)
        # Set the sm xs
        self.get_smxs_from_ws()
        self.set_smxs_xH(LatestBinning.obs_pth_xH.crosssection())
        self.exp_binning = [ 0., 15., 30., 45., 80., 120., 200., 350., 600., 800. ]
        self.select_decay_channel('hgg')
        self.tag = ''

    def get_histogram(self, ys):
        return differentials.plotting.pywrappers.Histogram(
            differentials.plotting.plotting_utils.get_unique_rootname(),
            '',
            self.exp_binning,
            ys,
            )

    def get_hist_ggH(self, ct, cb):
        xs = self.get_xs_exp_ggH(ct=ct, cb=cb)
        hist = self.get_histogram(xs)
        hist.title = 'ggH (#kappa_{{t}} = {0}, #kappa_{{b}} = {1})'.format(ct, cb)
        if len(self.tag) > 0: hist.title = self.tag + ':' + hist.title
        return hist

    def get_hist_xH(self, ct, cb):
        xs = self.get_xs_exp_xH(ct=ct, cb=cb)
        hist = self.get_histogram(xs)
        hist.title = 'xH (#kappa_{{t}} = {0}, #kappa_{{b}} = {1})'.format(ct, cb)
        if len(self.tag) > 0: hist.title = self.tag + ':' + hist.title
        return hist

    def get_hist_smH(self, ct, cb):
        xs_ggH = self.get_xs_exp_ggH(ct=ct, cb=cb)
        xs_xH = self.get_xs_exp_xH(ct=ct, cb=cb)
        xs = [ x1+x2 for x1, x2 in zip(xs_ggH, xs_xH) ]
        hist = self.get_histogram(xs)
        hist.title = 'smH (#kappa_{{t}} = {0}, #kappa_{{b}} = {1})'.format(ct, cb)
        if len(self.tag) > 0: hist.title = self.tag + ':' + hist.title
        return hist


@flag_as_option
def check_ctcb_xHBRscaling(args):

    ws_noBRscaling = 'out/workspaces_May28/combWithHbb_TopCtCb_reweighted_scalingbbHttH_couplingdependentBRs.root'
    ws_xHBRscaling = 'out/workspaces_May29/combWithHbb_TopCtCb_reweighted_scalingbbHttH_couplingdependentBRs.root'

    noBRscaling = WSParametrizationCtCb(ws_noBRscaling)
    xHBRscaling = WSParametrizationCtCb(ws_xHBRscaling)
    xHBRscaling.tag = 'xH_scaling'

    points = [
        ( 1.0, 1.0 ),
        ( 0.5, 0.0 ),
        ( 1.65, 1.4 ),
        ]

    plotname = 'xcheck_xHBRscaling'
    plot = differentials.plotting.plots.QuickPlot(
        plotname,
        x_min = 0., x_max= 350., y_min=0., y_max=20.
        )
    plot.x_title = 'pT'
    plot.y_title = '#sigma'

    for ct, cb in points:
        plot.clear()
        colors = iter([ 2, 4, 3, 8, 9, 42 ])
        plot.plotname = plotname + '_ct{0}_cb{1}'.format(differentials.core.float_to_str(ct), differentials.core.float_to_str(cb))

        color = next(colors)
        smH_noBRscaling = noBRscaling.get_hist_smH(ct, cb)
        smH_xHBRscaling = xHBRscaling.get_hist_smH(ct, cb)
        smH_noBRscaling.color = color
        smH_xHBRscaling.color = color

        color = next(colors)
        ggH_noBRscaling = noBRscaling.get_hist_ggH(ct, cb)
        ggH_xHBRscaling = xHBRscaling.get_hist_ggH(ct, cb)
        ggH_noBRscaling.color = color
        ggH_xHBRscaling.color = color

        color = next(colors)
        xH_noBRscaling = noBRscaling.get_hist_xH(ct, cb)
        xH_xHBRscaling = xHBRscaling.get_hist_xH(ct, cb)
        xH_noBRscaling.color = color
        xH_xHBRscaling.color = color

        plot.add(smH_noBRscaling, 'repr_basic_histogram')
        plot.add(smH_xHBRscaling, 'repr_basic_dashed_histogram')

        plot.add(ggH_noBRscaling, 'repr_basic_histogram')
        plot.add(ggH_xHBRscaling, 'repr_basic_dashed_histogram')

        plot.add(xH_noBRscaling, 'repr_basic_histogram')
        plot.add(xH_xHBRscaling, 'repr_basic_dashed_histogram')

        plot.draw()
        plot.wrapup()



persistence = []
new_tcolor_index = 5001
def make_lighter_color(color):
    tcolor = ROOT.gROOT.GetColor(color)
    r = tcolor.GetRed()
    b = tcolor.GetBlue()
    g = tcolor.GetGreen()
    new_r = min(1.0, r + 4./256.)
    new_b = min(1.0, b + 4./256.)
    new_g = min(1.0, g + 4./256.)

    print 'old r: {0}, new r: {1}'.format(r, new_r)
    print 'old b: {0}, new b: {1}'.format(b, new_b)
    print 'old g: {0}, new g: {1}'.format(g, new_g)

    global new_tcolor_index
    new_tcolor_index += 1
    new_color_index = new_tcolor_index

    new_tcolor = ROOT.TColor(new_color_index, new_r, new_b, new_g)
    ROOT.SetOwnership(new_tcolor, False)
    persistence.append(new_tcolor)
    return new_tcolor.GetNumber()

#____________________________________________________________________
class ShapeGetterFromWS(object):
    """Wrapper around WSParametrization"""
    def __init__(self, ws):
        super(ShapeGetterFromWS, self).__init__()
        self.ws = ws
        
        self.parametrization = differentials.parametrization.WSParametrization(ws)
        self.smxs = [ roovar.getVal() for roovar in differentials.core.read_set(ws, 'SMXS', return_names=False) ]
        self.parametrization.set_smxs(self.smxs)
        self.exp_binning = [ 0., 15., 30., 45., 80., 120. ]

    def get_shape_histogram(self, kb, kc):
        shape = self.parametrization.get_shape_exp(kappab=kb, kappac=kc)
        histogram = differentials.plotting.pywrappers.Histogram(
            differentials.plotting.plotting_utils.get_unique_rootname(),
            '#kappa_{{b}} = {0}, #kappa_{{c}} = {1}'.format(kb, kc),
            self.exp_binning,
            shape,
            # color = 2
            )
        return histogram


@flag_as_option
def check_shapes(args):

    ggH = ShapeGetterFromWS('out/workspaces_May24/combination_Yukawa_reweighted_floatingBRs.root')
    ggH_plus_bbH = ShapeGetterFromWS('out/workspaces_May22/combination_Yukawa_reweighted_scalingbbH_floatingBRs.root')

    # plot_parametrizations_in_workspace_for_some_points(
    #     ggH_plus_bbH,
    #     plotname = 'xcheck_shapes_ggH_plus_bbH'
    #     )

    # plot_parametrizations_in_workspace_for_some_points(
    #     ggH,
    #     plotname = 'xcheck_shapes_ggH'
    #     )
    
    # shapes for some specific points for both
    points = [
        ( 1.0, 1.0 ),
        ( 40.0, 1.0 ),
        ( -40.0, 1.0 ),
        ]

    plot = differentials.plotting.plots.QuickPlot(
        'xcheck_shapes_directcompare',
        x_min = 0., x_max= 120., y_min=0., y_max=0.55
        )
    plot.x_title = 'pT'
    plot.y_title = 'Shape'

    colors = iter([ 2, 4, 3, 8, 9, 42 ])
    for kb, kc in points:
        color = next(colors)

        hist_ggH = ggH.get_shape_histogram(kb, kc)
        hist_ggH.color = color

        hist_ggH_plus_bbH = ggH_plus_bbH.get_shape_histogram(kb, kc)
        hist_ggH_plus_bbH.color = color
        hist_ggH_plus_bbH.title += ' (incl. bbH)'
        hist_ggH_plus_bbH.line_width = 4

        plot.add(hist_ggH, 'repr_basic_histogram')
        plot.add(hist_ggH_plus_bbH, 'repr_basic_dashed_histogram')

    plot.draw()
    plot.wrapup()




def plot_parametrizations_in_workspace_for_some_points(
    shapegetter,
    plotname='xcheck_shapes'
    ):

    points = [
        ( 1., 1. ),
        # 
        ( 10., 1. ),
        ( 10., 10. ),
        ( 1., 10. ),
        ( -10., 1. ),
        ( -10., -10. ),
        ( 1., -10. ),
        ( 10., -10. ),
        ( -10., 10. ),
        # 
        ( 10000., 1. ),
        ( 10000., 10000. ),
        ( 1., 10000. ),
        ( -10000., 1. ),
        ( -10000., -10000. ),
        ( 1., -10000. ),
        ( 10000., -10000. ),
        ( -10000., 10000. ),
        ]

    histograms = [ shapegetter.get_shape_histogram(kb, kc) for kb, kc in points ]

    plot = differentials.plotting.plots.QuickPlot(
        plotname,
        x_min = 0., x_max= 120., y_min=0., y_max=1.0
        )
    plot.x_title = 'pT'
    plot.y_title = 'Shape'

    for histogram in histograms:
        plot.add(histogram, 'repr_basic_histogram')

    plot.draw()
    plot.wrapup()


#____________________________________________________________________
@flag_as_option
def check_fiducial_vs_inclusive_shape(args):
    # LatestBinning.shape_pth_smH
    # LatestBinning.shape_pth_ggH
    # LatestBinning.shape_pth_xH

    auc = AcceptanceUncertaintyCalculator('suppliedInput/fromVittorio/scaleWeightVariationFullPhaseSpaceCombination_Pt.npz')
    smH = auc.get_central_shape()

    print LatestBinning.shape_pth_smH
    print LatestBinning.shape_pth_ggH



#____________________________________________________________________
@flag_as_option
def xH_systematic_uncertainty(args):
    
    ggH = LatestBinning.obs_pth_ggH
    xH  = LatestBinning.obs_pth_xH

    ggH.Print()

    print

    xH.Print()

    print 'inclusive xH (pb):', LatestBinning.YR4_xH
    print 'Unc on inclusive xH (pb): {0} (fraction: {1:.4f} )'.format(
        LatestBinning.xH_unc_inclusive, LatestBinning.xH_unc_inclusive_fraction
        )


#____________________________________________________________________

class CouplingDependenceOfFunction(object):
    """docstring for CouplingDependenceOfFunction"""
    def __init__(self):
        super(CouplingDependenceOfFunction, self).__init__()
        self.n_points = 100

    def get_x_y(self, x_var, y_func, x_min=-10., x_max=10.):
        x_axis = self.get_axis(self.n_points, x_min, x_max)
        y_axis = []
        init_val_x_var = x_var.getVal()
        for x in x_axis:
            x_var.setVal(x)
            y_axis.append(y_func.getVal())
        x_var.setVal(init_val_x_var)
        return x_axis, y_axis

    def get_axis(self, n_points, x_min, x_max):
        dx = (x_max-x_min) / (n_points-1)
        return [ x_min + i*dx for i in xrange(n_points) ]

    def to_graph(self, x_axis, y_axis, title=''):
        g = differentials.plotting.pywrappers.Graph(
            differentials.plotting.plotting_utils.get_unique_rootname(),
            title,
            x_axis, y_axis
            )
        return g

    def get_x_y_graph(self, x_var, y_func, x_min=-10., x_max=10.):
        x_axis, y_axis = self.get_x_y(x_var, y_func, x_min, x_max)
        return self.to_graph(x_axis, y_axis, '{0}({1})'.format(y_func.GetTitle(), x_var.GetTitle()) )


@flag_as_option
def plot_mu_parabolas(args):
    ws = 'out/workspaces_May16/combination_Yukawa_reweighted_G1A_noTheoryUnc_scaledByMuTotalXS.root'
    w = differentials.core.get_ws(ws)

    mu_inc_unreweighted = w.function('totalXSmodifier')
    kb = w.var('kappab')
    kc = w.var('kappac')


    calc = CouplingDependenceOfFunction()
    g_kb = calc.get_x_y_graph(kb, mu_inc_unreweighted)
    g_kc = calc.get_x_y_graph(kc, mu_inc_unreweighted)

    plot = differentials.plotting.plots.QuickPlot(
        'xy_' + mu_inc_unreweighted.GetTitle(),
        x_min = -10., x_max = 10., y_min = 1., y_max = 3.
        )
    plot.x_title = '#kappa_{x}'
    plot.y_title = mu_inc_unreweighted.GetTitle()
    plot.add(g_kb, 'repr_basic_line')
    plot.add(g_kc, 'repr_basic_line')
    plot.draw()
    plot.wrapup()


@flag_as_option
def crosscheck_htt_kappat_scaling(args):
    exp_binning = LatestBinning.binning_pth
    n_bins = len(exp_binning)-1

    coupling_variations = differentials.theory.theory_utils.FileFinder(
        cb=1.0, muR=1.0, muF=1.0, Q=1.0, directory=LatestPaths.theory.top.filedir
        ).get()
    sm = [ v for v in coupling_variations if v.ct==1.0 and v.cg==0.0 ][0]
    coupling_variations.pop(coupling_variations.index(sm))

    parametrization = differentials.parametrization.Parametrization()
    parametrization.parametrize_by_matrix_inversion = True
    parametrization.do_linear_terms = False
    parametrization.c1_name = 'ct'
    parametrization.c2_name = 'cg'
    parametrization.c2_SM   = 0.0
    for v in coupling_variations[:6]:
        parametrization.add_variation(v.ct, v.cg, v.crosssection)
    parametrization.parametrize()
    parametrization.make_rebinner(sm.binBoundaries, exp_binning)


    # Assume ttH shape is the same
    # fraction_ttH_xH = LatestBinning.YR4_ttH / LatestBinning.YR4_xH
    # smxs_ttH = [ fraction_ttH_xH * xs for xs in parametrization.evaluate(1.0, 0.0) ]
    # smxs_ttH = [ s * LatestBinning.YR4_ttH for s in LatestBinning.shape_pth_xH ]
    smxs_ttH = LatestBinning.obs_pth_ttH.crosssection()

    # # Check at SM
    # smxs_ggH = parametrization.get_xs_exp_integrated_per_bin(1., 0.)
    # print sum(smxs_ttH), LatestBinning.YR4_ttH
    # print sum(smxs_ggH), LatestBinning.YR4_ggF_n3lo
    # for i in xrange(n_bins):
    #     r = (
    #         'bin {0}, left = {1:<6.1f}, smxs_ttH = {2:<7.4f}, smxs_ggH = {3:<7.4f}'
    #         .format(
    #             i, exp_binning[i], smxs_ttH[i], smxs_ggH[i]
    #             )
    #         )
    #     print r
    # sys.exit()

    
    # Function that calculates the fraction ttH / ggH per bin
    def fraction_ttH(ct, cg, verbose=True):
        r = []
        xs_ggH = parametrization.get_xs_exp_integrated_per_bin(ct, cg)
        if verbose: print '\nct = {0:5.2f}, cg={1:5.2f}'.format(ct, cg)
        for i in xrange(n_bins):
            ttH = ct*ct * smxs_ttH[i]
            ggH = xs_ggH[i]
            r.append(ttH/ggH)
            if verbose:
                print 'Bin {0}: ttH = {1:<8.4f}, ggH = {2:<8.4f}, fraction = {3}'.format(
                    i, ttH, ggH, ttH/ggH
                    )
        return r

    # observed:
    # 0.6 < ct < 3.4
    # -0.2 < cg < 0.04
    ct_min = 0.6
    ct_max = 3.4
    cg_min = -0.2
    cg_max = 0.04

    points = [
        ( 1.0, 0.0 ),
        ( ct_max, 0.0 ),
        ( ct_min, cg_max ),
        ( ct_max, cg_min )
        ]

    histograms = []
    for ct, cg in points:
        histogram = differentials.plotting.pywrappers.Histogram(
            differentials.plotting.plotting_utils.get_unique_rootname(),
            '#kappa_{{t}} = {0}, c_{{g}} = {1}'.format(ct, cg),
            exp_binning,
            fraction_ttH(ct, cg),
            # color = 2
            )
        histograms.append(histogram)

    y_max = 1.25 * max([ H.y_max() for H in histograms ])

    plot = differentials.plotting.plots.QuickPlot(
        'ttH_fractions',
        x_min = exp_binning[0], x_max=800.,
        y_min = 0.0,
        # y_max = 2.0,
        y_max = y_max,
        )
    plot.x_title = 'p_{T} (GeV)'
    plot.y_title = '(#kappa_{t}^{2} #sigma_{ttH}^{SM}) / #sigma_{ggH}(#kappa_{t},c_{g})'

    for H in histograms:
        plot.add(H, 'repr_basic_histogram')

    plot.draw()
    plot.wrapup()


@flag_as_option
def crosscheck_hbb_kappab_scaling(args):
    coupling_variations = differentials.theory.theory_utils.FileFinder(
        muR=1.0, muF=1.0, Q=1.0, directory=LatestPaths.theory.yukawa.filedir
        ).get()

    sm = [ v for v in coupling_variations if v.kappab==1.0 and v.kappac==1.0 ][0]
    coupling_variations.pop(coupling_variations.index(sm))

    parametrization = differentials.parametrization.Parametrization()
    parametrization.parametrize_by_matrix_inversion = True
    parametrization.c1_name = 'kappab'
    parametrization.c2_name = 'kappac'
    for v in coupling_variations[:6]:
        parametrization.add_variation(v.kappab, v.kappac, v.crosssection)
    parametrization.parametrize()

    # Assume bbH shape is the same
    fraction_bbH_SM = LatestBinning.YR4_bbH / LatestBinning.YR4_ggF_n3lo
    smxs_bbH = [ fraction_bbH_SM * xs for xs in parametrization.evaluate(1.0, 1.0) ]

    n_bins = len(sm.crosssection)
    # Function that calculates the fraction bbH / ggH per bin
    def fraction_bbH(kappab, kappac):
        r = []
        xs_ggH = parametrization.evaluate(kappab, kappac)
        for i in xrange(n_bins):
            r.append(
                kappab*kappab * smxs_bbH[i] / xs_ggH[i]
                )
        return r

    # ( -2.1 < kb < 3.8 expected) 
    # ( -8.6 < kc < 10.5 expected)
    kb_min = -2.1
    kb_max = 3.8
    kc_min = -8.6
    kc_max = 10.5

    points = [
        ( 1.0, 1.0 ),
        ( 1.0, kc_min ),
        ( 1.0, kc_max ),
        ( kb_min, 1.0 ),
        ( kb_max, 1.0 ),
        ( kb_min, kc_min ),
        ( kb_max, kc_max ),
        ]

    histograms = []
    for kappab, kappac in points:
        histogram = differentials.plotting.pywrappers.Histogram(
            differentials.plotting.plotting_utils.get_unique_rootname(),
            '#kappa_{{b}} = {0}, #kappa_{{c}} = {1}'.format(kappab, kappac),
            sm.binBoundaries,
            fraction_bbH(kappab, kappac),
            # color = 2
            )
        histograms.append(histogram)

    y_max = 1.25 * max([ H.y_max() for H in histograms ])

    plot = differentials.plotting.plots.QuickPlot(
        'bbH_fractions',
        x_min = sm.binBoundaries[0], x_max=sm.binBoundaries[-1], y_min=0.0, y_max = y_max
        )
    plot.x_title = 'p_{T} (GeV)'
    plot.y_title = '(#kappa_{b}^{2} #sigma_{bbH}^{SM}) / #sigma_{ggH}(#kappa_{b},#kappa_{c})'

    for H in histograms:
        plot.add(H, 'repr_basic_histogram')

    plot.draw()
    plot.wrapup()


@flag_as_option
def theory_uncertainty_on_inclusive_xs_yukawa(args):
    get_inc_xs_uncertainty(LatestPaths.theory.yukawa.filedir_gluoninduced)
    get_inc_xs_uncertainty(LatestPaths.theory.yukawa.filedir)


def get_inc_xs_uncertainty(theory_filedir):
    print '\nComputing inc xs uncertainty for {0}'.format(theory_filedir)

    coupling_variations = differentials.theory.theory_utils.FileFinder(
        kappab=1.0, kappac=1.0, directory=theory_filedir
        ).get()

    sm = [ v for v in coupling_variations if v.muR==1.0 and v.muF==1.0 and v.Q==1.0 ][0]
    coupling_variations.pop(coupling_variations.index(sm))

    def print_var(v, is_sm=False):
        r = '  muR = {0}, muF = {1}, Q = {2}:  {3}'.format(v.muR, v.muF, v.Q, v.inc_xs)
        if is_sm: r += ' (SM)'
        print r

    print 'Inc xs per scale variation:'
    inc_xs = []
    sm['inc_xs'] = sum(sm['crosssection_integrated'])
    print_var(sm, is_sm=True)
    for v in coupling_variations:
        v['inc_xs'] = sum(v['crosssection_integrated'])
        inc_xs.append(v.inc_xs)
        print_var(v)

    down = (sm.inc_xs - min(inc_xs)) / sm.inc_xs
    up   = (max(inc_xs) - sm.inc_xs) / sm.inc_xs
    symm = 0.5*(abs(up)+abs(down))

    print 'Unc down: {0}'.format(down)
    print 'Unc up  : {0}'.format(up)
    print 'Unc symm: {0}'.format(symm)


    print '  Double check for sm:'
    n_bins = len(sm.binBoundaries)-1
    smxs_inc = 0.0
    for i in xrange(n_bins):
        smxs_inc += sm.crosssection[i] * (sm.binBoundaries[i+1] - sm.binBoundaries[i])
    print '  smxs_inc = {0} (from xs/GeV * bin widths)'.format(smxs_inc)
    print '  smxs_inc = {0} (from sum(crosssection_integrated))'.format(sum(sm.crosssection_integrated))


class AcceptanceUncertaintyCalculator(object):
    """docstring for AcceptanceUncertaintyCalculator"""
    def __init__(self, npz_file):
        super(AcceptanceUncertaintyCalculator, self).__init__()
        self.npz_file = npz_file
        self.npz = numpy.load(self.npz_file)
        self.reinit_lists()

    def reinit_lists(self):
        self.central = []
        self.up = []
        self.down = []
        self.symm = []

    def get_index(self, binstr):
        return int(re.match(r'bin(\d+)', binstr).group(1))

    def get_keys(self):
        r = self.npz.keys()
        r.sort(key=lambda k: self.get_index(k))
        return r

    def binstr_to_proper_str(self, binstr):
        index = self.get_index(binstr)
        boundaries = [ 0., 15., 30., 45., 80., 120., 200., 350., 600. ]
        if index < 8:
            return '[{0:0d},{1:0d})'.format(boundaries[index], boundaries[index+1])
        else:
            return '[{0:0d},#infty)'.format(boundaries[index])

    def get_unc_for_bin(self, A):
        B = list(A[:].flatten())
        logging.debug(B)
        B.pop(7) # Remove ratio 4 scale variations (higher index first to not mess up the next pop)
        B.pop(5) # Remove ratio 4 scale variations
        central = B[0]
        up   = abs(max(B))/central - 1
        down = 1 - abs(min(B))/central
        symm = 0.5*(abs(down)+abs(up))
        return central, down, up, symm

    def get_cud(self):
        self.reinit_lists()
        for k in self.get_keys():
            logging.debug('Doing key {0}'.format(k))
            central, down, up, symm = self.get_unc_for_bin(self.npz[k])
            self.central.append(central)
            self.up.append(up)
            self.down.append(down)
            self.symm.append(symm)
        
    def get_central_shape(self):
        self.get_cud()
        return self.central

@flag_as_option
def check_acceptance_uncertainties(args):
    # Get acceptance uncertainties
    auc = AcceptanceUncertaintyCalculator('suppliedInput/fromVittorio/scaleWeightVariationFullPhaseSpaceCombination_Pt.npz')
    auc.get_cud()
    logging.info('Found following set of acceptance uncertainties per bin: {0}'.format(auc.symm))

    obstuple = LatestBinning.obstuple_pth_ggH
    scandict = LatestPaths.scan.pth_ggH.asimov if args.asimov else LatestPaths.scan.pth_ggH.observed

    systshapemaker = differentials.systshapemaker.SystShapeMaker()
    systshapemaker.set_sm(obstuple.combWithHbb.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=False))

    combWithHbb = differentials.scans.DifferentialSpectrum('combWithHbb', scandict.combWithHbb)
    combWithHbb.set_sm(obstuple.combWithHbb.crosssection_over_binwidth(normalize_by_second_to_last_bin_width=True))
    combWithHbb.read()

    systonly_histogram, systonly_histogram_xs = systshapemaker.get_systonly_histogram(combWithHbb, scandict.combWithHbb_statonly)

    # get symmetrized systematic uncertainty from histogram
    symm = [ 0.5*(abs(up)+abs(down)) for down, up in zip(systonly_histogram.errs_down, systonly_histogram.errs_up) ]

    logging.info('acceptance uncertainties: {0}'.format(auc.symm))
    logging.info('syst error before adding: {0}'.format(symm))

    symm_plus_AU = [ sqrt(s1**2+s2**2) for s1, s2 in zip(auc.symm, symm) ]
    change_in_syst = [ after/before-1. for before, after in zip(symm, symm_plus_AU) ]

    logging.info('syst error after adding:  {0}'.format(symm_plus_AU))
    logging.info('syst error rel change:    {0}'.format(change_in_syst))


    total_histogram = combWithHbb.to_hist()
    total = [ 0.5*(abs(up)+abs(down)) for down, up in zip(total_histogram.errs_down, total_histogram.errs_up) ]
    total_plus_AU = [ sqrt(s1**2+s2**2) for s1, s2 in zip(auc.symm, total) ]
    change_in_total = [ after/before-1. for before, after in zip(total, total_plus_AU) ]


    bounds = [ '0', '15', '30', '45', '80', '120', '200', '350', '600', 'infinity' ]
    bounds_line = [ '[{0},{1})'.format(left, right) for left, right in zip(bounds[:-1], bounds[1:]) ]

    table_py = [
        ['Bins'] + bounds_line,
        ['Acc. uncertainties'] + auc.symm,
        ['Rel. change in syst. unc.'] + change_in_syst,
        ['Rel. change in tot. unc.'] + change_in_total,
        ]
    table = differentials.plotting.tables.Table()
    table.from_list(table_py)

    formatter = differentials.plotting.tables.Formatter(
        include_sign = False,
        n_decimals = 1,
        is_percentage = True,
        )
    table.formatter = formatter
    table.max_col_width = 28

    print table.repr_terminal()
    print table.repr_twiki()


@flag_as_option
def check_reweighting_for_inc_xs(args):
    ws = LatestPaths.ws.yukawa.nominal.combination
    w = differentials.core.get_ws(ws)

    reweightor_names = [
        'reweightor_ggH_PTH_0_15',
        'reweightor_ggH_PTH_15_30',
        'reweightor_ggH_PTH_30_45',
        'reweightor_ggH_PTH_45_80',
        'reweightor_ggH_PTH_80_120',
        ]

    weights = []
    for name in reweightor_names:
        weights.append(w.function(name).getVal())

    smxs_Pier = [ 13.0137626356, 12.3785279848, 6.87788869197, 7.08329398674, 3.07720573252 ]
    smxs_Vitt = [ 12.158274514, 12.6947320234, 8.0889999211, 9.15091631806, 3.78165937525 ]

    mu = sum(smxs_Pier) / sum(smxs_Vitt)

    print 'Overall mu: {0}'.format(mu)
    print '  sum(smxs) Pier (up to 120 GeV): {0}'.format(sum(smxs_Pier))
    print '  sum(smxs) Vitt (up to 120 GeV): {0}'.format(sum(smxs_Vitt))



class PDFDrawer(object):
    """docstring for PDFDrawer"""
    def __init__(self, ws):
        super(PDFDrawer, self).__init__()
        self.ws = ws
        self.w = differentials.core.get_ws(self.ws)

        self.MH = self.w.var('MH')
        self.MH.setVal(125.)

        self.mH_hgg   = self.w.var('CMS_hgg_mass')
        ROOT.SetOwnership( self.mH_hgg, False )
        self.default_mH_hgg_range = [ self.mH_hgg.getMin(), self.mH_hgg.getMax() ]
        self.default_mH_hgg_nBins = self.mH_hgg.getBins()

        self.all_pdfs = []
        all_pdfs_arglist = ROOT.RooArgList(self.w.allPdfs())
        for i_pdf in xrange(all_pdfs_arglist.getSize()):
            self.all_pdfs.append( all_pdfs_arglist[i_pdf].GetName() )
        self.all_pdfs.sort()

        self.bin_width = 0.25

        self.hgg_cat = self.w.cat('CMS_channel')
        ROOT.SetOwnership(self.hgg_cat, False )


    def set_testable_yield_parameter(self, name):
        self.mu = self.w.var(name)
        if self.mu == None:
            raise ValueError
        self.mu.setRange(-100., 100.)

    def get_histogram(self, pdf_name):
        mH = self.mH_hgg
        default_mH_range = self.default_mH_hgg_range

        # Get bkg and sig pdfs
        mH.setRange( *default_mH_range )
        mH.setBins( int((default_mH_range[1]-default_mH_range[0])/self.bin_width) )

        pdf = self.w.pdf(pdf_name)
        if pdf == None:
            raise ValueError('Pdf {0} does not exist in {1}'.format(pdf_name, self.ws))

        H = pdf.createHistogram('H_'+pdf_name, mH)
        # Hbkg.Scale(self.bin_width) # default is /GeV
        ROOT.SetOwnership(H, False)

        return H

    def set_category(self, some_string):
        i_cat = int(re.search(r'SigmaMpTTag_(\d)', some_string).group(1))
        # cat_str = 'SigmaMpTTag_{0}_recoPt_600p0to10000p0_13TeV'.format(i_cat)
        # cat_str = 'recoPt_600p0_10000p0_SigmaMpTTag_0_13TeV'
        cat_str = 'recoPt_0p0_15p0_SigmaMpTTag_0_13TeV'

        # self.hgg_cat.setIndex(i_cat)
        self.hgg_cat.setLabel(cat_str)

    def get_data_histogram(self, pdf_name):
        self.set_category(pdf_name)
        dataset = self.w.data('data_obs')
        ROOT.SetOwnership( dataset, False )
        Hdata = self.roodataset_to_hist(dataset, self.mH_hgg, self.hgg_cat)
        return Hdata

    def roodataset_to_hist(self, dataset, xVar, category):
        ctemp = ROOT.TCanvas( 'ctemp', 'ctemp', 1000, 800 )
        ctemp.cd()
        ctemp.Clear()

        frame = xVar.frame()

        reduceStr = '{0}=={1}'.format( category.GetName(), category.getIndex() )

        dataset_reduced = dataset.reduce( reduceStr )
        dataset_reduced.plotOn( frame )
        frame.Draw()

        l = ctemp.GetListOfPrimitives()
        for i in xrange(l.GetEntries()):
            if isinstance( l.At(i), ROOT.RooHist ):
                H = l.At(i)
                break
        else:
            raise ValueError( 'ERROR: did not find a histogram', throwException=True )

        Hcopy = ROOT.RooHist( H )
        ROOT.SetOwnership( Hcopy, False )

        Hcopy.SetName( differentials.plotting.plotting_utils.get_unique_rootname() )

        # ctemp.SaveAs( 'plots_{0}_onetimeplots/roodatasetplottest.pdf'.format(datestr) )
        del ctemp
        del frame

        c.cd()

        return Hcopy

        

@flag_as_option
def hgg_pdf_xcheck(args):

    ws_smH = 'out/workspaces_Mar14/ws_pth_smH_hgg.root'
    ws_ggH = 'out/workspaces_Mar02/ws_pth_ggH_hgg.root'

    yp = {
        ws_smH : 'r_smH_PTH_GT600',
        ws_ggH : 'r_ggH_PTH_GT600'
        }

    for ws in [ ws_smH, ws_ggH ]:

        pdfdrawer = PDFDrawer(ws)
        pdfdrawer.set_testable_yield_parameter(yp[ws])

        c.Clear()
        c.set_margins()

        base = differentials.plotting.pywrappers.Base(
                x_min = 100.,
                x_max = 180.,
                y_min = 0.,
                # y_max = 0.20,
                y_max = 1.20,
                x_title = 'mH',
                y_title = '',
                )
        base.Draw()

        pdf_name = 'pdf_binrecoPt_600p0_10000p0_SigmaMpTTag_0_13TeV'

        for mu_val in [ 0.5, 1.0, 2.0, -1.0, -2.0, -3.0, -5.0, -10.0 ]:
            pdfdrawer.mu.setVal(mu_val)
            # H = pdfdrawer.get_histogram('pdf_binrecoPt_600p0_10000p0_SigmaMpTTag_0_13TeV_obsOnly')
            H = pdfdrawer.get_histogram(pdf_name)
            H.Draw('HISTSAME')

        H_data = pdfdrawer.get_data_histogram(pdf_name)
        H_data.Draw('HISTSAME')

        outname = 'hgg_pdf_xcheck_' + os.path.basename(ws).replace('.root','')
        c.save(outname)



#____________________________________________________________________
class Spectrum(object):
    """docstring for Spectrum"""
    def __init__(self, name, values, isSM, title=None, color=4, line_style=1):
        self.values = values
        self.name = name
        self.isSM = isSM
        self.color = color
        self.line_style = line_style
        if title is None:
            self.title = name
        else:
            self.title = title

class NormalizationCrossCheck(object):
    """docstring for NormalizationCrossCheck"""
    color_cycle = differentials.plotting.canvas.global_color_cycle
 
    def __init__(self):
        self.logscale = True
        self.last_bin_is_overflow = False
        self.spectra = []

    def add_spectrum(self, name, values, isSM=False, color=None, line_style=1):
        if len(values) < self.n_bins:
            raise ValueError(
                'The passed spectrum \'{0}\' has {1} bins, but the number of passed bins is {2}'.format(
                    name, len(values), self.n_bins) +
                '\n    values: {0}'.format(values) +
                '\n    bin_boundaries: {0}'.format(self.bin_boundaries)
                )
        elif len(values) > self.n_bins:
            logging.warning(
                'The passed spectrum \'{0}\' has {1} bins, but the number of passed bins is {2}'.format(
                    name, len(values), self.n_bins) +
                '\n    Will keep only the first {0} bins'.format(self.n_bins)
                )
            values = values[:self.n_bins]

        if color is None:
            color = next(self.color_cycle)
        self.spectra.append(Spectrum(name, values, isSM, color=color, line_style=line_style))

    def set_bin_boundaries(self, bin_boundaries, add_overflow=False):
        self.bin_boundaries = bin_boundaries
        if add_overflow:
            self.n_bins = len(self.bin_boundaries)
            self.last_bin_is_overflow = True
            self.bin_boundaries.append(10000.)
        else:
            self.n_bins = len(self.bin_boundaries)-1
    
    def get_SM(self):
        return [s for s in self.spectra if s.isSM][0]

    def get_extrema(self):
        SM = self.get_SM()
        x_min = self.bin_boundaries[0]
        if self.last_bin_is_overflow:
            x_max = self.bin_boundaries[-2] + (self.bin_boundaries[-2]-self.bin_boundaries[-3])
        else:
            x_max = self.bin_boundaries[-1]

        if self.logscale:
            y_minAbs = min([ y for y in SM.values if y >= 0.0 ])
            y_maxAbs = max(SM.values)
            y_min = 0.01 * y_minAbs
            y_max = 2.0 * y_maxAbs
        else:
            y_minAbs = min(SM.values)
            y_maxAbs = max(SM.values)
            y_min = y_minAbs - 0.1*(y_maxAbs-y_minAbs)
            y_max = y_maxAbs + 0.1*(y_maxAbs-y_minAbs)

        return x_min, y_min, x_max, y_max


    def plot_spectra(self, tag=None):
        outname = 'normalizationCrosscheck_reimpl'
        if not(tag is None):
            outname += '_' + tag
        plot = differentials.plotting.plots.BottomPanelPlot(outname)

        x_min, y_min, x_max, y_max = self.get_extrema()
        plot.top_x_min = x_min
        plot.top_y_min = y_min
        plot.top_x_max = x_max
        plot.top_y_max = y_max
        plot.bottom_x_max = x_max
        plot.bottom_y_max = y_max

        for spectrum in self.spectra:
            H = self.spectrum_to_hist(spectrum)
            # H.Draw('HISTSAME')
            # leg.AddEntry( H.GetName(), spectrum.title, 'l' )
            plot.add_top(H, 'HISTSAME')

        plot.draw()
        plot.wrapup()



    def plot_spectra_old(self, tag=None):
        SM = self.get_SM()

        x_min = self.bin_boundaries[0]
        if self.last_bin_is_overflow:
            x_max = self.bin_boundaries[-2] + (self.bin_boundaries[-2]-self.bin_boundaries[-3])
        else:
            x_max = self.bin_boundaries[-1]

        if self.logscale:
            y_minAbs = min([ y for y in SM.values if y >= 0.0 ])
            y_maxAbs = max(SM.values)
            y_min = 0.01 * y_minAbs
            y_max = 2.0 * y_maxAbs
        else:
            y_minAbs = min(SM.values)
            y_maxAbs = max(SM.values)
            y_min = y_minAbs - 0.1*(y_maxAbs-y_minAbs)
            y_max = y_maxAbs + 0.1*(y_maxAbs-y_minAbs)

        c.Clear()
        c.set_margins()

        base = differentials.plotting.pywrappers.Base(
            x_min = x_min,
            x_max = x_max,
            y_min = y_min,
            y_max = y_max,
            x_title = 'p_{T} (GeV)',
            y_title = '#sigma (pb/GeV)',
            )
        base.Draw()

        leg = ROOT.TLegend(
            1 - c.GetRightMargin() - 0.3,
            1 - c.GetTopMargin()   - 0.25,
            1 - c.GetRightMargin(),
            1 - c.GetTopMargin() 
            )
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)

        for spectrum in self.spectra:
            H = self.spectrum_to_hist(spectrum)
            H.Draw('HISTSAME')
            leg.AddEntry( H.GetName(), spectrum.title, 'l' )

        leg.Draw()

        if self.logscale: c.SetLogy()

        outname = 'normalizationCrosscheck_reimpl'
        if not(tag is None):
            outname += '_' + tag
        c.save(outname)


    def spectrum_to_hist(self, s):
        H = ROOT.TH1F(
            s.name+'_'+ differentials.plotting.plotting_utils.get_unique_rootname(), s.title,
            self.n_bins, array('f', self.bin_boundaries)
            )
        ROOT.SetOwnership( H, False )
        for i_bin in xrange(self.n_bins):
            H.SetBinContent(i_bin+1, s.values[i_bin])
        H.SetLineWidth(2)
        H.SetLineColor(s.color)
        H.SetLineStyle(s.line_style)
        return H



@flag_as_option
def normalization_xcheck_yukawa(args):
    normalizationCrossCheck = NormalizationCrossCheck()

    from LatestBinning import obs_pth_ggH as obs_pth
    obs_pth.drop_bins_up_to_value(120.)

    normalizationCrossCheck.set_bin_boundaries(obs_pth.binning, add_overflow=False)
    bin_boundaries = normalizationCrossCheck.bin_boundaries # Convenience

    # Add SM
    SMXSs = obs_pth.crosssection_over_binwidth()
    normalizationCrossCheck.add_spectrum('YR4, shape Vitt.', SMXSs, isSM=True, color=4)

    # Add raw from Pier, only rebinned
    SM_Pier = differentials.theory.theory_utils.FileFinder(
        kappab=1.0, kappac=1.0, muR=1.0, muF=1.0, Q=1.0,
        directory=LatestPaths.theory.yukawa.filedir
        ).get_one()
    SM_Pier_rebinned = differentials.theory.theory_utils.rebin_theory(SM_Pier, [0., 15., 30., 45., 80., 120.])
    normalizationCrossCheck.add_spectrum('Pier', SM_Pier_rebinned.crosssection, color=2)

    ws = LatestPaths.ws.yukawa.nominal.combination
    with differentials.core.openroot(ws) as ws_fp:
        w = ws_fp.Get('w')
        reweightors = []
        for left, right in zip(obs_pth.binning[:-1], obs_pth.binning[1:]):
            name = 'reweightor_ggH_PTH_{0:d}_{1:d}'.format(int(left), int(right))
            reweightor = w.function(name)
            if reweightor == None:
                raise ValueError('Workspace {0} does not have a function called {1}'.format(ws, name))
            reweightors.append(reweightor)
        reweighting_factors = [ rew.getVal() for rew in reweightors ]
        # ws_smxs = [ xspar.getVal() for xspar in differentials.core.read_set(w, 'SMXS', return_names=False) ]
        # widths = [ r-l for l,r in zip(bin_boundaries[:-1], bin_boundaries[1:]) ]
        # ws_xs = [ rew*xs/width for rew, xs, width in zip(reweighting_factors, ws_smxs, widths) ]
        # normalizationCrossCheck.add_spectrum('Pier reweighted', ws_xs, line_style=2, color=2)
        Vitt_reweighted = [ rew*xs for rew, xs in zip(reweighting_factors, SMXSs)]
        normalizationCrossCheck.add_spectrum('Vitt. reweighted', Vitt_reweighted, line_style=2, color=4)

    normalizationCrossCheck.plot_spectra('Yukawa')

