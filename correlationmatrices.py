import logging
from OptionHandler import flag_as_option

import differentials
import differentials.core as core
import differentials.combine.combine as combine

import LatestPaths

from vittinterpreter import RTransformer

import glob, re, sys, os


def base_config(args):
    config = combine.CombineConfig(args)
    config.subDirectory = 'out/corrmats_{0}'.format(core.datestr())
    config.queue = 'short.q'
    config.onBatch = False
    return config


@flag_as_option
def corrmats_cwr(args):
    scandirs = {
        'pth_ggH'  : LatestPaths.scan.pth_ggH.observed.combWithHbb,
        'pth_smH'  : LatestPaths.scan.pth_smH.observed.combWithHbb[0],
        'njets'    : LatestPaths.scan.njets.observed.combination,
        'ptjet'    : LatestPaths.scan.ptjet.observed.combination,
        'rapidity' : LatestPaths.scan.rapidity.observed.combination,
        }

    config = base_config(args)
    for obs in [
            # 'pth_ggH',
            # 'pth_smH',
            # 'njets',
            'ptjet',
            'rapidity'
            ]:
        config.datacard = glob.glob(
            os.path.join(scandirs[obs], 'postfit_and_fastscan/*.root')
            )[0]
        corrmat = combine.CombineCorrMat(config)
        corrmat.run()


@flag_as_option
def plot_corrmats_cwr(args):
    corrmat_files = {
        'njets'    : 'out/corrmats_Sep28/higgsCombine_CORRMAT_higgsCombine_POSTFIT_ws_njets_combination.MultiDimFit.mH125.MultiDimFit.mH125.root',
        'pth_ggH'  : 'out/corrmats_Sep28/higgsCombine_CORRMAT_higgsCombine_POSTFIT_ws_pth_ggH_combWithHbb.MultiDimFit.mH125.MultiDimFit.mH125.root',
        'pth_smH'  : 'out/corrmats_Sep28/higgsCombine_CORRMAT_higgsCombine_POSTFIT_ws_pth_smH_combWithHbb.MultiDimFit.mH125.MultiDimFit.mH125.root',
        'ptjet'    : 'out/corrmats_Sep28/higgsCombine_CORRMAT_higgsCombine_POSTFIT_ws_ptjet_combination.MultiDimFit.mH125.MultiDimFit.mH125.root',
        'rapidity' : 'out/corrmats_Sep28/higgsCombine_CORRMAT_higgsCombine_POSTFIT_ws_rapidity_combination.MultiDimFit.mH125.MultiDimFit.mH125.root',
        }

    for obs in [
            'pth_ggH',
            'pth_smH',
            'njets',
            'ptjet',
            'rapidity'
            ]:
        plot = differentials.plotting.plots_matrix.CorrelationMatrixFromCombinePlot(
            'corrmat_' + obs, corrmat_files[obs]
            )
        plot.get_pois_from_ws()
        plot.pois.sort(key=differentials.core.range_sorter)

        plot.x_title = differentials.core.standard_titles[obs]
        if obs.startswith('pt'): plot.x_title += ' (GeV)'

        n_decimals = 2 if obs == 'rapidity' else 0
        for poi in plot.pois:
            left, right = differentials.core.get_range_from_str(poi)
            if right == 'INF':
                label = '>{0:.{n_decimals}f}'.format(left, n_decimals=n_decimals)
            elif left == '-INF':
                label = '<{0:.{n_decimals}f}'.format(right, n_decimals=n_decimals)
            elif right == 'SINGLE':
                label = '{0:.{n_decimals}f}'.format(left, n_decimals=n_decimals)
            else:
                label = '{0:.{n_decimals}f}-{1:.{n_decimals}f}'.format(left, right, n_decimals=n_decimals)
            plot.binlabels.append(label)

        plot.disable_cms_text = True
        plot.draw()

        l = differentials.plotting.pywrappers.Latex(
            lambda c: 1.0 - 0.01,
            lambda c: 1.0 - c.GetTopMargin(),
            'Correlation'
            )
        l.SetTextFont(42)
        l.SetNDC()
        l.SetTextAngle(-90)
        l.SetTextSize(0.04)
        l.SetTextAlign(13)
        l.Draw()

        differentials.plotting.canvas.c.set_margins(
            LeftMargin  = 0.17,
            TopMargin   = 0.08,
            RightMargin = 0.14,
            BottomMargin = 0.16,
            )
        plot.H.GetXaxis().SetTitleOffset(1.5)
        plot.H.GetYaxis().SetTitleOffset(1.74)
        differentials.plotting.pywrappers.CMS_Latex_type().Draw()
        differentials.plotting.pywrappers.CMS_Latex_lumi().Draw()
        differentials.plotting.canvas.c.change_plotdir_temporarily('corrmats_cwr_{0}'.format(differentials.core.datestr()))
        plot.wrapup()


@flag_as_option
def corrmats_vittorio(args):
    config = base_config(args)
    for bestfit in glob.glob('out/bestfits_Vittorio_Mar12/*.root'):
        config.datacard = bestfit
        corrmat = combine.CombineCorrMat(config)
        corrmat.run()

@flag_as_option
def plot_corrmats_vittorio(args):
    for corrmat_file in glob.glob('out/corrmats_Mar12/*.root'):
        logging.info('Plotting corrmat for {0}'.format(corrmat_file))
        name = 'corrmat_' + re.search(r'differential_(.*)_DataBestFit', corrmat_file).group(1)
        plot = differentials.plotting.plots_matrix.CorrelationMatrixFromCombinePlot(name, corrmat_file)

        rtransformer = RTransformer()
        rtransformer.set_obsname_from_corrmat(corrmat_file)

        # print rtransformer.get_sorted_rs()
        # print rtransformer.get_sorted_binlabels()
        # print rtransformer.process_str_dict
        # sys.exit()

        plot.pois = rtransformer.get_sorted_rs()
        plot.n_pois = len(plot.pois)
        plot.binlabels = rtransformer.get_sorted_binlabels()
        plot.x_title = rtransformer.get_x_title()

        if rtransformer.obsname in ['PtNjets2p5NNLOPS_newBins_v2', 'PtNNLOPS_newBins' ]:
            plot.marker_size = 1.4

        plot.draw()

        if rtransformer.obsname == 'PtNjets2p5NNLOPS_newBins_v2':
            plot.H.GetXaxis().SetLabelSize(0.03)
            plot.H.GetYaxis().SetLabelSize(0.03)

        l = differentials.plotting.pywrappers.Latex(
            lambda c: 1.0 - 0.01,
            lambda c: 1.0 - c.GetTopMargin(),
            'Correlation'
            )
        l.SetTextFont(42)
        l.SetNDC()
        l.SetTextAngle(-90)
        l.SetTextSize(0.04)
        l.SetTextAlign(13)
        l.Draw()

        differentials.plotting.canvas.c.change_plotdir_temporarily('corrmatsvitt_{0}'.format(differentials.core.datestr()))
        plot.wrapup()



