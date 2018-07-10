import logging
from OptionHandler import flag_as_option

import differentials
import differentials.core as core
import differentials.combine.combine as combine

from vittinterpreter import RTransformer

import glob, re, sys


def base_config(args):
    config = combine.CombineConfig(args)
    config.subDirectory = 'out/corrmats_{0}'.format(core.datestr())
    config.queue = 'short.q'
    config.onBatch = False
    return config

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



