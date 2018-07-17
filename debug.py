#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

from OptionHandler import flag_as_option

import LatestPaths
import LatestBinning

import differentials
import differentialutils

import os.path
import logging
import copy
import random
import glob
random.seed(1002)

from time import strftime
datestr = strftime('%b%d')

########################################
# Main
########################################

@flag_as_option
def debug_print_hbb_pdfs(args):

    w = differentials.core.get_ws('projections/workspaces_Jun28/ws_pth_smH_hbb.root')

    pars = {
        'qcdeff' : 0.0130548,
        'r0p1' : 1.60581,
        'r1p0' : 1.0848,
        'r1p1' : 1.23162,
        'r2p0' : 1.45694,
        'r2p1' : 0.915351,
        'r3p0' : 0.98153,
        'r3p1' : 0.0387635,
        }

    def set_w(w):
        print '\n' + '-'*80 + '\nSetting some values\n\n'
        for key, value in pars.iteritems():
            w.var(key).setVal(value)
            print 'Setting {0} to {1}'.format(key, value)

    def print_w(w):
        print '\n' + '-'*80 + '\nPrinting some values\n\n'
        w.pdf('shapeBkg_qcd_cat2_pass_cat2').Print()

    print_w(w)
    set_w(w)
    print_w(w)




@flag_as_option
def testing_spline2d(args):
    x_coupling = 'ct'
    y_coupling = 'cg'

    x_min = 0.0
    x_max = 2.5
    y_min = -0.065
    y_max = 0.08

    x_min = -0.1
    x_max = 2.0
    y_min = -0.05
    y_max = 0.085


    args = differentialutils.set_one_decay_channel(args, 'combWithHbb')
    if args.asimov:
        scandir = 'out/Scan_May22_Top_combWithHbb_scalingttH_couplingdependentBRs_asimov'
    else:
        scandir = 'out/Scan_May31_Top_combWithHbb_scalingttH_couplingdependentBRs'

    scalingttH = differentials.scans.Scan2D(
        'scalingttH', x_coupling, y_coupling,
        scandir = scandir
        )
    scalingttH.title = 'Combination (incl. bbH / BR(#vec{#kappa}))'
    scalingttH.color = 1
    scalingttH.read()
    unsplined_hist = scalingttH.to_hist()

    spline = scalingttH.to_spline(
        x_min = x_min,
        x_max = x_max,
        y_min = y_min,
        y_max = y_max,
        )
    spline.add_noise_selector(
        lambda ct, cg: (cg  <  (1./12.)-0.02 - (1./12.)*ct)
        )
    spline.add_noise_selector(
        lambda ct, cg: (cg  >  (1./12.)+0.04 - (1./13.)*ct)
        )

    splined_hist = spline.to_hist()
    splined_hist.color = 1


    plot = differentials.plotting.plots.Single2DHistPlot(
        'debug_unsplined',
        unsplined_hist,
        x_min = x_min,
        x_max = x_max,
        y_min = y_min,
        y_max = y_max,
        )
    plot.set_ranges_by_contour = False
    plot.x_SM = 1.0
    plot.y_SM = 0.0
    plot.x_title = '#kappa_{t}'
    plot.y_title = 'c_{g}'
    plot.draw()
    plot.wrapup()

    plot = differentials.plotting.plots.Single2DHistPlot(
        'debug_splined',
        splined_hist,
        x_min = x_min,
        x_max = x_max,
        y_min = y_min,
        y_max = y_max,
        )
    plot.set_ranges_by_contour = False
    plot.x_SM = 1.0
    plot.y_SM = 0.0
    plot.x_title = '#kappa_{t}'
    plot.y_title = 'c_{g}'
    plot.draw()
    plot.wrapup()


    #____________________________________________________________________
    # Get 1d

    ct_graph = get_1d_x(unsplined_hist)
    get_unc_for_graph(ct_graph)
    ct_graph.draw_style = 'repr_smooth_line'
    ct_graph.title = '{0} observed; ({1:.2f} - {2:.2f}) @ 68% CL'.format(
        differentials.core.standard_titles['ct'],
        ct_graph.unc.left_bound, ct_graph.unc.right_bound
        )
    ct_graph.color = 1

    differentials.plotting.canvas.c.resize_temporarily(850, 800)
    plot = differentials.plotting.plots.MultiScanPlot('onekappascan_ct')
    plot.manual_graphs.append(ct_graph)
    # plot.manual_graphs.append(exp1D)
    plot.x_title = differentials.core.standard_titles['ct']
    plot.x_min = x_min
    plot.x_max = x_max
    plot.leg.SetNColumns(1)
    plot.leg._y1 = lambda c: 1. - c.GetTopMargin() - 0.20
    plot.draw()
    plot.wrapup()


    cg_graph = get_1d_y(unsplined_hist)
    get_unc_for_graph(cg_graph)
    cg_graph.draw_style = 'repr_smooth_line'
    cg_graph.title = '{0} observed; ({1:.2f} - {2:.2f}) @ 68% CL'.format(
        differentials.core.standard_titles['cg'],
        cg_graph.unc.left_bound, cg_graph.unc.right_bound
        )
    cg_graph.color = 1

    differentials.plotting.canvas.c.resize_temporarily(850, 800)
    plot = differentials.plotting.plots.MultiScanPlot('onekappascan_cg')
    plot.manual_graphs.append(cg_graph)
    # plot.manual_graphs.append(exp1D)
    plot.x_title = differentials.core.standard_titles['cg']
    plot.x_min = y_min
    plot.x_max = y_max
    plot.leg.SetNColumns(1)
    plot.leg._y1 = lambda c: 1. - c.GetTopMargin() - 0.20
    plot.draw()
    plot.wrapup()





def get_1d_x(histogram2D):
    # Histogram should be filled with 2*deltaNLL!!
    xs = []
    deltaNLLs = []

    x_bestfit = histogram2D.bestfit().x
    plugin_bestfit = True

    for i_center, center in enumerate(histogram2D.x_bin_centers):
        deltaNLL = min(histogram2D.H2_array[i_center][:])

        if plugin_bestfit and center > x_bestfit:
            xs.append(x_bestfit)
            deltaNLLs.append(0.0)
            plugin_bestfit = False

        xs.append(center)
        deltaNLLs.append(deltaNLL)

    graph = differentials.plotting.pywrappers.Graph(
        differentials.plotting.plotting_utils.get_unique_rootname(),
        'title',
        xs,
        deltaNLLs
        )
    return graph

def get_1d_y(histogram2D):
    # Histogram should be filled with 2*deltaNLL!!
    xs = []
    deltaNLLs = []

    x_bestfit = histogram2D.bestfit().y
    plugin_bestfit = True

    for i_center, center in enumerate(histogram2D.y_bin_centers):
        deltaNLL = min([ row[i_center] for row in histogram2D.H2_array ])

        if plugin_bestfit and center > x_bestfit:
            xs.append(x_bestfit)
            deltaNLLs.append(0.0)
            plugin_bestfit = False

        xs.append(center)
        deltaNLLs.append(deltaNLL)

    graph = differentials.plotting.pywrappers.Graph(
        differentials.plotting.plotting_utils.get_unique_rootname(),
        'title',
        xs,
        deltaNLLs
        )
    return graph




uncertaintycalculator = differentials.uncertaintycalculator.UncertaintyCalculator()
def get_unc_for_graph(graph):
    graph.unc = uncertaintycalculator.create_uncertainties(
        graph.xs,
        [ 0.5*y for y in graph.ys ] # Histogram has 2dNLL, but unc calculator expects dNLL
        )


@flag_as_option
def debug_test_ctcg_files(args):
    differentials.theory.ctcg_interpreter.create_all_ctcg()


@flag_as_option
def debug_scan_hgg_Mar15_pth_ggH_GT600_bestfit(args):
    cmd = [
        'combine',
        'out/workspaces_Mar02/ws_pth_ggH_hgg.root',
        '-n _debug_hgg_GT600_bestfit',
        '-M MultiDimFit',
        '-m 125.0',
        '--cminDefaultMinimizerType Minuit2',
        '--cminDefaultMinimizerAlgo migrad',
        '--saveNLL',
        '--saveInactivePOI 1',
        '--floatOtherPOIs=1',
        '-v 1',
        '-P r_ggH_PTH_GT600',
        '--setPhysicsModelParameterRanges r_ggH_PTH_0_15=-1.0,4.0:r_ggH_PTH_15_30=-1.0,4.0:r_ggH_PTH_30_45=-1.0,4.0:r_ggH_PTH_45_80=-1.0,4.0:r_ggH_PTH_80_120=-1.0,4.0:r_ggH_PTH_120_200=-1.0,4.0:r_ggH_PTH_200_350=-1.0,4.0:r_ggH_PTH_350_600=-15.0,10.0:r_ggH_PTH_GT600=-1000.0,10.0',
        '--setPhysicsModelParameters r_ggH_PTH_0_15=1.0,r_ggH_PTH_15_30=1.0,r_ggH_PTH_30_45=1.0,r_ggH_PTH_45_80=1.0,r_ggH_PTH_80_120=1.0,r_ggH_PTH_120_200=1.0,r_ggH_PTH_200_350=1.0,r_ggH_PTH_350_600=1.0,r_ggH_PTH_GT600=1.0',
        ]
    differentials.core.execute(cmd)


pdf_indices = [
    'pdfindex_recoPt_350p0_600p0_SigmaMpTTag_0_13TeV',
    'pdfindex_recoPt_350p0_600p0_SigmaMpTTag_1_13TeV',
    'pdfindex_recoPt_350p0_600p0_SigmaMpTTag_2_13TeV',
    'pdfindex_recoPt_600p0_10000p0_SigmaMpTTag_0_13TeV',
    'pdfindex_recoPt_600p0_10000p0_SigmaMpTTag_1_13TeV',
    'pdfindex_recoPt_600p0_10000p0_SigmaMpTTag_2_13TeV',
    'pdfindex_recoPt_45p0_80p0_SigmaMpTTag_1_13TeV',
    'pdfindex_recoPt_45p0_80p0_SigmaMpTTag_0_13TeV',
    'pdfindex_recoPt_45p0_80p0_SigmaMpTTag_2_13TeV',
    'pdfindex_recoPt_80p0_120p0_SigmaMpTTag_2_13TeV',
    'pdfindex_recoPt_80p0_120p0_SigmaMpTTag_1_13TeV',
    'pdfindex_recoPt_80p0_120p0_SigmaMpTTag_0_13TeV',
    'pdfindex_recoPt_200p0_350p0_SigmaMpTTag_1_13TeV',
    'pdfindex_recoPt_200p0_350p0_SigmaMpTTag_0_13TeV',
    'pdfindex_recoPt_200p0_350p0_SigmaMpTTag_2_13TeV',
    'pdfindex_recoPt_0p0_15p0_SigmaMpTTag_2_13TeV',
    'pdfindex_recoPt_0p0_15p0_SigmaMpTTag_0_13TeV',
    'pdfindex_recoPt_0p0_15p0_SigmaMpTTag_1_13TeV',
    'pdfindex_recoPt_15p0_30p0_SigmaMpTTag_0_13TeV',
    'pdfindex_recoPt_15p0_30p0_SigmaMpTTag_1_13TeV',
    'pdfindex_recoPt_15p0_30p0_SigmaMpTTag_2_13TeV',
    'pdfindex_recoPt_120p0_200p0_SigmaMpTTag_2_13TeV',
    'pdfindex_recoPt_120p0_200p0_SigmaMpTTag_1_13TeV',
    'pdfindex_recoPt_120p0_200p0_SigmaMpTTag_0_13TeV',
    'pdfindex_recoPt_30p0_45p0_SigmaMpTTag_1_13TeV',
    'pdfindex_recoPt_30p0_45p0_SigmaMpTTag_0_13TeV',
    'pdfindex_recoPt_30p0_45p0_SigmaMpTTag_2_13TeV',
    ]

@flag_as_option
def debug_scan_hgg_Mar15_pth_smH(args):
    cmd = [
        'combine',
        'out/workspaces_Mar14/ws_pth_smH_hgg.root',
        '-n _debug_hgg_Mar15_pth_smH',
        '-t -1',
        '--cminDefaultMinimizerType Minuit2',
        '--cminDefaultMinimizerAlgo migrad',
        '--saveNLL',
        '--saveInactivePOI 1',
        '--floatOtherPOIs=1',
        '-v 1',
        '-P r_smH_PTH_GT600',
        '--setPhysicsModelParameters r_smH_PTH_0_15=1.0,r_smH_PTH_15_30=1.0,r_smH_PTH_30_45=1.0,r_smH_PTH_45_80=1.0,r_smH_PTH_80_120=1.0,r_smH_PTH_120_200=1.0,r_smH_PTH_200_350=1.0,r_smH_PTH_350_600=1.0,r_smH_PTH_GT600=1.0',
        '--algo=grid',
        '-M MultiDimFit',
        '-m 125.0',
        '--setPhysicsModelParameterRanges r_smH_PTH_0_15=-1.0,4.0:r_smH_PTH_15_30=-1.0,4.0:r_smH_PTH_30_45=-1.0,4.0:r_smH_PTH_45_80=-1.0,4.0:r_smH_PTH_80_120=-1.0,4.0:r_smH_PTH_120_200=-1.0,4.0:r_smH_PTH_200_350=-1.0,4.0:r_smH_PTH_350_600=-15.0,10.0:r_smH_PTH_GT600=-10.0,20.0',
        '--points 40',
        '--firstPoint 5',
        '--lastPoint 7',
        '--saveSpecifiedIndex {0}'.format(','.join(pdf_indices))
        ]
    differentials.core.execute(cmd)

@flag_as_option
def debug_scan_hgg_Mar15_pth_ggH(args):
    cmd = [
        'combine',
        'out/workspaces_Mar02/ws_pth_ggH_hgg.root',
        '-n _debug_hgg_Mar15_pth_ggH',
        '-t -1',
        '--cminDefaultMinimizerType Minuit2',
        '--cminDefaultMinimizerAlgo migrad',
        '--saveNLL',
        '--saveInactivePOI 1',
        '--floatOtherPOIs=1',
        '-v 1',
        '-P r_ggH_PTH_GT600',
        '--setPhysicsModelParameters r_ggH_PTH_0_15=1.0,r_ggH_PTH_15_30=1.0,r_ggH_PTH_30_45=1.0,r_ggH_PTH_45_80=1.0,r_ggH_PTH_80_120=1.0,r_ggH_PTH_120_200=1.0,r_ggH_PTH_200_350=1.0,r_ggH_PTH_350_600=1.0,r_ggH_PTH_GT600=1.0',
        '--algo=grid',
        '-M MultiDimFit',
        '-m 125.0',
        '--setPhysicsModelParameterRanges r_ggH_PTH_0_15=-1.0,4.0:r_ggH_PTH_15_30=-1.0,4.0:r_ggH_PTH_30_45=-1.0,4.0:r_ggH_PTH_45_80=-1.0,4.0:r_ggH_PTH_80_120=-1.0,4.0:r_ggH_PTH_120_200=-1.0,4.0:r_ggH_PTH_200_350=-1.0,4.0:r_ggH_PTH_350_600=-15.0,10.0:r_ggH_PTH_GT600=-10.0,20.0',
        '--points 40',
        '--firstPoint 5',
        '--lastPoint 7',
        '--saveSpecifiedIndex {0}'.format(','.join(pdf_indices))
        ]
    differentials.core.execute(cmd)


@flag_as_option
def debug_scan_hzz_Mar15(args):
    cmd = [
        'combine',
        'out/workspaces_Feb28/ws_pth_ggH_hzz.root',
        '-n _debug_scan_hzz_Mar15',
        '-M MultiDimFit',
        '-m 125.0',
        '--cminDefaultMinimizerType Minuit2',
        '--cminDefaultMinimizerAlgo migrad',
        '--saveNLL',
        '--saveInactivePOI 1',
        '--floatOtherPOIs=1',
        # '-v -1',
        '-P r_ggH_PTH_GT200',
        '--setPhysicsModelParameterRanges r_ggH_PTH_0_15=0.0,4.0:r_ggH_PTH_15_30=0.0,4.0:r_ggH_PTH_30_80=0.0,4.0:r_ggH_PTH_80_200=0.0,4.0:r_ggH_PTH_GT200=-10.0,1.0',
        '--setPhysicsModelParameters r_ggH_PTH_0_15=1.0,r_ggH_PTH_15_30=1.0,r_ggH_PTH_30_80=1.0,r_ggH_PTH_80_200=1.0,r_ggH_PTH_GT200=1.0',
        '--algo=grid',
        '--points=40',
        ]
    differentials.core.execute(cmd)


@flag_as_option
def debug_draw_fastscans_Mar12(args):
    scandirs = [
        'out/Scan_Mar06_Top_combWithHbb_asimov',
        'out/Scan_Mar07_Top_combWithHbb_asimov',
        'out/Scan_Mar12_Top_combWithHbb_last2BinsDropped_asimov',
        'out/Scan_Mar06_Top_combination_asimov',
        'out/Scan_Mar07_Top_combination_asimov_5',
        'out/Scan_Mar12_Top_combination_last2BinsDropped_asimov',
        ]
    for scandir in scandirs:
        fastscan_file = glob.glob(scandir + '/postfit_and_fastscan/*FASTSCAN*.root')[0]
        name = os.path.basename(fastscan_file).replace('.root','')
        fastscan = differentials.scans.Scan2D(name, 'ct', 'cg')
        fastscan.root_files = [fastscan_file]
        fastscan.read()
        fastscan.plot(name, draw_style='repr_2D_rainbow_high_contours')

@flag_as_option
def debug_test_setreading(args):
    ws = LatestPaths.ws.yukawa.nominal.combination
    differentials.core.read_set(ws, 'SMXS')
    differentials.core.read_set(ws, 'reweightors')


@flag_as_option
def debug_test_integrator(args):
    bin_boundaries = [ 0., 5., 10., 15., 20., 25., 30. ]
    bin_values = [1.0 for i in xrange(len(bin_boundaries)-1)]
    integral = differentials.integral.Integrator(bin_boundaries, bin_values)
    integral.integral(7., 23.)


@flag_as_option
def debug_test_kappabkappac_files(args):
    interp = differentials.theory.kappabkappac_interpreter.KappabKappacInterpreter()
    # interp.dump_gluon_induced()
    interp.dump_quark_induced_scaled()
    # interp.dump_summed_quark_gluon_induced()


@flag_as_option
def debug_test_job_registering(args):

    output = """
Your job 8086673 ("job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_68_0.sh") has been submitted
Created job script: job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_68_0.sh
>> qsub -q all.q /mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/v4_NewBinning/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/out/Scan_Mar06_Top_hgg_asimov/job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_68_0.sh
Your job 8086674 ("job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_69_0.sh") has been submitted
Created job script: job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_69_0.sh
>> qsub -q all.q /mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/v4_NewBinning/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/out/Scan_Mar06_Top_hgg_asimov/job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_69_0.sh
Your job 8086675 ("job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_70_0.sh") has been submitted
Created job script: job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_70_0.sh
>> qsub -q all.q /mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/v4_NewBinning/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/out/Scan_Mar06_Top_hgg_asimov/job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_70_0.sh
Your job 8086677 ("job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_71_0.sh") has been submitted
Created job script: job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_71_0.sh
>> qsub -q all.q /mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/v4_NewBinning/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/out/Scan_Mar06_Top_hgg_asimov/job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_71_0.sh
Your job 8086685 ("job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_72_0.sh") has been submitted
Created job script: job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_72_0.sh
>> qsub -q all.q /mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/v4_NewBinning/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/out/Scan_Mar06_Top_hgg_asimov/job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_72_0.sh
Your job 8086691 ("job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_73_0.sh") has been submitted
Created job script: job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_73_0.sh
>> qsub -q all.q /mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/v4_NewBinning/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/out/Scan_Mar06_Top_hgg_asimov/job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_73_0.sh
Your job 8086692 ("job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_74_0.sh") has been submitted
Created job script: job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_74_0.sh
>> qsub -q all.q /mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/v4_NewBinning/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/out/Scan_Mar06_Top_hgg_asimov/job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_74_0.sh"""

    testconfig = CombineToolWrapper.CombineConfig(args)
    testscan = CombineToolWrapper.BaseCombineScan(testconfig)

    testscan.register_jobids_in_jobmanager(output)
    


@flag_as_option
def debug_draw_fastscans(args):
    
    all_fastscans = glob.glob('out/Scan*Mar12*Top*/postfit_and_fastscan/*FASTSCAN*.root')
    for fastscan_file in all_fastscans:
        name = os.path.basename(fastscan_file).replace('.root','')
        fastscan = differentials.scans.Scan2D(name, 'ct', 'cg')
        fastscan.root_files = [fastscan_file]
        fastscan.read()
        fastscan.plot(name, draw_style='repr_2D_rainbow_high_contours')


@flag_as_option
def debug_Scan_hbb(args):
    cmd = [
        'combine',
        # 'out/workspaces_Feb28/ws_pth_ggH_hbb.root',
        'out/postfits_Mar02/higgsCombine_POSTFIT_ws_pth_ggH_hbb.MultiDimFit.mH125.root',
        # 
        '-n _DEBUG_{0}_SCAN_bPOI_r_ggH_PTH_350_600_ePOI_ws_pth_ggH_hbb'.format(datestr),
        '-M MultiDimFit',
        '-m 125.0',
        # '--cminDefaultMinimizerType Minuit2',
        # '--cminDefaultMinimizerAlgo migrad',
        '--minimizerStrategy 2',
        '--minimizerTolerance 0.001',
        '--robustFit 1',
        '--minimizerAlgoForMinos Minuit2,Migrad',
        '--saveNLL',
        '--saveInactivePOI 1',
        '--floatOtherPOIs=1',
        '-P r_ggH_PTH_350_600',
        '--setPhysicsModelParameterRanges r_ggH_PTH_200_350=-10.0,10.0:r_ggH_PTH_350_600=-10.0,10.0:r_ggH_PTH_GT600=-10.0,10.0',
        # :qcdeff=0.001,8.0:r1p0=0.0,8.0:r2p0=0.0,8.0:r3p0=0.0,8.0:r0p1=0.0,8.0:r1p1=0.0,8.0:r2p1=0.0,8.0:r3p1=0.0,8.0',
        # '--setPhysicsModelParameters r_ggH_PTH_200_350=1.0,r_ggH_PTH_350_600=1.0,r_ggH_PTH_GT600=1.0',
        '--algo=grid',
        '--points=3',
        # 
        '--snapshotName MultiDimFit',
        '--skipInitialFit',
        '-v 2',
        ]
    differentials.core.execute(cmd)


@flag_as_option
def debug_combWithHbb_Scan_Mar02(args):
    # ws = 'out/workspaces_Mar01/ws_pth_ggH_combWithHbb.root'

    # ws = 'out/postfits_Mar02_v1/higgsCombine_POSTFIT_ws_pth_ggH_combination.MultiDimFit.mH125.root'
    # ws = 'out/postfits_Mar02_v1/higgsCombine_POSTFIT_ASIMOV_ws_pth_ggH_combination.MultiDimFit.mH125.root'
    # ws = 'out/postfits_Mar02_v1/higgsCombine_POSTFIT_ws_pth_ggH_combWithHbb.MultiDimFit.mH125.root'
    # ws = 'out/postfits_Mar02_v1/higgsCombine_POSTFIT_ASIMOV_ws_pth_ggH_combWithHbb.MultiDimFit.mH125.root'

    ws = 'out/postfits_Mar02/higgsCombine_POSTFIT_ASIMOV_ws_pth_ggH_combWithHbb.MultiDimFit.mH125.root'
    # ws = 'out/postfits_Mar02/higgsCombine_POSTFIT_ASIMOV_ws_pth_ggH_hbb.MultiDimFit.mH125.root'
    # ws = 'out/postfits_Mar02/higgsCombine_POSTFIT_ws_pth_ggH_combWithHbb.MultiDimFit.mH125.root'
    # ws = 'out/postfits_Mar02/higgsCombine_POSTFIT_ws_pth_ggH_hbb.MultiDimFit.mH125.root'

    cmd = [
        'combine',
        ws,
        '-n _DEBUG_{0}_SCAN'.format(datestr),
        '-M MultiDimFit',
        '-m 125.0',
        # 
        '--minimizerStrategy 2',
        '--minimizerTolerance 0.001',
        '--robustFit 1',
        '--minimizerAlgoForMinos Minuit2,Migrad',
        # 
        '--saveNLL',
        '--saveInactivePOI 1',
        '--floatOtherPOIs=1',
        # '--job-mode psi --task-name _SCAN_bPOI_r_ggH_PTH_600_10000_ePOI_ws_pth_ggH_combWithHbb --sub-opts='-q short.q' ',
        # 
        # '-P r_ggH_PTH_GT600',
        '-P r_ggH_PTH_120_200',
        # 
        '--setPhysicsModelParameterRanges r_ggH_PTH_0_15={0},{1}:r_ggH_PTH_15_30={0},{1}:r_ggH_PTH_30_45={0},{1}:r_ggH_PTH_45_80={0},{1}:r_ggH_PTH_80_120={0},{1}:r_ggH_PTH_120_200={0},{1}:r_ggH_PTH_200_350={0},{1}:r_ggH_PTH_350_600={0},{1}:r_ggH_PTH_GT600={0},{1}:qcdeff=0.001,8.0:r1p0=0.0,8.0:r2p0=0.0,8.0:r3p0=0.0,8.0:r0p1=0.0,8.0:r1p1=0.0,8.0:r2p1=0.0,8.0:r3p1=0.0,8.0'.format(
            0.0, 8.0
            ),
        '--setPhysicsModelParameters r_ggH_PTH_0_15=1.0,r_ggH_PTH_15_30=1.0,r_ggH_PTH_30_45=1.0,r_ggH_PTH_45_80=1.0,r_ggH_PTH_80_120=1.0,r_ggH_PTH_120_200=1.0,r_ggH_PTH_200_350=1.0,r_ggH_PTH_350_600=1.0,r_ggH_PTH_GT600=1.0',
        '--algo=grid',
        '--points=3',
        # '--split-points 2',
        '--snapshotName MultiDimFit',
        '--skipInitialFit',
        '-v 2',
        # 
        # '--freezeNuisances qcdeff,r1p0,r2p0,r3p0,r0p1,r1p1,r2p1,r3p1',
        # 
        '--cminDefaultMinimizerType Minuit2',
        '--cminDefaultMinimizerAlgo migrad',
        # 
        '-t -1'
        ]
    differentials.core.execute(cmd)

    # qcdeff r1p0 r2p0 r3p0 r0p1 r1p1 r2p1 r3p1

    # combine comb_2017_ggHbb.root
    # --cminDefaultMinimizerType Minuit2
    # --cminDefaultMinimizerAlgo migrad
    # --algo=grid
    # --floatOtherPOIs=1
    # -P r_ggH_PTH_350_600
    # --setPhysicsModelParameters r_ggH_PTH_GT600=1.0,r_ggH_PTH_350_600=1.0
    # --squareDistPoi
    # --saveNLL
    # --saveInactivePOI 1
    # -t -1
    # --minimizerStrategy 2
    # --minimizerTolerance 0.001
    # --robustFit 1
    # --minimizerAlgoForMinos Minuit2,Migrad
    # -M MultiDimFit
    # -m 125.00
    # --setPhysicsModelParameterRanges r_ggH_PTH_350_600=-10.000,10.000


@flag_as_option
def debug_hbb_bestfit(args):
    cmd = [
        'combine',
        'out/workspaces_Feb28/ws_pth_ggH_hbb.root',
        '-n _POSTFIT_ASIMOV_ws_pth_ggH_hbb',
        '-M MultiDimFit',
        '-m 125.0',
        '-t -1',
        # 
        '--cminDefaultMinimizerType Minuit2',
        '--cminDefaultMinimizerAlgo migrad',
        # 
        '--minimizerStrategy 2',
        '--minimizerTolerance 0.001',
        '--robustFit 1',
        '--minimizerAlgoForMinos Minuit2,Migrad',
        '--saveNLL',
        '--saveInactivePOI 1',
        '--floatOtherPOIs=1',
        '-P r_ggH_PTH_200_350 -P r_ggH_PTH_350_600 -P r_ggH_PTH_GT600',
        '--setPhysicsModelParameterRanges r_ggH_PTH_200_350=-10.0,10.0:r_ggH_PTH_350_600=-10.0,10.0:r_ggH_PTH_GT600=-10.0,10.0:qcdeff=0.001,8.0:r1p0=0.0,8.0:r2p0=0.0,8.0:r3p0=0.0,8.0:r0p1=0.0,8.0:r1p1=0.0,8.0:r2p1=0.0,8.0:r3p1=0.0,8.0',
        '--setPhysicsModelParameters r_ggH_PTH_200_350=1.0,r_ggH_PTH_350_600=1.0,r_ggH_PTH_GT600=1.0',
        '--saveWorkspace',
        ]
    differentials.core.execute(cmd)


@flag_as_option
def debug_combWithHbb_Mar01(args):
    ws = 'out/workspaces_Mar01/ws_pth_ggH_combWithHbb.root'
    cmd = [
        'combine',
        ws,
        '-n _DEBUG_SCAN_bPOI_r_ggH_PTH_600_10000_ePOI_ws_pth_ggH_combWithHbb',
        '-M MultiDimFit',
        '-m 125.0',
        '--minimizerStrategy 2',
        '--minimizerTolerance 0.001',
        '--robustFit 1',
        '--minimizerAlgoForMinos Minuit2,Migrad',
        '--saveNLL',
        '--saveInactivePOI 1',
        '--floatOtherPOIs=1',
        # '--job-mode psi --task-name _SCAN_bPOI_r_ggH_PTH_600_10000_ePOI_ws_pth_ggH_combWithHbb --sub-opts='-q short.q' ',
        '-P r_ggH_PTH_GT600',
        '--setPhysicsModelParameterRanges r_ggH_PTH_0_15=0.0,8.5:r_ggH_PTH_15_30=0.0,8.5:r_ggH_PTH_30_45=0.0,8.5:r_ggH_PTH_45_80=0.0,8.5:r_ggH_PTH_80_120=0.0,8.5:r_ggH_PTH_120_200=0.0,8.5:r_ggH_PTH_200_350=0.0,8.5:r_ggH_PTH_350_600=0.0,8.5:r_ggH_PTH_GT600=0.0,8.5:r_ggH_PTH_600_10000=0.0,8.5:qcdeff=0.001,8.0:r1p0=0.0,8.0:r2p0=0.0,8.0:r3p0=0.0,8.0:r0p1=0.0,8.0:r1p1=0.0,8.0:r2p1=0.0,8.0:r3p1=0.0,8.0',
        '--setPhysicsModelParameters r_ggH_PTH_0_15=1.0,r_ggH_PTH_15_30=1.0,r_ggH_PTH_30_45=1.0,r_ggH_PTH_45_80=1.0,r_ggH_PTH_80_120=1.0,r_ggH_PTH_120_200=1.0,r_ggH_PTH_200_350=1.0,r_ggH_PTH_350_600=1.0,r_ggH_PTH_GT600=1.0,r_ggH_PTH_600_10000=1.0',
        '--algo=grid',
        '--points=10',
        # '--split-points 2',
        ]
    differentials.core.execute(cmd)

@flag_as_option
def debug_combination_Mar01(args):
    # ws = 'out/workspaces_Mar01/ws_pth_ggH_combination.root'
    cmd = [
        'combine',
        'out/workspaces_Mar01/ws_pth_ggH_combination.root',
        '-n _DEBUG_SCAN_bPOI_r_ggH_PTH_GT600_ePOI_ws_pth_ggH_combination',
        '-M MultiDimFit',
        '-m 125.0',
        '--cminDefaultMinimizerType Minuit2',
        '--cminDefaultMinimizerAlgo migrad',
        '--saveNLL',
        '--saveInactivePOI 1',
        '--floatOtherPOIs=1',
        # '--job-mode psi --task-name _SCAN_bPOI_r_ggH_PTH_GT600_ePOI_ws_pth_ggH_combination --sub-opts='-q short.q' ',
        '-P r_ggH_PTH_GT600',
        '--setPhysicsModelParameterRanges r_ggH_PTH_0_15=0.0,4.0:r_ggH_PTH_15_30=0.0,4.0:r_ggH_PTH_30_45=0.0,4.0:r_ggH_PTH_45_80=0.0,4.0:r_ggH_PTH_80_120=0.0,4.0:r_ggH_PTH_120_200=0.0,4.0:r_ggH_PTH_200_350=0.0,4.0:r_ggH_PTH_350_600=0.0,4.0:r_ggH_PTH_GT600=0.0,4.0',
        '--setPhysicsModelParameters r_ggH_PTH_0_15=1.0,r_ggH_PTH_15_30=1.0,r_ggH_PTH_30_45=1.0,r_ggH_PTH_45_80=1.0,r_ggH_PTH_80_120=1.0,r_ggH_PTH_120_200=1.0,r_ggH_PTH_200_350=1.0,r_ggH_PTH_350_600=1.0,r_ggH_PTH_GT600=1.0',
        '--algo=grid',
        '--points=5',
        # '--split-points 5',
        ]
    differentials.core.execute(cmd)