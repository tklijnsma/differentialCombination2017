class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

card = AttrDict()
ws = AttrDict()
theory = AttrDict()
scan = AttrDict()

#____________________________________________________________________
# Text datacards

card.pth_ggH_noXHunc = AttrDict()
card.pth_ggH_noXHunc.hgg = 'suppliedInput/fromVittorio/pT_newBins_Feb28/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_newBins_renamedProcesses.txt'
card.pth_ggH_noXHunc.hzz = 'suppliedInput/fromDavid/PTH_Jan24_newBinning/ggH/hzz4l_comb_13TeV_xs_processesRenumbered.txt'
card.pth_ggH_noXHunc.hbb = 'suppliedInput/fromJavier/bernstein_r7428/comb_2017_ggHbb.txt'
card.pth_ggH_noXHunc.combination = 'suppliedInput/combination_pth_ggH_Mar01.txt'
card.pth_ggH_noXHunc.combWithHbb = 'suppliedInput/combWithHbb_pth_ggH_Mar02.txt'
card.pth_ggH = AttrDict()
card.pth_ggH.hgg = 'suppliedInput/fromVittorio/pT_newBins_Feb28/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_newBins_renamedProcesses_xHNuisPar.txt'
card.pth_ggH.hzz = 'suppliedInput/fromDavid/PTH_Jan24_newBinning/ggH/hzz4l_comb_13TeV_xs_processesRenumbered_xHNuisPar.txt'
card.pth_ggH.hbb = 'suppliedInput/fromJavier/bernstein_r7428/comb_2017_ggHbb_xHNuisPar.txt'
card.pth_ggH.combination = 'suppliedInput/combination_pth_ggH_Mar01_xHNuisPar.txt'
card.pth_ggH.combWithHbb = 'suppliedInput/combWithHbb_pth_ggH_Mar02_xHNuisPar.txt'

card.pth_smH = AttrDict()
card.pth_smH.hgg = card.pth_ggH_noXHunc.hgg # scale ggH/xH by smH!
card.pth_smH.hzz = 'suppliedInput/fromDavid/PTH_Jan24_newBinning/smH/hzz4l_comb_13TeV_xs_processesRenumbered.txt'
card.pth_smH.hbb = 'suppliedInput/fromJavier/bernstein_r7428/comb_2017_ggHbb.txt' # scale ggH/xH by smH!
card.pth_smH.combination = 'suppliedInput/combination_pth_smH_Mar14.txt' # scale ggH/xH by smH!
card.pth_smH.combWithHbb = 'suppliedInput/combWithHbb_pth_smH_Mar14.txt' # scale ggH/xH by smH!

card.njets = AttrDict()
card.njets.hgg = 'suppliedInput/fromVittorio/differential_Njets2p5NNLOPS_Nov10/Datacard_13TeV_differential_Njets2p5NNLOPS_renamedProcesses.txt'
card.njets.hzz = 'suppliedInput/fromDavid/NJ_NNLOPS_Nov01/smH/hzz4l_comb_13TeV_xs.txt'
card.njets.combination = 'suppliedInput/combination_njets_Mar19.txt'

card.rapidity = AttrDict()
card.rapidity.hgg = 'suppliedInput/fromVittorio/differential_AbsRapidityNNLOPS_newBins_Nov28/Datacard_13TeV_differential_AbsRapidityNNLOPS_newBins_combination_renamedProcesses.txt'
card.rapidity.hzz = 'suppliedInput/fromDavid/YH_NNLOPS_Nov12/smH/hzz4l_comb_13TeV_xs.txt' # No OutsideAcceptance in hzz, so no renumbering necessary
card.rapidity.combination = 'suppliedInput/combination_rapidity_Mar19.txt'

card.ptjet = AttrDict()
card.ptjet.hgg = 'suppliedInput/fromVittorio/differential_Jet2p5Pt0NNLOPS_newBins_Nov28/Datacard_13TeV_differential_Jet2p5Pt0NNLOPS_newBins_renamedProcesses.txt'
card.ptjet.hzz = 'suppliedInput/fromDavid/PTJET_NNLOPS_Nov28/smH/hzz4l_comb_13TeV_xs.txt'
card.ptjet.combination = 'suppliedInput/combination_ptjet_Mar19.txt'

card.inclusive = AttrDict()
# Had unrenamed process names... this probably lead to unsimultaneous scaling of what should have been the same process!!
# card.inclusive.hgg = 'suppliedInput/fromVittorio/inclusive_Nov27/Datacard_13TeV_differential_InclusiveNNLOPS.txt'
card.inclusive.hgg = 'suppliedInput/fromVittorio/inclusive_Nov27/Datacard_13TeV_differential_InclusiveNNLOPS_renamedProcesses.txt'
card.inclusive.hzz = 'suppliedInput/fromDavid/differential_Nov27/smH/hzz4l_comb_13TeV_xs.txt'
# card.inclusive.combination = 'suppliedInput/combination_inclusive_Mar19.txt'
card.inclusive.combination = 'suppliedInput/combination_inclusive_May09.txt'

card.yukawa_noXHunc = AttrDict()
card.yukawa_noXHunc.hgg = 'suppliedInput/Yukawa_hgg_pth_ggH_Mar08.txt'
card.yukawa_noXHunc.hzz = 'suppliedInput/Yukawa_hzz_pth_ggH_Mar08.txt'
card.yukawa_noXHunc.combination = 'suppliedInput/Yukawa_combination_pth_ggH_Mar08.txt'
card.yukawa = AttrDict()
card.yukawa.hgg = 'suppliedInput/Yukawa_hgg_pth_ggH_Mar08_xHNuisPar.txt'
card.yukawa.hzz = 'suppliedInput/Yukawa_hzz_pth_ggH_Mar08_xHNuisPar.txt'
card.yukawa.combination = 'suppliedInput/Yukawa_combination_pth_ggH_Mar08_xHNuisPar.txt'

card.top = AttrDict()
card.top.nominal = AttrDict()
card.top.nominal.combWithHbb = 'suppliedInput/Top_combWithHbb_pth_ggH_Mar13_noBinsDropped.txt'
card.top.nominal.combination = 'suppliedInput/Top_combination_pth_ggH_Mar13_noBinsDropped.txt'
card.top.nominal.hbb = 'suppliedInput/Top_hbb_pth_ggH_Mar14_noBinsDropped.txt'
card.top.nominal.hgg = 'suppliedInput/Top_hgg_pth_ggH_Mar14_noBinsDropped.txt'
card.top.nominal.hzz = 'suppliedInput/Top_hzz_pth_ggH_Mar14_noBinsDropped.txt'
card.top.combWithHbb_last2BinsDropped = 'suppliedInput/Top_combWithHbb_pth_ggH_Mar12_last2BinsDropped.txt'
card.top.combination_last2BinsDropped = 'suppliedInput/Top_combination_pth_ggH_Mar12_last2BinsDropped.txt'

#____________________________________________________________________
# Workspaces

ws.pth_ggH = AttrDict()
ws.pth_ggH.hgg = 'out/workspaces_Mar02/ws_pth_ggH_hgg.root'
ws.pth_ggH.hzz = 'out/workspaces_Feb28/ws_pth_ggH_hzz.root'
ws.pth_ggH.hbb = 'out/workspaces_Mar02/ws_pth_ggH_hbb.root'
ws.pth_ggH.combination = 'out/workspaces_Mar01/ws_pth_ggH_combination.root'
# ws.pth_ggH.combWithHbb = 'out/workspaces_Mar02/ws_pth_ggH_combWithHbb.root'
ws.pth_ggH.combWithHbb = 'out/workspaces_May17/ws_pth_ggH_combWithHbb.root' # Contains a nuisance for xH

ws.pth_smH = AttrDict()
ws.pth_smH.hgg = 'out/workspaces_Mar14/ws_pth_smH_hgg.root'
ws.pth_smH.hzz = 'out/workspaces_Mar14/ws_pth_smH_hzz.root'
ws.pth_smH.hbb = 'out/workspaces_Mar14/ws_pth_smH_hbb.root'
ws.pth_smH.combination = 'out/workspaces_Mar14/ws_pth_smH_combination.root'
ws.pth_smH.combWithHbb = 'out/workspaces_Mar14/ws_pth_smH_combWithHbb.root'

ws.njets = AttrDict()
ws.njets.combination = 'out/workspaces_Mar19/ws_njets_combination.root'
ws.njets.hgg = 'out/workspaces_Mar19/ws_njets_hgg.root'
ws.njets.hzz = 'out/workspaces_Mar19/ws_njets_hzz.root'
ws.ptjet = AttrDict()
ws.ptjet.combination = 'out/workspaces_Mar19/ws_ptjet_combination.root'
ws.ptjet.hgg = 'out/workspaces_Mar19/ws_ptjet_hgg.root'
ws.ptjet.hzz = 'out/workspaces_Mar19/ws_ptjet_hzz.root'
ws.rapidity = AttrDict()
ws.rapidity.combination = 'out/workspaces_Mar19/ws_rapidity_combination.root'
ws.rapidity.hgg = 'out/workspaces_Mar19/ws_rapidity_hgg.root'
ws.rapidity.hzz = 'out/workspaces_Mar19/ws_rapidity_hzz.root'
ws.totalXS = AttrDict()
ws.totalXS.combination = 'out/workspaces_Apr12/combination_inclusive_Mar19_multiSignalModel.root'
ws.totalXS.hgg = 'out/workspaces_Apr12/Datacard_13TeV_differential_InclusiveNNLOPS_multiSignalModel.root'
ws.totalXS.hzz = 'out/workspaces_Apr12/hzz4l_comb_13TeV_xs_multiSignalModel.root'
ws.ratioOfBRs = 'out/workspaces_Apr12/combination_inclusive_Mar19_extendedMultiSignalModel_ratioOfBRs.root'

ws.top = AttrDict()
ws.top.nominal = AttrDict()
# ws.top.nominal.combination = 'out/workspaces_Mar03/combination_Top_reweighted_nominal.root'
# ws.top.nominal.combWithHbb = 'out/workspaces_Mar03/combWithHbb_Top_reweighted_nominal.root'
# ws.top.nominal.hgg = 'out/workspaces_Mar03/hgg_Top_reweighted_nominal.root'
# ws.top.nominal.hzz = 'out/workspaces_Mar03/hzz_Top_reweighted_nominal.root'
ws.top.nominal.combination = 'out/workspaces_Mar13/combination_Top_reweighted_noBinsDropped.root'
ws.top.nominal.combWithHbb = 'out/workspaces_Mar13/combWithHbb_Top_reweighted_noBinsDropped.root'
ws.top.nominal.hbb = 'out/workspaces_Mar14/hbb_Top_reweighted_noBinsDropped.root'
ws.top.nominal.hgg = 'out/workspaces_Mar14/hgg_Top_reweighted_noBinsDropped.root'
ws.top.nominal.hzz = 'out/workspaces_Mar14/hzz_Top_reweighted_noBinsDropped.root'

ws.top.lumiScale = 'out/workspaces_Apr11/combWithHbb_Top_reweighted_noBinsDropped_lumiScale.root'
ws.top.last2BinsDropped = AttrDict()
ws.top.last2BinsDropped.combWithHbb = 'out/workspaces_Mar12/combWithHbb_Top_reweighted_last2BinsDropped.root'
ws.top.last2BinsDropped.combination = 'out/workspaces_Mar12/combination_Top_reweighted_last2BinsDropped.root'
ws.top.BRcouplingDependency = 'combWithHbb_Top_reweighted_BRcouplingDependency.root'

ws.yukawa = AttrDict()
ws.yukawa.nominal = AttrDict()
ws.yukawa.nominal.combination = 'out/workspaces_Mar09/combination_Yukawa_reweighted_nominal.root'
ws.yukawa.nominal.hgg = 'out/workspaces_Mar09/hgg_Yukawa_reweighted_nominal.root'
ws.yukawa.nominal.hzz = 'out/workspaces_Mar15/hzz_Yukawa_reweighted_nominal.root'

ws.yukawa.couplingdependentBRs = AttrDict()
ws.yukawa.couplingdependentBRs.combination = 'out/workspaces_May22/combination_Yukawa_reweighted_scalingbbH_couplingdependentBRs.root'
ws.yukawa.couplingdependentBRs.hgg = 'out/workspaces_May26/hgg_Yukawa_reweighted_scalingbbH_couplingdependentBRs.root'
ws.yukawa.couplingdependentBRs.hzz = 'out/workspaces_May26/hzz_Yukawa_reweighted_scalingbbH_couplingdependentBRs.root'


ws.yukawa.unreweighted = AttrDict()
ws.yukawa.unreweighted.combination = 'out/workspaces_Mar09/combination_Yukawa_nominal.root'
ws.yukawa.unreweighted.hgg = 'out/workspaces_Mar09/hgg_Yukawa_nominal.root'
ws.yukawa.unreweighted.hzz = 'out/workspaces_Mar15/hzz_Yukawa_nominal.root'

ws.yukawa.lumiScale = 'out/workspaces_Mar12/combination_Yukawa_reweighted_lumiScale.root'
ws.yukawa.noTheoryUnc = 'out/workspaces_Mar12/combination_Yukawa_reweighted_noTheoryUnc.root'
ws.yukawa.profiledTotalXS = 'out/workspaces_Mar12/combination_Yukawa_reweighted_profiledTotalXS.root'
ws.yukawa.uncorrelatedTheoryUnc = 'out/workspaces_Mar12/combination_Yukawa_reweighted_uncorrelatedTheoryUnc.root'
ws.yukawa.BRcouplingDependency = 'out/workspaces_Mar20/combination_Yukawa_reweighted_BRcouplingDependency.root'

ws.topctcb = AttrDict()
ws.topctcb.nominal = AttrDict()
ws.topctcb.nominal.combWithHbb = 'out/workspaces_Mar20/combWithHbb_TopCtCb_reweighted_noBinsDropped.root'
ws.topctcb.nominal.hgg = 'out/workspaces_Mar20/hgg_TopCtCb_reweighted_noBinsDropped.root'
ws.topctcb.nominal.hzz = 'out/workspaces_Mar20/hzz_TopCtCb_reweighted_noBinsDropped.root'
ws.topctcb.lumiScale = 'out/workspaces_Apr11/combWithHbb_TopCtCb_reweighted_noBinsDropped_lumiScale.root'


#____________________________________________________________________
# Theory

theory.top = AttrDict()
# theory.top.filedir = 'out/derivedTheoryFiles_Dec07_TopHighPt'
# theory.top.correlation_matrix = 'out/correlationMatrices_Nov24_Top/corrMat_exp.txt'
# theory.top.uncertainties = 'out/correlationMatrices_Nov24_Top/errors_for_corrMat_exp.txt'
theory.top.filedir = 'out/theories_May17_tophighpt'
theory.top.correlation_matrix = ''
theory.top.uncertainties = ''


theory.yukawa = AttrDict()
theory.yukawa.filedir = 'out/theories_Mar09_yukawa_summed'
theory.yukawa.correlation_matrix = 'out/scalecorrelations_Mar09/corrMat_yukawa.txt'
theory.yukawa.correlation_matrix_uncorrelated = 'out/scalecorrelations_Mar09/corrMat_yukawa_uncorrelated.txt'
theory.yukawa.uncertainties = 'out/scalecorrelations_Mar09/errors_yukawa.txt'

theory.yukawa.filedir_gluoninduced = 'out/theories_Mar09_yukawa_gluoninduced'
theory.yukawa.filedir_quarkinduced = 'out/theories_Mar09_yukawa_quarkinduced'
theory.yukawa.filedir_quarkinduced_scaled = 'out/theories_Mar09_yukawa_quarkinduced_scaled'

#____________________________________________________________________
# Scans

scan.pth_ggH = AttrDict()
scan.pth_ggH.observed = AttrDict()
scan.pth_ggH.observed.hgg = [ 'out/Scan_Mar02_pth_ggH_hgg', 'out/Scan_Mar14_pth_ggH_hgg_rescan', 'out/Scan_Mar14_pth_ggH_hgg_rescan_0', 'out/Scan_Mar15_pth_ggH_hgg_rescan', 'out/Scan_Mar15_pth_ggH_hgg_rescan_0', 'out/Scan_Mar15_pth_ggH_hgg_rescan_1' ]
scan.pth_ggH.observed.hzz = [ 'out/Scan_Mar02_pth_ggH_hzz', 'out/Scan_Mar14_pth_ggH_hzz_rescan' ]
scan.pth_ggH.observed.hbb = 'out/Scan_Mar02_pth_ggH_hbb'
scan.pth_ggH.observed.combination = 'out/Scan_Mar02_pth_ggH_combination'
scan.pth_ggH.observed.combination_statonly = 'out/Scan_Mar20_pth_ggH_combination_statonly'
scan.pth_ggH.observed.combWithHbb_noxHunc = [ 'out/Scan_Mar02_pth_ggH_combWithHbb', 'out/Scan_Mar14_pth_ggH_combWithHbb_rescan', 'out/Scan_Mar14_pth_ggH_combWithHbb_rescan_0' ]
scan.pth_ggH.observed.combWithHbb_statonly_noxHunc = 'out/Scan_Mar20_pth_ggH_combWithHbb_statonly'
scan.pth_ggH.observed.combWithHbb = 'out/Scan_May17_pth_ggH_combWithHbb'
scan.pth_ggH.observed.combWithHbb_statonly = 'out/Scan_May17_pth_ggH_combWithHbb_statonly'
scan.pth_ggH.asimov = AttrDict()
scan.pth_ggH.asimov.hgg = [ 'out/Scan_Mar02_pth_ggH_hgg_asimov' ]
scan.pth_ggH.asimov.hzz = [ 'out/Scan_Mar02_pth_ggH_hzz_asimov', 'out/Scan_Mar14_pth_ggH_hzz_rescan_asimov' ]
scan.pth_ggH.asimov.hbb = 'out/Scan_Mar02_pth_ggH_hbb_asimov'
scan.pth_ggH.asimov.combination = 'out/Scan_Mar02_pth_ggH_combination_asimov'
scan.pth_ggH.asimov.combination_statonly = 'out/Scan_Mar20_pth_ggH_combination_statonly_asimov'
scan.pth_ggH.asimov.combWithHbb_noxHunc = [ 'out/Scan_Mar02_pth_ggH_combWithHbb_asimov', 'out/Scan_Mar14_pth_ggH_combWithHbb_rescan_asimov', 'out/Scan_Mar07_pth_ggH_combWithHbb_asimov' ]
scan.pth_ggH.asimov.combWithHbb_statonly_noxHunc = 'out/Scan_Mar20_pth_ggH_combWithHbb_statonly_asimov'
# scan.pth_ggH.asimov.combWithHbb = [ 'out/Scan_Mar02_pth_ggH_combWithHbb_asimov', 'out/Scan_Mar14_pth_ggH_combWithHbb_rescan_asimov', 'out/Scan_Mar07_pth_ggH_combWithHbb_asimov' ]
# scan.pth_ggH.asimov.combWithHbb_statonly = 'out/Scan_Mar20_pth_ggH_combWithHbb_statonly_asimov'

scan.pth_smH = AttrDict()
scan.pth_smH.observed = AttrDict()
scan.pth_smH.observed.hgg = [ 'out/Scan_Mar14_pth_smH_hgg', 'out/Scan_Mar14_pth_smH_hgg_rescan' ]
scan.pth_smH.observed.hzz = [ 'out/Scan_Mar14_pth_smH_hzz', 'out/Scan_Mar14_pth_smH_hzz_rescan' ]
scan.pth_smH.observed.hbb = 'out/Scan_Mar14_pth_smH_hbb'
scan.pth_smH.observed.combination = 'out/Scan_Mar14_pth_smH_combination'
scan.pth_smH.observed.combWithHbb = [ 'out/Scan_Mar14_pth_smH_combWithHbb', 'out/Scan_Mar14_pth_smH_combWithHbb_rescan', 'out/Scan_Mar14_pth_smH_combWithHbb_rescan_0' ]
scan.pth_smH.observed.combWithHbb_statonly = 'out/Scan_Mar20_pth_smH_combWithHbb_statonly'
scan.pth_smH.observed.combination_statonly = 'out/Scan_Mar20_pth_smH_combination_statonly'
scan.pth_smH.asimov = AttrDict()
scan.pth_smH.asimov.hgg = [ 'out/Scan_Mar14_pth_smH_hgg_asimov', 'out/Scan_Mar14_pth_smH_hgg_rescan_asimov' ]
scan.pth_smH.asimov.hzz = 'out/Scan_Mar14_pth_smH_hzz_asimov'
scan.pth_smH.asimov.hbb = 'out/Scan_Mar14_pth_smH_hbb_asimov'
scan.pth_smH.asimov.combination = 'out/Scan_Mar14_pth_smH_combination_asimov'
scan.pth_smH.asimov.combWithHbb = [ 'out/Scan_Mar14_pth_smH_combWithHbb_asimov', 'out/Scan_Mar14_pth_smH_combWithHbb_rescan_asimov' ]
scan.pth_smH.asimov.combWithHbb_statonly = 'out/Scan_Mar20_pth_smH_combWithHbb_statonly_asimov'
scan.pth_smH.asimov.combination_statonly = 'out/Scan_Mar20_pth_smH_combination_statonly_asimov'

scan.njets = AttrDict()
scan.njets.observed = AttrDict()
# scan.njets.observed.hgg = 'out/Scan_njets_Feb06_hgg'
# scan.njets.observed.hzz = 'out/Scan_njets_Feb06_hzz'
# scan.njets.observed.combination = 'out/Scan_njets_Feb06_combination'
scan.njets.observed.combination = 'out/Scan_Mar19_njets_combination'
scan.njets.observed.hgg = 'out/Scan_Mar19_njets_hgg'
scan.njets.observed.hzz = 'out/Scan_Mar19_njets_hzz'
scan.njets.observed.combination_statonly = 'out/Scan_Mar19_njets_combination_statonly'
scan.njets.asimov = AttrDict()
# scan.njets.asimov.hgg = 'out/Scan_njets_Feb06_hgg_asimov'
# scan.njets.asimov.hzz = 'out/Scan_njets_Feb06_hzz_asimov'
# scan.njets.asimov.combination = 'out/Scan_njets_Feb06_combination_asimov'
scan.njets.asimov.combination = 'out/Scan_Mar19_njets_combination_asimov'
scan.njets.asimov.hgg = 'out/Scan_Mar19_njets_hgg_asimov'
scan.njets.asimov.hzz = 'out/Scan_Mar19_njets_hzz_asimov'
scan.njets.asimov.combination_statonly = 'out/Scan_Mar19_njets_combination_statonly_asimov'

scan.ptjet = AttrDict()
scan.ptjet.observed = AttrDict()
# scan.ptjet.observed.hgg = 'out/Scan_ptjet_Feb06_hgg'
# scan.ptjet.observed.hzz = 'out/Scan_ptjet_Feb06_hzz'
# scan.ptjet.observed.combination = 'out/Scan_ptjet_Feb06_combination'
scan.ptjet.observed.combination = 'out/Scan_Mar19_ptjet_combination'
scan.ptjet.observed.hgg = 'out/Scan_Mar19_ptjet_hgg'
scan.ptjet.observed.hzz = 'out/Scan_Mar19_ptjet_hzz'
scan.ptjet.observed.combination_statonly = 'out/Scan_Mar19_ptjet_combination_statonly'
scan.ptjet.asimov = AttrDict()
# scan.ptjet.asimov.hgg = 'out/Scan_ptjet_Feb06_hgg_asimov'
# scan.ptjet.asimov.hzz = 'out/Scan_ptjet_Feb06_hzz_asimov'
# scan.ptjet.asimov.combination = 'out/Scan_ptjet_Feb06_combination_asimov'
scan.ptjet.asimov.combination = 'out/Scan_Mar19_ptjet_combination_asimov'
scan.ptjet.asimov.hgg = 'out/Scan_Mar19_ptjet_hgg_asimov'
scan.ptjet.asimov.hzz = 'out/Scan_Mar19_ptjet_hzz_asimov'
scan.ptjet.asimov.combination_statonly = 'out/Scan_Mar19_ptjet_combination_statonly_asimov'

scan.rapidity = AttrDict()
scan.rapidity.observed = AttrDict()
# scan.rapidity.observed.hgg = 'out/Scan_rapidity_Feb06_hgg'
# scan.rapidity.observed.hzz = 'out/Scan_rapidity_Feb06_hzz'
# scan.rapidity.observed.combination = 'out/Scan_rapidity_Feb06_combination'
scan.rapidity.observed.combination = 'out/Scan_Mar19_rapidity_combination'
scan.rapidity.observed.hgg = 'out/Scan_Mar19_rapidity_hgg'
scan.rapidity.observed.hzz = 'out/Scan_Mar19_rapidity_hzz'
scan.rapidity.observed.combination_statonly = 'out/Scan_Mar19_rapidity_combination_statonly'
scan.rapidity.asimov = AttrDict()
# scan.rapidity.asimov.hgg = 'out/Scan_rapidity_Feb06_hgg_asimov'
# scan.rapidity.asimov.hzz = 'out/Scan_rapidity_Feb06_hzz_asimov'
# scan.rapidity.asimov.combination = 'out/Scan_rapidity_Feb06_combination_asimov'
scan.rapidity.asimov.combination = 'out/Scan_Mar19_rapidity_combination_asimov'
scan.rapidity.asimov.hgg = 'out/Scan_Mar19_rapidity_hgg_asimov'
scan.rapidity.asimov.hzz = 'out/Scan_Mar19_rapidity_hzz_asimov'
scan.rapidity.asimov.combination_statonly = 'out/Scan_Mar19_rapidity_combination_statonly_asimov'

scan.totalXS = AttrDict()
scan.totalXS.combination = 'out/Scan_Apr12_totalXS_combination'
scan.totalXS.hgg = 'out/Scan_Apr12_totalXS_hgg'
scan.totalXS.hzz = 'out/Scan_Apr12_totalXS_hzz'
scan.totalXS.statonly = AttrDict()
scan.totalXS.statonly.combination = 'out/Scan_Apr12_totalXS_combination_statonly'
scan.totalXS.statonly.hgg = 'out/Scan_Apr12_totalXS_hgg_statonly'
scan.totalXS.statonly.hzz = 'out/Scan_Apr12_totalXS_hzz_statonly'

scan.ratioOfBRs = 'out/Scan_Apr12_ratioOfBRs_1'
scan.ratioOfBRs_statonly = 'out/Scan_Apr12_ratioOfBRs_1_statonly'

scan.top = AttrDict()
scan.top.reweighted = AttrDict()
scan.top.reweighted.observed = AttrDict()
scan.top.reweighted.observed.combWithHbb = 'out/Scan_Mar15_Top_combWithHbb'
scan.top.reweighted.observed.combination = 'out/Scan_Mar15_Top_combination'
scan.top.reweighted.observed.hgg = 'out/Scan_Mar15_Top_hgg'
scan.top.reweighted.observed.hzz = 'out/Scan_Mar15_Top_hzz'
scan.top.reweighted.asimov = AttrDict()
# scan.top.reweighted.asimov.hgg = 'out/Scan_Mar06_Top_hgg_asimov' # Redo properly for reweighted
# scan.top.reweighted.asimov.hzz = 'out/Scan_Mar06_Top_hzz_asimov' # Redo properly for reweighted
scan.top.reweighted.asimov.combWithHbb = 'out/Scan_Mar13_Top_combWithHbb_noBinsDropped_asimov'

# scan.top.lumi300fb = 'out/Scan_Mar13_Top_combWithHbb_lumiStudy_asimov' # Not sure what reweighting scheme is implemented here
scan.top.lumi300fb = 'out/Scan_Apr11_Top_combWithHbb_lumiStudy_asimov'
scan.top.profiledTotalXS = 'out/Scan_Mar13_Top_combWithHbb_profiledTotalXS_asimov'

scan.yukawa = AttrDict()
scan.yukawa.reweighted = AttrDict()
scan.yukawa.reweighted.observed = AttrDict()
scan.yukawa.reweighted.observed.combination = 'out/Scan_Yukawa_Mar09_combination'
scan.yukawa.reweighted.observed.hgg = 'out/Scan_Yukawa_Mar09_hgg'
scan.yukawa.reweighted.observed.hzz = 'out/Scan_Yukawa_Mar15_hzz'
scan.yukawa.reweighted.asimov = AttrDict()
scan.yukawa.reweighted.asimov.combination = 'out/Scan_Yukawa_Mar09_combination_asimov'
scan.yukawa.reweighted.asimov.hgg = 'out/Scan_Yukawa_Mar09_hgg_asimov'
scan.yukawa.reweighted.asimov.hzz = 'out/Scan_Yukawa_Mar15_hzz_asimov'

scan.yukawa.unreweighted = AttrDict()
scan.yukawa.unreweighted.observed = AttrDict()
scan.yukawa.unreweighted.observed.combination = 'out/Scan_Yukawa_Mar09_combination_0'
scan.yukawa.unreweighted.observed.hgg = 'out/Scan_Yukawa_Mar09_hgg_0'
scan.yukawa.unreweighted.observed.hzz = 'out/Scan_Yukawa_Mar15_hzz_0'
scan.yukawa.unreweighted.asimov = AttrDict()
scan.yukawa.unreweighted.asimov.combination = 'out/Scan_Yukawa_Mar09_combination_asimov_0'
scan.yukawa.unreweighted.asimov.hgg = 'out/Scan_Yukawa_Mar09_hgg_asimov_0'
scan.yukawa.unreweighted.asimov.hzz = 'out/Scan_Yukawa_Mar15_hzz_asimov_0'

scan.yukawa.lumi300fb = 'out/Scan_Yukawa_Mar12_combination_lumi300fb_asimov'
scan.yukawa.noTheoryUnc = 'out/Scan_Yukawa_Mar12_combination_noTheoryUnc_asimov'
scan.yukawa.profiledTotalXS = 'out/Scan_Yukawa_Mar12_combination_profiledTotalXS_asimov'
scan.yukawa.uncorrelatedTheoryUnc = 'out/Scan_Yukawa_Mar12_combination_uncorrelatedTheoryUnc_asimov'

scan.yukawa.onekappa = AttrDict()
scan.yukawa.onekappa.observed = AttrDict()
scan.yukawa.onekappa.observed.kappab = 'out/Scan_Yukawa_Mar16_combination_oneKappa_kappab_0'
scan.yukawa.onekappa.observed.kappac = 'out/Scan_Yukawa_Mar16_combination_oneKappa_kappac'
scan.yukawa.onekappa.asimov = AttrDict()
scan.yukawa.onekappa.asimov.kappab = 'out/Scan_Yukawa_Mar16_combination_oneKappa_kappab_asimov'
scan.yukawa.onekappa.asimov.kappac = 'out/Scan_Yukawa_Mar16_combination_oneKappa_kappac_asimov'

scan.topctcb = AttrDict()
scan.topctcb.reweighted = AttrDict()
scan.topctcb.reweighted.observed = AttrDict()
scan.topctcb.reweighted.observed.combWithHbb =  'out/Scan_Mar29_TopCtCb_combWithHbb'
scan.topctcb.reweighted.observed.hgg =  'out/Scan_Mar29_TopCtCb_hgg'
scan.topctcb.reweighted.observed.hzz =  'out/Scan_Mar29_TopCtCb_hzz'
scan.topctcb.reweighted.asimov = AttrDict()
scan.topctcb.reweighted.asimov.combWithHbb = 'out/Scan_Mar29_TopCtCb_combWithHbb_asimov'
scan.topctcb.reweighted.asimov.hgg = 'out/Scan_Mar29_TopCtCb_hgg_asimov'
scan.topctcb.reweighted.asimov.hzz = 'out/Scan_Mar29_TopCtCb_hzz_asimov'
scan.topctcb.lumi300fb = 'out/Scan_Apr11_TopCtCb_combWithHbb_lumiStudy_asimov'
