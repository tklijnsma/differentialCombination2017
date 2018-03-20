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

card.pth_ggH = AttrDict()
card.pth_ggH.hgg = 'suppliedInput/fromVittorio/pT_newBins_Feb28/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_newBins_renamedProcesses.txt'
card.pth_ggH.hzz = 'suppliedInput/fromDavid/PTH_Jan24_newBinning/ggH/hzz4l_comb_13TeV_xs_processesRenumbered.txt'
card.pth_ggH.hbb = 'suppliedInput/fromJavier/bernstein_r7428/comb_2017_ggHbb.txt'
card.pth_ggH.combination = 'suppliedInput/combination_pth_ggH_Mar01.txt'
card.pth_ggH.combWithHbb = 'suppliedInput/combWithHbb_pth_ggH_Mar02.txt'

card.pth_smH = AttrDict()
card.pth_smH.hgg = card.pth_ggH.hgg # scale ggH/xH by smH!
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
card.inclusive.hgg = 'suppliedInput/fromVittorio/inclusive_Nov27/Datacard_13TeV_differential_InclusiveNNLOPS.txt'
card.inclusive.hzz = 'suppliedInput/fromDavid/differential_Nov27/smH/hzz4l_comb_13TeV_xs.txt'
card.inclusive.combination = 'suppliedInput/combination_inclusive_Mar19.txt'

card.yukawa = AttrDict()
card.yukawa.hgg = 'suppliedInput/Yukawa_hgg_pth_ggH_Mar08.txt'
card.yukawa.hzz = 'suppliedInput/Yukawa_hzz_pth_ggH_Mar08.txt'
card.yukawa.combination = 'suppliedInput/Yukawa_combination_pth_ggH_Mar08.txt'

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
ws.pth_ggH.combWithHbb = 'out/workspaces_Mar02/ws_pth_ggH_combWithHbb.root'

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

ws.top.last2BinsDropped = AttrDict()
ws.top.last2BinsDropped.combWithHbb = 'out/workspaces_Mar12/combWithHbb_Top_reweighted_last2BinsDropped.root'
ws.top.last2BinsDropped.combination = 'out/workspaces_Mar12/combination_Top_reweighted_last2BinsDropped.root'

ws.yukawa = AttrDict()
ws.yukawa.nominal = AttrDict()
ws.yukawa.nominal.combination = 'out/workspaces_Mar09/combination_Yukawa_reweighted_nominal.root'
ws.yukawa.nominal.hgg = 'out/workspaces_Mar09/hgg_Yukawa_reweighted_nominal.root'
ws.yukawa.nominal.hzz = 'out/workspaces_Mar15/hzz_Yukawa_reweighted_nominal.root'

ws.yukawa.unreweighted = AttrDict()
ws.yukawa.unreweighted.combination = 'out/workspaces_Mar09/combination_Yukawa_nominal.root'
ws.yukawa.unreweighted.hgg = 'out/workspaces_Mar09/hgg_Yukawa_nominal.root'
ws.yukawa.unreweighted.hzz = 'out/workspaces_Mar15/hzz_Yukawa_nominal.root'

ws.yukawa.lumiScale = 'out/workspaces_Mar12/combination_Yukawa_reweighted_lumiScale.root'
ws.yukawa.noTheoryUnc = 'out/workspaces_Mar12/combination_Yukawa_reweighted_noTheoryUnc.root'
ws.yukawa.profiledTotalXS = 'out/workspaces_Mar12/combination_Yukawa_reweighted_profiledTotalXS.root'
ws.yukawa.uncorrelatedTheoryUnc = 'out/workspaces_Mar12/combination_Yukawa_reweighted_uncorrelatedTheoryUnc.root'


#____________________________________________________________________
# Theory

theory.top = AttrDict()
theory.top.filedir = 'out/derivedTheoryFiles_Dec07_TopHighPt'
theory.top.correlation_matrix = 'out/correlationMatrices_Nov24_Top/corrMat_exp.txt'
theory.top.uncertainties = 'out/correlationMatrices_Nov24_Top/errors_for_corrMat_exp.txt'

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
scan.pth_ggH.observed.combWithHbb = [ 'out/Scan_Mar02_pth_ggH_combWithHbb', 'out/Scan_Mar14_pth_ggH_combWithHbb_rescan', 'out/Scan_Mar14_pth_ggH_combWithHbb_rescan_0' ]
scan.pth_ggH.observed.combWithHbb_statonly = 'out/Scan_Mar20_pth_ggH_combWithHbb_statonly'
scan.pth_ggH.observed.combination_statonly = 'out/Scan_Mar20_pth_ggH_combination_statonly'
scan.pth_ggH.asimov = AttrDict()
scan.pth_ggH.asimov.hgg = [ 'out/Scan_Mar02_pth_ggH_hgg_asimov' ]
scan.pth_ggH.asimov.hzz = [ 'out/Scan_Mar02_pth_ggH_hzz_asimov', 'out/Scan_Mar14_pth_ggH_hzz_rescan_asimov' ]
scan.pth_ggH.asimov.hbb = 'out/Scan_Mar02_pth_ggH_hbb_asimov'
scan.pth_ggH.asimov.combination = 'out/Scan_Mar02_pth_ggH_combination_asimov'
scan.pth_ggH.asimov.combWithHbb = [ 'out/Scan_Mar02_pth_ggH_combWithHbb_asimov', 'out/Scan_Mar14_pth_ggH_combWithHbb_rescan_asimov', 'out/Scan_Mar07_pth_ggH_combWithHbb_asimov' ]
scan.pth_ggH.asimov.combWithHbb_statonly = 'out/Scan_Mar20_pth_ggH_combWithHbb_statonly_asimov'
scan.pth_ggH.asimov.combination_statonly = 'out/Scan_Mar20_pth_ggH_combination_statonly_asimov'

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

scan.top.lumi300fb = 'out/Scan_Mar13_Top_combWithHbb_lumiStudy_asimov'
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


# # ----- PTH -----
# card_hgg_smH_PTH_unprocessed   = 'suppliedInput/fromVittorio/pT_NNLOPS_Nov01/Datacard_13TeV_differential_PtNNLOPS_systs.txt'
# card_hzz_smH_PTH_unprocessed   = 'suppliedInput/fromDavid/PTH_NNLOPS_Nov01/smH/hzz4l_comb_13TeV_xs.txt'

# card_hgg_ggHxH_PTH_unprocessed = 'suppliedInput/fromVittorio/pT_NNLOPS_ggHxH_Nov03/Datacard_13TeV_differential_PtGghPlusHxNNLOPS.txt'
# card_hzz_ggHxH_PTH_unprocessed = 'suppliedInput/fromDavid/PTH_NNLOPS_Nov01/ggH/hzz4l_comb_13TeV_xs.txt'

# # card_hbb_ggHxH_PTH             = 'suppliedInput/fromJavier/cards_PTH_diff_r7402_Dec14/comb_2017_ggHbb.txt'
# # card_hbb_ggHxH_PTH_debuggingTest = 'suppliedInput/combinedCard_hbb_debuggingTest_Dec19.txt'
# # card_combinedWithHbb_ggHxH_PTH = 'suppliedInput/combinedCard_hgg_hzz_hbb_ggHxH_Dec15.txt'

# # card_hbb_ggHxH_PTH             = 'suppliedInput/fromJavier/bernstein_r7428/comb_2017_ggHbb.txt'
# # card_combinedWithHbb_ggHxH_PTH = 'suppliedInput/combinedCard_hgg_hzz_hbb_ggHxH_Dec21.txt'
# card_hbb_ggHxH_PTH             = 'suppliedInput/fromJavier/negativity_fixed_Jan19/comb_2017_ggHbb.txt'
# card_combinedWithHbb_ggHxH_PTH = 'suppliedInput/combinedCard_hgg_hzz_hbb_ggHxH_Jan19.txt'
# card_hgg_ggHxH_PTH             = 'suppliedInput/fromVittorio/pT_NNLOPS_ggHxH_Nov03/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses.txt'
# card_hzz_ggHxH_PTH             = 'suppliedInput/fromDavid/PTH_NNLOPS_Nov01/ggH/hzz4l_comb_13TeV_xs_processesShifted.txt'
# card_combined_ggHxH_PTH        = 'suppliedInput/combinedCard_Nov03.txt'

# card_hgg_smH_PTH               = 'suppliedInput/fromVittorio/pT_NNLOPS_Nov01/Datacard_13TeV_differential_PtNNLOPS_systs_renamedProcesses.txt'
# card_hzz_smH_PTH               = 'suppliedInput/fromDavid/PTH_NNLOPS_Nov01/smH/hzz4l_comb_13TeV_xs_processesShifted.txt'
# card_combined_smH_PTH          = 'suppliedInput/combinedCard_smH_Nov07.txt'

# # ----- NJ -----
# card_hgg_smH_NJ_unprocessed    = 'suppliedInput/fromVittorio/differential_Njets2p5NNLOPS_Nov10/Datacard_13TeV_differential_Njets2p5NNLOPS.txt'
# card_hzz_smH_NJ                = 'suppliedInput/fromDavid/NJ_NNLOPS_Nov01/smH/hzz4l_comb_13TeV_xs.txt'
# card_hgg_smH_NJ                = 'suppliedInput/fromVittorio/differential_Njets2p5NNLOPS_Nov10/Datacard_13TeV_differential_Njets2p5NNLOPS_renamedProcesses.txt'
# card_combined_smH_NJ           = 'suppliedInput/combinedCard_NJ_smH_Nov12.txt'

# # ----- YH -----
# card_hgg_smH_YH_unprocessed    = 'suppliedInput/fromVittorio/differential_AbsRapidityNNLOPS_newBins_Nov28/Datacard_13TeV_differential_AbsRapidityNNLOPS_newBins_combination.txt'
# card_hgg_smH_YH                = 'suppliedInput/fromVittorio/differential_AbsRapidityNNLOPS_newBins_Nov28/Datacard_13TeV_differential_AbsRapidityNNLOPS_newBins_combination_renamedProcesses.txt'
# card_hzz_smH_YH                = 'suppliedInput/fromDavid/YH_NNLOPS_Nov12/smH/hzz4l_comb_13TeV_xs.txt' # No OutsideAcceptance in hzz, so no renumbering necessary
# card_combined_smH_YH           = 'suppliedInput/combinedCard_YH_smH_Nov28.txt'

# # ----- PTJ -----
# card_hzz_smH_PTJ               = 'suppliedInput/fromDavid/PTJET_NNLOPS_Nov28/smH/hzz4l_comb_13TeV_xs.txt'
# card_hgg_smH_PTJ_unprocessed   = 'suppliedInput/fromVittorio/differential_Jet2p5Pt0NNLOPS_newBins_Nov28/Datacard_13TeV_differential_Jet2p5Pt0NNLOPS_newBins.txt'
# card_hgg_smH_PTJ               = 'suppliedInput/fromVittorio/differential_Jet2p5Pt0NNLOPS_newBins_Nov28/Datacard_13TeV_differential_Jet2p5Pt0NNLOPS_newBins_renamedProcesses.txt'
# card_combined_smH_PTJ          = 'suppliedInput/combinedCard_PTJ_smH_Nov28.txt'

# # ----- INCLUSIVE -----
# card_hgg_INC_unprocessed       = 'suppliedInput/fromVittorio/inclusive_Nov27/Datacard_13TeV_differential_InclusiveNNLOPS.txt'
# card_hzz_INC_unprocessed       = 'suppliedInput/fromDavid/differential_Nov27/smH/hzz4l_comb_13TeV_xs.txt'
# card_combined_INC              = 'suppliedInput/combinedCard_smH_Nov27_INCLUSIVE.txt'

# # ----- Premature implementation of new binning -----
# card_hzz_ggHxH_PTH_newBins_unprocessed = 'suppliedInput/fromDavid/PTH_Jan24_newBinning/ggH/hzz4l_comb_13TeV_xs.txt'
# card_hzz_ggHxH_PTH_newBins      = 'suppliedInput/fromDavid/PTH_Jan24_newBinning/ggH/hzz4l_comb_13TeV_xs_processesShifted.txt'
# card_combined_ggHxH_PTH_newBins = 'suppliedInput/combinedCard_newBins_hzz_hbb_ggHxH_Jan24.txt'


# ########################################
# # Derived theory files
# ########################################

# derivedTheoryFiles_YukawaQuarkInduced       = 'out/derivedTheoryFiles_Nov03_YukawaQuarkInduced'
# derivedTheoryFiles_YukawaGluonInduced       = 'out/derivedTheoryFiles_Nov03_YukawaGluonInduced'
# derivedTheoryFiles_YukawaQuarkInducedScaled = 'out/derivedTheoryFiles_Nov03_YukawaQuarkInducedScaled'
# derivedTheoryFiles_YukawaSummed             = 'out/derivedTheoryFiles_Nov03_YukawaSummed'

# # Problem with SM file
# # derivedTheoryFiles_Top                      = 'out/derivedTheoryFiles_Nov06_Top'
# derivedTheoryFiles_Top                      = 'out/derivedTheoryFiles_Nov24_Top'

# derivedTheoryFiles_TopHighPt                = 'out/derivedTheoryFiles_Dec07_TopHighPt'


# ########################################
# # Correlation matrices and uncertainties
# ########################################

# correlationMatrix_Yukawa              = 'out/correlationMatrices_Nov03/corrMat_exp.txt'
# theoryUncertainties_Yukawa            = 'out/correlationMatrices_Nov03/errors_for_corrMat_exp.txt'
# correlationMatrix_Yukawa_Uncorrelated = 'out/correlationMatrices_Nov03/corrMat_exp_UNCORRELATED.txt'

# # Problem with SM file
# # correlationMatrix_Top                 = 'out/correlationMatrices_Nov06_Top/corrMat_exp.txt'
# # theoryUncertainties_Top               = 'out/correlationMatrices_Nov06_Top/errors_for_corrMat_exp.txt'
# correlationMatrix_Top                 = 'out/correlationMatrices_Nov24_Top/corrMat_exp.txt'
# theoryUncertainties_Top               = 'out/correlationMatrices_Nov24_Top/errors_for_corrMat_exp.txt'


# # Correlation matrices for differential observables
# # correlationMatrix_PTH                 = 'out/corrMat_Nov15_combinedCard_smH_Nov07/higgsCombine_CORRMAT_combinedCard_smH_Nov07.MultiDimFit.mH125.root'
# # correlationMatrix_NJ                  = 'out/corrMat_Nov15_combinedCard_NJ_smH_Nov12/higgsCombine_CORRMAT_combinedCard_NJ_smH_Nov12.MultiDimFit.mH125.root'
# correlationMatrix_PTH                   = 'out/corrMat_Nov29_combinedCard_smH_Nov07/higgsCombine_CORRMAT_combinedCard_smH_Nov07.MultiDimFit.mH125.root'
# correlationMatrix_PTH_ggH               = 'out/corrMat_Nov29_combinedCard_Nov03_xHfixed_1/higgsCombine_CORRMAT_combinedCard_Nov03_xHfixed.MultiDimFit.mH125.root'
# correlationMatrix_NJ                    = 'out/corrMat_Nov29_combinedCard_NJ_smH_Nov12/higgsCombine_CORRMAT_combinedCard_NJ_smH_Nov12.MultiDimFit.mH125.root'
# correlationMatrix_YH                    = 'out/corrMat_Nov29_combinedCard_YH_smH_Nov28/higgsCombine_CORRMAT_combinedCard_YH_smH_Nov28.MultiDimFit.mH125.root'
# correlationMatrix_PTJ                   = 'out/corrMat_Nov30_combinedCard_PTJ_smH_Nov28/higgsCombine_CORRMAT_combinedCard_PTJ_smH_Nov28.MultiDimFit.mH125.root'


# ########################################
# # Workspaces after running t2ws
# ########################################

# # ======================================
# # Unsplit workspaces for non-coupling combinations

# # ---------------------
# # Workspaces for extra studies

# # ws_ratio_of_BRs              = 'out/workspaces_Nov08/combinedCard_smH_Nov07_FitBRModel.root' # Probably not wrong but safer to use new result
# # ws_ratio_of_BRs              = 'out/workspaces_Dec11/combinedCard_smH_Nov07_FitBRModel.root'
# # ws_ratio_of_BRs_globalScales = 'out/workspaces_Nov20/combinedCard_smH_Nov07_FitBRModel_globalScales.root'
# ws_totalXS                   = 'out/workspaces_Nov08/combinedCard_smH_Nov07_FitBRModel_fitTotalXS.root'
# ws_combined_ratioOfBRs       = 'out/workspaces_Nov14/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_ratioOfBRs.root'
# ws_combined_totalXS          = 'out/workspaces_Nov27/combinedCard_smH_Nov27_INCLUSIVE_FitBRModel_fitTotalXS.root'

# ws_ratio_of_BRs_NJ           = 'out/workspaces_Dec20/combinedCard_NJ_smH_Nov12_FitBRModel.root'
# ws_ratio_of_BRs_PTJ          = 'out/workspaces_Dec20/combinedCard_PTJ_smH_Nov28_FitBRModel.root'
# ws_ratio_of_BRs_YH           = 'out/workspaces_Dec20/combinedCard_YH_smH_Nov28_FitBRModel.root'
# ws_ratio_of_BRs_PTH          = 'out/workspaces_Dec20/combinedCard_smH_Nov07_FitBRModel.root'
# ws_ratio_of_BRs_globalScales = 'out/workspaces_Dec20/combinedCard_smH_Nov27_INCLUSIVE_FitBRModel_globalScales.root'

# # ======================================
# # kappab kappac

# # With old model, still valid though. Unreweighted.
# # ws_combined_Yukawa = 'out/workspaces_Nov03/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties.root'
# # ws_hgg_Yukawa      = 'out/workspaces_Nov03/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_Yukawa_withTheoryUncertainties.root'
# # ws_hzz_Yukawa      = 'out/workspaces_Nov03/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_Yukawa_withTheoryUncertainties.root'

# # New model, reweighted
# # ws_combined_Yukawa = 'out/workspaces_Jan19/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties.root'
# # ws_combined_Yukawa = 'out/workspaces_Feb01/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_reweighted.root'
# # ws_hgg_Yukawa      = 'out/workspaces_Jan19/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_Yukawa_withTheoryUncertainties.root'
# # ws_hzz_Yukawa      = 'out/workspaces_Jan19/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_Yukawa_withTheoryUncertainties.root'

# # Redid t2ws to be sure; unreweighted:
# ws_combined_Yukawa = 'out/workspaces_Feb01/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties.root'
# ws_hgg_Yukawa      = 'out/workspaces_Feb01/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_Yukawa_withTheoryUncertainties.root'
# ws_hzz_Yukawa      = 'out/workspaces_Feb01/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_Yukawa_withTheoryUncertainties.root'

# # Reweighted:
# ws_combined_Yukawa_reweighted = 'out/workspaces_Feb01/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_reweighted.root'

# # ws_combined_Yukawa_profiledTotalXS                     = 'out/workspaces_Nov09/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_profiledTotalXS.root'
# # ws_combined_Yukawa_profiledTotalXS                     = 'out/workspaces_Feb05/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_profiledTotalXS_reweighted.root'
# ws_combined_Yukawa_profiledTotalXS                     = 'out/workspaces_Feb07/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_profiledTotalXS_reweighted.root'

# ws_combined_Yukawa_noTheoryUncertainties               = 'out/workspaces_Nov06/combinedCard_Nov03_CouplingModel_Yukawa_noTheoryUncertainties.root'
# ws_combined_Yukawa_withUncorrelatedTheoryUncertainties = 'out/workspaces_Nov06/combinedCard_Nov03_CouplingModel_Yukawa_withUncorrelatedTheoryUncertainties.root'
# ws_combined_Yukawa_lumiScalable                        = 'out/workspaces_Nov06/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_lumiScale.root'
# # ws_combined_Yukawa_couplingDependentBR                 = 'out/workspaces_Nov10/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_couplingDependentBR.root'
# # ws_combined_Yukawa_profiledTotalXS_fitOnlyNormalization = 'out/workspaces_Nov17/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_profiledTotalXS_fitOnlyNormalization.root'
# # ws_combined_Yukawa_profiledTotalXS_fitOnlyNormalization = 'out/workspaces_Nov20/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_profiledTotalXS_fitOnlyNormalization.root'
# ws_combined_Yukawa_profiledTotalXS_fitOnlyNormalization = 'out/workspaces_Nov28/combinedCard_smH_Nov27_INCLUSIVE_CouplingModel_Yukawa_profiledTotalXS_fitOnlyNormalization.root'


# # Using the new model
# ws_combined_Yukawa_couplingDependentBR_reweighted = 'out/workspaces_Jan25/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_couplingDependentBR_reweighted.root'
# ws_combined_Yukawa_couplingDependentBR = 'out/workspaces_Feb07/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_couplingDependentBR.root'
# ws_combined_Yukawa_couplingDependentBR_profiledTotalXS = 'out/workspaces_Feb07/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_couplingDependentBR_profiledTotalXS.root'


# # ======================================
# # kappat kappag

# # Problem with last theory uncertainty:
# # ws_combined_Top    = 'out/workspaces_Nov06/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties.root'
# # ws_hgg_Top         = 'out/workspaces_Nov06/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_Top_withTheoryUncertainties.root'
# # ws_hzz_Top         = 'out/workspaces_Nov06/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_Top_withTheoryUncertainties.root'

# # Last theory uncertainty skipped:  # Had problem with SM file
# # ws_hgg_Top         = 'out/workspaces_Nov08/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_Top_withTheoryUncertainties.root'
# # ws_combined_Top    = 'out/workspaces_Nov08/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties.root'
# # ws_hzz_Top         = 'out/workspaces_Nov08/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_Top_withTheoryUncertainties.root'

# # Problem with SM file resolved
# ws_combined_Top    = 'out/workspaces_Nov24/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties.root'
# ws_hgg_Top         = 'out/workspaces_Nov24/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_Top_withTheoryUncertainties.root'
# ws_hzz_Top         = 'out/workspaces_Nov24/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_Top_withTheoryUncertainties.root'

# # ws_combined_TopHighPt = 'out/workspaces_Dec11/combinedCard_Nov03_CouplingModel_TopHighPt_withTheoryUncertainties.root'
# # ws_hgg_TopHighPt      = 'out/workspaces_Dec11/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_TopHighPt_withTheoryUncertainties.root'
# # ws_hzz_TopHighPt      = 'out/workspaces_Dec11/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_TopHighPt_withTheoryUncertainties.root'

# # Reweighted
# ws_combined_TopHighPt = 'out/workspaces_Jan19/combinedCard_Nov03_CouplingModel_TopHighPt_withTheoryUncertainties_reweighted.root'
# ws_hgg_TopHighPt      = 'out/workspaces_Jan23/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_TopHighPt_withTheoryUncertainties_reweighted.root'
# ws_hzz_TopHighPt      = 'out/workspaces_Jan23/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_TopHighPt_withTheoryUncertainties_reweighted.root'
# ws_combWithHbb_TopHighPt = 'out/workspaces_Jan19/combinedCard_hgg_hzz_hbb_ggHxH_Jan19_CouplingModel_TopHighPt_withTheoryUncertainties_reweighted.root'

# # ws_combined_Top_lumiScalable           = 'out/workspaces_Nov09/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_lumiScale.root'
# # ws_combined_Top_profiledTotalXS        = 'out/workspaces_Nov09/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_profiledTotalXS.root'

# ws_combined_Top_couplingDependentBR    = 'out/workspaces_Nov10/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_couplingDependentBR.root'
# ws_combined_Top_couplingDependentBR_profiledTotalXS = 'out/workspaces_Feb07/combinedCard_Nov03_CouplingModel_TopHighPt_withTheoryUncertainties_profiledTotalXS_couplingDependentBR.root'

# ws_combined_Top_skippedLastBin         = 'out/workspaces_Nov17/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_skippedLastBin.root'

# # ws_combined_Top_profiledTotalXS_fitOnlyNormalization = 'out/workspaces_Nov17/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_profiledTotalXS_fitOnlyNormalization.root'
# ws_combined_Top_profiledTotalXS_fitOnlyNormalization = 'out/workspaces_Nov28/combinedCard_smH_Nov27_INCLUSIVE_CouplingModel_Top_profiledTotalXS_fitOnlyNormalization.root'
# ws_combined_TopHighPt_profiledTotalXS_fitOnlyNormalization = 'out/workspaces_Dec12/combinedCard_smH_Nov27_INCLUSIVE_CouplingModel_Top_profiledTotalXS_fitOnlyNormalization.root'

# # Without problematic SM norm
# ws_combined_Top_lumiScalable           = 'out/workspaces_Nov27/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_lumiScale.root'
# ws_combined_Top_profiledTotalXS        = 'out/workspaces_Nov27/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_profiledTotalXS.root'


# ws_combined_TopCtCb = 'out/workspaces_Nov15/combinedCard_Nov03_CouplingModel_TopCtCb_withTheoryUncertainties.root'
# ws_hgg_TopCtCb      = 'out/workspaces_Nov15/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_TopCtCb_withTheoryUncertainties.root'
# ws_hzz_TopCtCb      = 'out/workspaces_Nov15/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_TopCtCb_withTheoryUncertainties.root'

# ws_combined_TopCtCbHighPt = 'out/workspaces_Dec12/combinedCard_Nov03_CouplingModel_TopCtCbHighPt_withTheoryUncertainties.root'
# ws_hgg_TopCtCbHighPt      = 'out/workspaces_Dec12/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_TopCtCbHighPt_withTheoryUncertainties.root'
# ws_hzz_TopCtCbHighPt      = 'out/workspaces_Dec12/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_TopCtCbHighPt_withTheoryUncertainties.root'


# ########################################
# # Scans
# ########################################

# # Non-coupling combinations Moved to EOF

# # ======================================
# # Extra studies

# # scan_ratioOfBRs                  = 'out/Scan_ratioOfBRs_Nov08' # Probably not wrong, but renewed to be safe
# scan_ratioOfBRs                  = 'out/Scan_ratioOfBRs_Dec11'
# scan_ratioOfBRs_globalScales     = 'out/Scan_ratioOfBRs_Nov20_globalScales'
# # scan_combined_totalXS            = 'out/Scan_TotalXS_Nov08'
# # scan_combined_PTH_xHfixed        = 'out/Scan_PTH_Nov08_xHfixed'
# # scan_combined_PTH_xHfixed_asimov = 'out/Scan_PTH_Nov17_xHfixed_asimov'

# scan_combined_totalXS            = 'out/Scan_TotalXS_Nov27'

# scan_ratioOfBRs_NJ   = 'out/Scan_ratioOfBRs_Dec20'
# scan_ratioOfBRs_PTJ  = 'out/Scan_ratioOfBRs_Dec20_0'
# scan_ratioOfBRs_YH   = 'out/Scan_ratioOfBRs_Dec20_1'
# scan_ratioOfBRs_PTH  = 'out/Scan_ratioOfBRs_Dec20_2'
# scan_ratioOfBRs_INC  = 'out/Scan_ratioOfBRs_Dec20_3'

# scan_ratioOfBRs_NJ_asimov   = 'out/Scan_ratioOfBRs_Dec20_asimov'
# scan_ratioOfBRs_PTJ_asimov  = 'out/Scan_ratioOfBRs_Dec20_asimov_0'
# scan_ratioOfBRs_YH_asimov   = 'out/Scan_ratioOfBRs_Dec20_asimov_1'
# scan_ratioOfBRs_PTH_asimov  = 'out/Scan_ratioOfBRs_Dec20_asimov_2'
# scan_ratioOfBRs_INC_asimov  = 'out/Scan_ratioOfBRs_Dec20_asimov_3'



# # ======================================
# # kappab kappac scans

# scan_combined_Yukawa_old                                        = 'out/Scan_Yukawa_Nov03_0'
# scan_hgg_Yukawa_old                                             = 'out/Scan_Yukawa_Nov03_hgg'
# scan_hzz_Yukawa_old                                             = 'out/Scan_Yukawa_Nov03_hzz'
# scan_combined_Yukawa_asimov_old                                 = 'out/Scan_Yukawa_Nov03_asimov'
# scan_hgg_Yukawa_asimov_old                                      = 'out/Scan_Yukawa_Nov03_hgg_asimov'
# scan_hzz_Yukawa_asimov_old                                      = 'out/Scan_Yukawa_Nov03_hzz_asimov'

# scan_combined_Yukawa                                            = 'out/Scan_Yukawa_Feb01_0'
# scan_hgg_Yukawa                                                 = 'out/Scan_Yukawa_Feb01_hgg'
# scan_hzz_Yukawa                                                 = 'out/Scan_Yukawa_Feb01_hzz'
# scan_combined_Yukawa_asimov                                     = 'out/Scan_Yukawa_Feb01_asimov_0'
# scan_hgg_Yukawa_asimov                                          = 'out/Scan_Yukawa_Feb01_hgg_asimov'
# scan_hzz_Yukawa_asimov                                          = 'out/Scan_Yukawa_Feb01_hzz_asimov'

# scan_combined_Yukawa_reweighted                                 = 'out/Scan_Yukawa_Feb01'
# scan_combined_Yukawa_reweighted_asimov                          = 'out/Scan_Yukawa_Feb01_asimov'


# # scan_combined_Yukawa_fitOnlyNormalization_asimov            = 'out/Scan_Yukawa_Nov17_profiledTotalXS_fitOnlyNormalization_asimov'
# scan_combined_Yukawa_fitOnlyNormalization                   = 'out/Scan_Yukawa_Nov20_profiledTotalXS_fitOnlyNormalization'
# # scan_combined_Yukawa_fitOnlyNormalization_asimov            = 'out/Scan_Yukawa_Nov20_profiledTotalXS_fitOnlyNormalization_asimov'
# scan_combined_Yukawa_fitOnlyNormalization_asimov            = 'out/Scan_Yukawa_Nov28_fitOnlyNormalization_asimov_1'

# scan_combined_Yukawa_lumiStudy_asimov                       = 'out/Scan_Yukawa_Nov06_lumiStudy_asimov'
# scan_combined_Yukawa_noTheoryUncertainties_asimov           = 'out/Scan_Yukawa_Nov06_noTheoryUncertainties_asimov'
# scan_combined_Yukawa_uncorrelatedTheoryUncertainties_asimov = 'out/Scan_Yukawa_Nov06_uncorrelatedTheoryUncertainties_asimov'

# # Now reweighted
# # scan_combined_Yukawa_profiledTotalXS_asimov                 = 'out/Scan_Yukawa_Nov09_profiledTotalXS_asimov'
# scan_combined_Yukawa_profiledTotalXS_asimov                 = 'out/Scan_Yukawa_Feb05_combined_profiledTotalXS_asimov'

# # Reweighting wrong way around, and only half applied
# # scan_combined_Yukawa_oneKappa_kappac                        = 'out/Scan_Yukawa_Nov06_oneKappa_kappac'
# # scan_combined_Yukawa_oneKappa_kappac_asimov                 = 'out/Scan_Yukawa_Nov06_oneKappa_kappac_asimov'
# # scan_combined_Yukawa_oneKappa_kappab                        = 'out/Scan_Yukawa_Jan30_oneKappa_kappab_0'
# # scan_combined_Yukawa_oneKappa_kappab_asimov                 = 'out/Scan_Yukawa_Nov06_oneKappa_kappab_asimov'
# # Corrected
# scan_combined_Yukawa_oneKappa_kappac                        = 'out/Scan_Yukawa_Feb07_combined_oneKappa_kappac'
# scan_combined_Yukawa_oneKappa_kappac_asimov                 = 'out/Scan_Yukawa_Feb07_combined_oneKappa_kappac_asimov'
# scan_combined_Yukawa_oneKappa_kappab                        = 'out/Scan_Yukawa_Feb07_combined_oneKappa_kappab'
# scan_combined_Yukawa_oneKappa_kappab_asimov                 = 'out/Scan_Yukawa_Feb07_combined_oneKappa_kappab_asimov'

# scan_combined_Yukawa_couplingDependentBR_asimov              = 'out/Scan_Yukawa_Nov10_couplingDependentBR_asimov'
# scan_combined_Yukawa_couplingDependentBR_fixedKappaV_asimov  = 'out/Scan_Yukawa_Nov10_couplingDependentBR_fixedKappaV_asimov'
# # scan_combined_Yukawa_couplingDependentBR_kappaVMaxOne_asimov = 'out/Scan_Yukawa_Nov16_couplingDependentBR_kappaVMaxOne_asimov'
# scan_combined_Yukawa_couplingDependentBR_kappaVMaxOne_asimov = 'out/Scan_Yukawa_Nov17_couplingDependentBR_kappaVMaxOne_asimov'
# scan_combined_Yukawa_ratioOfBRs_asimov                       = 'out/Scan_Yukawa_Nov14_ratioOfBRs_asimov'

# scan_combined_Yukawa_ratioOfBRs_onedimRatioScan              = 'out/Scan_Yukawa_Nov15_ratioOfBRs_onedimRatioScan_0'
# scan_combined_Yukawa_ratioOfBRs_onedimRatioScan_asimov       = 'out/Scan_Yukawa_Nov15_ratioOfBRs_onedimRatioScan_asimov_0'

# scan_combined_Yukawa_profiledTotalXS_onedimTotalXSScan        = 'out/Scan_Yukawa_Nov16_profiledTotalXS_onedimTotalXSScan'
# scan_combined_Yukawa_profiledTotalXS_onedimTotalXSScan_asimov = 'out/Scan_Yukawa_Nov16_profiledTotalXS_onedimTotalXSScan_asimov'

# scan_combined_Yukawa_couplingDependentBR_profiledTotalXS_asimov = 'out/Scan_Yukawa_Feb07_combination_couplingDependentBR_profiledTotalXS'


# # ======================================
# # kappat kappag scans

# # Problems with last theory uncertainty
# # scan_combined_Top        = 'out/Scan_Top_Nov06'
# # scan_hzz_Top             = 'out/Scan_Top_Nov06_0'
# # scan_hgg_Top             = 'out/Scan_Top_Nov06_1'
# # scan_combined_Top_asimov = 'out/Scan_Top_Nov07_asimov'
# # scan_hzz_Top_asimov      = 'out/Scan_Top_Nov07_hzz_asimov'
# # scan_hgg_Top_asimov      = 'out/Scan_Top_Nov07_hgg_asimov'

# # scan_combined_Top                                         = 'out/Scan_Top_Nov08'
# # scan_hgg_Top                                              = 'out/Scan_Top_Nov08_hgg'
# # scan_hzz_Top                                              = 'out/Scan_Top_Nov08_hzz'
# # scan_hgg_Top                                              = 'out/Scan_Top_Nov27_hgg_0'

# # scan_combined_Top                                         = 'out/Scan_Top_Nov27_0'
# # scan_hgg_Top                                              = 'out/Scan_Top_Nov28_hgg_0'
# # scan_hzz_Top                                              = 'out/Scan_Top_Nov26_hzz'

# # scan_hzz_Top                                              = 'out/Scan_Top_Dec03_hzz'

# scan_combined_Top                                         = 'out/Scan_Top_Dec03_2'
# scan_hgg_Top                                              = 'out/Scan_Top_Dec03_hgg_0'
# scan_hzz_Top                                              = 'out/Scan_Top_Dec03_hzz_0'

# # scan_combined_Top_asimov                                  = 'out/Scan_Top_Nov08_asimov'
# # scan_hgg_Top_asimov                                       = 'out/Scan_Top_Nov08_hgg_asimov'
# # scan_hzz_Top_asimov                                       = 'out/Scan_Top_Nov08_hzz_asimov'
# scan_combined_Top_asimov                                  = 'out/Scan_Top_Nov27_asimov_1'
# scan_hgg_Top_asimov                                       = 'out/Scan_Top_Nov26_hgg_asimov'
# scan_hzz_Top_asimov                                       = 'out/Scan_Top_Nov26_hzz_asimov'

# scan_combined_TopHighPt                                   = 'out/Scan_TopHighPt_Dec11'
# scan_hgg_TopHighPt                                        = 'out/Scan_TopHighPt_Dec11_hgg'
# scan_hzz_TopHighPt                                        = 'out/Scan_TopHighPt_Dec11_hzz_0'


# # Implementing reweighting and hbb

# scan_combWithHbb_TopHighPt                                = 'out/Scan_TopHighPt_Jan22_combWithHbb'
# scan_hgg_TopHighPt                                        = 'out/Scan_TopHighPt_Jan22_hgg'
# scan_hzz_TopHighPt                                        = 'out/Scan_TopHighPt_Jan22_hzz'

# # scan_combined_TopHighPt_asimov                            = 'out/Scan_TopHighPt_Dec12_asimov_2'
# # scan_combined_TopHighPt_asimov                            = 'out/Scan_TopHighPt_Jan19_asimov_2'
# # scan_hzz_TopHighPt_asimov                                 = 'out/Scan_TopHighPt_Jan22_hzz_asimov_1'
# # scan_combWithHbb_TopHighPt_asimov                         = 'out/Scan_TopHighPt_Jan19_asimov_0'
# # scan_hzz_TopHighPt_asimov                                 = 'out/Scan_TopHighPt_Jan22_hzz_asimov_0'
# # 
# # scan_combined_TopHighPt_asimov                            = 'out/Scan_TopHighPt_Jan22_asimov_0'
# # scan_combWithHbb_TopHighPt_asimov                         = 'out/Scan_TopHighPt_Jan22_combWithHbb_asimov_0'
# # scan_hgg_TopHighPt_asimov                                 = 'out/Scan_TopHighPt_Jan22_hgg_asimov'
# # scan_hzz_TopHighPt_asimov                                 = 'out/Scan_TopHighPt_Jan22_hzz_asimov_1'
# # 
# # combWithHbb not trustworthy! Last few bins are not processed correctly
# scan_combined_TopHighPt_asimov                            = 'out/Scan_TopHighPt_Jan23_asimov'
# scan_combWithHbb_TopHighPt_asimov                         = 'out/Scan_TopHighPt_Jan23_combWithHbb_asimov'
# scan_hgg_TopHighPt_asimov                                 = 'out/Scan_TopHighPt_Jan23_hgg_asimov'
# scan_hzz_TopHighPt_asimov                                 = 'out/Scan_TopHighPt_Jan23_hzz_asimov'


# # TODO:
# # Scan_TotalXS_Nov27

# scan_combined_Top_extendedRange_asimov                    = 'out/Scan_Top_Nov20_asimov'
# scan_hgg_Top_extendedRange_asimov                         = 'out/Scan_Top_Nov20_hzz_asimov'
# scan_hzz_Top_extendedRange_asimov                         = 'out/Scan_Top_Nov20_hgg_asimov'

# scan_combined_Top_fitOnlyNormalization                    = 'out/Scan_Top_Nov17_profiledTotalXS_fitOnlyNormalization'
# # scan_combined_Top_fitOnlyNormalization_asimov             = 'out/Scan_Top_Nov17_profiledTotalXS_fitOnlyNormalization_asimov'
# scan_combined_Top_bigRange                                = 'out/Scan_Top_Nov17'
# scan_combined_Top_bigRange2                               = 'out/Scan_Top_Nov19'
# scan_combined_Top_skippedLastBin                          = 'out/Scan_Top_Nov17_skippedLastBin'
# scan_combined_Top_skippedLastBin_asimov                   = 'out/Scan_Top_Nov17_skippedLastBin_asimov'
# scan_combined_Top_lumiStudy_asimov                        = 'out/Scan_Top_Nov09_lumiStudy_asimov'
# # scan_combined_Top_profiledTotalXS_asimov                  = 'out/Scan_Top_Nov09_profiledTotalXS_asimov'
# scan_combined_Top_couplingDependentBR_asimov              = 'out/Scan_Top_Nov10_couplingDependentBR_asimov'
# scan_combined_Top_couplingDependentBR_fixedKappaV_asimov  = 'out/Scan_Top_Nov10_couplingDependentBR_fixedKappaV_asimov'
# # scan_combined_Top_couplingDependentBR_kappaVMaxOne_asimov = 'out/Scan_Top_Nov16_couplingDependentBR_kappaVMaxOne_asimov'
# scan_combined_Top_couplingDependentBR_kappaVMaxOne_asimov = 'out/Scan_Top_Nov17_couplingDependentBR_kappaVMaxOne_asimov'
# scan_combined_Top_couplingDependentBR_bigRange_asimov     = 'out/Scan_Top_Nov19_couplingDependentBR_asimov'


# scan_combined_Top_fitOnlyNormalization_asimov             = 'out/Scan_Top_Nov28_fitOnlyNormalization_asimov_0'
# scan_combined_Top_profiledTotalXS_asimov                  = 'out/Scan_Top_Nov27_profiledTotalXS_asimov_0'
# scan_combined_Top_lumiStudy_asimov                        = 'out/Scan_Top_Nov27lumiStudy_asimov_1'


# scan_combined_TopCtCb                                    = 'out/Scan_TopCtCb_Nov15'
# scan_hgg_TopCtCb                                         = 'out/Scan_TopCtCb_Nov15_hgg'
# scan_hzz_TopCtCb                                         = 'out/Scan_TopCtCb_Nov15_hzz'
# scan_combined_TopCtCb_asimov                             = 'out/Scan_TopCtCb_Nov15_asimov'
# scan_hgg_TopCtCb_asimov                                  = 'out/Scan_TopCtCb_Nov15_hgg_asimov'
# scan_hzz_TopCtCb_asimov                                  = 'out/Scan_TopCtCb_Nov15_hzz_asimov'

# # For pre-app
# scan_combined_TopHighPt_couplingDependentBR_profiledTotalXS_asimov = 'out/Scan_TopHighPt_Feb07_combination_couplingDependentBR_profiledTotalXS'


# ########################################
# # Differential combination
# ########################################

# # The up-to-date naming shouldn't be the legacy naming... but fix later
# card_hgg_pth_ggH           = card_hgg_ggHxH_PTH
# card_hzz_pth_ggH           = card_hzz_ggHxH_PTH
# card_combination_pth_ggH   = card_combined_ggHxH_PTH
# card_hgg_pth_smH           = card_hgg_smH_PTH
# card_hzz_pth_smH           = card_hzz_smH_PTH
# card_combination_pth_smH   = card_combined_smH_PTH
# card_hgg_ptjet             = card_hgg_smH_PTJ
# card_hzz_ptjet             = card_hzz_smH_PTJ
# card_combination_ptjet     = card_combined_smH_PTJ
# card_hgg_njets             = card_hgg_smH_NJ
# card_hzz_njets             = card_hzz_smH_NJ
# card_combination_njets     = card_combined_smH_NJ
# card_hgg_rapidity          = card_hgg_smH_YH
# card_hzz_rapidity          = card_hzz_smH_YH
# card_combination_rapidity  = card_combined_smH_YH

# ws_hgg_pth_ggH           = 'out/workspaces_Nov29/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_xHfixed.root'
# ws_hgg_pth_smH           = 'out/workspaces_Nov07/Datacard_13TeV_differential_PtNNLOPS_systs_renamedProcesses.root'
# ws_hgg_ptjet             = 'out/workspaces_Nov30/Datacard_13TeV_differential_Jet2p5Pt0NNLOPS_newBins_renamedProcesses.root'
# ws_hgg_njets             = 'out/workspaces_Nov12/Datacard_13TeV_differential_Njets2p5NNLOPS_renamedProcesses.root'
# ws_hgg_rapidity          = 'out/workspaces_Nov28/Datacard_13TeV_differential_AbsRapidityNNLOPS_newBins_combination_renamedProcesses.root'
# ws_hzz_pth_ggH           = 'out/workspaces_Nov29/hzz4l_comb_13TeV_xs_processesShifted_xHfixed.root'
# ws_hzz_pth_smH           = 'out/workspaces_Nov08/hzz4l_comb_13TeV_xs_processesShifted.root'
# ws_hzz_ptjet             = 'out/workspaces_Nov28/hzz4l_comb_13TeV_xs_ptjet.root'
# ws_hzz_njets             = 'out/workspaces_Nov12/hzz4l_comb_13TeV_xs.root'
# ws_hzz_rapidity          = 'out/workspaces_Nov28/hzz4l_comb_13TeV_xs.root'
# ws_combination_pth_ggH   = 'out/workspaces_Nov08/combinedCard_Nov03_xHfixed.root'
# ws_combination_pth_smH   = 'out/workspaces_Nov07/combinedCard_smH_Nov07.root'
# ws_combination_ptjet     = 'out/workspaces_Nov30/combinedCard_PTJ_smH_Nov28.root'
# ws_combination_njets     = 'out/workspaces_Nov12/combinedCard_NJ_smH_Nov12.root'
# ws_combination_rapidity  = 'out/workspaces_Nov28/combinedCard_YH_smH_Nov28.root'

# ws_combination_pth_ggH_lumiScale   = 'out/workspaces_Feb12/combinedCard_Nov03_DifferentialModel_lumiScale.root'
# ws_combination_pth_smH_lumiScale   = 'out/workspaces_Feb12/combinedCard_smH_Nov07_DifferentialModel_lumiScale.root'
# ws_combination_ptjet_lumiScale     = 'out/workspaces_Feb12/combinedCard_PTJ_smH_Nov28_DifferentialModel_lumiScale.root'
# ws_combination_njets_lumiScale     = 'out/workspaces_Feb12/combinedCard_NJ_smH_Nov12_DifferentialModel_lumiScale.root'
# ws_combination_rapidity_lumiScale  = 'out/workspaces_Feb12/combinedCard_YH_smH_Nov28_DifferentialModel_lumiScale.root'


# # Legacy in order not to break old code
# ws_hgg_ggH_xHfixed       = ws_hgg_pth_ggH
# ws_hgg_smH               = ws_hgg_pth_smH
# ws_hgg_smH_PTJ           = ws_hgg_ptjet
# ws_hgg_smH_NJ            = ws_hgg_njets
# ws_hgg_smH_YH            = ws_hgg_rapidity
# ws_hzz_ggH_xHfixed       = ws_hzz_pth_ggH
# ws_hzz_smH               = ws_hzz_pth_smH
# ws_hzz_smH_PTJ           = ws_hzz_ptjet
# ws_hzz_smH_NJ            = ws_hzz_njets
# ws_hzz_smH_YH            = ws_hzz_rapidity
# ws_combined_ggH_xHfixed  = ws_combination_pth_ggH
# ws_combined_smH          = ws_combination_pth_smH
# ws_combined_smH_PTJ      = ws_combination_ptjet
# ws_combined_smH_NJ       = ws_combination_njets
# ws_combined_smH_YH       = ws_combination_rapidity

# postfit_hgg_pth_ggH                 = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_xHfixed_pth_ggH.MultiDimFit.mH125.root'
# postfit_hgg_pth_smH                 = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_Datacard_13TeV_differential_PtNNLOPS_systs_renamedProcesses_pth_smH.MultiDimFit.mH125.root'
# postfit_hgg_ptjet                   = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_Datacard_13TeV_differential_Jet2p5Pt0NNLOPS_newBins_renamedProcesses_ptjet.MultiDimFit.mH125.root'
# postfit_hgg_njets                   = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_Datacard_13TeV_differential_Njets2p5NNLOPS_renamedProcesses_njets.MultiDimFit.mH125.root'
# postfit_hgg_rapidity                = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_Datacard_13TeV_differential_AbsRapidityNNLOPS_newBins_combination_renamedProcesses_rapidity.MultiDimFit.mH125.root'
# postfit_hgg_pth_ggH_asimov          = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_xHfixed_pth_ggH_asimov.MultiDimFit.mH125.root'
# postfit_hgg_pth_smH_asimov          = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_Datacard_13TeV_differential_PtNNLOPS_systs_renamedProcesses_pth_smH_asimov.MultiDimFit.mH125.root'
# postfit_hgg_ptjet_asimov            = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_Datacard_13TeV_differential_Jet2p5Pt0NNLOPS_newBins_renamedProcesses_ptjet_asimov.MultiDimFit.mH125.root'
# postfit_hgg_njets_asimov            = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_Datacard_13TeV_differential_Njets2p5NNLOPS_renamedProcesses_njets_asimov.MultiDimFit.mH125.root'
# postfit_hgg_rapidity_asimov         = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_Datacard_13TeV_differential_AbsRapidityNNLOPS_newBins_combination_renamedProcesses_rapidity_asimov.MultiDimFit.mH125.root'
# postfit_hzz_pth_ggH                 = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_hzz4l_comb_13TeV_xs_processesShifted_xHfixed_pth_ggH.MultiDimFit.mH125.root'
# postfit_hzz_pth_smH                 = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_hzz4l_comb_13TeV_xs_processesShifted_pth_smH.MultiDimFit.mH125.root'
# postfit_hzz_ptjet                   = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_hzz4l_comb_13TeV_xs_ptjet_ptjet.MultiDimFit.mH125.root'
# postfit_hzz_njets                   = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_hzz4l_comb_13TeV_xs_njets.MultiDimFit.mH125.root'
# postfit_hzz_rapidity                = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_hzz4l_comb_13TeV_xs_rapidity.MultiDimFit.mH125.root'
# postfit_hzz_pth_ggH_asimov          = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_hzz4l_comb_13TeV_xs_processesShifted_xHfixed_pth_ggH_asimov.MultiDimFit.mH125.root'
# postfit_hzz_pth_smH_asimov          = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_hzz4l_comb_13TeV_xs_processesShifted_pth_smH_asimov.MultiDimFit.mH125.root'
# postfit_hzz_ptjet_asimov            = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_hzz4l_comb_13TeV_xs_ptjet_ptjet_asimov.MultiDimFit.mH125.root'
# postfit_hzz_njets_asimov            = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_hzz4l_comb_13TeV_xs_njets_asimov.MultiDimFit.mH125.root'
# postfit_hzz_rapidity_asimov         = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_hzz4l_comb_13TeV_xs_rapidity_asimov.MultiDimFit.mH125.root'
# postfit_combination_pth_ggH         = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_combinedCard_Nov03_xHfixed_pth_ggH.MultiDimFit.mH125.root'
# postfit_combination_pth_smH         = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_combinedCard_smH_Nov07_pth_smH.MultiDimFit.mH125.root'
# postfit_combination_ptjet           = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_combinedCard_PTJ_smH_Nov28_ptjet.MultiDimFit.mH125.root'
# postfit_combination_njets           = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_combinedCard_NJ_smH_Nov12_njets.MultiDimFit.mH125.root'
# postfit_combination_rapidity        = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_combinedCard_YH_smH_Nov28_rapidity.MultiDimFit.mH125.root'
# postfit_combination_pth_ggH_asimov  = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_combinedCard_Nov03_xHfixed_pth_ggH_asimov.MultiDimFit.mH125.root'
# postfit_combination_pth_smH_asimov  = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_combinedCard_smH_Nov07_pth_smH_asimov.MultiDimFit.mH125.root'
# postfit_combination_ptjet_asimov    = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_combinedCard_PTJ_smH_Nov28_ptjet_asimov.MultiDimFit.mH125.root'
# postfit_combination_njets_asimov    = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_combinedCard_NJ_smH_Nov12_njets_asimov.MultiDimFit.mH125.root'
# postfit_combination_rapidity_asimov = 'out/postfitWss_Feb05/higgsCombine_POSTFIT_combinedCard_YH_smH_Nov28_rapidity_asimov.MultiDimFit.mH125.root'

# # scan_combination_njets_statonly     = 'out/Scan_njets_Feb06_combination_statonly'
# # scan_hgg_njets_statonly             = 'out/Scan_njets_Feb06_hgg_statonly'
# # scan_hzz_njets_statonly             = 'out/Scan_njets_Feb06_hzz_statonly'
# # scan_combination_pth_ggH_statonly   = 'out/Scan_pth_ggH_Feb06_combination_statonly'
# # scan_hgg_pth_ggH_statonly           = 'out/Scan_pth_ggH_Feb06_hgg_statonly'
# # scan_hzz_pth_ggH_statonly           = 'out/Scan_pth_ggH_Feb06_hzz_statonly'
# # scan_combination_pth_smH_statonly   = 'out/Scan_pth_smH_Feb06_combination_statonly'
# # scan_hgg_pth_smH_statonly           = 'out/Scan_pth_smH_Feb06_hgg_statonly'
# # scan_hzz_pth_smH_statonly           = 'out/Scan_pth_smH_Feb06_hzz_statonly'
# # scan_combination_ptjet_statonly     = 'out/Scan_ptjet_Feb06_combination_statonly'
# # scan_hgg_ptjet_statonly             = 'out/Scan_ptjet_Feb06_hgg_statonly'
# # scan_hzz_ptjet_statonly             = 'out/Scan_ptjet_Feb06_hzz_statonly'
# # scan_combination_rapidity_statonly  = 'out/Scan_rapidity_Feb06_combination_statonly'
# # scan_hgg_rapidity_statonly          = 'out/Scan_rapidity_Feb06_hgg_statonly'
# # scan_hzz_rapidity_statonly          = 'out/Scan_rapidity_Feb06_hzz_statonly'

# # scan_combination_njets    = 'out/Scan_nJets_Nov12'
# # scan_hgg_njets            = 'out/Scan_nJets_Nov12_0'
# # scan_hzz_njets            = 'out/Scan_nJets_Nov12_1'
# # scan_combination_pth_ggH  = 'out/Scan_PTH_Nov08_xHfixed'
# # scan_hgg_pth_ggH          = 'out/Scan_PTH_Nov29_xHfixed_hgg'
# # scan_hzz_pth_ggH          = 'out/Scan_PTH_Nov29_xHfixed_hzz'
# # scan_combination_pth_smH  = 'out/Scan_PTH_Nov07'
# # scan_hgg_pth_smH          = 'out/Scan_PTH_Nov07_hgg'
# # scan_hzz_pth_smH          = 'out/Scan_PTH_Nov08_hzz'
# # scan_combination_ptjet    = 'out/Scan_PTJ_Nov28_1'
# # scan_hgg_ptjet            = 'out/Scan_PTJ_Nov28_0'
# # scan_hzz_ptjet            = 'out/Scan_PTJ_Nov28'
# # scan_combination_rapidity = 'out/Scan_YH_Nov28'
# # scan_hgg_rapidity         = 'out/Scan_YH_Nov28_0'
# # scan_hzz_rapidity         = 'out/Scan_YH_Nov28_1'

# scan_combination_njets                    = 'out/Scan_njets_Feb06_combination'
# scan_combination_njets_asimov             = 'out/Scan_njets_Feb06_combination_asimov'
# scan_combination_njets_statonly           = 'out/Scan_njets_Feb06_combination_statonly'
# scan_combination_njets_statonly_asimov    = 'out/Scan_njets_Feb06_combination_statonly_asimov'
# scan_hgg_njets                            = 'out/Scan_njets_Feb06_hgg'
# scan_hgg_njets_asimov                     = 'out/Scan_njets_Feb06_hgg_asimov'
# scan_hgg_njets_statonly                   = 'out/Scan_njets_Feb06_hgg_statonly'
# scan_hgg_njets_statonly_asimov            = 'out/Scan_njets_Feb06_hgg_statonly_asimov'
# scan_hzz_njets                            = 'out/Scan_njets_Feb06_hzz'
# scan_hzz_njets_asimov                     = 'out/Scan_njets_Feb06_hzz_asimov'
# scan_hzz_njets_statonly                   = 'out/Scan_njets_Feb06_hzz_statonly'
# scan_hzz_njets_statonly_asimov            = 'out/Scan_njets_Feb06_hzz_statonly_asimov'
# scan_combination_pth_ggH                  = 'out/Scan_pth_ggH_Feb06_combination'
# scan_combination_pth_ggH_asimov           = 'out/Scan_pth_ggH_Feb06_combination_asimov'
# scan_combination_pth_ggH_statonly         = 'out/Scan_pth_ggH_Feb06_combination_statonly'
# scan_combination_pth_ggH_statonly_asimov  = 'out/Scan_pth_ggH_Feb06_combination_statonly_asimov'
# scan_hgg_pth_ggH                          = 'out/Scan_pth_ggH_Feb06_hgg'
# scan_hgg_pth_ggH_asimov                   = 'out/Scan_pth_ggH_Feb06_hgg_asimov'
# scan_hgg_pth_ggH_statonly                 = 'out/Scan_pth_ggH_Feb06_hgg_statonly'
# scan_hgg_pth_ggH_statonly_asimov          = 'out/Scan_pth_ggH_Feb06_hgg_statonly_asimov'
# scan_hzz_pth_ggH                          = 'out/Scan_pth_ggH_Feb06_hzz'
# scan_hzz_pth_ggH_asimov                   = 'out/Scan_pth_ggH_Feb06_hzz_asimov'
# scan_hzz_pth_ggH_statonly                 = 'out/Scan_pth_ggH_Feb06_hzz_statonly'
# scan_hzz_pth_ggH_statonly_asimov          = 'out/Scan_pth_ggH_Feb06_hzz_statonly_asimov'
# scan_combination_pth_smH                  = 'out/Scan_pth_smH_Feb06_combination'
# scan_combination_pth_smH_asimov           = 'out/Scan_pth_smH_Feb06_combination_asimov'
# scan_combination_pth_smH_statonly         = 'out/Scan_pth_smH_Feb06_combination_statonly'
# scan_combination_pth_smH_statonly_asimov  = 'out/Scan_pth_smH_Feb06_combination_statonly_asimov'
# scan_hgg_pth_smH                          = 'out/Scan_pth_smH_Feb06_hgg'
# scan_hgg_pth_smH_asimov                   = 'out/Scan_pth_smH_Feb06_hgg_asimov'
# scan_hgg_pth_smH_statonly                 = 'out/Scan_pth_smH_Feb06_hgg_statonly'
# scan_hgg_pth_smH_statonly_asimov          = 'out/Scan_pth_smH_Feb06_hgg_statonly_asimov'
# scan_hzz_pth_smH                          = 'out/Scan_pth_smH_Feb06_hzz'
# scan_hzz_pth_smH_asimov                   = 'out/Scan_pth_smH_Feb06_hzz_asimov'
# scan_hzz_pth_smH_statonly                 = 'out/Scan_pth_smH_Feb06_hzz_statonly'
# scan_hzz_pth_smH_statonly_asimov          = 'out/Scan_pth_smH_Feb06_hzz_statonly_asimov'
# scan_combination_ptjet                    = 'out/Scan_ptjet_Feb06_combination'
# scan_combination_ptjet_asimov             = 'out/Scan_ptjet_Feb06_combination_asimov'
# scan_combination_ptjet_statonly           = 'out/Scan_ptjet_Feb06_combination_statonly'
# scan_combination_ptjet_statonly_asimov    = 'out/Scan_ptjet_Feb06_combination_statonly_asimov'
# scan_hgg_ptjet                            = 'out/Scan_ptjet_Feb06_hgg'
# scan_hgg_ptjet_asimov                     = 'out/Scan_ptjet_Feb06_hgg_asimov'
# scan_hgg_ptjet_statonly                   = 'out/Scan_ptjet_Feb06_hgg_statonly'
# scan_hgg_ptjet_statonly_asimov            = 'out/Scan_ptjet_Feb06_hgg_statonly_asimov'
# scan_hzz_ptjet                            = 'out/Scan_ptjet_Feb06_hzz'
# scan_hzz_ptjet_asimov                     = 'out/Scan_ptjet_Feb06_hzz_asimov'
# scan_hzz_ptjet_statonly                   = 'out/Scan_ptjet_Feb06_hzz_statonly'
# scan_hzz_ptjet_statonly_asimov            = 'out/Scan_ptjet_Feb06_hzz_statonly_asimov'
# scan_combination_rapidity                 = 'out/Scan_rapidity_Feb06_combination'
# scan_combination_rapidity_asimov          = 'out/Scan_rapidity_Feb06_combination_asimov'
# scan_combination_rapidity_statonly        = 'out/Scan_rapidity_Feb06_combination_statonly'
# scan_combination_rapidity_statonly_asimov = 'out/Scan_rapidity_Feb06_combination_statonly_asimov'
# scan_hgg_rapidity                         = 'out/Scan_rapidity_Feb06_hgg'
# scan_hgg_rapidity_asimov                  = 'out/Scan_rapidity_Feb06_hgg_asimov'
# scan_hgg_rapidity_statonly                = 'out/Scan_rapidity_Feb06_hgg_statonly'
# scan_hgg_rapidity_statonly_asimov         = 'out/Scan_rapidity_Feb06_hgg_statonly_asimov'
# scan_hzz_rapidity                         = 'out/Scan_rapidity_Feb06_hzz'
# scan_hzz_rapidity_asimov                  = 'out/Scan_rapidity_Feb06_hzz_asimov'
# scan_hzz_rapidity_statonly                = 'out/Scan_rapidity_Feb06_hzz_statonly'
# scan_hzz_rapidity_statonly_asimov         = 'out/Scan_rapidity_Feb06_hzz_statonly_asimov'

# scan_combination_njets_lumiScale_asimov    = 'out/Scan_njets_Feb12_combination_lumiScale_asimov'
# scan_combination_pth_ggH_lumiScale_asimov  = 'out/Scan_pth_ggH_Feb12_combination_lumiScale_asimov'
# scan_combination_pth_smH_lumiScale_asimov  = 'out/Scan_pth_smH_Feb12_combination_lumiScale_asimov'
# scan_combination_ptjet_lumiScale_asimov    = 'out/Scan_ptjet_Feb12_combination_lumiScale_asimov'
# scan_combination_rapidity_lumiScale_asimov = 'out/Scan_rapidity_Feb12_combination_lumiScale_asimov'

# # Legacy keys, so old code doesn't break
# scan_combined_NJ      = scan_combination_njets
# scan_hgg_NJ           = scan_hgg_njets
# scan_hzz_NJ           = scan_hzz_njets
# scan_combined_PTH_ggH = scan_combination_pth_ggH
# scan_hgg_PTH_ggH      = scan_hgg_pth_ggH
# scan_hzz_PTH_ggH      = scan_hzz_pth_ggH
# scan_combined_PTH     = scan_combination_pth_smH
# scan_hgg_PTH          = scan_hgg_pth_smH
# scan_hzz_PTH          = scan_hzz_pth_smH
# scan_combined_PTJ     = scan_combination_ptjet
# scan_hgg_PTJ          = scan_hgg_ptjet
# scan_hzz_PTJ          = scan_hzz_ptjet
# scan_combined_YH      = scan_combination_rapidity
# scan_hgg_YH           = scan_hgg_rapidity
# scan_hzz_YH           = scan_hzz_rapidity

# #____________________________________________________________________
# # To implement
# scan_hbb_PTH_ggH         = 'out/Scan_PTH_Jan23_xHfixed_hbb'
# scan_combWithHbb_PTH_ggH = 'out/Scan_PTH_Jan23_xHfixed_combWithHbb_2'

# # Asimovs
# scan_combined_PTH_ggH_asimov    = 'out/Scan_PTH_Dec15_xHfixed_asimov_asimov'
# scan_hgg_PTH_ggH_asimov         = 'out/Scan_PTH_Dec15_xHfixed_hgg_asimov_asimov'
# scan_hzz_PTH_ggH_asimov         = 'out/Scan_PTH_Dec15_xHfixed_hzz_asimov_asimov'
# scan_hbb_PTH_ggH_asimov         = 'out/Scan_PTH_Jan23_xHfixed_hbb_asimov'
# scan_combWithHbb_PTH_ggH_asimov = 'out/Scan_PTH_Jan23_xHfixed_combWithHbb_asimov'

# scan_combined_YH_asimov  = 'out/Scan_YH_Nov28_asimov'
# scan_hgg_YH_asimov       = 'out/Scan_YH_Nov28_asimov_0'
# scan_hzz_YH_asimov       = 'out/Scan_YH_Nov28_asimov_1'

# scan_combined_PTJ_asimov = 'out/Scan_PTJ_Nov30_asimov'
# scan_hgg_PTJ_asimov      = 'out/Scan_PTJ_Nov30_asimov_0'
# scan_hzz_PTJ_asimov      = 'out/Scan_PTJ_Nov28_asimov_1'


# ########################################
# # End
# ########################################
# if __name__ == "__main__":
#     vardict = vars()
#     for key in vardict.keys():
#         if key.startswith('__'): continue
#         print '{0}: {1}'.format( key, vardict[key] )
