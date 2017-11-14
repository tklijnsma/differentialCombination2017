########################################
# Variables
########################################

YR4_totalXS = 55.70628722 # pb


########################################
# Text datacards
########################################

# ----- PTH -----
card_hgg_smH_PTH_unprocessed   = 'suppliedInput/fromVittorio/pT_NNLOPS_Nov01/Datacard_13TeV_differential_PtNNLOPS_systs.txt'
card_hzz_smH_PTH_unprocessed   = 'suppliedInput/fromDavid/PTH_NNLOPS_Nov01/smH/hzz4l_comb_13TeV_xs.txt'

card_hgg_ggHxH_PTH_unprocessed = 'suppliedInput/fromVittorio/pT_NNLOPS_ggHxH_Nov03/Datacard_13TeV_differential_PtGghPlusHxNNLOPS.txt'
card_hzz_ggHxH_PTH_unprocessed = 'suppliedInput/fromDavid/PTH_NNLOPS_Nov01/ggH/hzz4l_comb_13TeV_xs.txt'

card_hgg_ggHxH_PTH             = 'suppliedInput/fromVittorio/pT_NNLOPS_ggHxH_Nov03/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses.txt'
card_hzz_ggHxH_PTH             = 'suppliedInput/fromDavid/PTH_NNLOPS_Nov01/ggH/hzz4l_comb_13TeV_xs_processesShifted.txt'

card_combined_ggHxH_PTH        = 'suppliedInput/combinedCard_Nov03.txt'

card_hgg_smH_PTH               = 'suppliedInput/fromVittorio/pT_NNLOPS_Nov01/Datacard_13TeV_differential_PtNNLOPS_systs_renamedProcesses.txt'
card_hzz_smH_PTH               = 'suppliedInput/fromDavid/PTH_NNLOPS_Nov01/smH/hzz4l_comb_13TeV_xs_processesShifted.txt'
card_combined_smH_PTH          = 'suppliedInput/combinedCard_smH_Nov07.txt'


# ----- NJ -----
card_hgg_smH_NJ_unprocessed    = 'suppliedInput/fromVittorio/differential_Njets2p5NNLOPS_Nov10/Datacard_13TeV_differential_Njets2p5NNLOPS.txt'

card_hzz_smH_NJ                = 'suppliedInput/fromDavid/NJ_NNLOPS_Nov01/smH/hzz4l_comb_13TeV_xs.txt'
card_hgg_smH_NJ                = 'suppliedInput/fromVittorio/differential_Njets2p5NNLOPS_Nov10/Datacard_13TeV_differential_Njets2p5NNLOPS_renamedProcesses.txt'

# card_combined_smH_NJ           = 'suppliedInput/combinedCard_NJ_smH_Nov10.txt' # used JER as scale instead of JEC
card_combined_smH_NJ           = 'suppliedInput/combinedCard_NJ_smH_Nov12.txt'


# ----- YH -----

card_hgg_smH_YH_unprocessed    = 'suppliedInput/fromVittorio/differential_AbsRapidityNNLOPS_newBins_Nov11/Datacard_13TeV_differential_AbsRapidityNNLOPS_newBins.txt'
card_hgg_smH_YH                = 'suppliedInput/fromVittorio/differential_AbsRapidityNNLOPS_newBins_Nov11/Datacard_13TeV_differential_AbsRapidityNNLOPS_newBins_renamedProcesses.txt'

# No OutsideAcceptance in hzz, so no renumbering necessary
card_hzz_smH_YH                = 'suppliedInput/fromDavid/YH_NNLOPS_Nov12/smH/hzz4l_comb_13TeV_xs.txt'

card_combined_smH_YH           = 'suppliedInput/combinedCard_YH_smH_Nov12.txt'


########################################
# Derived theory files
########################################

derivedTheoryFiles_YukawaQuarkInduced       = 'derivedTheoryFiles_Nov03_YukawaQuarkInduced'
derivedTheoryFiles_YukawaGluonInduced       = 'derivedTheoryFiles_Nov03_YukawaGluonInduced'
derivedTheoryFiles_YukawaQuarkInducedScaled = 'derivedTheoryFiles_Nov03_YukawaQuarkInducedScaled'
derivedTheoryFiles_YukawaSummed             = 'derivedTheoryFiles_Nov03_YukawaSummed'
derivedTheoryFiles_Top                      = 'derivedTheoryFiles_Nov06_Top'


########################################
# Correlation matrices and uncertainties
########################################

correlationMatrix_Yukawa              = 'correlationMatrices_Nov03/corrMat_exp.txt'
theoryUncertainties_Yukawa            = 'correlationMatrices_Nov03/errors_for_corrMat_exp.txt'
correlationMatrix_Yukawa_Uncorrelated = 'correlationMatrices_Nov03/corrMat_exp_UNCORRELATED.txt'

correlationMatrix_Top                 = 'correlationMatrices_Nov06_Top/corrMat_exp.txt'
theoryUncertainties_Top               = 'correlationMatrices_Nov06_Top/errors_for_corrMat_exp.txt'


########################################
# Workspaces after running t2ws
########################################

# ======================================
# Unsplit workspaces for non-coupling combinations

ws_hgg_smH         = 'workspaces_Nov07/Datacard_13TeV_differential_PtNNLOPS_systs_renamedProcesses.root'
ws_hzz_smH         = 'workspaces_Nov08/hzz4l_comb_13TeV_xs_processesShifted.root'
ws_combined_smH    = 'workspaces_Nov07/combinedCard_smH_Nov07.root'

ws_hgg_smH_NJ      = 'workspaces_Nov12/Datacard_13TeV_differential_Njets2p5NNLOPS_renamedProcesses.root'
ws_hzz_smH_NJ      = 'workspaces_Nov12/hzz4l_comb_13TeV_xs.root'
ws_combined_smH_NJ = 'workspaces_Nov12/combinedCard_NJ_smH_Nov12.root'

# ---------------------
# Workspaces for extra studies

ws_ratio_of_BRs         = 'workspaces_Nov08/combinedCard_smH_Nov07_FitBRModel.root'
ws_totalXS              = 'workspaces_Nov08/combinedCard_smH_Nov07_FitBRModel_fitTotalXS.root'
ws_combined_ggH_xHfixed = 'workspaces_Nov08/combinedCard_Nov03_xHfixed.root'


# ======================================
# kappab kappac

ws_combined_Yukawa = 'workspaces_Nov03/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties.root'
ws_hgg_Yukawa      = 'workspaces_Nov03/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_Yukawa_withTheoryUncertainties.root'
ws_hzz_Yukawa      = 'workspaces_Nov03/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_Yukawa_withTheoryUncertainties.root'

ws_combined_Yukawa_noTheoryUncertainties               = 'workspaces_Nov06/combinedCard_Nov03_CouplingModel_Yukawa_noTheoryUncertainties.root'
ws_combined_Yukawa_withUncorrelatedTheoryUncertainties = 'workspaces_Nov06/combinedCard_Nov03_CouplingModel_Yukawa_withUncorrelatedTheoryUncertainties.root'
ws_combined_Yukawa_lumiScalable                        = 'workspaces_Nov06/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_lumiScale.root'
ws_combined_Yukawa_profiledTotalXS                     = 'workspaces_Nov09/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_profiledTotalXS.root'
ws_combined_Yukawa_couplingDependentBR                 = 'workspaces_Nov10/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_couplingDependentBR.root'


# ======================================
# kappat kappag

# Problem with last theory uncertainty:
# ws_combined_Top    = 'workspaces_Nov06/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties.root'
# ws_hgg_Top         = 'workspaces_Nov06/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_Top_withTheoryUncertainties.root'
# ws_hzz_Top         = 'workspaces_Nov06/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_Top_withTheoryUncertainties.root'

# Last theory uncertainty skipped:
ws_hgg_Top         = 'workspaces_Nov08/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_Top_withTheoryUncertainties.root'
ws_combined_Top    = 'workspaces_Nov08/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties.root'
ws_hzz_Top         = 'workspaces_Nov08/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_Top_withTheoryUncertainties.root'

ws_combined_Top_lumiScalable           = 'workspaces_Nov09/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_lumiScale.root'
ws_combined_Top_profiledTotalXS        = 'workspaces_Nov09/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_profiledTotalXS.root'
ws_combined_Top_couplingDependentBR    = 'workspaces_Nov10/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_couplingDependentBR.root'


########################################
# Scans
########################################

# ======================================
# Non-coupling combinations

scan_combined_PTH = 'Scan_PTH_Nov07'
scan_hgg_PTH      = 'Scan_PTH_Nov07_hgg'
scan_hzz_PTH      = 'Scan_PTH_Nov08_hzz'

scan_combined_NJ  = 'Scan_nJets_Nov12'
scan_hgg_NJ       = 'Scan_nJets_Nov12_0'
scan_hzz_NJ       = 'Scan_nJets_Nov12_1'


# ======================================
# Extra studies

scan_ratioOfBRs           = 'Scan_ratioOfBRs_Nov08'
scan_combined_totalXS     = 'Scan_TotalXS_Nov08'
scan_combined_PTH_xHfixed = 'Scan_PTH_Nov08_xHfixed'

# ======================================
# kappab kappac scans

scan_combined_Yukawa                                        = 'Scan_Yukawa_Nov03_0'
scan_hgg_Yukawa                                             = 'Scan_Yukawa_Nov03_hgg'
scan_hzz_Yukawa                                             = 'Scan_Yukawa_Nov03_hzz'
scan_combined_Yukawa_asimov                                 = 'Scan_Yukawa_Nov03_asimov'
scan_hgg_Yukawa_asimov                                      = 'Scan_Yukawa_Nov03_hgg_asimov'
scan_hzz_Yukawa_asimov                                      = 'Scan_Yukawa_Nov03_hzz_asimov'

scan_combined_Yukawa_lumiStudy_asimov                       = 'Scan_Yukawa_Nov06_lumiStudy_asimov'
scan_combined_Yukawa_noTheoryUncertainties_asimov           = 'Scan_Yukawa_Nov06_noTheoryUncertainties_asimov'
scan_combined_Yukawa_uncorrelatedTheoryUncertainties_asimov = 'Scan_Yukawa_Nov06_uncorrelatedTheoryUncertainties_asimov'
scan_combined_Yukawa_profiledTotalXS_asimov                 = 'Scan_Yukawa_Nov09_profiledTotalXS_asimov'

scan_combined_Yukawa_oneKappa_kappac                        = 'Scan_Yukawa_Nov06_oneKappa_kappac'
scan_combined_Yukawa_oneKappa_kappac_asimov                 = 'Scan_Yukawa_Nov06_oneKappa_kappac_asimov'
scan_combined_Yukawa_oneKappa_kappab                        = 'Scan_Yukawa_Nov06_oneKappa_kappab'
scan_combined_Yukawa_oneKappa_kappab_asimov                 = 'Scan_Yukawa_Nov06_oneKappa_kappab_asimov'

scan_combined_Yukawa_couplingDependentBR_asimov             = 'Scan_Yukawa_Nov10_couplingDependentBR_asimov'
scan_combined_Yukawa_couplingDependentBR_fixedKappaV_asimov = 'Scan_Yukawa_Nov10_couplingDependentBR_fixedKappaV_asimov'
 

# ======================================
# kappat kappag scans

# Problems with last theory uncertainty
# scan_combined_Top        = 'Scan_Top_Nov06'
# scan_hzz_Top             = 'Scan_Top_Nov06_0'
# scan_hgg_Top             = 'Scan_Top_Nov06_1'
# scan_combined_Top_asimov = 'Scan_Top_Nov07_asimov'
# scan_hzz_Top_asimov      = 'Scan_Top_Nov07_hzz_asimov'
# scan_hgg_Top_asimov      = 'Scan_Top_Nov07_hgg_asimov'

scan_combined_Top                                        = 'Scan_Top_Nov08'
scan_hgg_Top                                             = 'Scan_Top_Nov08_hgg'
scan_hzz_Top                                             = 'Scan_Top_Nov08_hzz'
scan_combined_Top_asimov                                 = 'Scan_Top_Nov08_asimov'
scan_hgg_Top_asimov                                      = 'Scan_Top_Nov08_hgg_asimov'
scan_hzz_Top_asimov                                      = 'Scan_Top_Nov08_hzz_asimov'

scan_combined_Top_profiledTotalXS_asimov                 = 'Scan_Top_Nov09_profiledTotalXS_asimov'
scan_combined_Top_couplingDependentBR_asimov             = 'Scan_Top_Nov10_couplingDependentBR_asimov'
scan_combined_Top_couplingDependentBR_fixedKappaV_asimov = 'Scan_Top_Nov10_couplingDependentBR_fixedKappaV_asimov'



########################################
# End
########################################
if __name__ == "__main__":
    keys = [ k for k in vars().keys() ]
    for key in keys:
        if key.startswith('__'): continue
        print '{0}: {1}'.format( key, vars()[key] )