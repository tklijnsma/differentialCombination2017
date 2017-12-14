########################################
# Text datacards
########################################

# ----- PTH -----
card_hgg_smH_PTH_unprocessed   = 'out/suppliedInput/fromVittorio/pT_NNLOPS_Nov01/Datacard_13TeV_differential_PtNNLOPS_systs.txt'
card_hzz_smH_PTH_unprocessed   = 'out/suppliedInput/fromDavid/PTH_NNLOPS_Nov01/smH/hzz4l_comb_13TeV_xs.txt'

card_hgg_ggHxH_PTH_unprocessed = 'out/suppliedInput/fromVittorio/pT_NNLOPS_ggHxH_Nov03/Datacard_13TeV_differential_PtGghPlusHxNNLOPS.txt'
card_hzz_ggHxH_PTH_unprocessed = 'out/suppliedInput/fromDavid/PTH_NNLOPS_Nov01/ggH/hzz4l_comb_13TeV_xs.txt'

card_hgg_ggHxH_PTH             = 'out/suppliedInput/fromVittorio/pT_NNLOPS_ggHxH_Nov03/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses.txt'
card_hzz_ggHxH_PTH             = 'out/suppliedInput/fromDavid/PTH_NNLOPS_Nov01/ggH/hzz4l_comb_13TeV_xs_processesShifted.txt'
card_combined_ggHxH_PTH        = 'out/suppliedInput/combinedCard_Nov03.txt'

card_hgg_smH_PTH               = 'out/suppliedInput/fromVittorio/pT_NNLOPS_Nov01/Datacard_13TeV_differential_PtNNLOPS_systs_renamedProcesses.txt'
card_hzz_smH_PTH               = 'out/suppliedInput/fromDavid/PTH_NNLOPS_Nov01/smH/hzz4l_comb_13TeV_xs_processesShifted.txt'
card_combined_smH_PTH          = 'out/suppliedInput/combinedCard_smH_Nov07.txt'

# ----- NJ -----
card_hgg_smH_NJ_unprocessed    = 'out/suppliedInput/fromVittorio/differential_Njets2p5NNLOPS_Nov10/Datacard_13TeV_differential_Njets2p5NNLOPS.txt'
card_hzz_smH_NJ                = 'out/suppliedInput/fromDavid/NJ_NNLOPS_Nov01/smH/hzz4l_comb_13TeV_xs.txt'
card_hgg_smH_NJ                = 'out/suppliedInput/fromVittorio/differential_Njets2p5NNLOPS_Nov10/Datacard_13TeV_differential_Njets2p5NNLOPS_renamedProcesses.txt'
card_combined_smH_NJ           = 'out/suppliedInput/combinedCard_NJ_smH_Nov12.txt'

# ----- YH -----
card_hgg_smH_YH_unprocessed      = 'out/suppliedInput/fromVittorio/differential_AbsRapidityNNLOPS_newBins_Nov28/Datacard_13TeV_differential_AbsRapidityNNLOPS_newBins_combination.txt'
card_hgg_smH_YH                  = 'out/suppliedInput/fromVittorio/differential_AbsRapidityNNLOPS_newBins_Nov28/Datacard_13TeV_differential_AbsRapidityNNLOPS_newBins_combination_renamedProcesses.txt'
card_hzz_smH_YH                = 'out/suppliedInput/fromDavid/YH_NNLOPS_Nov12/smH/hzz4l_comb_13TeV_xs.txt' # No OutsideAcceptance in hzz, so no renumbering necessary
card_combined_smH_YH           = 'out/suppliedInput/combinedCard_YH_smH_Nov28.txt'

# ----- PTJ -----
card_hzz_smH_PTJ               = 'out/suppliedInput/fromDavid/PTJET_NNLOPS_Nov28/smH/hzz4l_comb_13TeV_xs.txt'
card_hgg_smH_PTJ_unprocessed   = 'out/suppliedInput/fromVittorio/differential_Jet2p5Pt0NNLOPS_newBins_Nov28/Datacard_13TeV_differential_Jet2p5Pt0NNLOPS_newBins.txt'
card_hgg_smH_PTJ               = 'out/suppliedInput/fromVittorio/differential_Jet2p5Pt0NNLOPS_newBins_Nov28/Datacard_13TeV_differential_Jet2p5Pt0NNLOPS_newBins_renamedProcesses.txt'
card_combined_smH_PTJ          = 'out/suppliedInput/combinedCard_PTJ_smH_Nov28.txt'

# ----- INCLUSIVE -----
card_hgg_INC_unprocessed       = 'out/suppliedInput/fromVittorio/inclusive_Nov27/Datacard_13TeV_differential_InclusiveNNLOPS.txt'
card_hzz_INC_unprocessed       = 'out/suppliedInput/fromDavid/differential_Nov27/smH/hzz4l_comb_13TeV_xs.txt'
card_combined_INC              = 'out/suppliedInput/combinedCard_smH_Nov27_INCLUSIVE.txt'


########################################
# Derived theory files
########################################

derivedTheoryFiles_YukawaQuarkInduced       = 'out/derivedTheoryFiles_Nov03_YukawaQuarkInduced'
derivedTheoryFiles_YukawaGluonInduced       = 'out/derivedTheoryFiles_Nov03_YukawaGluonInduced'
derivedTheoryFiles_YukawaQuarkInducedScaled = 'out/derivedTheoryFiles_Nov03_YukawaQuarkInducedScaled'
derivedTheoryFiles_YukawaSummed             = 'out/derivedTheoryFiles_Nov03_YukawaSummed'

# Problem with SM file
# derivedTheoryFiles_Top                      = 'out/derivedTheoryFiles_Nov06_Top'
derivedTheoryFiles_Top                      = 'out/derivedTheoryFiles_Nov24_Top'

derivedTheoryFiles_TopHighPt                = 'out/derivedTheoryFiles_Dec07_TopHighPt'


########################################
# Correlation matrices and uncertainties
########################################

correlationMatrix_Yukawa              = 'out/correlationMatrices_Nov03/corrMat_exp.txt'
theoryUncertainties_Yukawa            = 'out/correlationMatrices_Nov03/errors_for_corrMat_exp.txt'
correlationMatrix_Yukawa_Uncorrelated = 'out/correlationMatrices_Nov03/corrMat_exp_UNCORRELATED.txt'

# Problem with SM file
# correlationMatrix_Top                 = 'out/correlationMatrices_Nov06_Top/corrMat_exp.txt'
# theoryUncertainties_Top               = 'out/correlationMatrices_Nov06_Top/errors_for_corrMat_exp.txt'
correlationMatrix_Top                 = 'out/correlationMatrices_Nov24_Top/corrMat_exp.txt'
theoryUncertainties_Top               = 'out/correlationMatrices_Nov24_Top/errors_for_corrMat_exp.txt'


# Correlation matrices for differential observables
# correlationMatrix_PTH                 = 'out/corrMat_Nov15_combinedCard_smH_Nov07/higgsCombine_CORRMAT_combinedCard_smH_Nov07.MultiDimFit.mH125.root'
# correlationMatrix_NJ                  = 'out/corrMat_Nov15_combinedCard_NJ_smH_Nov12/higgsCombine_CORRMAT_combinedCard_NJ_smH_Nov12.MultiDimFit.mH125.root'
correlationMatrix_PTH                   = 'out/corrMat_Nov29_combinedCard_smH_Nov07/higgsCombine_CORRMAT_combinedCard_smH_Nov07.MultiDimFit.mH125.root'
correlationMatrix_PTH_ggH               = 'out/corrMat_Nov29_combinedCard_Nov03_xHfixed_1/higgsCombine_CORRMAT_combinedCard_Nov03_xHfixed.MultiDimFit.mH125.root'
correlationMatrix_NJ                    = 'out/corrMat_Nov29_combinedCard_NJ_smH_Nov12/higgsCombine_CORRMAT_combinedCard_NJ_smH_Nov12.MultiDimFit.mH125.root'
correlationMatrix_YH                    = 'out/corrMat_Nov29_combinedCard_YH_smH_Nov28/higgsCombine_CORRMAT_combinedCard_YH_smH_Nov28.MultiDimFit.mH125.root'
correlationMatrix_PTJ                   = 'out/corrMat_Nov30_combinedCard_PTJ_smH_Nov28/higgsCombine_CORRMAT_combinedCard_PTJ_smH_Nov28.MultiDimFit.mH125.root'


########################################
# Workspaces after running t2ws
########################################

# ======================================
# Unsplit workspaces for non-coupling combinations

ws_hgg_smH                   = 'out/workspaces_Nov07/Datacard_13TeV_differential_PtNNLOPS_systs_renamedProcesses.root'
ws_hzz_smH                   = 'out/workspaces_Nov08/hzz4l_comb_13TeV_xs_processesShifted.root'
ws_combined_smH              = 'out/workspaces_Nov07/combinedCard_smH_Nov07.root'

ws_hgg_smH_NJ                = 'out/workspaces_Nov12/Datacard_13TeV_differential_Njets2p5NNLOPS_renamedProcesses.root'
ws_hzz_smH_NJ                = 'out/workspaces_Nov12/hzz4l_comb_13TeV_xs.root'
ws_combined_smH_NJ           = 'out/workspaces_Nov12/combinedCard_NJ_smH_Nov12.root'

ws_hgg_smH_YH                = 'out/workspaces_Nov28/Datacard_13TeV_differential_AbsRapidityNNLOPS_newBins_combination_renamedProcesses.root'
ws_hzz_smH_YH                = 'out/workspaces_Nov28/hzz4l_comb_13TeV_xs.root'
ws_combined_smH_YH           = 'out/workspaces_Nov28/combinedCard_YH_smH_Nov28.root'

ws_hzz_smH_PTJ               = 'out/workspaces_Nov28/hzz4l_comb_13TeV_xs_ptjet.root'

# Wrong variable name (causes problems with file-naming, although physics is correct)
# ws_hgg_smH_PTJ               = 'out/workspaces_Nov28/Datacard_13TeV_differential_Jet2p5Pt0NNLOPS_newBins_renamedProcesses.root'
# ws_combined_smH_PTJ          = 'out/workspaces_Nov28/combinedCard_PTJ_smH_Nov28.root'
ws_hgg_smH_PTJ               = 'out/workspaces_Nov30/Datacard_13TeV_differential_Jet2p5Pt0NNLOPS_newBins_renamedProcesses.root'
ws_combined_smH_PTJ          = 'out/workspaces_Nov30/combinedCard_PTJ_smH_Nov28.root'

ws_combined_ggH_xHfixed      = 'out/workspaces_Nov08/combinedCard_Nov03_xHfixed.root'
ws_hgg_ggH_xHfixed           = 'out/workspaces_Nov29/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_xHfixed.root'
ws_hzz_ggH_xHfixed           = 'out/workspaces_Nov29/hzz4l_comb_13TeV_xs_processesShifted_xHfixed.root'

# ---------------------
# Workspaces for extra studies

# ws_ratio_of_BRs              = 'out/workspaces_Nov08/combinedCard_smH_Nov07_FitBRModel.root' # Probably not wrong but safer to use new result
ws_ratio_of_BRs              = 'out/workspaces_Dec11/combinedCard_smH_Nov07_FitBRModel.root'
ws_ratio_of_BRs_globalScales = 'out/workspaces_Nov20/combinedCard_smH_Nov07_FitBRModel_globalScales.root'
ws_totalXS                   = 'out/workspaces_Nov08/combinedCard_smH_Nov07_FitBRModel_fitTotalXS.root'
ws_combined_ratioOfBRs       = 'out/workspaces_Nov14/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_ratioOfBRs.root'
ws_combined_totalXS          = 'out/workspaces_Nov27/combinedCard_smH_Nov27_INCLUSIVE_FitBRModel_fitTotalXS.root'


# ======================================
# kappab kappac

ws_combined_Yukawa = 'out/workspaces_Nov03/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties.root'
ws_hgg_Yukawa      = 'out/workspaces_Nov03/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_Yukawa_withTheoryUncertainties.root'
ws_hzz_Yukawa      = 'out/workspaces_Nov03/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_Yukawa_withTheoryUncertainties.root'

ws_combined_Yukawa_noTheoryUncertainties               = 'out/workspaces_Nov06/combinedCard_Nov03_CouplingModel_Yukawa_noTheoryUncertainties.root'
ws_combined_Yukawa_withUncorrelatedTheoryUncertainties = 'out/workspaces_Nov06/combinedCard_Nov03_CouplingModel_Yukawa_withUncorrelatedTheoryUncertainties.root'
ws_combined_Yukawa_lumiScalable                        = 'out/workspaces_Nov06/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_lumiScale.root'
ws_combined_Yukawa_profiledTotalXS                     = 'out/workspaces_Nov09/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_profiledTotalXS.root'
ws_combined_Yukawa_couplingDependentBR                 = 'out/workspaces_Nov10/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_couplingDependentBR.root'
# ws_combined_Yukawa_profiledTotalXS_fitOnlyNormalization = 'out/workspaces_Nov17/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_profiledTotalXS_fitOnlyNormalization.root'
# ws_combined_Yukawa_profiledTotalXS_fitOnlyNormalization = 'out/workspaces_Nov20/combinedCard_Nov03_CouplingModel_Yukawa_withTheoryUncertainties_profiledTotalXS_fitOnlyNormalization.root'
ws_combined_Yukawa_profiledTotalXS_fitOnlyNormalization = 'out/workspaces_Nov28/combinedCard_smH_Nov27_INCLUSIVE_CouplingModel_Yukawa_profiledTotalXS_fitOnlyNormalization.root'


# ======================================
# kappat kappag

# Problem with last theory uncertainty:
# ws_combined_Top    = 'out/workspaces_Nov06/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties.root'
# ws_hgg_Top         = 'out/workspaces_Nov06/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_Top_withTheoryUncertainties.root'
# ws_hzz_Top         = 'out/workspaces_Nov06/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_Top_withTheoryUncertainties.root'

# Last theory uncertainty skipped:  # Had problem with SM file
# ws_hgg_Top         = 'out/workspaces_Nov08/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_Top_withTheoryUncertainties.root'
# ws_combined_Top    = 'out/workspaces_Nov08/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties.root'
# ws_hzz_Top         = 'out/workspaces_Nov08/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_Top_withTheoryUncertainties.root'

# Problem with SM file resolved, last theory bin still skipped
ws_combined_Top    = 'out/workspaces_Nov24/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties.root'
ws_hgg_Top         = 'out/workspaces_Nov24/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_Top_withTheoryUncertainties.root'
ws_hzz_Top         = 'out/workspaces_Nov24/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_Top_withTheoryUncertainties.root'

ws_combined_TopHighPt = 'out/workspaces_Dec11/combinedCard_Nov03_CouplingModel_TopHighPt_withTheoryUncertainties.root'
ws_hgg_TopHighPt      = 'out/workspaces_Dec11/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_TopHighPt_withTheoryUncertainties.root'
ws_hzz_TopHighPt      = 'out/workspaces_Dec11/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_TopHighPt_withTheoryUncertainties.root'

# ws_combined_Top_lumiScalable           = 'out/workspaces_Nov09/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_lumiScale.root'
# ws_combined_Top_profiledTotalXS        = 'out/workspaces_Nov09/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_profiledTotalXS.root'
ws_combined_Top_couplingDependentBR    = 'out/workspaces_Nov10/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_couplingDependentBR.root'
ws_combined_Top_skippedLastBin         = 'out/workspaces_Nov17/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_skippedLastBin.root'

# ws_combined_Top_profiledTotalXS_fitOnlyNormalization = 'out/workspaces_Nov17/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_profiledTotalXS_fitOnlyNormalization.root'
ws_combined_Top_profiledTotalXS_fitOnlyNormalization = 'out/workspaces_Nov28/combinedCard_smH_Nov27_INCLUSIVE_CouplingModel_Top_profiledTotalXS_fitOnlyNormalization.root'
ws_combined_TopHighPt_profiledTotalXS_fitOnlyNormalization = 'out/workspaces_Dec12/combinedCard_smH_Nov27_INCLUSIVE_CouplingModel_Top_profiledTotalXS_fitOnlyNormalization.root'

# Without problematic SM norm
ws_combined_Top_lumiScalable           = 'out/workspaces_Nov27/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_lumiScale.root'
ws_combined_Top_profiledTotalXS        = 'out/workspaces_Nov27/combinedCard_Nov03_CouplingModel_Top_withTheoryUncertainties_profiledTotalXS.root'


ws_combined_TopCtCb = 'out/workspaces_Nov15/combinedCard_Nov03_CouplingModel_TopCtCb_withTheoryUncertainties.root'
ws_hgg_TopCtCb      = 'out/workspaces_Nov15/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_TopCtCb_withTheoryUncertainties.root'
ws_hzz_TopCtCb      = 'out/workspaces_Nov15/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_TopCtCb_withTheoryUncertainties.root'

ws_combined_TopCtCbHighPt = 'out/workspaces_Dec12/combinedCard_Nov03_CouplingModel_TopCtCbHighPt_withTheoryUncertainties.root'
ws_hgg_TopCtCbHighPt      = 'out/workspaces_Dec12/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_renamedProcesses_CouplingModel_TopCtCbHighPt_withTheoryUncertainties.root'
ws_hzz_TopCtCbHighPt      = 'out/workspaces_Dec12/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_TopCtCbHighPt_withTheoryUncertainties.root'


########################################
# Scans
########################################

# ======================================
# Non-coupling combinations

scan_combined_PTH = 'out/Scan_PTH_Nov07'
scan_hgg_PTH      = 'out/Scan_PTH_Nov07_hgg'
scan_hzz_PTH      = 'out/Scan_PTH_Nov08_hzz'

scan_combined_PTH_ggH = 'out/Scan_PTH_Nov08_xHfixed'
scan_hgg_PTH_ggH      = 'out/Scan_PTH_Nov29_xHfixed_hgg'
scan_hzz_PTH_ggH      = 'out/Scan_PTH_Nov29_xHfixed_hzz'

scan_combined_NJ  = 'out/Scan_nJets_Nov12'
scan_hgg_NJ       = 'out/Scan_nJets_Nov12_0'
scan_hzz_NJ       = 'out/Scan_nJets_Nov12_1'

scan_combined_YH  = 'out/Scan_YH_Nov28'
scan_hgg_YH       = 'out/Scan_YH_Nov28_0'
scan_hzz_YH       = 'out/Scan_YH_Nov28_1'

scan_combined_PTJ = 'out/Scan_PTJ_Nov28_1'
scan_hgg_PTJ      = 'out/Scan_PTJ_Nov28_0'
scan_hzz_PTJ      = 'out/Scan_PTJ_Nov28'

scan_combined_PTH_ggH_asimov = 'out/Scan_PTH_Nov17_xHfixed_asimov'

scan_combined_YH_asimov  = 'out/Scan_YH_Nov28_asimov'
scan_hgg_YH_asimov       = 'out/Scan_YH_Nov28_asimov_0'
scan_hzz_YH_asimov       = 'out/Scan_YH_Nov28_asimov_1'

# scan_combined_PTJ_asimov = 'out/Scan_PTJ_Nov28_asimov'
# scan_hgg_PTJ_asimov      = 'out/Scan_PTJ_Nov28_asimov_0'
scan_combined_PTJ_asimov = 'out/Scan_PTJ_Nov30_asimov'
scan_hgg_PTJ_asimov      = 'out/Scan_PTJ_Nov30_asimov_0'
scan_hzz_PTJ_asimov      = 'out/Scan_PTJ_Nov28_asimov_1'


# ======================================
# Extra studies

# scan_ratioOfBRs                  = 'out/Scan_ratioOfBRs_Nov08' # Probably not wrong, but renewed to be safe
scan_ratioOfBRs                  = 'out/Scan_ratioOfBRs_Dec11'
scan_ratioOfBRs_globalScales     = 'out/Scan_ratioOfBRs_Nov20_globalScales'
# scan_combined_totalXS            = 'out/Scan_TotalXS_Nov08'
# scan_combined_PTH_xHfixed        = 'out/Scan_PTH_Nov08_xHfixed'
# scan_combined_PTH_xHfixed_asimov = 'out/Scan_PTH_Nov17_xHfixed_asimov'

scan_combined_totalXS            = 'out/Scan_TotalXS_Nov27'


# ======================================
# kappab kappac scans

scan_combined_Yukawa                                        = 'out/Scan_Yukawa_Nov03_0'
scan_hgg_Yukawa                                             = 'out/Scan_Yukawa_Nov03_hgg'
scan_hzz_Yukawa                                             = 'out/Scan_Yukawa_Nov03_hzz'

scan_combined_Yukawa_asimov                                 = 'out/Scan_Yukawa_Nov03_asimov'
scan_hgg_Yukawa_asimov                                      = 'out/Scan_Yukawa_Nov03_hgg_asimov'
scan_hzz_Yukawa_asimov                                      = 'out/Scan_Yukawa_Nov03_hzz_asimov'

# scan_combined_Yukawa_fitOnlyNormalization_asimov            = 'out/Scan_Yukawa_Nov17_profiledTotalXS_fitOnlyNormalization_asimov'
scan_combined_Yukawa_fitOnlyNormalization                   = 'out/Scan_Yukawa_Nov20_profiledTotalXS_fitOnlyNormalization'
# scan_combined_Yukawa_fitOnlyNormalization_asimov            = 'out/Scan_Yukawa_Nov20_profiledTotalXS_fitOnlyNormalization_asimov'
scan_combined_Yukawa_fitOnlyNormalization_asimov            = 'out/Scan_Yukawa_Nov28_fitOnlyNormalization_asimov_1'

scan_combined_Yukawa_lumiStudy_asimov                       = 'out/Scan_Yukawa_Nov06_lumiStudy_asimov'
scan_combined_Yukawa_noTheoryUncertainties_asimov           = 'out/Scan_Yukawa_Nov06_noTheoryUncertainties_asimov'
scan_combined_Yukawa_uncorrelatedTheoryUncertainties_asimov = 'out/Scan_Yukawa_Nov06_uncorrelatedTheoryUncertainties_asimov'
scan_combined_Yukawa_profiledTotalXS_asimov                 = 'out/Scan_Yukawa_Nov09_profiledTotalXS_asimov'

scan_combined_Yukawa_oneKappa_kappac                        = 'out/Scan_Yukawa_Nov06_oneKappa_kappac'
scan_combined_Yukawa_oneKappa_kappac_asimov                 = 'out/Scan_Yukawa_Nov06_oneKappa_kappac_asimov'
scan_combined_Yukawa_oneKappa_kappab                        = 'out/Scan_Yukawa_Nov06_oneKappa_kappab'
scan_combined_Yukawa_oneKappa_kappab_asimov                 = 'out/Scan_Yukawa_Nov06_oneKappa_kappab_asimov'

scan_combined_Yukawa_couplingDependentBR_asimov              = 'out/Scan_Yukawa_Nov10_couplingDependentBR_asimov'
scan_combined_Yukawa_couplingDependentBR_fixedKappaV_asimov  = 'out/Scan_Yukawa_Nov10_couplingDependentBR_fixedKappaV_asimov'
# scan_combined_Yukawa_couplingDependentBR_kappaVMaxOne_asimov = 'out/Scan_Yukawa_Nov16_couplingDependentBR_kappaVMaxOne_asimov'
scan_combined_Yukawa_couplingDependentBR_kappaVMaxOne_asimov = 'out/Scan_Yukawa_Nov17_couplingDependentBR_kappaVMaxOne_asimov'
scan_combined_Yukawa_ratioOfBRs_asimov                       = 'out/Scan_Yukawa_Nov14_ratioOfBRs_asimov'

scan_combined_Yukawa_ratioOfBRs_onedimRatioScan              = 'out/Scan_Yukawa_Nov15_ratioOfBRs_onedimRatioScan_0'
scan_combined_Yukawa_ratioOfBRs_onedimRatioScan_asimov       = 'out/Scan_Yukawa_Nov15_ratioOfBRs_onedimRatioScan_asimov_0'

scan_combined_Yukawa_profiledTotalXS_onedimTotalXSScan        = 'out/Scan_Yukawa_Nov16_profiledTotalXS_onedimTotalXSScan'
scan_combined_Yukawa_profiledTotalXS_onedimTotalXSScan_asimov = 'out/Scan_Yukawa_Nov16_profiledTotalXS_onedimTotalXSScan_asimov'


# ======================================
# kappat kappag scans

# Problems with last theory uncertainty
# scan_combined_Top        = 'out/Scan_Top_Nov06'
# scan_hzz_Top             = 'out/Scan_Top_Nov06_0'
# scan_hgg_Top             = 'out/Scan_Top_Nov06_1'
# scan_combined_Top_asimov = 'out/Scan_Top_Nov07_asimov'
# scan_hzz_Top_asimov      = 'out/Scan_Top_Nov07_hzz_asimov'
# scan_hgg_Top_asimov      = 'out/Scan_Top_Nov07_hgg_asimov'

# scan_combined_Top                                         = 'out/Scan_Top_Nov08'
# scan_hgg_Top                                              = 'out/Scan_Top_Nov08_hgg'
# scan_hzz_Top                                              = 'out/Scan_Top_Nov08_hzz'
# scan_hgg_Top                                              = 'out/Scan_Top_Nov27_hgg_0'

# scan_combined_Top                                         = 'out/Scan_Top_Nov27_0'
# scan_hgg_Top                                              = 'out/Scan_Top_Nov28_hgg_0'
# scan_hzz_Top                                              = 'out/Scan_Top_Nov26_hzz'

# scan_hzz_Top                                              = 'out/Scan_Top_Dec03_hzz'

scan_combined_Top                                         = 'out/Scan_Top_Dec03_2'
scan_hgg_Top                                              = 'out/Scan_Top_Dec03_hgg_0'
scan_hzz_Top                                              = 'out/Scan_Top_Dec03_hzz_0'

# scan_combined_Top_asimov                                  = 'out/Scan_Top_Nov08_asimov'
# scan_hgg_Top_asimov                                       = 'out/Scan_Top_Nov08_hgg_asimov'
# scan_hzz_Top_asimov                                       = 'out/Scan_Top_Nov08_hzz_asimov'
scan_combined_Top_asimov                                  = 'out/Scan_Top_Nov27_asimov_1'
scan_hgg_Top_asimov                                       = 'out/Scan_Top_Nov26_hgg_asimov'
scan_hzz_Top_asimov                                       = 'out/Scan_Top_Nov26_hzz_asimov'

scan_combined_TopHighPt                                   = 'out/Scan_TopHighPt_Dec11'
scan_hgg_TopHighPt                                        = 'out/Scan_TopHighPt_Dec11_hgg'
scan_hzz_TopHighPt                                        = 'out/Scan_TopHighPt_Dec11_hzz_0'

scan_combined_TopHighPt_asimov                            = 'out/Scan_TopHighPt_Dec12_asimov_2'

# TODO:
# Scan_TotalXS_Nov27

scan_combined_Top_extendedRange_asimov                    = 'out/Scan_Top_Nov20_asimov'
scan_hgg_Top_extendedRange_asimov                         = 'out/Scan_Top_Nov20_hzz_asimov'
scan_hzz_Top_extendedRange_asimov                         = 'out/Scan_Top_Nov20_hgg_asimov'

scan_combined_Top_fitOnlyNormalization                    = 'out/Scan_Top_Nov17_profiledTotalXS_fitOnlyNormalization'
# scan_combined_Top_fitOnlyNormalization_asimov             = 'out/Scan_Top_Nov17_profiledTotalXS_fitOnlyNormalization_asimov'
scan_combined_Top_bigRange                                = 'out/Scan_Top_Nov17'
scan_combined_Top_bigRange2                               = 'out/Scan_Top_Nov19'
scan_combined_Top_skippedLastBin                          = 'out/Scan_Top_Nov17_skippedLastBin'
scan_combined_Top_skippedLastBin_asimov                   = 'out/Scan_Top_Nov17_skippedLastBin_asimov'
scan_combined_Top_lumiStudy_asimov                        = 'out/Scan_Top_Nov09_lumiStudy_asimov'
# scan_combined_Top_profiledTotalXS_asimov                  = 'out/Scan_Top_Nov09_profiledTotalXS_asimov'
scan_combined_Top_couplingDependentBR_asimov              = 'out/Scan_Top_Nov10_couplingDependentBR_asimov'
scan_combined_Top_couplingDependentBR_fixedKappaV_asimov  = 'out/Scan_Top_Nov10_couplingDependentBR_fixedKappaV_asimov'
# scan_combined_Top_couplingDependentBR_kappaVMaxOne_asimov = 'out/Scan_Top_Nov16_couplingDependentBR_kappaVMaxOne_asimov'
scan_combined_Top_couplingDependentBR_kappaVMaxOne_asimov = 'out/Scan_Top_Nov17_couplingDependentBR_kappaVMaxOne_asimov'
scan_combined_Top_couplingDependentBR_bigRange_asimov     = 'out/Scan_Top_Nov19_couplingDependentBR_asimov'


scan_combined_Top_fitOnlyNormalization_asimov             = 'out/Scan_Top_Nov28_fitOnlyNormalization_asimov_0'
scan_combined_Top_profiledTotalXS_asimov                  = 'out/Scan_Top_Nov27_profiledTotalXS_asimov_0'
scan_combined_Top_lumiStudy_asimov                        = 'out/Scan_Top_Nov27lumiStudy_asimov_1'


scan_combined_TopCtCb                                    = 'out/Scan_TopCtCb_Nov15'
scan_hgg_TopCtCb                                         = 'out/Scan_TopCtCb_Nov15_hgg'
scan_hzz_TopCtCb                                         = 'out/Scan_TopCtCb_Nov15_hzz'
scan_combined_TopCtCb_asimov                             = 'out/Scan_TopCtCb_Nov15_asimov'
scan_hgg_TopCtCb_asimov                                  = 'out/Scan_TopCtCb_Nov15_hgg_asimov'
scan_hzz_TopCtCb_asimov                                  = 'out/Scan_TopCtCb_Nov15_hzz_asimov'



########################################
# End
########################################
if __name__ == "__main__":
    keys = [ k for k in vars().keys() ]
    for key in keys:
        if key.startswith('__'): continue
        print '{0}: {1}'.format( key, vars()[key] )