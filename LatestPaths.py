YR4_totalXS = 55.70628722 # pb

# ======================================
# Text card paths

card_onlyhgg_unsplit_unrenamed       = 'suppliedInput/fromVittorio/pT_unsplit_Mar21/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5.txt'
card_onlyhgg_unsplit_renamed         = 'suppliedInput/fromVittorio/pT_unsplit_Mar21/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamedProcesses_Jul24.txt'

card_onlyhgg_split_ggHonly_unrenamed = 'suppliedInput/fromVittorio/pT_ggHonly_Jun26/Datacard_13TeV_differential_pT_moriond17_ggHonly_v2.txt'
card_onlyhgg_split_ggHonly_reparsed  = 'suppliedInput/fromVittorio/Datacard_13TeV_differential_pT_moriond17_ggHonly_v2_reparsed.txt'
card_onlyhgg_split_xHonly_unrenamed  = 'suppliedInput/fromVittorio/differential_pT_moriond17_HxOnly_July19/Datacard_13TeV_differential_pT_moriond17_HxOnly.txt'
card_onlyhgg_split_xHonly_reparsed   = 'suppliedInput/fromVittorio/Datacard_13TeV_differential_pT_moriond17_HxOnly_reparsed.txt'
card_onlyhgg_split_ggHxHmerged       = 'suppliedInput/fromVittorio/Merged_Aug03.txt'

card_onlyhgg_split_both_unrenamed    = 'suppliedInput/fromVittorio/differential_PtGghPlusHx_Aug21/Datacard_13TeV_differential_PtGghPlusHx.txt'
card_onlyhgg_split_both_renamed      = 'suppliedInput/fromVittorio/differential_PtGghPlusHx_Aug21/Datacard_13TeV_differential_PtGghPlusHx_renamedProcesses.txt'

card_onlyhzz_unsplit_unrenamed       = 'suppliedInput/fromDavid/PTH_May15/hzz4l_comb_13TeV_xs.txt'
card_onlyhzz_unsplit_OAsignal        = 'suppliedInput/fromDavid/PTH_May15/hzz4l_comb_13TeV_xs_processesShifted.txt'

card_onlyhzz_split_unrenamed         = 'suppliedInput/fromDavid/ggHonly_Jun26/hzz4l_all_13TeV_xs.txt'
card_onlyhzz_split_renamed           = 'suppliedInput/fromDavid/ggHonly_Jun26/hzz4l_all_13TeV_xs_processesShifted.txt'

card_combined_unsplit                = 'suppliedInput/combinedCard_Jul26.txt'
card_combined_split                  = 'suppliedInput/combinedCard_Aug21.txt'

# -----------------
# nJets

card_njets_hgg_unsplit_unrenamed     = 'suppliedInput/fromVittorio/nJets2p5_v5_Mar15/Datacard_13TeV_differential_Njets_moriond17_skipAndDebug_reminiaod_corrections_genJetID_v5.txt'
card_njets_hgg_unsplit_renamed       = 'suppliedInput/fromVittorio/nJets2p5_v5_Mar15/Datacard_13TeV_differential_Njets_moriond17_skipAndDebug_reminiaod_corrections_genJetID_v5_renamedProcesses.txt'

card_njets_hzz_unsplit               = 'suppliedInput/fromDavid/nJets_Sep19/smH/hzz4l_comb_13TeV_xs.txt'

card_njets_combined_unsplit          = 'suppliedInput/combinedCard_nJets_Sep19.txt'



# ======================================
# Workspaces after text2workspace

ws_onlyhgg_unsplit         = 'workspaces_Aug12/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamedProcesses_Jul24.root'
ws_onlyhzz_unsplit         = 'workspaces_Aug13/hzz4l_comb_13TeV_xs_processesShifted.root'
ws_combined_unsplit        = 'workspaces_Aug13/combinedCard_Jul26.root'
# ws_combined_unsplit        = 'workspaces_May30/combinedCard_May15.root'

ws_onlyhzz_unsplit_yukawa             = 'workspaces_Aug13/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_Yukawa_withTheoryUncertainties.root'
ws_onlyhgg_unsplit_yukawa_fulllength  = 'workspaces_Aug13/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamedProcesses_Jul24_CouplingModel_Yukawa_withTheoryUncertainties.root'
ws_onlyhgg_unsplit_yukawa             = 'workspaces_Aug13/hgg_renamedProcesses_Jul24_CouplingModel_Yukawa_withTheoryUncertainties.root'

ws_onlyhzz_unsplit_yukawa             = 'workspaces_Aug13/hzz4l_comb_13TeV_xs_processesShifted_CouplingModel_Yukawa_withTheoryUncertainties.root'

ws_combined_unsplit_yukawa            = 'workspaces_Aug10/combinedCard_Jul26_CouplingModel_Yukawa_withTheoryUncertainties.root'
ws_combined_unsplit_yukawa_onlyGluonInduced = 'workspaces_Aug17/combinedCard_Jul26_CouplingModel_Yukawa_withTheoryUncertainties.root'

# Has proper theory uncertainties as well now
ws_combined_split_yukawa              = 'workspaces_Sep14/combinedCard_Aug21_CouplingModel_Yukawa_withTheoryUncertainties.root'
ws_onlyhgg_split_yukawa               = 'workspaces_Sep14/Datacard_13TeV_differential_PtGghPlusHx_renamedProcesses_CouplingModel_Yukawa_withTheoryUncertainties.root'
# ws_onlyhzz_split_yukawa               = 'workspaces_Sep15/hzz4l_all_13TeV_xs_processesShifted_CouplingModel_Yukawa_withTheoryUncertainties.root'
# ws_onlyhzz_split_yukawa               = 'workspaces_Sep18/hzz4l_all_13TeV_xs_processesShifted_CouplingModel_Yukawa_withTheoryUncertainties.root'
ws_onlyhzz_split_yukawa               = 'workspaces_Sep19/hzz4l_all_13TeV_xs_processesShifted_CouplingModel_Yukawa_withTheoryUncertainties.root'

# Lumi scalable
ws_combined_unsplit_lumiScalableWS    = 'workspaces_Aug18/combinedCard_Jul26_CouplingModel_Yukawa_withTheoryUncertainties_lumiScale.root'
ws_combined_split_lumiScalableWS      = 'workspaces_Sep19/combinedCard_Aug21_CouplingModel_Yukawa_withTheoryUncertainties_lumiScale.root'

ws_onlyhzz_split         = 'workspaces_Sep15/hzz4l_all_13TeV_xs_processesShifted.root'
ws_onlyhgg_split           = 'workspaces_Aug21/Datacard_13TeV_differential_PtGghPlusHx_renamedProcesses.root'
ws_combined_split          = 'workspaces_Aug21/combinedCard_Aug21.root'


ws_combined_split_top      = 'workspaces_Aug21/combinedCard_Aug21_CouplingModel_Top_withTheoryUncertainties.root'
ws_onlyhgg_split_top       = 'workspaces_Aug22/Datacard_13TeV_differential_PtGghPlusHx_renamedProcesses_CouplingModel_Top_withTheoryUncertainties.root'
ws_onlyhzz_split_top       = 'workspaces_Aug22/hzz4l_all_13TeV_xs_processesShifted_CouplingModel_Top_withTheoryUncertainties.root'

ws_combined_split_top_lumiScalableWS   = 'workspaces_Sep20/combinedCard_Aug21_CouplingModel_Top_withTheoryUncertainties_lumiScale.root'

# NOT ACTUALLY TOP!! CHANGE TO YUKAWA!!
ws_combined_split_top_notheoryunc  = 'workspaces_Sep19/combinedCard_Aug21_CouplingModel_Yukawa.root'
ws_combined_split_top_uncorrelated = 'workspaces_Sep19/combinedCard_Aug21_CouplingModel_Yukawa_withUncorrelatedTheoryUncertainties.root'


ws_njets_combined_unsplit = 'workspaces_Sep19/combinedCard_nJets_Sep19.root'
ws_njets_hgg_unsplit      = 'workspaces_Sep19/Datacard_13TeV_differential_Njets_moriond17_skipAndDebug_reminiaod_corrections_genJetID_v5_renamedProcesses.root'
ws_njets_hzz_unsplit      = 'workspaces_Sep19/hzz4l_comb_13TeV_xs.root'

ws_FitBR_combined_unsplit = 'workspaces_Sep25/combinedCard_Jul26_FitBRModel.root'

ws_combined_split_top_couplingDependentBR = 'workspaces_Sep27/combinedCard_Aug21_CouplingModel_Top_withTheoryUncertainties_couplingDependentBR.root'


# Also scaling xH with the BR modifier, and better binning decisions w.r.t. underflow and overflow
ws_combined_split_betterTop                        = 'workspaces_Sep28/combinedCard_Aug21_CouplingModel_Top_withTheoryUncertainties.root'
# ws_combined_split_betterTop_couplingDependentBR    = 'workspaces_Sep28/combinedCard_Aug21_CouplingModel_Top_withTheoryUncertainties_couplingDependentBR.root'
ws_combined_split_betterTop_couplingDependentBR    = 'workspaces_Sep29/combinedCard_Aug21_CouplingModel_Top_withTheoryUncertainties_couplingDependentBR.root'

ws_combined_split_betterYukawa                     = 'workspaces_Sep28/combinedCard_Aug21_CouplingModel_Yukawa_withTheoryUncertainties.root'
# ws_combined_split_betterYukawa_couplingDependentBR = 'workspaces_Sep28/combinedCard_Aug21_CouplingModel_Yukawa_withTheoryUncertainties_couplingDependentBR.root'
# ws_combined_split_betterYukawa_couplingDependentBR = 'workspaces_Sep29/combinedCard_Aug21_CouplingModel_Yukawa_withTheoryUncertainties_couplingDependentBR.root'
ws_combined_split_betterYukawa_couplingDependentBR = 'workspaces_Oct02/combinedCard_Aug21_CouplingModel_Yukawa_withTheoryUncertainties_couplingDependentBR.root'



# ======================================
# Derived theory files

derivedTheoryFilesDirectory_Top                      = 'derivedTheoryFiles_Aug11_Top'

derivedTheoryFilesDirectory_YukawaGluonInduced       = 'derivedTheoryFiles_Aug09_YukawaGluonInduced'
derivedTheoryFilesDirectory_YukawaQuarkInduced       = 'derivedTheoryFiles_Aug09_YukawaQuarkInduced'
derivedTheoryFilesDirectory_YukawaQuarkInducedScaled = 'derivedTheoryFiles_Sep14_YukawaQuarkInducedScaled'
# derivedTheoryFilesDirectory_YukawaSummed             = 'derivedTheoryFiles_Aug09_YukawaSummed'
derivedTheoryFilesDirectory_YukawaSummed             = 'derivedTheoryFiles_Sep14_YukawaSummed'


# ======================================
# Correlation matrices and uncertainties

correlationMatrix_Yukawa_uncorrelated = 'plots_CorrelationMatrices_Sep14/corrMat_exp_uncorrelated.txt'
correlationMatrix_Yukawa              = 'plots_CorrelationMatrices_Sep14/corrMat_exp.txt'
theoryUncertainties_Yukawa            = 'plots_CorrelationMatrices_Sep14/errors_for_corrMat_exp.txt'
correlationMatrix_Top                 = 'plots_CorrelationMatrices_Aug11_Top/corrMat_exp.txt'
theoryUncertainties_Top               = 'plots_CorrelationMatrices_Aug11_Top/errors_for_corrMat_exp.txt'

correlationMatrix_YukawaGluon         = 'plots_CorrelationMatrices_Aug09/corrMat_exp.txt'
theoryUncertainties_YukawaGluon       = 'plots_CorrelationMatrices_Aug09/errors_for_corrMat_exp.txt'

# ======================================
# Scans

# ---------------------------
# Not-coupling results (yieldParameter per bin)

scan_ptcombination_combined_profiled_asimov = 'Scan_Aug22_1'
scan_ptcombination_hgg_profiled_asimov      = 'Scan_Aug22_2'
scan_ptcombination_hzz_profiled_asimov      = 'Scan_Aug22_3'

scan_ptcombination_combined_profiled        = 'Scan_May15'
scan_ptcombination_hgg_profiled             = 'Scan_May15'
scan_ptcombination_hzz_profiled             = 'Scan_May15'

# ---------------------------
# Yukawa results

# hgg:
scan_hgg_fastscan      = 'Scan_couplings_Aug13_0'
scan_hgg_profiled      = 'Scan_couplings_Aug14_0'

# hzz:
scan_hzz_fastscan      = 'Scan_couplings_Aug13'
scan_hzz_profiled      = 'Scan_couplings_Aug14'

# combination
scan_combined_fastscan = 'Scan_couplings_Aug09_0'
scan_combined_profiled = 'Scan_couplings_Aug09_1'
scan_combined_profiled_addition = 'manual_Scan_coupling_Aug09_1' # Note the missing 's'... now too late

scan_combined_fastscan_asimov        = 'Scan_yukawa_Aug18_asimov'
scan_combined_fastscan_asimov_lumi10 = 'Scan_yukawa_Aug18_asimov_5'
scan_combined_profiled_asimov        = 'Scan_yukawa_Aug18_asimov_1'
scan_combined_profiled_asimov_lum10  = 'Scan_yukawa_Aug18_asimov_6'
scan_combined_profiled_asimov_lum8   = 'Scan_yukawa_Aug21_asimov_0'

scan_yukawa_hzz_profiled_asimov      = 'Scan_yukawa_Aug22_asimov_0'
scan_yukawa_hgg_profiled_asimov      = 'Scan_yukawa_Aug22_asimov_1'

scan_combined_split_fastscan         = 'Scan_yukawa_Sep14'


scan_yukawa_combined_split_profiled         = 'Scan_yukawa_Sep14_0'
scan_yukawa_hgg_split_profiled              = 'Scan_yukawa_Sep14_1'
# scan_yukawa_hzz_split_profiled              = 'Scan_yukawa_Sep18_0'
scan_yukawa_hzz_split_profiled              = 'Scan_yukawa_Sep19'

scan_yukawa_combined_split_profiled_asimov  = 'Scan_yukawa_Sep18_asimov'
scan_yukawa_hgg_split_profiled_asimov       = 'Scan_yukawa_Sep18_asimov_0'
# scan_yukawa_hzz_split_profiled_asimov       = 'Scan_yukawa_Sep18_asimov_1'
scan_yukawa_hzz_split_profiled_asimov       = 'Scan_yukawa_Sep19_asimov'


scan_yukawa_combined_split_profiled_uncorrelated  = 'Scan_yukawa_Sep19_0'
scan_yukawa_combined_split_profiled_notheoryunc   = 'Scan_yukawa_Sep19_1'

scan_yukawa_combined_split_profiled_uncorrelated_asimov  = 'Scan_yukawa_Sep19_asimov_0'
scan_yukawa_combined_split_profiled_notheoryunc_asimov   = 'Scan_yukawa_Sep19_asimov_1'
scan_yukawa_combined_split_profiled_asimov_lum8          = 'Scan_yukawa_Sep19_asimov_2'

# ---------------------------
# Top results

scan_top_combined_profiled           = 'Scan_Top_Aug21_0'
scan_top_hzz_profiled                = 'Scan_Top_Aug22'
scan_top_hgg_profiled                = 'Scan_Top_Aug22_0'

scan_top_combined_profiled_asimov    = 'Scan_Top_Sep20_asimov'
scan_top_hzz_profiled_asimov         = 'Scan_Top_Sep20_asimov_0'
scan_top_hgg_profiled_asimov         = 'Scan_Top_Sep20_asimov_1'
scan_top_combined_profiled_asimov_lum8 = 'Scan_Top_Sep20_asimov_2'
# scan_top_combined_profiled_asimov_lum8 = 'Scan_Top_Sep22_asimov' # New one, wait till finished


scan_top_combined_profiled_asimov_couplingDependentBR = 'Scan_Top_Sep27_asimov'


scan_betterTop_combined_asimov                        = 'Scan_Top_Sep28_asimov'
# scan_betterTop_combined_asimov_couplingDependentBR    = 'Scan_Top_Sep28_asimov_0'
scan_betterTop_combined_asimov_couplingDependentBR    = 'Scan_Top_Sep29_asimov'
scan_betterYukawa_combined_asimov                     = 'Scan_yukawa_Sep28_asimov'
# scan_betterYukawa_combined_asimov_couplingDependentBR = 'Scan_yukawa_Sep28_asimov_0'
scan_betterYukawa_combined_asimov_couplingDependentBR = 'Scan_yukawa_Sep29_asimov'
scan_betterYukawa_combined_asimov_couplingDependentBR_fast = 'Scan_yukawa_Oct02_asimov'


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument( 'varName', type=str, help='Variable to print' )
    args = parser.parse_args()

    import sys
    current_module = sys.modules[__name__]

    if hasattr( current_module, args.varName ):
        print getattr( current_module, args.varName )
    else:
        print 'Variable \'{0}\' was not found in module \'{1}\''.format( args.varName, __name__ )


########################################
# End of Main
########################################
if __name__ == "__main__":
    main()