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

# Lumi scalable
ws_combined_unsplit_lumiScalableWS    = 'workspaces_Aug18/combinedCard_Jul26_CouplingModel_Yukawa_withTheoryUncertainties_lumiScale.root'


ws_onlyhgg_split           = 'workspaces_Aug21/Datacard_13TeV_differential_PtGghPlusHx_renamedProcesses.root'
ws_combined_split          = 'workspaces_Aug21/combinedCard_Aug21.root'


ws_combined_split_top      = 'workspaces_Aug21/combinedCard_Aug21_CouplingModel_Top_withTheoryUncertainties.root'
ws_onlyhgg_split_top       = 'workspaces_Aug22/Datacard_13TeV_differential_PtGghPlusHx_renamedProcesses_CouplingModel_Top_withTheoryUncertainties.root'
ws_onlyhzz_split_top       = 'workspaces_Aug22/hzz4l_all_13TeV_xs_processesShifted_CouplingModel_Top_withTheoryUncertainties.root'


# ======================================
# Scans

# ---------------------------
# Not-coupling results (yieldParameter per bin)

scan_ptcombination_combined_profiled_asimov = 'Scan_Aug22_1'
scan_ptcombination_hgg_profiled_asimov      = 'Scan_Aug22_2'
scan_ptcombination_hzz_profiled_asimov      = 'Scan_Aug22_3'

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

# ---------------------------
# Top results

scan_top_combined_profiled           = 'Scan_Top_Aug21_0'
scan_top_hzz_profiled                = 'Scan_Top_Aug22'


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