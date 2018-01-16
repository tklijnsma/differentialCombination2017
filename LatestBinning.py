import sys
sys.path.append('src')
from Observable import Observable

########################################
# Physical quantities
########################################

YR4_ggF_n3lo        = 4.852E+01
YR4_VBF             = 3.779E+00
YR4_WH              = 1.369E+00
YR4_ZH              = 8.824E-01
YR4_ttH             = 5.065E-01
YR4_bbH             = 4.863E-01
YR4_tH_t_ch         = 7.426E-02
YR4_tH_s_ch         = 2.875E-03
YR4_tH_W_associated = 0.000E+00

YR4_totalXS         = YR4_ggF_n3lo + YR4_VBF + YR4_WH + YR4_ZH + YR4_ttH + YR4_bbH + YR4_tH_t_ch + YR4_tH_s_ch + YR4_tH_W_associated
YR4_xH              = YR4_totalXS - YR4_ggF_n3lo

# YR4_totalXS = 55.70628722 # pb
# YR4_ggHXS   = 48.58 # pb

SM_BR_hgg = 2.270E-03
SM_BR_hzz = 0.02641
SM_ratio_of_BRs = SM_BR_hgg / SM_BR_hzz

Observable.YR4_totalXS = YR4_totalXS


########################################
# Binning and SM differentials
########################################

binning_ptjet = [ -1000.0, 30.0, 55.0, 95.0, 120.0, 200.0, 13000.0 ]
binning_yh    = [ 0.0, 0.15, 0.3, 0.6, 0.9, 1.2, 2.5 ]
binning_pth   = [ 0.0, 15.0, 30.0, 45.0, 85.0, 125.0, 200.0, 350.0, 10000.0 ]
binning_njets = [ -0.5, 0.5, 1.5, 2.5, 3.5, 100.0 ]

shape_ptjet   = [ 78.53302819,  20.52882251,  13.40260032,   3.79521629, 4.44346963, 1.62691503]
shape_yh      = [ 8.88641092, 8.86086463, 17.34621208, 16.38134718, 15.14497828, 51.87741436 ]
shape_pth     = [ 27.46607434,  29.44961617,  19.6934684,   26.59903008,  10.46543842, 6.31602595,   2.0582112, 0.28218741 ]
shape_njets   = [ 78.53302819, 30.8512632, 9.1331588, 2.41446376, 1.39813801 ]

# Normalize to 1.0
shape_njets = [ i/sum(shape_njets) for i in shape_njets ]
shape_pth   = [ i/sum(shape_pth) for i in shape_pth ]
shape_ptjet = [ i/sum(shape_ptjet) for i in shape_ptjet ]
shape_yh    = [ i/sum(shape_yh) for i in shape_yh ]

obs_pth   = Observable( name = 'pth', title = 'p_{T}^{H}', shape = shape_pth, binning = binning_pth )
obs_njets = Observable( name = 'njets', title = 'N_{jets}', shape = shape_njets, binning = binning_njets )
obs_yh    = Observable( name = 'yh', title = '|y_{H}|', shape = shape_yh, binning = binning_yh, lastBinIsOverflow=False )
obs_ptjet = Observable( name = 'ptjet', title = 'p_{T}^{jet}', shape = shape_ptjet, binning = binning_ptjet )

hzz_binMerging_pth   = [ 0, 1, [2,3], [4,5], [6,7] ]
hzz_binMerging_ptjet = [ 0, 1, 2, [3,4,5] ]
hzz_binMerging_njets = [ 0, 1, 2, [3,4] ]
obs_pth_hzzBinning = obs_pth.mergeBins(hzz_binMerging_pth)
obs_ptjet_hzzBinning = obs_ptjet.mergeBins(hzz_binMerging_ptjet)
obs_njets_hzzBinning = obs_njets.mergeBins(hzz_binMerging_njets)



binning_pth_hbb = [ 350., 600., 10000. ]
shape_pth_hbb   = [ 0.5*shape_pth[-1], 0.5*shape_pth[-1] ]
obs_pth_hbbBinning = Observable( name = 'pth', title = 'p_{T}^{H}', shape = shape_pth_hbb, binning = binning_pth_hbb )

binning_pth_combWithHbb = [ 0.0, 15.0, 30.0, 45.0, 85.0, 125.0, 200.0, 350.0, 600.0, 10000.0 ]
shape_pth_combWithHbb   = shape_pth[:-1] + [ 0.9*shape_pth[-1], 0.1*shape_pth[-1] ]
obs_pth_combWithHbbBinning = Observable( name = 'pth', title = 'p_{T}^{H}', shape = shape_pth_combWithHbb, binning = binning_pth_combWithHbb )