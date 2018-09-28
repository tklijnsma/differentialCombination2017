import sys
from differentials.observable import Observable
from collections import namedtuple
from copy import deepcopy
from math import sqrt
import numpy

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

def normalize(l, normalization=1.0):
    s = float(sum(l))
    return [ e/s*normalization for e in l ]

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
YR4_totalXS_uncertainty = 2.5

# YR4_totalXS = 55.70628722 # pb
# YR4_ggHXS   = 48.58 # pb

SM_BR_hgg = 2.270E-03
SM_BR_hzz = 0.02641
SM_ratio_of_BRs = SM_BR_hgg / SM_BR_hzz

SM_ratio_of_BRs_unc_up = +0.00223
SM_ratio_of_BRs_unc_down = -0.00220

Observable.YR4_totalXS = YR4_totalXS


########################################
# Binning and SM differentials
########################################

# ----------------
# Observables other than pth
binning_ptjet = [ -1000.0, 30.0, 55.0, 95.0, 120.0, 200.0, 13000.0 ]
binning_yh    = [ 0.0, 0.15, 0.3, 0.6, 0.9, 1.2, 2.5 ]
binning_njets = [ 0., 1.0, 2.0, 3.0, 4.0, 100.0 ]

shape_ptjet   = normalize([ 78.53302819,  20.52882251,  13.40260032,   3.79521629, 4.44346963, 1.62691503])
shape_yh      = normalize([ 8.88641092, 8.86086463, 17.34621208, 16.38134718, 15.14497828, 51.87741436 ])
shape_njets   = normalize([ 78.53302819, 30.8512632, 9.1331588, 2.41446376, 1.39813801 ])

obs_njets = Observable( name = 'njets', title = 'N_{jets}', shape = shape_njets, binning = binning_njets )
obs_yh    = Observable( name = 'yh', title = '|y_{H}|', shape = shape_yh, binning = binning_yh, lastBinIsOverflow=False )
obs_ptjet = Observable( name = 'ptjet', title = 'p_{T}^{jet}', shape = shape_ptjet, binning = binning_ptjet )

hzz_binMerging_ptjet = [ 0, 1, 2, [3,4,5] ]
hzz_binMerging_njets = [ 0, 1, 2, [3,4] ]
obs_ptjet_hzzBinning = obs_ptjet.mergeBins(hzz_binMerging_ptjet)
obs_njets_hzzBinning = obs_njets.mergeBins(hzz_binMerging_njets)


# ----------------
# pth
binning_pth   = [ 0.0, 15.0, 30.0, 45.0, 80.0, 120.0, 200.0, 350.0, 600., 10000.0 ]

shape_pth_smH = [
    2.74660743e+01, 2.94496162e+01, 1.96934684e+01, 2.44876793e+01,
    1.17284547e+01, 7.16436043e+00, 2.05821120e+00, 2.63714145e-01,
    1.84732619e-02]
shape_pth_ggH = [
    2.70051430e+01, 2.81966864e+01, 1.79667435e+01, 2.03254009e+01,
    8.39956789e+00, 4.53454528e+00, 1.19397255e+00, 1.39795879e-01,
    7.51079245e-03
    ]
shape_pth_xH  = [ s_tot-s_ggH for s_tot, s_ggH in zip(shape_pth_smH, shape_pth_ggH) ]

shape_pth_smH = normalize(shape_pth_smH)
shape_pth_ggH = normalize(shape_pth_ggH)
shape_pth_xH  = normalize(shape_pth_xH)

obs_pth_smH = Observable( name = 'pth_smH', title = 'p_{T}^{H}', shape = shape_pth_smH, binning = binning_pth )
obs_pth_ggH = Observable( name = 'pth_ggH', title = 'p_{T}^{H} (ggH)', shape = shape_pth_ggH, binning = binning_pth )
obs_pth_ggH.YR4_totalXS = YR4_ggF_n3lo
obs_pth_xH = Observable( name = 'pth_xH', title = 'p_{T}^{H} (xH)', shape = shape_pth_xH, binning = binning_pth )
obs_pth_xH.YR4_totalXS = YR4_xH

# xH individual shapes
ttH_npz = numpy.load('suppliedInput/fromVittorio/FullPhaseSpace_spectrumNNLOPS_ttH_Pt.npz')
shape_pth_ttH = normalize(ttH_npz['spectrum'])
obs_pth_ttH = Observable( name = 'pth_ttH', title = 'p_{T}^{H} (ttH)', shape = shape_pth_ttH, binning = binning_pth )
obs_pth_ttH.YR4_totalXS = YR4_ttH

# Different binning schemes for HZZ and Hbb
hzz_binMerging_pth   = [ 0, 1, [2,3], [4,5], [6,7,8] ]
obs_pth_smH_hzzBinning = obs_pth_smH.mergeBins(hzz_binMerging_pth)
obs_pth_ggH_hzzBinning = obs_pth_ggH.mergeBins(hzz_binMerging_pth)

obs_pth_smH_hbbBinning = deepcopy(obs_pth_smH)
obs_pth_ggH_hbbBinning = deepcopy(obs_pth_ggH)
# !!!: Dropping up to 350 now!! 200-350 no longer included; careful when multiplying this with a scan!
for i in xrange(7):
    obs_pth_smH_hbbBinning.drop_first_bin()
    obs_pth_ggH_hbbBinning.drop_first_bin()


#____________________________________________________________________
# Convenience tuples

obstuple_pth_smH  = AttrDict(
    hgg         = obs_pth_smH,
    hzz         = obs_pth_smH_hzzBinning,
    combination = obs_pth_smH,
    combWithHbb = obs_pth_smH,
    hbb         = obs_pth_smH_hbbBinning
    )
obstuple_pth_ggH  = AttrDict(
    hgg         = obs_pth_ggH,
    hzz         = obs_pth_ggH_hzzBinning,
    combination = obs_pth_ggH,
    combWithHbb = obs_pth_ggH,
    hbb         = obs_pth_ggH_hbbBinning
    )
obstuple_njets    = AttrDict(
    hgg         = obs_njets,
    hzz         = obs_njets_hzzBinning,
    combination = obs_njets,
    combWithHbb = obs_njets,
    )
obstuple_ptjet    = AttrDict(
    hgg         = obs_ptjet,
    hzz         = obs_ptjet_hzzBinning,
    combination = obs_ptjet,
    combWithHbb = obs_ptjet,
    )
obstuple_rapidity = AttrDict(
    hgg         = obs_yh,
    hzz         = obs_yh,
    combination = obs_yh,
    combWithHbb = obs_yh,
    )


#____________________________________________________________________
# xH uncertainty, inclusive

# Uncertainties per mode, all in percentages of total XS
# first one is scale, second PDF, third alpha_s
uncs_VBF             = [ 0.35, 2.1, 0.05 ]
uncs_WH              = [ 0.6, 1.7, 0.9 ]
uncs_ZH              = [ 3.4, 1.3, 0.9 ]
uncs_ttH             = [ 7.5, 3.0, 2.0 ]
uncs_bbH             = [ 22.0 ]
uncs_tH_t_ch         = [ 10.6, 3.5, 1.2 ]
uncs_tH_s_ch         = [ 2.1, 2.2, 0.2 ]
uncs_tH_W_associated = [ 5.8, 6.1, 1.5 ]

unc_squared_per_mode = lambda uncs, xs: sum([ (0.01*unc * xs)**2 for unc in uncs ])

tot_unc_squared = 0.0
tot_unc_squared += unc_squared_per_mode(uncs_VBF, YR4_VBF)
tot_unc_squared += unc_squared_per_mode(uncs_WH, YR4_WH)
tot_unc_squared += unc_squared_per_mode(uncs_ZH, YR4_ZH)
tot_unc_squared += unc_squared_per_mode(uncs_ttH, YR4_ttH)
tot_unc_squared += unc_squared_per_mode(uncs_bbH, YR4_bbH)
tot_unc_squared += unc_squared_per_mode(uncs_tH_t_ch, YR4_tH_t_ch)
tot_unc_squared += unc_squared_per_mode(uncs_tH_s_ch, YR4_tH_s_ch)
tot_unc_squared += unc_squared_per_mode(uncs_tH_W_associated, YR4_tH_W_associated)

xH_unc_inclusive = sqrt(tot_unc_squared)
xH_unc_inclusive_fraction = xH_unc_inclusive / YR4_xH

uncs_ggF = [ 5.65, 3.2 ]
tot_unc_squared += unc_squared_per_mode(uncs_ggF, YR4_ggF_n3lo)
smH_unc_inclusive = sqrt(tot_unc_squared)
smH_unc_inclusive_fraction = smH_unc_inclusive / YR4_totalXS

#____________________________________________________________________
# smH uncertainties

add_quad = lambda *args: sqrt(sum([ x**2 for x in args ]))
def add_unc(obs, unc_file):
    xs = obs.crosssection()
    shape_unc_perc = list(numpy.load(unc_file)['uncetainty'])
    # Add in quadrature with incl xs uncertainty from YR4
    # Do this all in 'percent-space'
    unc_perc = [ add_quad(err, smH_unc_inclusive_fraction) for err in shape_unc_perc ]
    obs.unc_fraction = unc_perc

add_unc(obs_yh, 'suppliedInput/fromVittorio/sm_uncertainties/uncertaintyFullPhaseSpaceCombination_AbsRapidity.npz')
add_unc(obs_ptjet, 'suppliedInput/fromVittorio/sm_uncertainties/uncertaintyFullPhaseSpaceCombination_Jet2p5Pt0.npz')
add_unc(obs_njets, 'suppliedInput/fromVittorio/sm_uncertainties/uncertaintyFullPhaseSpaceCombination_Njets2p5.npz')
add_unc(obs_pth_smH, 'suppliedInput/fromVittorio/sm_uncertainties/uncertaintyFullPhaseSpaceCombination_Pt.npz')
add_unc(obs_pth_ggH, 'suppliedInput/fromVittorio/sm_uncertainties/uncertaintyFullPhaseSpaceCombination_ggH_Pt.npz')


if __name__ == "__main__":
    # print 'Unc from VBF:             ', sqrt(unc_squared_per_mode(uncs_VBF, YR4_VBF))
    # print 'Unc from WH:              ', sqrt(unc_squared_per_mode(uncs_WH, YR4_WH))
    # print 'Unc from ZH:              ', sqrt(unc_squared_per_mode(uncs_ZH, YR4_ZH))
    # print 'Unc from ttH:             ', sqrt(unc_squared_per_mode(uncs_ttH, YR4_ttH))
    # print 'Unc from bbH:             ', sqrt(unc_squared_per_mode(uncs_bbH, YR4_bbH))
    # print 'Unc from tH_t_ch:         ', sqrt(unc_squared_per_mode(uncs_tH_t_ch, YR4_tH_t_ch))
    # print 'Unc from tH_s_ch:         ', sqrt(unc_squared_per_mode(uncs_tH_s_ch, YR4_tH_s_ch))
    # print 'Unc from tH_W_associated: ', sqrt(unc_squared_per_mode(uncs_tH_W_associated, YR4_tH_W_associated))
    # print
    print 'YR4_xH:                   ', YR4_xH
    print 'xH_unc_inclusive:         ', xH_unc_inclusive
    print 'xH_unc_inclusive_fraction:', xH_unc_inclusive_fraction

    print 'smH_unc_inclusive:         ', smH_unc_inclusive
    print 'smH_unc_inclusive_fraction:', smH_unc_inclusive_fraction

    # print list(numpy.load('suppliedInput/fromVittorio/sm_uncertainties/uncertaintyFullPhaseSpaceCombination_AbsRapidity.npz')['uncetainty'])
    # print obs_yh.unc_fraction