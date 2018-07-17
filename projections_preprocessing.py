#!/usr/bin/env python
"""
Thomas Klijnsma
"""

#____________________________________________________________________
# Imports

import logging
import os, sys, re, copy
from OptionHandler import flag_as_option
import differentials
import differentialutils
import LatestBinning
import LatestPaths

# differentials.combine.preprocessing.Datacard.outdir = 'projections/cards_{0}'.format(differentials.core.datestr())

#____________________________________________________________________
# Notes:
# Remember to let CMS_xH_incxs scale with lumi

@flag_as_option
def print_nuisances(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    card = LatestPaths.card.pth_smH[decay_channel]
    datacard = differentials.combine.preprocessing.Datacard(card)
    datacard.print_nuisances()

@flag_as_option
def print_nuisances_caterina(args):
    card = 'projections/caterina_Jun12/comb_2017_ggHbb.txt'
    datacard = differentials.combine.preprocessing.Datacard(card)
    datacard.print_nuisances()

@flag_as_option
def print_nuisances_ftr_hgg(args):
    card = 'projections/ftr_hgg/simple_datacards/hgg_card_13TeV_projection.simple.txt'
    datacard = differentials.combine.preprocessing.Datacard(card)
    datacard.print_nuisances()

@flag_as_option
def create_projection_datacards(args):
    create_hgg_projection_datacard(args)
    create_hbb_projection_datacard(args)
    create_hzz_projection_datacard(args)

@flag_as_option
def create_combWithHbb_projection_datacard(args):
    differentials.combine.preprocessing.combine_cards(
        os.path.join('suppliedInput/projection_combWithHbb'),
        [ 'hgg', 'suppliedInput/fromVittorio/pT_newBins_Feb28/projection_hgg_Jun28.txt' ],
        [ 'hzz', 'suppliedInput/fromJavier/bernstein_r7428/projection_hbb_Jun28.txt' ],
        [ 'hbb', 'suppliedInput/fromDavid/PTH_Jan24_newBinning/smH/projection_hzz_Jun28.txt' ],
        )

@flag_as_option
def create_kbkc_projection_datacard(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    input_card = LatestPaths.card.yukawa[decay_channel]
    datacard = differentials.combine.preprocessing.Datacard(input_card)
    datacard.add_lumiscale_rateparam()
    datacard.out('projection_yukawa_{0}'.format(decay_channel))

@flag_as_option
def create_ktcgkb_projection_datacard(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    input_card = LatestPaths.card.pth_ggH[decay_channel]
    datacard = differentials.combine.preprocessing.Datacard(input_card)
    datacard.add_lumiscale_rateparam()
    datacard.out('projection_ktcgkb_{0}'.format(decay_channel))

@flag_as_option
def create_hgg_projection_datacard(args):
    card = LatestPaths.card.pth_smH.hgg
    datacard = differentials.combine.preprocessing.Datacard(card)

    pho_energy = [
        'CMS_hgg_nuisance_deltafracright',
        'CMS_hgg_nuisance_HighR9EE_13TeVscale',
        'CMS_hgg_nuisance_LowR9EE_13TeVscale',
        'CMS_hgg_nuisance_HighR9EB_13TeVscale',
        'CMS_hgg_nuisance_LowR9EB_13TeVscale',
        'CMS_hgg_nuisance_HighR9EE_13TeVsmear',
        'CMS_hgg_nuisance_LowR9EE_13TeVsmear',
        'CMS_hgg_nuisance_HighR9EB_13TeVsmear',
        'CMS_hgg_nuisance_LowR9EB_13TeVsmear',
        'CMS_hgg_nuisance_MaterialCentral_scale',
        'CMS_hgg_nuisance_MaterialForward_scale',
        'CMS_hgg_nuisance_NonLinearity_13TeVscale',
        'CMS_hgg_nuisance_Geant4_13TeVscale',
        'CMS_hgg_nuisance_Gain1EB_13TeVscale',
        'CMS_hgg_nuisance_Gain6EB_13TeVscale',
        'CMS_hgg_nuisance_LightColl_13TeVscale',
        'CMS_hgg_nuisance_Absolute_13TeVscale',
        ]
    datacard.add_group('pho_energy', pho_energy)

    scale_with_lumi = [
        'CMS_hgg_JER',
        'CMS_hgg_TriggerWeight',
        'CMS_hgg_phoIdMva',
        'CMS_hgg_LooseMvaSF',
        'CMS_hgg_SigmaEOverEShift',
        'CMS_hgg_JEC',
        'CMS_hgg_PreselSF',
        'CMS_hgg_electronVetoSF',
        # 'CMS_hgg_tth_parton_shower',
        # 'CMS_hgg_tth_gluon_splitting',
        # 'CMS_hgg_tth_mc_low_stat',
        'lumi_13TeV',
        'CMS_hgg_PUJIDShift',
        # 'CMS_hgg_eff_e',
        ]
    datacard.add_group('scale_with_lumi', scale_with_lumi)

    datacard.add_lumiscale_rateparam()
    datacard.out('projection_hgg')


@flag_as_option
def create_hbb_projection_datacard(args):
    card = LatestPaths.card.pth_smH.hbb
    datacard = differentials.combine.preprocessing.Datacard(card)

    mcstat_group = [
        'hqq125GenpT2failcat1mcstat',
        'hqq125GenpT2passcat1mcstat',
        'hqq125GenpT3failcat1mcstat',
        'hqq125GenpT3failcat2mcstat',
        'hqq125GenpT3passcat1mcstat',
        'hqq125GenpT3passcat2mcstat',
        'hqq125GenpT4failcat1mcstat',
        'hqq125GenpT4failcat2mcstat',
        'hqq125GenpT4passcat1mcstat',
        'hqq125GenpT4passcat2mcstat',
        'tqqfailcat1mcstat',
        'tqqfailcat2mcstat',
        'tqqpasscat1mcstat',
        'tqqpasscat2mcstat',
        'wqqfailcat1mcstat',
        'wqqfailcat2mcstat',
        'wqqpasscat1mcstat',
        'wqqpasscat2mcstat',
        'xhqq125GenpT1failcat1mcstat',
        'xhqq125GenpT1failcat2mcstat',
        'xhqq125GenpT1passcat1mcstat',
        'xhqq125GenpT1passcat2mcstat',
        'xhqq125GenpT2failcat1mcstat',
        'xhqq125GenpT2failcat2mcstat',
        'xhqq125GenpT2passcat1mcstat',
        'xhqq125GenpT2passcat2mcstat',
        'xhqq125GenpT3failcat1mcstat',
        'xhqq125GenpT3failcat2mcstat',
        'xhqq125GenpT3passcat1mcstat',
        'xhqq125GenpT3passcat2mcstat',
        'xhqq125GenpT4failcat1mcstat',
        'xhqq125GenpT4failcat2mcstat',
        'xhqq125GenpT4passcat1mcstat',
        'xhqq125GenpT4passcat2mcstat',
        'zqqfailcat1mcstat',
        'zqqfailcat2mcstat',
        'zqqpasscat1mcstat',
        'zqqpasscat2mcstat',
        ]
    # Group is already in the datacard!
    # datacard.add_group('mcstat', mcstat_group)

    theory_group = [
        'hqq125pt',
        'hqq125ptShape',
        ]
    datacard.add_group('theory', theory_group)

    tagging_group = [ 
        'bbeff',
        'scale',
        'scalept',
        'smear',
        'veff',
        ]
    datacard.add_group('tagging', tagging_group)

    ttbarcr_group = [
        'muid',
        'muiso',
        'mutrigger',
        'qcdpassmuonCRmcstat',
        'qcdfailmuonCRmcstat',
        'stqqfailmuonCRmcstat',
        'stqqpassmuonCRmcstat',
        'tqqfailmuonCRmcstat',
        'tqqpassmuonCRmcstat',
        'zqqfailmuonCRmcstat',
        'hqq125failmuonCRmcstat',
        'zllfailmuonCRmcstat',
        'wlnupassmuonCRmcstat',
        'wlnufailmuonCRmcstat',
        'vvqqfailmuonCRmcstat ',
        ]
    datacard.add_group('ttbarcr', ttbarcr_group)

    datacard.add_lumiscale_rateparam()
    datacard.out('projection_hbb')


@flag_as_option
def create_hzz_projection_datacard(args):
    card = LatestPaths.card.pth_smH.hzz
    datacard = differentials.combine.preprocessing.Datacard(card)
    datacard.add_lumiscale_rateparam()
    datacard.out('projection_hzz')


