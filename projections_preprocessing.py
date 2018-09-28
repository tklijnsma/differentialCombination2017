#!/usr/bin/env python
"""
Thomas Klijnsma
"""

#____________________________________________________________________
# Imports

import logging
import os, sys, re, copy, glob, json
from OptionHandler import flag_as_option
import differentials
import differentialutils
import LatestBinning
import LatestPaths

import projections_paths as paths

# differentials.combine.preprocessing.Datacard.outdir = 'projections/cards_{0}'.format(differentials.core.datestr())

#____________________________________________________________________
pt_cards = differentials.core.AttrDict()
pt_cards.hgg         = paths.pt_cards.hgg
pt_cards.hzz         = paths.pt_cards.hzz
pt_cards.hbb         = paths.pt_cards.hbb
pt_cards.combWithHbb = paths.pt_cards.combWithHbb

@flag_as_option
def list_nuisances_s2grouping(args):
    card_comb = differentials.combine.preprocessing.Datacard(pt_cards.combWithHbb)
    card_hgg  = differentials.combine.preprocessing.Datacard(pt_cards.hgg)
    card_hzz  = differentials.combine.preprocessing.Datacard(pt_cards.hzz)
    card_hbb  = differentials.combine.preprocessing.Datacard(pt_cards.hbb)

    nuis_comb = card_comb.get_nuisance_names()
    nuis_hgg  = card_hgg.get_nuisance_names()
    nuis_hzz  = card_hzz.get_nuisance_names()
    nuis_hbb  = card_hbb.get_nuisance_names()

    nuis_name_length = max(map(len, nuis_comb))
    group_name_length = max(map(len, card_comb.nuis_group_dict.keys()))

    print '{0:{nuis_name_length}} {1:{group_name_length}}  {2}'.format(
        'Name', 'Group', 'Decay channels',
        nuis_name_length = nuis_name_length,
        group_name_length = group_name_length,
        )

    for nuis in card_comb.nuisances:

        is_in = []
        if nuis.name in nuis_hgg: is_in.append('hgg')
        if nuis.name in nuis_hzz: is_in.append('hzz')
        if nuis.name in nuis_hbb: is_in.append('hbb')

        group = nuis.group if not(nuis.group is None) else 'pOther'

        r = '{0:{nuis_name_length}} {1:{group_name_length}}  {2}'.format(
            nuis.name,
            group,
            ', '.join(is_in),
            nuis_name_length = nuis_name_length,
            group_name_length = group_name_length,            
            )
        print r


@flag_as_option
def create_projection_datacard_s2grouping(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    card = LatestPaths.card.pth_smH[decay_channel]
    datacard = differentials.combine.preprocessing.Datacard(card)
    datacard.delete_current_nuisance_groups()
    nuisgroupfinder = NuisanceGroupFinder()
    datacard.lines.extend(
        nuisgroupfinder.compile_datacard_lines(datacard.get_nuisance_names())
        )
    datacard.add_lumiscale_rateparam()
    datacard.out('projection_{0}_s2groups'.format(decay_channel))

@flag_as_option
def create_projection_datacards_s2grouping(args):
    for dc in [ 'hbb', 'hgg', 'hzz', 'combWithHbb' ]:
        args = differentialutils.set_one_decay_channel(args, dc)
        create_projection_datacard_s2grouping(args)

@flag_as_option
def create_kbkc_projection_datacard_s2grouping(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    input_card = LatestPaths.card.yukawa[decay_channel]
    datacard = differentials.combine.preprocessing.Datacard(input_card)
    NuisanceGroupFinder().process_datacard(datacard)
    datacard.add_lumiscale_rateparam()
    datacard.out('projection_yukawa_{0}_s2groups'.format(decay_channel))

@flag_as_option
def create_ktcgkb_projection_datacard_s2grouping(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    input_card = LatestPaths.card.pth_ggH[decay_channel]
    datacard = differentials.combine.preprocessing.Datacard(input_card)
    NuisanceGroupFinder().process_datacard(datacard)
    datacard.add_lumiscale_rateparam()
    datacard.out('projection_ktcgkb_{0}_s2groups'.format(decay_channel))

@flag_as_option
def create_coupling_projection_datacards_s2grouping(args):
    for dc in [ 'hzz', 'hgg', 'combWithHbb' ]:
        newargs = differentialutils.set_one_decay_channel(args, dc)
        create_ktcgkb_projection_datacard_s2grouping(newargs)

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
def create_hbb_projection_datacard_befores2(args):
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


@flag_as_option
def test_nuisgroups(args):
    nuisgroupfinder = NuisanceGroupFinder()
    nuisgroupfinder.print_groups()


class NuisanceGroupFinder(object):
    """docstring for NuisanceGroupFinder"""
    grouping_path = 'projections/andrew/np-groups-July17'
    unassigned_means_pother = False
    return_value_for_not_found = 'unassigned'

    def __init__(self):
        super(NuisanceGroupFinder, self).__init__()
        self.grouping = {}
        self.read_grouping()

        # hbb card
        self.add_grouping({'mcstat' : [
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
            # Originally in ttbarcr group
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
            # Originally no group for:
            'vvqqfailmuonCRmcstat',
            'xhqq125failmuonCRmcstat',
            'xhqq125passmuonCRmcstat',
            ]})

        # The following are not in a hbb group, and also not in Andrew's groups
        # Looks very much like it belongs to the bkgTheory group
        # --> Actually these are rateParams, and t2ws crashes if you put them in
        #     a group. Instead make them unassigned.
        self.add_grouping({
            'unassigned' : [
                'tqqfailcat1norm',
                'tqqfailcat2norm',
                'tqqfailmuonCRnorm',
                'tqqpasscat1norm',
                'tqqpasscat2norm',
                'tqqpassmuonCRnorm',
                ]
            })

        # Fixing hgg
        self.add_grouping({
            'pLumi' : [
                'lumi_13TeV',
                ],
            'pResJ' : [
                'CMS_hgg_JER',
                ],
            'pScaleJ' : [
                'CMS_hgg_JEC',
                ],
            'pScaleJPileup' : [
                'CMS_hgg_PUJIDShift',
                ],
            'unassigned' : [
                # Nuisances in Hgg datacard that are not grouped by Andrews grouping
                'CMS_hgg_nuisance_MaterialCentral_scale',
                'CMS_hgg_nuisance_LightColl_13TeVscale',
                'CMS_hgg_nuisance_Absolute_13TeVscale',
                'CMS_hgg_nuisance_HighR9EE_13TeVsmear',
                'CMS_hgg_nuisance_LowR9EE_13TeVsmear',
                'CMS_hgg_nuisance_HighR9EB_13TeVsmear',
                'CMS_hgg_nuisance_LowR9EB_13TeVsmear',
                'CMS_hgg_LooseMvaSF',
                # These ones are simply empty in the datacard! That's why they crash...
                'CMS_hgg_tth_parton_shower',
                'CMS_hgg_tth_gluon_splitting',
                'CMS_hgg_tth_mc_low_stat',
                'CMS_hgg_eff_e',
                ]
            })

        # Fixing hzz
        self.add_grouping({
            'pMuonID' : [
                'CMS_eff_m',
                'CMS_zz4l_mean_m_sig',
                'CMS_zz4l_sigma_m_sig',
                ],
            'bkgTheory' : [ 'QCDscale_ggVV', 'kfactor_ggzz', 'norm_nonResH' ],
            'unassigned' : [
                'CMS_eff_e',
                'CMS_hzz2e2mu_Zjets',
                'CMS_hzz4e_Zjets',
                'CMS_hzz4mu_Zjets',
                'CMS_zjets_bkgdcompo',
                # No idea for these:
                'CMS_zz4l_n_sig_1_8',
                'CMS_zz4l_n_sig_3_8',
                'CMS_zz4l_mean_e_sig',
                'CMS_zz4l_sigma_e_sig',
                'CMS_zz4l_n_sig_2_8',
                ],
            })

        # The nuis due to uncertainty on xH should scale with lumi
        self.add_grouping({'pOther' : ['CMS_xH_incxs']})


    def read_grouping(self):
        json_files = glob.glob(os.path.join(self.grouping_path, '*.json'))
        for json_file in json_files:
            with open(json_file, 'rb') as json_fp:
                grouping = json.load(json_fp)
                self.add_grouping(grouping, skip_deletion=True)

    def add_grouping(self, grouping, skip_deletion=False):
        for key, values in grouping.iteritems():
            if not key in self.grouping: self.grouping[key] = []
            if not skip_deletion:
                for val in values: self.delete_from_already_grouped(val)
            self.grouping[key].extend(values)

    def delete_from_already_grouped(self, nuis):
        group = self.check_already_grouped(nuis)
        if group is None: return
        logging.info('Removing {0} from group {1}'.format(nuis, group))
        self.grouping[group].pop(self.grouping[group].index(nuis))

    def check_already_grouped(self, nuis):
        for key, values in self.grouping.iteritems():
            if nuis in values:
                return key
        else:
            return None

    def lookup(self, nuis):
        if nuis.startswith('pdfindex_'): return None
        for key, values in self.grouping.iteritems():
            if nuis in values:
                if self.unassigned_means_pother and key == 'unassigned':
                    logging.info('Found \'unassigned\' for \'{0}\'; returning pOther'.format(nuis))
                    return 'pOther'
                return key
        else:
            logging.info('Found no group for \'{0}\'; returning {1}'.format(nuis, self.return_value_for_not_found))
            return self.return_value_for_not_found


    def print_groups(self):
        for key, values in self.grouping.iteritems():
            print 'group {0}'.format(key)
            print '    ' + '\n    '.join(values)


    def compile_datacard_lines(self, nuisances, verbose=True):
        # Gather all the groups that should be created
        groups = {}
        for nuis in nuisances:
            group = self.lookup(nuis)
            if group is None or group == 'unassigned': continue
            if not group in groups.keys(): groups[group] = []
            groups[group].append(nuis)

        # Create lines
        lines = []
        for groupname, nuisnames in groups.iteritems():
            line = '{0} group = {1}'.format(
                groupname, ' '.join(nuisnames)
                )
            lines.append(line)

        if verbose:
            logging.info(
                'Compiled the following new group lines:\n{0}'
                .format('\n'.join(lines))
                )

        return lines

    def process_datacard(self, datacard):
        datacard.delete_current_nuisance_groups()
        datacard.lines.extend(
            self.compile_datacard_lines(datacard.get_nuisance_names())
            )

