#!/usr/bin/env python
"""
Thomas Klijnsma
"""

import os, re, copy
from os.path import *

import logging
import differentials
import differentials.core as core
from OptionHandler import flag_as_option, flag_as_parser_options
import LatestPaths
import LatestBinning
import differentialutils

#____________________________________________________________________
# Text changes to datacards

# Fixes to PTH workspaces
@flag_as_option
def rename_processes_hgg_pth_ggH(args):
    rename_processes_hgg_pth(
        'suppliedInput/fromVittorio/pT_newBins_Feb28/Datacard_13TeV_differential_PtGghPlusHxNNLOPS_newBins.txt'
        )
@flag_as_option
def renumber_processes_hzz_pth_ggH(args):
    renumber_processes_hzz_pth(
        'suppliedInput/fromDavid/PTH_Jan24_newBinning/ggH/hzz4l_comb_13TeV_xs.txt'
        )
@flag_as_option
def rename_processes_hgg_pth_smH(args):
    rename_processes_hgg_pth(
        'suppliedInput/fromVittorio/pT_newBins_Mar13_smH/Datacard_13TeV_differential_PtNNLOPS_newBins.txt'
        )
@flag_as_option
def renumber_processes_hzz_pth_smH(args):
    renumber_processes_hzz_pth(
        'suppliedInput/fromDavid/PTH_Jan24_newBinning/smH/hzz4l_comb_13TeV_xs.txt'
        )

# Hgg differentials: njets, ptjet, rapidity
@flag_as_option
def rename_processes_hgg_njets(args):
    rename_processes_hgg_differentials(
        'suppliedInput/fromVittorio/differential_Njets2p5NNLOPS_Nov10/Datacard_13TeV_differential_Njets2p5NNLOPS.txt',
        global_replace = [( 'CMS_hgg_JEC', 'CMS_scale_j' )]
        )
@flag_as_option
def rename_processes_hgg_ptjet(args):
    rename_processes_hgg_differentials(
        'suppliedInput/fromVittorio/differential_Jet2p5Pt0NNLOPS_newBins_Nov28/Datacard_13TeV_differential_Jet2p5Pt0NNLOPS_newBins.txt'
        )
@flag_as_option
def rename_processes_hgg_rapidity(args):
    rename_processes_hgg_differentials(
        'suppliedInput/fromVittorio/differential_AbsRapidityNNLOPS_newBins_Nov28/Datacard_13TeV_differential_AbsRapidityNNLOPS_newBins_combination.txt'
        )

# Fixes to INC workspaces
@flag_as_option
def rename_processes_hgg_INC(args):
    rename_processes_hgg_differentials(
        'suppliedInput/fromVittorio/inclusive_Nov27/Datacard_13TeV_differential_InclusiveNNLOPS.txt',
        )

# Adding an xH nuisance parameter to the combined ggH datacard
@flag_as_option
def add_xH_nuisance_parameter_to_ggH_card(args):
    # args = differentialutils.set_one_decay_channel(args, 'combWithHbb')
    # add_xH_nuisance_parameter('suppliedInput/combWithHbb_pth_ggH_Mar02.txt')
    add_xH_nuisance_parameter('suppliedInput/Yukawa_combination_pth_ggH_Mar08.txt', 24, 20, 43)

@flag_as_option
def add_xH_nuisance_parameter_to_ggH_card(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    card = LatestPaths.card.pth_ggH_noXHunc[decay_channel]
    add_xH_nuisance_parameter(card, 24, 20, 43)


#____________________________________________________________________
# combineCards

@flag_as_option
def combine_all_cards(args):
    args = differentialutils.set_one_decay_channel(args, 'combination')
    # combine_pth_ggH(args)
    # combine_pth_smH(args)
    combine_njets(args)
    combine_ptjet(args)
    combine_rapidity(args)
    # args = differentialutils.set_one_decay_channel(args, 'combWithHbb')
    # combine_pth_ggH(args)
    # combine_pth_smH(args)

@flag_as_option
def combine_pth_ggH(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    out_card = 'suppliedInput/{0}_pth_ggH_{1}.txt'.format(
            decay_channel, core.datestr()
            )
    cmd = []
    if args.combWithHbb or args.combination or args.hgg:
        cmd.append('hgg=' + LatestPaths.card.pth_ggH.hgg)
    if args.combWithHbb or args.hbb:
        cmd.append('hbb=' + LatestPaths.card.pth_ggH.hbb)
    if args.combWithHbb or args.combination or args.hzz:
        cmd.append('hzz=' + LatestPaths.card.pth_ggH.hzz)
    combine_cards(out_card, *cmd)

@flag_as_option
def combine_pth_smH(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    out_card = 'suppliedInput/{0}_pth_smH_{1}.txt'.format(
            decay_channel, core.datestr()
            )
    cmd = []
    if args.combWithHbb or args.combination or args.hgg:
        cmd.append('hgg=' + LatestPaths.card.pth_ggH.hgg) # use ggH, but use special scaling in T2WS
    if args.combWithHbb or args.hbb:
        cmd.append('hbb=' + LatestPaths.card.pth_ggH.hbb) # use ggH, but use special scaling in T2WS
    if args.combWithHbb or args.combination or args.hzz:
        cmd.append('hzz=' + LatestPaths.card.pth_smH.hzz)
    combine_cards(out_card, *cmd)

@flag_as_option
def combine_njets(args):
    args = differentialutils.set_one_decay_channel(args, 'combination')
    combine_cards_for_observable(args, 'njets')

@flag_as_option
def combine_ptjet(args):
    args = differentialutils.set_one_decay_channel(args, 'combination')
    combine_cards_for_observable(args, 'ptjet')

@flag_as_option
def combine_rapidity(args):
    args = differentialutils.set_one_decay_channel(args, 'combination')
    combine_cards_for_observable(args, 'rapidity')

@flag_as_option
def combine_inclusive(args):
    args = differentialutils.set_one_decay_channel(args, 'combination')
    combine_cards_for_observable(args, 'inclusive')

def combine_cards_for_observable(args, obsname):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    out_card = 'suppliedInput/{0}_{1}_{2}.txt'.format(
            decay_channel, obsname, core.datestr()
            )
    cmd = []
    if args.combWithHbb or args.combination or args.hgg:
        cmd.append('hgg=' + LatestPaths.card[obsname].hgg)
    if args.combWithHbb or args.combination or args.hzz:
        cmd.append('hzz=' + LatestPaths.card[obsname].hzz)
    if args.combWithHbb or args.hbb:
        cmd.append('hbb=' + LatestPaths.card[obsname].hbb)
    combine_cards(out_card, *cmd)


#____________________________________________________________________
# Yukawa

@flag_as_option
def add_xH_Yukawa(real_args):
    args = copy.deepcopy(real_args)
    for dc in ['hgg', 'hzz', 'combination']:
        if dc == 'combination': continue
        differentialutils.set_one_decay_channel(args, dc)
        card = LatestPaths.card.yukawa_noXHunc[dc]
        add_xH_nuisance_parameter(card, 24, 20, 43)

@flag_as_option
def combine_all_cards_Yukawa(real_args):
    args = copy.deepcopy(real_args)
    for dc in ['hgg', 'hzz', 'combination']:
        differentialutils.set_one_decay_channel(args, dc)
        combine_cards_Yukawa(args)

@flag_as_option
def combine_cards_Yukawa(args):
    hgg_cat_pats = [
        'recoPt_120p0_200p0',
        'recoPt_200p0_350p0',
        'recoPt_350p0_600p0',
        'recoPt_600p0_10000p0',
        ]
    hzz_cat_pats = [
        'hzz_PTH_GT200_cat2e2mu',
        'hzz_PTH_GT200_cat4e',
        'hzz_PTH_GT200_cat4mu',
        ]

    decay_channel = differentialutils.get_decay_channel_tag(args)
    out_card = 'suppliedInput/Yukawa_{0}_pth_ggH_{1}.txt'.format(
            decay_channel, core.datestr()
            )

    cmd = []
    if args.combination or args.hzz:
        cmd.append('hzz=' + LatestPaths.card.pth_ggH.hzz)
        for cat_pat in hzz_cat_pats:
            cmd.append('--xc={0}.*'.format(cat_pat))
    if args.combination or args.hgg:
        cmd.append('hgg=' + LatestPaths.card.pth_ggH.hgg)
        for cat_pat in hgg_cat_pats:
            cmd.append('--xc={0}.*'.format(cat_pat))

    combine_cards(out_card, *cmd)

    # Have to manually remove the pdf indices
    if args.combWithHbb or args.combination or args.hgg:
        drop_pdfindices(out_card, hgg_cat_pats)

#____________________________________________________________________
# Top

@flag_as_option
def combine_cards_Top_noBinsDropped(args):
    decay_channel = differentialutils.get_decay_channel_tag(args)
    out_card = 'suppliedInput/Top_{0}_pth_ggH_{1}_noBinsDropped.txt'.format(
            decay_channel, core.datestr()
            )
    cmd = []
    if args.combWithHbb or args.combination or args.hgg:
        cmd.append('hgg=' + LatestPaths.card.pth_ggH.hgg)
    if args.combWithHbb or args.hbb:
        cmd.append('hbb=' + LatestPaths.card.pth_ggH.hbb)
    if args.combWithHbb or args.combination or args.hzz:
        cmd.append('hzz=' + LatestPaths.card.pth_ggH.hzz)
    combine_cards(out_card, *cmd)

@flag_as_option
def combine_cards_pth_ggH_lastBinDroppedHgg(args):
    out_card = 'suppliedInput/combination_pth_ggH_lastBinDroppedHgg_{0}.txt'.format(core.datestr())
    combine_cards(
        out_card,
        'hgg=' + LatestPaths.card.pth_ggH.hgg,
        'hzz=' + LatestPaths.card.pth_ggH.hzz,
        '--xc=recoPt_600p0_10000p0_.*'
        )
    drop_pdfindices(out_card)
    
@flag_as_option
def combine_cards_pth_ggH_lastBinDroppedHgg_hbb(args):
    out_card = 'suppliedInput/combWithHbb_pth_ggH_lastBinDroppedHgg_{0}.txt'.format(core.datestr())
    combine_cards(
        out_card,
        'hgg=' + LatestPaths.card.pth_ggH.hgg,
        'hzz=' + LatestPaths.card.pth_ggH.hzz,
        'hbb=' + LatestPaths.card.pth_ggH.hbb,
        '--xc=recoPt_600p0_10000p0_.*'
        )
    drop_pdfindices(out_card)

@flag_as_option
def combine_cards_Top_last2BinsDropped(args):
    hgg_cat_pats = [
        'recoPt_350p0_600p0',
        'recoPt_600p0_10000p0',
        ]

    decay_channel = differentialutils.get_decay_channel_tag(args)
    out_card = 'suppliedInput/Top_{0}_pth_ggH_{1}_last2BinsDropped.txt'.format(
            decay_channel, core.datestr()
            )

    cmd = []
    if args.combWithHbb or args.combination or args.hgg:
        cmd.append('hgg=' + LatestPaths.card.pth_ggH.hgg)
        for cat_pat in hgg_cat_pats:
            cmd.append('--xc={0}.*'.format(cat_pat))
    if args.combWithHbb or args.hbb:
        cmd.append('hbb=' + LatestPaths.card.pth_ggH.hbb)
    if args.combWithHbb or args.combination or args.hzz:
        cmd.append('hzz=' + LatestPaths.card.pth_ggH.hzz)

    combine_cards(out_card, *cmd)

    # Have to manually remove the 
    if args.combWithHbb or args.combination or args.hgg:
        drop_pdfindices(out_card, hgg_cat_pats)



#____________________________________________________________________
# Helper functions

def combine_cards(
        output_file,
        *input_list
        ):
    cmd = [ 'combineCards.py' ]
    for datacard in input_list:
        cmd.append(datacard)
    cmd.append( '> {0}'.format( output_file ) )
    core.execute(cmd)

def drop_pdfindices(card_file, category_pats=None):
    if differentials.core.is_testmode():
        return
    with open(card_file, 'r') as card_fp:
        card = card_fp.read()

    if category_pats is None:
        category_pats = ['recoPt_600p0_10000p0']

    lines = []
    for line in card.split('\n'):
        for category_pat in category_pats:
            if re.match(r'pdfindex_.*{0}'.format(category_pat), line):
                logging.debug(
                    'Dropping following line from {0} (matched to {2}):\n{1}'
                    .format(card_file, line, category_pat)
                    )
                break
        else:
            lines.append(line)

    new_card = '\n'.join(lines)
    logging.trace('Datacard after removing lines:\n{0}'.format(new_card))
    logging.info('Writing new card after deleting lines to {0}'.format(card_file))
    if not core.is_testmode():
        with open(card_file, 'w') as card_fp:
            card_fp.write(new_card)


def sub(pat, repl, text, ntimes=3):
    """Does a replacement 3 times to process matches next to each other as well"""
    for i in xrange(ntimes):
        text = re.sub(pat, repl, text)
    return text

def write_to_file(out_file, contents):
    logging.debug('Contents of datacard {0}:\n{1}'.format(out_file, contents))
    logging.info('Opening {0} and dumping contents'.format(out_file))
    if not core.is_testmode():
        with open(out_file, 'w') as out_fp:
            out_fp.write(contents)


def disable_200_350_process_hbb_pth(
        datacard_file,
        rename_OutsideAcceptance = True,
        global_replace = None,
        tag = '_disabled200350'
        ):
    with open(datacard_file, 'r') as datacard_fp:
        datacard = datacard_fp.read()
    out = datacard

    # ggH
    out = sub(
        r'(\W)ggH_PTH_200_350(\W)',
        r'\1DISABLEDggH_PTH_200_350\2',
        out
        )

    # Process simple global replacements
    if not global_replace == None:
        for string, replacement in global_replace:
            out = out.replace(string, replacement)

    # Write to file
    out_file = datacard_file.replace('.txt', '{0}.txt'.format(tag))
    write_to_file(out_file, out)


def make_hbb_pth_smH_fn(
        datacard_file,
        tag = '_smH'
        ):
    with open(datacard_file, 'r') as datacard_fp:
        datacard = datacard_fp.read()
    out = datacard

    # ggH
    # out = sub(
    #     r'(\W)ggH_PTH_200_350(\W)',
    #     r'\1DISABLEDggH_PTH_200_350\2',
    #     out
    #     )

    out = sub(
        r'(\W)[xg]+H_PTH_350_600(\W)',
        r'\1smH_PTH_350_600\2',
        out
        )

    out = sub(
        r'(\W)[xg]+H_PTH_GT600(\W)',
        r'\1smH_PTH_GT600\2',
        out
        )

    # Write to file
    out_file = datacard_file.replace('.txt', '{0}.txt'.format(tag))
    write_to_file(out_file, out)


def rename_processes_hgg_pth(
        datacard_file,
        rename_OutsideAcceptance = True,
        global_replace = None,
        tag = '_renamedProcesses'
        ):
    with open(datacard_file, 'r') as datacard_fp:
        datacard = datacard_fp.read()
    out = datacard

    # ggH
    out = sub(
        r'(\W)gghInsideAcceptance_genPt_350p0_10000p0(\W)',
        r'\1ggH_PTH_GT350\2',
        out
        )
    out = sub(
        r'(\W)gghInsideAcceptance_genPt_600p0_10000p0(\W)',
        r'\1ggH_PTH_GT600\2',
        out
        )
    out = sub(
        r'(\W)gghInsideAcceptance_genPt_(\d+)p0_(\d+)p0(\W)',
        r'\1ggH_PTH_\2_\3\4',
        out
        )
    if rename_OutsideAcceptance:
        out = re.sub(
            r'(\W)gghOutsideAcceptance(\W)',
            r'\1ggH_OutsideAcceptance\2',
            out
            )

    # xH
    out = sub(
        r'(\W)hxInsideAcceptance_genPt_350p0_10000p0(\W)',
        r'\1xH_PTH_GT350\2',
        out
        )
    out = sub(
        r'(\W)hxInsideAcceptance_genPt_600p0_10000p0(\W)',
        r'\1xH_PTH_GT600\2',
        out
        )
    out = sub(
        r'(\W)hxInsideAcceptance_genPt_(\d+)p0_(\d+)p0(\W)',
        r'\1xH_PTH_\2_\3\4',
        out
        )
    if rename_OutsideAcceptance:
        out = re.sub(
            r'(\W)hxOutsideAcceptance(\W)',
            r'\1xH_OutsideAcceptance\2',
            out
            )

    # smH
    out = sub(
        r'(\W)InsideAcceptance_genPt_350p0_10000p0(\W)',
        r'\1smH_PTH_GT350\2',
        out
        )
    out = sub(
        r'(\W)InsideAcceptance_genPt_600p0_10000p0(\W)',
        r'\1smH_PTH_GT600\2',
        out
        )
    out = sub(
        r'(\W)InsideAcceptance_genPt_(\d+)p0_(\d+)p0(\W)',
        r'\1smH_PTH_\2_\3\4',
        out
        )

    # Process simple global replacements
    if not global_replace == None:
        for string, replacement in global_replace:
            out = out.replace(string, replacement)

    # Write to file
    out_file = datacard_file.replace('.txt', '{0}.txt'.format(tag))
    write_to_file(out_file, out)


def renumber_processes_hzz_pth(datacard_file):
    with open( datacard_file, 'r' ) as datacard_fp:
        lines = datacard_fp.readlines()

    first_match = True
    for i_line, line in enumerate(lines):
        if line.startswith('process '):
            if first_match:
                names = line.split()[1:]
                first_match = False
                i_name_line = i_line
            else:
                numbers = line.split()[1:]
                i_number_line = i_line
                break
    else:
        raise RuntimeError(
            'Reached end of file while searching two lines '
            'that start with \'process\'; datacard {0} is faulty.'
            .format(datacard_file)
            )

    signals = []
    bkgs    = []
    for process in list(set(names)):
        if process.startswith('ggH') or process.startswith('xH') or process.startswith('smH'):
            signals.append(process)
        else:
            bkgs.append(process)

    signals.sort( key=process_sorter )
    bkgs.sort( key=process_sorter )

    number_dict = {}
    for iSignal, signal in enumerate( signals ):
        number_dict[signal] = -iSignal
    for iBkg, bkg in enumerate( bkgs ):
        number_dict[bkg] = iBkg + 1

    new_number_line = '{0:39}'.format('process')
    for name in names:
        new_number_line += '{0:<25}'.format( number_dict[name] )
    lines[i_number_line] = new_number_line + '\n'

    # Write to file
    out_file = datacard_file.replace( '.txt', '_processesRenumbered.txt' )
    write_to_file(out_file, ''.join(lines))


def add_xH_nuisance_parameter(datacard_file, name_width=24, lnN_width=23, col_width=46):
    with open( datacard_file, 'r' ) as datacard_fp:
        lines = datacard_fp.readlines()

    for i_line, line in enumerate(lines):
        if line.startswith('process '):
            names = line.split()[1:]
            i_name_line = i_line
            break
    else:
        raise RuntimeError(
            'Reached end of file while searching two lines '
            'that start with \'process\'; datacard {0} is faulty.'
            .format(datacard_file)
            )

    nuis_name = 'CMS_xH_incxs'
    nuis_line = [ '{0:{nwidth}}{1:{lwidth}}'.format(nuis_name, 'lnN', nwidth=name_width, lwidth=lnN_width) ]

    nuis_size = 1. + LatestBinning.xH_unc_inclusive_fraction
    nuis_size_str = '{0:<{cwidth}.5f}'.format(nuis_size, cwidth=col_width)
    nuis_no_influence_str = '{0:{cwidth}}'.format('-', cwidth=col_width)

    for name in names:
        if name.startswith('xH'):
            nuis_line.append(nuis_size_str)
        else:
            nuis_line.append(nuis_no_influence_str)

    lines.append(''.join(nuis_line) + '\n')


    # Change kmax, the number of nuisance parameters
    for i_line, line in enumerate(lines):
        if line.startswith('kmax '):
            components = line.split()
            if components[1] == '*': break # * means no reason to track number of nuisances
            kmax_old = int(components[1])
            kmax_new = kmax_old + 1
            new_line = line.replace('kmax {0}'.format(kmax_old), 'kmax {0}'.format(kmax_new))
            lines[i_line] = new_line
            break
    else:
        raise RuntimeError(
            'Reached end of file while searching a line that starts with \'kmax\'; file {0} is faulty'
            .format(datacard_file)
            )

    
    out_file = datacard_file.replace( '.txt', '_xHNuisPar.txt' )
    write_to_file(out_file, ''.join(lines))




def process_sorter(process):
    if not isinstance(process, basestring):
        process = process[0]

    ret = 100000
    if process.startswith( 'ggH' ):
        ret -= 10000
    elif process.startswith( 'xH' ):
        ret -= 20000
    elif process.startswith( 'smH' ):
        ret -= 30000

    if process.startswith( 'ggH' ) or process.startswith( 'xH' ) or process.startswith( 'smH' ):
        components = process.split('_')
        if len(components) >= 3:
            leftBound = int( components[2].replace('GT','') )
            ret += leftBound
    else:
        ret += 1000
    return ret


def rename_processes_hgg_differentials(
        datacard_file,
        global_replace = None,
        tag = '_renamedProcesses'
        ):

    with open(datacard_file, 'r') as datacard_fp:
        datacard = datacard_fp.read()
    out = datacard

    # njets
    out = sub(
        r'(\W)InsideAcceptance_myGenNjets2p5_3p5to100p0(\W)',
        r'\1smH_NJ_GE4\2',
        out
        )
    out = sub(
        r'(\W)InsideAcceptance_myGenNjets2p5_([m\d]+)p5to(\d+)p5(\W)',
        r'\1smH_NJ_\3\4',
        out
        )
    # Slightly updated conventions for the NNLOPS datacards (from Nov 10)
    out = sub(
        r'(\W)InsideAcceptance_genNjets2p5_3p5_100p0(\W)',
        r'\1smH_NJ_GE4\2',
        out
        )
    out = sub(
        r'(\W)InsideAcceptance_genNjets2p5_([m\d]+)p5_(\d+)p5(\W)',
        r'\1smH_NJ_\3\4',
        out
        )

    # Rapidity renaming Nov12
    # Has to be done mosty manually... David used "0p30", Vittorio used "0p3"
    out = sub(
        r'(\W)InsideAcceptance_genAbsRapidity_0p0_0p15(\W)',
        r'\1smH_YH_0p0_0p15\2',
        out
        )
    out = sub(
        r'(\W)InsideAcceptance_genAbsRapidity_0p15_0p3(\W)',
        r'\1smH_YH_0p15_0p30\2',
        out
        )
    out = sub(
        r'(\W)InsideAcceptance_genAbsRapidity_0p3_0p6(\W)',
        r'\1smH_YH_0p30_0p60\2',
        out
        )
    out = sub(
        r'(\W)InsideAcceptance_genAbsRapidity_0p6_0p9(\W)',
        r'\1smH_YH_0p60_0p90\2',
        out
        )
    out = sub(
        r'(\W)InsideAcceptance_genAbsRapidity_0p9_3p0(\W)',
        r'\1smH_YH_0p90_2p50\2',
        out
        )
    # New bin boundaries
    out = sub(
        r'(\W)InsideAcceptance_genAbsRapidity_0p9_1p2(\W)',
        r'\1smH_YH_0p90_1p20\2',
        out
        )
    out = sub(
        r'(\W)InsideAcceptance_genAbsRapidity_1p2_3p0(\W)',
        r'\1smH_YH_1p20_2p50\2',
        out
        )

    # ptjet renaming Nov28
    out = sub(
        r'(\W)InsideAcceptance_genJet2p5Pt0_m1000p0_30p0(\W)',
        r'\1smH_PTJ_LT30\2',
        out
        )
    out = sub(
        r'(\W)InsideAcceptance_genJet2p5Pt0_200p0_13000p0(\W)',
        r'\1smH_PTJ_GT200\2',
        out
        )
    out = sub(
        r'(\W)InsideAcceptance_genJet2p5Pt0_([m\d]+)p0_(\d+)p0(\W)',
        r'\1smH_PTJ_\2_\3\4',
        out
        )

    # INC renaming May09
    out = sub(
        r'(\W)InsideAcceptance_genInclusive_0p0_2p0(\W)',
        r'\1smH_INC_INC\2',
        out
        )

    # Process simple global replacements
    if not global_replace == None:
        for string, replacement in global_replace:
            out = out.replace( string, replacement )

    # Write to file
    out_file = datacard_file.replace('.txt', '{0}.txt'.format(tag))
    write_to_file(out_file, out)
