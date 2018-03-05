#!/usr/bin/env python
"""
Thomas Klijnsma
"""

import os, re
from os.path import *
from copy import deepcopy

import logging
import differentials.core as core
from OptionHandler import flag_as_option, flag_as_parser_options
import LatestPaths

#____________________________________________________________________
# Text changes to datacards

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
    raise NotImplementedError
    rename_processes_hgg_pth(
        ''
        )

@flag_as_option
def renumber_processes_hzz_pth_smH(args):
    renumber_processes_hzz_pth(
        'suppliedInput/fromDavid/PTH_Jan24_newBinning/smH/hzz4l_comb_13TeV_xs.txt'
        )

@flag_as_option
def disable_200_350_process_hbb_pth_ggH(args):
    disable_200_350_process_hbb_pth(
        'suppliedInput/fromJavier/bernstein_r7428/comb_2017_ggHbb.txt'
        )

#____________________________________________________________________
# combineCards

@flag_as_option
def combine_all_cards(args):
    for obsname in [
            # 'pth_smH',
            'pth_ggH',
            # 'njets',
            # 'rapidity',
            # 'ptjet'
            ]:
        combine_cards_for_observable(obsname)
    for obsname in [
            # 'pth_smH',
            'pth_ggH',
            ]:
        combine_cards_for_observable_with_hbb(obsname)

@flag_as_option
def combine_pth_ggH_hbb(args):
    combine_cards(
        'suppliedInput/combWithHbb_pth_ggH_{0}.txt'.format(core.datestr()),
        'hgg=' + LatestPaths.card.pth_ggH.hgg,
        'hzz=' + LatestPaths.card.pth_ggH.hzz,
        'hbb=' + LatestPaths.card.pth_ggH.hbb
        )

def combine_cards_for_observable(obsname):
    combine_cards(
        'suppliedInput/combination_{0}_{1}.txt'.format(obsname, core.datestr()),
        'hgg=' + LatestPaths.card[obsname].hgg,
        'hzz=' + LatestPaths.card[obsname].hzz
        )
def combine_cards_for_observable_with_hbb(obsname):
    combine_cards(
        'suppliedInput/combWithHbb_{0}_{1}.txt'.format(obsname, core.datestr()),
        'hgg=' + LatestPaths.card[obsname].hgg,
        'hzz=' + LatestPaths.card[obsname].hzz,
        'hbb=' + LatestPaths.card[obsname].hbb
        )

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

def sub(pat, repl, text, ntimes=3):
    """Does a replacement 3 times to process matches next to each other as well"""
    for i in xrange(ntimes):
        text = re.sub(pat, repl, text)
    return text

def write_to_file(out_file, contents):
    logging.debug('Contents of datacard {0}:\n{1}'.format(out_file, contents))
    logging.info('Opening {0} and dumping contents'.format(out_file))
    if core.is_testmode():
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

    # Process simple global replacements
    if not global_replace == None:
        for string, replacement in global_replace:
            out = out.replace( string, replacement )

    # Write to file
    out_file = datacard_file.replace('.txt', '{0}.txt'.format(tag))
    write_to_file(out_file, out)
