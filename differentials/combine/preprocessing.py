import os, shutil
import logging
from os.path import *
import glob, re, copy
from collections import namedtuple
import sys

import differentials
import differentials.core as core

import combine_utils as utils



class Nuisance(differentials.core.AttrDict):
    """docstring for Nuisance"""
    def __init__(self, name, type, **kwargs):
        super(Nuisance, self).__init__(
            name=name,
            type=type,
            **kwargs
            )


class NuisanceGroup(object):
    """docstring for NuisanceGroup"""
    def __init__(self, name, elements):
        super(NuisanceGroup, self).__init__()
        self.name = name
        self.elements = elements

class Datacard(object):
    """docstring for Datacard"""

    outdir = ''

    def __init__(self, card_file):
        super(Datacard, self).__init__()
        self.card_file = card_file
        self.get_lines()

        self.nuisance_groups = []
        self.nuisances = []


    def get_lines(self):
        with open(self.card_file, 'r') as card_fp:
            lines = card_fp.readlines()
        self.lines = []
        for i, line in enumerate(lines):
            line = line.strip()
            if len(line) == 0: continue
            if line.startswith('-----------------'): continue
            self.lines.append(line)

    def print_nuisances(self):
        if len(self.nuisances) == 0:
            self.get_nuisances()
        col_width = max([ len(nuis.name) for nuis in self.nuisances ])
        for nuis in self.nuisances:
            print '{0:{width}} - {1:10} {2}'.format(
                nuis.name, nuis.type,
                '- ' + nuis.group if not(nuis.group is None) and not(nuis.is_group) else '',
                width=col_width,
                )

    def get_nuisance_lines(self, return_indices=False):
        fill_in_nuisances = False
        nuis_lines = []
        for i_line, line in enumerate(self.lines):
            if fill_in_nuisances:
                if return_indices:
                    nuis_lines.append(i_line)
                else:
                    nuis_lines.append(line)
                continue
            if line.startswith('rate '):
                # Start filling the next lines
                fill_in_nuisances = True
        return nuis_lines
        
    def get_nuisances(self):
        nuis_lines = self.get_nuisance_lines()
        nuisances = [ self.interpret_nuis(line) for line in nuis_lines ]
        self.nuisance_groups, self.nuisances = self.assign_nuisances_to_group(nuisances)
        self.get_nuisance_dict()
        self.get_group_dict()

    def get_nuisance_dict(self):
        self.nuisance_dict = {}
        for name in self.get_nuisance_names():
            for nuis in self.nuisances:
                if nuis.name == name:
                    self.nuisance_dict[name] = nuis
                    break
            else:
                print 'Could not find nuis for name \'{0}\''.format(name)

    def get_group_dict(self):
        self.nuis_group_dict = { g['name'] : g for g in self.nuisance_groups }

    def get_nuisance_names(self):
        if len(self.nuisances) == 0: self.get_nuisances()
        return [ nuis.name for nuis in self.nuisances ]

    def interpret_nuis(self, line):
        components = line.split()
        name = components[0]
        nuistype = components[1]
        r = Nuisance(
            name = name,
            type = nuistype,
            is_group = False,
            elements = [],
            group = None,
            )
        if nuistype == 'group':
            r.is_group = True
            r.elements = components[2:]
        return r

    def assign_nuisances_to_group(self, nuisances):
        groups = [ n for n in nuisances if n.is_group ]
        nuisances = [ n for n in nuisances if not n.is_group ]
        for group in groups:
            for nuis in nuisances:
                if nuis.name in group.elements:
                    nuis.group = group.name
        return groups, nuisances

    def add_lumiscale_rateparam(self):
        self.lines.append('lumiscale rateParam * * 1')

    def add_group(self, group_name, containing_nuisances):
        self.nuisance_groups.append(
            NuisanceGroup(group_name, containing_nuisances)
            )

    def delete_current_nuisance_groups(self):
        to_delete = []
        nuis_lines = self.get_nuisance_lines(return_indices=True)
        for i_line in nuis_lines:
            line = self.lines[i_line]
            components = line.split()
            if len(components) >= 2 and components[1] == 'group':
                to_delete.append(i_line)

        if len(to_delete) == 0:
            logging.info('Found no groups to delete')
            return

        logging.info('Deleting the following lines because they are groups:')
        for i_line in to_delete:
            line = self.lines[i_line]
            if len(line) > 40: line = line[:40] + ' ...'
            logging.info(line)

        # High to low so pop() doesn't mess up future pops
        to_delete.sort(reverse=True)
        for i_line in to_delete:
            self.lines.pop(i_line)

    def parse(self):
        lines = self.lines[:]
        for group in self.nuisance_groups:
            line = '{0} group = {1}'.format(
                group.name,
                ' '.join(group.elements)
                )
            lines.append(line)
        return '\n'.join(lines)

    def out_to_temp(self):
        card_dir = os.path.dirname(self.card_file)
        out_temp = os.path.join(card_dir, '__temp.txt')
        dump_txt_to_file(self.parse(), out_temp)
        return out_temp

    def out_temp_to_combineCards(self, final_output):
        """
        Stage out in 2 steps:
        first to a temporary datacard in the same dir as the input card file
        second using combineCards to the final desired output path
        combineCards has the nice property that it takes care of the paths
        in the datacard to the root files
        """
        temp_out_card = self.out_to_temp()
        combine_cards(
            final_output,
            temp_out_card
            )

    def get_output_path(self, outname):
        if len(self.outdir) > 0:
            outname = os.path.join(self.outdir, outname)
        else:
            outname = os.path.join(os.path.dirname(self.card_file), outname)
        if not outname.endswith('.txt'): outname += '.txt'
        return outname

    def out(self, outname):
        full_output_path = self.get_output_path(outname)

        requested_outdir = os.path.dirname(full_output_path)
        input_card_dir   = os.path.dirname(self.card_file)
        if os.path.normpath(requested_outdir) != os.path.normpath(input_card_dir):
            logging.info(
                'Requested outdir {0} is not the input card dir {1}; will use a trick involving combineCards.py'
                .format(requested_outdir, input_card_dir)
                )
            self.out_temp_to_combineCards(full_output_path)
        else:
            dump_txt_to_file(self.parse(), full_output_path)


def dump_txt_to_file(text, out):
    actual_out_dir = os.path.dirname(out)
    if not os.path.isdir(actual_out_dir):
        if core.is_testmode():
            logging.info('Would now create directory {0}'.format(actual_out_dir))
        else:
            logging.info('Creating directory {0}'.format(actual_out_dir))
            os.makedirs(actual_out_dir)
    out = fix_extension_for_txt(out)
    if core.is_testmode():
        logging.info('Would now dump text to {0}'.format(out))
    else:
        logging.info('Dumping text to {0}'.format(out))
        with open(out, 'w') as out_fp:
            out_fp.write(text)

def combine_cards(
        output_file,
        *input_list
        ):
    cmd = [ 'combineCards.py' ]
    for datacard in input_list:
        if not isinstance(datacard, basestring) and len(datacard) == 2:
            cmd.append('{0}={1}'.format(datacard[0], datacard[1]))            
        else:
            cmd.append(datacard)
    output_file = fix_extension_for_txt(output_file)
    cmd.append( '> {0}'.format( output_file ) )
    core.execute(cmd)

def fix_extension_for_txt(out):
    if not out.endswith('.txt'): out += '.txt'
    out = out.replace('.txt', '_{0}.txt'.format(core.datestr()))
    return out




