import os, shutil
import logging
from os.path import *
import glob, re, copy
from collections import namedtuple
import sys

import differentials
import differentials.core as core

import combine_utils as utils



class T2WS(object):
    """docstring for T2WS"""

    default_model_file = 'HiggsAnalysis.CombinedLimit.PhysicsModel'

    def __init__(self, card=None, model_file=None, model_name=None, name=None):
        super(T2WS, self).__init__()
        self.name = name
        self.card = card
        self.ignore_xH = True

        self.model_file = self.default_model_file
        if not(model_file is None):
            self.model_file = model_file

        if model_name is None:
            if self.model_file == T2WS.default_model_file:
                self.model_name = 'multiSignalModel'
            else:
                base = basename(self.model_file).replace('.py','')
                self.model_name = base[0].lower() + base[1:]
        else:
            self.model_name = model_name

        self.outdir = 'workspaces_{0}'.format(core.datestr())
        self.output_ws = None
        self.tags = []
        self.extra_options = []

    def add_variable(self, name, val, x_min=None, x_max=None, is_POI=False):
        """Use only for variables (not expressions)"""
        if x_min is None and x_max is None:
            factory_str = '{0}[{1}]'.format(name, val)
        elif not(x_min is None) and not(x_max is None):
            factory_str = '{0}[{1},{2},{3}]'.format(name, val, x_min, x_max)
        else:
            raise ValueError('Pass both x_min and x_max, or neither')

        if is_POI:
            self.extra_options.append('--PO \'poi={0}\''.format(factory_str))
        else:
            self.extra_options.append('--PO \'variable={0}\''.format(factory_str))
        
    def add_expr(self, expr):
        """Use only for expressions (not variables)"""
        if self.model_file == T2WS.default_model_file:
            dummyname = 'dummy_{0}'.format(core.__uniqueid__().next())
            mapstr = '--PO \'map={0}:{1}\''.format(dummyname, expr)
            self.extra_options.insert(0, mapstr)
        else:
            self.extra_options.append('--PO \'expr={0}\''.format(expr))

    def add_map(self, mapstr):
        self.extra_options.append('--PO \'map={0}\''.format(mapstr))

    def add_maps(self, *args):
        for l in args:
            self.add_map(l)
        

    def get_processes_from_card(self):
        with open(self.card, 'r') as card_fp:
            lines = [ i.strip() for i in card_fp.readlines() ]
        process_lines = []
        for line in lines:
            if line.startswith('process'):
                process_lines.append(line)
                if len(process_lines) == 2: break
        else:
            raise RuntimeError('Could not find two lines in {0} that start with \'process\''.format(self.card))

        processes = []
        for process, integer in zip(process_lines[0].split()[1:], process_lines[1].split()[1:]):
            if int(integer) <= 0 and not 'OutsideAcceptance' in process:
                if self.ignore_xH and process.startswith('xH'): continue
                processes.append(process)
        processes = list(set(processes))
        logging.debug('Determined list of processes from {0}: {1}'.format(self.card, processes))
        return processes

    def make_maps_from_processes(self, binning=None, add_overflow=False, add_underflow=False, scale_ggH_xH_with_smH=False):
        processes = self.get_processes_from_card()
        self.processinterpreter = differentials.processinterpreter.ProcessInterpreter(processes, binning, scale_ggH_xH_with_smH)
        self.processinterpreter.make_yield_parameters(add_underflow=add_underflow, add_overflow=add_overflow)
        self.processinterpreter.link_processes_to_yield_parameters()
        self.extra_options.extend(self.processinterpreter.make_maps())

    def get_outdir(self):
        outdir = abspath(join(utils.get_global_outdir(), self.outdir))
        logging.debug('Output directory: {0}'.format(outdir))
        return outdir

    def get_output_ws(self):
        if self.name is None:
            if self.output_ws is None:
                output_ws = basename(self.card)
            else:
                output_ws = self.output_ws
            output_ws = output_ws.replace('.txt', '')
            output_ws += '_' + self.model_name
        else:
            output_ws = self.name

        if len(self.tags) > 0:
            output_ws += '_' + '_'.join(self.tags)
        output_ws += '.root'
        output_ws = join(self.get_outdir(), output_ws)
        logging.debug('Output ws: {0}'.format(output_ws))
        return output_ws

    def get_model_string(self):
        if self.model_file == self.default_model_file:
            model_file_string = self.model_file
        else:
            copy_physics_model_dir()
            model_file = relpath(self.model_file, os.getcwd())
            if model_file.startswith('/'): model_file = model_file[1:]
            if model_file.endswith('/'): model_file = model_file[:-1]
            model_file = model_file.replace('/','.').replace('.py','')
            model_file_string = model_file

        model_string = '{0}:{1}'.format(model_file_string, self.model_name)
        logging.debug('Model string: {0}'.format(model_string))
        return model_string

    def get_cmd(self):
        cmd = []
        cmd.append('text2workspace.py')
        cmd.append(self.card)
        cmd.append('-o {0}'.format(self.get_output_ws()))
        cmd.append('-P {0}'.format(self.get_model_string()))
        cmd.append('--PO verbose=2')
        cmd.append('--PO \'higgsMassRange=123,127\'')
        for line in self.extra_options:
            cmd.append(line)
        return cmd

    def run(self):
        logging.info('Creating {0} if not yet existing'.format(self.get_outdir()))
        if not core.is_testmode():
            if not isdir(self.get_outdir()):
                os.makedirs(self.get_outdir())

        cmd = self.get_cmd()
        core.execute(cmd)


def copy_physics_model_dir():
    """
    Copies models to compiled directory
    (scram b takes unnecessarily long)
    """
    physics_models_dir = 'physicsModels'
    dst = join(os.environ['CMSSW_BASE'], 'bin', os.environ['SCRAM_ARCH'], basename(physics_models_dir))
    logging.info('Copying {0} to {1}'.format(physics_models_dir, dst))
    if not core.is_testmode():
        if isdir(dst):
            shutil.rmtree(dst)
        shutil.copytree(physics_models_dir, dst)
