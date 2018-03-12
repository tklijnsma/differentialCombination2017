import os, shutil
import logging
from os.path import *
import glob, re, copy
from collections import namedtuple

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

        self.model_file = self.default_model_file
        if not(model_file is None):
            self.model_file = model_file

        if model_name is None:
            self.model_name = 'multiSignalModel'
        else:
            self.model_name = model_name

        self.outdir = 'workspaces_{0}'.format(core.datestr())
        self.output_ws = None

        self.tags = []

        self.extra_options = []
        self.tags = []


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
