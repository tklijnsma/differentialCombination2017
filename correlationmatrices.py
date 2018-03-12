import logging
from OptionHandler import flag_as_option

import differentials
import differentials.core as core
import differentials.combine.combine as combine

import glob


def base_config(args):
    config = combine.CombineConfig(args)
    config.subDirectory = 'out/corrmats_{0}'.format(core.datestr())
    config.queue = 'short.q'
    config.onBatch = False
    return config

@flag_as_option
def corrmats_vittorio(args):
    config = base_config(args)
    for bestfit in glob.glob('out/bestfits_Vittorio_Mar12/*.root'):
        config.datacard = bestfit
        corrmat = combine.CombineCorrMat(config)
        corrmat.run()
