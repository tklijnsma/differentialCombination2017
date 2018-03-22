import LatestPaths
import differentialTools

# import sys
# sys.path.append('src')
# import CombineToolWrapper
import differentials.combine.combine as combine

import os, copy, traceback

from os.path import *

#____________________________________________________________________
# Some scanning macros

def run_postfit_fastscan_scan(config):
    # Make sure no previous run directory is overwritten
    config.make_unique_directory()

    postfit = combine.CombinePostfit(config)
    postfit.run()
    postfit_file = postfit.get_output()

    fastscan = combine.CombineFastScan(config)
    fastscan.run(postfit_file)
    fastscan_file = fastscan.get_output()

    pointwisescan = combine.CombinePointwiseScan(config)
    pointwisescan.run(postfit_file, fastscan_file)


def run_postfit_scan(config):
    # Make sure no previous run directory is overwritten
    config.make_unique_directory()

    postfit = combine.CombinePostfit(config)
    postfit.run()
    postfit_file = postfit.get_output()

    scan = combine.CombineScanFromPostFit(config)
    scan.run(postfit_file)

def scan_directly(config):
    # Make sure no previous run directory is overwritten
    config.make_unique_directory()
    scan = combine.CombineScan(config)
    scan.run()

#____________________________________________________________________
# Argument interpretation

def assert_decay_channel(args, allowed_decay_channels=['combination']):
    decay_channel = differentials.core.get_decay_channel_tag(args)
    if not decay_channel in allowed_decay_channels:
        raise NotImplementedError(
            'Decay channel \'{0}\' is not allowed for this option '
            '(allowed: {1})'
            .format(decay_channel, allowed_decay_channels)
            )

def assert_asimov(args):
    if not args.asimov:
        raise RuntimeError(
            'Option is only allowed in combination with the flag --asimov'
            )

class DecayChannelNotFoundError(Exception):
    def __init__(self):
        Exception.__init__(self, 'Could not determine the decay channel') 

def get_decay_channel_tag(args, allow_default=False):
    if args.combination:
        tag = 'combination'
    elif args.hgg:
        tag = 'hgg'
    elif args.hzz:
        tag = 'hzz'
    elif args.hbb:
        tag = 'hbb'
    elif args.combWithHbb:
        tag = 'combWithHbb'
    else:
        if allow_default:
            tag = 'combination'
        else:
            raise DecayChannelNotFoundError()
    return tag

def set_one_decay_channel(real_args, decay_channel, asimov=None):
    args = copy.deepcopy(real_args)
    for dc in [ 'hgg', 'hzz', 'hbb', 'combination', 'combWithHbb' ]:
        setattr(args, dc, False)
    setattr(args, decay_channel, True)
    if not(asimov is None):
        args.asimov = asimov
    return args

def force_asimov(real_args):
    args = copy.deepcopy(real_args)
    args.asimov = True
    return args

def try_call_function_with_args(fn, args):
    try:
        fn(args)
    except Exception as exc:
        print traceback.format_exc()
        print exc
        pass
