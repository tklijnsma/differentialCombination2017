import LatestPaths
import differentials
import differentials.combine.combine as combine
import os, copy, traceback, shutil, logging
from os.path import *

#____________________________________________________________________
# Some scanning macros



def run_fastscan_scan_reused_postfit(config, postfit):
    config.make_unique_directory()
    postfit_and_fastscan_dir = join(config.subDirectory, 'postfit_and_fastscan')
    copied_postfit = join(postfit_and_fastscan_dir, basename(postfit))
    if differentials.core.is_testmode():
        logging.info('Would now create {0}'.format(postfit_and_fastscan_dir))
        logging.info('Would now copy {0} to {1}'.format(postfit, copied_postfit))
    else:
        logging.info('Creating {0}'.format(postfit_and_fastscan_dir))
        os.makedirs(postfit_and_fastscan_dir)
        logging.info('Copying {0} to {1}'.format(postfit, copied_postfit))
        shutil.copyfile(postfit, copied_postfit)

    fastscan = combine.CombineFastScan(config)
    fastscan.run(copied_postfit)
    fastscan_file = fastscan.get_output()

    pointwisescan = combine.CombinePointwiseScan(config)
    pointwisescan.run(copied_postfit, fastscan_file)



def run_postfit_fastscan_scan(config, point_minimizer_settings=None):
    # Make sure no previous run directory is overwritten
    config.make_unique_directory()

    postfit = combine.CombinePostfit(config)
    postfit.run()
    postfit_file = postfit.get_output()

    fastscan = combine.CombineFastScan(config)
    fastscan.run(postfit_file)
    fastscan_file = fastscan.get_output()

    if not(point_minimizer_settings is None):
        # Change the minimizer settings only for the scan, but not for the best fit
        config.minimizer_settings = point_minimizer_settings

    pointwisescan = combine.CombinePointwiseScan(config)
    pointwisescan.run(postfit_file, fastscan_file)


def run_postfit_scan(config, postfit_file=None):
    # Make sure no previous run directory is overwritten
    config.make_unique_directory()

    if postfit_file is None:
        postfit = combine.CombinePostfit(config)
        postfit.run()
        postfit_file = postfit.get_output()
    else:
        postfit_file = os.path.abspath(postfit_file)

    scan = combine.CombineScanFromPostFit(config)
    scan.run(postfit_file)

def run_postfit_fastscan(config):
    # Make sure no previous run directory is overwritten
    config.make_unique_directory()

    postfit = combine.CombinePostfit(config)
    postfit.run()
    postfit_file = postfit.get_output()

    fastscan = combine.CombineFastScan(config)
    fastscan.run(postfit_file)

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

def force_asimov(real_args, asimov=True):
    args = copy.deepcopy(real_args)
    args.asimov = asimov
    return args

def try_call_function_with_args(fn, args):
    try:
        fn(args)
    except Exception as exc:
        print traceback.format_exc()
        print exc
        pass
