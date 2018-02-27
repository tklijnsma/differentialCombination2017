import os, re, random, copy
from os.path import *
from datetime import datetime

import logging
import logger

import ROOT

from time import strftime
GLOBAL_DATESTR = strftime( '%b%d' )

def datestr():
    return GLOBAL_DATESTR

TESTMODE = False
def testmode(flag=True):
    global TESTMODE
    TESTMODE = flag
    logger.enable_testmode()
    logging.info('Test mode enabled')

def is_testmode():
    global TESTMODE
    return TESTMODE

standard_titles = {
    'hgg' : 'H#rightarrow#gamma#gamma',
    'hzz' : 'H#rightarrowZZ',
    'combination' : 'Combination',
    'hbb' : 'H#rightarrowbb',
    'combWithHbb' : 'Comb. with H#rightarrowbb',
    # 
    'kappac' : '#kappa_{c}',
    'kappab' : '#kappa_{b}',
    'kappat' : '#kappa_{t}',
    'ct' : '#kappa_{t}',
    'cg' : 'c_{g}',
    # 
    'pth_smH'  : 'p_{T}^{H}',
    'pth_ggH'  : 'p_{T}^{H (ggH)}',
    'njets'    : 'N_{jets}',
    'ptjet'    : 'p_{T}^{jet}',
    'rapidity' : '|y_{H}|',
    }

def execute(cmd, capture_output=False, ignore_testmode=False):
    # Allow both lists and strings to be passed as the cmd
    if not isinstance(cmd, basestring):
        cmd = [ l for l in cmd if not len(l.strip()) == 0 ]
        cmd_str = '\n    '.join( cmd )
        cmd_exec = ' '.join(cmd)
    else:
        cmd_str = cmd
        cmd_exec = cmd

    logging.info('Executing the following command:\n{0}'.format(cmd_str))
    logging.trace('Actual passed command: {0}'.format(cmd_exec))
    if not(is_testmode()) and not(ignore_testmode):
        if capture_output:
            raise NotImplementedError('Re-implement using the capture output context manager')
        else:
            os.system(cmd_exec)


def float_to_str(number, nDecimals=None):
    number = float(number)
    if not nDecimals is None:
        string = '{:.{nDecimals}f}'.format(number, nDecimals=nDecimals).replace('-','m').replace('.','p')
        return string
    if number.is_integer():
        number = int(number)
    string = str(number).replace('-','m').replace('.','p')
    return string

def str_to_float(string):
    string = str(string)
    number = string.replace('m','-').replace('p','.')
    number = float(number)

    return number

def get_range_from_str(text):
    regular_match = re.search(r'([\dpm\.\-]+)_([\dpm\.\-]+)', text)
    overflow_match = re.search(r'(GE|GT)([\dpm\.\-]+)', text)
    underflow_match = re.search(r'(LE|LT)([\dpm\.\-]+)', text)
    single_match = re.search(r'_([\dpm\.\-]+)', text)

    if regular_match:
        left = str_to_float(regular_match.group(1))
        right = str_to_float(regular_match.group(2))
    elif overflow_match:
        left = str_to_float(overflow_match.group(2))
        right = 'INF'
    elif underflow_match:
        left = '-INF'
        right = str_to_float(underflow_match.group(2))
    elif single_match:
        logging.debug('single_match for {0}; matched text is {1}'.format(text, single_match.group(1)))
        left = str_to_float(single_match.group(1))
        right = 'SINGLE'
    else:
        left = 'UNDEFINED'
        right = 'UNDEFINED'

    return left, right

def range_sorter(text):
    left, right = get_range_from_str(text)
    if left == 'UNDEFINED':
        return 900000
    elif right == 'SINGLE':
        return left
    elif right == 'INF':
        return 800000
    elif left == '-INF':
        return -800000
    else:
        return left


def __uniqueid__():
    mynow=datetime.now
    sft=datetime.strftime
    # store old datetime each time in order to check if we generate during same microsecond (glucky wallet !)
    # or if daylight savings event occurs (when clocks are adjusted backward) [rarely detected at this level]
    old_time=mynow() # fake init - on very speed machine it could increase your seed to seed + 1... but we have our contingency :)
    # manage seed
    seed_range_bits=14 # max range for seed
    seed_max_value=2**seed_range_bits - 1 # seed could not exceed 2**nbbits - 1
    # get random seed
    seed=random.getrandbits(seed_range_bits)
    current_seed=str(seed)
    # producing new ids
    while True:
        # get current time 
        current_time=mynow()
        if current_time <= old_time:
            # previous id generated in the same microsecond or Daylight saving time event occurs (when clocks are adjusted backward)
            seed = max(1,(seed + 1) % seed_max_value)
            current_seed=str(seed)
        # generate new id (concatenate seed and timestamp as numbers)
        #newid=hex(int(''.join([sft(current_time,'%f%S%M%H%d%m%Y'),current_seed])))[2:-1]
        newid=int(''.join([sft(current_time,'%f%S%M%H%d%m%Y'),current_seed]))
        # save current time
        old_time=current_time
        # return a new id
        yield newid

class openroot():
    """Context manager to safely open and close root files"""
    def __init__(self, root_file):
        self._root_file = root_file

    def __enter__(self):
        if not isfile(self._root_file):
            raise IOError('File {0} does not exist'.format(self._root_file))
        self._root_fp = ROOT.TFile.Open(self._root_file)
        return self._root_fp

    def __exit__(self, *args):
        self._root_fp.Close()


def list_POIs(root_file, only_r_=True):
    with openroot(root_file) as root_fp:
        POI_list = ROOT.RooArgList(root_fp.Get('w').set('POI'))

    par_names = []
    for i in xrange(POI_list.getSize()):
        par_name = POI_list[i].GetName()
        if only_r_:
            if par_name.startswith('r_'):
                par_names.append(par_name)
            else:
                continue
        else:
            par_names.append(par_name)
    return par_names


def last_bin_is_overflow(POIs):
    """Checks if the last bin is an overflow bin. Assumes POIs are pre-sorted"""
    left, right = get_range_from_str(POIs[-1])
    is_overflow = False
    if right == 'INF':
        is_overflow = True
    logging.debug('Checking if POI \'{0}\' is overflow: found right {1}; overflow={2}'.format(POIs[-1], right, is_overflow))
    return is_overflow

def first_bin_is_underflow(POIs):
    """Checks if the first bin is an underflow bin. Assumes POIs are pre-sorted"""
    left, right = get_range_from_str(POIs[0])
    if left == '-INF':
        return True
    return False

def binning_from_POIs(POIs_original):
    POIs = copy.copy(POIs_original)
    POIs.sort(key=range_sorter)
    logging.debug(
        'Determining bin boundaries of the following (sorted) POIs:\n    '
        + '\n    '.join(POIs)
        )
    is_underflow = first_bin_is_underflow(POIs)
    is_overflow  = last_bin_is_overflow(POIs)
    binning = []
    if is_underflow:
        POIs = POIs[1:]
        binning.append(-10000)
    for POI in POIs:
        left, right = get_range_from_str(POI)
        binning.append(left)
    if is_overflow:
        binning.append(10000)
    else:
        binning.append(right)
    logging.debug('Determined the following binning: {0}'.format(binning))
    return binning


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

