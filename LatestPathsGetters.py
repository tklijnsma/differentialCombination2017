#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, sys
import LatestPaths

sys.path.append('src')
import differentialTools


########################################
# Main
########################################

vardict = vars(LatestPaths)
variables = vardict.keys()

def get_generic(prefix, obs_name, args=None, decay_channel=None, asimov=None, statonly=None):
    if decay_channel is None:
        if args is None:
            decay_channel = False
        else:
            decay_channel = differentialTools.get_decay_channel_tag(args)

    key = '{0}_{1}_{2}'.format(prefix, decay_channel, obs_name)

    if asimov is None:
        if args is None:
            asimov = False
        else:
            asimov = args.asimov
    if statonly is None:
        if args is None:
            statonly = False
        else:
            statonly = args.statonly

    if args.lumiScale and not(prefix == 'card'):
        key += '_lumiScale'
    if statonly:
        key += '_statonly'
    if asimov and not(prefix == 'ws' or prefix == 'card'):
        key += '_asimov'
    if not key in vardict:
        raise NotImplementedError(
            'Error getting {0} for decay_channel={1}, obs_name={2}'
            '\n    LatestPaths.{3} does not exist'
            .format(prefix, decay_channel, obs_name, key)
            )
    return vardict[key]

def get_ws(*args, **kwargs):
    return get_generic('ws', *args, **kwargs)
def get_postfit(*args, **kwargs):
    return get_generic('postfit', *args, **kwargs)
def get_scan(*args, **kwargs):
    return get_generic('scan', *args, **kwargs)
def get_card(*args, **kwargs):
    return get_generic('card', *args, **kwargs)
