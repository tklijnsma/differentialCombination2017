from OptionHandler import flag_as_option, flag_as_parser_options

import LatestPaths
import sys
sys.path.append('src')
import Commands
import CombineToolWrapper
import differentialutils

from time import strftime
datestr = strftime( '%b%d' )

import os
from os.path import *


@flag_as_option
def totalXS_scan(args):
    if args.asimov:
        Commands.Warning('This is probably not what I want')
        return

    config = CombineToolWrapper.CombineConfig(args)
    config.onBatch       = True
    config.queue         = 'short.q'
    config.nPoints       = 55
    config.nPointsPerJob = config.nPoints

    r_ranges = [ 0.1, 2.0 ]
    config.POIs = [ 'r' ]
    config.PhysicsModelParameterRanges = [
        'r={0},{1}'.format( r_ranges[0], r_ranges[1] ),
        ]
    config.subDirectory = 'out/Scan_TotalXS_{0}'.format(datestr)

    config.datacard = LatestPaths.ws_combined_totalXS

    if args.asimov:
        config.subDirectory += '_asimov'
    config.make_unique_directory()

    postfit = CombineToolWrapper.CombinePostfit(config)
    postfit.run()
    postfit_file = postfit.get_output()

    # Stat+syst scan (regular)
    scan = CombineToolWrapper.CombineScanFromPostFit(config)
    scan.run(postfit_file)

    # Stat-only scan
    scan_stat_only = CombineToolWrapper.CombineScanFromPostFit(config)
    scan_stat_only.subDirectory += '_statonly'
    scan_stat_only.freezeNuisances.append('rgx{.*}')
    scan_stat_only.run(postfit_file)
