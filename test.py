#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import argparse

# New style option handling
from OptionHandler import OptionHandler


########################################
# Main
########################################

def main():
    parser = argparse.ArgumentParser()

    group_decay_channel = parser.add_mutually_exclusive_group(required=False)
    group_decay_channel.add_argument( '--combination',   action='store_true' )
    group_decay_channel.add_argument( '--hzz',           action='store_true' )
    group_decay_channel.add_argument( '--hgg',           action='store_true' )
    group_decay_channel.add_argument( '--hbb',           action='store_true' )
    group_decay_channel.add_argument( '--combWithHbb',   action='store_true' )

    parser.add_argument( '--bkg',                             action='store_true' )
    parser.add_argument( '--test',                            action='store_true' )
    parser.add_argument( '--debug',                           action='store_true' )
    parser.add_argument( '--trace',                           action='store_true' )

    parser.add_argument( '--fastscan',                        action='store_true' )
    parser.add_argument( '--asimov',                          action='store_true' )
    parser.add_argument( '--notAsimov',                       action='store_true' )

    parser.add_argument( '--table',   action='store_true' )
    parser.add_argument( '--saveroot',   action='store_true' )
    parser.add_argument( '--savepng',   action='store_true' )
    parser.add_argument( '--savepng_convert',   action='store_true' )
    parser.add_argument( '--savegray',   action='store_true' )

    parser.add_argument( '--statonly', action='store_true' )
    parser.add_argument( '--statsyst', action='store_true' )
    parser.add_argument( '--lumiScale', action='store_true' )
    parser.add_argument( '--no-preliminary-tag', action='store_true' )

    #____________________________________________________________________
    # New style imports
    optionHandler = OptionHandler()
    optionHandler.set_parser(parser)
    optionHandler.process_modules([
        'datacard_preprocessing',
        # 
        'differentials_scans',
        'differentials_plots',
        # 
        'yukawa_scans',
        'yukawa_t2ws',
        'yukawa_plots',
        'top_scans',
        'top_plots',
        'top_t2ws',
        # 
        'debug',
        'inclusive_scans',
        # 'onetimeplotsCommands', # Not optimized for new code base
        # 
        'scalecorrelationmatrices',
        'crosschecks',
        'correlationmatrices',
        # 
        'fermilab',
        'projections_preprocessing',
        'projections_t2ws',
        'projections_scans',
        'projections_scans_kbkc',
        'projections_scans_ktcgkb',
        'projections_plots',
        'projections_plots_kbkc',
        'projections_plots_ktcgkb',
        # 
        'parametrization_plots',
        # 
        'thesis_plots'
        ])

    args = parser.parse_args()

    import differentials
    differentials.logger.set_basic_format()
    if args.test:
        differentials.core.testmode()
    if args.debug:
        differentials.logger.set_level_debug()
    if args.trace:
        differentials.logger.set_level_trace()

    import logging
    # args.no_preliminary_tag = True
    # args.no_preliminary_tag = False
    if args.no_preliminary_tag:
        logging.info('Remove the default \'Preliminary\' tag from plots')
        differentials.plotting.pywrappers.CMS_Latex_type.CMS_type_str = ''

    # projections_plotting = True
    projections_plotting = False

    if projections_plotting:
        differentials.plotting.pywrappers.CMS_Latex_type.CMS_type_str = 'Projection'
        differentials.plotting.plots.MultiContourPlot.z_axis_title = 'q'

    if args.saveroot:
        differentials.core.save_root()
    if args.savepng:
        differentials.core.save_png()
    if args.savepng_convert:
        differentials.core.save_png_through_convert()
    if args.savegray:
        differentials.core.save_gray()


    ########################################
    # New style options
    ########################################

    optionHandler.args = args
    optionHandler.execute_functions()


########################################
# End of Main
########################################
if __name__ == "__main__":
    main()