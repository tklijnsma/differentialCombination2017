#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import argparse

import combineCommands
import plotCommands
import yukawaCommands
import topCommands
import onetimeplotsCommands
import extrastudyCommands
import crosscheckCommands
import differentialCombinations

import sys
sys.path.append('src')
import Commands
import TheoryCommands

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

    parser.add_argument( '--fastscan',                        action='store_true' )
    parser.add_argument( '--asimov',                          action='store_true' )
    parser.add_argument( '--notAsimov',                       action='store_true' )

    parser.add_argument( '--saveroot',   action='store_true' )
    parser.add_argument( '--savepng',   action='store_true' )
    parser.add_argument( '--savepng_convert',   action='store_true' )

    parser.add_argument( '--statonly', action='store_true' )
    parser.add_argument( '--statsyst', action='store_true' )
    parser.add_argument( '--lumiScale', action='store_true' )

    #____________________________________________________________________
    # Old style imports
    yukawaCommands.AppendParserOptions(parser)
    topCommands.AppendParserOptions(parser)
    extrastudyCommands.AppendParserOptions(parser)
    crosscheckCommands.AppendParserOptions(parser)

    #____________________________________________________________________
    # New style imports
    optionHandler = OptionHandler()
    optionHandler.set_parser(parser)
    optionHandler.process_modules([
        'scans_yukawa',
        'scans_top',
        'scans_other',
        'differentialCombinations',
        'differentialPlots',
        'lumiStudyPlots',
        'onetimeplotsCommands',
        ])

    args = parser.parse_args()

    if args.bkg:
        pass
    if args.test: Commands.test_mode()

    if args.saveroot:
        TheoryCommands.save_as_root()
    if args.savepng:
        TheoryCommands.save_as_png()
    if args.savepng_convert:
        TheoryCommands.save_as_png_through_convert()


    ########################################
    # New style options
    ########################################

    optionHandler.args = args
    optionHandler.execute_functions()


    ########################################
    # Old style options
    ########################################

    if args.yukawaCommands:
        yukawaCommands.main(args)

    if args.topCommands:
        topCommands.main(args)

    if args.extrastudyCommands:
        extrastudyCommands.main(args)

    if args.crosscheckCommands:
        crosscheckCommands.main(args)


########################################
# End of Main
########################################
if __name__ == "__main__":
    main()