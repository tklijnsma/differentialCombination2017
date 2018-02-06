#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

# import os, itertools, operator, re, argparse, sys
# from math import isnan, isinf
# from os.path import *
# from glob import glob
# from copy import deepcopy

import combineCommands
import plotCommands
import yukawaCommands
import topCommands
import onetimeplotsCommands
# import highLumiStudyCommands
import extrastudyCommands
import crosscheckCommands
import scanCommands
# import rapidityCommands
# import njetsCommands
import differentialCombinations

import sys
sys.path.append('src')
import Commands
# import PhysicsCommands
# import OneOfCommands
import TheoryCommands
# import CorrelationMatrices
# import MergeHGGWDatacards
# import TheoryFileInterface

# from time import strftime
# datestr = strftime( '%b%d' )

# New style option handling
from OptionHandler import OptionHandler


########################################
# Main
########################################

def main():

    # ======================================
    # Parser

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument( '--bkg',                             action='store_true' )

    parser.add_argument( '--test',                            action='store_true' )

    parser.add_argument( '--fastscan',                        action='store_true' )
    parser.add_argument( '--notFastscan',                     action='store_true' )

    parser.add_argument( '--asimov',                          action='store_true' )
    parser.add_argument( '--notAsimov',                       action='store_true' )

    # combineCommands.AppendParserOptions(parser)
    # njetsCommands.AppendParserOptions(parser)
    plotCommands.AppendParserOptions(parser)
    yukawaCommands.AppendParserOptions(parser)
    topCommands.AppendParserOptions(parser)
    # highLumiStudyCommands.AppendParserOptions(parser)
    onetimeplotsCommands.AppendParserOptions(parser)
    extrastudyCommands.AppendParserOptions(parser)
    # rapidityCommands.AppendParserOptions(parser)
    scanCommands.AppendParserOptions(parser)
    crosscheckCommands.AppendParserOptions(parser)

    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument( '--latest', dest='latest', action='store_true', default=True )
    group.add_argument( '--older',  dest='latest', action='store_false' )

    group_decay = parser.add_mutually_exclusive_group(required=False)
    group_decay.add_argument( '--combination',   action='store_true' )
    group_decay.add_argument( '--hzz',           action='store_true' )
    group_decay.add_argument( '--hgg',           action='store_true' )
    group_decay.add_argument( '--hbb',           action='store_true' )
    group_decay.add_argument( '--combWithHbb',   action='store_true' )

    parser.add_argument( '--saveroot',   action='store_true' )
    parser.add_argument( '--savepng',   action='store_true' )

    optionHandler = OptionHandler()
    optionHandler.set_parser(parser)
    optionHandler.process_modules([
        'scanCommands',
        'differentialCombinations',
        'differentialPlots'
        ])

    parser.add_argument( '--statonly', action='store_true' )
    parser.add_argument( '--statsyst', action='store_true' )


    args = parser.parse_args()
    # print args
    # print ''

    if args.bkg:
        pass

    if args.test: Commands.TestMode()

    if args.saveroot:
        TheoryCommands.SaveAsRoot()
    if args.savepng:
        TheoryCommands.SaveAsPng()


    ########################################
    # New style options
    ########################################

    optionHandler.args = args
    optionHandler.execute_functions()


    ########################################
    # Encompassing module for submitting scans from yukawa and top
    ########################################

    # --> Moved to new style
    # if args.scanCommands:
    #     scanCommands.main(args)


    ########################################
    # Encompassing module for stuff related to the kappac kappab fits
    ########################################

    if args.yukawaCommands:
        yukawaCommands.main(args)


    ########################################
    # Encompassing module for stuff related to the ct, cg, cb fits
    ########################################

    if args.topCommands:
        topCommands.main(args)


    ########################################
    # Stuff dealing with combine (datacard merging/combining, t2ws, bestfits, scans, etc.)
    ########################################

    # if args.combineCommands:
    #     combineCommands.main(args)

    # if args.njetsCommands:
    #     njetsCommands.main(args)

    # if args.rapidityCommands:
    #     rapidityCommands.main(args)

    if args.extrastudyCommands:
        extrastudyCommands.main(args)


    ########################################
    # Result and Test Plotting
    ########################################

    if args.plotCommands:
        plotCommands.main(args)

    if args.onetimeplotsCommands:
        onetimeplotsCommands.main(args)

    if args.crosscheckCommands:
        crosscheckCommands.main(args)


    # ########################################
    # # High lumi study
    # ########################################

    # if args.highLumiStudyCommands:
    #     highLumiStudyCommands.main(args)


########################################
# End of Main
########################################
if __name__ == "__main__":
    main()