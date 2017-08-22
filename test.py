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
import highLumiStudyCommands

import sys
sys.path.append('src')
import Commands
# import PhysicsCommands
# import OneOfCommands
# import TheoryCommands
# import CorrelationMatrices
# import MergeHGGWDatacards
# import TheoryFileInterface

# from time import strftime
# datestr = strftime( '%b%d' )


########################################
# Main
########################################

def main():

    # ======================================
    # Parser

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument( '--test',                            action='store_true' )
    parser.add_argument( '--notFastscan',                     action='store_true' )

    combineCommands.AppendParserOptions(parser)
    plotCommands.AppendParserOptions(parser)
    yukawaCommands.AppendParserOptions(parser)
    topCommands.AppendParserOptions(parser)
    highLumiStudyCommands.AppendParserOptions(parser)

    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument( '--latest', dest='latest', action='store_true', default=True )
    group.add_argument( '--older',  dest='latest', action='store_false' )

    parser.add_argument( '--hzz',   action='store_true' )
    parser.add_argument( '--hgg',   action='store_true' )

    args = parser.parse_args()

    print args
    print ''

    if args.test: Commands.TestMode()


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

    if args.combineCommands:
        combineCommands.main(args)


    ########################################
    # Result and Test Plotting
    ########################################

    if args.plotCommands:
        plotCommands.main(args)


    ########################################
    # High lumi study
    ########################################

    if args.highLumiStudyCommands:
        highLumiStudyCommands.main(args)


########################################
# End of Main
########################################
if __name__ == "__main__":
    main()