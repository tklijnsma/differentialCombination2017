#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, itertools, operator, re, argparse, sys
from math import isnan, isinf
from os.path import *
from glob import glob
from copy import deepcopy

import combineCommands
import plotCommands

sys.path.append('src')
import Commands
import PhysicsCommands
import OneOfCommands
import TheoryCommands
import CorrelationMatrices

from time import strftime
datestr = strftime( '%b%d' )


########################################
# Main
########################################

def main():

    # ======================================
    # Parser

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument( '--test',                            action='store_true' )

    parser.add_argument( '--newCommand',                      action='store_true' )

    combineCommands.AppendParserOptions(parser)
    plotCommands.AppendParserOptions(parser)

    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument( '--latest', dest='latest', action='store_true', default=True )
    group.add_argument( '--older',  dest='latest', action='store_false' )

    args = parser.parse_args()

    print args
    print ''



    if args.test:
        Commands.TestMode()



    ########################################
    # Stuff dealing with combine (datacard merging/combining, t2ws, bestfits, scans, etc.)
    ########################################

    # Moved to separate file
    if args.combineCommands:
        combineCommands.main(args)


    ########################################
    # Result and Test Plotting
    ########################################

    # Moved to separate file
    if args.plotCommands:
        plotCommands.main(args)



    ########################################
    # New commands
    ########################################


    if args.newCommand:

        print 'test'


























########################################
# End of Main
########################################
if __name__ == "__main__":
    main()
