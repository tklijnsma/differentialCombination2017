
#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

from OptionHandler import flag_as_option

import LatestPaths
import LatestPathsGetters
import LatestBinning

import differentials

import logging
import copy
import re
import random
random.seed(1002)

import sys
sys.path.append('src')
import TheoryFileInterface

from time import strftime
datestr = strftime('%b%d')

########################################
# Main
########################################

top_exp_binning = [ 0., 15., 30., 45., 80., 120., 200., 350., 600., 10000 ]

@flag_as_option
def top_scalecorrelations(args):

    # HIER VERDER
    # Variation files moeten nog gemaakt worden...

    variation_files = TheoryFileInterface.FileFinder(
        directory = LatestPaths.derivedTheoryFiles_YukawaSummed,
        kappab = 1, kappac = 1
        )
    variations = [
        TheoryFileInterface.ReadDerivedTheoryFile( variation_file, returnContainer=True )
            for variation_file in variation_files ]


    scaleCorrelation = process_variations(variations)
    # scaleCorrelation.make_scatter_plots(subdir='theory')
    scaleCorrelation.plot_correlation_matrix('theory')
    scaleCorrelation.write_correlation_matrix_to_file('theory')
    scaleCorrelation.write_errors_to_file('theory')

    variations_expbinning = deepcopy(variations)
    for variation in variations_expbinning:
        TheoryCommands.RebinDerivedTheoryContainer(variation, top_exp_binning)
    scaleCorrelation_exp = process_variations(variations_expbinning)
    scaleCorrelation_exp.make_scatter_plots(subdir='exp')
    scaleCorrelation_exp.plot_correlation_matrix('exp')
    scaleCorrelation_exp.write_correlation_matrix_to_file('exp')
    scaleCorrelation_exp.write_errors_to_file('exp')
