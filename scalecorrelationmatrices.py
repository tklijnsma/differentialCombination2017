
#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

from OptionHandler import flag_as_option

# import LatestPaths
# import LatestPathsGetters
# import LatestBinning

import differentials

import logging
import copy
import re
import random
random.seed(1002)

import sys
sys.path.append('src')
# import TheoryFileInterface
import TheoryCommands

# from time import strftime
# datestr = strftime('%b%d')

########################################
# Main
########################################

top_exp_binning = [ 0., 15., 30., 45., 80., 120., 200., 350., 600., 10000 ]
yukawa_exp_binning = [0.0, 15., 30., 45., 80., 120.]


@flag_as_option
def top_scalecorrelations(args):

    variations = differentials.theory.theory_utils.FileFinder(
        ct=1.0, cg=0.0, cb=1.0,
        directory='out/theories_Mar05_tophighpt/'
        ).get()
    sm = [v for v in variations if v.muR==1.0 and v.muF==1.0 and v.Q==1.0][0]
    # variations.pop(variations.index(sm))

    scalecorrelation = differentials.theory.scalecorrelation.ScaleCorrelation()

    bin_boundaries = top_exp_binning[:]
    if 10000. in bin_boundaries:
        logging.warning('Treating 10000. as overflow in {0}'.format(bin_boundaries))
        bin_boundaries.pop(bin_boundaries.index(10000.))
    logging.warning('Limiting binning to the maximum of the theory calculations: {0}'.format(sm.binBoundaries[-1]))
    bin_boundaries.append(sm.binBoundaries[-1])

    scalecorrelation.set_bin_boundaries(bin_boundaries, add_overflow=False)

    for variation in variations:
        scalecorrelation.add_variation(
            TheoryCommands.Rebin(variation.binBoundaries, variation.crosssection, top_exp_binning),
            {},
            is_central=(variation is sm),
            )

    scalecorrelation.plot_correlation_matrix('tophighpt')
    scalecorrelation.make_scatter_plots(subdir='scatterplots_tophighpt')
    scalecorrelation.write_correlation_matrix_to_file('tophighpt')
    scalecorrelation.write_errors_to_file('tophighpt')

@flag_as_option
def yukawa_scalecorrelations(args):

    variations = differentials.theory.theory_utils.FileFinder(
        kappab=1.0, kappac=1.0,
        directory='out/theories_Mar09_yukawa_summed/'
        ).get()
    sm = [v for v in variations if v.muR==1.0 and v.muF==1.0 and v.Q==1.0][0]
    # variations.pop(variations.index(sm))

    scalecorrelation = differentials.theory.scalecorrelation.ScaleCorrelation()
    scalecorrelation.set_bin_boundaries(yukawa_exp_binning, add_overflow=False)

    for variation in variations:
        scalecorrelation.add_variation(
            TheoryCommands.Rebin(variation.binBoundaries, variation.crosssection, yukawa_exp_binning),
            {},
            is_central=(variation is sm),
            )

    scalecorrelation.plot_correlation_matrix('yukawa')
    scalecorrelation.make_scatter_plots(subdir='scatterplots_yukawa')
    scalecorrelation.write_correlation_matrix_to_file('yukawa')
    scalecorrelation.write_errors_to_file('yukawa')
