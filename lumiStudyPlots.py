#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, sys, re
from os.path import *
from glob import glob
from copy import deepcopy

from OptionHandler import flag_as_option

sys.path.append('src')
import Commands
import PhysicsCommands
import TheoryCommands
import LatestPaths
import LatestPathsGetters
import LatestBinning
from Container import Container
import PlotCommands
from differentialTools import *
import CombineToolWrapper

from time import strftime
datestr = strftime( '%b%d' )


########################################
# Main
########################################

lumiMultiplier300 = 8.356546
lumiMultiplier3000 = 83.56546


@flag_as_option
def oneKappa_scans(args):
    TheoryCommands.set_plot_dir( 'plots_{0}_lumiStudy'.format(datestr) )

    scans = [
        # ( 'kappac', LatestPaths.scan_combined_Yukawa_oneKappa_kappac ),
        ( 'kappac', LatestPaths.scan_combined_Yukawa_oneKappa_kappac_asimov ),
        # ( 'kappab', LatestPaths.scan_combined_Yukawa_oneKappa_kappab ),
        ( 'kappab', LatestPaths.scan_combined_Yukawa_oneKappa_kappab_asimov ),
        ]

    for kappa, scanDir in scans:

        container = Container()
        container.kappa = kappa
        container.x, container.y = TheoryCommands.basic_read_scan( glob( scanDir + '/*.root' ), kappa )

        container_300 = deepcopy(container)
        container_300.y = [ y * lumiMultiplier300 for y in container_300.y ]

        container_3000 = deepcopy(container)
        container_3000.y = [ y * lumiMultiplier3000 for y in container_3000.y ]

        if kappa == 'kappab':
            xMax = 7.
            xMin = -15.
        if kappa == 'kappac':
            xMax = 18.
            xMin = -45.

        PlotCommands.plot_multiple_scans(
            [ container, container_300, container_3000 ],
            xTitle   = kappa,
            yTitle   = '2#DeltaNLL',
            yMax     = 5.,
            xMin     = xMin, xMax     = xMax,
            plotname = 'oneKappaScan_{0}_{1}'.format( kappa, basename(scanDir).replace('/','') ),
            draw1sigmaline = True,
            draw2sigmaline = True,
            translateToChi2 = True,
            printUncertainties = True,
            )

