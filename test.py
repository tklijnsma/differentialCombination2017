#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os


import Commands

from time import strftime
datestr = strftime( '%b%d' )


########################################
# Main
########################################

def main():

    # Commands.TestMode()

    Commands.BasicT2WS( 'suppliedInput/pT_v5_renamed/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.txt' )
    Commands.BasicBestfit( 'workspaces_Apr27/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.root' )
    
    # Commands.BasicCombineTool(
    #     'workspaces_Apr27/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.root',
    #     POIpattern    = '30_45',
    #     nPoints       = 10,
    #     nPointsPerJob = 2,
    #     )



    # Commands.BasicCombineTool(
    #     'workspaces_Apr27/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.root',
    #     POIpattern    = '30_45',
    #     nPoints       = 32,
    #     nPointsPerJob = 4,
    #     # notOnBatch    = True,
    #     jobDirectory  = 'test002_{0}'.format(datestr),
    #     queue         = '8nh',
    #     )













########################################
# End of Main
########################################
if __name__ == "__main__":
    main()