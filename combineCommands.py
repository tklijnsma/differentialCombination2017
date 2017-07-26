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

sys.path.append('src')
import Commands
import PhysicsCommands
import OneOfCommands
import TheoryCommands
import CorrelationMatrices
import MergeHGGWDatacards


from time import strftime
datestr = strftime( '%b%d' )


########################################
# Main
########################################

def AppendParserOptions( parser ):

    parser.add_argument( '--combineCommands', default=False )
    class CustomAction(argparse.Action):
        def __init__(self, option_strings, dest, **kwargs):
            # Only 1 argument allowed (this is basically the "store_true" action)
            super(CustomAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            # Make sure the parser knows a combineCommand was entered
            setattr( namespace, 'combineCommands', True )
            setattr( namespace, self.dest, True )
            # super(FooAction, self).__call__( parser, namespace, values, option_string=None )
            # setattr( namespace, self.dest, values )

    parser.add_argument( '--renamehgg',                       action=CustomAction )
    parser.add_argument( '--mergehgg',                        action=CustomAction )
    parser.add_argument( '--t2ws',                            action=CustomAction )
    parser.add_argument( '--bestfit',                         action=CustomAction )
    parser.add_argument( '--couplingT2WS',                    action=CustomAction )
    parser.add_argument( '--combineCards',                    action=CustomAction )
    parser.add_argument( '--couplingBestfit',                 action=CustomAction )
    parser.add_argument( '--couplingImportanceHistogram',     action=CustomAction )
    parser.add_argument( '--doFastscan',                      action=CustomAction )
    parser.add_argument( '--OutsideAcceptancetosignal',       action=CustomAction )

    # group = parser.add_mutually_exclusive_group(required=False)
    # group.add_argument( '--latest', dest='latest', action='store_true', default=True )
    # group.add_argument( '--older',  dest='latest', action='store_false' )



########################################
# Methods
########################################    

def main( args ):

    #____________________________________________________________________
    if args.OutsideAcceptancetosignal:
        OneOfCommands.ChangeOutsideAcceptanceToSignalProcess(
            'suppliedInput/fromDavid/PTH_May15/hzz4l_comb_13TeV_xs.txt',
            )

    #____________________________________________________________________
    if args.renamehgg:

        MergeHGGWDatacards.Rename_fea(
            'ggH',
            'suppliedInput/fromVittorio/pT_ggHonly_Jun26/CMS-HGG_sigfit_differential_pT_moriond17_ggHonly_v2.root'
            )

        MergeHGGWDatacards.Rename_fea(
            'xH',
            'suppliedInput/fromVittorio/differential_pT_moriond17_HxOnly_July19/CMS-HGG_sigfit_differential_pT_moriond17_HxOnly.root'
            )

    #____________________________________________________________________
    if args.mergehgg:

        xHCard  = 'suppliedInput/fromVittorio/differential_pT_moriond17_HxOnly_July19/Datacard_13TeV_differential_pT_moriond17_HxOnly.txt'
        ggHCard = 'suppliedInput/fromVittorio/pT_ggHonly_Jun26/Datacard_13TeV_differential_pT_moriond17_ggHonly_v2.txt'
        # xHCard  = 'suppliedInput/fromVittorio/differential_pT_moriond17_HxOnly_July19/CMS-HGG_sigfit_differential_pT_moriond17_HxOnly_feaRenamed.root'
        # ggHCard = 'suppliedInput/fromVittorio/pT_ggHonly_Jun26/CMS-HGG_sigfit_differential_pT_moriond17_ggHonly_v2_feaRenamed.root'

        MergeHGGWDatacards.Merge_xH_ggH_hgg_cards( xHCard, ggHCard )

        MergeHGGWDatacards.RenameProcesses( 'xH',  xHCard,  outdatacard='auto',
            globalReplace = [(
                'CMS-HGG_sigfit_differential_pT_moriond17_HxOnly.root',
                'CMS-HGG_sigfit_differential_pT_moriond17_HxOnly_feaRenamed.root'
                )]
            )
        MergeHGGWDatacards.RenameProcesses( 'ggH', ggHCard, outdatacard='auto',
            globalReplace = [(
                'CMS-HGG_sigfit_differential_pT_moriond17_ggHonly_v2.root',
                'CMS-HGG_sigfit_differential_pT_moriond17_ggHonly_v2_feaRenamed.root'
                )]
            )

        MergeHGGWDatacards.RenameProcesses(
            'smH',
            'suppliedInput/fromVittorio/pT_unsplit_Mar21/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5.txt',
            outdatacard='auto' )



    #____________________________________________________________________
    if args.combineCards:

        if args.latest:

            Commands.BasicCombineCards(
                'suppliedInput/combinedCard_{0}.txt'.format(datestr),
                'suppliedInput/fromVittorio/pT_unsplit_Mar21/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamedProcesses_Jul21.txt',
                'suppliedInput/fromDavid/PTH_May15/hzz4l_comb_13TeV_xs_processesShifted.txt'
                )

        else:
            print 'Nothing for args.combineCards'


    #____________________________________________________________________
    if args.t2ws:

        # ======================================
        # July 24

        if args.latest:
            print 'Running latest'

            Commands.BasicT2WS(
                'suppliedInput/fromVittorio/hggMerged_Jul24.txt',
                manualMaps=[
                    '--PO \'map=.*/ggH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_30_45:r_smH_PTH_30_45[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_45_85:r_smH_PTH_45_85[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_85_125:r_smH_PTH_85_125[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_125_200:r_smH_PTH_125_200[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_200_350:r_smH_PTH_200_350[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_GT350:r_smH_PTH_GT350[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_30_45:r_smH_PTH_30_45[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_45_85:r_smH_PTH_45_85[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_85_125:r_smH_PTH_85_125[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_125_200:r_smH_PTH_125_200[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_200_350:r_smH_PTH_200_350[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_GT350:r_smH_PTH_GT350[1.0,0.0,3.0]\'',
                    ],
                )

        else:
            print 'Running older'

            # unsplitCard = 'suppliedInput/combinedCard_Jul21.txt'
            unsplitCard = 'suppliedInput/combinedCard_Jul24.txt'

            Commands.BasicT2WS(
                unsplitCard,
                manualMaps=[
                    '--PO \'map=.*/smH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/smH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/smH_PTH_30_45:r_smH_PTH_30_45[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/smH_PTH_45_85:r_smH_PTH_45_85[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/smH_PTH_85_125:r_smH_PTH_85_125[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/smH_PTH_125_200:r_smH_PTH_125_200[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/smH_PTH_200_350:r_smH_PTH_200_350[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/smH_PTH_GT350:r_smH_PTH_GT350[1.0,0.0,3.0]\'',
                    ],
                )

            sys.exit()


            # ======================================
            # July 21

            Commands.BasicT2WS(
                'suppliedInput/fromVittorio/hggMerged_Jul20.txt',
                manualMaps=[
                    '--PO \'map=.*/ggH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_30_45:r_smH_PTH_30_45[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_45_85:r_smH_PTH_45_85[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_85_125:r_smH_PTH_85_125[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_125_200:r_smH_PTH_125_200[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_200_350:r_smH_PTH_200_350[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_GT350:r_smH_PTH_GT350[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_30_45:r_smH_PTH_30_45[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_45_85:r_smH_PTH_45_85[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_85_125:r_smH_PTH_85_125[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_125_200:r_smH_PTH_125_200[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_200_350:r_smH_PTH_200_350[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_GT350:r_smH_PTH_GT350[1.0,0.0,3.0]\'',
                    ],
                )


            renamed_xHCard  = 'suppliedInput/fromVittorio/differential_pT_moriond17_HxOnly_July19/Datacard_13TeV_differential_pT_moriond17_HxOnly_renamedProcesses_Jul21.txt'
            renamed_ggHCard = 'suppliedInput/fromVittorio/pT_ggHonly_Jun26/Datacard_13TeV_differential_pT_moriond17_ggHonly_v2_renamedProcesses_Jul21.txt'

            Commands.BasicT2WS(
                renamed_xHCard,
                manualMaps=[
                    # '--PO \'map=.*/ggH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
                    # '--PO \'map=.*/ggH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
                    # '--PO \'map=.*/ggH_PTH_30_45:r_smH_PTH_30_45[1.0,0.0,3.0]\'',
                    # '--PO \'map=.*/ggH_PTH_45_85:r_smH_PTH_45_85[1.0,0.0,3.0]\'',
                    # '--PO \'map=.*/ggH_PTH_85_125:r_smH_PTH_85_125[1.0,0.0,3.0]\'',
                    # '--PO \'map=.*/ggH_PTH_125_200:r_smH_PTH_125_200[1.0,0.0,3.0]\'',
                    # '--PO \'map=.*/ggH_PTH_200_350:r_smH_PTH_200_350[1.0,0.0,3.0]\'',
                    # '--PO \'map=.*/ggH_PTH_GT350:r_smH_PTH_GT350[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_0_15:r_xH_PTH_0_15[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_15_30:r_xH_PTH_15_30[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_30_45:r_xH_PTH_30_45[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_45_85:r_xH_PTH_45_85[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_85_125:r_xH_PTH_85_125[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_125_200:r_xH_PTH_125_200[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_200_350:r_xH_PTH_200_350[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/xH_PTH_GT350:r_xH_PTH_GT350[1.0,0.0,3.0]\'',
                    ],
                )

            Commands.BasicT2WS(
                renamed_ggHCard,
                manualMaps=[
                    '--PO \'map=.*/ggH_PTH_0_15:r_ggH_PTH_0_15[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_15_30:r_ggH_PTH_15_30[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_30_45:r_ggH_PTH_30_45[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_45_85:r_ggH_PTH_45_85[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_85_125:r_ggH_PTH_85_125[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_125_200:r_ggH_PTH_125_200[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_200_350:r_ggH_PTH_200_350[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/ggH_PTH_GT350:r_ggH_PTH_GT350[1.0,0.0,3.0]\'',
                    # '--PO \'map=.*/xH_PTH_0_15:r_xH_PTH_0_15[1.0,0.0,3.0]\'',
                    # '--PO \'map=.*/xH_PTH_15_30:r_xH_PTH_15_30[1.0,0.0,3.0]\'',
                    # '--PO \'map=.*/xH_PTH_30_45:r_xH_PTH_30_45[1.0,0.0,3.0]\'',
                    # '--PO \'map=.*/xH_PTH_45_85:r_xH_PTH_45_85[1.0,0.0,3.0]\'',
                    # '--PO \'map=.*/xH_PTH_85_125:r_xH_PTH_85_125[1.0,0.0,3.0]\'',
                    # '--PO \'map=.*/xH_PTH_125_200:r_xH_PTH_125_200[1.0,0.0,3.0]\'',
                    # '--PO \'map=.*/xH_PTH_200_350:r_xH_PTH_200_350[1.0,0.0,3.0]\'',
                    # '--PO \'map=.*/xH_PTH_GT350:r_xH_PTH_GT350[1.0,0.0,3.0]\'',
                    ],
                )

            unsplitCard = 'suppliedInput/fromVittorio/pT_unsplit_Mar21/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamedProcesses_Jul21.txt'
            Commands.BasicT2WS(
                unsplitCard,
                manualMaps=[
                    '--PO \'map=.*/smH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/smH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/smH_PTH_30_45:r_smH_PTH_30_45[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/smH_PTH_45_85:r_smH_PTH_45_85[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/smH_PTH_85_125:r_smH_PTH_85_125[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/smH_PTH_125_200:r_smH_PTH_125_200[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/smH_PTH_200_350:r_smH_PTH_200_350[1.0,0.0,3.0]\'',
                    '--PO \'map=.*/smH_PTH_GT350:r_smH_PTH_GT350[1.0,0.0,3.0]\'',
                    ],
                )


            # # Specific for HZZ
            # Commands.BasicT2WS(
            #     'suppliedInput/hzz_ggH_xH_split_Jun26/hzz4l_all_13TeV_xs.txt',
            #     manualMaps=[
            #         '--PO \'map=.*/ggH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/xH_PTH_0_15:r_smH_PTH_0_15[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/ggH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/xH_PTH_15_30:r_smH_PTH_15_30[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/ggH_PTH_30_45:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/xH_PTH_30_45:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/ggH_PTH_45_85:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/xH_PTH_45_85:r_smH_PTH_30_85[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/ggH_PTH_85_125:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/xH_PTH_85_125:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/ggH_PTH_125_200:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/xH_PTH_125_200:r_smH_PTH_85_200[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/ggH_PTH_200_350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/xH_PTH_200_350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/ggH_PTH_GT350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
            #         '--PO \'map=.*/xH_PTH_GT350:r_smH_PTH_GT200[1.0,0.0,3.0]\'',
            #         ],
            #     )


            # Commands.BasicT2WS( 'suppliedInput/pT_v5_renamed/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.txt' )
            # Commands.BasicT2WS( 'suppliedInput/PTH/hzz4l_comb_13TeV_xs.txt' )
            # Commands.BasicT2WS( 'suppliedInput/fromDavid/PTH_May09/hzz4l_comb_13TeV_xs.txt' )
            # Commands.BasicT2WS( 'suppliedInput/combinedCard_{0}.txt'.format(datestr) )
            # Commands.BasicT2WS( 'suppliedInput/pT_v5_renamed/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.txt' )



    #____________________________________________________________________
    if args.bestfit:

        if args.latest:

            Commands.BasicBestfit(
                'workspaces_Jul24/hggMerged_Jul24.root',
                # 'workspaces_Jul21/combinedCard_Jul21.root',
                # 
                # onBatch = True
                extraOptions = [
                    '--minimizerStrategy 2',
                    '-v 2',
                    '-m 125',
                    '--floatOtherPOIs=1',
                    ]
                )


        else:

            # merged hgg July 20
            Commands.BasicBestfit(
                'workspaces_Jul20/hggMerged_Jul20.root',
                # onBatch = True
                extraOptions = [
                    '--minimizerStrategy 2',
                    '-v 2',
                    '-m 125',
                    '--floatOtherPOIs=1',
                    ]
                )



            # hgg
            # Commands.BasicBestfit(
            #     'workspaces_May15/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.root',
            #     # onBatch = True
            #     )
            
            # hzz
            # Commands.BasicBestfit(
            #     'workspaces_May10/hzz4l_comb_13TeV_xs.root',
            #     # setPOIs = False,
            #     # onBatch = True
            #     )



            # May 15: workspaces with included but faulty _norm

            # extraVars = []
            # nProcesses = len(Commands.ListProcesses( 'suppliedInput/PTH/hzz4l_comb_13TeV_xs.txt' )[0])
            # extraVars.extend( [ 'K1Bin{0}'.format(i) for i in xrange(nProcesses) ] )
            # extraVars.extend( [ 'K2Bin{0}'.format(i) for i in xrange(nProcesses) ] )

            # # hzz
            # Commands.BasicBestfit(
            #     'workspaces_May15/hzz4l_comb_13TeV_xs.root',
            #     setPOIs = False,
            #     # onBatch = True,
            #     extraOptions = (
            #         '--saveSpecifiedFunc {0}'.format( ','.join(extraVars) )
            #         + ' -m 125.09 --floatOtherPOIs=1 --freezeNuisances MH'
            #         )
            #     )



            # Commands.BasicBestfit(
            #     'workspaces_May15/combinedCard_May15.root',
            #     # setPOIs = False,
            #     # onBatch = True
            #     extraOptions = (
            #         ' --floatOtherPOIs=1'
            #         # + ' -m 125.09 --freezeNuisances MH'
            #         )
            #     )



    #____________________________________________________________________
    if args.couplingT2WS:

        # Commands.TestMode()

        if args.latest:

            fullpath = lambda path: abspath('derivedTheoryFiles_Jul25/{0}'.format( path ) )

            Commands.BasicT2WSwithModel(
                # 'suppliedInput/combinedCard_May15.txt',
                # 'suppliedInput/combinedCard_Jul21.txt',
                'suppliedInput/combinedCard_Jul25.txt',
                'CouplingModel.py',
                extraOptions = [
                    '--PO verbose=2',
                    '--PO \'higgsMassRange=123,127\'',
                    '--PO linearTerms=True',
                    # Theory uncertainties
                    '--PO correlationMatrix=plots_CorrelationMatrices_Jul25/corrMat_exp.txt',
                    '--PO theoryUncertainties=plots_CorrelationMatrices_Jul25/errors_for_corrMat_exp.txt'
                    # SM
                    '--PO SM=[kappab=1,kappac=1,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_1_kappac_1.txt') ),
                    # For parametrization
                    '--PO theory=[kappab=-2,kappac=5,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_m2_kappac_5.txt') ),
                    '--PO theory=[kappab=1,kappac=0,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_1_kappac_0.txt') ),
                    '--PO theory=[kappab=2,kappac=-10,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_2_kappac_m10.txt') ),
                    '--PO theory=[kappab=-1,kappac=-5,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_m1_kappac_m5.txt') ),
                    '--PO theory=[kappab=2,kappac=10,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_2_kappac_10.txt') ),
                    # Rest
                    '--PO theory=[kappab=0,kappac=0,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_0_kappac_0.txt') ),
                    '--PO theory=[kappab=0,kappac=1,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_0_kappac_1.txt') ),
                    '--PO theory=[kappab=0,kappac=10,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_0_kappac_10.txt') ),
                    '--PO theory=[kappab=0,kappac=5,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_0_kappac_5.txt') ),
                    '--PO theory=[kappab=0,kappac=-10,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_0_kappac_m10.txt') ),
                    '--PO theory=[kappab=0,kappac=-5,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_0_kappac_m5.txt') ),
                    '--PO theory=[kappab=1,kappac=10,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_1_kappac_10.txt') ),
                    '--PO theory=[kappab=1,kappac=5,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_1_kappac_5.txt') ),
                    '--PO theory=[kappab=1,kappac=-10,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_1_kappac_m10.txt') ),
                    '--PO theory=[kappab=1,kappac=-5,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_1_kappac_m5.txt') ),
                    '--PO theory=[kappab=2,kappac=0,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_2_kappac_0.txt') ),
                    '--PO theory=[kappab=2,kappac=1,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_2_kappac_1.txt') ),
                    '--PO theory=[kappab=2,kappac=5,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_2_kappac_5.txt') ),
                    '--PO theory=[kappab=2,kappac=-5,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_2_kappac_m5.txt') ),
                    '--PO theory=[kappab=-1,kappac=0,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_m1_kappac_0.txt') ),
                    '--PO theory=[kappab=-1,kappac=1,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_m1_kappac_1.txt') ),
                    '--PO theory=[kappab=-1,kappac=10,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_m1_kappac_10.txt') ),
                    '--PO theory=[kappab=-1,kappac=5,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_m1_kappac_5.txt') ),
                    '--PO theory=[kappab=-1,kappac=-10,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_m1_kappac_m10.txt') ),
                    '--PO theory=[kappab=-2,kappac=0,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_m2_kappac_0.txt') ),
                    '--PO theory=[kappab=-2,kappac=1,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_m2_kappac_1.txt') ),
                    '--PO theory=[kappab=-2,kappac=10,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_m2_kappac_10.txt') ),
                    '--PO theory=[kappab=-2,kappac=-10,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_m2_kappac_m10.txt') ),
                    '--PO theory=[kappab=-2,kappac=-5,file={0}]'.format( fullpath('muR_1_muF_1_Q_1_kappab_m2_kappac_m5.txt') ),

                    # ======================================
                    # These do not contain the ratios

                    # '--PO SM=[kappab=+1,kappac=1,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_1_kappac_1.txt')),
                    # # Used for the parametrization
                    # '--PO theory=[kappab=-2,kappac=5,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_m2_kappac_5.txt')),
                    # '--PO theory=[kappab=+1,kappac=0,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_1_kappac_0.txt')),
                    # '--PO theory=[kappab=+2,kappac=-10,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_2_kappac_m10.txt')),
                    # # The rest
                    # '--PO theory=[kappab=+2,kappac=-5,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_2_kappac_m5.txt')),
                    # '--PO theory=[kappab=+0,kappac=0,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_0_kappac_0.txt')),
                    # '--PO theory=[kappab=+0,kappac=1,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_0_kappac_1.txt')),
                    # '--PO theory=[kappab=+0,kappac=10,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_0_kappac_10.txt')),
                    # '--PO theory=[kappab=+0,kappac=5,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_0_kappac_5.txt')),
                    # '--PO theory=[kappab=+0,kappac=-10,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_0_kappac_m10.txt')),
                    # '--PO theory=[kappab=+0,kappac=-5,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_0_kappac_m5.txt')),
                    # '--PO theory=[kappab=+1,kappac=10,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_1_kappac_10.txt')),
                    # '--PO theory=[kappab=+1,kappac=5,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_1_kappac_5.txt')),
                    # '--PO theory=[kappab=+1,kappac=-10,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_1_kappac_m10.txt')),
                    # '--PO theory=[kappab=+1,kappac=-5,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_1_kappac_m5.txt')),
                    # '--PO theory=[kappab=+2,kappac=0,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_2_kappac_0.txt')),
                    # '--PO theory=[kappab=+2,kappac=1,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_2_kappac_1.txt')),
                    # '--PO theory=[kappab=+2,kappac=10,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_2_kappac_10.txt')),
                    # '--PO theory=[kappab=+2,kappac=5,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_2_kappac_5.txt')),
                    # '--PO theory=[kappab=-1,kappac=0,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_m1_kappac_0.txt')),
                    # '--PO theory=[kappab=-1,kappac=1,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_m1_kappac_1.txt')),
                    # '--PO theory=[kappab=-1,kappac=10,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_m1_kappac_10.txt')),
                    # '--PO theory=[kappab=-1,kappac=5,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_m1_kappac_5.txt')),
                    # '--PO theory=[kappab=-1,kappac=-10,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_m1_kappac_m10.txt')),
                    # '--PO theory=[kappab=-1,kappac=-5,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_m1_kappac_m5.txt')),
                    # '--PO theory=[kappab=-2,kappac=0,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_m2_kappac_0.txt')),
                    # '--PO theory=[kappab=-2,kappac=1,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_m2_kappac_1.txt')),
                    # '--PO theory=[kappab=-2,kappac=10,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_m2_kappac_10.txt')),
                    # '--PO theory=[kappab=-2,kappac=-10,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_m2_kappac_m10.txt')),
                    # '--PO theory=[kappab=-2,kappac=-5,file={0}]'.format( abspath('derivedTheoryFiles_Jul07/muR_100_muF_100_Q_50_kappab_m2_kappac_m5.txt')),
                    ]
                )



        else:

            Commands.BasicT2WSwithModel(
                # 'suppliedInput/combinedCard_May15.txt',
                'suppliedInput/combinedCard_May15_theoryUncertainties.txt',
                'CouplingModel.py',
                extraOptions = [
                    '--PO verbose=1',
                    '--PO \'higgsMassRange=123,127\'',
                    '--PO SM=[ct=1.0,cg=0.0,file={0}]'.format( abspath('derivedTheoryFiles_May23/SM_NNLO.txt') ),
                    '--PO theory=[ct=0.1,cg=0.075,file={0}]'.format( abspath('derivedTheoryFiles_May23/ct_0p1_cg_0p075.txt') ),
                    '--PO theory=[ct=0.5,cg=0.042,file={0}]'.format( abspath('derivedTheoryFiles_May23/ct_0p5_cg_0p042.txt') ),
                    '--PO theory=[ct=1.5,cg=-0.042,file={0}]'.format( abspath('derivedTheoryFiles_May23/ct_1p5_cg_m0p042.txt') ),
                    '--PO theory=[ct=2.0,cg=-0.083,file={0}]'.format( abspath('derivedTheoryFiles_May23/ct_2p0_cg_m0p083.txt') ),
                    ]
                )



    #____________________________________________________________________
    if args.couplingImportanceHistogram:

        if args.doFastscan:

            datacard = 'workspaces_May30/combinedCard_May15.root'
            Commands.MultiDimCombineTool(
                datacard,
                nPoints       = 10000,
                nPointsPerJob = 10,
                queue         = 'all.q',
                # notOnBatch    = False,
                notOnBatch    = True,
                jobDirectory  = 'Fastscan_couplings_{0}'.format( datestr ),
                extraOptions  = [
                    '-P ct -P cg',
                    '--setPhysicsModelParameters ct=1.0,cg=0.0',
                    '--setPhysicsModelParameterRanges ct=-0.5,4.2:cg=-0.32,0.1',
                    '--fastScan',
                    '--saveSpecifiedFunc {0}'.format( ','.join(Commands.ListSet( datacard, 'yieldParameters' )) ),
                    ]
                )
            TheoryCommands.WriteTH2DToFile(
                glob( 'Fastscan_couplings_{0}/*.root'.format( datestr ) )
                )

        else:

            TheoryCommands.WriteTH2DToFile(
                # glob( 'Scan_couplings_May31_3/*.root' )
                glob( 'Fastscan_couplings_{0}/*.root'.format( datestr ) )
                )


    #____________________________________________________________________
    if args.couplingBestfit:

        TESTFIT = True

        # datacard = 'workspaces_May30/combinedCard_May15.root'
        # datacard = 'workspaces_Jun22/combinedCard_May15_theoryUncertainties.root'
        # datacard = 'workspaces_Jul24/combinedCard_Jul21.root'
        datacard = 'workspaces_Jul25/combinedCard_Jul25_CouplingModel.root'

        if TESTFIT:
            Commands.BasicBestfit(
                datacard,
                setPOIs = False,
                extraOptions = [
                    # '--setPhysicsModelParameters ct=1.0,cg=0.0'
                    '--setPhysicsModelParameters kappab=1.0,kappac=1.0'
                    ]
                )

        else:
            
            Commands.MultiDimCombineTool(
                datacard,
                nPoints       = 1600,
                # nPoints       = 100,
                nPointsPerJob = 10,
                queue         = 'all.q',
                notOnBatch    = False,
                # notOnBatch    = True,
                jobDirectory  = 'Scan_couplings_{0}'.format( datestr ),
                extraOptions  = [
                    '--importanceSampling={0}:couplingScan'.format( abspath('scanTH2D_Jun01.root') ),
                    '-P ct -P cg',
                    '--setPhysicsModelParameters ct=1.0,cg=0.0',
                    '--setPhysicsModelParameterRanges ct=-0.5,4.2:cg=-0.32,0.1',
                    # '--fastScan', # TURN THIS OFF FOR REAL RUN!!!!
                    '--saveSpecifiedFunc {0}'.format(','.join(
                        Commands.ListSet( datacard, 'yieldParameters' ) + [ i for i in Commands.ListSet( datacard, 'ModelConfig_NuisParams' ) if i.startswith('theoryUnc') ]  ) ),
                    '--squareDistPoiStep',
                    ]
                )

        # datacard = 'workspaces_May15/hzz4l_comb_13TeV_xs.root'
        # Commands.MultiDimCombineTool(
        #     datacard,
        #     # nPoints       = 1600,
        #     nPoints       = 100,
        #     nPointsPerJob = 10,
        #     queue         = 'all.q',
        #     # notOnBatch    = False,
        #     notOnBatch    = True,
        #     jobDirectory  = 'Scan_couplings_{0}_hzz1'.format( datestr ),
        #     extraOptions  = [
        #         '--importanceSampling={0}:couplingScan'.format( abspath('scanTH2D_Jun01.root') ),
        #         # '-P ct -P cg',
        #         '-P r_smH_PTH_0_15 -P r_smH_PTH_15_30',
        #         # '--setPhysicsModelParameters ct=1.0,cg=0.0',
        #         # '--setPhysicsModelParameterRanges ct=-0.5,4.2:cg=-0.32,0.1',
        #         '--fastScan', # TURN THIS OFF FOR REAL RUN!!!!
        #         # '--saveSpecifiedFunc {0}'.format( ','.join(Commands.ListSet( datacard, 'yieldParameters' )) ),
        #         # '--squareDistPoiStep',
        #         ]
        #     )








    # ======================================
    # combineTool.py

    # Commands.BasicCombineTool(
    #     'workspaces_Apr27/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.root',
    #     POIpattern    = '30_45',
    #     nPoints       = 10,
    #     nPointsPerJob = 2,
    #     )



    # # Ran on 15 May

    # Commands.BasicCombineTool(
    #     'workspaces_May15/usingDavidsCommand_OAshifted.root',
    #     POIpattern    = '*',
    #     nPoints       = 40,
    #     nPointsPerJob = 8,
    #     # notOnBatch    = True,
    #     jobDirectory  = 'Scan_{0}'.format(datestr),
    #     queue         = 'short.q',
    #     )

    # Commands.BasicCombineTool(
    #     'workspaces_May15/Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5_renamed.root',
    #     POIpattern    = '*',
    #     nPoints       = 40,
    #     nPointsPerJob = 4,
    #     # notOnBatch    = True,
    #     jobDirectory  = 'Scan_{0}'.format(datestr),
    #     queue         = 'short.q',
    #     )

    # Commands.BasicCombineTool(
    #     'workspaces_May15/combinedCard_May15.root',
    #     POIpattern    = '*',
    #     nPoints       = 40,
    #     nPointsPerJob = 4,
    #     # notOnBatch    = True,
    #     jobDirectory  = 'Scan_{0}'.format(datestr),
    #     queue         = 'short.q',
    #     )








########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'