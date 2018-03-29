from OptionHandler import flag_as_option, flag_as_parser_options

import LatestPaths
import sys

# sys.path.append('src')
# import Commands
# import CombineToolWrapper

import differentials
import differentialutils

from time import strftime
datestr = strftime( '%b%d' )

import os
from os.path import *


@flag_as_option
def totalXS_scan(args):
    if args.asimov:
        Commands.Warning('This is probably not what I want')
        return

    config = differentials.combine.CombineConfig(args)
    config.onBatch       = True
    config.queue         = 'short.q'
    config.nPoints       = 55
    config.nPointsPerJob = config.nPoints

    r_ranges = [ 0.1, 2.0 ]
    config.POIs = [ 'r' ]
    config.PhysicsModelParameterRanges = [
        'r={0},{1}'.format( r_ranges[0], r_ranges[1] ),
        ]
    config.subDirectory = 'out/Scan_TotalXS_{0}'.format(datestr)

    config.datacard = LatestPaths.ws_combined_totalXS

    if args.asimov:
        config.subDirectory += '_asimov'
    config.make_unique_directory()

    postfit = differentials.combine.CombinePostfit(config)
    postfit.run()
    postfit_file = postfit.get_output()

    # Stat+syst scan (regular)
    scan = differentials.combine.CombineScanFromPostFit(config)
    scan.run(postfit_file)

    # Stat-only scan
    scan_stat_only = differentials.combine.CombinePostfitScanFromPostFit(config)
    scan_stat_only.subDirectory += '_statonly'
    scan_stat_only.freezeNuisances.append('rgx{.*}')
    scan_stat_only.run(postfit_file)



@flag_as_option
def ratioBR_t2ws(args):
    card = LatestPaths.card.pth_smH.combination
    t2ws = differentials.combine.t2ws.T2WS(card, 'physicsModels/ExtendedMultiSignalModel.py')

    t2ws.extra_options.append('--PO \'get_splines\'')
    t2ws.add_expr('hgg_BRmodifier[1.0,-2.0,4.0]')
    t2ws.add_expr('ratio_BR_hgg_hzz[0.086,0.0,0.5]')

    hzz_modifier_expr = (
        'expr::hzz_BRmodifier("@0*(@1/@2)/@3",'
        'hgg_BRmodifier,'
        '#spline_hgg,' # Will be replaced by model
        '#spline_hzz,' # Will be replaced by model
        'ratio_BR_hgg_hzz'
        ')'
        )
    t2ws.add_expr(hzz_modifier_expr)

    t2ws.add_map('hzz.*/.*:sum::r_hzz(hzz_BRmodifier)')
    t2ws.add_map('hgg.*/.*:sum::r_hgg(hgg_BRmodifier)')

    t2ws.run()

# Gives:
# text2workspace.py
#     suppliedInput/combination_pth_smH_Mar14.txt
#     -o /mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/v4_NewBinning/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/differentialCombination2017/out/workspaces_Mar29/combination_pth_smH_Mar14_extendedMultiSignalModel.root
#     -P physicsModels.ExtendedMultiSignalModel:extendedMultiSignalModel
#     --PO verbose=2
#     --PO 'higgsMassRange=123,127'
#     --PO 'get_splines'
#     --PO 'factory=hgg_BRmodifier[1.0,-2.0,4.0]'
#     --PO 'factory=ratio_BR_hgg_hzz[0.086,0.0,0.5]'
#     --PO 'factory=expr::hzz_BRmodifier("@0*(@1/@2)/@3",hgg_BRmodifier,#spline_hgg,#spline_hzz,ratio_BR_hgg_hzz)'
#     --PO 'map=hzz.*/.*:sum::r_hzz(hzz_BRmodifier)'
#     --PO 'map=hgg.*/.*:sum::r_hgg(hgg_BRmodifier)'
# Will create a POI  sum::r_hzz(hzz_BRmodifier)  with factory  sum::r_hzz(hzz_BRmodifier)
# Mapping  sum::r_hzz(hzz_BRmodifier)  to  ['hzz.*/.*']  patterns
# Will create a POI  sum::r_hgg(hgg_BRmodifier)  with factory  sum::r_hgg(hgg_BRmodifier)
# Mapping  sum::r_hgg(hgg_BRmodifier)  to  ['hgg.*/.*']  patterns
# Making splines
# Processing factory hgg_BRmodifier[1.0,-2.0,4.0]
# Processing factory ratio_BR_hgg_hzz[0.086,0.0,0.5]
# Processing factory expr::hzz_BRmodifier("@0*(@1/@2)/@3",hgg_BRmodifier,fbr_13TeV,BR_hzz,ratio_BR_hgg_hzz)
# MH will be left floating within 123 and 127
# [#0] ERROR:InputArguments -- RooWorkspace::defineSet(w) ERROR proposed set constituent "sum::r_hgg(hgg_BRmodifier)" is not in workspace
# Traceback (most recent call last):
#   File "/mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/v4_NewBinning/CMSSW_7_4_7/bin/slc6_amd64_gcc491/text2workspace.py", line 56, in <module>
#     MB.doModel()
#   File "/mnt/t3nfs01/data01/shome/tklijnsm/differentialCombination2017/v4_NewBinning/CMSSW_7_4_7/python/HiggsAnalysis/CombinedLimit/ModelTools.py", line 97, in doModel
#     poiIter = self.out.set('POI').createIterator()
# ReferenceError: attempt to access a null-pointer



#     self.modelBuilder.doVar( 'hgg_BRmodifier[1.0,-2.0,4.0]' )
#     # self.modelBuilder.doVar( 'hzz_BRmodifier[1.0,-2.0,4.0]' )


#     # Load spline for HZZ BR (function of MH)
#     self.SMH = SMHiggsBuilder(self.modelBuilder)
#     datadir = os.environ['CMSSW_BASE'] + '/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg'
#     self.SMH.textToSpline( 'BR_hzz', os.path.join( datadir, 'sm/br/BR4.txt' ), ycol=11 );
#     spline_SMBR_hzz = self.modelBuilder.out.function('BR_hzz')

#     # Seems to work:
#     # self.modelBuilder.out.var('MH').setVal(125.09)
#     # spline_SMBR_hzz.Print()

#     # Load spling for Hgg BR
#     spline_SMBR_hgg = self.modelBuilder.out.function('fbr_13TeV')

#     #  --> Can't do this, has to be a spline
#     # # HZZ BR from YR4 excel sheet at mH = 125.09
#     # BR_hzz = 0.02641
#     # self.modelBuilder.doVar( 'hzz_SMBR[{0}]'.format(BR_hzz) )
#     # self.modelBuilder.out.var('hzz_SMBR').setConstant()

#     # BR_hgg @ 125.09 GeV = 2.270E-03 = 0.00227
#     # So BR_hgg/BR_hzz = 0.086
#     # --> Varying ratio between 0.0 and 0.5 should suffice

#     self.modelBuilder.doVar( 'ratio_BR_hgg_hzz[0.086,0.0,0.5]' )


#     # ======================================
#     # Create 'hzz_BRmodifier', as a function of 'hgg_BRmodifier' and 'ratio_BR_hgg_hzz'

#     hzz_modifier_expr = 'expr::{name}("{formula}",{commaSeparatedParameters})'.format(
#         name                     = 'hzz_BRmodifier',
#         formula                  = '@0*(@1/@2)/@3',
#         commaSeparatedParameters = ','.join([
#             'hgg_BRmodifier',
#             spline_SMBR_hgg.GetName(),
#             spline_SMBR_hzz.GetName(),
#             'ratio_BR_hgg_hzz'
#             ])
#         )
#     print 'Processing expr:'
#     print '    ',hzz_modifier_expr
#     self.modelBuilder.factory_( hzz_modifier_expr )


#     # for poiname in self.pois:
#     #     if not poiname.startswith( 'r_' ): continue

#     #     for channel in [ 'hzz', 'hgg' ]:
#     #         newexpr = 'expr::{name}("{formula}",{commaSeparatedParameters})'.format(
#     #             name                     = poiname.replace( 'r_', 'r_{0}_'.format(channel) ),
#     #             formula                  = '@0*@1',
#     #             commaSeparatedParameters = '{channel}_BRmodifier,{poiname}'.format( channel=channel, poiname=poiname )
#     #             )
#     #         print 'Processing expr:'
#     #         print '    ',newexpr
#     #         self.modelBuilder.factory_( newexpr )


#     POIset = self.modelBuilder.out.set('POI')
#     POIset.add( self.modelBuilder.out.var('hgg_BRmodifier') )
#     POIset.add( self.modelBuilder.out.var('ratio_BR_hgg_hzz') )
#     # POIset.add( self.modelBuilder.out.var('hzz_BRmodifier') )


#     self.chapter( 'Starting model.getYieldScale()' )


# def getYieldScale( self, bin, process ):

#     if process in self.DC.signals and not 'OutsideAcceptance' in bin:
#         if self.BinIsHzz( bin ):
#             poi = 'hzz_BRmodifier'
#         else:
#             poi = 'hgg_BRmodifier'
#     else:
#         poi = 1.0


#     string = "%s/%s" % (bin,process)
#     print "Will scale ", string, " by ", poi
#     if poi in ["1","0"]: return int(poi)
#     return poi;




@flag_as_option
def ratioBR_scan(args):
    pass







