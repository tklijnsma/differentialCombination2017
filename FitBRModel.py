
from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder

import os, sys, numpy, itertools, re
import ROOT
from math import sqrt
from copy import deepcopy
from array import array

class FitBRModelError(Exception):
    pass



class GeneralDifferentialModel( MultiSignalModel ):

    def __init__(self):
        MultiSignalModel.__init__(self)

        self.mHRange=[]
        self.mass = 125.
        self.verbose = 1


    def chapter( self, text, indent=0 ):
        if self.verbose:
            print '\n{tabs}{line}\n{tabs}{text}\n'.format(
                tabs = '    ' * indent,
                line = '-' * 70,
                text = text
                )
    def setPhysicsOptions( self, physOptions ):
        MultiSignalModel.setPhysicsOptions( self, physOptions )


    def BinIsHzz( self, bin ):
        # This is not robust coding AT ALL
        if bin.startswith('hzz'):
            if self.verbose: print '    Bin \'{0}\' is classified as a hzz bin!'.format( bin )
            return True
        elif bin.startswith('hgg'):
            if self.verbose: print '    Bin \'{0}\' is classified as a hgg bin!'.format( bin )
            return False
        else:
            raise FitBRModelError( 'Bin \'{0}\' raises an error when determining the channel'.format( bin ) )



class FitBRPerPTHBinModel(GeneralDifferentialModel):
    """docstring for FitBRPerPTHBinModel

    DROP THIS IDEA

    This can only work for the coarsest binning scheme.

    The idea is to let kappa_hzz_j modify the hzz BR, kappa_hgg_i modify the hgg BR
    Make kappa_hzz_j a function of kappa_hgg_i and R_j
    R_j cannot be R_i, that would be unconstrained.

    Would have been possible to do this for the hzz bins (and scaling 2 hgg bins by the same modifier: kappa_hgg_i = kappa_hgg_i+1 = kappa_hgg_j)
    But this is not worth the effort.
    """
    
    def __init__(self):
        super(FitBRPerPTHBinModel, self).__init__()

        self.hgg_bins = [ 0., 15., 30., 45., 85., 125., 200., 350. ]
        self.hzz_bins = [ 0., 15., 30., 85., 200. ]

        self.map_hgg_to_hzz = {
            0 : 0,
            1 : 1
            2 : 3
            3 : 3
            4 : 4
            5 : 4
            6 : 6
            7 : 6
            }

        self.n_bins_hgg = len(self.hgg_bins)-1
        self.n_bins_hzz = len(self.hzz_bins)-1

    def get_corresponding_hzz_bin(self, i_hgg):
        return self.map_hgg_to_hzz[i_hgg]

    def get_range_str_hgg(self, i):
        return '{0}_{1}'.format(int(self.hgg_bins[i]), int(self.hgg_bins[i+1]))

    def get_range_str_hzz(self, i):
        return '{0}_{1}'.format(int(self.hzz_bins[i]), int(self.hzz_bins[i+1]))


    def doParametersOfInterest(self):
        GeneralDifferentialModel.doParametersOfInterest( self )


        # Load spline for HZZ BR (function of MH)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        datadir = os.environ['CMSSW_BASE'] + '/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg'
        self.SMH.textToSpline( 'BR_hzz', os.path.join( datadir, 'sm/br/BR4.txt' ), ycol=11 );
        spline_SMBR_hzz = self.modelBuilder.out.function('BR_hzz')
        spline_SMBR_hgg = self.modelBuilder.out.function('fbr_13TeV')


        # self.modelBuilder.doVar( 'ratio_BR_hgg_hzz[0.086,0.0,0.5]' )
        # self.modelBuilder.doVar( 'hgg_BRmodifier[1.0,-2.0,4.0]' )


        for i_hgg in xrange(self.n_bins_hgg):
            self.modelBuilder.doVar('kappaBR_{0}[1.0,-2.0,4.0]'.format(self.get_range_str_hgg(i_hgg)))





        


class FitTotalXSModel( GeneralDifferentialModel ):

    def __init__(self):
        GeneralDifferentialModel.__init__(self)


    def doParametersOfInterest(self):
        GeneralDifferentialModel.doParametersOfInterest( self )
        self.modelBuilder.doVar( 'r[1.0,-2.0,4.0]' )


    def getYieldScale( self, bin, process ):

        if not self.DC.isSignal[process]:
            # If it's a background, always return the 1.0 yieldParameter
            yieldParameter = 1.0

        else:
            yieldParameter = 'r'

        print 'Scaling {0}/{1} by {2}'.format( bin, process, yieldParameter )

        return yieldParameter





class FitGlobalBRModel( GeneralDifferentialModel ):

    def __init__(self):
        GeneralDifferentialModel.__init__(self)


    def doParametersOfInterest(self):
        GeneralDifferentialModel.doParametersOfInterest( self )
        # self.chapter( 'Starting model.doParametersOfInterest()' )

        self.modelBuilder.doVar( 'hgg_BRmodifier[1.0,-2.0,4.0]' )
        # self.modelBuilder.doVar( 'hzz_BRmodifier[1.0,-2.0,4.0]' )


        # Load spline for HZZ BR (function of MH)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        datadir = os.environ['CMSSW_BASE'] + '/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg'
        self.SMH.textToSpline( 'BR_hzz', os.path.join( datadir, 'sm/br/BR4.txt' ), ycol=11 );
        spline_SMBR_hzz = self.modelBuilder.out.function('BR_hzz')

        # Seems to work:
        # self.modelBuilder.out.var('MH').setVal(125.09)
        # spline_SMBR_hzz.Print()

        # Load spling for Hgg BR
        spline_SMBR_hgg = self.modelBuilder.out.function('fbr_13TeV')

        #  --> Can't do this, has to be a spline
        # # HZZ BR from YR4 excel sheet at mH = 125.09
        # BR_hzz = 0.02641
        # self.modelBuilder.doVar( 'hzz_SMBR[{0}]'.format(BR_hzz) )
        # self.modelBuilder.out.var('hzz_SMBR').setConstant()

        # BR_hgg @ 125.09 GeV = 2.270E-03 = 0.00227
        # So BR_hgg/BR_hzz = 0.086
        # --> Varying ratio between 0.0 and 0.5 should suffice

        self.modelBuilder.doVar( 'ratio_BR_hgg_hzz[0.086,0.0,0.5]' )


        # ======================================
        # Create 'hzz_BRmodifier', as a function of 'hgg_BRmodifier' and 'ratio_BR_hgg_hzz'

        hzz_modifier_expr = 'expr::{name}("{formula}",{commaSeparatedParameters})'.format(
            name                     = 'hzz_BRmodifier',
            formula                  = '@0*(@1/@2)/@3',
            commaSeparatedParameters = ','.join([
                'hgg_BRmodifier',
                spline_SMBR_hgg.GetName(),
                spline_SMBR_hzz.GetName(),
                'ratio_BR_hgg_hzz'
                ])
            )
        print 'Processing expr:'
        print '    ',hzz_modifier_expr
        self.modelBuilder.factory_( hzz_modifier_expr )


        # for poiname in self.pois:
        #     if not poiname.startswith( 'r_' ): continue

        #     for channel in [ 'hzz', 'hgg' ]:
        #         newexpr = 'expr::{name}("{formula}",{commaSeparatedParameters})'.format(
        #             name                     = poiname.replace( 'r_', 'r_{0}_'.format(channel) ),
        #             formula                  = '@0*@1',
        #             commaSeparatedParameters = '{channel}_BRmodifier,{poiname}'.format( channel=channel, poiname=poiname )
        #             )
        #         print 'Processing expr:'
        #         print '    ',newexpr
        #         self.modelBuilder.factory_( newexpr )


        POIset = self.modelBuilder.out.set('POI')
        POIset.add( self.modelBuilder.out.var('hgg_BRmodifier') )
        POIset.add( self.modelBuilder.out.var('ratio_BR_hgg_hzz') )
        # POIset.add( self.modelBuilder.out.var('hzz_BRmodifier') )


        self.chapter( 'Starting model.getYieldScale()' )


    def getYieldScale( self, bin, process ):

        if process in self.DC.signals and not 'OutsideAcceptance' in bin:
            if self.BinIsHzz( bin ):
                poi = 'hzz_BRmodifier'
            else:
                poi = 'hgg_BRmodifier'
        else:
            poi = 1.0


        string = "%s/%s" % (bin,process)
        print "Will scale ", string, " by ", poi
        if poi in ["1","0"]: return int(poi)
        return poi;



class FitBRModel( GeneralDifferentialModel ):

    def __init__(self):
        GeneralDifferentialModel.__init__(self)


    def doParametersOfInterest(self):
        GeneralDifferentialModel.doParametersOfInterest( self )
        # self.chapter( 'Starting model.doParametersOfInterest()' )

        self.modelBuilder.doVar( 'hgg_BRmodifier[1.0,-2.0,4.0]' )
        # self.modelBuilder.doVar( 'hzz_BRmodifier[1.0,-2.0,4.0]' )


        # Load spline for HZZ BR (function of MH)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        datadir = os.environ['CMSSW_BASE'] + '/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg'
        self.SMH.textToSpline( 'BR_hzz', os.path.join( datadir, 'sm/br/BR4.txt' ), ycol=11 );
        spline_SMBR_hzz = self.modelBuilder.out.function('BR_hzz')

        # Seems to work:
        # self.modelBuilder.out.var('MH').setVal(125.09)
        # spline_SMBR_hzz.Print()

        # Load spling for Hgg BR
        spline_SMBR_hgg = self.modelBuilder.out.function('fbr_13TeV')

        #  --> Can't do this, has to be a spline
        # # HZZ BR from YR4 excel sheet at mH = 125.09
        # BR_hzz = 0.02641
        # self.modelBuilder.doVar( 'hzz_SMBR[{0}]'.format(BR_hzz) )
        # self.modelBuilder.out.var('hzz_SMBR').setConstant()

        # BR_hgg @ 125.09 GeV = 2.270E-03 = 0.00227
        # So BR_hgg/BR_hzz = 0.086
        # --> Varying ratio between 0.0 and 0.5 should suffice

        self.modelBuilder.doVar( 'ratio_BR_hgg_hzz[0.086,0.0,0.5]' )


        # ======================================
        # Create 'hzz_BRmodifier', as a function of 'hgg_BRmodifier' and 'ratio_BR_hgg_hzz'

        hzz_modifier_expr = 'expr::{name}("{formula}",{commaSeparatedParameters})'.format(
            name                     = 'hzz_BRmodifier',
            formula                  = '@0*(@1/@2)/@3',
            commaSeparatedParameters = ','.join([
                'hgg_BRmodifier',
                spline_SMBR_hgg.GetName(),
                spline_SMBR_hzz.GetName(),
                'ratio_BR_hgg_hzz'
                ])
            )
        print 'Processing expr:'
        print '    ',hzz_modifier_expr
        self.modelBuilder.factory_( hzz_modifier_expr )


        for poiname in self.pois:
            if not poiname.startswith( 'r_' ): continue

            for channel in [ 'hzz', 'hgg' ]:
                newexpr = 'expr::{name}("{formula}",{commaSeparatedParameters})'.format(
                    name                     = poiname.replace( 'r_', 'r_{0}_'.format(channel) ),
                    formula                  = '@0*@1',
                    commaSeparatedParameters = '{channel}_BRmodifier,{poiname}'.format( channel=channel, poiname=poiname )
                    )
                print 'Processing expr:'
                print '    ',newexpr
                self.modelBuilder.factory_( newexpr )


        POIset = self.modelBuilder.out.set('POI')
        POIset.add( self.modelBuilder.out.var('hgg_BRmodifier') )
        POIset.add( self.modelBuilder.out.var('ratio_BR_hgg_hzz') )
        # POIset.add( self.modelBuilder.out.var('hzz_BRmodifier') )


        self.chapter( 'Starting model.getYieldScale()' )


    def getYieldScale( self, bin, process ):

        string = "%s/%s" % (bin,process)
        poi = 1
        for p, list in self.poiMap:
            for l in list:
                if re.match(l, string): poi = p

        if isinstance( poi, basestring ) and poi.startswith( 'r_' ):
            if self.BinIsHzz( bin ):
                poi = poi.replace( 'r_', 'r_hzz_' )
            else:
                poi = poi.replace( 'r_', 'r_hgg_' )

        print "Will scale ", string, " by ", poi
        if poi in ["1","0"]: return int(poi)
        return poi;


fitBRModel=FitBRModel()
fitTotalXSModel=FitTotalXSModel()
fitGlobalBRModel=FitGlobalBRModel()