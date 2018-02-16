#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, sys, re
from os.path import *
from copy import deepcopy

import Commands
import TheoryCommands
from Container import Container

from time import strftime
datestr = strftime( '%b%d' )

import ROOT

from TheoryCommands import c
from TheoryCommands import SaveC
from TheoryCommands import GetPlotBase
from TheoryCommands import SetCMargins


########################################
# Main
########################################


class MuShapeDrawer(object):

    def __init__( self, w ):
        
        if isinstance( w, basestring ):
            wFile = w
            wFp = ROOT.TFile.Open( wFile )
            w_temp = wFp.Get('w')
            # w = ROOT.RooWorkspace(w_temp)
            # ROOT.SetOwnership( w, False )
            self.w = deepcopy(w_temp)
            wFp.Close()

        elif isinstance( w, ROOT.RooWorkspace ):
            self.w = deepcopy(w)
        else:
            Commands.throw_error( 'Pass either a path to a root file or a RooWorkspace', throwException=True )

        self.w.loadSnapshot('MultiDimFit')

        self._category = self.w.cat('CMS_channel')
        ROOT.SetOwnership( self._category, False )

        self.mH_hgg   = self.w.var('CMS_hgg_mass')
        ROOT.SetOwnership( self.mH_hgg, False )
        self.default_mH_hgg_range = [ self.mH_hgg.getMin(), self.mH_hgg.getMax() ]
        self.default_mH_hgg_nBins = self.mH_hgg.getBins()

        self.mH_hzz   = self.w.var('CMS_zz4l_mass')
        ROOT.SetOwnership( self.mH_hzz, False )
        self.default_mH_hzz_range = [ self.mH_hzz.getMin(), self.mH_hzz.getMax() ]
        self.default_mH_hzz_nBins = self.mH_hzz.getBins()


        self.allPdfsInWorkspace = []
        allPdfsArgList = ROOT.RooArgList( self.w.allPdfs() )
        for iPdf in xrange( allPdfsArgList.getSize() ):
            self.allPdfsInWorkspace.append( allPdfsArgList[iPdf].GetName() )
        self.allPdfsInWorkspace.sort()

        self.yieldParameterNames = [ y for y in Commands.list_set( self.w, 'POI' ) if y.startswith('r_') ]
        getLeftBound = lambda y: Commands.convert_str_to_float( re.search( r'_[GLET]*([\dpm]+)(_|$)', y ).group(1) )
        self.yieldParameterNames.sort( key = getLeftBound )
        self.yieldParameterDefaultValues = [ self.w.var(yP).getVal() for yP in self.yieldParameterNames ]


        self.binWidth = 0.25 # GeV

        self.hgg_binBoundaries = [ 0., 15., 30., 45., 85., 125., 200., 350 ]
        self.hzz_binBoundaries = [ 0., 15., 30., 85., 200. ]

        self.do_hgg = False
        self.do_hzz = False

        self.drawBestfit = False


    # Easy setter and getter for the CMS_channel index
    def get_category( self ):
        return self._category.getIndex()
    def set_category( self, val ):
        if isinstance( val, int ):
            self._category.setIndex(val)
        elif isinstance( val, basestring ):
            self._category.setLabel(val)
            # print '\nSet category to new value {0}:'.format(val)
            # self._category.Print()
        else:
            Commands.throw_error( 'MuShapeDrawer.category needs an integer or a string', throwException=True )
    category = property( get_category, set_category )


    # def get_mH( self ):
    #     return self._mH_hgg.getVal()
    # def set_mH( self, val ):
    #     self._mH_hgg.setVal(val)
    #     self._mH_hzz.setVal(val)
    # mH = property( get_mH, set_mH )

    # def set_mH_range( self, xMin, xMax ):
    #     self._mH_hgg.setRange( xMin, xMax )
    #     self._mH_hzz.setRange( xMin, xMax )



    def get_bin_str(
            self,
            hggBin = None,
            hzzBin = None,
            hggCat = None,
            hzzCat = None,
            ):

        PdfFound = lambda pdf: 'pdf_bin{0}'.format(pdf) in self.allPdfsInWorkspace

        if hggBin != None and hggCat != None:

            if not hggCat in [ 0, 1, 2 ]:
                Commands.throw_error( 'hggCat should be 0, 1 or 2', throwException=True )

            hggBinStr = {
                0 : 'ch1_recoPt_0p0_15p0_SigmaMpTTag_{0}_13TeV'.format(hggCat),
                1 : 'ch1_recoPt_15p0_30p0_SigmaMpTTag_{0}_13TeV'.format(hggCat),
                2 : 'ch1_recoPt_30p0_45p0_SigmaMpTTag_{0}_13TeV'.format(hggCat),
                3 : 'ch1_recoPt_45p0_85p0_SigmaMpTTag_{0}_13TeV'.format(hggCat),
                4 : 'ch1_recoPt_85p0_125p0_SigmaMpTTag_{0}_13TeV'.format(hggCat),
                5 : 'ch1_recoPt_125p0_200p0_SigmaMpTTag_{0}_13TeV'.format(hggCat),
                6 : 'ch1_recoPt_200p0_350p0_SigmaMpTTag_{0}_13TeV'.format(hggCat),
                7 : 'ch1_recoPt_GT350p0_SigmaMpTTag_{0}_13TeV'.format(hggCat),
                }[hggBin]
            if PdfFound(hggBinStr):
                return hggBinStr
            else:
                print 'Attempt \'{0}\' is not a valid pdf'.format(hggBinStr)

            hggBinStr = {
                0 : 'ch1_SigmaMpTTag_{0}_recoPt_0p0_15p0_13TeV'.format(hggCat),
                1 : 'ch1_SigmaMpTTag_{0}_recoPt_15p0_30p0_13TeV'.format(hggCat),
                2 : 'ch1_SigmaMpTTag_{0}_recoPt_30p0_45p0_13TeV'.format(hggCat),
                3 : 'ch1_SigmaMpTTag_{0}_recoPt_45p0_85p0_13TeV'.format(hggCat),
                4 : 'ch1_SigmaMpTTag_{0}_recoPt_85p0_125p0_13TeV'.format(hggCat),
                5 : 'ch1_SigmaMpTTag_{0}_recoPt_125p0_200p0_13TeV'.format(hggCat),
                6 : 'ch1_SigmaMpTTag_{0}_recoPt_200p0_350p0_13TeV'.format(hggCat),
                7 : 'ch1_SigmaMpTTag_{0}_recoPt_GT350p0_13TeV'.format(hggCat),
                }[hggBin]
            if PdfFound(hggBinStr):
                return hggBinStr
            else:
                print 'Attempt \'{0}\' is not a valid pdf'.format(hggBinStr)

            hggBinStr = {
                0 : 'ch1_SigmaMpTTag_{0}_recoPt_0p0to15p0_13TeV'.format(hggCat),
                1 : 'ch1_SigmaMpTTag_{0}_recoPt_15p0to30p0_13TeV'.format(hggCat),
                2 : 'ch1_SigmaMpTTag_{0}_recoPt_30p0to45p0_13TeV'.format(hggCat),
                3 : 'ch1_SigmaMpTTag_{0}_recoPt_45p0to85p0_13TeV'.format(hggCat),
                4 : 'ch1_SigmaMpTTag_{0}_recoPt_85p0to125p0_13TeV'.format(hggCat),
                5 : 'ch1_SigmaMpTTag_{0}_recoPt_125p0to200p0_13TeV'.format(hggCat),
                6 : 'ch1_SigmaMpTTag_{0}_recoPt_200p0to350p0_13TeV'.format(hggCat),
                7 : 'ch1_SigmaMpTTag_{0}_recoPt_350p0to10000p0_13TeV'.format(hggCat),
                }[hggBin]
            if PdfFound(hggBinStr):
                return hggBinStr
            else:
                print 'Attempt \'{0}\' is not a valid pdf'.format(hggBinStr)

            Commands.throw_error( 'No pdf available for hggBin={0}, hggCat={1}\nList of all PDFs in workspace:\n  {2}'.format( hggBin, hggCat, '\n  '.join(self.allPdfsInWorkspace) ), throwException=True )


        elif hzzBin != None and hzzCat != None:

            if not hzzCat in [ '4mu', '4e', '2e2mu' ]:
                Commands.throw_error( 'hzzCat should be 4mu, 4e or 2e2mu', throwException=True )

            hzzBinStr = {
                0 : 'ch2_hzz_PTH_0_15_cat{0}'.format(hzzCat),
                1 : 'ch2_hzz_PTH_15_30_cat{0}'.format(hzzCat),
                2 : 'ch2_hzz_PTH_30_85_cat{0}'.format(hzzCat),
                3 : 'ch2_hzz_PTH_85_200_cat{0}'.format(hzzCat),
                4 : 'ch2_hzz_PTH_GE200_cat{0}'.format(hzzCat),
                }[hzzBin]

            if PdfFound(hzzBinStr):
                return hzzBinStr
            else:
                print 'Attempt \'{0}\' is not a valid pdf'.format(hzzBinStr)

            Commands.throw_error( 'No pdf available for hzzBin={0}, hzzCat={1}\nList of all PDFs in workspace:\n  {2}'.format( hzzBin, hzzCat, '\n  '.join(self.allPdfsInWorkspace) ), throwException=True )
            
        else:
            Commands.throw_error( 'Not implemented', throwException=True )


    def get_yield_parameter(
            self,
            iBin,
            ):

        if self.do_hzz:
            iBin = {
                0 : 0,
                1 : 1,
                2 : 2,
                # 2 : 3,
                3 : 4,
                # 3 : 5,
                4 : 6,
                # 4 : 7,
                }[iBin]

        return self.w.var( self.yieldParameterNames[iBin] )


    def draw_shapes( self, **kwargs ):

        if 'drawBestfit' in kwargs:
            self.drawBestfit = kwargs['drawBestfit']
            del kwargs['drawBestfit']

        res = self.get_shapes( **kwargs )
        self.draw_shape_result(res)


    def get_shapes(
            self,
            hggBin   = None,
            hggCat   = None,
            hzzBin   = None,
            hzzCat   = None,
            muValues = [ 0.5, 1.0, 2.0 ],
            ):

        if hggBin != None and hggCat != None and hzzBin==None and hzzCat==None:
            self.do_hgg = True
            self.do_hzz = False
        elif hzzBin != None and hzzCat != None and hggBin==None and hggCat==None:
            self.do_hgg = False
            self.do_hzz = True
        else:
            Commands.throw_error( 'Can\'t decide between hgg and hzz:\n    hggBin={0}, hggCat={1}, hzzBin={2}, hzzCat={3}'.format( hggBin, hggCat, hzzBin, hzzCat ), throwException=True )

        res = Container()

        if self.do_hgg:
            mH = self.mH_hgg
            default_mH_range = self.default_mH_hgg_range
            binStr = self.get_bin_str( hggBin=hggBin, hggCat=hggCat )
            iBin = hggBin
            iCat = hggCat
        elif self.do_hzz:
            mH = self.mH_hzz
            default_mH_range = self.default_mH_hzz_range
            binStr = self.get_bin_str( hzzBin=hzzBin, hzzCat=hzzCat )
            iBin = hzzBin
            iCat = hzzCat

        # Set the proper category
        self.category = binStr

        # Get bkg and sig pdfs
        mH.setRange( *default_mH_range )
        mH.setBins( int((default_mH_range[1]-default_mH_range[0])/self.binWidth) )

        bkgPdf = self.w.pdf( 'pdf_bin{0}_bonly'.format(binStr) )
        Hbkg = bkgPdf.createHistogram( 'bkg_'+binStr, mH )
        Hbkg.Scale( self.binWidth ) # default is /GeV

        Hsigs = []
        yieldParameter = self.get_yield_parameter(iBin)
        for r in muValues:
            if self.drawBestfit:
                yieldParameter.setVal( r * self.yieldParameterDefaultValues[ self.yieldParameterNames.index(yieldParameter.GetName()) ] )
            else:
                yieldParameter.setVal(r)
            sigPdf = self.w.pdf( 'pdf_bin{0}'.format(binStr) )
            Hsig = sigPdf.createHistogram( 'sig_'+binStr, mH )
            Hsig.Scale( self.binWidth ) # default is /GeV
            if self.drawBestfit and r != 1.:
                Hsig.SetLineStyle(2)
            Hsigs.append( Hsig )

        dataset = self.w.data('data_obs')
        ROOT.SetOwnership( dataset, False )
        Hdata = get_hist_from_RooDataSet( dataset, mH, self._category )


        res.do_hgg              = self.do_hgg
        res.do_hzz              = self.do_hzz

        res.iBin                = iBin
        res.iCat                = iCat
        res.hggBin              = hggBin
        res.hggCat              = hggCat
        res.hzzBin              = hzzBin
        res.hzzCat              = hzzCat

        res.muValues            = muValues

        res.yieldParameterName  = yieldParameter.GetName()
        res.binStr              = binStr
        res.Hbkg                = Hbkg
        res.Hdata               = Hdata
        res.Hsigs               = Hsigs

        return res



    # def get_shapes_hzz(
    #         self,
    #         hzzBin   = None,
    #         hzzCat   = None,
    #         muValues = [ 0.5, 1.0, 2.0 ],
    #         ):

    #     res = Container()


    #     binStr = self.get_bin_str( hzzBin=hzzBin, hzzCat=hzzCat )

    #     # Set the proper category
    #     self.category = binStr

    #     # Get bkg and sig pdfs
    #     self.mH_hzz.setRange( *self.default_mH_hzz_range )
    #     self.mH_hzz.setBins( int((self.default_mH_hzz_range[1]-self.default_mH_hzz_range[0])/self.binWidth) )

    #     bkgPdf = self.w.pdf( 'pdf_bin{0}_bonly'.format(binStr) )
    #     Hbkg = bkgPdf.createHistogram( 'bkg_'+binStr, self.mH_hzz )
    #     Hbkg.Scale( self.binWidth ) # default is /GeV

    #     Hsigs = []
    #     yieldParameter = self.get_yield_parameter(hzzBin)
    #     for r in muValues:
    #         yieldParameter.setVal(r)
    #         sigPdf = self.w.pdf( 'pdf_bin{0}'.format(binStr) )
    #         Hsig = sigPdf.createHistogram( 'sig_'+binStr, self.mH_hzz )
    #         Hsig.Scale( self.binWidth ) # default is /GeV
    #         Hsigs.append( Hsig )


    #     dataset = self.w.data('data_obs')

    #     ROOT.SetOwnership( dataset, False )

    #     Hdata = get_hist_from_RooDataSet( dataset, self.mH_hzz, self._category )

    #     res.hzzBin              = hzzBin
    #     res.hzzCat              = hzzCat
    #     res.muValues            = muValues

    #     res.yieldParameterName  = yieldParameter.GetName()
    #     res.binStr              = binStr
    #     res.Hbkg                = Hbkg
    #     res.Hdata               = Hdata
    #     res.Hsigs               = Hsigs

    #     return res





    def draw_shape_result( self, res ):

        c.Clear()
        set_cmargins( RightMargin=0.05 )

        xMin = 120.
        xMax = 130.

        # yMinAbs = min( get_TH1_minimum_in_range( Hbkg, xMin, xMax ), get_TH1_minimum_in_range( Hsig, xMin, xMax ) )
        yMinAbs = 0.
        # yMaxAbs = max(
        #     get_TH1_maximum_in_range( res.Hbkg, xMin, xMax ),
        #     max([ get_TH1_maximum_in_range( Hsig, xMin, xMax ) for Hsig in res.Hsigs ])
        #     )

        if res.do_hgg:
            yMaxAbs = res.Hbkg.GetBinContent( res.Hbkg.FindBin(xMin) )
    
            yMax = yMaxAbs + 0.6*(yMaxAbs-yMinAbs)

        elif res.do_hzz:
            yMaxAbs = get_TGraph_maximum_in_range( res.Hdata, xMin, xMax )
            yMax = yMaxAbs + 1
            if yMax < 3.: yMax = 3.

        yMin = yMinAbs # - 0.6*(yMaxAbs-yMinAbs)


        base = get_plot_base(
            xMin = xMin,
            xMax = xMax,
            yMin = yMin,
            yMax = yMax,
            xTitle = 'm_{#gamma#gamma}' if res.do_hgg else 'm_{4l}',
            yTitle = '# of events / {0:.2f} GeV'.format(self.binWidth),
            )
        base.Draw('P')


        res.Hdata.Draw('PSAME')

        res.Hbkg.SetLineColor(4)
        res.Hbkg.SetLineWidth(2)
        res.Hbkg.Draw('C SAME')

        for Hsig in res.Hsigs:
            Hsig.SetLineWidth(2)
            Hsig.SetLineColor(2)
            Hsig.Draw('C SAME')


        l = ROOT.TLatex()
        l.SetNDC()
        l.SetTextAlign(13)
        l.SetTextAngle(90)
        l.DrawLatex( 1-c.GetRightMargin()+0.01, c.GetBottomMargin()+0.02,
            ( '#sigma_{{m}}/m cat {0}'.format(res.hggCat) if res.do_hgg else
              'H #rightarrow {0}'.format( res.iCat.replace('mu','#mu') ) )
            )


        binBoundaries = self.hgg_binBoundaries if res.do_hgg else self.hzz_binBoundaries
        if res.iBin < len(binBoundaries)-1:
            binStr = 'p_{{T}} [ {0}, {1} ] GeV'.format( int(binBoundaries[res.iBin]), int(binBoundaries[res.iBin+1]) )
        else:
            binStr = 'p_{{T}} >{0} GeV'.format( int(binBoundaries[res.iBin]) )

        l.SetTextAngle(0)
        l.SetTextAlign(23)
        l.DrawLatex(
            0.5*( c.GetLeftMargin() + 1-c.GetRightMargin() ),
            1-c.GetTopMargin()-0.01,
            binStr
            )


        save_c( res.binStr )




def get_hist_from_RooDataSet( dataset, xVar, category ):

    ctemp = ROOT.TCanvas( 'ctemp', 'ctemp', 1000, 800 )
    ctemp.cd()
    ctemp.Clear()

    frame = xVar.frame()

    reduceStr = '{0}=={1}'.format( category.GetName(), category.getIndex() )

    dataset_reduced = dataset.reduce( reduceStr )
    dataset_reduced.plotOn( frame )
    frame.Draw()

    l = ctemp.GetListOfPrimitives()
    for i in xrange(l.GetEntries()):
        if isinstance( l.At(i), ROOT.RooHist ):
            H = l.At(i)
            break
    else:
        Commands.throw_error( 'ERROR: did not find a histogram', throwException=True )

    Hcopy = ROOT.RooHist( H )
    ROOT.SetOwnership( Hcopy, False )

    Hcopy.SetName( TheoryCommands.get_unique_root_name() )

    # ctemp.SaveAs( 'plots_{0}_onetimeplots/roodatasetplottest.pdf'.format(datestr) )
    del ctemp
    del frame

    c.cd()

    return Hcopy



def get_TH1_extremum_in_range( H, xMin, xMax, findMinimum ):
    
    if findMinimum:
        res = 10e12
    else:
        res = -10e12

    iMinimum = -99

    for iBin in xrange( H.GetNbinsX() ):
        x = H.GetBinCenter(iBin+1)
        if x < xMin or x > xMax: continue
        y = H.GetBinContent(iBin+1)
        if findMinimum:
            if y < res:
                res = y
                iMinimum = iBin
        else:
            if y > res:
                res = y
                iMinimum = iBin
    return res

def get_TH1_minimum_in_range( H, xMin, xMax ):
    return get_TH1_extremum_in_range( H, xMin, xMax, True )
def get_TH1_maximum_in_range( H, xMin, xMax ):
    return get_TH1_extremum_in_range( H, xMin, xMax, False )


def get_TGraph_extremum_in_range( Tg, xMin, xMax, findMinimum ):

    if findMinimum:
        res = 10e12
    else:
        res = -10e12

    iMinimum = -99
    
    x_Double = ROOT.Double(0)
    y_Double = ROOT.Double(0)
    for iPoint in xrange( Tg.GetN() ):
        Tg.GetPoint( iPoint, x_Double, y_Double )
        x = float(x_Double)
        y = float(y_Double)
        if x < xMin or x > xMax: continue

        if findMinimum:
            if y < res:
                res = y
                iMinimum = iPoint
        else:
            if y > res:
                res = y
                iMinimum = iPoint

    return res

def get_TGraph_minimum_in_range( Tg, xMin, xMax ):
    return get_TGraph_extremum_in_range( Tg, xMin, xMax, True )
def get_TGraph_maximum_in_range( Tg, xMin, xMax ):
    return get_TGraph_extremum_in_range( Tg, xMin, xMax, False )



