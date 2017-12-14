#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, itertools, operator, re, argparse, sys, random
from math import isnan, isinf
from os.path import *
from glob import glob
from copy import deepcopy

import LatestPaths

sys.path.append('src')
import Commands
import PhysicsCommands
import OneOfCommands
import TheoryCommands
import CorrelationMatrices
import MergeHGGWDatacards
import TheoryFileInterface
from Container import Container
from Parametrization import Parametrization, WSParametrization
import MuShapeDrawer
import PlotCommands

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

def AppendParserOptions( parser ):

    parser.add_argument( '--onetimeplotsCommands', default=False )
    class CustomAction(argparse.Action):
        def __init__(self, option_strings, dest, **kwargs):
            # Only 1 argument allowed (this is basically the "store_true" action)
            super(CustomAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            # Make sure the parser knows a combineCommand was entered
            setattr( namespace, 'onetimeplotsCommands', True )
            setattr( namespace, self.dest, True )
            # super(FooAction, self).__call__( parser, namespace, values, option_string=None )
            # setattr( namespace, self.dest, values )

    parser.add_argument( '--derivedTheoryFilePlots',                  action=CustomAction )
    parser.add_argument( '--FittingMuIllustration',                   action=CustomAction )
    parser.add_argument( '--DrawAllFits',                             action=CustomAction )


########################################
# Methods
########################################    


def GetTH1FExtremumInRange( H, xMin, xMax, findMinimum ):
    
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

def GetTH1FMinimumInRange( H, xMin, xMax ):
    return GetTH1FExtremumInRange( H, xMin, xMax, True )
def GetTH1FMaximumInRange( H, xMin, xMax ):
    return GetTH1FExtremumInRange( H, xMin, xMax, False )


def GetHistFromRooDataSet( dataset, xVar, category ):

    ctemp = ROOT.TCanvas( 'ctemp', 'ctemp', 1000, 800 )
    ctemp.cd()
    ctemp.Clear()

    frame = xVar.frame()

    dataset_reduced = dataset.reduce( '{0}=={1}'.format( category.GetName(), category.getIndex() ) )
    dataset_reduced.plotOn( frame )
    frame.Draw()

    l = ctemp.GetListOfPrimitives()
    for i in xrange(l.GetEntries()):
        if isinstance( l.At(i), ROOT.RooHist ):
            H = l.At(i)
            break
    else:
        print 'ERROR: did not find a histogram'

    Hcopy = ROOT.RooHist( H )
    ROOT.SetOwnership( Hcopy, False )

    Hcopy.SetName( TheoryCommands.GetUniqueRootName() )

    ctemp.SaveAs( 'plots_{0}_onetimeplots/roodatasetplottest.pdf'.format(datestr) )
    del ctemp
    del frame

    c.cd()

    return Hcopy



def main( args ):

    TheoryCommands.SetPlotDir( 'plots_{0}_onetimeplots'.format(datestr) )


    #____________________________________________________________________
    if args.DrawAllFits:

        TheoryCommands.SetPlotDir( 'plots_{0}_muShapes'.format(datestr) )

        # postfitFilename = 'corrMat_Oct17/higgsCombine_POSTFIT_combinedCard_Aug21.MultiDimFit.mH125.root'
        postfitFilename = 'corrMat_Oct17/higgsCombine_POSTFIT_combinedCard_Jul26.MultiDimFit.mH125.root'

        muShapeDrawer = MuShapeDrawer.MuShapeDrawer( postfitFilename )

        for hggCat in xrange(3):
            for hggBin in xrange(8):

                muShapeDrawer.DrawShapes(
                    hggBin   = hggBin,
                    hggCat   = hggCat,
                    muValues = [ 0.5, 1.0, 2.0 ],
                    )

        for hzzCat in [ '4e', '4mu', '2e2mu' ]:
            for hzzBin in xrange(5):

                muShapeDrawer.DrawShapes(
                    hzzBin   = hzzBin,
                    hzzCat   = hzzCat,
                    muValues = [ 0.5, 1.0, 2.0 ],
                    drawBestfit = True,
                    )


    #____________________________________________________________________
    if args.FittingMuIllustration:
        TheoryCommands.SetPlotDir( 'plots_{0}_muShapes'.format(datestr) )

        # postfitFilename = 'corrMat_Oct17/higgsCombine_POSTFIT_combinedCard_Aug21.MultiDimFit.mH125.root'
        postfitFilename = 'corrMat_Oct17/higgsCombine_POSTFIT_combinedCard_Jul26.MultiDimFit.mH125.root'

        postfitFp = ROOT.TFile.Open( postfitFilename )
        w = postfitFp.Get('w')
        w.loadSnapshot('MultiDimFit')

        mH = w.var('CMS_hgg_mass')


        # PDFs

        # mH.setBins( mH.getBinning().numBins() * 1 ) # Set per 0.25 GeV - default
        mH.setMin(100.)
        mH.setMax(180.)

        # Get the best fit bkg
        bkgPdfName = 'pdf_binch1_recoPt_0p0_15p0_SigmaMpTTag_0_13TeV_bonly'
        bkgPdf = w.pdf(bkgPdfName)
        Hbkg = bkgPdf.createHistogram( 'bkg_PTH_0_15_cat0', mH )
        Hbkg.Scale(0.25)

        Hsigs = []
        yieldParameter = w.var('r_smH_PTH_0_15')
        for r in [ 0.5, 1.0, 2.0 ]:
            yieldParameter.setVal(r)

            # Signal shape
            sigPdfName = 'pdf_binch1_recoPt_0p0_15p0_SigmaMpTTag_0_13TeV'
            sigPdf = w.pdf(sigPdfName)
            Hsig = sigPdf.createHistogram( 'sig_PTH_0_15_cat0', mH )
            Hsig.Scale(0.25)

            Hsigs.append( Hsig )


        # Get the data histogram
        dataset = w.data('data_obs')

        print dataset

        CMS_channel = w.cat('CMS_channel')
        CMS_channel.setIndex(0)
        CMS_channel.Print()

        mH.setRange( 120., 130. )
        mH.setBins( 40 )

        Hdata = GetHistFromRooDataSet( dataset, mH, CMS_channel )

        # This works, but the errors aren't right and getting them so is a pain
        # dataset_thisCat = dataset.reduce( '{0}=={1}'.format( CMS_channel.GetName(), CMS_channel.getIndex() ) )
        # Hdata = ROOT.RooAbsData.createHistogram(
        #     dataset_thisCat, 'data_PTH_0_15_cat0', mH
        #     )

        c.Clear()
        SetCMargins()

        xMin = 120.
        xMax = 130.

        # yMinAbs = min( GetTH1FMinimumInRange( Hbkg, xMin, xMax ), GetTH1FMinimumInRange( Hsig, xMin, xMax ) )
        yMinAbs = 0.
        yMaxAbs = max( GetTH1FMaximumInRange( Hbkg, xMin, xMax ), GetTH1FMaximumInRange( Hsig, xMin, xMax ) )
        yMin = yMinAbs # - 0.6*(yMaxAbs-yMinAbs)
        yMax = yMaxAbs + 0.4*(yMaxAbs-yMinAbs)


        base = GetPlotBase(
            xMin = xMin,
            xMax = xMax,
            yMin = yMin,
            yMax = yMax,
            xTitle = 'm_{#gamma#gamma}',
            yTitle = '# of events / 0.25 GeV',
            )
        base.Draw('P')

        Hbkg.SetLineColor(4)

        Hdata.Draw('PSAME')
        Hbkg.Draw('C SAME')

        for Hsig in Hsigs:
            Hsig.SetLineColor(2)
            Hsig.Draw('C SAME')

        SaveC( 'muFitIllustration1' )
        postfitFp.Close()



    #____________________________________________________________________
    if args.derivedTheoryFilePlots:

        DO_PARAMETRIZATION = True
        # DO_PARAMETRIZATION = False

        PLOT_CROSSSECTIONS = True
        # PLOT_CROSSSECTIONS = False


        # ======================================
        # Input

        readList = lambda theoryFiles: [ TheoryFileInterface.ReadDerivedTheoryFile( theoryFile ) for theoryFile in theoryFiles ]


        # ----------------------------
        TheoryFileInterface.SetFileFinderDir( LatestPaths.derivedTheoryFiles_Top )

        theoryFiles_ct_cg = TheoryFileInterface.FileFinder(
            ct='*', cg='*', muR=1, muF=1, Q=1, filter='cb'
            )
        containers_ct_cg = readList( theoryFiles_ct_cg )
        SM_ct_cg = TheoryFileInterface.ReadDerivedTheoryFile( TheoryFileInterface.FileFinder( ct=1, cg=0, cb=1, muR=1, muF=1, Q=1, expectOneFile=True ) )

        theoryFiles_ct_cb = TheoryFileInterface.FileFinder(
            ct='*', cb='*', muR=1, muF=1, Q=1, filter='cg'
            )
        containers_ct_cb = readList( theoryFiles_ct_cb )
        SM_ct_cb = SM_ct_cg # Same SM

        # ----------------------------
        containers_ct_cg_HighPt = TheoryFileInterface.FileFinder(
            ct='*', cg='*', cb=1, muR=1, muF=1, Q=1, filter='ct_1_cg_0', loadImmediately=True,
            directory = LatestPaths.derivedTheoryFiles_TopHighPt
            )
        SM_ct_cg_HighPt = TheoryFileInterface.FileFinder(
            ct=1, cg=0, cb=1, muR=1, muF=1, Q=1, expectOneFile=True, loadImmediately=True,
            directory = LatestPaths.derivedTheoryFiles_TopHighPt
            )

        # ----------------------------
        containers_ct_cb_HighPt = TheoryFileInterface.FileFinder(
            ct='*', cb='*', cg=0, muR=1, muF=1, Q=1, filter='ct_1_cg_0', loadImmediately=True,
            directory = LatestPaths.derivedTheoryFiles_TopHighPt
            )            
        SM_ct_cb_HighPt = SM_ct_cg_HighPt # Same SM


        # ----------------------------
        theoryFiles_kappab_kappac = TheoryFileInterface.FileFinder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFiles_YukawaSummed
            )
        containers_kappab_kappac = readList( theoryFiles_kappab_kappac )
        SM_kappab_kappac = [ container for container in containers_kappab_kappac if container.kappab == 1. and container.kappac == 1. ][0]

        theoryFiles_kappab_kappac_Gluon = TheoryFileInterface.FileFinder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFiles_YukawaGluonInduced
            )
        containers_kappab_kappac_Gluon = readList( theoryFiles_kappab_kappac_Gluon )
        SM_kappab_kappac_Gluon = [ container for container in containers_kappab_kappac_Gluon if container.kappab == 1. and container.kappac == 1. ][0]


        theoryFiles_kappab_kappac_Quark = TheoryFileInterface.FileFinder(
            kappab='*', kappac='*', filter='muR',
            directory = LatestPaths.derivedTheoryFiles_YukawaQuarkInduced
            )
        containers_kappab_kappac_Quark = readList( theoryFiles_kappab_kappac_Quark )
        SM_kappab_kappac_Quark = [ container for container in containers_kappab_kappac_Quark if container.kappab == 1. and container.kappac == 1. ][0]


        theoryFiles_kappab_kappac_QuarkScaled = TheoryFileInterface.FileFinder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFiles_YukawaQuarkInducedScaled
            )
        containers_kappab_kappac_QuarkScaled = readList( theoryFiles_kappab_kappac_QuarkScaled )
        SM_kappab_kappac_QuarkScaled = [ container for container in containers_kappab_kappac_QuarkScaled if container.kappab == 1. and container.kappac == 1. ][0]


        def kappaSorter( container ):
            if hasattr( container, 'kappab' ) and hasattr( container, 'kappac' ):
                return ( container.kappab, container.kappac )
            elif hasattr( container, 'ct' ) and hasattr( container, 'cg' ):
                return ( container.ct, container.cg )
            elif hasattr( container, 'ct' ) and hasattr( container, 'cb' ):
                return ( container.ct, container.cb )
            else:
                Commands.ThrowError( 'Can\'t sort; not a known set of couplings found' )

        containers_kappab_kappac.sort( key=kappaSorter )
        containers_kappab_kappac_Gluon.sort( key=kappaSorter )
        containers_kappab_kappac_Quark.sort( key=kappaSorter )
        containers_kappab_kappac_QuarkScaled.sort( key=kappaSorter )
        containers_ct_cg.sort( key=kappaSorter )
        containers_ct_cb.sort( key=kappaSorter )
        containers_ct_cg_HighPt.sort( key=kappaSorter )


        ct_cg = Container()
        ct_cg.name = 'kappat_kappag'
        ct_cg.couplings = [ 'ct', 'cg' ]
        ct_cg.containers = containers_ct_cg
        ct_cg.ws = LatestPaths.ws_combined_Top
        ct_cg.SM = SM_ct_cg

        ct_cg_HighPt = Container()
        ct_cg_HighPt.name = 'kappat_kappag_HighPt'
        ct_cg_HighPt.couplings = [ 'ct', 'cg' ]
        ct_cg_HighPt.containers = containers_ct_cg_HighPt
        ct_cg_HighPt.ws = LatestPaths.ws_combined_TopHighPt
        # ct_cg_HighPt.parametrizeByFitting = True
        # ct_cg_HighPt.includeLinearTerms = False
        ct_cg_HighPt.SM = SM_ct_cg_HighPt

        ct_cb = Container()
        ct_cb.name = 'kappat_kappab'
        ct_cb.couplings = [ 'ct', 'cb' ]
        ct_cb.containers = containers_ct_cb
        ct_cb.ws = LatestPaths.ws_combined_TopCtCb
        ct_cb.SM = SM_ct_cb

        ct_cb_HighPt = Container()
        ct_cb_HighPt.name = 'kappat_kappab_HighPt'
        ct_cb_HighPt.couplings = [ 'ct', 'cb' ]
        ct_cb_HighPt.containers = containers_ct_cb_HighPt
        ct_cb_HighPt.ws = LatestPaths.ws_combined_TopCtCbHighPt
        ct_cb_HighPt.SM = SM_ct_cb_HighPt

        kappab_kappac_Gluon = Container()
        kappab_kappac_Gluon.name = 'kappab_kappac_Gluon'
        kappab_kappac_Gluon.couplings = [ 'kappab', 'kappac' ]
        kappab_kappac_Gluon.containers = containers_kappab_kappac_Gluon
        kappab_kappac_Gluon.parametrizeByFitting = True
        kappab_kappac_Gluon.SM = SM_kappab_kappac_Gluon

        kappab_kappac_Quark = Container()
        kappab_kappac_Quark.name = 'kappab_kappac_Quark'
        kappab_kappac_Quark.couplings = [ 'kappab', 'kappac' ]
        kappab_kappac_Quark.containers = containers_kappab_kappac_Quark
        kappab_kappac_Quark.parametrizeByFitting = True
        kappab_kappac_Quark.SM = SM_kappab_kappac_Quark

        kappab_kappac_QuarkScaled = Container()
        kappab_kappac_QuarkScaled.name = 'kappab_kappac_QuarkScaled'
        kappab_kappac_QuarkScaled.couplings = [ 'kappab', 'kappac' ]
        kappab_kappac_QuarkScaled.containers = containers_kappab_kappac_QuarkScaled
        kappab_kappac_QuarkScaled.SM = SM_kappab_kappac_QuarkScaled
        kappab_kappac_QuarkScaled.DrawFileContentsAsLine = True


        couplingVariationContainers = [
            # ct_cg,
            ct_cb,
            # ct_cg_HighPt,
            ct_cb_HighPt,
            # kappab_kappac_Gluon,
            # kappab_kappac_Quark,
            # kappab_kappac_QuarkScaled,
            ]

        # n_theoryFiles_kappab_kappac = len( theoryFiles_kappab_kappac )
        # for iChunk, chunk in enumerate( chunks( containers_kappab_kappac, n_theoryFiles_kappab_kappac/4+1 ) ):
        #     couplingVariationContainer = Container()
        #     couplingVariationContainer.name = 'kappab_kappac_{0}'.format(iChunk)
        #     couplingVariationContainer.couplings = [ 'kappab', 'kappac' ]
        #     couplingVariationContainer.containers = list(chunk)
        #     # couplingVariationContainer.ws = LatestPaths.ws_combined_split_yukawa
        #     couplingVariationContainer.ws = LatestPaths.ws_combined_Yukawa
        #     # couplingVariationContainers.append( couplingVariationContainer )
        #     couplingVariationContainer.SM = SM_kappab_kappac

        for iVal, kappabVal in enumerate([ -2.0, -1.0, 0.0, 1.0, 2.0 ]):
            couplingVariationContainer = Container()
            couplingVariationContainer.name = 'kappab_kappac_{0}'.format(iVal)
            couplingVariationContainer.couplings = [ 'kappab', 'kappac' ]
            couplingVariationContainer.containers = [ container for container in containers_kappab_kappac if container.kappab == kappabVal ]
            couplingVariationContainer.ws = LatestPaths.ws_combined_Yukawa
            couplingVariationContainer.SM = SM_kappab_kappac
            # couplingVariationContainers.append( couplingVariationContainer )



    if args.derivedTheoryFilePlots:

        DO_PARAMETRIZATION = True

        TOP_PANEL_LOGSCALE = True

        SKIP_KAPPAB_0_KAPPAC_1 = True
        # SKIP_KAPPAB_0_KAPPAC_1 = False


        def SetParametrizedTgs( couplingVariationContainer ):

            if DO_PARAMETRIZATION and hasattr( couplingVariationContainer, 'ws' ):
                wsParametrization = WSParametrization( couplingVariationContainer.ws, verbose=True )

                for container in couplingVariationContainer.containers:
                    # for coupling in couplingVariationContainer.couplings:
                    #     print 'Setting \'{0}\' to \'{1}\''.format( coupling, float(getattr( container, coupling )) )
                    #     setattr( wsParametrization, coupling, getattr( container, coupling ) )
                    kwargs = { coupling : float(getattr( container, coupling )) for coupling in couplingVariationContainer.couplings }
                    kwargs['returnWhat'] = 'theory'

                    outputContainer = wsParametrization.GetOutputContainer( **kwargs )
                    container.Tg_parametrization_ratio = outputContainer.Tg

                    outputContainer.crosssection = [ SMxs * ratio for SMxs, ratio in zip( couplingVariationContainer.SM.crosssection, outputContainer.mus ) ]
                    container.Tg_parametrization_crosssection = outputContainer.GetTGraph( yAttr='crosssection' )

                couplingVariationContainer.isParametrized = True


            elif DO_PARAMETRIZATION and hasattr( couplingVariationContainer, 'parametrizeByFitting' ):
                parametrization = Parametrization()
                parametrization.ParametrizeByFitting(
                    couplingVariationContainer.containers,
                    fitWithScipy=True,
                    includeLinearTerms = getattr( couplingVariationContainer, 'includeLinearTerms', True ),
                    couplingsToParametrize = couplingVariationContainer.couplings
                    )

                for container in couplingVariationContainer.containers:
                    # kwargs = { coupling : float(getattr( container, coupling )) for coupling in couplingVariationContainer.couplings }
                    for coupling in couplingVariationContainer.couplings:
                        setattr( parametrization, coupling, float(getattr( container, coupling )) )
                    parametrizationResult = parametrization.GetOutputContainer()
                    parametrizationResult.ratios = [ xs / SMxs if not SMxs == 0. else 0. for xs, SMxs in zip( parametrizationResult.binValues, couplingVariationContainer.SM.crosssection ) ]

                    container.Tg_parametrization_crosssection = parametrizationResult.GetTGraph( yAttr = 'crosssection' )
                    container.Tg_parametrization_ratio = parametrizationResult.GetTGraph( yAttr = 'ratios' )

                couplingVariationContainer.isParametrized = True

            else:
                couplingVariationContainer.isParametrized = False

        def CouplingTitle( coupling ):
            if re.match( r'c\w', coupling ):
                particle = coupling.replace('c','')
            elif re.match( r'kappa\w', coupling ):
                particle = coupling.replace('kappa','')
            return '#kappa_{{{0}}}'.format( particle )


        newColorCycle = lambda : itertools.cycle( range(2,10) + range(27,50) )



        # ======================================
        # Loop

        for DO_PARAMETRIZATION in [ False, True ]:
            for couplingVariationContainer in couplingVariationContainers:
                SetParametrizedTgs( couplingVariationContainer )

                if SKIP_KAPPAB_0_KAPPAC_1:
                    containers = [ container for container in couplingVariationContainer.containers if not ( getattr( container, 'kappab', 999. ) == 0.0 and getattr( container, 'kappac', 999. ) == 1.0 ) ]
                else:
                    containers = couplingVariationContainer.containers

                c.Clear()


                # ======================================
                # Gather panel objects

                topPanelObjects = []
                bottomPanelObjects = []


                # -----------------------
                # Make base for top panel

                yMinAbs = min( [ min(container.crosssection) for container in containers ] )
                yMaxAbs = max( [ max(container.crosssection) for container in containers ] )
                yMin = yMinAbs - 0.1*(yMaxAbs-yMinAbs)
                yMax = yMaxAbs + 0.5*(yMaxAbs-yMinAbs)

                if TOP_PANEL_LOGSCALE:
                    yMin = max( 0.9*yMinAbs, 0.0001 )
                    yMax = 20*yMaxAbs

                xMin = 0.
                xMax = max( containers[0].binBoundaries )

                baseTop = GetPlotBase(
                    xMin = xMin,
                    xMax = xMax,
                    yMin = yMin,
                    yMax = yMax,
                    # xTitle = 'p_{T}^{H}',
                    xTitle = '',
                    yTitle = '#sigma_{ggH} (pb/GeV)',
                    )
                topPanelObjects.append( ( baseTop, 'P' ) )
                baseTop.SetTitleOffset(1.0)


                leg = PlotCommands.TLegendMultiPanel(
                    lambda c: c.GetLeftMargin() + 0.01,
                    lambda c: 1 - c.GetTopMargin() - 0.23,
                    lambda c: 1 - c.GetRightMargin() - 0.01,
                    lambda c: 1 - c.GetTopMargin()
                    )

                leg.SetNColumns(2)
                if len(containers) > 10:
                    leg.SetNColumns(3)
                if len(containers) > 20:
                    leg.SetNColumns(4)
                leg.SetBorderSize(0)
                leg.SetFillStyle(0)


                # -----------------------
                # Make base for bottom panel

                yMinAbs = min( [ min(container.ratios) for container in containers ] )
                yMaxAbs = max( [ max(container.ratios) for container in containers ] )
                yMin = yMinAbs - 0.1*(yMaxAbs-yMinAbs)
                yMax = yMaxAbs + 0.1*(yMaxAbs-yMinAbs)

                xMin = 0.
                xMax = max( containers[0].binBoundaries )

                baseBottom = GetPlotBase(
                    xMin = xMin,
                    xMax = xMax,
                    yMin = yMin,
                    yMax = yMax,
                    xTitle = 'p_{T}^{H}',
                    yTitle = 'ratio',
                    )
                bottomPanelObjects.append( ( baseBottom, 'P' ) )
                baseBottom.SetTitleOffset(1.0)

                lineAtOne = ROOT.TLine( xMin, 1.0, xMax, 1.0 )
                lineAtOne.SetLineColor(14)
                bottomPanelObjects.append( ( lineAtOne, '' ) )


                # -----------------------
                # Fill in objects

                colorCycle = newColorCycle()
                for container in containers:
                    color = next(colorCycle)
                    container.color = color
                    
                    variationName  = []
                    variationTitle = []
                    for coupling in [ 'ct', 'cg', 'cb', 'kappab', 'kappac' ]:
                        if hasattr( container, coupling ):
                            variationName.append( '{0}_{1}'.format( coupling, Commands.ConvertFloatToStr( getattr( container, coupling ) ) ))
                            variationTitle.append( '{0} = {1}'.format( CouplingTitle(coupling), getattr( container, coupling ) ))
                    variationName  = '_'.join(variationName)
                    variationTitle = ', '.join(variationTitle)

                    Tg_xs = TheoryFileInterface.ReadDerivedTheoryContainerToTGraph(
                        container,
                        name = variationName,
                        yAttr = 'crosssection'
                        )
                    Tg_xs.SetLineColor(color)
                    Tg_xs.SetLineWidth(2)
                    Tg_xs.SetMarkerColor(color)
                    Tg_xs.SetMarkerStyle(8)
                    Tg_xs.SetMarkerSize(0.9)

                    Tg_ratio = TheoryFileInterface.ReadDerivedTheoryContainerToTGraph(
                        container,
                        name = variationName + '_ratios',
                        yAttr = 'ratios'
                        )
                    Tg_ratio.SetLineColor(color)
                    Tg_ratio.SetLineWidth(2)
                    Tg_ratio.SetMarkerColor(color)
                    Tg_ratio.SetMarkerStyle(8)
                    Tg_ratio.SetMarkerSize(0.9)


                    if hasattr( couplingVariationContainer, 'DrawFileContentsAsLine' ) and couplingVariationContainer.DrawFileContentsAsLine:
                        leg.AddEntry( Tg_xs.GetName(), variationTitle, 'l' )
                        topPanelObjects.append( ( Tg_xs, 'LX' ) )
                        bottomPanelObjects.append( ( Tg_ratio, 'LX' ) )
                    else:
                        topPanelObjects.append( ( Tg_xs, 'PX' ) )
                        bottomPanelObjects.append( ( Tg_ratio, 'PX' ) )
                        leg.AddEntry( Tg_xs.GetName(), variationTitle, 'p' )

                    if couplingVariationContainer.isParametrized:
                        container.Tg_parametrization_crosssection.SetLineColor(color)
                        container.Tg_parametrization_crosssection.SetLineWidth(2)
                        container.Tg_parametrization_ratio.SetLineColor(color)
                        container.Tg_parametrization_ratio.SetLineWidth(2)
                        topPanelObjects.append( ( container.Tg_parametrization_crosssection, 'LX' ) )
                        bottomPanelObjects.append( ( container.Tg_parametrization_ratio, 'LX' ) )


                topPanelObjects.append( ( leg, '' ) )

                PlotCommands.PlotWithBottomPanel(
                    couplingVariationContainer.name + '_twoPanel' + ( '_parametrized' if DO_PARAMETRIZATION else '' ),
                    topPanelObjects,
                    bottomPanelObjects,
                    xTitle       = 'p_{T}^{H}',
                    yTitleTop    = '#sigma_{ggH} (pb/GeV)',
                    yTitleBottom = 'ratio w.r.t. SM',
                    SetTopPanelLogScale = TOP_PANEL_LOGSCALE,
                    topPadLeftMargin = 0.14,
                    bottomPadLeftMargin = 0.14,
                    disableCMSText = True,
                    )



    if False and args.derivedTheoryFilePlots:
        # This is the old style, only doing xs or ratio at one go

        # ======================================
        # Loop

        for couplingVariationContainer in couplingVariationContainers:
            colorCycle = itertools.cycle( range(2,10) + range(27,50) )

            containers = couplingVariationContainer.containers
            name       = couplingVariationContainer.name

            # Try to get the ws parametrizations
            isParametrized = False

            if DO_PARAMETRIZATION and hasattr( couplingVariationContainer, 'ws' ):
                wsParametrization = WSParametrization( couplingVariationContainer.ws, verbose=True )

                for container in containers:
                    # for coupling in couplingVariationContainer.couplings:
                    #     print 'Setting \'{0}\' to \'{1}\''.format( coupling, float(getattr( container, coupling )) )
                    #     setattr( wsParametrization, coupling, getattr( container, coupling ) )
                    kwargs = { coupling : float(getattr( container, coupling )) for coupling in couplingVariationContainer.couplings }
                    kwargs['returnWhat'] = 'theory'

                    outputContainer = wsParametrization.GetOutputContainer( **kwargs )
                    container.Tg_parametrization = outputContainer.Tg

                    if PLOT_CROSSSECTIONS:
                        outputContainer.crosssection = [ SMxs * ratio for SMxs, ratio in zip( couplingVariationContainer.SM.crosssection, outputContainer.mus ) ]
                        container.Tg_parametrization = outputContainer.GetTGraph( yAttr='crosssection' )

                isParametrized = True


            elif DO_PARAMETRIZATION and hasattr( couplingVariationContainer, 'parametrizeByFitting' ):

                parametrization = Parametrization()
                parametrization.ParametrizeByFitting( containers, fitWithScipy=True )

                for container in containers:
                    # kwargs = { coupling : float(getattr( container, coupling )) for coupling in couplingVariationContainer.couplings }
                    for coupling in couplingVariationContainer.couplings:
                        setattr( parametrization, coupling, float(getattr( container, coupling )) )
                    parametrizationResult = parametrization.GetOutputContainer()
                    parametrizationResult.ratios = [ xs / SMxs if not SMxs == 0. else 0. for xs, SMxs in zip( parametrizationResult.binValues, couplingVariationContainer.SM.crosssection ) ]
                    # container.Tg_parametrization = parametrizationResult.GetTGraph( yAttr='ratios' )
                    container.Tg_parametrization = parametrizationResult.GetTGraph( yAttr = 'crosssection' if PLOT_CROSSSECTIONS else 'ratios' )

                isParametrized = True


            c.Clear()
            SetCMargins()

            yMinAbs = min( [ min( container.crosssection if PLOT_CROSSSECTIONS else container.ratios ) for container in containers ] )
            yMaxAbs = max( [ max( container.crosssection if PLOT_CROSSSECTIONS else container.ratios ) for container in containers ] )
            yMin = yMinAbs - 0.1*(yMaxAbs-yMinAbs)
            yMax = yMaxAbs + 0.5*(yMaxAbs-yMinAbs)

            xMin = 0.
            xMax = max( containers[0].binBoundaries )

            base = GetPlotBase(
                xMin = xMin,
                xMax = xMax,
                yMin = yMin,
                yMax = yMax,
                xTitle = 'p_{T}^{H}',
                yTitle = '#sigma_{ggH} (pb/GeV)' if PLOT_CROSSSECTIONS else 'Ratio w.r.t. ggH^{SM}',
                )
            base.Draw('P')

            leg = ROOT.TLegend(
                c.GetLeftMargin() + 0.01,
                1 - c.GetTopMargin() - 0.23,
                1 - c.GetRightMargin() - 0.01,
                1 - c.GetTopMargin()
                )
            leg.SetNColumns(2)
            if len(containers) > 10:
                leg.SetNColumns(3)
            if len(containers) > 20:
                leg.SetNColumns(4)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)


            def c_to_kappa( coupling ):
                if re.match( r'c\w', coupling ):
                    particle = coupling.replace('c','')
                elif re.match( r'kappa\w', coupling ):
                    particle = coupling.replace('kappa','')
                return '#kappa_{{{0}}}'.format( particle )

            for container in containers:
                color = next(colorCycle)
                container.color = color
                
                variationName  = []
                variationTitle = []
                for coupling in [ 'ct', 'cg', 'cb', 'kappab', 'kappac' ]:
                    if hasattr( container, coupling ):
                        variationName.append( '{0}_{1}'.format(
                            coupling,
                            Commands.ConvertFloatToStr( getattr( container, coupling ) )
                            ))
                        variationTitle.append( '{0} = {1}'.format(
                            c_to_kappa(coupling),
                            getattr( container, coupling )
                            ))

                variationName  = '_'.join(variationName)
                variationTitle = ', '.join(variationTitle)


                Tg = TheoryFileInterface.ReadDerivedTheoryContainerToTGraph(
                    container,
                    name = variationName,
                    yAttr = 'crosssection' if PLOT_CROSSSECTIONS else 'ratios'
                    )
                Tg.SetLineColor(color)
                Tg.SetLineWidth(2)
                Tg.SetMarkerColor(color)
                Tg.SetMarkerStyle(8)
                Tg.SetMarkerSize(0.9)
                
                if hasattr( couplingVariationContainer, 'DrawFileContentsAsLine' ) and couplingVariationContainer.DrawFileContentsAsLine:
                    Tg.Draw('LX')
                    leg.AddEntry( Tg.GetName(), variationTitle, 'l' )
                else:
                    Tg.Draw('PX')
                    leg.AddEntry( Tg.GetName(), variationTitle, 'p' )


                if isParametrized:
                    container.Tg_parametrization.SetLineColor(color)
                    container.Tg_parametrization.SetLineWidth(2)
                    container.Tg_parametrization.Draw('LX')

                

            leg.Draw()


            if isParametrized:
                name += '_parametrized'

            SaveC( name )






def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'