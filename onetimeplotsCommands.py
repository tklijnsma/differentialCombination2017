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

from OptionHandler import flag_as_option, flag_as_parser_options

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

#____________________________________________________________________
@flag_as_option
def DrawAllFits(args):
    TheoryCommands.SetPlotDir( 'plots_{0}_muShapes'.format(datestr) )
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
@flag_as_option
def derivedTheoryFilePlots(args):
    Commands.Warning('Do yukawa with a python anaconda environment')
    Commands.Warning(
        '[Feb 14] Causes segmentation violation; '
        'Check if the set-naming that WSParametrization relies on '
        '(yieldParameter vs all_ggH_yieldParameters or something)'
        )
    functions = [
        # read_yukawa_full,
        # read_yukawa_Gluon,
        # read_yukawa_Quark,
        # read_yukawa_QuarkScaled,
        read_top_highpt,
        read_top_highpt_ctcb,
        # read_yukawa_chunked,
        ]
    for function in functions:
        list_of_containers = function()
        if isinstance(list_of_containers, Container):
            list_of_containers = [list_of_containers]

        for container in list_of_containers:
            theory_plot_container(container)


def readList(theoryFiles):
    return [ TheoryFileInterface.ReadDerivedTheoryFile( theoryFile ) for theoryFile in theoryFiles ]

def only_load_kappab_kappac():
    theoryFiles_kappab_kappac = TheoryFileInterface.FileFinder(
        kappab='*', kappac='*', muR=1, muF=1, Q=1,
        directory = LatestPaths.derivedTheoryFiles_YukawaSummed
        )
    containers_kappab_kappac = readList( theoryFiles_kappab_kappac )
    SM_kappab_kappac = [ container for container in containers_kappab_kappac if container.kappab == 1. and container.kappac == 1. ][0]
    containers_kappab_kappac.sort(key=kappaSorter)
    return containers_kappab_kappac, SM_kappab_kappac

def read_yukawa_chunked():
    containers_kappab_kappac, SM_kappab_kappac = only_load_kappab_kappac()
    couplingVariationContainers = []
    for iVal, kappabVal in enumerate([ -2.0, -1.0, 0.0, 1.0, 2.0 ]):
        couplingVariationContainer = Container()
        couplingVariationContainer.name = 'kappab_kappac_{0}'.format(iVal)
        couplingVariationContainer.couplings = [ 'kappab', 'kappac' ]
        couplingVariationContainer.containers = [ container for container in containers_kappab_kappac if container.kappab == kappabVal ]
        couplingVariationContainer.ws = LatestPaths.ws_combined_Yukawa
        couplingVariationContainer.SM = SM_kappab_kappac
        couplingVariationContainers.append( couplingVariationContainer )
    return couplingVariationContainers

def read_yukawa_full():
    containers_kappab_kappac, SM_kappab_kappac = only_load_kappab_kappac()
    kappab_kappac = Container()
    kappab_kappac.name = 'kappab_kappac'
    kappab_kappac.couplings = [ 'kappab', 'kappac' ]
    kappab_kappac.containers = containers_kappab_kappac
    kappab_kappac.parametrizeByFitting = True
    kappab_kappac.SM = SM_kappab_kappac
    return kappab_kappac

def read_yukawa_Gluon():
    theoryFiles_kappab_kappac_Gluon = TheoryFileInterface.FileFinder(
        kappab='*', kappac='*', muR=1, muF=1, Q=1,
        directory = LatestPaths.derivedTheoryFiles_YukawaGluonInduced
        )
    containers_kappab_kappac_Gluon = readList( theoryFiles_kappab_kappac_Gluon )
    SM_kappab_kappac_Gluon = [ container for container in containers_kappab_kappac_Gluon if container.kappab == 1. and container.kappac == 1. ][0]
    containers_kappab_kappac_Gluon.sort(key=kappaSorter)
    kappab_kappac_Gluon = Container()
    kappab_kappac_Gluon.name = 'kappab_kappac_Gluon'
    kappab_kappac_Gluon.couplings = [ 'kappab', 'kappac' ]
    kappab_kappac_Gluon.containers = containers_kappab_kappac_Gluon
    kappab_kappac_Gluon.parametrizeByFitting = True
    kappab_kappac_Gluon.SM = SM_kappab_kappac_Gluon
    return kappab_kappac_Gluon

def read_yukawa_Quark():
    theoryFiles_kappab_kappac_Quark = TheoryFileInterface.FileFinder(
        kappab='*', kappac='*', filter='muR',
        directory = LatestPaths.derivedTheoryFiles_YukawaQuarkInduced
        )
    containers_kappab_kappac_Quark = readList( theoryFiles_kappab_kappac_Quark )
    SM_kappab_kappac_Quark = [ container for container in containers_kappab_kappac_Quark if container.kappab == 1. and container.kappac == 1. ][0]
    containers_kappab_kappac_Quark.sort(key=kappaSorter)
    kappab_kappac_Quark = Container()
    kappab_kappac_Quark.name = 'kappab_kappac_Quark'
    kappab_kappac_Quark.couplings = [ 'kappab', 'kappac' ]
    kappab_kappac_Quark.containers = containers_kappab_kappac_Quark
    kappab_kappac_Quark.parametrizeByFitting = True
    kappab_kappac_Quark.SM = SM_kappab_kappac_Quark
    return kappab_kappac_Quark

def read_yukawa_QuarkScaled():
    theoryFiles_kappab_kappac_QuarkScaled = TheoryFileInterface.FileFinder(
        kappab='*', kappac='*', muR=1, muF=1, Q=1,
        directory = LatestPaths.derivedTheoryFiles_YukawaQuarkInducedScaled
        )
    containers_kappab_kappac_QuarkScaled = readList( theoryFiles_kappab_kappac_QuarkScaled )
    SM_kappab_kappac_QuarkScaled = [ container for container in containers_kappab_kappac_QuarkScaled if container.kappab == 1. and container.kappac == 1. ][0]
    containers_kappab_kappac_QuarkScaled.sort(key=kappaSorter)
    kappab_kappac_QuarkScaled = Container()
    kappab_kappac_QuarkScaled.name = 'kappab_kappac_QuarkScaled'
    kappab_kappac_QuarkScaled.couplings = [ 'kappab', 'kappac' ]
    kappab_kappac_QuarkScaled.containers = containers_kappab_kappac_QuarkScaled
    kappab_kappac_QuarkScaled.SM = SM_kappab_kappac_QuarkScaled
    kappab_kappac_QuarkScaled.DrawFileContentsAsLine = True
    return kappab_kappac_QuarkScaled

def read_top_highpt():
    containers_ct_cg_HighPt = TheoryFileInterface.FileFinder(
        ct='*', cg='*', cb=1, muR=1, muF=1, Q=1, filter='ct_1_cg_0', loadImmediately=True,
        directory = LatestPaths.derivedTheoryFiles_TopHighPt
        )
    SM_ct_cg_HighPt = TheoryFileInterface.FileFinder(
        ct=1, cg=0, cb=1, muR=1, muF=1, Q=1, expectOneFile=True, loadImmediately=True,
        directory = LatestPaths.derivedTheoryFiles_TopHighPt
        )
    containers_ct_cg_HighPt.sort(key=kappaSorter)
    ct_cg_HighPt = Container()
    ct_cg_HighPt.name = 'kappat_kappag_HighPt'
    ct_cg_HighPt.couplings = [ 'ct', 'cg' ]
    ct_cg_HighPt.containers = containers_ct_cg_HighPt
    ct_cg_HighPt.ws = LatestPaths.ws_combined_TopHighPt
    # ct_cg_HighPt.parametrizeByFitting = True
    # ct_cg_HighPt.includeLinearTerms = False
    ct_cg_HighPt.SM = SM_ct_cg_HighPt
    return ct_cg_HighPt

def read_top_highpt_ctcb():
    containers_ct_cb_HighPt = TheoryFileInterface.FileFinder(
        ct='*', cb='*', cg=0, muR=1, muF=1, Q=1, filter='ct_1_cg_0', loadImmediately=True,
        directory = LatestPaths.derivedTheoryFiles_TopHighPt
        )            
    SM_ct_cb_HighPt = TheoryFileInterface.FileFinder(
        ct=1, cg=0, cb=1, muR=1, muF=1, Q=1, expectOneFile=True, loadImmediately=True,
        directory = LatestPaths.derivedTheoryFiles_TopHighPt
        )
    containers_ct_cb_HighPt.sort(key=kappaSorter)
    ct_cb_HighPt = Container()
    ct_cb_HighPt.name = 'kappat_kappab_HighPt'
    ct_cb_HighPt.couplings = [ 'ct', 'cb' ]
    ct_cb_HighPt.containers = containers_ct_cb_HighPt
    ct_cb_HighPt.ws = LatestPaths.ws_combined_TopCtCbHighPt
    ct_cb_HighPt.SM = SM_ct_cb_HighPt
    return ct_cb_HighPt

def kappaSorter(container):
    if hasattr( container, 'kappab' ) and hasattr( container, 'kappac' ):
        return ( container.kappab, container.kappac )
    elif hasattr( container, 'ct' ) and hasattr( container, 'cg' ):
        return ( container.ct, container.cg )
    elif hasattr( container, 'ct' ) and hasattr( container, 'cb' ):
        return ( container.ct, container.cb )
    else:
        Commands.ThrowError( 'Can\'t sort; not a known set of couplings found' )


def theory_plot_container(couplingVariationContainer):

    DO_PARAMETRIZATION = True
    TOP_PANEL_LOGSCALE = True

    SKIP_KAPPAB_0_KAPPAC_1 = True
    # SKIP_KAPPAB_0_KAPPAC_1 = False

    PLOT_CROSSSECTIONS = True
    # PLOT_CROSSSECTIONS = False

    newColorCycle = lambda : itertools.cycle( range(2,10) + range(27,50) )


    # ======================================
    # Loop

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


def SetParametrizedTgs(couplingVariationContainer):

    if hasattr( couplingVariationContainer, 'ws' ):
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


    elif hasattr( couplingVariationContainer, 'parametrizeByFitting' ):
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


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

########################################
# End of Main
########################################
if __name__ == "__main__":
    print 'Don\'t run this program directly'