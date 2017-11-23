#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import Commands

import os, itertools, operator, re, argparse, sys, random
from math import isnan, isinf
from os.path import *
from glob import glob
from copy import deepcopy
from array import array


from time import strftime
datestr = strftime( '%b%d' )

import ROOT


import TheoryCommands
from TheoryCommands import c
from TheoryCommands import SetCMargins
from TheoryCommands import SetPlotDir
from TheoryCommands import SaveAsRoot
from TheoryCommands import SaveC
from TheoryCommands import GetUniqueRootName
from TheoryCommands import GetPlotBase

from Container import Container

import PhysicsCommands
import CorrelationMatrices
from Parametrization import Parametrization, WSParametrization

from time import strftime
datestr = strftime( '%b%d' )


########################################
# Functions
########################################


#____________________________________________________________________
def PlotParametrizationsOnCombination(
        container,
        ):

    expBinBoundaries    = container.expBinBoundaries
    ws_combination      = container.ws_combination
    scanDir_combination = container.scanDir_combination
    ws_coupling         = container.ws_coupling
    scanDir_coupling    = container.scanDir_coupling
    xCoupling           = container.xCoupling
    yCoupling           = container.yCoupling
    xCouplingTitle      = container.xCouplingTitle
    yCouplingTitle      = container.yCouplingTitle
    plotTitle           = container.plotTitle


    DO_X_MAXIMA = False
    if hasattr( container, 'OnlyXMaximaOfContour' ) and container.OnlyXMaximaOfContour:
        DO_X_MAXIMA = True

    DO_Y_MAXIMA = False
    if hasattr( container, 'OnlyYMaximaOfContour' ) and container.OnlyYMaximaOfContour:
        DO_Y_MAXIMA = True

    if hasattr( container, 'MaximaOfContour' ) and container.MaximaOfContour:
        DO_X_MAXIMA = True
        DO_Y_MAXIMA = True


    DO_BESTFIT = False
    if hasattr( container, 'BestFit' ) and container.BestFit:
        DO_BESTFIT = True

    STRAIGHT_LINE_TO_SM = False
    if hasattr( container, 'StraightLineToSM' ) and container.StraightLineToSM:
        STRAIGHT_LINE_TO_SM = True
        if not hasattr( container, 'xSM' ):
            Commands.ThrowError( 'Requested mode \'StraightLineToSM\', but xSM is not specified' )
        if not hasattr( container, 'ySM' ):
            Commands.ThrowError( 'Requested mode \'StraightLineToSM\', but ySM is not specified' )


    # ======================================
    # Begin plotting

    newColorCycle = lambda: itertools.cycle( [ 2, 4, 6, 41, 46, 30, 43, 3, 5, 8, 9 ] )

    c.Clear()
    SetCMargins()

    xMin = expBinBoundaries[0]
    xMax = expBinBoundaries[-2] + ( expBinBoundaries[-2] - expBinBoundaries[-3] ) # Overflow will screw up plot
    yMin = 0.0
    yMax = 3.0

    base = GetPlotBase(
        xMin = xMin,
        xMax = xMax,
        yMin = yMin,
        yMax = yMax,
        xTitle = 'p_{T} [GeV]', yTitle = '#mu'
        )
    base.Draw('P')

    leg = ROOT.TLegend(
        c.GetLeftMargin(),
        1 - c.GetTopMargin() - 0.2,
        1 - c.GetRightMargin(),
        1 - c.GetTopMargin() 
        )
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetNColumns(2)


    # ======================================
    # Draw the combination as blocks

    combinationPOIs = Commands.ListPOIs( ws_combination )
    combinationscans = PhysicsCommands.GetScanResults( combinationPOIs, scanDir_combination, pattern = '' )

    Tg = PhysicsCommands.GetTGraphForSpectrum(
        combinationPOIs,
        combinationscans,
        name = GetUniqueRootName()
        )

    Tg.SetLineColor( 1 )
    # Tg.SetMarkerStyle( 2 )
    Tg.SetFillColorAlpha( 1, 0.2 )
    # Tg.SetFillColor( 13 )
    # Tg.SetFillStyle( 3544 )
    # Tg.SetFillStyle( 3345 )

    Tg.title = 'Combination'

    CorrelationMatrices.ConvertTGraphToLinesAndBoxes(
        Tg,
        drawImmediately=True,
        legendObject=leg,
        verbose=False,
        noBoxes=False,
        # xMinExternal=xMin,
        xMaxExternal=xMax,
        yMinExternal=yMin,
        yMaxExternal=yMax,
        )


    # ======================================
    # Get a contour

    combined_rootfiles  = glob( '{0}/*.root'.format( scanDir_coupling ) )
    combined = TheoryCommands.GetTH2FromListOfRootFiles(
        combined_rootfiles, xCoupling, yCoupling, verbose = False,
        )
    combined.color = 1
    combined.name = 'regular'
    combined.title = 'Nominal'

    allcontours = TheoryCommands.GetContoursFromTH2( combined.H2, 2.30 )
    candidatecontours = []

    xBestfit = combined.xBestfit
    yBestfit = combined.yBestfit
    for Tg in allcontours:
        Tg.x, Tg.y = TheoryCommands.GetXYfromTGraph(Tg)
        Tg.xMin = min(Tg.x)
        Tg.xMax = max(Tg.x)
        Tg.yMin = min(Tg.y)
        Tg.yMax = max(Tg.y)

        # Check if bestfit is at least inside minima and maxima
        if not (
                Tg.xMin < xBestfit
                and Tg.xMax > xBestfit
                and Tg.yMin < yBestfit
                and Tg.yMax > yBestfit
                ):
            continue

        # Compute some numbers that may help in selection: minimum distance to bestfit
        Tg.minDist = min([ ( Tg.x[i] - xBestfit )**2 + ( Tg.y[i] - yBestfit )**2 for i in xrange(Tg.GetN()) ])

        # Distance to bestfit should be close to half the difference between xMax-xMin (or yMax-yMin)
        # Compute ratio of this, minus 1 - result *should* be close to zero
        Tg.distRatio = Tg.minDist / ( 0.5 * min( Tg.yMax-Tg.yMin, Tg.xMax-Tg.xMin ) )  -  1.0
        # if abs(Tg.distRatio) > 0.8: continue

        candidatecontours.append( Tg )


    if len(candidatecontours) == 0:
        Commands.ThrowError( 'Can\'t find contour' )
    elif len(candidatecontours) > 1:
        candidatecontours.sort( key = lambda Tg: Tg.minDist )
        candidatecontours = candidatecontours[:2]

    # # Pick the contour with the highest ratio (more likely to be the 'outer shell')
    # candidatecontours.sort( key = lambda Tg: Tg.distRatio, reverse=True )

    # Actually, pick the 'inner' shell, outer shell is too likely to be a misfit
    candidatecontours.sort( key = lambda Tg: Tg.distRatio )    

    contour = candidatecontours[0]

    # ======================================
    # Determine the points to analyze:

    points = []

    if DO_BESTFIT:
        points.append( ( xBestfit, yBestfit ) )

    if DO_X_MAXIMA:
        points.extend([
            ( contour.xMin, contour.y[ contour.x.index(contour.xMin) ] ),
            ( contour.xMax, contour.y[ contour.x.index(contour.xMax) ] ),
            ])

    if DO_Y_MAXIMA:
        points.extend([
            ( contour.x[ contour.y.index(contour.yMin) ], contour.yMin ),
            ( contour.x[ contour.y.index(contour.yMax) ], contour.yMax ),
            ])

    if STRAIGHT_LINE_TO_SM:
        dxSM = xBestfit - container.xSM
        dySM = yBestfit - container.ySM

        points.extend([
            ( xBestfit - 2.0*dxSM, yBestfit - 2.0*dySM ),
            ( xBestfit - 1.5*dxSM, yBestfit - 1.5*dySM ),
            ( xBestfit - 0.5*dxSM, yBestfit - 0.5*dySM ),
            ( xBestfit + 0.5*dxSM, yBestfit + 0.5*dySM ),
            # ( xBestfit - 1.5*dxSM, yBestfit - 1.5*dySM ),
            ])


    if hasattr( container, 'ManualPoints' ):
        for xRaw, yRaw in container.ManualPoints:

            if callable(xRaw):
                try:
                    x = xRaw( xBestfit )
                except TypeError:
                    x = xRaw( xBestfit, container.xSM )
            else:
                x = xRaw

            if callable(yRaw):
                try:
                    y = yRaw( yBestfit )
                except TypeError:
                    y = yRaw( yBestfit, container.ySM )
            else:
                y = xRaw


            points.append( ( x, y ) )


    # ======================================
    # Get the parametrization

    wsParametrization = WSParametrization( ws_coupling )

    colorCycle = newColorCycle()
    for xPoint, yPoint in points:
        color = next(colorCycle)

        kwargs = { xCoupling : xPoint, yCoupling : yPoint }
        Tg_param = wsParametrization.GetOutputContainer( returnWhat='exp', **kwargs ).Tg

        Tg_param.SetLineColor(color)
        Tg_param.SetMarkerColor(color)
        Tg_param.SetLineStyle(1)
        Tg_param.SetLineWidth(4)

        Tg_param.title = '{0} = {1:.1f}, {2} = {3:.1f}'.format(
            xCouplingTitle, xPoint, yCouplingTitle, yPoint
            )

        CorrelationMatrices.ConvertTGraphToLinesAndBoxes(
            Tg_param,
            drawImmediately=True,
            legendObject=leg,
            noBoxes=True,
            xMaxExternal=xMax
            )

        Tg_param.RemovePoint( Tg_param.GetN()-1 )
        Tg_param.SetLineWidth(1)
        Tg_param.Draw('SAMEL')

    leg.Draw()

    SaveC( '{0}_onCombination'.format(plotTitle) )


    # ======================================
    # Second plot: The contour with the points on it

    c.Clear()
    SetCMargins( for2Dhist = True )

    combined.H2.GetXaxis().SetLabelSize(0.045)
    combined.H2.GetYaxis().SetLabelSize(0.045)

    combined.H2.GetXaxis().SetTitleSize(0.055)
    combined.H2.GetYaxis().SetTitleSize(0.055)

    # combined.H2.GetXaxis().SetTitleOffset(1.2)
    # combined.H2.GetYaxis().SetTitleOffset(1.2)

    combined.H2.GetXaxis().SetTitle( xCouplingTitle )
    combined.H2.GetYaxis().SetTitle( yCouplingTitle )

    combined.H2.Draw('COLZ')
    combined.H2.SetMaximum(7.0)

    contour.SetLineColor(1)
    contour.SetLineWidth(3)
    contour.Draw('SAMEL')


    colorCycle = newColorCycle()
    for xPoint, yPoint in points:
        color = next(colorCycle)

        point = ROOT.TGraph( 1, array( 'f', [xPoint] ), array( 'f', [yPoint] ) )
        ROOT.SetOwnership( point, False )

        point.SetMarkerStyle(33)
        point.SetMarkerColor(color)
        point.SetMarkerSize( 4.5 if not hasattr( container, 'MarkerSize' ) else container.MarkerSize )
        point.Draw('SAMEP')

        point2 = ROOT.TGraph( point )
        ROOT.SetOwnership( point2, False )
        point2.SetMarkerStyle( 27 )
        point2.SetMarkerColor(1)
        point2.Draw('SAMEP')

    c.Update()
    SaveC( '{0}_pointsOnContour'.format(plotTitle) )






#____________________________________________________________________
def PlotMultipleScans(
        containers,
        xMin      = None,
        xMax      = None,
        yMin      = None,
        yMax      = None,
        xTitle    = 'x',
        yTitle    = 'y',
        plotname  = 'unnamedscans',
        draw1sigmaline = True,
        translateToChi2 = False,
        printUncertainties = False,
        ):

    if xMin is None: xMin = min([ min(container.x) for container in containers if not hasattr( container, 'line' ) ])
    if xMax is None: xMax = max([ max(container.x) for container in containers if not hasattr( container, 'line' ) ])
    if yMin is None: yMin = min([ min(container.y) for container in containers if not hasattr( container, 'line' ) ])
    if yMax is None: yMax = max([ max(container.y) for container in containers if not hasattr( container, 'line' ) ])


    # ======================================
    # Make plot

    c.cd()
    c.Clear()
    SetCMargins()


    base = GetPlotBase(
        xMin = xMin,
        xMax = xMax,
        yMin = yMin,
        yMax = yMax,
        xTitle = xTitle,
        yTitle = yTitle,
        )
    base.Draw('P')

    base.GetXaxis().SetTitleSize(0.06)
    base.GetXaxis().SetLabelSize(0.05)
    base.GetYaxis().SetTitleSize(0.06)
    base.GetYaxis().SetLabelSize(0.05)


    leg = ROOT.TLegend(
        c.GetLeftMargin() + 0.01,
        c.GetBottomMargin() + 0.02,
        1 - c.GetRightMargin() - 0.01,
        c.GetBottomMargin() + 0.09
        )
    leg.SetNColumns(3)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    colorCycle = itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )
    for container in containers:

        if hasattr( container, 'line' ):
            if container.line.GetX1() == -999: container.line.SetX1( xMin )
            if container.line.GetX2() ==  999: container.line.SetX2( xMax )
            if container.line.GetY1() == -999: container.line.SetY1( yMin )
            if container.line.GetY2() ==  999: container.line.SetY2( yMax )
            container.line.Draw()
            continue

        if printUncertainties:
            container.unc = PhysicsCommands.FindMinimaAndErrors( container.x, container.y, returnContainer=True )

        if translateToChi2:
            container.y = [ 2.*y for y in container.y ]

        if not hasattr( container, 'Tg' ):
            container.Tg = ROOT.TGraph(
                len(container.x), array( 'f', container.x ), array( 'f', container.y )
                )


        if not hasattr( container, 'color' ):
            container.color = next(colorCycle)

        container.Tg.SetLineColor(container.color)
        container.Tg.SetMarkerColor(container.color)
        container.Tg.SetLineWidth(2)
        container.Tg.SetMarkerSize(0.9)
        container.Tg.SetMarkerStyle(5)

        container.Tg.Draw('PLSAME')


        if printUncertainties:

            withPercentages = False

            bestfitLine = ROOT.TLine( container.unc.min, yMin, container.unc.min, yMax )
            bestfitLine.SetLineWidth(2)
            bestfitLine.SetLineColor(container.color)
            bestfitLine.Draw()


            l = ROOT.TLatex()
            l.SetTextColor(container.color)
            l.SetTextSize(0.04)

            if withPercentages:
                l.SetTextAlign(31)
                l.DrawLatex(
                    container.unc.leftBound, 1.0 + 0.013*(yMax-yMin),
                    '-{0:.2f} ({1:d}%)'.format(
                        abs(container.unc.leftError),
                        int( abs(container.unc.leftError) / container.unc.min * 100. )
                        ))

                l.SetTextAlign(11)
                l.DrawLatex(
                    container.unc.rightBound, 1.0 + 0.013*(yMax-yMin),
                    '+{0:.2f} ({1:d}%)'.format(
                        abs(container.unc.rightError),
                        int( abs(container.unc.rightError) / container.unc.min * 100. )
                        ))
            else:
                l.SetTextAlign(31)
                l.DrawLatex(
                    container.unc.leftBound - 0.013*(xMax-xMin), 1.0 + 0.013*(yMax-yMin),
                    '{0:+.2f}'.format(
                        abs(container.unc.leftBound),
                        ))
                l.SetTextAlign(11)
                l.DrawLatex(
                    container.unc.rightBound + 0.013*(xMax-xMin), 1.0 + 0.013*(yMax-yMin),
                    '{0:+.2f}'.format(
                        abs(container.unc.rightBound),
                        ))


            l.SetTextAlign(21)
            l.DrawLatex(
                container.unc.min, 1.0 + 0.013*(yMax-yMin),
                '{0:.3f}'.format( container.unc.min )
                )


            TgPoints = ROOT.TGraph( 2,
                array( 'f', [ container.unc.leftBound, container.unc.rightBound ] ),
                array( 'f', [ 1.0, 1.0 ] ),
                )
            TgPoints.SetMarkerSize(1.2)
            TgPoints.SetMarkerStyle(8)
            TgPoints.SetMarkerColor(container.color)
            TgPoints.Draw('PSAME')


    if draw1sigmaline:
        y1sigma = 0.5
        if translateToChi2: y1sigma = 1.0
        line1sigma = ROOT.TLine( xMin, y1sigma, xMax, y1sigma )
        line1sigma.SetLineColor(14)
        line1sigma.Draw()

    SaveC( plotname )



#____________________________________________________________________
def FilterContourHeuristic(
        contours,
        xBestfit,
        yBestfit,
        ):

    filteredContours = []

    for Tg in contours:
        Tg.x, Tg.y = TheoryCommands.GetXYfromTGraph(Tg)
        Tg.xMin = min(Tg.x)
        Tg.xMax = max(Tg.x)
        Tg.yMin = min(Tg.y)
        Tg.yMax = max(Tg.y)

        # Check if bestfit is at least inside minima and maxima
        if not (
                Tg.xMin < xBestfit
                and Tg.xMax > xBestfit
                and Tg.yMin < yBestfit
                and Tg.yMax > yBestfit
                ):
            continue

        # Compute some numbers that may help in selection: minimum distance to bestfit
        Tg.minDist = min([ ( Tg.x[i] - xBestfit )**2 + ( Tg.y[i] - yBestfit )**2 for i in xrange(Tg.GetN()) ])

        # Distance to bestfit should be close to half the difference between xMax-xMin (or yMax-yMin)
        # Compute ratio of this, minus 1 - result *should* be close to zero
        Tg.distRatio = Tg.minDist / ( 0.5 * min( Tg.yMax-Tg.yMin, Tg.xMax-Tg.xMin ) )  -  1.0
        # if abs(Tg.distRatio) > 0.8: continue

        filteredContours.append( Tg )


    if len(filteredContours) == 0:
        # Heuristic failed
        filteredContours = contours
        return filteredContours
    elif len(filteredContours) > 1:
        filteredContours.sort( key = lambda Tg: Tg.minDist )
        filteredContours = filteredContours[:2]

    # Pick the contour with the highest ratio (more likely to be the 'outer shell')
    filteredContours.sort( key = lambda Tg: Tg.distRatio, reverse=True )

    return filteredContours


#____________________________________________________________________
def BasicMixedContourPlot(
        containers,
        xMin      = 0.,
        xMax      = 1.,
        yMin      = 0.,
        yMax      = 1.,
        xTitle    = 'x',
        yTitle    = 'y',
        plotname  = 'contours',
        x_SM      = 1.,
        y_SM      = 1.,
        plotIndividualH2s = False,
        filterContours = True,
        ):

    print '\nRunning BasicMixedContourPlot for \'{0}\''.format( plotname )

    # ======================================
    # Check whether the passed containers fulfill requirements

    for container in containers:
        attrs = container.ListAttributes()

        for expectedAttr in [ 'H2', 'name' ]:
            if not expectedAttr in attrs:
                Commands.ThrowError(
                    'Container misses mandatory attribute \'{0}\' (defined attributes: {1})'.format( expectedAttr, ', '.join(attrs) ),
                    throwException = True
                    )

        if not hasattr( container, 'color' ):
            container.color = 1


    # ======================================
    # Calculate contours

    for container in containers:
        print '\nGetting contours for {0}'.format( container.name )
        container.contours_1sigma = TheoryCommands.GetContoursFromTH2( container.H2, 2.30 )
        container.contours_2sigma = TheoryCommands.GetContoursFromTH2( container.H2, 6.18 )

        if filterContours:
            container.contours_1sigma = FilterContourHeuristic( container.contours_1sigma, container.xBestfit, container.yBestfit )
            container.contours_2sigma = FilterContourHeuristic( container.contours_2sigma, container.xBestfit, container.yBestfit )

        print '  Bestfit:'
        print '  {0} = {1}'.format( xTitle, container.xBestfit )
        print '  {0} = {1}'.format( yTitle, container.yBestfit )


    # ======================================
    # Make plot

    c.cd()
    c.Clear()
    SetCMargins()


    base = GetPlotBase(
        xMin = xMin,
        xMax = xMax,
        yMin = yMin,
        yMax = yMax,
        xTitle = xTitle,
        yTitle = yTitle,
        )
    base.Draw('P')

    base.GetXaxis().SetTitleSize(0.06)
    base.GetXaxis().SetLabelSize(0.05)
    base.GetYaxis().SetTitleSize(0.06)
    base.GetYaxis().SetLabelSize(0.05)


    leg = ROOT.TLegend(
        c.GetLeftMargin() + 0.01,
        c.GetBottomMargin() + 0.02,
        1 - c.GetRightMargin() - 0.01,
        c.GetBottomMargin() + 0.09
        )
    leg.SetNColumns( min( 3, len(containers) ) )
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)


    for container in containers:

        for Tg in container.contours_1sigma:

            Tg.SetLineWidth(2)
            Tg.SetLineColor( container.color )
            Tg.SetLineStyle(1)
            Tg.Draw('CSAME')
            if Tg == container.contours_1sigma[0]:
                Tg.SetName( '{0}_contour_1sigma'.format(container.name) )
                leg.AddEntry(
                    Tg.GetName(),
                    container.name if not hasattr( container, 'title' ) else container.title,
                    'l' )

        for Tg in container.contours_2sigma:
            Tg.SetLineWidth(2)
            Tg.SetLineColor( container.color )
            Tg.SetLineStyle(2)
            Tg.Draw('CSAME')

        Tpoint = ROOT.TGraph( 1, array( 'd', [container.xBestfit] ), array( 'd', [container.yBestfit] ) )
        ROOT.SetOwnership( Tpoint, False )
        Tpoint.SetMarkerSize(2)
        Tpoint.SetMarkerStyle(34)
        Tpoint.SetMarkerColor( container.color )
        Tpoint.Draw('PSAME')
        Tpoint.SetName( '{0}_bestfitpoint'.format( container.name ) )
        container.bestfitPoint = Tpoint

    TpointSM = ROOT.TGraph( 1, array( 'd', [x_SM] ), array( 'd', [y_SM] ) )
    ROOT.SetOwnership( TpointSM, False )
    TpointSM.SetMarkerSize(2)
    TpointSM.SetMarkerStyle(21)
    TpointSM.SetMarkerColor( 12 )
    TpointSM.Draw('PSAME')

    leg.Draw()

    SaveC( plotname )


    if plotIndividualH2s:

        for container in containers:

            c.Clear()
            SetCMargins(
                LeftMargin   = 0.12,
                RightMargin  = 0.10,
                BottomMargin = 0.12,
                TopMargin    = 0.09,
                )

            container.H2.Draw('COLZ')

            container.H2.GetXaxis().SetRangeUser( xMin, xMax )
            container.H2.GetYaxis().SetRangeUser( yMin, yMax )
            container.H2.SetMaximum( 7. )

            container.H2.GetXaxis().SetTitle( xTitle )
            container.H2.GetYaxis().SetTitle( yTitle )
            container.H2.GetXaxis().SetTitleSize(0.06)
            container.H2.GetXaxis().SetLabelSize(0.05)
            container.H2.GetYaxis().SetTitleSize(0.06)
            container.H2.GetYaxis().SetLabelSize(0.05)

            for Tg in container.contours_1sigma: Tg.Draw('CSAME')
            for Tg in container.contours_2sigma: Tg.Draw('CSAME')
            container.bestfitPoint.Draw('PSAME')

            c.Update()

            SaveC( plotname + '_' + container.name )



########################################
# Execute this if called directly
########################################

if __name__ == '__main__':
    print 'Import functions from this file'
    print 'Available functions:'
    functions = [ name for (name, thing) in locals().items() if callable(thing) ]
    print '\n'.join( functions )