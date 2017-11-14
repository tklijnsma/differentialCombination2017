#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, re, itertools, sys
from os.path import *
from glob import glob
from sys import exit
from copy import deepcopy
from numpy import corrcoef, var, std

from time import strftime
datestr = strftime( '%b%d' )

datestr_detailed = strftime( '%y-%m-%d %H:%M:%S' )
src = os.path.basename(__file__)
from subprocess import check_output
currentcommit = check_output(['git', 'log', '-1', '--oneline' ])


import Commands
import TheoryCommands
from Container import Container

import ROOT
from array import array

ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
ROOT.gStyle.SetOptStat(0)
c = ROOT.TCanvas( 'c', 'c', 1000, 800 )

CLeftMargin   = 0.15
CRightMargin  = 0.03
CBottomMargin = 0.15
CTopMargin    = 0.03
def SetCMargins(
    LeftMargin   = CLeftMargin,
    RightMargin  = CRightMargin,
    BottomMargin = CBottomMargin,
    TopMargin    = CTopMargin,
    ):
    c.SetLeftMargin( LeftMargin )
    c.SetRightMargin( RightMargin )
    c.SetBottomMargin( BottomMargin )
    c.SetTopMargin( TopMargin )

def GetPlotBase(
    xMin = 0, xMax = 1,
    yMin = 0, yMax = 1,
    xTitle = 'x', yTitle = 'y',
    SetTitleSizes = True,
    ):

    base = ROOT.TH1F()
    ROOT.SetOwnership( base, False )
    base.SetName( 'genericbase' )
    base.GetXaxis().SetLimits( xMin, xMax )
    base.SetMinimum( yMin )
    base.SetMaximum( yMax )
    base.SetMarkerColor(0)
    base.GetXaxis().SetTitle( xTitle )
    base.GetYaxis().SetTitle( yTitle )

    if SetTitleSizes:
        base.GetXaxis().SetTitleSize( 0.06 )
        base.GetYaxis().SetTitleSize( 0.06 )

    return base

CPlotDir = 'correlationMatrices_{0}'.format(datestr)
def SetPlotDir( newdir ):
    global CPlotDir
    CPlotDir = newdir
def SaveC( outname, subdir=None, asPDF=True, asPNG=False, asROOT=False ):

    if subdir == None:
        outdir = CPlotDir
    else:
        outdir = join( CPlotDir, subdir )

    if not isdir(outdir): os.makedirs( outdir )
    outname = join( outdir, outname.replace('.pdf','').replace('.png','').replace('.root','') )

    if asPDF: c.SaveAs( outname + '.pdf' )
    if asPNG: c.SaveAs( outname + '.png' )
    if asROOT or TheoryCommands.SAVEROOT: c.SaveAs( outname + '.root' )


########################################
# Main
########################################

def main():

    # variationFiles = glob( '../suppliedInput/fromAgnieszka/ScaleVarNNLO_Jul17/*.top' )

    # # Add central SM point as well - Actually, not necessary, HRes_mR1mF1.top is the central SM
    # # variationFiles.append( '../suppliedInput/fromAgnieszka/HRes_SMEFT_May16/SM_NNLO' )

    # # ======================================
    # # In theory binning

    # variations = []
    # for variationFile in variationFiles:
    #     variation = ReadVariationFile( variationFile )
    #     variations.append( variation )

    # GetCorrelationMatrix(
    #     variations,
    #     makeScatterPlots          = True,
    #     makeCorrelationMatrixPlot = True,
    #     outname                   = 'corrMat_theory',
    #     verbose                   = True,
    #     )

    print 'Don\'t execute this directly anymore.'
    return



def mean( l ):
    return sum(l) / len(l)

def GetCorrelationMatrix(
        variations,
        makeScatterPlots = False,
        makeCorrelationMatrixPlot = True,
        outname = 'corrMat',
        verbose = False,
        halfNumberOfPlots = False,
        withRespectToCentralScale = False,
        ):

    if isinstance( variations[0], Container ):
        # Probably a returned object by TheoryCommands.ReadDerivedTheoryFile;
        # make sure all the proper attributes are filled
        for variation in variations:
            variation.binValues  = variation.crosssection
            variation.binCenters = [
                0.5*(variation.binBoundaries[i]+variation.binBoundaries[i+1]) for i in xrange(len(variation.binValues))
                ]

    INCLUDE_LABELS = True

    if INCLUDE_LABELS:
        labels = []
        for variation in variations:
            labels.append(
                '(#mu_{{R}}={0}, #mu_{{F}}={1})'.format( variation.muR, variation.muF )
                )


    if not makeScatterPlots:
        halfNumberOfPlots = True

    # ======================================
    # Try to get correlations

    nPoints = len( variations[0].binValues )
    binBoundaries = variations[0].binBoundaries

    if makeScatterPlots:
        base = GetPlotBase()

    corrMatrix = [ [ 999 for jPoint in xrange(nPoints) ] for iPoint in xrange(nPoints) ]
    errors     = []

    jStartPoint = 0
    for iPoint in xrange(nPoints):
        if verbose: print 'Processing point {0}/{1}'.format( iPoint, nPoints-1 )
        iBinValues = [ variation.binValues[iPoint] for variation in variations ]
        iPt = variations[0].binCenters[iPoint]

        binWidth   = binBoundaries[iPoint+1] - binBoundaries[iPoint]
        maxVarUp   = max( iBinValues )
        maxVarDown = min( iBinValues )
        if withRespectToCentralScale:
            print 'Not yet implemented'
        else:
            errors.append((
                binWidth * abs( maxVarUp - mean(iBinValues) ),
                binWidth * abs( maxVarDown - mean(iBinValues) ),
                ))

        # Essentially doing this loop from 0 to nPoints is double work, but for the plots
        # it's nicer to do it anyway
        if halfNumberOfPlots: jStartPoint = iPoint
        for jPoint in xrange(jStartPoint,nPoints):
            jBinValues = [ variation.binValues[jPoint] for variation in variations ]
            jPt = variations[0].binCenters[jPoint]

            corr = corrcoef( iBinValues, jBinValues )[0][1]

            try:
                corrMatrix[iPoint][jPoint] = corr
            except IndexError:
                print 'Problem with the dimensions of corrMatrix'
                print '  nRows = {0} ;  nCols = {1}'.format( len(corrMatrix), len(corrMatrix[0]) )
                print '  Trying to access point ( {0} , {1} )'.format( iPoint, jPoint )
                sys.exit()
            
            if halfNumberOfPlots:
                corrMatrix[jPoint][iPoint] = corr

            if makeScatterPlots:
                ########################################
                # Make plot
                ########################################

                c.Clear()
                SetCMargins()

                xMin = min(iBinValues) - 0.3*( max(iBinValues) - min(iBinValues) )
                xMax = max(iBinValues) + 0.3*( max(iBinValues) - min(iBinValues) )
                yMin = min(jBinValues) - 0.3*( max(jBinValues) - min(jBinValues) )
                yMax = max(jBinValues) + 0.3*( max(jBinValues) - min(jBinValues) )

                base.GetXaxis().SetLimits( xMin, xMax )
                base.SetMinimum( yMin )
                base.SetMaximum( yMax )

                base.GetXaxis().SetTitle( 'p_{{T}} = {0:.1f} (Bin {1})'.format( iPt, iPoint ) )
                base.GetYaxis().SetTitle( 'p_{{T}} = {0:.1f} (Bin {1})'.format( jPt, jPoint ) )

                base.Draw('P')


                Tg = ROOT.TGraph( len( iBinValues ), array( 'd', iBinValues ), array( 'd', jBinValues ) )
                ROOT.SetOwnership( Tg, False )
                Tg.SetMarkerColor(9)
                Tg.SetMarkerStyle(8)
                Tg.SetMarkerSize(0.8)
                Tg.Draw('PSAME')

                if INCLUDE_LABELS:
                    l = ROOT.TLatex()
                    l.SetTextSize( 0.03 )
                    # l.SetTextAlign(21)
                    for iVariation in xrange( len( iBinValues ) ):
                        l.DrawLatex( iBinValues[iVariation], jBinValues[iVariation], labels[iVariation] )


                # ======================================
                # Draw helpful lines

                iMean = sum(iBinValues)/len(iBinValues)
                iMeanLine = ROOT.TLine( iMean, yMin, iMean, yMax )
                ROOT.SetOwnership( iMeanLine, False )
                iMeanLine.Draw()

                jMean = sum(jBinValues)/len(jBinValues)
                jMeanLine = ROOT.TLine( xMin, jMean, xMax, jMean )
                ROOT.SetOwnership( jMeanLine, False )
                jMeanLine.Draw()

                slope  = corr * std(jBinValues)/std(iBinValues)
                offset = jMean - slope * iMean
                fnCorrLine = lambda x: slope*x + offset
                corrLine = ROOT.TLine( xMin, fnCorrLine(xMin), xMax, fnCorrLine(xMax) )
                ROOT.SetOwnership( corrLine, False )
                corrLine.SetLineWidth(2)
                corrLine.SetLineColor(2)
                corrLine.Draw()

                # Labels
                lbl = ROOT.TLatex()
                lbl.SetNDC()
                lbl.SetTextAlign(33)
                lbl.SetTextSize(0.06)
                lbl.DrawLatex( 1-CRightMargin-0.01, 1-CTopMargin-0.01, '#rho = {0:.2f}'.format(corr) )

                SaveC( 'Correlation_Bin{0}_Bin{1}'.format( iPoint, jPoint ), subdir = '{0}/Bin{1}'.format( outname, iPoint ), asPNG=True )
                if halfNumberOfPlots:
                    SaveC( 'Correlation_Bin{0}_Bin{1}'.format( jPoint, iPoint ), subdir = '{0}/Bin{1}'.format( outname, jPoint ), asPNG=True )

                del Tg



    if not isdir(CPlotDir): os.makedirs( CPlotDir )

    with open( join( CPlotDir, '{0}.txt'.format(outname) ), 'w' ) as corrMatFp:
        for iPoint in xrange(nPoints):
            corrMatFp.write(
                ' '.join([ '{0:+.4f}'.format(corrMatrix[iPoint][jPoint]) for jPoint in xrange(nPoints) ])
                + '\n'
                )

    with open( join( CPlotDir, 'errors_for_{0}.txt'.format(outname) ), 'w' ) as errorsFp:
        for up, down in errors:
            errorsFp.write(
                '{0:+.8f} {1:+.8f}\n'.format( up, down )
                )

    with open( join( CPlotDir, 'errors_for_{0}.tex'.format(outname) ), 'w' ) as errorsFp:

        header    = [ 'Bin boundaries (GeV)' ]
        scaleDown = [ '$\\Delta^\\text{scale}_-$ (pb)' ]
        scaleUp   = [ '$\\Delta^\\text{scale}_+$ (pb)' ]

        for iBin in xrange(len(binBoundaries)-1):
            up, down = errors[iBin]
            down = -abs(down)

            binBoundLeft  = int(binBoundaries[iBin])
            binBoundRight = int(binBoundaries[iBin+1])

            header.append( '{0} to {1}'.format( binBoundLeft, binBoundRight ) )
            scaleDown.append( '{0:+.2f}'.format(down) )
            scaleUp.append(   '{0:+.2f}'.format(up) )

        errorsFp.write( '% File generated on {0} by the script {1}\n'.format( datestr_detailed, src ) )
        errorsFp.write( '% Current git commit: {0}'.format( currentcommit ) )
        errorsFp.write( ' & '.join(header) + ' \\\\\n' )
        errorsFp.write( ' & '.join(scaleDown) + ' \\\\\n' )
        errorsFp.write( ' & '.join(scaleUp) + ' \\\\' )


    if makeCorrelationMatrixPlot:

        # ======================================
        # Make plot of correlation matrix

        c.Clear()
        SetCMargins(
            LeftMargin   = 0.21,
            RightMargin  = 0.12,
            TopMargin    = 0.12,
            BottomMargin = 0.19,
            )

        T = ROOT.TH2D(
            'corrMat', '#scale[0.85]{Correlation between p_{T} bins}',
            nPoints, 0., nPoints,
            nPoints, 0., nPoints
            )
        ROOT.SetOwnership( T, False )
        T.SetContour(100)

        for iRow in xrange(nPoints):
            for iCol in xrange(nPoints):
                T.SetBinContent( iCol+1, iRow+1, corrMatrix[iRow][iCol] )

        T.GetZaxis().SetRangeUser(-1.0,1.0)


        toString = lambda number: str(int(number)) if number.is_integer() else '{0:.1f}'.format(number)
        for iPoint in xrange(nPoints):
            if nPoints < 20 or iPoint % int(0.1*nPoints) == 0:
                if not hasattr( variations[0], 'binBoundaries' ):
                    binLabel = toString(variations[0].binCenters[iPoint])
                else:
                    binLabel = '{0} - {1}'.format(
                        toString(variations[0].binBoundaries[iPoint] ),
                        toString(variations[0].binBoundaries[iPoint+1] )
                        )
                T.GetXaxis().SetBinLabel( iPoint+1, binLabel )
                T.GetYaxis().SetBinLabel( iPoint+1, binLabel )


        n_stops = 3
        stops  = [ 0.0, 0.5, 1.0 ]
        reds   = [ 0.0, 1.0, 1.0 ]
        blues  = [ 1.0, 1.0, 0.0 ]
        greens = [ 0.0, 1.0, 0.0 ]

        ROOT.TColor.CreateGradientColorTable(
            n_stops,
            array('d', stops ),
            array('d', reds ),
            array('d', greens ),
            array('d', blues ),
            255 )

        T.GetXaxis().SetTitle( 'p_{T}' )
        T.GetXaxis().SetTitleOffset( 1.7 )
        T.GetXaxis().SetTitleSize(0.045)

        T.GetYaxis().SetTitle( 'p_{T}' )
        T.GetYaxis().SetTitleOffset( 2.45 )
        T.GetYaxis().SetTitleSize(0.045)

        T.GetXaxis().SetLabelSize(0.045)
        T.GetYaxis().SetLabelSize(0.045)

        ROOT.gStyle.SetHistMinimumZero() # To draw the "0", otherwise ROOT leaves it empty
        ROOT.gStyle.SetPaintTextFormat('1.2g')

        T.SetMarkerSize(1.35)
        
        if nPoints < 20:
            T.Draw('COLZ TEXT')
        else:
            T.Draw('COLZ')


        c.cd()
        c.Update()

        SaveC( outname, asPNG=True )

    # Return the found correlation matrix
    return corrMatrix



def ReadVariationFile(
        variationFile,
        fromAgnieszka=False,
        fromHqT=False,
        ):

    if not isfile( variationFile ):
        print 'File \'{0}\' does not exist'.format( variationFile )
        return

    with open( variationFile, 'r' ) as variationFp:
        lines = variationFp.readlines()

    res = Container()
    res.variationFileFull = variationFile
    res.variationFile = basename(variationFile)
    res.binCenters   = []
    res.binValues    = []
    res.binValErrors = []

    for line in lines:
        line = line.strip()
        if len(line) == 0: continue

        if fromAgnieszka:
            if line.startswith('#'): continue
            components = line.split()

            if not len(components) == 3:
                print 'Not understanding the following line:'
                print line
                print '    Continuing...'
                continue

            res.binCenters.append(   float(components[0]) )
            res.binValues.append(    float(components[1]) )
            res.binValErrors.append( float(components[2]) )

        elif fromHqT:
            if line.startswith('('): continue
            components = line.split()

            if not len(components) == 6:
                print 'Not understanding the following line:'
                print line
                print '    Continuing...'
                continue

            res.binCenters.append(   float(components[0]) )
            res.binValues.append(    float(components[5]) )


    if fromAgnieszka:
        print '[warning] Dividing by 1/2.27 because of Agnieszka\'s normalization'
        res.binValues = [ i / 2.27 for i in res.binValues ]

    # Try also to get the coupling values from the file name
    variationFile = basename(variationFile)
    if variationFile == 'SM_NNLO':
        res.muR = 1.0
        res.muF = 1.0
        res.Q   = 1.0
    else:

        if fromAgnieszka:
            match = re.search( r'HRes_mR([12h])mF([12h]).*\.top', variationFile )
            if not match:
                print 'Could not determine couplings from \'{0}\''.format( variationFile )
            else:
                res.muR = float( match.group(1).replace( 'h', '0.5' ) )
                res.muF = float( match.group(2).replace( 'h', '0.5' ) )

            match = re.search( r'Q([12h])\.top', variationFile )
            if not match:
                res.Q   = 1.0
            else:
                res.Q   = float( match.group(1).replace( 'h', '0.5' ) )

        elif fromHqT:
            # LO_muR_1_muF_1_Q_1p0_minPt_1_maxPt_51.out
            match = re.search( r'([NLO]+)_muR_([\dpm]+)_muF_([\dpm]+)_Q_([\dpm]+)_minPt_([\dpm]+)_maxPt_([\dpm]+).out', variationFile )
            if not match:
                print 'Could not determine couplings from \'{0}\''.format( variationFile )
            else:
                res.order = match.group(1)
                res.muR   = float( match.group(2).replace('p','.').replace('m','-') )
                res.muF   = float( match.group(3).replace('p','.').replace('m','-') )
                res.Q     = float( match.group(4).replace('p','.').replace('m','-') )
                res.minPt = float( match.group(5).replace('p','.').replace('m','-') )
                res.maxPt = float( match.group(6).replace('p','.').replace('m','-') )


    # Also get some sort of bin boundaries using the heuristic for theory binning
    heuristicBinCenters, heuristicBinBoundaries, heuristicBinWidths = TheoryCommands.BinningHeuristic(
        res.binCenters,
        manualSwitchAt50=( fromAgnieszka or fromHqT )
        )
    res.binBoundaries = heuristicBinBoundaries

    return res



def PlotVariationSpectra(
        variations,
        suffix=None,
        ):

    c.Clear()
    SetCMargins()
    c.SetLogy()

    binBoundaries = variations[0].binBoundaries

    xMin = binBoundaries[0]
    xMax = binBoundaries[-1]

    yMin = max( 0.5 * min([ min(variation.binValues) for variation in variations ]), 0.0000001 )
    yMax = 2.0 * max([ max(variation.binValues) for variation in variations ])


    base = GetPlotBase(
        xMin = xMin, xMax = xMax, yMin = yMin, yMax = yMax,
        xTitle = 'pT [GeV]', yTitle='d#sigma/dp_{T} [pb/GeV]'
        )
    base.Draw('P')



    leg = ROOT.TLegend( 1-CRightMargin-0.4, 1-CTopMargin-0.5, 1-CRightMargin, 1-CTopMargin )
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    colorCycle = itertools.cycle( range(2,9+1) + [ 30, 38, 40, 41, 42 ] + range( 45, 48+1 ) )
    for variation in variations:
        color = next(colorCycle)
        Tg = TheoryCommands.GetTheoryTGraph(
            variation.variationFile,
            variation.binBoundaries,
            variation.binValues,
            boundaries = True,
            )
        Tg.SetLineColor(color)
        Tg.SetMarkerColor(color)
        Tg.Draw('PSAME')

        leg.AddEntry( Tg.GetName(), variation.variationFile, 'l' )


    leg.Draw()


    outname = 'scalevariations'
    if not suffix == None:
        outname += '_' + suffix

    SaveC( outname, asPNG=True )

    c.SetLogy(False)
    yMax = yMax/2.0 * 1.1
    base.SetMaximum(yMax)
    SaveC( outname + '_linear', asPNG=True )


def PlotWithEnvelop(
        *args,
        **kwargs
        ):

    Tgs = []
    for name, variations in args:

        binBoundaries = variations[0].binBoundaries

        if 'ptMax' in kwargs:
            ptMax = kwargs['ptMax']
            for iBinBoundary, binBoundary in enumerate(binBoundaries):
                if ptMax < binBoundary:
                    break
            else:
                print '[error] ptMax {0} is beyond the limit of the bin boundaries ({1})'.format( ptMax, binBoundaries[-1] )
            binBoundaries = binBoundaries[:iBinBoundary+1]


        # Mostly for plotting purposes
        nPoints       = len(binBoundaries)-1
        binCenters    = [ 0.5*(binBoundaries[i]+binBoundaries[i+1]) for i in xrange(nPoints) ]
        halfBinWidths = [ 0.5*(binBoundaries[i+1]-binBoundaries[i]) for i in xrange(nPoints) ]


        # Get pointer to the SM variation
        for variation in variations:
            if variation.muF == 1 and variation.muR == 1 and variation.Q == 1 :
                SMvariation = variation
                break
        else:
            print '[error] Could not find a SM variation in the supplied list'
            return


        minUncs = []
        maxUncs = []

        for iPoint in xrange( nPoints ):

            # Min and max variations
            minVar = min([ variation.binValues[iPoint] for variation in variations ])
            maxVar = max([ variation.binValues[iPoint] for variation in variations ])

            minUncs.append( abs( minVar - SMvariation.binValues[iPoint] ) )
            maxUncs.append( abs( maxVar - SMvariation.binValues[iPoint] ) )


        Tg = ROOT.TGraphAsymmErrors(
            nPoints,
            array( 'd', binCenters ),
            array( 'd', SMvariation.binValues ),
            array( 'd', halfBinWidths ),
            array( 'd', halfBinWidths ),
            array( 'd', minUncs ),
            array( 'd', maxUncs ),
            )
        ROOT.SetOwnership( Tg, False )

        Tg.xMin = binBoundaries[0]
        Tg.xMax = binBoundaries[-1]
        Tg.yMin = min([ center-abs(err) for center, err in zip(SMvariation.binValues, minUncs) ])
        Tg.yMax = max([ center+abs(err) for center, err in zip(SMvariation.binValues, maxUncs) ])
        Tg.yMinPositive = min([ center-abs(err) for center, err in zip(SMvariation.binValues, minUncs) if center-abs(err) > 0. ])

        Tg.SetName(name)
        Tgs.append(Tg)


    # ======================================
    # Make a plot

    xMin = min([ Tg.xMin for Tg in Tgs ])
    if 'ptMax' in kwargs:
        xMax = kwargs['ptMax']
    else:
        xMax = max([ Tg.xMax for Tg in Tgs ])

    # yMinAbs = min([ Tg.yMin for Tg in Tgs ])
    yMinAbs = min([ Tg.yMinPositive for Tg in Tgs ]) # Minimal value that is still positive
    yMaxAbs = max([ Tg.yMax for Tg in Tgs ])
    yMin = yMinAbs - 0.1*( yMaxAbs - yMinAbs )
    yMax = yMaxAbs + 0.1*( yMaxAbs - yMinAbs )

    base = GetPlotBase(
        xMin = xMin, xMax = xMax, yMin = yMin, yMax = yMax,
        xTitle = 'pT [GeV]', yTitle='d#sigma/dp_{T} [pb/GeV]'
        )
    base.Draw('P')

    leg = ROOT.TLegend( 1-CRightMargin-0.4, 1-CTopMargin-0.5, 1-CRightMargin, 1-CTopMargin )
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    colorCycle = itertools.cycle( range(2,9+1) + [ 30, 38, 40, 41, 42 ] + range( 45, 48+1 ) )
    for Tg in Tgs:
        color = next(colorCycle)

        Tg.SetLineColor(color)
        Tg.SetMarkerColor(color)
        Tg.SetFillStyle(1001)
        Tg.SetFillColorAlpha( color, 0.4 )

        # lines, legendDummy = PhysicsCommands.ConvertToLinesAndBoxes( Tg )
        # for line in lines:
        #     line.Draw()
        # legendDummy.Draw('SAME')
        # leg.AddEntry( legendDummy.GetName(), Tg.name, 'l' )

        ConvertTGraphToLinesAndBoxes( Tg, drawImmediately=True, legendObject=leg )

        # Tg.Draw('PSAME E3')

        # leg.AddEntry( Tg.GetName(), Tg.GetName(), 'fpl' )

    leg.Draw()

    outname = 'multipleSpectra'
    if 'suffix' in kwargs:
        outname += '_' + kwargs['suffix']

    SaveC( outname, asPNG=True )

    yMin = 0.5*yMinAbs
    yMax = 2.*yMaxAbs
    base.SetMinimum(yMin)
    base.SetMaximum(yMax)

    c.SetLogy()
    SaveC( outname + '_logscale', asPNG=True )
    c.SetLogy(False)



def PlotRelativeUncertainty(
        variations,
        suffix=None,
        ):
    
    binBoundaries = variations[0].binBoundaries

    # Mostly for plotting purposes
    nPoints       = len(binBoundaries)-1
    binCenters    = [ 0.5*(binBoundaries[i]+binBoundaries[i+1]) for i in xrange(nPoints) ]
    halfBinWidths = [ 0.5*(binBoundaries[i+1]-binBoundaries[i]) for i in xrange(nPoints) ]


    # Get pointer to the SM variation
    for variation in variations:
        if variation.muF == 1 and variation.muR == 1 and variation.Q == 1 :
            SMvariation = variation
            break
    else:
        print '[error] Could not find a SM variation in the supplied list'
        return


    minRelUncs = []
    maxRelUncs = []

    for iPoint in xrange( nPoints ):

        # Min and max variations
        minVar = min([ variation.binValues[iPoint] for variation in variations ])
        maxVar = max([ variation.binValues[iPoint] for variation in variations ])

        minRelUncs.append( -abs( minVar - SMvariation.binValues[iPoint] ) / SMvariation.binValues[iPoint] )
        maxRelUncs.append(  abs( maxVar - SMvariation.binValues[iPoint] ) / SMvariation.binValues[iPoint] )


    c.Clear()
    SetCMargins()

    xMin = SMvariation.binBoundaries[0]
    xMax = SMvariation.binBoundaries[-1]

    yMin = min(minRelUncs) - 0.1*( max(maxRelUncs) - min(minRelUncs) )
    yMax = max(maxRelUncs) + 0.1*( max(maxRelUncs) - min(minRelUncs) )


    base = GetPlotBase(
        xMin = xMin, xMax = xMax, yMin = yMin, yMax = yMax,
        xTitle = 'pT [GeV]', yTitle='#Delta d#sigma/dp_{T} [pb/GeV]'
        )
    base.Draw('P')


    Tg = ROOT.TGraphAsymmErrors(
        nPoints,
        array( 'd', binCenters ),
        array( 'd', [ 0. for i in xrange(nPoints) ] ),
        array( 'd', halfBinWidths ),
        array( 'd', halfBinWidths ),
        array( 'd', [ -i for i in minRelUncs ] ),
        array( 'd', maxRelUncs ),
        )

    Tg.Draw('PSAME')


    outname = 'scaleRelativeUncertainty'
    if not suffix == None:
        outname += '_' + suffix

    SaveC( outname, asPNG=True )



def MergeLowHighPtFilesFromHqT(
        variations_lowPt,
        variations_highPt,
        ):

    HqTvariations = []
    for variation_lowPt in variations_lowPt:
        firstPartOfFilename_lowPt = variation_lowPt.variationFile.split('minPt')[0]

        # Find corresponding highPt variation
        for variation_highPt in variations_highPt:
            firstPartOfFilename_highPt = variation_highPt.variationFile.split('minPt')[0]
            if firstPartOfFilename_highPt == firstPartOfFilename_lowPt:
                break
        else:
            print 'Could not find matching high Pt file for \'{0}\''.format( variation_lowPt.variationFile )
            continue

        variation = deepcopy( variation_lowPt )
        variation.binCenters = variation_lowPt.binCenters[:-1] + variation_highPt.binCenters
        variation.binValues  = variation_lowPt.binValues[:-1] + variation_highPt.binValues

        # Also get some sort of bin boundaries using the heuristic for theory binning
        heuristicBinCenters, heuristicBinBoundaries, heuristicBinWidths = TheoryCommands.BinningHeuristic(
            variation.binCenters,
            manualSwitchAt50=True
            )
        variation.binBoundaries = heuristicBinBoundaries

        HqTvariations.append( variation )

    return HqTvariations



def ConvertTGraphToLinesAndBoxes(
        Tg,
        drawImmediately=False,
        legendObject=None,
        verbose=False,
        noBoxes=False,
        xMaxExternal=None,
        yMinExternal=None,
        ):

    yBand = ( Tg.GetErrorYhigh(0) != -1 and Tg.GetErrorYlow(0) != -1 )
    xBand = ( Tg.GetErrorXhigh(0) != -1 and Tg.GetErrorXlow(0) != -1 )

    if not xBand or not yBand:
        Commands.ThrowError( 'Make sure all errors are filled' )

    nPoints = Tg.GetN()

    if verbose: print '[debug] Converting \'{0}\' to lines and boxes ({1} points)'.format( Tg.GetName(), nPoints )

    lineStyle = Tg.GetLineStyle()
    lineWidth = Tg.GetLineWidth()
    lineColor = Tg.GetLineColor()

    fillStyle = Tg.GetFillStyle()
    fillColor = Tg.GetFillColor()

    lines = []
    boxes = []

    x_Double = ROOT.Double(0)
    y_Double = ROOT.Double(0)
    for iPoint in xrange( nPoints ):
        Tg.GetPoint( iPoint, x_Double, y_Double )
        x = float(x_Double)
        y = float(y_Double)

        xMin = abs(Tg.GetErrorXlow(iPoint))
        xMax = abs(Tg.GetErrorXhigh(iPoint))
        yMin = abs(Tg.GetErrorYlow(iPoint))
        yMax = abs(Tg.GetErrorYhigh(iPoint))

        if verbose:
            print '[debug] Point {0:<3}:'.format( iPoint )
            print '        x = {0:+8.3f}, xMin = {1:+8.3f}, xMax = {2:+8.3f}'.format( x, xMin, xMax )
            print '        y = {0:+8.3f}, yMin = {1:+8.3f}, yMax = {2:+8.3f}'.format( y, yMin, yMax )

        if xMaxExternal != None:
            if x-xMin > xMaxExternal:
                if verbose:
                    print 'x-xMin {0} > xMaxExternal {1}; Continuing'.format( x-xMin, xMaxExternal )
                continue
            elif x+xMax > xMaxExternal:
                if verbose:
                    print 'x+xMax {0} > xMaxExternal {1}; Limiting xMax'.format( x+xMax, xMaxExternal )
                xMax = xMaxExternal-x

        # if yMinExternal != None:
        #     if x-xMin > yMinExternal:
        #         if verbose:
        #             print 'x-xMin {0} > yMinExternal {1}; Continuing'.format( x-xMin, yMinExternal )
        #         continue
        #     elif x+yMin > yMinExternal:
        #         if verbose:
        #             print 'x+yMin {0} > yMinExternal {1}; Limiting yMin'.format( x+yMin, yMinExternal )
        #         yMin = yMinExternal-x

        line = ROOT.TLine( x-xMin, y, x+xMax, y )
        ROOT.SetOwnership( line, False )
        line.SetLineStyle( lineStyle )
        line.SetLineWidth( lineWidth )
        line.SetLineColor( lineColor )
        lines.append(line)

        box  = ROOT.TBox( x-xMin, y-yMin, x+xMax, y+yMax )
        ROOT.SetOwnership( box, False )
        box.SetLineWidth(0)
        if not fillStyle == 0: box.SetFillStyle( fillStyle )
        box.SetFillColor( fillColor )
        boxes.append(box)

        if drawImmediately:

            if yMinExternal is None:
                if not noBoxes: box.Draw()
                line.Draw()
            else:
                if y < yMinExternal:
                    if y+yMax > yMinExternal and not noBoxes:
                        box.SetY1( yMinExternal )
                        box.Draw()
                elif y-yMin < yMinExternal:
                    if not noBoxes:
                        box.SetY1( yMinExternal )
                        box.Draw()
                    line.Draw()
                else:
                    if not noBoxes: box.Draw()
                    line.Draw()


    # Create a dummy for the legend
    legendDummy = ROOT.TGraph( 1, array( 'd', [-999.]), array( 'd', [-999.]) )
    ROOT.SetOwnership( legendDummy, False )
    legendDummy.SetLineColor( lineColor )
    legendDummy.SetLineWidth( Tg.GetLineWidth() )
    legendDummy.SetFillColor( fillColor )
    if not fillStyle == 0: legendDummy.SetFillStyle( fillStyle )
    legendDummy.SetMarkerStyle(8)
    legendDummy.SetMarkerColor( lineColor )
    legendDummy.SetMarkerSize(0)
    legendDummy.SetName( Tg.GetName() + '_dummy' )

    if drawImmediately:
        legendDummy.Draw('PSAME')

        if legendObject:
            legendObject.AddEntry(
                legendDummy.GetName(), Tg.GetName(),
                'lf' if not noBoxes else 'l'
                )



def ConvertTGraphToHistogram(
        Tg,
        drawImmediately=False,
        legendObject=None,
        verbose=False,
        ):

    yBand = ( Tg.GetErrorYhigh(0) != -1 and Tg.GetErrorYlow(0) != -1 )
    xBand = ( Tg.GetErrorXhigh(0) != -1 and Tg.GetErrorXlow(0) != -1 )

    if not xBand or not yBand:
        Commands.ThrowError( 'Make sure all errors are filled' )

    nPoints = Tg.GetN()

    if verbose: print '[debug] Converting \'{0}\' to lines and boxes ({1} points)'.format( Tg.GetName(), nPoints )

    lineStyle = Tg.GetLineStyle()
    lineWidth = Tg.GetLineWidth()
    lineColor = Tg.GetLineColor()

    fillStyle = Tg.GetFillStyle()
    fillColor = Tg.GetFillColor()

    binBoundaries = []
    yValues       = []

    x_Double = ROOT.Double(0)
    y_Double = ROOT.Double(0)
    for iPoint in xrange( nPoints ):
        Tg.GetPoint( iPoint, x_Double, y_Double )
        x = float(x_Double)
        y = float(y_Double)

        xMin = abs(Tg.GetErrorXlow(iPoint))
        xMax = abs(Tg.GetErrorXhigh(iPoint))
        yMin = abs(Tg.GetErrorYlow(iPoint))
        yMax = abs(Tg.GetErrorYhigh(iPoint))

        if verbose:
            print '[debug] Point {0:<3}:'.format( iPoint )
            print '        x = {0:+8.3f}, xMin = {1:+8.3f}, xMax = {2:+8.3f}'.format( x, xMin, xMax )
            print '        y = {0:+8.3f}, yMin = {1:+8.3f}, yMax = {2:+8.3f}'.format( y, yMin, yMax )

        binBoundaries.append( x - xMin )
        rightBound = x + xMax

        yValues.append( y )


    binBoundaries.append( rightBound )


    Hname = TheoryCommands.GetUniqueRootName()
    H = ROOT.TH1F(
        Hname, Hname,
        nPoints,
        array( 'f', binBoundaries )
        )
    ROOT.SetOwnership( H, False )

    H.SetLineStyle( lineStyle )
    H.SetLineWidth( lineWidth )
    H.SetLineColor( lineColor )

    print '[debug] Filling the histogram'
    for iBin in xrange(nPoints):
        H.SetBinContent( iBin+1, yValues[iBin] )

    # H.Scale( 1./H.Integral() )

    if drawImmediately:
        H.Draw('SAME')

        if legendObject:
            legendObject.AddEntry( H.GetName(), Tg.GetName(), 'l' )







########################################
# End of Main
########################################
if __name__ == "__main__":
    main()