#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, re, itertools
from os.path import *
from glob import glob
from sys import exit
from numpy import corrcoef, var, std

from time import strftime
datestr = strftime( '%b%d' )


import TheoryCommands


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

CPlotDir = 'plots_CorrelationMatrices_{0}'.format(datestr)
def SaveC( outname, subdir=None, asPDF=True, asPNG=False, asROOT=False ):

    if subdir == None:
        outdir = CPlotDir
    else:
        outdir = join( CPlotDir, subdir )

    if not isdir(outdir): os.makedirs( outdir )
    outname = join( outdir, outname.replace('.pdf','').replace('.png','').replace('.root','') )

    if asPDF: c.SaveAs( outname + '.pdf' )
    if asPNG: c.SaveAs( outname + '.png' )
    if asROOT: c.SaveAs( outname + '.proot' )

class Container:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


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



def GetCorrelationMatrix(
        variations,
        makeScatterPlots = False,
        makeCorrelationMatrixPlot = True,
        outname = 'corrMat',
        verbose = False,
        halfNumberOfPlots = False,
        ):

    if not makeScatterPlots:
        halfNumberOfPlots = True

    # ======================================
    # Try to get correlations

    nPoints = len( variations[0].binValues )

    if makeScatterPlots:
        base = GetPlotBase()

    corrMatrix = [ [ 999 for jPoint in xrange(nPoints) ] for iPoint in xrange(nPoints) ]

    jStartPoint = 0
    for iPoint in xrange(nPoints):
        if verbose: print 'Processing point {0}/{1}'.format( iPoint, nPoints-1 )
        iBinValues = [ variation.binValues[iPoint] for variation in variations ]
        iPt = variations[0].binCenters[iPoint]

        # Essentially doing this loop from 0 to nPoints is double work, but for the plots
        # it's nicer to do it anyway
        if halfNumberOfPlots: jStartPoint = iPoint
        for jPoint in xrange(jStartPoint,nPoints):
            jBinValues = [ variation.binValues[jPoint] for variation in variations ]
            jPt = variations[0].binCenters[jPoint]

            corr = corrcoef( iBinValues, jBinValues )[0][1]

            corrMatrix[iPoint][jPoint] = corr
            
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
                ' '.join([ '{0:+.2f}'.format(corrMatrix[iPoint][jPoint]) for jPoint in xrange(nPoints) ])
                + '\n'
                )

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

        for iPoint in xrange(nPoints):
            if nPoints < 20 or iPoint % int(0.1*nPoints) == 0:
                if not hasattr( variations[0], 'binBoundaries' ):
                    binLabel = '{0:.1f}'.format(variations[0].binCenters[iPoint])
                else:
                    binLabel = '{0:.1f} - {1:.1f}'.format(
                        variations[0].binBoundaries[iPoint],
                        variations[0].binBoundaries[iPoint+1]
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




def ReadVariationFile( variationFile ):

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

    # Try also to get the coupling values from the file name
    variationFile = basename(variationFile)
    if variationFile == 'SM_NNLO':
        res.muR = 1.0
        res.muF = 1.0
        res.Q   = 1.0
    else:
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

    return res



def PlotVariationSpectra(
        variations,
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




    SaveC( 'scalevariations', asPNG=True )

    c.SetLogy(False)




########################################
# End of Main
########################################
if __name__ == "__main__":
    main()