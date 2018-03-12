import ROOT
import plotting_utils as utils
import pywrappers
from canvas import c

import differentials
import differentials.core as core
import plots

import logging
from math import isnan, isinf, log10, sqrt
from array import array
from collections import namedtuple


class CorrelationMatrixPlot(plots.PlotBase):

    def __init__(self, plotname, ws=None):
        super(CorrelationMatrixPlot, self).__init__(plotname)
        self.ws = ws

    def read_ws(self):
        self.w = core.get_ws(self.ws)
        self.POIs = core.read_set(w, 'POI')

    def correlation_matrix(self):
        pass
        # HIER VERDER





# #____________________________________________________________________
# def PlotCorrelationMatrix(
#         container
#         ):
#     numpy.set_printoptions( precision=2, linewidth=100 )

#     POIs = Commands.ListPOIs(container.ws)
#     Commands.SortPOIs(POIs)
#     nBins = len(POIs)

#     # Obtain correlation matrix
#     corrMat = [ [ 999 for j in xrange(nBins) ] for i in xrange(nBins) ]
#     with Commands.OpenRootFile( container.corrRootFile ) as rootFp:
#         fit = rootFp.Get('fit')
#         for iBin1, POI1 in enumerate( POIs ):
#             for iBin2, POI2 in enumerate( POIs ):
#                 corrMat[iBin1][iBin2] = fit.correlation( POI1, POI2 )

#     print 'Found the following corrMat from', container.corrRootFile
#     print numpy.array(corrMat)
    

#     # ======================================
#     # Make plot

#     c.Clear()
#     SetCMargins(
#         TopMargin   = 0.08,
#         RightMargin = 0.14,
#         BottomMargin = 0.17,
#         )

#     titleDict = {
#         'PTH' : 'p_{T}^{H} (GeV)',
#         }
#     productionMode, observableName, _ = Commands.InterpretPOI(POIs[0])
#     observableName = titleDict.get( observableName, observableName )
#     xTitle = getattr( container, 'xTitle', observableName )

#     # Construct the binning labels
#     def toStr( number ):
#         if number == '-INF':
#             string = '-#infty'
#         elif number == 'INF':
#             string = '#infty'
#         elif number.is_integer():
#             string = '{0:d}'.format(int(number))
#         else:
#             string = '{0:0.2f}'.format(number)
#         return string

#     binningLabels = []
#     for POI in POIs:
#         _1, _2, binBoundaries = Commands.InterpretPOI(POI)
#         binBoundariesAsStrs = [ toStr(i) for i in binBoundaries ]
#         if len(binBoundaries) == 1:
#             binningLabels.append( binBoundariesAsStrs[0] )
#         elif len(binBoundaries) == 2:
#             binningLabels.append( '(' + ', '.join(binBoundariesAsStrs) + ')'  )


#     H = ROOT.TH2D(
#         GetUniqueRootName(),
#         # '#scale[0.85]{{Bin-to-bin correlation matrix for {0}}}'.format(observableName),
#         '',
#         nBins, 0., nBins,
#         nBins, 0., nBins
#         )
#     ROOT.SetOwnership( H, False )
#     H.SetContour(100)

#     for iRow in xrange(nBins):
#         for iCol in xrange(nBins):
#             H.SetBinContent( iCol+1, iRow+1, corrMat[iRow][iCol] )
#     H.GetZaxis().SetRangeUser(-1.0,1.0)

#     for iBin in xrange(nBins):
#         H.GetXaxis().SetBinLabel( iBin+1, binningLabels[iBin] )
#         H.GetYaxis().SetBinLabel( iBin+1, binningLabels[iBin] )
#     H.GetXaxis().SetTitle( xTitle )
#     H.GetXaxis().SetTitleSize(0.05)
#     H.GetXaxis().SetTitleOffset( 1.6 )
#     H.GetXaxis().SetLabelSize(0.045)
#     H.GetYaxis().SetLabelSize(0.045)

#     H.Draw('COLZ TEXT')


#     # ======================================
#     # Set some style

#     ROOT.gStyle.SetHistMinimumZero() # To draw the "0", otherwise ROOT leaves it empty
#     ROOT.gStyle.SetPaintTextFormat('1.2g')

#     n_stops = 3
#     stops  = [ 0.0, 0.5, 1.0 ]
#     reds   = [ 0.0, 1.0, 1.0 ]
#     blues  = [ 1.0, 1.0, 0.0 ]
#     greens = [ 0.0, 1.0, 0.0 ]

#     ROOT.TColor.CreateGradientColorTable(
#         n_stops,
#         array('d', stops ),
#         array('d', reds ),
#         array('d', greens ),
#         array('d', blues ),
#         255 )

#     Commands.GetCMSLabel()
#     Commands.GetCMSLumi()

#     plotname = 'corrMat_' + basename(container.corrRootFile).replace('/','').replace('higgsCombine_','').replace('higgsCombine','').replace('.root','').replace('.','_')
#     SaveC( plotname )

#     # Set back to default
#     numpy.set_printoptions( precision=8, linewidth=75 )

