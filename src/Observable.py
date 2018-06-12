import Commands
from copy import deepcopy
from array import array
import ROOT
import logging

#____________________________________________________________________
class Observable(object):

    YR4_totalXS = 0.

    """Simple class that returns shape, cross section, or cross section over bin width"""
    def __init__(self,
            name,
            title,
            shape,
            binning,
            lastBinIsOverflow = True,
            ):
        if not len(binning) == len(shape)+1:
            raise RuntimeError( 'len(binning)[={0}] =/= len(shape)+1[={1}]'.format( len(binning), len(shape) ) )
        self.name = name
        self.title = title
        self.shape = shape
        self.binning = binning
        self.lastBinIsOverflow = lastBinIsOverflow
        self.nBins = len(self.shape)

        # Set empty uncertainties
        self.unc_fraction = [ 0. for s in shape ]

    #____________________________________________________________________
    def drop_first_bin(self):
        self.shape = self.shape[1:]
        self.binning = self.binning[1:]
        self.unc_fraction = self.unc_fraction[1:]
        self.nBins -= 1
        logging.debug(
            'Obs {0}: First bin dropped, remaining bins: {1}'
            .format(self.name, self.binning)
            )

    #____________________________________________________________________
    def drop_last_bin(self):
        self.shape = self.shape[:-1]
        self.binning = self.binning[:-1]
        self.unc_fraction[:-1]
        self.nBins -= 1
        self.lastBinIsOverflow=False
        logging.debug(
            'Obs {0}: Last bin dropped, remaining bins: {1}'
            .format(self.name, self.binning)
            )

    #____________________________________________________________________
    def drop_bins_up_to_value(self, value):
        logging.info(
            'Obs {0}: Dropping bins up to value {1}; current binning is {2}'
            .format(self.name, value, self.binning)
            )
        while self.binning[-1] > value:
            self.drop_last_bin()
        logging.info('Obs {0}: Final binning: {1}'.format(self.name, self.binning))

    #____________________________________________________________________
    def crosssection(self):
        return [ s * self.YR4_totalXS for s in self.shape ]

    #____________________________________________________________________
    def inclusive_crosssection(self):
        return sum(self.crosssection())

    #____________________________________________________________________
    def crosssection_over_binwidth(self, normalize_by_second_to_last_bin_width=False):
        xs = self.crosssection()
        if self.lastBinIsOverflow:
            xs_o_binwidth = [ xs[i] / ( self.binning[i+1]-self.binning[i] ) for i in xrange(self.nBins-1) ] + [ xs[-1] ]
        else:
            xs_o_binwidth = [ xs[i] / ( self.binning[i+1]-self.binning[i] ) for i in xrange(self.nBins) ]

        if normalize_by_second_to_last_bin_width:
            xs_o_binwidth[-1] = xs[-1] / ( self.binning[-2]-self.binning[-3] )

        return xs_o_binwidth

    #____________________________________________________________________
    def unc_xs_over_binwidth(self, normalize_by_second_to_last_bin_width=False):
        xss = self.crosssection_over_binwidth(normalize_by_second_to_last_bin_width)
        return [ xs * unc for xs, unc in zip(xss, self.unc_fraction) ]

    #____________________________________________________________________
    def mergeBins( self, mergeList ):
        newshape   = []
        newbinning = [ self.binning[0] ]
        for bins in mergeList:
            if isinstance( bins, int ):
                newshape.append( self.shape[bins] )
                newbinning.append( self.binning[bins+1] ) # Always add right bound of bin
            else:
                shapesum = 0.
                for iBin in bins:
                    shapesum += self.shape[iBin]
                newshape.append(shapesum)
                newbinning.append( self.binning[bins[-1]+1] ) # Always add right bound of bin

        if not len(newbinning) == len(newshape)+1:
            raise RuntimeError( 'len(newbinning)[={0}] =/= len(newshape)+1[={1}]'.format( len(newbinning), len(newshape) ) )

        copy = deepcopy(self)
        copy.shape = newshape
        copy.binning = newbinning
        copy.nBins = len(newshape)
        return copy

    #____________________________________________________________________
    def Print( self ):
        xs = self.crosssection()
        xs_over_binwidth = self.crosssection_over_binwidth()
        strList = lambda L: ', '.join([ '{0:<9.3f}'.format(f) for f in L ])
        print '\n Cross sections for observable {0} ({1})'.format( self.name, self.title )
        print 'shape:            ' + strList(self.shape)
        print 'binning:          ' + strList(self.binning)
        print 'xs:               ' + strList(xs)
        print 'xs_over_binwidth: ' + strList(xs_over_binwidth)
        print 'inclusive xs:     ' + str(self.inclusive_crosssection())

    #____________________________________________________________________
    def BasicHistogram( self, what=None ):

        if what is None:
            what = 'crosssection_over_binwidth'
        if not what in [ 'crosssection', 'crosssection_over_binwidth', 'shape' ]:
            Commands.ThrowError( 'Choose from {0}'.format([ 'crosssection', 'crosssection_over_binwidth', 'shape' ]) )

        if what == 'crosssection':
            ys = self.crosssection()
        elif what == 'crosssection_over_binwidth':
            ys = self.crosssection_over_binwidth()
            # if self.lastBinIsOverflow:
            #     ys[-1] /= self.binning[-2] - self.binning[-3]

        elif what == 'shape':
            ys = self.shape

        uniqueName = 'ObsHist_{0}'.format( Commands.__uniqueid__().next() )
        H = ROOT.TH1F(
            uniqueName, uniqueName,
            len(ys), array( 'f', self.binning )
            )
        ROOT.SetOwnership( H, False )
        H.SetLineWidth(2)

        for iBin in xrange(len(ys)):
            H.SetBinContent( iBin+1, ys[iBin] )

        return H
