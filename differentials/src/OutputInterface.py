#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, tempfile, shutil, re, glob, itertools, sys, numpy, operator, pprint, re, random
from os.path import *
from operator import itemgetter
from array import array
from math import log, exp, sqrt, copysign
from copy import deepcopy
from glob import glob

from time import strftime
datestr = strftime( '%b%d' )

from Container import Container
import Commands
import TheoryCommands

import ROOT


########################################
# General functions
########################################

class OutputContainer(Container):
    def __init__( self, *args ):
        super(OutputContainer, self).__init__()
        
        if len(args) == 1:
            if isinstance( args[0], Container ):
                self.read_derived_theory_file_container( args[0] )




    def get_TGraph(
            self,
            name=None,
            xAttr=None,
            yAttr=None,
            xAttrErrDown=None,
            xAttrErrUp=None,
            yAttrErrDown=None,
            yAttrErrUp=None,
            xAreBinBoundaries=True,
            ):
        self.xAreBinBoundaries = xAreBinBoundaries


        if xAttr is None:
            for attempt in [ 'binBoundaries' ]:
                if hasattr( self, attempt ):
                    xAttr = attempt
                    break
            else:
                Commands.throw_error( 'Can not determine an attribute to use for the x axis' )
                sys.exit()

        if yAttr is None:
            for attempt in [ 'binValues', 'mus', 'ratios', 'crosssection' ]:
                if hasattr( self, attempt ):
                    yAttr = attempt
                    break
            else:
                Commands.throw_error( 'Can not determine an attribute to use for the x axis' )
                sys.exit()

        xValues = getattr( self, xAttr )
        yValues = getattr( self, yAttr )

        if self.xAreBinBoundaries:
            if len(xValues)-1 > len(yValues):
                Commands.warning( 'Found len(xValues) = {0}, len(yValues) = {1}; need to limit number of bins to make a TGraph'.format( len(xValues), len(yValues) ) )
                xValues = xValues[:len(yValues)+1]
            elif len(xValues)-1 < len(yValues):
                Commands.throw_error(
                    'len(xValues)-1 = {0}, len(yValues) = {1} ; should be the same'.format( len(xValues)-1, len(yValues) ),
                    throwException = True
                    )

            nBins = len(xValues)-1
            binBoundaries = xValues
            binCenters    = [ 0.5*(binBoundaries[i]+binBoundaries[i+1]) for i in xrange(nBins) ]
            binWidths     = [ (binBoundaries[i+1]-binBoundaries[i]) for i in xrange(nBins) ]
            halfBinWidths = [ 0.5*(binBoundaries[i+1]-binBoundaries[i]) for i in xrange(nBins) ]

        else:
            if len(xValues) != len(yValues):
                Commands.throw_error(
                    'len(xValues) = {0}, len(yValues) = {1} ; should be the same'.format( len(xValues), len(yValues) ),
                    throwException = True
                    )

            nBins = len(xValues)
            binCenters    = xValues
            
            binWidths     = [ xValues[i+1] - xValues[i] for i in xrange(nBins-1) ]
            binWidths     = [ binWidths[0] ] + binWidths + [ binWidths[-1] ]
            halfBinWidths = [ 0.5*binWidth for binWidth in binWidths ]

            binBoundaries = [ binCenters[0]-halfBinWidths[0] ] + [ center + halfBinWidth for center, halfBinWidth in zip( binCenters, halfBinWidths ) ]


        if xAttrErrDown is None:
            xErrDown = [ 0 for i in xrange(nBins) ]
        else:
            xErrDown = getattr( self, xAttrErrDown )
        if xAttrErrUp is None:
            xErrUp   = [ 0 for i in xrange(nBins) ]
        else:
            xErrUp = getattr( self, xAttrErrUp )
        if yAttrErrDown is None:
            yErrDown = [ 0 for i in xrange(nBins) ]
        else:
            yErrDown = getattr( self, yAttrErrDown )
        if yAttrErrUp is None:
            yErrUp   = [ 0 for i in xrange(nBins) ]
        else:
            yErrUp = getattr( self, yAttrErrUp )


        Tg = ROOT.TGraphAsymmErrors(
            len(binCenters),
            array( 'd', binCenters ),
            array( 'd', yValues ),
            array( 'd', halfBinWidths ),
            array( 'd', halfBinWidths ),
            array( 'd', [ abs(l) for l in yErrDown ] ),
            array( 'd', [ abs(r) for r in yErrUp ] ),
            )

        ROOT.SetOwnership( Tg, False )
        Tg.SetLineWidth(2)

        # if name is None:
        #     if not hasattr( self, 'name' ):
        #         Commands.throw_error( 'Need to specify a name.' )
        #         sys.exit()
        #     name = self.name
        # Tg.SetName( name )

        if name is None:
            if hasattr( self, 'name' ):
                Tg.SetName( self.name )
        else:
            Tg.SetName( name )


        # ======================================
        # Set additional attributes for convenience

        Tg.name          = name
        Tg.xValues       = xValues
        Tg.binValues     = yValues
        Tg.yValues       = yValues
        Tg.binCenters    = binCenters
        Tg.binWidths     = binWidths
        Tg.binBoundaries = binBoundaries
        # Tg.binErrLeft    = binErrLeft
        # Tg.binErrRight   = binErrRight

        Tg.xMin = min( binBoundaries )
        Tg.xMax = max( binBoundaries )
        if len(binBoundaries) >= 5:
            Tg.fourth_xMin = sorted( binBoundaries )[4]
            Tg.fourth_xMax = sorted( binBoundaries, reverse=True )[4]
        else:
            Tg.fourth_xMin = sorted( binBoundaries )[-1]
            Tg.fourth_xMax = sorted( binBoundaries, reverse=True )[-1]

        if yAttrErrUp is None and yAttrErrDown is None:
            Tg.yMin = min( yValues )
            Tg.yMax = max( yValues )
            if len(yValues) >= 5:
                Tg.fourth_yMin = sorted( yValues )[4]
                Tg.fourth_yMax = sorted( yValues, reverse=True )[4]
            else:
                Tg.fourth_yMin = sorted( yValues )[-1]
                Tg.fourth_yMax = sorted( yValues, reverse=True )[-1]

        self.Tg = Tg
        return Tg


    def read_derived_theory_file_container(
            self,
            theoryContainer,
            ):
        for attr in theoryContainer.list_attributes():
            setattr( self, attr, getattr( theoryContainer, attr ) )





# def output_container_from_wsparametrization(
#         wsFile, returnWhat='theory', **kwargs
#         ):

    
    










########################################
# Execute this if called directly
########################################

if __name__ == '__main__':
    print 'Import functions from this file'
    print 'Available functions:'
    functions = [ name for (name, thing) in locals().items() if callable(thing) ]
    print '\n'.join( functions )