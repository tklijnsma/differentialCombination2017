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

    parser.add_argument( '--derivedTheoryFilePlots_Top',              action=CustomAction )


########################################
# Methods
########################################    

def main( args ):

    TheoryCommands.SetPlotDir( 'plots_{0}_onetimeplots'.format(datestr) )


    #____________________________________________________________________
    if args.derivedTheoryFilePlots_Top:

        DO_PARAMETRIZATION = True
        # DO_PARAMETRIZATION = False

        # ======================================
        # Input

        readList = lambda theoryFiles: [ TheoryFileInterface.ReadDerivedTheoryFile( theoryFile ) for theoryFile in theoryFiles ]


        TheoryFileInterface.SetFileFinderDir( LatestPaths.derivedTheoryFilesDirectory_Top )

        theoryFiles_ct_cg = TheoryFileInterface.FileFinder(
            ct='*', cg='*', muR=1, muF=1, Q=1, filter='cb'
            )
        containers_ct_cg = readList( theoryFiles_ct_cg )


        theoryFiles_ct_cb = TheoryFileInterface.FileFinder(
            ct='*', cb='*', muR=1, muF=1, Q=1, filter='cg'
            )
        containers_ct_cb = readList( theoryFiles_ct_cb )


        theoryFiles_kappab_kappac = TheoryFileInterface.FileFinder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFilesDirectory_YukawaSummed
            )
        containers_kappab_kappac = readList( theoryFiles_kappab_kappac )


        theoryFiles_kappab_kappac_Gluon = TheoryFileInterface.FileFinder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFilesDirectory_YukawaGluonInduced
            )
        containers_kappab_kappac_Gluon = readList( theoryFiles_kappab_kappac_Gluon )


        theoryFiles_kappab_kappac_Quark = TheoryFileInterface.FileFinder(
            kappab='*', kappac='*', filter='muR',
            directory = LatestPaths.derivedTheoryFilesDirectory_YukawaQuarkInduced
            )
        containers_kappab_kappac_Quark = readList( theoryFiles_kappab_kappac_Quark )


        theoryFiles_kappab_kappac_QuarkScaled = TheoryFileInterface.FileFinder(
            kappab='*', kappac='*', muR=1, muF=1, Q=1,
            directory = LatestPaths.derivedTheoryFilesDirectory_YukawaQuarkInducedScaled
            )
        containers_kappab_kappac_QuarkScaled = readList( theoryFiles_kappab_kappac_QuarkScaled )


        def kappaSorter( container ):
            return ( container.kappac, container.kappab )
        containers_kappab_kappac.sort( key=kappaSorter )



        ct_cg = Container()
        ct_cg.name = 'kappat_kappag'
        ct_cg.couplings = [ 'ct', 'cg' ]
        ct_cg.containers = containers_ct_cg
        ct_cg.ws = LatestPaths.ws_combined_split_top

        ct_cb = Container()
        ct_cb.name = 'kappat_kappab'
        ct_cb.couplings = [ 'ct', 'cb' ]
        ct_cb.containers = containers_ct_cb

        kappab_kappac_Gluon = Container()
        kappab_kappac_Gluon.name = 'kappab_kappac_Gluon'
        kappab_kappac_Gluon.couplings = [ 'kappab', 'kappac' ]
        kappab_kappac_Gluon.containers = containers_kappab_kappac_Gluon

        kappab_kappac_Quark = Container()
        kappab_kappac_Quark.name = 'kappab_kappac_Quark'
        kappab_kappac_Quark.couplings = [ 'kappab', 'kappac' ]
        kappab_kappac_Quark.containers = containers_kappab_kappac_Quark
        kappab_kappac_Quark.parametrizeByFitting = True

        kappab_kappac_QuarkScaled = Container()
        kappab_kappac_QuarkScaled.name = 'kappab_kappac_QuarkScaled'
        kappab_kappac_QuarkScaled.couplings = [ 'kappab', 'kappac' ]
        kappab_kappac_QuarkScaled.containers = containers_kappab_kappac_QuarkScaled


        couplingVariationContainers = [
            # ct_cg,
            # ct_cb,
            # kappab_kappac_Gluon,
            kappab_kappac_Quark,
            # kappab_kappac_QuarkScaled,
            ]

        n_theoryFiles_kappab_kappac = len( theoryFiles_kappab_kappac )
        for iChunk, chunk in enumerate( chunks( containers_kappab_kappac, n_theoryFiles_kappab_kappac/4+1 ) ):
            couplingVariationContainer = Container()
            couplingVariationContainer.name = 'kappab_kappac_{0}'.format(iChunk)
            couplingVariationContainer.couplings = [ 'kappab', 'kappac' ]
            couplingVariationContainer.containers = list(chunk)
            couplingVariationContainer.ws = LatestPaths.ws_combined_split_yukawa
            # couplingVariationContainers.append( couplingVariationContainer )


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
                    container.Tg_parametrization = wsParametrization.GetOutputContainer( **kwargs ).Tg

                isParametrized = True

            if DO_PARAMETRIZATION and hasattr( couplingVariationContainer, 'parametrizeByFitting' ):

                for container in containers:
                    couplingVals = [ float(getattr( container, coupling )) for coupling in couplingVariationContainer.couplings ]
                    if all( [ val == 1.0 for val in couplingVals ] ):
                        SM = container
                        break
                else:
                    Commands.ThrowError( 'Could not find standard model' )
                    sys.exit()

                parametrization = Parametrization()
                parametrization.ParametrizeByFitting( containers )

                for container in containers:
                    # kwargs = { coupling : float(getattr( container, coupling )) for coupling in couplingVariationContainer.couplings }
                    for coupling in couplingVariationContainer.couplings:
                        setattr( parametrization, coupling, float(getattr( container, coupling )) )
                    parametrizationResult = parametrization.GetOutputContainer()
                    parametrizationResult.ratios = [ xs / SMxs if not SMxs == 0. else 0. for xs, SMxs in zip( parametrizationResult.binValues, SM.crosssection ) ]
                    container.Tg_parametrization = parametrizationResult.GetTGraph( yAttr='ratios' )

                isParametrized = True


            c.Clear()
            SetCMargins()

            yMinAbs = min( [ min(container.ratios) for container in containers ] )
            yMaxAbs = max( [ max(container.ratios) for container in containers ] )
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
                yTitle = 'Ratio w.r.t. ggH^{SM}',
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
                    name = variationName
                    )
                Tg.SetMarkerColor(color)
                Tg.SetMarkerStyle(8)
                Tg.SetMarkerSize(0.9)
                Tg.Draw('PX')

                if isParametrized:
                    container.Tg_parametrization.SetLineColor(color)
                    container.Tg_parametrization.SetLineWidth(2)
                    container.Tg_parametrization.Draw('LX')

                leg.AddEntry( Tg.GetName(), variationTitle, 'p' )

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