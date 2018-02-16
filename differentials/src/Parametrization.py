import os, tempfile, shutil, re, glob, itertools, sys, numpy, operator, pprint, re, random
from os.path import *
from operator import itemgetter
from array import array
from math import log, exp, sqrt, copysign
from copy import deepcopy
from glob import glob

from time import strftime
datestr = strftime( '%b%d' )

import Commands
import TheoryCommands
from OutputInterface import OutputContainer

import ROOT
ROOT.gROOT.SetBatch(True)
# ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
# ROOT.RooMsgService.instance().setSilentMode(True)


class WSParametrization():
    def __init__( self, wsFile=None, returnWhat='theory', verbose=False, newStyle=False ):
        self.verbose = verbose
        self.Persistence = []
        self.newStyle = newStyle

        if not returnWhat in [ 'both', 'theory', 'exp' ]:
            Commands.ThrowError( 'Argument \'returnWhat\' should be \'both\', \'theory\' or \'exp\'' )
            sys.exit()
        self.returnWhat = returnWhat

        if self.returnWhat == 'theory':
            self.returnTheory = True
            self.returnExp    = False
        elif self.returnWhat == 'exp':
            self.returnTheory = False
            self.returnExp    = True
        elif self.returnWhat == 'both':
            self.returnTheory = True
            self.returnExp    = True

        if wsFile != None:
            self.Parametrize( wsFile, verbose = self.verbose )


    def Parametrize( self, wsFile, verbose=False ):

        self.file = wsFile

        wsFp = ROOT.TFile.Open( wsFile )
        w = wsFp.Get('w')

        couplings = []
        couplingsList = ROOT.RooArgList( w.set('POI') )
        for i in xrange( couplingsList.getSize() ):
            couplings.append( couplingsList[i] )

        yieldParameters = []
        yieldParameterList = ROOT.RooArgList( w.set('yieldParameters' if not self.newStyle else 'all_ggH_yieldParameters') )
        for i in xrange( yieldParameterList.getSize() ):
            yieldParameters.append( yieldParameterList[i] )

        parametrizations = []
        parametrizationList = ROOT.RooArgList( w.set('parametrizations') )
        for i in xrange( parametrizationList.getSize() ):
            parametrizations.append( parametrizationList[i] )

        self.w = w
        ROOT.SetOwnership( self.w, False )

        self.couplings = couplings
        for coupling in self.couplings:
            ROOT.SetOwnership( coupling, False )

        self.yieldParameters = yieldParameters
        for yieldParameter in self.yieldParameters:
            ROOT.SetOwnership( yieldParameter, False )

        self.parametrizations = parametrizations
        for parametrization in self.parametrizations:
            ROOT.SetOwnership( parametrization, False )

        wsFp.Close()


        if verbose:
            print '\n--------------'
            print 'Parametrized workspace \'{0}\''.format( wsFile )
            print '   couplings :           len {0:<3}, '.format(len(self.couplings)), [ i.GetName() for i in self.couplings ]
            print '   yieldParameters :     len {0:<3}, '.format(len(self.yieldParameters)), [ i.GetName() for i in self.yieldParameters ]
            print '   parametrizations :    len {0:<3}, '.format(len(self.parametrizations)), [ i.GetName() for i in self.parametrizations ]

            theoryBinBoundaries = Commands.ReadTheoryBinBoundariesFromWS( self.file )
            expBinBoundaries    = Commands.ReadExpBinBoundariesFromWS( self.file )
            print '   theoryBinBoundaries : len {0:<3}, '.format(len(theoryBinBoundaries)), theoryBinBoundaries
            print '   expBinBoundaries :    len {0:<3}, '.format(len(expBinBoundaries)), expBinBoundaries
            print


    def Evaluate( self, **kwargs ):

        if 'returnWhat' in kwargs:
            self.returnWhat = kwargs['returnWhat']
            del kwargs['returnWhat']

        if self.verbose:
            print 'Setting coupling values:'
            print '    '+ '\n    '.join([ '{0} = {1}'.format( key, value ) for key, value in kwargs.iteritems() ])

        for couplingName, couplingValue in kwargs.iteritems():

            for coupling in self.couplings:
                if couplingName == coupling.GetName():
                    coupling.setVal( couplingValue )
                    if self.verbose: print 'Set \'{0}\' to {1}'.format( coupling.GetName(), coupling.getVal() )
                    break
            else:
                print '[warning] There was no variable found for \'{0}\'; continuing'.format( couplingName )


        if self.returnWhat == 'theory':
            return [ parametrization.getVal() for parametrization in self.parametrizations ]

        elif self.returnWhat == 'exp':
            return [ yieldParameter.getVal() for yieldParameter in self.yieldParameters ]

        elif self.returnWhat == 'both':
            return [
                [ parametrization.getVal() for parametrization in self.parametrizations ],
                [ yieldParameter.getVal() for yieldParameter in self.yieldParameters ],
                ]


    def GetOutputContainer( self, **kwargs ):

        self.theoryBinBoundaries = Commands.ReadTheoryBinBoundariesFromWS( self.file )
        self.expBinBoundaries    = Commands.ReadExpBinBoundariesFromWS( self.file )


        mus = self.Evaluate( **kwargs )
        if self.returnWhat == 'theory': mus = mus[1:]

        container = OutputContainer()

        container.mus       = mus
        container.ratios    = mus
        container.binValues = mus

        if self.returnWhat == 'theory':
            container.binBoundaries = deepcopy( self.theoryBinBoundaries )
        elif self.returnWhat == 'exp':
            container.binBoundaries = deepcopy( self.expBinBoundaries )
            if 'xMax' in kwargs:
                container.binBoundaries[-1] = kwargs['xMax']
                del kwargs['xMax']
        elif self.returnWhat == 'both':
            container.binBoundaries = [ deepcopy(self.theoryBinBoundaries), deepcopy(self.expBinBoundaries) ]

        container.name = 'param'
        for coupling in self.couplings:
            container.name += '_{0}_{1:.2f}'.format( coupling.GetName(), coupling.getVal() )

        container.GetTGraph()

        return container



class Parametrization():

    THROW_EVALUATION_WARNINGS = True

    def __init__( self, verbose=True, SM=None ):
        self.verbose = verbose
        self.Persistence = []

        self.parametrizedByMatrixInversion = False
        self.parametrizedByFitting         = False

        self.hasSM = False
        if SM:
            self.SetSM(SM)

    def SetSM( self, SM ):
        self.SM = SM
        self.hasSM = True


    def Parametrize(
            self,
            containers,
            couplingsToParametrize = [ 'kappab', 'kappac' ],
            includeLinearTerms = True,
            ):
        self.parametrizedByMatrixInversion = True

        self.couplings = couplingsToParametrize

        couplingCombinations = []
        couplingCombinations += [ list(couplingTuple) for couplingTuple in itertools.combinations( self.couplings, 2 ) ]
        couplingCombinations += [ [ coupling, coupling ] for coupling in self.couplings ]
        if includeLinearTerms:
            couplingCombinations += [ [coupling] for coupling in self.couplings ] + [[]]
        self.couplingCombinations = couplingCombinations

        nComponents = len(couplingCombinations)
        self.nComponents = nComponents

        if len(containers) < nComponents:
            Commands.ThrowError( 'Need at least as much input ({0}) as desired components ({1})'.format( len(containers), nComponents ) )
            sys.exit()
        elif len(containers) > nComponents:
            print '[info] Need only {0} containers; limiting number of containers'.format( nComponents )
            print '[fixme] Now taking the first {0} (may lead to degenerate matrix...)'.format( nComponents )
            # containers = containers[:nComponents]
            random.seed(1001)
            containers = random.sample( containers, nComponents )


        couplingMatrix = []
        for container in containers:
            row = []
            for couplingCombination in couplingCombinations:
                product = 1.0
                for coupling in couplingCombination:
                    product *= getattr( container, coupling )
                row.append( product )
            couplingMatrix.append( row )

        couplingMatrix = numpy.array( couplingMatrix )
        if self.verbose:
            print '\nSquared coupling terms:'
            print '  ', couplingCombinations
            print 'Found the following coupling matrix:'
            print couplingMatrix
        couplingMatrixInv = numpy.linalg.inv( couplingMatrix )

        nBins = len(containers[0].binBoundaries) - 1
        self.nBins = nBins
        self.binBoundaries = containers[0].binBoundaries

        # ======================================
        # Create parametrization per bin

        # def make_evaluateParametrization( components ):
        #     if not len(self.couplingCombinations) == len(components):
        #         Commands.ThrowError('Number of coupling combinations ({0}) does not agree with number of components ({1})'.format( len(couplingCombinations), len(components) ) )
        #         sys.exit()

        #     def evaluateParametrization():
        #         res = 0.
        #         for couplingCombination, component in zip( self.couplingCombinations, components ):
        #             product = 1.0
        #             for couping in couplingCombination:
        #                 product *= getattr( self, coupling )
        #             res += product * component
        #         return res

        #     return 

        componentsPerBin = []
        for iBin in xrange(nBins):
            yValues = numpy.array([ [container.crosssection[iBin]] for container in containers ])
            components = list(itertools.chain.from_iterable( couplingMatrixInv.dot( yValues ) ))
            componentsPerBin.append( components )

            if self.verbose:
                print 'Components for bin {0}'.format(iBin)
                for couplingCombination, component in zip( couplingCombinations, components ):
                    print '    comp. {0:30} = {1}'.format( ', '.join(couplingCombination), component )

        if self.verbose: print 'Maximal variations per component over all bins:'
        for iComponent in xrange(nComponents):
            maxVar = max([ abs(components[iComponent]) for components in componentsPerBin ])
            couplingCombination = couplingCombinations[iComponent]
            if self.verbose: print '    comp. {0:30} = {1}'.format( ', '.join(couplingCombination), maxVar )

        self.componentsPerBin = componentsPerBin


        # ======================================
        # If a SM was set, calculate the SM cross section based on this parametrization

        if self.hasSM:
            self.SMXS_fromParametrization = deepcopy( self.Evaluate(**{
                'kappab' : 1, 'kappac' : 1, 'ct' : 1, 'cg' : 0, 'kappab' : 1,
                }) )


    def Evaluate( self, **kwargs ):

        self.verbose = kwargs.get( 'verbose', False )
        for key, value in kwargs.iteritems():
            if self.verbose:
                if self.THROW_EVALUATION_WARNINGS:print '[in Parametrization] Setting {0} to {1}'.format( key, value )
            setattr( self, key, value )


        if self.parametrizedByMatrixInversion:

            xsFromParametrization = []
            for iBin in xrange(self.nBins):
                components = self.componentsPerBin[iBin]

                res = 0.
                for couplingCombination, component in zip( self.couplingCombinations, components ):
                    product = 1.0
                    for coupling in couplingCombination:
                        product *= getattr( self, coupling )
                    res += product * component

                xsFromParametrization.append( res )

            self.THROW_EVALUATION_WARNINGS = False
            return xsFromParametrization


        elif self.ParametrizeByFitting:

            if self.fitWithScipy:

                xs = []
                for coefficients in self.coefficientsPerBin:
                    C = tuple([ getattr( self, coupling ) for coupling in self.couplings ])
                    xs.append(
                        self.parabolicParametrizationFunction( C, *coefficients )
                        )
                self.THROW_EVALUATION_WARNINGS = False
                return xs

            else:

                # Set RooRealVars to the right couplings
                for couplingRooRealVar in self.couplingRooRealVars:
                    couplingRooRealVar.setVal(
                        getattr( self, couplingRooRealVar.GetName() )
                        )

                xs = []
                for fitEval in self.fitEvals:
                    xs.append( fitEval.getVal() )

                self.THROW_EVALUATION_WARNINGS = False
                return xs


    def EvaluateForBinning( self, theory_binning, exp_binning, returnRatios=True, lastBinIsOverflow=False, **kwargs ):
        theory_xs = self.Evaluate(**kwargs)

        if len(theory_binning)-1 > len(theory_xs):
            if self.THROW_EVALUATION_WARNINGS: Commands.Warning( 'Supplied binning has {0} bins, whereas the parametrization only yields {1} values; assuming the last few bins can be thrown away'.format( len(theory_binning), len(theory_xs) ) )
            theory_binning = theory_binning[:len(theory_xs)+1]
        elif len(theory_binning)-1 < len(theory_xs):
            Commands.ThrowError( 'Supplied binning has {0} bins, whereas the parametrization yields {1} values; This cannot be integrated over'.format( len(theory_binning), len(theory_xs) ) )

        integral = TheoryCommands.GetIntegral( theory_binning, theory_xs )

        exp_xs = []
        for left, right in zip( exp_binning[:-1], exp_binning[1:] ):
            if lastBinIsOverflow and right == exp_binning[-1]:
                exp_xs.append( integral( left, right ) )
            else:
                exp_xs.append( integral( left, right ) / ( right - left ) )

        if returnRatios:
            if not self.hasSM:
                Commands.ThrowError( 'Use \'SetSM(container)\' method before calling \'Parametrize()\'', throwException=True )
            exp_xs_SM_integral = TheoryCommands.GetIntegral( self.SM.binBoundaries, self.SM.crosssection )
            exp_xs_SM = []
            for left, right in zip( exp_binning[:-1], exp_binning[1:] ):
                exp_xs_SM.append( exp_xs_SM_integral( left, right ) / ( right - left ) )

            ratios = [ xs / SMxs if SMxs != 0. else 0. for xs, SMxs in zip( exp_xs, exp_xs_SM ) ]

            self.THROW_EVALUATION_WARNINGS = False
            return ratios

        else:
            self.THROW_EVALUATION_WARNINGS = False
            return exp_xs


    def ParametrizeByFitting( self, *args, **kwargs ):
        self.fitWithScipy = False
        if 'fitWithScipy' in kwargs:
            self.fitWithScipy = kwargs['fitWithScipy']
            del kwargs['fitWithScipy']

        if self.fitWithScipy:
            self.ParametrizeByFitting_Scipy( *args, **kwargs )
        else:
            self.ParametrizeByFitting_RooFit( *args, **kwargs )

        # ======================================
        # If a SM was set, calculate the SM cross section based on this parametrization

        if self.hasSM:
            self.SMXS_fromParametrization = deepcopy( self.Evaluate(**{
                'kappab' : 1, 'kappac' : 1, 'ct' : 1, 'cg' : 0, 'kappab' : 1,
                }) )



    def ParametrizeByFitting_Scipy(
            self,
            containers,
            couplingsToParametrize = [ 'kappab', 'kappac' ],
            includeLinearTerms = True,
            ):

        print 'Trying to import scipy.optimize.curve_fit'
        from scipy.optimize import curve_fit

        self.parametrizedByFitting = True

        self.couplings = couplingsToParametrize
        self.nCouplings = len(self.couplings)

        couplingCombinations = []
        couplingCombinations += [ list(couplingTuple) for couplingTuple in itertools.combinations( self.couplings, 2 ) ]
        couplingCombinations += [ [ coupling, coupling ] for coupling in self.couplings ]
        if includeLinearTerms:
            couplingCombinations += [ [coupling] for coupling in self.couplings ] + [[]]
        self.couplingCombinations = couplingCombinations

        nComponents = len(couplingCombinations)
        self.nComponents = nComponents

        if len(containers) < nComponents:
            Commands.ThrowError( 'Need at least as much input ({0}) as desired components ({1})'.format( len(containers), nComponents ) )
            sys.exit()

        nBins = len(containers[0].binBoundaries) - 1
        self.nBins = nBins
        self.binBoundaries = containers[0].binBoundaries

        self.nContainers = len(containers)

        # ======================================
        # Fitting


        # Construct function to be minimized

        couplingNameToIndex = { coupling : iCoupling for iCoupling, coupling in enumerate(self.couplings) }
        couplingIndexToName = { iCoupling : coupling for iCoupling, coupling in enumerate(self.couplings) }

        def parabolicParametrization( C, *coefficients ):
            paramValue = 0.
            for couplingList, coefficient in zip( self.couplingCombinations, coefficients ):
                component = coefficient
                for coupling in couplingList:
                    component *= C[ couplingNameToIndex[coupling] ]
                paramValue += component
            return paramValue
        self.parabolicParametrizationFunction = parabolicParametrization


        self.coefficientsPerBin = []

        for iBin in xrange(self.nBins):
            # if iBin == 0: continue

            print '\n' + '='*50 + '\nFitting bin {0}'.format(iBin)

            # Construct input for least_sq: ( c1, c2 ) per input container
            couplingTuple = [ [] for coupling in self.couplings ]
            for iCoupling, coupling in enumerate(self.couplings):
                for container in containers:
                    couplingTuple[iCoupling].append( getattr(container, coupling) )
            couplingTuple = tuple(couplingTuple)

            # Construct input for least_sq: xs[iBin] per input container
            xsTuple = tuple([ container.crosssection[iBin] for container in containers ])

            # initial guess for coefficients
            p0 = tuple([ 1. for i in xrange(self.nComponents) ])


            print '\n    Overview of given input to minimizer:'
            print '    couplingTuple :', couplingTuple
            print '    xsTuple       :', xsTuple
            print '    p0            :', p0

            print '\n    Fitting...'
            fittedCoefficients, covMat = curve_fit(
                parabolicParametrization,
                couplingTuple,
                xsTuple,
                p0
                )
            fittedCoefficients = list(fittedCoefficients)
            print '\n    fittedCoefficients:', fittedCoefficients


            print '\n    Test evaluations'

            for i in xrange(len(xsTuple)):

                couplings = tuple([ couplingTuple[iCoupling][i] for iCoupling in xrange(len(self.couplings)) ])
                xs        = xsTuple[i]

                xsParametrization = parabolicParametrization( couplings, *fittedCoefficients )

                line = []
                for iCoupling, couplingValue in enumerate(couplings):
                    couplingName = couplingIndexToName[iCoupling]
                    line.append( '{0:7} = {1:+5.1f}'.format( couplingName, couplingValue ) )

                line.append( 'xs_file = {0:10.4f}'.format(xs) )
                line.append( 'xs_param = {0:10.4f}'.format(xsParametrization) )

                print '    ' + ' | '.join(line)

            self.coefficientsPerBin.append( fittedCoefficients )








    def ParametrizeByFitting_RooFit(
            self,
            containers,
            couplingsToParametrize = [ 'kappab', 'kappac' ],
            includeLinearTerms = True,
            ):
        self.parametrizedByFitting = True

        self.couplings = couplingsToParametrize

        couplingCombinations = []
        couplingCombinations += [ list(couplingTuple) for couplingTuple in itertools.combinations( self.couplings, 2 ) ]
        couplingCombinations += [ [ coupling, coupling ] for coupling in self.couplings ]
        if includeLinearTerms:
            couplingCombinations += [ [coupling] for coupling in self.couplings ] + [[]]
        self.couplingCombinations = couplingCombinations

        nComponents = len(couplingCombinations)
        self.nComponents = nComponents

        if len(containers) < nComponents:
            Commands.ThrowError( 'Need at least as much input ({0}) as desired components ({1})'.format( len(containers), nComponents ) )
            sys.exit()

        nBins = len(containers[0].binBoundaries) - 1
        self.nBins = nBins
        self.binBoundaries = containers[0].binBoundaries

        self.nContainers = len(containers)


        # ======================================
        # Start fit

        # Open up RooRealVars per coupling
        couplingRooRealVars = []
        for coupling in self.couplings:
            couplingRooRealVar = ROOT.RooRealVar(
                coupling, coupling, 1., -50., 50.
                )
            ROOT.SetOwnership( couplingRooRealVar, False )
            couplingRooRealVars.append( couplingRooRealVar )
        self.couplingRooRealVars = couplingRooRealVars

        self.fitEvals = []
        for iBin in xrange(self.nBins):

            # Open up RooRealVars per component
            componentRooRealVars = []
            for couplingCombination in self.couplingCombinations:

                if len(couplingCombination) == 0:
                    name = 'bin{0}_component_CONSTANT'.format( iBin )
                else:
                    name = 'bin{0}_component_{1}'.format( iBin, '_'.join(couplingCombination) )

                # comp. kappab, kappac                 = 0.001285759
                # comp. kappab, kappab                 = 0.006308878
                # comp. kappac, kappac                 = 6.67915333333e-05
                # comp. kappab                         = 0.116259917
                # comp. kappac                         = 0.0142519758
                # comp.                                = 1.1356275955

                componentRooRealVar = ROOT.RooRealVar(
                    name, name, 0.0001
                    )
                ROOT.SetOwnership( componentRooRealVar, False )
                componentRooRealVar.setConstant(False)

                if len(couplingCombination) == 0:
                    componentRooRealVar.setRange( -3., 3. )
                elif len(couplingCombination) == 1:
                    componentRooRealVar.setRange( -.5, .5 )
                elif len(couplingCombination) == 2:
                    componentRooRealVar.setRange( -.01, .01 )
                else:
                    componentRooRealVar.setRange( -3., 3. )

                # Append the list of couplings as a python list
                componentRooRealVar.couplings = couplingCombination

                componentRooRealVars.append( componentRooRealVar )


            # Dict to get a consistent number back per component
            argumentIndices = {}
            fitArgList = ROOT.RooArgList()
            ROOT.SetOwnership( fitArgList, False )
            for var in componentRooRealVars:
                argumentIndices[ var.GetName() ] = '@{0}'.format( len(argumentIndices) )
                fitArgList.add(var)
            self.Persistence.append( fitArgList )

            # Construct chi2 components
            chi2Components = []
            chi2ComponentsArgList = ROOT.RooArgList()
            ROOT.SetOwnership( chi2ComponentsArgList, False )
            for iContainer in xrange(self.nContainers):
                container = containers[iContainer]

                if self.verbose:
                    print '\n' + '-'*50
                    print 'Processing container {0}: {1}'.format(
                        iContainer,
                        ', '.join([ '{0} = {1}'.format( coupling, getattr( container, coupling ) ) for coupling in self.couplings ])
                        )

                # ======================================
                # Build the fit function

                fitName = 'fit_bin{0}_container{1}'.format( iBin, iContainer )

                fitFunction = []
                for componentRooRealVar in componentRooRealVars:
                    functionPart = argumentIndices[componentRooRealVar.GetName()]
                    for coupling in componentRooRealVar.couplings:
                        functionPart += '*{0}'.format( getattr( container, coupling ) )
                    fitFunction.append( functionPart )
                fitFunction = '+'.join(fitFunction)

                fit = ROOT.RooFormulaVar( fitName, fitName, fitFunction, fitArgList )
                ROOT.SetOwnership( fit, False )

                if self.verbose:
                    print '\nDoing \'fit.Print():\''
                    fit.Print()

                chi2ComponentName = 'chi2Component_bin{0}_container{1}'.format( iBin, iContainer )
                chi2ComponentFunction = '(@0-{0})**2'.format( container.crosssection[iBin] )
                chi2ComponentArgList = ROOT.RooArgList( fit )
                ROOT.SetOwnership( chi2ComponentArgList, False )
                chi2Component = ROOT.RooFormulaVar( chi2ComponentName, chi2ComponentName, chi2ComponentFunction, chi2ComponentArgList )
                ROOT.SetOwnership( chi2Component, False )

                if self.verbose:
                    print '\nDoing \'chi2Component.Print():\''
                    chi2Component.Print()

                chi2Components.append( chi2Component )
                chi2ComponentsArgList.add( chi2Component )


            chi2Name = 'chi2_bin{0}'.format(iBin)
            chi2 = ROOT.RooAddition( chi2Name, chi2Name, chi2ComponentsArgList )
            ROOT.SetOwnership( chi2, False )

            if self.verbose:
                print '\nDoing \'chi2.Print():\''
                chi2.Print()


            # ======================================
            # Minimize

            print
            Minimizer = ROOT.RooMinimizer( chi2 )
            if self.verbose:
                print '\nSuccesfully initialized the Minimizer, now trying MIGRAD\n'
            Minimizer.migrad()


            # ======================================
            # Post processing

            if self.verbose:

                print '\nMinimization done; best chi2 = {0}'.format( chi2.getVal() )

                print '\nBest fit of components in bin {0}:'.format(iBin)
                for componentRooRealVar in componentRooRealVars:
                    print '    {0} = {1}'.format(
                        componentRooRealVar.GetName(),
                        componentRooRealVar.getVal()
                        )


            # Create a fit function that can be evaluated for a random coupling

            # Dict to get a consistent number back per component
            argumentIndices = {}
            fitEvalArgList = ROOT.RooArgList()
            ROOT.SetOwnership( fitEvalArgList, False )

            allVars = couplingRooRealVars + componentRooRealVars
            for var in couplingRooRealVars + componentRooRealVars:
                argumentIndices[ var.GetName() ] = '@{0}'.format( len(argumentIndices) )
                fitEvalArgList.add(var)

            fitEvalName = 'fit_bin{0}'.format(iBin)

            fitEvalFunction = []
            for componentRooRealVar in componentRooRealVars:
                functionPart = argumentIndices[componentRooRealVar.GetName()]
                for coupling in componentRooRealVar.couplings:
                    functionPart += '*{0}'.format( argumentIndices[coupling] )
                fitEvalFunction.append( functionPart )
            fitEvalFunction = '+'.join(fitEvalFunction)

            fitEval = ROOT.RooFormulaVar( fitEvalName, fitEvalName, fitEvalFunction, fitEvalArgList )
            ROOT.SetOwnership( fitEval, False )

            if self.verbose:

                print '\nDoing \'fitEval.Print():\''
                fitEval.Print()

                for container in containers:

                    for couplingRooRealVar in couplingRooRealVars:
                        couplingRooRealVar.setVal(
                            getattr( container, couplingRooRealVar.GetName() )
                            )

                    line = []
                    for coupling in self.couplings:
                        line.append( '{0} = {1:5.1f}'.format( coupling, getattr( container, coupling ) ) )
                    line.append( 'exp. xs = {0:8.3f}'.format( container.crosssection[iBin] ) )
                    line.append( 'param. xs = {0:8.3f}'.format( fitEval.getVal() ) )

                    if self.hasSM and container == self.SM:
                        print '*   ' + ' | '.join(line)
                    else:
                        print '    ' + ' | '.join(line)

            self.fitEvals.append( fitEval )



    def GetOutputContainer( self, **kwargs ):

        # self.theoryBinBoundaries = Commands.ReadTheoryBinBoundariesFromWS( self.file )
        # self.expBinBoundaries    = Commands.ReadExpBinBoundariesFromWS( self.file )

        crosssection = self.Evaluate( **kwargs )


        container = OutputContainer()

        container.crosssection = crosssection
        container.binBoundaries = self.binBoundaries

        if self.hasSM:
            mus = [ xs / SMxs if not SMxs == 0. else 0. for xs, SMxs in zip( crosssection, self.SMXS_fromParametrization ) ]
            binWidths = [ 0.5*(container.binBoundaries[i]+container.binBoundaries[i+1]) for i in xrange(len(container.binBoundaries)-1) ]
            container.integral = [ width * xs for width, crosssection in zip( binWidths, crosssection ) ]
        else:
            mus = crosssection

        container.mus       = mus
        container.ratios    = mus
        container.binValues = mus

        container.name = 'param'
        for coupling in self.couplings:
            setattr( container, coupling, getattr( self, coupling ) )
            container.name += '_{0}_{1:.2f}'.format( coupling, getattr( self, coupling ) )

        container.GetTGraph()

        return container
