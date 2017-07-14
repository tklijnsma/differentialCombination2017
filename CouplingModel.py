
from HiggsAnalysis.CombinedLimit.PhysicsModel import *

import os, sys, numpy, itertools, re
import ROOT

class CouplingModel( PhysicsModel ):
    ''' Model used to unfold differential distributions '''

    def __init__(self):
        PhysicsModel.__init__(self)

        self.mHRange=[]
        self.mass = 125.
        self.verbose = 1

        self.theories = []

        self.SMbinBoundaries = []
        self.SMXS = []


    def chapter( self, text, indent=0 ):
        if self.verbose:
            print '\n{tabs}{line}\n{tabs}{text}\n'.format(
                tabs = '    ' * indent,
                line = '-' * 70,
                text = text
                )

    def readTheoryFile( self, theoryFile ):
        with open( theoryFile, 'r' ) as theoryFp:
            lines = [ l.strip() for l in theoryFp.readlines() if len(l.strip()) > 0 and not l.strip().startswith('#') ]

        couplings    = {}
        ratios       = []
        crosssection = []

        for line in lines:
            key, value = line.split('=',1)
            if key == 'binBoundaries':
                binBoundaries = [ float(v) for v in value.split(',') ]
            elif key == 'crosssection':
                crosssection = [ float(v) for v in value.split(',') ]
            elif key == 'ratios':
                ratios = [ float(v) for v in value.split(',') ]
            else:
                couplings[key] = float(value)

        return binBoundaries, ratios, crosssection



    def setPhysicsOptions( self, physOptions ):
        self.chapter( 'Starting model.setPhysicsOptions()' )

        for physOption in physOptions:
            optionName, optionValue = physOption.split('=',1)

            if False:
                pass

            elif optionName == 'verbose':
                self.verbose = int(optionValue)

            elif optionName == 'theory':
                # Syntax: --PO theory=[ct=1,cg=1,file=some/path/.../] , brackets optional

                # Delete enclosing brackets if present
                if optionValue.startswith('['): optionValue = optionValue[1:]
                if optionValue.endswith(']'): optionValue = optionValue[:-1]

                theoryDict = {                
                    'couplings'      : {},
                    'binBoundaries'  : [],
                    'crosssection'   : [],
                    'ratios'         : [],
                    }
                for component in optionValue.split(','):
                    couplingName, couplingValue = component.split('=',1)
                    # Only exception is file, which is the file to read the shape from
                    if couplingName == 'file':
                        theoryDict['binBoundaries'], theoryDict['ratios'], theoryDict['crosssection'] = self.readTheoryFile( couplingValue )
                        continue
                    theoryDict['couplings'][couplingName] = float(couplingValue)

                self.theories.append( theoryDict )
                
            elif optionName == 'SM':

                # Delete enclosing brackets if present
                if optionValue.startswith('['): optionValue = optionValue[1:]
                if optionValue.endswith(']'): optionValue = optionValue[:-1]

                self.SMDict = {                
                    'couplings'      : {},
                    'binBoundaries'  : [],
                    'crosssection'   : [],
                    }
                for component in optionValue.split(','):
                    couplingName, couplingValue = component.split('=',1)
                    # Only exception is file, which is the file to read the shape from
                    if couplingName == 'file':
                        self.SMDict['binBoundaries'], dummy, self.SMDict['crosssection'] = self.readTheoryFile( couplingValue )
                        continue
                    self.SMDict['couplings'][couplingName] = float(couplingValue)

                # SM cross section extended with underflow and overflow bin
                self.SMXS = [ self.SMDict['crosssection'][0] ] + self.SMDict['crosssection'] + [ self.SMDict['crosssection'][-1] ]


            elif optionName == 'higgsMassRange':
                self.mHRange = [ float(v) for v in optionValue.split(',') ]


            else:
                print 'Unknown option \'{0}\'; skipping'.format(optionName)
                continue



    def doParametersOfInterest(self):
        self.chapter( 'Starting model.doParametersOfInterest()' )

        bins    = self.DC.bins
        signals = self.DC.signals
        signals = [ s for s in signals if not s == 'OutsideAcceptance' ]
        signals.sort( key = lambda n: float(
            n.split('_')[2].replace('p','.').replace('m','-').replace('GT','').replace('GE','').replace('LT','').replace('LE','') ) )

        if self.modelBuilder.out.var("MH"):
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                print 'MH will be assumed to be', self.mass
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.mass)
        else:
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
            else:
                print 'MH (not there before) will be assumed to be', self.mass
                self.modelBuilder.doVar("MH[%g]" % self.mass)


        # ======================================
        # Figure out experimental binning from signals

        productionMode, observableName = signals[0].split('_')[:2]

        expBinBoundaries = []
        for signal in signals:
            if signal == 'OutsideAcceptance': continue
            bounds = signal.split('_')[2:]
            for bound in bounds:
                bound = bound.replace('p','.').replace('m','-').replace('GT','').replace('GE','').replace('LT','').replace('LE','')
                expBinBoundaries.append( float(bound) )
        expBinBoundaries = list(set(expBinBoundaries))
        expBinBoundaries.sort()
        expBinBoundaries.append( 500. )  # Overflow bin
        nExpBins = len(expBinBoundaries) - 1

        print 'Determined the following experimental bin boundaries:', expBinBoundaries
        print '   from the following signals:'
        print '   ' + '\n   '.join( signals )


        # ======================================
        # Handling theory

        # Sort according to ct value
        try:
            self.theories.sort( key=lambda t: t['ct'] )
        except:
            pass

        # Make sure list lengths do not differ between theories
        #   (Could also check whether binBoundaries matches exactly)
        def checkListLengths( key ):
            elements = list(set([ len(theory[key]) for theory in self.theories ]))
            if len(elements) > 1:
                print 'WARNING: Found different list lengths for key \'{0}\''.format(key)
        checkListLengths( 'couplings' )
        checkListLengths( 'binBoundaries' )
        checkListLengths( 'crosssection' )


        nCouplings                = len(self.theories[0]['couplings'])
        couplings                 = self.theories[0]['couplings'].keys()
        nTheoryBins               = len(self.theories[0]['binBoundaries']) - 1
        theoryBinBoundaries       = self.theories[0]['binBoundaries']

        # Get list of squared terms and unique combinations (works as expected also for 1 couplings cases)
        couplingCombinations = []
        couplingCombinations.extend( [ ( coupling, coupling ) for coupling in couplings ] )
        couplingCombinations.extend( list(itertools.combinations( couplings, 2 )) )

        # Number of theories necessary to perform a parametrization
        nParameters = sum(range(nCouplings+1))
        if not len(couplingCombinations) == nParameters:
            print 'ERROR: amount of coupling combinations should be equal to number of expected parameters'
            sys.exit()

        if len(self.theories) > nParameters:
            print '\n{0} theories supplied but only {1} needed for parametrization; taking only the first {1} theories:'.format(
                len(self.theories), nParameters )
            self.theories = self.theories[:nParameters]
            for theory in self.theories:
                print '    ' + ', '.join([ '{0:10} = {1:10}'.format( cName, cValue ) for cName, cValue in theory['couplings'].iteritems() ])
        elif len(self.theories) < nParameters:
            print 'ERROR: cannot parametrize because number of supplied theories is too small (need {0} theories, found {1})'.format(
                nParameters, len(self.theories) )
            sys.exit()


        # ======================================
        # Set variables for all couplings

        for coupling in couplings:
            self.modelBuilder.doVar(
                '{coupling}[{default},{down},{up}]'.format(
                    coupling = coupling,
                    default  = self.SMDict['couplings'][coupling],
                    down     = self.SMDict['couplings'][coupling] - 2.0,
                    up       = self.SMDict['couplings'][coupling] + 2.0,
                    )
                )


        # ======================================
        # Find the parametrizations for all nTheorybins, plus 1 underflow and 1 overflow parametrization

        # Calculate the coupling matrix
        couplingMat = []
        for theory in self.theories:
            couplingMat.append(
                [ theory['couplings'][c1]*theory['couplings'][c2] for c1, c2 in couplingCombinations ]
                )
        couplingMat = numpy.array( couplingMat )
        couplingMatInv = numpy.linalg.inv( couplingMat )

        if self.verbose:
            print '\nSquared coupling terms:'
            print '  ', couplingCombinations
            print 'Found the following coupling matrix:'
            print couplingMat


        parametrizations = []

        # Underflow (left extrapolation)
        parametrizations.append( [ 0. for i in xrange(nParameters) ] )

        # Add parametrization for each bin
        for iTheoryBin in xrange(nTheoryBins):
            crosssections = numpy.array([ [theory['ratios'][iTheoryBin]] for theory in self.theories ])
            parametrization = list(itertools.chain.from_iterable( couplingMatInv.dot( crosssections ) ))
            parametrizations.append( parametrization )

        # Overflow (right extrapolation)
        # parametrizations.append( [ 0. for i in xrange(nParameters) ] )
        parametrizations.append( parametrizations[-1][:] )

        if self.verbose:
            print '\nParametrizations per theory bin:'
            for i, parametrization in enumerate(parametrizations):
                print '{0:4}: {1}'.format( i, parametrization )
            print


        # Assign coupling index for RooFormulas
        iCoupling = 0
        argCoupling = {}
        for coupling in couplings:
            argCoupling[coupling] = '@{0}'.format(iCoupling)
            iCoupling += 1
        # argString = ','.join([ argCoupling[coupling] for coupling in couplings ])

        parametrizationNames = []
        for iParametrization, parametrization in enumerate( parametrizations ):

            # Also import it to the workspace
            parametrizationName = 'parametrization{0}'.format( iParametrization )
            parametrizationNames.append( parametrizationName )

            parametrizationString = ''
            for parameter, (coupling1, coupling2) in zip( parametrization, couplingCombinations ):
                parametrizationString += '+{0}*{1}*{2}'.format(
                    parameter, argCoupling[coupling1], argCoupling[coupling2]
                    )
            parametrizationString = parametrizationString[1:]

            parametrizationExpression = (
                'expr::{name}('
                '"({string})", '
                '{arguments} )'
                ).format(
                    name      = parametrizationName,
                    string    = parametrizationString,
                    arguments = ','.join([ coupling for coupling in couplings ]),
                    )

            if self.verbose > 1:
                print 'Importing parametrization {0}:'.format( iParametrization )
                print parametrizationExpression

            self.modelBuilder.factory_( parametrizationExpression )

        self.modelBuilder.out.defineSet( 'parametrizations', ','.join(parametrizationNames) )



        # ======================================
        # Find the theory bin widths and parametrization indices inside the experimental bins

        for iExpBin in xrange(nExpBins):

            expBoundLeft  = expBinBoundaries[iExpBin]
            expBoundRight = expBinBoundaries[iExpBin+1]
            expBinStr     = signals[iExpBin]

            theoryBinBoundariesInsideExperimentalBin = []
            parametrizationIndices                   = []

            # Extrapolation on left of theory spectrum
            if expBoundLeft < theoryBinBoundaries[0]:
                theoryBinBoundariesInsideExperimentalBin.append( expBoundLeft )
                parametrizationIndices.append( 0 )

            for iTheoryBin in xrange(nTheoryBins):

                # Left-most bin
                if (
                    expBoundLeft >= theoryBinBoundaries[iTheoryBin] and
                    expBoundLeft < theoryBinBoundaries[iTheoryBin+1]
                    ):
                    theoryBinBoundariesInsideExperimentalBin.append( expBoundLeft )
                    parametrizationIndices.append( iTheoryBin+1 )

                # Bins in between; only stricly 'between' (and not 'on') the left and right bounds
                if (
                    expBoundLeft < theoryBinBoundaries[iTheoryBin] and
                    theoryBinBoundaries[iTheoryBin] < expBoundRight
                    ):
                    theoryBinBoundariesInsideExperimentalBin.append( theoryBinBoundaries[iTheoryBin] )
                    parametrizationIndices.append( iTheoryBin+1 )

                # Right-most bin
                if (
                    expBoundRight > theoryBinBoundaries[iTheoryBin] and
                    expBoundRight <= theoryBinBoundaries[iTheoryBin+1]
                    ):
                    theoryBinBoundariesInsideExperimentalBin.append( expBoundRight )
                    break

            # Extrapolation on right of theory spectrum
            if expBoundRight > theoryBinBoundaries[-1]:
                theoryBinBoundariesInsideExperimentalBin.append( expBoundRight )
                parametrizationIndices.append( nTheoryBins+1 )


            nTheoryBinsInsideExperimentalBin = len(theoryBinBoundariesInsideExperimentalBin)-1
            theoryBinWidthsInsideExperimentalBin = [
                theoryBinBoundariesInsideExperimentalBin[i+1] - theoryBinBoundariesInsideExperimentalBin[i] for i in xrange(nTheoryBinsInsideExperimentalBin)
                ]


            if self.verbose:
                print '\nProcessing experimental bin {0}'.format(expBinStr)
                print 'Found the following theory bin boundaries inside experimental bin {0} to {1}:'.format(
                    expBoundLeft, expBoundRight )
                print '  ', theoryBinBoundariesInsideExperimentalBin
                print '   The following parametrizations will be used to evaluate the cross section:'
                print '  ', parametrizationIndices


            # Build the analytic function that described the XS in the experimental bin
            #   as a function of the supplied couplings

            # sumweights = sum(theoryBinWidthsInsideExperimentalBin)
            # weights = [ w/sumweights for w in theoryBinWidthsInsideExperimentalBin ]
            weights = theoryBinWidthsInsideExperimentalBin[:]
            # sumweights = sum(weights)
            # weights = [ w / sumweights for w in weights ]

            if self.verbose:
                print '   Which will carry the following weights:'
                print '   [ ' + ', '.join([ '{0:.4f}'.format(w) for w in weights ]) + ']'

            weightedParameters = []
            for iParameter in xrange(nParameters):
                weightedParameter = 0.
                for iParametrization, weight in zip( parametrizationIndices, weights ):
                    weightedParameter += weight * parametrizations[iParametrization][iParameter] * self.SMXS[iParametrization]

                    if self.verbose > 1:
                        print '\n    Adding contribution from parametrization {0} to component {1}'.format( iParametrization, iParameter )
                        print '      weight                                         = ', weight
                        print '      parametrizations[iParametrization][iParameter] = ', parametrizations[iParametrization][iParameter]
                        print '      self.SMXS[iParametrization]                    = ', self.SMXS[iParametrization]
                        print '      contribution                                   = ', weight * parametrizations[iParametrization][iParameter] * self.SMXS[iParametrization]

                weightedParameters.append( weightedParameter )

            # Divide by SM cross section so we get a ratio back
            # SMXSInsideExperimentalBin = sum([
            #     self.SMXS[iParametrization] for iParametrization in parametrizationIndices
            #     ])
            SMXSInsideExperimentalBin = 0.
            for iParametrization, binWidth in zip( parametrizationIndices, theoryBinWidthsInsideExperimentalBin ):
                SMXSInsideExperimentalBin += self.SMXS[iParametrization] * binWidth

                if self.verbose > 1:
                    print '\n    Adding contribution of theory bin {0} to SM cross section in exp bin {1}'.format( iParametrization, expBinStr )
                    print '      self.SMXS[iParametrization] = ', self.SMXS[iParametrization]
                    print '      binWidth                    = ', binWidth
                    print '      contribution                = ', self.SMXS[iParametrization] * binWidth


            weightedParameters = [ p / SMXSInsideExperimentalBin for p in weightedParameters ]

            if self.verbose:
                print '   Weighted parameters:'
                print '  ', weightedParameters


            # Importing into WS is somewhat delicate
            weightedParameterNames = []
            productNames = []
            for iWeightedParameter, weightedParameter in enumerate(weightedParameters):

                coupling1, coupling2 = couplingCombinations[iWeightedParameter]
                componentName = coupling1 + coupling2


                weightedParameterName = 'weight_{signal}_component_{componentName}'.format(
                    signal = expBinStr,
                    componentName = componentName,
                    )

                self.modelBuilder.doVar(
                    '{0}[{1}]'.format(
                        weightedParameterName,
                        weightedParameter,
                        )
                    )
                weightedParameterNames.append( weightedParameterName )

                productName = 'r_{signal}_component_{componentName}'.format(
                    signal     = expBinStr,
                    componentName = componentName,
                    )
                productExpression = (
                    'expr::{productName}('
                    '"(@0*@1*@2)", '
                    '{weight},{coupling1},{coupling2} )'
                    ).format(
                        productName = productName,
                        weight      = weightedParameterName,
                        coupling1   = coupling1,
                        coupling2   = coupling2,
                        )
                if self.verbose: print productExpression
                self.modelBuilder.factory_( productExpression )
                productNames.append( productName )


            # Compile expression string for final yieldParameter
            yieldParameterExpression = 'expr::r_{signal}( "({formulaString})", {commaSeparatedParameters} )'.format(
                signal                   = expBinStr,
                formulaString            = '+'.join([ '@' + str(iComponent) for iComponent in xrange(nParameters) ]),
                commaSeparatedParameters = ','.join( productNames )
                )
            if self.verbose: print yieldParameterExpression
            self.modelBuilder.factory_( yieldParameterExpression )



            # ======================================
            # Much neater, but can't figure out the import without making the WS disfunctional

            # # Build yieldParameter expression
            # yieldParameterComponents = ROOT.RooArgList()
            # iComponent = 0
            # for weightedParameter, couplingCombination in zip( weightedParameters, couplingCombinations ):
            #     w =  ROOT.RooFit.RooConst( weightedParameter )
            #     ROOT.SetOwnership( w, False )
            #     c1 = self.modelBuilder.out.var( couplingCombination[0] )
            #     c2 = self.modelBuilder.out.var( couplingCombination[1] )
            #     product = ROOT.RooProduct(
            #             'mu_expbin{0}_component{1}'.format( iExpBin, iComponent ),
            #             'mu_expbin{0}_component{1}'.format( iExpBin, iComponent ),
            #             ROOT.RooArgList( w, c1, c2 )
            #             )
            #     ROOT.SetOwnership( product, False )
            #     yieldParameterComponents.add( product )
            #     iComponent += 1

            #     if self.verbose:
            #         print 'Created component {0}; doing product.Print():'.format( iComponent )
            #         product.Print()


            # yieldParameter = ROOT.RooAddition(
            #     'r_' + expBinStr,
            #     'r_' + expBinStr,
            #     yieldParameterComponents,
            #     )
            # ROOT.SetOwnership( yieldParameter, False )

            # if self.verbose:
            #     print 'Created yieldParameter; doing yieldParameter.Print():'
            #     yieldParameter.Print()

            # Only import todo - Neither work
            # getattr( self.modelBuilder.out, 'import' )( yieldParameter )
            # self.modelBuilder.doVar( yieldParameter.GetName() )


        self.modelBuilder.out.defineSet( 'POI', ','.join(couplings) )

        # Define 2 extra sets for plotting convenience
        self.modelBuilder.out.defineSet( 'couplings', ','.join(couplings) )
        self.modelBuilder.out.defineSet( 'yieldParameters', ','.join([ 'r_' + s for s in signals ]) )

        self.chapter( 'Starting model.getYieldScale()' )




    def getYieldScale( self, bin, process ):

        def p( yieldScale ):
            if self.verbose > 0:
                print 'Getting yield scale: bin: {0:30} | process: {1:16} | yieldScale: {2}'.format( bin, process, yieldScale )

        if not self.DC.isSignal[process]:
            p( 1 )
            return 1
        else:
            if process == 'OutsideAcceptance':
                p( 1 )
                return 1
            else:
                yieldParameter = 'r_' + process
                p( yieldParameter )
                return yieldParameter




couplingModel=CouplingModel()

