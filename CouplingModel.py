
from HiggsAnalysis.CombinedLimit.PhysicsModel import *

import os, sys, numpy, itertools, re
import ROOT
from math import sqrt
from copy import deepcopy
from array import array


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

        self.includeLinearTerms = False
        self.splitggH           = False
        self.MakeLumiScalable   = False

        self.theoryUncertaintiesPassed = False
        self.correlationMatrixPassed   = False
        self.covarianceMatrixPassed    = False

        self.manualExpBinBoundaries = []
        self.skipBins = []
        


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
            key = key.strip()
            if key == 'binBoundaries':
                binBoundaries = [ float(v) for v in value.split(',') ]
            elif key == 'crosssection':
                crosssection = [ float(v) for v in value.split(',') ]
            elif key == 'ratios':
                ratios = [ float(v) for v in value.split(',') ]
            else:
                try:
                    couplings[key] = float(value)
                except ValueError:
                    # print '[error]: Could not call float(value) for key/value pair:'
                    # print '    key   : {0}'.format( key )
                    # print '    value :'
                    # print value
                    # sys.exit()
                    continue

        return binBoundaries, ratios, crosssection


    def readCorrelationMatrixFile( self, correlationMatrixFile ):
        with open( correlationMatrixFile, 'r' ) as correlationMatrixFp:
            lines = [ l.strip() for l in correlationMatrixFp.readlines() if len(l.strip()) > 0 and not l.strip().startswith('#') ]

        corrMat = [ [ float(number) for number in line.split() ] for line in lines ]

        # Check if it is square
        if not all([ len(row) == len(corrMat) for row in corrMat ]):
            print '[ERROR] inputted matrix is not square - Found:'
            print corrMat

        return corrMat


    def readErrorFile( self, errorFile ):
        with open( errorFile, 'r' ) as errorFp:
            lines = [ l.strip() for l in errorFp.readlines() if len(l.strip()) > 0 and not l.strip().startswith('#') ]

        symmErrors = []
        for line in lines:
            line = line.split()

            if len(line) == 1:
                symmErrors.append( abs(float(line[0])) )
            elif len(line) == 2:
                symmErrors.append( 0.5*( abs(float(line[0])) + abs(float(line[1])) ) )
            else:
                print '[ERROR] Found {0} elements on line in \'{1}\''.format( len(line), errorFile )
                return
            
        return symmErrors

    # Literal example from Pasquale:

    # def Decorrelate( correlationMatrix, errors ):

    #     ndim = 2
    #     covariance = ROOT.TMatrixDSym(ndim)
    #     errors = [ 0.7, 1 ]
    #     correlations = [ [1, 0.5], [0.5, 1] ]
    #     transform = [ [1, 0], [0, 1] ]

    #     for ix in xrange(ndim):
    #         for iy in xrange(ndim):
    #             covariance[ix][iy] = correlations[ix][iy]*errors[ix]*errors[iy]
    #     covariance.Print()

    #     eigen = ROOT.TMatrixDSymEigen(covariance)
    #     vectors = eigen.GetEigenVectors()
    #     values  = eigen.GetEigenValues()

    #     vectors.Print()
    #     values.Print()


    #     for ix in xrange(ndim):
    #         for iy in xrange(ndim):
    #             transform[ix][iy] = vectors(ix,iy) * sqrt(values(iy))

    #     print transform

    def Decorrelate( self, covarianceMatrix ):

        N = len(covarianceMatrix)

        covarianceTMatrix = ROOT.TMatrixDSym(N)
        for i in xrange(N):
            for j in xrange(N):
                # covarianceTMatrix[i][j] = covarianceMatrix[i][j]
                covarianceTMatrix[i,j] = covarianceMatrix[i][j]

        eigenObject  = ROOT.TMatrixDSymEigen(covarianceTMatrix)
        eigenVectors = eigenObject.GetEigenVectors()
        eigenValues  = eigenObject.GetEigenValues()

        decorrelatedMatrix = [ [ 999 for j in xrange(N) ] for i in xrange(N) ]
        for i in xrange(N):
            for j in xrange(N):
                decorrelatedMatrix[i][j] = eigenVectors(i,j) * sqrt( eigenValues(j) )

        return decorrelatedMatrix


    def setPhysicsOptions( self, physOptions ):
        self.chapter( 'Starting model.setPhysicsOptions()' )

        for physOption in physOptions:
            optionName, optionValue = physOption.split('=',1)

            if False:
                pass

            elif optionName == 'lumiScale':
                self.MakeLumiScalable = eval(optionValue)

            elif optionName == 'splitggH':
                self.splitggH = eval(optionValue)

            elif optionName == 'binBoundaries':
                self.manualExpBinBoundaries = [ float(i) for i in optionValue.split(',') ]

            elif optionName == 'verbose':
                self.verbose = int(optionValue)

            elif optionName == 'skipBins':
                self.skipBins = optionValue.split(',')

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


            elif optionName == 'linearTerms':
                self.includeLinearTerms = eval(optionValue)


            # ======================================
            # Options to pass the theory uncertainties

            elif optionName == 'theoryUncertainties':
                self.theoryUncertainties = self.readErrorFile( optionValue )
                self.theoryUncertaintiesPassed = True

            elif optionName == 'correlationMatrix':
                self.correlationMatrix = self.readCorrelationMatrixFile( optionValue )
                self.correlationMatrixPassed = True

            elif optionName == 'covarianceMatrix':
                self.covarianceMatrix = self.readCorrelationMatrixFile( optionValue )
                self.covarianceMatrixPassed = True

            else:
                print 'Unknown option \'{0}\'; skipping'.format(optionName)
                continue



    def doParametersOfInterest(self):
        self.chapter( 'Starting model.doParametersOfInterest()' )

        bins    = self.DC.bins
        signals = self.DC.signals
        signals = [ s for s in signals if not 'OutsideAcceptance' in s ]
        
        if self.splitggH:
            signals = [ s for s in signals if 'ggH' in s ]

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

        if len(self.manualExpBinBoundaries) > 0:
            expBinBoundaries = self.manualExpBinBoundaries
            expBinBoundaries.append( 500. )
            nExpBins = len(expBinBoundaries) - 1
            print 'Taking manually specified bin boundaries:', expBinBoundaries

        else:
            expBinBoundaries = []
            for signal in signals:
                if 'OutsideAcceptance' in signal: continue
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

        # Also attach to class so they are accesible in other methods
        self.expBinBoundaries = expBinBoundaries
        self.nExpBins         = nExpBins


        ########################################
        # Handling theory
        ########################################

        # Make sure list lengths do not differ between theories
        #   (Could also check whether binBoundaries matches exactly)
        for key in [ 'couplings', 'binBoundaries', 'crosssection' ]:
            lengthsOfLists = list(set([ len(theory[key]) for theory in self.theories ]))
            if len(lengthsOfLists) > 1:
                print 'WARNING: Found different list lengths for theory attribute \'{0}\': {1}'.format( key, lengthsOfLists )

        nCouplings                = len(self.theories[0]['couplings'])
        couplings                 = self.theories[0]['couplings'].keys()
        nTheoryBins               = len(self.theories[0]['binBoundaries']) - 1
        theoryBinBoundaries       = self.theories[0]['binBoundaries']

        # Get list of squared terms and unique combinations (works as expected also for 1 couplings cases)
        couplingCombinations = []
        couplingCombinations.extend( [ [ coupling, coupling ] for coupling in couplings ] )
        couplingCombinations.extend( [ list(couplingTuple) for couplingTuple in itertools.combinations( couplings, 2 ) ] )
        if self.includeLinearTerms:
            # Include also singular coupling terms and a constant
            couplingCombinations.extend( [ [coupling] for coupling in couplings ] )
            couplingCombinations.append( [] )


        # Number of theories necessary to perform a parametrization
        nComponents = sum(range(nCouplings+1))
        if self.includeLinearTerms:
            nComponents += len(couplings) + 1

        if not len(couplingCombinations) == nComponents:
            print 'ERROR: amount of coupling combinations should be equal to number of expected parameters'
            sys.exit()

        if len(self.theories) > nComponents:
            print '[FIXME!] Need to choose {0} theories with DIFFERENT COUPLINGS, otherwise matrix is singular!'.format( nComponents )
            print '\n{0} theories supplied but only {1} needed for parametrization; taking only the first {1} theories:'.format(
                len(self.theories), nComponents )
            self.theories = self.theories[:nComponents]
            for theory in self.theories:
                print '    ' + ', '.join([ '{0:10} = {1:10}'.format( cName, cValue ) for cName, cValue in theory['couplings'].iteritems() ])
        elif len(self.theories) < nComponents:
            print 'ERROR: cannot parametrize because number of supplied theories is too small (need {0} theories, found {1})'.format(
                nComponents, len(self.theories) )
            sys.exit()


        # ======================================
        # Find the parametrizations for all nTheorybins, plus 1 underflow and 1 overflow parametrization

        # First make sure there are RooRealVars for all the couplings in the workspace
        for coupling in couplings:
            self.modelBuilder.doVar(
                '{coupling}[{default},{down},{up}]'.format(
                    coupling = coupling,
                    default  = self.SMDict['couplings'][coupling],
                    down     = self.SMDict['couplings'][coupling] - 1000.0,
                    up       = self.SMDict['couplings'][coupling] + 1000.0,
                    )
                )

        # Calculate the coupling matrix
        couplingMat = []
        for theory in self.theories:
            if self.includeLinearTerms:
                row = []
                for couplingList in couplingCombinations:
                    product = 1.0
                    for coupling in couplingList:
                        product *= theory['couplings'][coupling]
                    row.append(product)
                couplingMat.append(row)
            else:
                couplingMat.append(
                    [ theory['couplings'][c1]*theory['couplings'][c2] for c1, c2 in couplingCombinations ]
                    )

        couplingMat = numpy.array( couplingMat )
        if self.verbose:
            print '\nSquared coupling terms:'
            print '  ', couplingCombinations
            print 'Found the following coupling matrix:'
            print couplingMat
        couplingMatInv = numpy.linalg.inv( couplingMat )


        # Create the actual parametrization for each theory bin
        parametrizations = []

        # Underflow (left extrapolation)
        parametrizations.append( [ 0. for i in xrange(nComponents) ] )

        # Add parametrization for each bin
        for iTheoryBin in xrange(nTheoryBins):
            ratios = numpy.array([ [theory['ratios'][iTheoryBin]] for theory in self.theories ])
            parametrization = list(itertools.chain.from_iterable( couplingMatInv.dot( ratios ) ))
            parametrizations.append( parametrization )

        # Use the SAME parametrization for the underflow bin as the first actually defined bin
        # Otherwise the zero contribution pulls down the integral
        parametrizations[0] = deepcopy( parametrizations[1] )

        # Overflow (right extrapolation)
        # parametrizations.append( [ 0. for i in xrange(nComponents) ] )
        parametrizations.append( parametrizations[-1][:] )

        if self.verbose:
            print '\nParametrizations per theory bin:'
            for i, parametrization in enumerate(parametrizations):
                print '{0:4}: {1}'.format( i, parametrization )
            print

        # Create a "@number" string per coupling
        argCoupling = { coupling : '@{0}'.format(iCoupling) for iCoupling, coupling in enumerate(couplings) }

        parametrizationNames = []
        for iParametrization, parametrization in enumerate( parametrizations ):

            # Also import it to the workspace
            parametrizationName = 'parametrization{0}'.format( iParametrization )
            parametrizationNames.append( parametrizationName )

            parametrizationString = ''
            for parameter, couplingList in zip( parametrization, couplingCombinations ):
                if len(couplingList) == 2:
                    coupling1, coupling2 = couplingList
                    parametrizationString += '+{0}*{1}*{2}'.format(
                        parameter, argCoupling[coupling1], argCoupling[coupling2]
                        )
                elif len(couplingList) == 1:
                    parametrizationString += '+{0}*{1}'.format(
                        parameter, argCoupling[couplingList[0]]
                        )
                elif len(couplingList) == 0:
                    parametrizationString += '+{0}'.format( parameter )

            # Strip off the '+' at the beginning
            parametrizationString = parametrizationString[1:]

            parametrizationExpression = (
                'expr::{name}('
                '"({string})", '
                '{arguments} )'
                ).format(
                    name      = parametrizationName,
                    string    = parametrizationString,
                    arguments = ','.join(couplings),
                    )

            if self.verbose > 1:
                print 'Importing parametrization {0}:'.format( iParametrization )
                print parametrizationExpression

            self.modelBuilder.factory_( parametrizationExpression )

        self.modelBuilder.out.defineSet( 'parametrizations', ','.join(parametrizationNames) )


        ########################################
        # Find the theory bin widths and parametrization indices inside the experimental bins
        ########################################

        if self.MakeLumiScalable:
            self.modelBuilder.doVar( 'lumiScale[1.0]' )
            self.modelBuilder.out.var('lumiScale').setConstant()


        SMXSInsideExperimentalBins = []
        self.yieldParameterNames = []
        for iExpBin in xrange(nExpBins):

            expBoundLeft  = expBinBoundaries[iExpBin]
            expBoundRight = expBinBoundaries[iExpBin+1]
            # expBinStr     = signals[iExpBin]

            if iExpBin == nExpBins-1:
                expBinStr     =  '{0}_{1}_GT{2}'.format(
                    productionMode, observableName,
                    expBoundLeft if not expBoundLeft.is_integer() else int(expBoundLeft)
                    )
            else:
                expBinStr     =  '{0}_{1}_{2}_{3}'.format(
                    productionMode, observableName,
                    expBoundLeft  if not expBoundLeft.is_integer() else int(expBoundLeft),
                    expBoundRight if not expBoundRight.is_integer() else int(expBoundRight)
                    )

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


            # Calculate total cross section (*not* /GeV) in experimental bin
            SMXSInsideExperimentalBin = 0.
            for iParametrization, binWidth in zip( parametrizationIndices, theoryBinWidthsInsideExperimentalBin ):
                SMXSInsideExperimentalBin += self.SMXS[iParametrization] * binWidth

                if self.verbose > 1:
                    print '\n    Adding contribution of theory bin {0} to SM cross section in exp bin {1}'.format( iParametrization, expBinStr )
                    print '      self.SMXS[iParametrization] = ', self.SMXS[iParametrization]
                    print '      binWidth                    = ', binWidth
                    print '      contribution                = ', self.SMXS[iParametrization] * binWidth

            if self.verbose:
                print '\nTotal cross section in {0} is {1}'.format( SMXSInsideExperimentalBin, expBinStr )
            SMXSInsideExperimentalBins.append( SMXSInsideExperimentalBin )

            # if self.verbose:
            #     print '   Which will carry the following weights:'
            #     print '   [ ' + ', '.join([ '{0:.4f}'.format(w) for w in theoryBinWidthsInsideExperimentalBin ]) + ']'

            if len(theoryBinWidthsInsideExperimentalBin) == 0:
                print '  Did not find any theoretical bin boundaries inside this experimental bin'
                print '   (i.e. can not do a parametrization here)'
                print '   Yield parameter will be 1.0'
                self.modelBuilder.doVar( 'r_{0}[1.0]'.format( expBinStr ))
                continue


            componentWeights = []
            for iComponent in xrange(nComponents):
                parametrizationWeights = []
                for iParametrization, binWidth in zip( parametrizationIndices, theoryBinWidthsInsideExperimentalBin ):
                    weight = (
                        ( binWidth * self.SMXS[iParametrization] )
                        /
                        SMXSInsideExperimentalBin
                        )
                    parametrizationWeights.append(weight)

                    if self.verbose > 1:
                        print '\n  Weight for component {0}, parametrization {1}'.format( iComponent, iParametrization )
                        print '    binWidth (theoryBin)    = ', binWidth
                        print '    SMXS/GeV (theoryBin)    = ', self.SMXS[iParametrization]
                        print '    SMXS (expBin, not /GeV) = ', SMXSInsideExperimentalBin
                        print '    ( binWidth * SMXS/GeV ) / SMXS_expBin = ', weight

                componentWeights.append( parametrizationWeights )


            # ======================================
            # Calculate the weighted average of components

            if self.verbose:
                print '\nCalculating weighted average components'

            averageComponents = []
            for iComponent in xrange(nComponents):

                if self.verbose > 1:
                    print '\n  Component {0}:'.format(iComponent)

                averageComponent = 0.
                parametrizationWeights = componentWeights[iComponent]
                for iParametrization, weight in zip( parametrizationIndices, parametrizationWeights ):
                    averageComponent += weight * parametrizations[iParametrization][iComponent]

                    if self.verbose > 1:
                        print '    Parametrization {0}'.format(iParametrization)
                        print '      weight          = ', weight
                        print '      parameter value = ', parametrizations[iParametrization][iComponent]
                        print '      product         = ', weight * parametrizations[iParametrization][iComponent]

                if self.verbose > 1:
                    print '  average Compontent {0} = {1}'.format( iComponent, averageComponent )

                averageComponents.append( averageComponent )


            # ======================================
            # Importing into WS is somewhat delicate

            argumentIndices = { coupling : '@{0}'.format(iCoupling) for iCoupling, coupling in enumerate(couplings) }

            averageComponentNames = []
            productNames = []
            yieldParameterFormula = []
            for iAverageComponent, averageComponent in enumerate(averageComponents):

                couplingList = couplingCombinations[iAverageComponent]
                nCouplingsForThisParameter = len(couplingList)

                averageComponentName = 'averageComponent_{signal}_{whichCouplings}'.format(
                    signal = expBinStr,
                    whichCouplings = ''.join(couplingList) if nCouplingsForThisParameter > 0 else 'CONSTANT',
                    )

                self.modelBuilder.doVar(
                    '{0}[{1}]'.format(
                        averageComponentName,
                        averageComponent,
                        )
                    )
                averageComponentNames.append( averageComponentName )


                # First add a "@number" entry in the argumentIndices dict, then use that in the formula
                argumentIndices[averageComponentName] = '@{0}'.format( len(argumentIndices) )

                yieldParameterFormulaComponent = argumentIndices[averageComponentName]
                for coupling in couplingList:
                    yieldParameterFormulaComponent += '*{0}'.format( argumentIndices[coupling] )
                yieldParameterFormula.append( yieldParameterFormulaComponent )

            # Compile expression string for final yieldParameter
            if not self.MakeLumiScalable:
                yieldParameterExpression = 'expr::r_{signal}( "({formulaString})", {commaSeparatedParameters} )'.format(
                    signal                   = expBinStr,
                    formulaString            = '+'.join(yieldParameterFormula),
                    commaSeparatedParameters = ','.join( couplings + averageComponentNames )
                    )
            else:
                argumentIndices['lumiScale'] = '@{0}'.format( len(argumentIndices) )
                yieldParameterExpression = 'expr::r_{signal}( "({formulaString})", {commaSeparatedParameters} )'.format(
                    signal                   = expBinStr,
                    formulaString            =  '{0}*({1})'.format( argumentIndices['lumiScale'], '+'.join(yieldParameterFormula) ),
                    commaSeparatedParameters = ','.join( couplings + averageComponentNames + ['lumiScale'] )
                    )

            if self.verbose:
                print 'Final yield parameter expression: {0}'.format(yieldParameterExpression)
                print '  Overview:'
                for key, value in argumentIndices.iteritems():
                    print '    {0:4} = {1}'.format( value, key )

            self.modelBuilder.factory_( yieldParameterExpression )

            if self.verbose:
                print '\nTest evaluation of yieldParameter:'
                for coupling in couplings:
                    self.modelBuilder.out.var(coupling).Print()
                self.modelBuilder.out.function( 'r_{0}'.format(expBinStr) ).Print()
                print ''

            self.yieldParameterNames.append( 'r_{0}'.format(expBinStr) )


        self.modelBuilder.out.defineSet( 'POI', ','.join(couplings) )

        # Define 2 extra sets for plotting convenience
        self.modelBuilder.out.defineSet( 'couplings', ','.join(couplings) )
        # self.modelBuilder.out.defineSet( 'yieldParameters', ','.join([ 'r_' + s for s in signals ]) )
        self.modelBuilder.out.defineSet( 'yieldParameters', ','.join(self.yieldParameterNames) )


        ########################################
        # Handle theory uncertainties
        ########################################

        def printMatrix( matrix, indent = '    ', scientific=True ):
            for row in matrix:
                print indent + ' '.join([ '{0:<+10.2{1}}'.format( number, 'E' if scientific else 'f' ) for number in row ])

        # First check if given input makes sense
        doDecorrelation = False
        if (
            not self.theoryUncertaintiesPassed
            and not self.correlationMatrixPassed
            and not self.covarianceMatrixPassed
            ):
            print 'No theory uncertainties are specified; Running without theory uncertainties'
            pass

        elif (
            self.theoryUncertaintiesPassed
            and self.correlationMatrixPassed
            and not self.covarianceMatrixPassed
            ):
            if not len(self.theoryUncertainties) == len(self.correlationMatrix):
                print '[ERROR] Cannot build a covariance matrix out of given input'
                print '        len(self.theoryUncertainties) = {0}'.format( len(self.theoryUncertainties) )
                print '        len(self.correlationMatrix)   = {0}'.format( len(self.correlationMatrix) )
                return

            self.covarianceMatrix = []
            self.nTheoryUncertainties = len(self.theoryUncertainties)
            for i in xrange(self.nTheoryUncertainties):
                self.covarianceMatrix.append(
                    [ self.theoryUncertainties[i] * self.theoryUncertainties[j] * self.correlationMatrix[i][j] for j in xrange(self.nTheoryUncertainties) ]
                    )
            doDecorrelation = True
            print '\nApplying theory uncertainties using the passed correlationMatrix and theoryUncertainties'

            print '  Using the following theory uncertainties:'
            for unc in self.theoryUncertainties:
                print '    {0:<+20.8f}'.format( unc )

            print '  Using the following correlation matrix:'
            printMatrix( self.correlationMatrix, scientific=False )


        elif (
            not self.theoryUncertaintiesPassed
            and not self.correlationMatrixPassed
            and self.covarianceMatrixPassed
            ):
            # self.covarianceMatrix should be filled already
            doDecorrelation = True
            self.nTheoryUncertainties = len(self.covarianceMatrix)
            print '\nApplying theory uncertainties using the passed covarianceMatrix'

        else:
            print '[ERROR] Given input makes no sense:'
            print '        self.theoryUncertaintiesPassed = ', self.theoryUncertaintiesPassed
            print '        self.correlationMatrixPassed   = ', self.correlationMatrixPassed
            print '        self.covarianceMatrixPassed    = ', self.covarianceMatrixPassed
            return


        if doDecorrelation:

            print '  Using the following covariance matrix:'
            printMatrix( self.covarianceMatrix )

            decorrelatedMatrix = self.Decorrelate( self.covarianceMatrix )

            print '  Found the following decorrelatedMatrix:'
            printMatrix( decorrelatedMatrix )

            print '  Divided by SM cross section:'
            decorrelatedMatrixNormalized = deepcopy( decorrelatedMatrix )
            for i in xrange( min( self.nTheoryUncertainties, len(SMXSInsideExperimentalBins) ) ):
                for j in xrange( min( self.nTheoryUncertainties, len(SMXSInsideExperimentalBins) ) ):
                    decorrelatedMatrixNormalized[i][j] = 1 + decorrelatedMatrix[i][j] / SMXSInsideExperimentalBins[i]
                    # # Skip some uncertainties if the effect is too small
                    # if abs(decorrelatedMatrixNormalized[i][j] - 1) < 0.0005:
                    #     decorrelatedMatrixNormalized[i][j] = 0.
            printMatrix( decorrelatedMatrixNormalized )

            # for nuis in self.DC.systs:
            #     print nuis[0], nuis[1], nuis[2], nuis[3]

            for iTheoryUncertainty in xrange(self.nTheoryUncertainties):

                errDict = {}
                for binName in self.DC.bins:
                    errDict[binName] = {}
                    for processName in self.DC.processes:

                        errDict[binName][processName] = 0.

                        if processName in signals:
                            iProcess = signals.index(processName)
                            if iProcess < self.nTheoryUncertainties:
                                errDict[binName][processName] = decorrelatedMatrixNormalized[iTheoryUncertainty][iProcess]
                        
                systematicName = 'theoryUncertainty_{0}'.format( signals[iTheoryUncertainty] )

                self.DC.systs.append(
                    ( systematicName, False, 'lnN', [], errDict )
                    )

                if self.verbose:
                    print 'Added nuisance \'{0}\''.format( systematicName )
                        


        # Import also the theory and exp binBoundaries into workspace
        # theoryBinBoundariesTVectorD = ROOT.TVectorD( len(theoryBinBoundaries) )
        # for i, theoryBinBoundary in enumerate(theoryBinBoundaries):
        #     theoryBinBoundariesTVectorD[i] = theoryBinBoundary


        theoryBinBoundarySet = []
        for i, theoryBinBoundary in enumerate(theoryBinBoundaries):
            self.modelBuilder.doVar( 'theoryBinBound{0}[{1}]'.format( i, theoryBinBoundary ) )
            theoryBinBoundarySet.append( 'theoryBinBound{0}'.format(i) )
        self.modelBuilder.out.defineSet( 'theoryBinBoundaries', ','.join(theoryBinBoundarySet) )

        expBinBoundarySet = []
        for i, expBinBoundary in enumerate(expBinBoundaries):
            self.modelBuilder.doVar( 'expBinBound{0}[{1}]'.format( i, expBinBoundary ) )
            expBinBoundarySet.append( 'expBinBound{0}'.format(i) )
        self.modelBuilder.out.defineSet( 'expBinBoundaries', ','.join(expBinBoundarySet) )

        self.chapter( 'Starting model.getYieldScale()' )


    def getYieldScale( self, bin, process ):

        def p( yieldScale ):
            if self.verbose > 0:
                print 'Getting yield scale: bin: {0:30} | process: {1:16} | yieldScale: {2}'.format( bin, process, yieldScale )

        one = 1 if not self.MakeLumiScalable else 'lumiScale'

        if not self.DC.isSignal[process]:
            p( one )
            return one
        else:
            if process == 'OutsideAcceptance':
                p( one )
                return one
            elif self.splitggH and not 'ggH' in process:
                p( one )
                return one
            else:

                if not( self.manualExpBinBoundaries is None ):

                    # Try to determine bin boundaries from process name
                    match = re.search( r'([\dmp]+)_([\dmp]+)', process )

                    if not match:
                        p( one )
                        return one
                    else:
                        leftBound  = float(match.group(1))
                        rightBound = float(match.group(2))

                        for iBin in xrange(self.nExpBins):
                            if leftBound >= self.expBinBoundaries[iBin] and rightBound <= self.expBinBoundaries[iBin+1]:
                                yieldParameter = self.yieldParameterNames[iBin]

                                # Still check if the bin is supposed to be skipped
                                for skipBin in self.skipBins:
                                    if skipBin in process:
                                        yieldParameter = one

                                p(yieldParameter)
                                return yieldParameter

                        else:
                            print 'ERROR: Could not find an appropriate yield parameter for \'{0}\''.format( process )
                            sys.exit()

                else:
                    # Case for when no manual bin boundaries were specified

                    yieldParameter = 'r_' + process

                    for skipBin in self.skipBins:
                        if skipBin in process:
                            yieldParameter = one

                    p( yieldParameter )
                    return yieldParameter



    def getYieldScaleOLD( self, bin, process ):

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

                if not( self.manualExpBinBoundaries is None ):

                    # Try to determine bin boundaries from process name
                    match = re.search( r'([\dmp]+)_([\dmp]+)', process )

                    if not match:
                        p( 1 )
                        return 1
                    else:
                        leftBound  = float(match.group(1))
                        rightBound = float(match.group(2))

                        for iBin in xrange(self.nExpBins):
                            if leftBound >= self.expBinBoundaries[iBin] and rightBound <= self.expBinBoundaries[iBin+1]:
                                yieldParameter = self.yieldParameterNames[iBin]

                                # Still check if the bin is supposed to be skipped
                                for skipBin in self.skipBins:
                                    if skipBin in process:
                                        yieldParameter = 1

                                p(yieldParameter)
                                return yieldParameter

                        else:
                            print 'ERROR: Could not find an appropriate yield parameter for \'{0}\''.format( process )
                            sys.exit()

                else:
                    # Case for when no manual bin boundaries were specified

                    yieldParameter = 'r_' + process

                    for skipBin in self.skipBins:
                        if skipBin in process:
                            yieldParameter = 1

                    p( yieldParameter )
                    return yieldParameter




couplingModel=CouplingModel()

