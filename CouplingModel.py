
from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder

import os, sys, numpy, itertools, re
import ROOT
from math import sqrt
from copy import deepcopy
from array import array

from time import sleep

class CouplingModelError(Exception):
    pass

class CouplingModel( PhysicsModel ):
    ''' Model used to unfold differential distributions '''

    def __init__(self):
        PhysicsModel.__init__(self)

        self.mHRange=[]
        self.mass = 125.
        self.verbose = 1

        self.Persistence = []

        self.theories = []

        self.SMbinBoundaries = []
        self.SMXS = []

        self.isOnlyHgg = False
        self.isOnlyHZZ = False

        self.includeLinearTerms   = False
        self.splitggH             = False
        self.MakeLumiScalable     = False
        self.FitBR                = False
        self.ProfileTotalXS       = False
        self.FitOnlyNormalization = False
        self.FitRatioOfBRs        = False
        self.DoBRUncertainties    = False

        self.theoryUncertaintiesPassed = False
        self.correlationMatrixPassed   = False
        self.covarianceMatrixPassed    = False
        self.skipOverflowBinTheoryUncertainty = False

        self.manualExpBinBoundaries = []
        self.skipBins = []
        self.skipOverflowBin = False
        

        self.SM_HIGG_DECAYS = [ "hww", "hzz", "hgg", "htt", "hbb", 'hzg', 'hmm', 'hcc', 'hgluglu' ]
        self.SM_HIGG_DECAYS.append( 'hss' )


    ################################################################################
    # Standard methods
    ################################################################################

    #____________________________________________________________________
    def setPhysicsOptions( self, physOptions ):
        self.chapter( 'Starting model.setPhysicsOptions()' )

        for physOption in physOptions:
            optionName, optionValue = physOption.split('=',1)

            if False:
                pass

            elif optionName == 'lumiScale':
                self.MakeLumiScalable = eval(optionValue)

            elif optionName == 'FitBR':
                self.FitBR = eval(optionValue)

            elif optionName == 'ProfileTotalXS':
                self.ProfileTotalXS = eval(optionValue)

            elif optionName == 'FitOnlyNormalization':
                self.FitOnlyNormalization = eval(optionValue)

            elif optionName == 'FitRatioOfBRs':
                self.FitRatioOfBRs = eval(optionValue)

            elif optionName == 'splitggH':
                self.splitggH = eval(optionValue)

            elif optionName == 'isOnlyHZZ':
                self.isOnlyHZZ = eval(optionValue)

            elif optionName == 'isOnlyHgg':
                self.isOnlyHgg = eval(optionValue)

            elif optionName == 'DoBRUncertainties':
                self.DoBRUncertainties = eval(optionValue)

            elif optionName == 'binBoundaries':
                self.manualExpBinBoundaries = [ float(i) for i in optionValue.split(',') ]

            elif optionName == 'verbose':
                self.verbose = int(optionValue)

            elif optionName == 'skipBins':
                self.skipBins = optionValue.split(',')

                for skipBin in self.skipBins:
                    if 'GT' in skipBin:
                        self.skipOverflowBin = True
                        break


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
                # self.SMXS = [ self.SMDict['crosssection'][0] ] + self.SMDict['crosssection'] + [ self.SMDict['crosssection'][-1] ]
                self.SMXS = [ 0.0 ] + self.SMDict['crosssection']


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

            elif optionName == 'skipOverflowBinTheoryUncertainty':
                self.skipOverflowBinTheoryUncertainty = eval(optionValue)

            else:
                print 'Unknown option \'{0}\'; skipping'.format(optionName)
                continue


        if self.FitBR and self.FitRatioOfBRs:
            raise CouplingModelError( 'Options FitBR and FitRatioOfBRs are exclusive!' )


    #____________________________________________________________________
    def doParametersOfInterest(self):
        self.chapter( 'Starting model.doParametersOfInterest()' )

        self.setMH()

        if self.FitOnlyNormalization:
            self.makeParametrizationsFromTheory()
            self.MakeTotalXSExpressions()
            self.modelBuilder.doVar('bkg_modifier[1.0]')
            self.modelBuilder.out.var('bkg_modifier').setConstant(True)
            self.modelBuilder.doVar('hgg_xH_modifier[1.0]')
            self.modelBuilder.out.var('hgg_xH_modifier').setConstant(True)
            self.modelBuilder.doVar('hzz_xH_modifier[1.0]')
            self.modelBuilder.out.var('hzz_xH_modifier').setConstant(True)

            self.modelBuilder.out.defineSet( 'POI', ','.join(self.couplings) )

        else:
            self.figureOutBinning()

            self.makeParametrizationsFromTheory()

            if self.FitBR:
                self.MakeWidthExpressions()
        
            if self.ProfileTotalXS:
                self.MakeTotalXSExpressions()

            if self.FitRatioOfBRs:
                self.MakeRatioOfBRsExpressions()

            self.defineYieldParameters()

            self.addTheoryUncertaintyNuisances()

        self.chapter( 'Starting model.getYieldScale()' )


    #____________________________________________________________________
    def getYieldScale( self, bin, process ):

        if self.verbose:
            print 'Getting scale for process = {0:21}, bin = {1}'.format( process, bin )

        # 'hgg_xH_modifier',
        # 'hzz_xH_modifier',
        # 'hgg_OutsideAcceptance_modifier',
        # 'hzz_OutsideAcceptance_modifier',
        # 'bkg_modifier',

        if not self.DC.isSignal[process]:
            yieldParameter = 'bkg_modifier'

        elif 'OutsideAcceptance' in process:
            yieldParameter = 'hgg_OutsideAcceptance_modifier' if self.BinIsHgg(bin) else 'hzz_OutsideAcceptance_modifier'

            if self.FitOnlyNormalization:
                yieldParameter = 'globalTotalXSmodifier'

        else:
            if self.splitggH and 'xH' in process:
                yieldParameter = 'hgg_xH_modifier' if self.BinIsHgg(bin) else 'hzz_xH_modifier'

            elif self.FitOnlyNormalization:
                yieldParameter = 'globalTotalXSmodifier'

            elif ( self.splitggH and 'ggH' in process ) or ( 'smH' in process ):
                # If it's a ggH process

                isOverflowBin = False
                isInsideUserBinBoundaries = False

                matchRegularBin  = re.search( r'([\dmp]+)_([\dmp]+)', process )
                matchOverflowBin = re.search( r'([GTLE]+)([\dmp]+)', process )

                if matchRegularBin:
                    leftBound  = float(matchRegularBin.group(1).replace('m','-').replace('p','.'))
                    rightBound = float(matchRegularBin.group(2).replace('m','-').replace('p','.'))
                    for iBin in xrange(self.nExpBins):
                        if leftBound >= self.expBinBoundaries[iBin] and rightBound <= self.expBinBoundaries[iBin+1]:
                            yieldParameter = self.hgg_yieldParameterNames[iBin] if self.BinIsHgg(bin) else self.hzz_yieldParameterNames[iBin]
                            break
                    else:
                        yieldParameter = 'hgg_OutsideAcceptance_modifier' if self.BinIsHgg(bin) else 'hzz_OutsideAcceptance_modifier'

                elif matchOverflowBin:
                    # bound = float(matchOverflowBin.group(1).replace('m','-').replace('p','.'))
                    lastYieldParameter = self.hgg_yieldParameterNames[-1] if self.BinIsHgg(bin) else self.hzz_yieldParameterNames[-1]
                    if matchOverflowBin.group() in lastYieldParameter:
                        yieldParameter = lastYieldParameter
                    else:
                        yieldParameter = 'hgg_OutsideAcceptance_modifier' if self.BinIsHgg(bin) else 'hzz_OutsideAcceptance_modifier'

                else:
                    raise CouplingModelError( 'Failure for process \'{0}\': Could not extract any bin boundary information'.format(process) )

            else:
                raise CouplingModelError( 'Failure for process \'{0}\': Production process is not \'xH\', \'ggH\' or \'smH\''.format(process) )


        if self.verbose:
            print '    --> Scaling with \'{0}\''.format( yieldParameter )
            
            # print '          test print:'
            # try:
            #     self.modelBuilder.out.var(yieldParameter).Print()
            # except ReferenceError:
            #     try:
            #         self.modelBuilder.out.function(yieldParameter).Print()
            #     except ReferenceError:
            #         print '          yieldParameter \'{0}\' does not seem to be in the ws!!'.format(yieldParameter)

        return yieldParameter



    ################################################################################
    # Split parts of doParametersOfInterest
    ################################################################################

    #____________________________________________________________________
    def setMH( self ):

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


    #____________________________________________________________________
    def makeParametrizationsFromTheory( self ):

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
            raise CouplingModelError( 'ERROR: amount of coupling combinations should be equal to number of expected parameters' )
            sys.exit()

        if len(self.theories) > nComponents:
            print '[FIXME!] Need to choose {0} theories with DIFFERENT COUPLINGS, otherwise matrix is singular!'.format( nComponents )
            print '\n{0} theories supplied but only {1} needed for parametrization; taking only the first {1} theories:'.format(
                len(self.theories), nComponents )
            self.theories = self.theories[:nComponents]
            for theory in self.theories:
                print '    ' + ', '.join([ '{0:10} = {1:10}'.format( cName, cValue ) for cName, cValue in theory['couplings'].iteritems() ])
        elif len(self.theories) < nComponents:
            raise CouplingModelError( 'ERROR: cannot parametrize because number of supplied theories is too small (need {0} theories, found {1})'.format(
                nComponents, len(self.theories) ) )


        # ======================================
        # Find the parametrizations for all nTheoryBins, plus 1 underflow and 1 overflow parametrization

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
        #   Sep28: Actually just put 0.0. The underflow bin is considered to have no contribution
        # parametrizations[0] = deepcopy( parametrizations[1] )

        # Overflow (right extrapolation)
        #   Sep28: Stop adding the overflow parametrization.
        #          Just use whatever parametrizations there are in the overflow bin, or return 1.0
        # parametrizations.append( parametrizations[-1][:] )

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

        # Attach to class

        self.parametrizations     = parametrizations

        self.couplingCombinations = couplingCombinations
        self.nComponents          = nComponents

        self.nCouplings           = nCouplings
        self.couplings            = couplings
        self.nTheoryBins          = nTheoryBins
        self.theoryBinBoundaries  = theoryBinBoundaries


        theoryBinBoundarySet = []
        for i, theoryBinBoundary in enumerate(self.theoryBinBoundaries):
            self.modelBuilder.doVar( 'theoryBinBound{0}[{1}]'.format( i, theoryBinBoundary ) )
            theoryBinBoundarySet.append( 'theoryBinBound{0}'.format(i) )
        self.modelBuilder.out.defineSet( 'theoryBinBoundaries', ','.join(theoryBinBoundarySet) )



    #____________________________________________________________________
    def figureOutBinning( self ):

        bins    = self.DC.bins
        signals = self.DC.signals
        signals = [ s for s in signals if not 'OutsideAcceptance' in s ]
        
        # Determine all available production modes
        productionModes = list(set([ s.split('_')[0] for s in signals ]))

        if self.splitggH:
            signals = [ s for s in signals if 'ggH' in s ]

        signals.sort( key = lambda n: float(
            n.split('_')[2].replace('p','.').replace('m','-').replace('GT','').replace('GE','').replace('LT','').replace('LE','') ) )


        # ======================================
        # Figure out experimental binning from signals

        productionMode, observableName = signals[0].split('_')[:2]

        if len(self.manualExpBinBoundaries) > 0:
            expBinBoundaries = self.manualExpBinBoundaries
            # Overflow bin should be added on command line, is safer
            # expBinBoundaries.append( 500. )
            nExpBins = len(expBinBoundaries) - 1
            if self.verbose:
                print 'Taking manually specified bin boundaries:', expBinBoundaries

        else:
            raise CouplingModelError(
                'Specify \'--PO binBoundaries=0,15,...\' on the command line. Automatic boundary determination no longer supported.'
                )

        self.allProcessBinBoundaries = []
        for signal in signals:
            if 'OutsideAcceptance' in signal: continue
            for bound in signal.split('_')[2:]:
                bound = bound.replace('p','.').replace('m','-').replace('GT','').replace('GE','').replace('LT','').replace('LE','')
                self.allProcessBinBoundaries.append( float(bound) )

        # Attach to class so they are accesible in other methods
        self.expBinBoundaries = expBinBoundaries
        self.nExpBins         = nExpBins

        self.signals          = signals
        self.productionMode   = productionMode
        self.observableName   = observableName


        expBinBoundarySet = []
        for i, expBinBoundary in enumerate(self.expBinBoundaries):
            if self.skipOverflowBin and expBinBoundary == expBinBoundaries[-1]: continue
            self.modelBuilder.doVar( 'expBinBound{0}[{1}]'.format( i, expBinBoundary ) )
            expBinBoundarySet.append( 'expBinBound{0}'.format(i) )
        self.modelBuilder.out.defineSet( 'expBinBoundaries', ','.join(expBinBoundarySet) )



    #____________________________________________________________________
    def defineYieldParameters( self ):
        self.chapter( 'Starting model.defineYieldParameters()' )

        # Define some variables in local for convenience
        nExpBins            = self.nExpBins
        expBinBoundaries    = self.expBinBoundaries
        nTheoryBins         = self.nTheoryBins
        theoryBinBoundaries = self.theoryBinBoundaries
        parametrizations    = self.parametrizations


        # ======================================
        # Extra variables for further studies

        self.hgg_yieldParameterNames = []
        self.hzz_yieldParameterNames = []

        if self.MakeLumiScalable:
            self.modelBuilder.doVar( 'lumiScale[1.0]' )
            self.modelBuilder.out.var('lumiScale').setConstant()


        # ======================================
        # Start of loop over experimental bins

        SMXSInsideExperimentalBins = []
        self.yieldParameterNames = []
        for iExpBin in xrange(nExpBins):

            expBoundLeft  = expBinBoundaries[iExpBin]
            expBoundRight = expBinBoundaries[iExpBin+1]

            # Determine if bin is the overflow bin
            isOverflowBin = False
            if expBoundLeft == self.allProcessBinBoundaries[-1]:
                isOverflowBin = True

            # if iExpBin == nExpBins-1:
            #     expBinStr     =  '{0}_{1}_GT{2}'.format(
            #         self.productionMode, self.observableName,
            #         expBoundLeft if not expBoundLeft.is_integer() else int(expBoundLeft)
            #         )
            if isOverflowBin:
                expBinStr     = '{0}_{1}_GT{2}'.format(
                    self.productionMode, self.observableName,
                    expBoundLeft if not expBoundLeft.is_integer() else int(expBoundLeft)
                    )
            else:
                expBinStr     =  '{0}_{1}_{2}_{3}'.format(
                    self.productionMode, self.observableName,
                    expBoundLeft  if not expBoundLeft.is_integer() else int(expBoundLeft),
                    expBoundRight if not expBoundRight.is_integer() else int(expBoundRight)
                    )


            if self.verbose:
                print '\n' + '='*80
                print 'Processing bin {0}: {1}'.format( iExpBin, expBinStr )

            theoryBinBoundariesInsideExperimentalBin = []
            parametrizationIndices                   = []

            # Extrapolation on left of theory spectrum
            if expBoundLeft < theoryBinBoundaries[0]:
                theoryBinBoundariesInsideExperimentalBin.append( expBoundLeft )
                parametrizationIndices.append( 0 )


            # ======================================
            # Determine which parametrizations enter into the experimental bin

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
            # Sep28: Stop doing this
            #        Use only the parametrizations defined by theorists
            if expBoundRight > theoryBinBoundaries[-1]:
                # theoryBinBoundariesInsideExperimentalBin.append( expBoundRight )
                # parametrizationIndices.append( nTheoryBins+1 )
                theoryBinBoundariesInsideExperimentalBin.append( theoryBinBoundaries[-1] )


            nTheoryBinsInsideExperimentalBin = len(theoryBinBoundariesInsideExperimentalBin)-1
            theoryBinWidthsInsideExperimentalBin = [
                theoryBinBoundariesInsideExperimentalBin[i+1] - theoryBinBoundariesInsideExperimentalBin[i] for i in xrange(nTheoryBinsInsideExperimentalBin)
                ]


            if self.verbose:
                # print '\nProcessing experimental bin {0}'.format(expBinStr)
                print 'Found the following theory bin boundaries inside experimental bin {0} to {1}:'.format(
                    expBoundLeft, expBoundRight )
                print '  ', theoryBinBoundariesInsideExperimentalBin
                print '   The following parametrizations will be used to evaluate the cross section:'
                print '  ', parametrizationIndices


            if self.verbose:
                print '\n' + '-'*60
                print 'Calculating cross section'

            # Calculate total cross section (*not* /GeV) in experimental bin
            SMXSInsideExperimentalBin = 0.
            for iParametrization, binWidth in zip( parametrizationIndices, theoryBinWidthsInsideExperimentalBin ):
                SMXSInsideExperimentalBin += self.SMXS[iParametrization] * binWidth

                if self.verbose > 1:
                    print '\n    Adding contribution of theory bin {0} to SM cross section in exp bin {1}'.format( iParametrization, expBinStr )
                    print '    contribution = {2:+8.3f} ( self.SMXS[iParametrization] = {0:+8.3f} * binWidth = {1:+8.3f} )'.format(
                        self.SMXS[iParametrization] ,
                        binWidth ,
                        self.SMXS[iParametrization] * binWidth ,
                        )

            if self.verbose:
                print '\nTotal cross section in {1} is {0}'.format( SMXSInsideExperimentalBin, expBinStr )


            if len(theoryBinWidthsInsideExperimentalBin) == 0:
                print '  Did not find any theoretical bin boundaries inside this experimental bin'
                print '   (i.e. can not do a parametrization here)'
                print '   Yield parameter will be 1.0'
                self.modelBuilder.doVar( 'r_{0}[1.0]'.format( expBinStr ))
                continue

            SMXSInsideExperimentalBins.append( SMXSInsideExperimentalBin )


            if self.verbose:
                print '\n' + '-'*60
                print 'Creating expression for yield parameter'

            componentWeights = []
            for iComponent in xrange(self.nComponents):
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
                        print '    weight = {0:+8.3f} ( ( binWidth[{1:+8.3f}] * SMXS/GeV[{2:+8.3f}] ) / SMXS_expBin[{3:+8.3f}] )'.format(
                            weight, binWidth, self.SMXS[iParametrization], SMXSInsideExperimentalBin
                            )

                componentWeights.append( parametrizationWeights )


            # ======================================
            # Calculate the weighted average of components

            if self.verbose:
                print '\nCalculating weighted average components'

            averageComponents = []
            for iComponent in xrange(self.nComponents):

                if self.verbose > 1:
                    print '\n  Component {0} ('.format(iComponent), self.couplingCombinations[iComponent], '):'

                averageComponent = 0.
                parametrizationWeights = componentWeights[iComponent]
                for iParametrization, weight in zip( parametrizationIndices, parametrizationWeights ):
                    averageComponent += weight * parametrizations[iParametrization][iComponent]

                    if self.verbose > 1:
                        print '    Parametrization {0:3}: product = {1:+8.3f} ( weight[{2:+8.3f}] * parameterValue[{3:+8.3f}]'.format(
                            iParametrization, 
                            weight * parametrizations[iParametrization][iComponent],
                            weight,
                            parametrizations[iParametrization][iComponent]
                            )

                if self.verbose > 1:
                    print '  average Component {0} = {1}'.format( iComponent, averageComponent )

                averageComponents.append( averageComponent )


            if self.verbose:
                print '\nThe determined polynomial for bin {0} is:\n'.format(expBinStr)

                line = 'r_{0} = '.format(expBinStr)
                for iComponent in xrange(self.nComponents):

                    line += '{0}*{1} + '.format(
                        averageComponents[iComponent],
                        '*'.join(self.couplingCombinations[iComponent])
                        )
                print line [:-2]


            # ======================================
            # Importing into WS is somewhat delicate

            argumentIndices = { coupling : '@{0}'.format(iCoupling) for iCoupling, coupling in enumerate(self.couplings) }

            averageComponentNames = []
            productNames = []
            yieldParameterFormula = []
            for iAverageComponent, averageComponent in enumerate(averageComponents):

                couplingList = self.couplingCombinations[iAverageComponent]
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
            yieldParameterExpression = 'expr::r_{signal}( "({formulaString})", {commaSeparatedParameters} )'.format(
                signal                   = expBinStr,
                formulaString            = '+'.join(yieldParameterFormula),
                commaSeparatedParameters = ','.join( self.couplings + averageComponentNames )
                )

            if self.verbose:
                print '\nFinal yield parameter expression: {0}'.format(yieldParameterExpression)
                print '  Overview:'
                for key, value in argumentIndices.iteritems():
                    print '    {0:4} = {1}'.format( value, key )

            self.modelBuilder.factory_( yieldParameterExpression )

            if self.verbose:
                print '\nTest evaluation of yieldParameter:'
                for coupling in self.couplings:
                    self.modelBuilder.out.var(coupling).Print()
                self.modelBuilder.out.function( 'r_{0}'.format(expBinStr) ).Print()
                print ''

            self.yieldParameterNames.append( 'r_{0}'.format(expBinStr) )


            # ======================================
            # Add other modifiers

            hgg_modifiers = [ 'r_{0}'.format(expBinStr) ]
            hzz_modifiers = [ 'r_{0}'.format(expBinStr) ]

            if self.MakeLumiScalable:
                hgg_modifiers.append( 'lumiScale' )
                hzz_modifiers.append( 'lumiScale' )

            if self.FitBR:
                hgg_modifiers.append( 'hggBRmodifier' )
                hzz_modifiers.append( 'hzzBRmodifier' )

            if self.ProfileTotalXS:
                hgg_modifiers.append( 'totalXSmodifier' )
                hzz_modifiers.append( 'totalXSmodifier' )

            if self.FitRatioOfBRs:
                hgg_modifiers.append( 'hgg_ratioBRmodifier' )
                hzz_modifiers.append( 'hzz_ratioBRmodifier' )

            self.modelBuilder.factory_( 'prod::r_hgg_{0}({1})'.format( expBinStr, ','.join(hgg_modifiers) ) )
            self.modelBuilder.factory_( 'prod::r_hzz_{0}({1})'.format( expBinStr, ','.join(hzz_modifiers) ) )

            self.hgg_yieldParameterNames.append( 'r_hgg_{0}'.format(expBinStr) )
            self.hzz_yieldParameterNames.append( 'r_hzz_{0}'.format(expBinStr) )


        # ======================================
        # Other yield parameters: xH, OutsideAcceptance, and background

        self.modelBuilder.doVar( 'one[1]' )
        self.modelBuilder.out.var('one').setConstant(True)

        # xH
        hgg_xH_modifier = [ 'one' ]
        hzz_xH_modifier = [ 'one' ]

        if self.MakeLumiScalable:
            hgg_xH_modifier.append( 'lumiScale' )
            hzz_xH_modifier.append( 'lumiScale' )
        if self.FitBR:
            hgg_xH_modifier.extend([ 'xH_modifier', 'hggBRmodifier' ])
            hzz_xH_modifier.extend([ 'xH_modifier', 'hzzBRmodifier' ])

        self.modelBuilder.factory_( 'prod::hgg_xH_modifier({0})'.format( ','.join(hgg_xH_modifier) ) )
        self.modelBuilder.factory_( 'prod::hzz_xH_modifier({0})'.format( ','.join(hzz_xH_modifier) ) )

        # OutsideAcceptance
        hgg_OutsideAcceptance_modifier = [ 'one' ]
        hzz_OutsideAcceptance_modifier = [ 'one' ]

        if self.FitBR:
            hgg_OutsideAcceptance_modifier.append( 'hggBRmodifier' )
            hzz_OutsideAcceptance_modifier.append( 'hzzBRmodifier' )
        if self.MakeLumiScalable:
            hgg_OutsideAcceptance_modifier.append( 'lumiScale' )
            hzz_OutsideAcceptance_modifier.append( 'lumiScale' )
        
        self.modelBuilder.factory_( 'prod::hgg_OutsideAcceptance_modifier({0})'.format( ','.join(hgg_OutsideAcceptance_modifier) ) )
        self.modelBuilder.factory_( 'prod::hzz_OutsideAcceptance_modifier({0})'.format( ','.join(hzz_OutsideAcceptance_modifier) ) )

        # bkg
        bkg_modifier = [ 'one' ]
        if self.MakeLumiScalable:
            bkg_modifier.append( 'lumiScale' )
        self.modelBuilder.factory_( 'prod::bkg_modifier({0})'.format( ','.join(bkg_modifier) ) )


        # ======================================
        # Save some sets to ws

        self.modelBuilder.out.defineSet( 'POI', ','.join(self.couplings) )
        self.modelBuilder.out.defineSet( 'couplings', ','.join(self.couplings) )
        self.modelBuilder.out.defineSet( 'yieldParameters', ','.join(self.yieldParameterNames) )
        self.modelBuilder.out.defineSet( 'hgg_yieldParameters', ','.join(self.hgg_yieldParameterNames) )
        self.modelBuilder.out.defineSet( 'hzz_yieldParameters', ','.join(self.hgg_yieldParameterNames) )

        self.modelBuilder.out.defineSet( 'modifiers', ','.join([
            'hgg_xH_modifier',
            'hzz_xH_modifier',
            'hgg_OutsideAcceptance_modifier',
            'hzz_OutsideAcceptance_modifier',
            'bkg_modifier',
            ]))

        self.SMXSInsideExperimentalBins = SMXSInsideExperimentalBins

        if self.FitRatioOfBRs:
            self.modelBuilder.out.set('POI').add( self.modelBuilder.out.var('hgg_ratioBRmodifier') )
            self.modelBuilder.out.set('POI').add( self.modelBuilder.out.var('ratio_BR_hgg_hzz') )


        if self.verbose:

            print '\n\n*** Done defining all yield parameters; control printout:'

            allVarsToPrint = self.hgg_yieldParameterNames + self.hzz_yieldParameterNames + [
                'hgg_xH_modifier',
                'hzz_xH_modifier',
                'hgg_OutsideAcceptance_modifier',
                'hzz_OutsideAcceptance_modifier',
                'bkg_modifier',
                ]

            for name in allVarsToPrint:
                # print ''
                self.modelBuilder.out.function(name).Print()

            print '\n'


    #____________________________________________________________________
    def addTheoryUncertaintyNuisances( self ):
    
        if self.verbose:
            print '\n' + '='*80
            print 'Inserting theory uncertainties\n'

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
                raise CouplingModelError()

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
            raise CouplingModelError()


        # ======================================
        # Align input matrix with the chosen binning for the yield parameters

        # Avoid dimension mismatch: The theory covariance matrix may have more bins then specified for the fit
        # This is a problem because the covariance matrix has to be normalized by the SM cross section,
        # which is not calculated for every bin if the theory covariance matrix has more bins then specified for the fit
        self.nTheoryUncertainties = min( self.nTheoryUncertainties, len(self.SMXSInsideExperimentalBins) )

        # In some cases the overflow bin has a huge (~400%) uncertainty; causes fit instabilities, needs to be skipped
        if self.skipOverflowBinTheoryUncertainty:
            self.nTheoryUncertainties -= 1
            print '\n[WARNING]: Not registering the nuisance parameter related to the last theory bin (self.skipOverflowBinTheoryUncertainty==True)\n'

        self.covarianceMatrix = [ line[:self.nTheoryUncertainties] for line in self.covarianceMatrix[:self.nTheoryUncertainties] ]


        # ======================================
        # Run ICA and write nuisances to output workspace

        if doDecorrelation:

            print '  Using the following covariance matrix:'
            printMatrix( self.covarianceMatrix )

            decorrelatedMatrix = self.Decorrelate( self.covarianceMatrix )

            print '  Found the following decorrelatedMatrix:'
            printMatrix( decorrelatedMatrix )

            print '  Divided by SM cross section:'
            decorrelatedMatrixNormalized = deepcopy( decorrelatedMatrix )
            for i in xrange( self.nTheoryUncertainties ):
                for j in xrange( self.nTheoryUncertainties ):
                    decorrelatedMatrixNormalized[i][j] = 1 + decorrelatedMatrix[i][j] / self.SMXSInsideExperimentalBins[i]
                    # # Skip some uncertainties if the effect is too small
                    # if abs(decorrelatedMatrixNormalized[i][j] - 1) < 0.0005:
                    #     decorrelatedMatrixNormalized[i][j] = 0.
            printMatrix( decorrelatedMatrixNormalized )


            sortedSignals = deepcopy(self.signals)
            sortedSignals.sort( key = lambda signal: float( re.search( r'_[GTLE]*([pm\d]+)', signal ).group(1).replace('p','.').replace('m','-') ))

            for iTheoryUncertainty in xrange(self.nTheoryUncertainties):

                # if self.verbose:
                #     print 'Doing theoryUncertainty_{0}'.format( iTheoryUncertainty )

                errDict = {}
                for binName in self.DC.bins:
                    errDict[binName] = {}
                    for processName in self.DC.processes:

                        errDict[binName][processName] = 0.

                        # for skipBin in self.skipBins:
                        #     if skipBin in processName:
                        #         continue

                        if processName in self.signals:
                            iProcess = self.signals.index(processName)
                            if iProcess < self.nTheoryUncertainties:
                                errDict[binName][processName] = decorrelatedMatrixNormalized[iTheoryUncertainty][iProcess]
                                # if self.verbose:
                                #     print '  Registering errDict[binName=\'{0}\'][processName=\'{1}\'] = '.format( binName, processName )
                                #     print '    {0:+9.4f}'.format( decorrelatedMatrixNormalized[iTheoryUncertainty][iProcess] )
                        
                        
                systematicName = 'theoryUncertainty_{0}'.format(
                    # signals[iTheoryUncertainty] # This is even incorrect I think
                    iTheoryUncertainty
                    )

                self.DC.systs.append(
                    ( systematicName, False, 'lnN', [], errDict )
                    )

                if self.verbose:
                    print 'Added nuisance \'{0}\''.format( systematicName )
                        

    ################################################################################
    # Helper functions
    ################################################################################

    #____________________________________________________________________
    def MakeRatioOfBRsExpressions( self ):
        self.chapter( 'Starting model.MakeRatioOfBRsExpressions()' )

        # Define the two floating pars
        self.modelBuilder.doVar( 'hgg_ratioBRmodifier[1.0,-2.0,4.0]' )
        self.modelBuilder.doVar( 'ratio_BR_hgg_hzz[0.086,0.0,0.5]' )

        # Load spline for HZZ BR (function of MH)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        datadir = os.environ['CMSSW_BASE'] + '/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg'
        self.SMH.textToSpline( 'BR_hzz', os.path.join( datadir, 'sm/br/BR4.txt' ), ycol=11 );
        spline_SMBR_hzz = self.modelBuilder.out.function('BR_hzz')

        # Seems to work:
        # self.modelBuilder.out.var('MH').setVal(125.09)
        # spline_SMBR_hzz.Print()

        # Load spling for Hgg BR
        spline_SMBR_hgg = self.modelBuilder.out.function('fbr_13TeV')


        # ======================================
        # Create 'hzz_BRmodifier', as a function of 'hgg_BRmodifier' and 'ratio_BR_hgg_hzz'

        hzz_modifier_expr = 'expr::{name}("{formula}",{commaSeparatedParameters})'.format(
            name                     = 'hzz_ratioBRmodifier',
            formula                  = '@0*(@1/@2)/@3',
            commaSeparatedParameters = ','.join([
                'hgg_ratioBRmodifier',
                spline_SMBR_hgg.GetName(),
                spline_SMBR_hzz.GetName(),
                'ratio_BR_hgg_hzz'
                ])
            )
        print 'Processing expr:'
        print '    ',hzz_modifier_expr
        self.modelBuilder.factory_( hzz_modifier_expr )


    #____________________________________________________________________
    def MakeTotalXSExpressions( self ):
        self.chapter( 'Starting model.MakeTotalXSExpressions()' )

        self.modelBuilder.doVar( 'r_totalXS[1.0,0.0,3.0]' )

        thingsToSum = []

        for iTheoryBin in xrange( self.nTheoryBins ):

            binWidth = self.theoryBinBoundaries[iTheoryBin+1] - self.theoryBinBoundaries[iTheoryBin]

            self.modelBuilder.factory_(
                'expr::IncXS{0}( "{1}*{2}*@0", parametrization{0} )'.format(
                    iTheoryBin, binWidth, self.SMXS[iTheoryBin]
                    )
                )

            thingsToSum.append( 'IncXS{0}'.format(iTheoryBin) )

        self.modelBuilder.factory_(
            'sum::totalXS( {0} )'.format(
                ','.join(thingsToSum)
                )
            )

        # Get the SM xs by evaluating the parametrization at SM

        # Make sure couplings are at SM
        for coupling in self.couplings:
            self.modelBuilder.out.var(coupling).setVal( self.SMDict['couplings'][coupling] )

        # self.modelBuilder.doVar( 'totalXS_SM[55.70628722]' )
        # Set SM spectrum integral to the SM xs
        self.modelBuilder.doVar( 'totalXS_SM[{0}]'.format(
            self.modelBuilder.out.function('totalXS').getVal()
            ))
        self.modelBuilder.out.var('totalXS_SM').setConstant(True)


        self.modelBuilder.factory_( 'expr::totalXSmodifier( "@0*@1/@2", r_totalXS, totalXS_SM, totalXS )' )


        if self.FitOnlyNormalization:
            # Make also modifier if only fitting normalization
            self.modelBuilder.out.var('r_totalXS').setVal(1.0)
            self.modelBuilder.out.var('r_totalXS').setConstant(True)
            self.modelBuilder.factory_( 'expr::globalTotalXSmodifier( "@0/@1", totalXS, totalXS_SM )' )

            if self.verbose:
                self.modelBuilder.out.function('globalTotalXSmodifier').Print()



    #____________________________________________________________________
    def GetCouplingOrDefineNew( self, couplingName, possibyExistingCouplings=None, defaultValue=1. ):
        # Look for possibly existing couplings that also define the desired coupling

        if possibyExistingCouplings is None:
            possibyExistingCouplings = [ couplingName ]
        if not couplingName in possibyExistingCouplings:
            possibyExistingCouplings.append( couplingName )

        # ======================================
        # Get list of all components in workspace

        componentArgset   = self.modelBuilder.out.components()
        nComponents       = componentArgset.getSize()
        componentIterator = self.modelBuilder.out.componentIterator()
        allComponentsList = []
        for i in xrange(nComponents):
            element = componentIterator.Next()
            allComponentsList.append( element.GetName() )

        # ======================================
        # Check for possibly existing couplings

        for possibyExistingCoupling in possibyExistingCouplings:
            if possibyExistingCoupling in allComponentsList:
                if self.verbose:
                    print 'Found pre-existing coupling \'{0}\', which will be used for \'{1}\''.format( possibyExistingCoupling, couplingName )
                usedCouplingName = possibyExistingCoupling
                break
        else:
            if self.verbose:
                print 'Creating new variable \'{0}\''.format( couplingName )
            self.modelBuilder.doVar( '{0}[{1}]'.format( couplingName, defaultValue ) )
            self.modelBuilder.out.var(couplingName).setConstant(True)
            usedCouplingName = couplingName

        # If the possibly existing variable was actually a function, this will break...
        couplingVar = self.modelBuilder.out.var(usedCouplingName)

        return couplingVar


    #____________________________________________________________________
    def MakeWidthExpressions( self ):
        self.chapter( 'Starting model.MakeWidthExpressions()' )

        # ======================================
        # Prepare some variables in expressions in WS

        # Make sure couplings are defined
        kappa_t   = self.GetCouplingOrDefineNew( 'kappa_t', [ 'kappat', 'ct' ] )
        kappa_b   = self.GetCouplingOrDefineNew( 'kappa_b', [ 'kappab', 'cb' ] )
        kappa_c   = self.GetCouplingOrDefineNew( 'kappa_c', [ 'kappac', 'cc' ] )
        kappa_glu = self.GetCouplingOrDefineNew(
            'kappa_glu',
            [ 'kappaglu', 'cglu', 'cg', 'kappag', 'kappa_g' ],
            defaultValue = 0.0
            )

        kappa_V   = self.GetCouplingOrDefineNew( 'kappa_V' )
        kappa_W   = self.GetCouplingOrDefineNew( 'kappa_V' )
        kappa_Z   = self.GetCouplingOrDefineNew( 'kappa_V' )

        kappa_tau = self.GetCouplingOrDefineNew( 'kappa_tau' )
        kappa_mu  = self.GetCouplingOrDefineNew( 'kappa_mu' )

        if self.verbose:
            print '\nPrintout of couplings:'
            kappa_t.Print()
            kappa_b.Print()
            kappa_c.Print()
            kappa_W.Print()
            kappa_Z.Print()
            kappa_tau.Print()
            kappa_mu.Print()
            print ''

        # Other needed expressions
        self.modelBuilder.factory_("expr::kappa_mu_expr(\"@0*@1+(1-@0)*@2\", CMS_use_kmu[0], kappa_mu, kappa_tau)")

        # Invisible BR, needed for the other expressions but fixing to 0. now
        self.modelBuilder.doVar("BRinv[0.,0.,1.]")
        self.modelBuilder.out.var("BRinv").setConstant(True)


        # ======================================
        # Application of SMHiggsBuilder

        if not hasattr( self, 'SMH' ): self.SMH = SMHiggsBuilder(self.modelBuilder)

        # SM BR's, called 'SM_BR_(decayChannel)'
        for decayChannel in self.SM_HIGG_DECAYS:
            self.SMH.makeBR( decayChannel )

        # BR uncertainties
        # doBRU = False
        if self.DoBRUncertainties:
            self.SMH.makePartialWidthUncertainties()
            if self.verbose:
                print '\nPrintout of HiggsDecayWidth_UncertaintyScaling_{decaychannel}:'
                for decayChannel in self.SM_HIGG_DECAYS:
                    if decayChannel == 'hss': continue
                    print 'HiggsDecayWidth_UncertaintyScaling_{0}:'.format( decayChannel )
                    self.modelBuilder.out.function( 'HiggsDecayWidth_UncertaintyScaling_{0}'.format( decayChannel ) ).Print()

                for otherVar in [
                        'HiggsDecayWidthTHU_hvv', 'HiggsDecayWidthTHU_hgg', 'HiggsDecayWidthTHU_hll', 'HiggsDecayWidthTHU_hqq', 'HiggsDecayWidthTHU_hzg', 'HiggsDecayWidthTHU_hll', 'HiggsDecayWidthTHU_hgluglu',
                        'param_mt', 'param_alphaS', 'param_mB', 'param_mC',
                        ]:
                    try:
                        self.modelBuilder.out.function( otherVar ).Print()
                    except AttributeError:
                        self.modelBuilder.out.var( otherVar ).Print()

        else:
            for decayChannel in self.SM_HIGG_DECAYS: 
                self.modelBuilder.factory_( 'HiggsDecayWidth_UncertaintyScaling_%s[1.0]' % decayChannel )

        # makeScaling functions copied from https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/74x-root6/python/LHCHCGModels.py#L355-L358

        self.modelBuilder.factory_( 'expr::kappa_t_plus_kappa_glu( "@0+@1", {0}, {1} )'.format( kappa_t.GetName(), kappa_glu.GetName() ) )
        self.SMH.makeScaling(
            'hgluglu',
            Cb      = kappa_b.GetName(),
            Ctop    = 'kappa_t_plus_kappa_glu'
            )
        self.SMH.makeScaling(
            'hgg',
            Cb      = kappa_b.GetName(),
            Ctop    = kappa_t.GetName(),
            CW      = kappa_W.GetName(),
            Ctau    = kappa_tau.GetName()
            )
        self.SMH.makeScaling(
            'hzg',
            Cb      = kappa_b.GetName(),
            Ctop    = kappa_t.GetName() ,
            CW      = kappa_W.GetName(),
            Ctau    = kappa_tau.GetName()
            )

        ## partial witdhs, normalized to the SM one
        self.modelBuilder.factory_(
            'expr::c7_Gscal_Z('
            '"@0*@0*@1*@2",'
            ' {0}, SM_BR_hzz, HiggsDecayWidth_UncertaintyScaling_hzz)'.format( kappa_Z.GetName() )
            )
        self.modelBuilder.factory_(
            'expr::c7_Gscal_W('
            '"@0*@0*@1*@2",'
            ' {0}, SM_BR_hww, HiggsDecayWidth_UncertaintyScaling_hww)'.format( kappa_W.GetName() )
            )
        self.modelBuilder.factory_(
            'expr::c7_Gscal_tau('
            '"@0*@0*@1*@4+@2*@2*@3*@5",'
            ' {0}, SM_BR_htt, kappa_mu_expr, SM_BR_hmm, HiggsDecayWidth_UncertaintyScaling_htt, HiggsDecayWidth_UncertaintyScaling_hmm)'.format( kappa_tau.GetName() )
            )
        self.modelBuilder.factory_(
            'expr::c7_Gscal_bottom('
            '"@0*@0 * (@1*@3+@2)",'
            ' {0}, SM_BR_hbb, SM_BR_hss, HiggsDecayWidth_UncertaintyScaling_hbb)'.format( kappa_b.GetName() )
            )
        self.modelBuilder.factory_(
            'expr::c7_Gscal_gluon('
            '"  @0  * @1 * @2",'
            ' Scaling_hgluglu, SM_BR_hgluglu, HiggsDecayWidth_UncertaintyScaling_hgluglu)'
            )
        self.modelBuilder.factory_(
            'expr::c7_Gscal_gamma('
            '"@0*@1*@4 + @2*@3*@5",'
            '  Scaling_hgg, SM_BR_hgg, Scaling_hzg, SM_BR_hzg, HiggsDecayWidth_UncertaintyScaling_hgg, HiggsDecayWidth_UncertaintyScaling_hzg)'
            )
        kappa_c_ForBR = kappa_c if kappa_c.GetName() in self.couplings else kappa_t
        self.modelBuilder.factory_(
            'expr::c7_Gscal_top('
            '"@0*@0 * @1*@2",'
            ' {0}, SM_BR_hcc, HiggsDecayWidth_UncertaintyScaling_hcc)'.format( kappa_c_ForBR.GetName() )
            )


        ## fix to have all BRs add up to unity
        self.modelBuilder.factory_(
            'sum::c7_SMBRs({0})'.format(
                ','.join([ 'SM_BR_{0}'.format(decayChannel) for decayChannel in self.SM_HIGG_DECAYS ])
                # "SM_BR_"+X for X in "hzz hww htt hmm hcc hbb hss hgluglu hgg hzg".split())
                ))
        self.modelBuilder.out.function("c7_SMBRs").Print("")      

        ## total witdh, normalized to the SM one
        self.modelBuilder.factory_(
            'expr::c7_Gscal_tot('
            '"(@1+@2+@3+@4+@5+@6+@7)/@8/(1-@0)", BRinv, c7_Gscal_Z, c7_Gscal_W, c7_Gscal_tau, c7_Gscal_top, c7_Gscal_bottom, c7_Gscal_gluon, c7_Gscal_gamma, c7_SMBRs)'
            )

        ## BRs, normalized to the SM ones: they scale as (partial/partial_SM) / (total/total_SM) 
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hww('
            '"@0*@0*@2/@1", {0}, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hww)'.format( kappa_W.GetName() )
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hzz('
            '"@0*@0*@2/@1", {0}, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hzz)'.format( kappa_Z.GetName() )
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_htt('
            '"@0*@0*@2/@1", {0}, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_htt)'.format( kappa_tau.GetName() )
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hmm('
            '"@0*@0*@2/@1", kappa_mu_expr, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hmm)'
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hbb('
            '"@0*@0*@2/@1", {0}, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hbb)'.format( kappa_b.GetName() )
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hcc('
            '"@0*@0*@2/@1", {0}, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hcc)'.format( kappa_c_ForBR.GetName() )
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hgg('
            '"@0*@2/@1", Scaling_hgg, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hgg)'
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hzg('
            '"@0*@2/@1", Scaling_hzg, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hzg)'
            )
        self.modelBuilder.factory_(
            'expr::c7_BRscal_hgluglu('
            '"@0*@2/@1", Scaling_hgluglu, c7_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hgluglu)'
            )

        if self.verbose:
            print '\nPrintout of c7_Gscal_{decayChannel}:'
            for Gscal in [
                    'c7_Gscal_Z',
                    'c7_Gscal_W',
                    'c7_Gscal_tau',
                    'c7_Gscal_bottom',
                    'c7_Gscal_gluon',
                    'c7_Gscal_gamma',
                    'c7_Gscal_top',
                    'c7_Gscal_tot',
                    ]:
                self.modelBuilder.out.function(Gscal).Print()
            print '\nPrintout of c7_BRscal_{decayChannel}:'
            for BRscal in [
                    'c7_BRscal_hww',
                    'c7_BRscal_hzz',
                    'c7_BRscal_htt',
                    'c7_BRscal_hmm',
                    'c7_BRscal_hbb',
                    'c7_BRscal_hcc',
                    'c7_BRscal_hgg',
                    'c7_BRscal_hzg',
                    'c7_BRscal_hgluglu',
                    ]:
                self.modelBuilder.out.function(BRscal).Print()


        # ======================================
        # The final scaling parameters used will be called 'hggBRmodifier' and 'hzzBRmodifier'

        self.modelBuilder.factory_(
            'expr::hggBRmodifier('
            '"@0", c7_BRscal_hgg )'
            )

        self.modelBuilder.factory_(
            'expr::hzzBRmodifier('
            '"@0", c7_BRscal_hzz )'
            )

        # xH modifier
        self.modelBuilder.factory_(
            'expr::xH_modifier('
            '"@0*@0", {0} )'.format( kappa_V.GetName() )
            )


        # self.modelBuilder.factory_(
        #     'expr::hggBRmodifier_times_xH_modifier('
        #     '"@0*@1", hggBRmodifier, xH_modifier )'
        #     )
        # self.modelBuilder.factory_(
        #     'expr::hzzBRmodifier_times_xH_modifier('
        #     '"@0*@1", hzzBRmodifier, xH_modifier )'
        #     )


        # if self.MakeLumiScalable:

        #     self.modelBuilder.factory_(
        #         'expr::lumiScale_times_hggBRmodifier('
        #         '"@0*@1", lumiScale, hggBRmodifier )'
        #         )

        #     self.modelBuilder.factory_(
        #         'expr::lumiScale_times_hzzBRmodifier('
        #         '"@0*@1", lumiScale, hzzBRmodifier )'
        #         )

        #     # xH
        #     self.modelBuilder.factory_(
        #         'expr::lumiScale_times_hggBRmodifier_times_xH_modifier('
        #         '"@0*@1*@2", lumiScale, hggBRmodifier, xH_modifier )'
        #         )
        #     self.modelBuilder.factory_(
        #         'expr::lumiScale_times_hzzBRmodifier_times_xH_modifier('
        #         '"@0*@1*@2", lumiScale, hzzBRmodifier, xH_modifier )'
        #         )

        self.modelBuilder.out.defineSet( 'BRvariables', ','.join([
            kappa_t.GetName(),
            kappa_b.GetName(),
            kappa_c.GetName(),
            kappa_W.GetName(),
            kappa_Z.GetName(),
            kappa_tau.GetName(),
            kappa_mu.GetName(),
            'hggBRmodifier',
            'hzzBRmodifier',
            'xH_modifier',
            'Scaling_hgg',
            'Scaling_hzg',
            'Scaling_hgluglu',
            ] ) )


    #____________________________________________________________________
    def chapter( self, text, indent=0 ):
        if self.verbose:
            print '\n{tabs}{line}\n{tabs}{text}\n'.format(
                tabs = '    ' * indent,
                line = '-' * 70,
                text = text
                )


    #____________________________________________________________________
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


    #____________________________________________________________________
    def readCorrelationMatrixFile( self, correlationMatrixFile ):
        with open( correlationMatrixFile, 'r' ) as correlationMatrixFp:
            lines = [ l.strip() for l in correlationMatrixFp.readlines() if len(l.strip()) > 0 and not l.strip().startswith('#') ]

        corrMat = [ [ float(number) for number in line.split() ] for line in lines ]

        # Check if it is square
        if not all([ len(row) == len(corrMat) for row in corrMat ]):
            print corrMat
            raise CouplingModelError( '[ERROR] inputted matrix is not square - Found ^ ' )

        return corrMat


    #____________________________________________________________________
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
                raise CouplingModelError(
                    '[ERROR] Found {0} elements on line in \'{1}\''.format( len(line), errorFile )
                    )
                return
            
        return symmErrors


    #____________________________________________________________________
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


    #____________________________________________________________________
    def BinIsHgg( self, bin ):

        # If the input is hardcoded to be only hzz or hgg:
        if self.isOnlyHgg:
            return True
        if self.isOnlyHZZ:
            return False
        
        # This is not very robust coding
        if bin.startswith('hgg_'):
            # if self.verbose: print '    Bin \'{0}\' is classified as a hgg bin!'.format( bin )
            return True
        elif bin.startswith('hzz_'):
            # if self.verbose: print '    Bin \'{0}\' is classified as a hzz bin!'.format( bin )
            return False
        else:
            raise CouplingModelError( 'Bin \'{0}\' can not be classified as either a hgg or hzz bin'.format(bin) )



couplingModel=CouplingModel()

