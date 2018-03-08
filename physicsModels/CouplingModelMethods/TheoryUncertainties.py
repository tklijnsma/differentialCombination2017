from physicsModels.MethodHandler import flag_as_method

import os, sys, numpy, itertools, re
import ROOT
from math import sqrt
from copy import deepcopy


@flag_as_method
def addTheoryUncertaintyNuisances( self ):

    if self.verbose:
        print '\n' + '='*80
        print 'Inserting theory uncertainties\n'

    # First check if given input makes sense
    doDecorrelation = False
    if (
        not self.theoryUncertaintiesPassed
        and not self.correlationMatrixPassed
        and not self.covarianceMatrixPassed
        ):
        print 'No theory uncertainties are specified; Running without theory uncertainties'
        return

    elif (
        self.theoryUncertaintiesPassed
        and self.correlationMatrixPassed
        and not self.covarianceMatrixPassed
        ):
        if not len(self.theoryUncertainties) == len(self.correlationMatrix):
            print '[ERROR] Cannot build a covariance matrix out of given input'
            print '        len(self.theoryUncertainties) = {0}'.format( len(self.theoryUncertainties) )
            print '        len(self.correlationMatrix)   = {0}'.format( len(self.correlationMatrix) )
            raise self.CouplingModelError()

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
        raise self.CouplingModelError()


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
                    

@flag_as_method
def readTheoryFile(self, theoryFile):
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

@flag_as_method
def readCorrelationMatrixFile( self, correlationMatrixFile ):
    with open( correlationMatrixFile, 'r' ) as correlationMatrixFp:
        lines = [ l.strip() for l in correlationMatrixFp.readlines() if len(l.strip()) > 0 and not l.strip().startswith('#') ]

    corrMat = [ [ float(number) for number in line.split() ] for line in lines ]

    # Check if it is square
    if not all([ len(row) == len(corrMat) for row in corrMat ]):
        print corrMat
        raise self.CouplingModelError( '[ERROR] inputted matrix is not square - Found ^ ' )

    # N = len(corrMat)
    # print '[WARNING] Adding 1e-12 to the diagonal of the correlation matrix in order to protect against non-positive-definiteness'
    # for i in xrange(N):
    #     corrMat[i][i] += 1e-12

    return corrMat

@flag_as_method
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
            raise self.CouplingModelError(
                '[ERROR] Found {0} elements on line in \'{1}\''.format( len(line), errorFile )
                )
            return
        
    return symmErrors

@flag_as_method
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

    print '[Decorrelating]: Found Eigenvalues: {0}'.format([eigenValues(i) for  i in xrange(N)])
    print '[Decorrelating]: Found Eigenvectors:'
    printMatrix([[eigenVectors(i,j) for j in xrange(N)] for i in xrange(N)])

    decorrelatedMatrix = [ [ 999 for j in xrange(N) ] for i in xrange(N) ]
    for i in xrange(N):
        for j in xrange(N):
            eigen_value = eigenValues(j)
            if abs(eigen_value) < 1e-9:
                print '[WARNING] Found abs(eigen_value) = {0} < 1e-6; rounding to 0.0'.format(eigen_value)
                eigen_value = 0.0
            decorrelatedMatrix[i][j] = eigenVectors(i,j) * sqrt(eigen_value)

    return decorrelatedMatrix


def printMatrix( matrix, indent = '    ', scientific=True ):
    for row in matrix:
        print indent + ' '.join([ '{0:<+10.2{1}}'.format( number, 'E' if scientific else 'f' ) for number in row ])
