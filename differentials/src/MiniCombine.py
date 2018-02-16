#!/usr/bin/env python
"""
Thomas Klijnsma
"""

import os, sys, re, numpy
from os.path import *
from glob import glob
from array import array

import Commands
import TheoryCommands
import PhysicsCommands
from Container import Container
from Parametrization import Parametrization
import PlotCommands

import scipy
scipyversion = scipy.__version__.split('.')[1]
if int(scipyversion) < 15:
    Commands.ThrowError( 'Need the python environment for MiniCombine!' )
from scipy.optimize import minimize

import ROOT


class MiniCombine(object):
    """docstring for MiniCombine"""

    TITLES = {
        'ct' : '#kappa_{t}',
        'cg' : '#kappa_{g}',
        'kappac' : '#kappa_{c}',
        'kappab' : '#kappa_{b}',
        'kappat' : '#kappa_{t}',
        }

    #____________________________________________________________________
    def __init__(
            self,
            couplings = [ 'kappab', 'kappac' ],
            lastBinIsOverflow = False,
            ):

        self.couplings = couplings
        self.SMIsSet = False
        self.verbose = True
        self.lastBinIsOverflow = lastBinIsOverflow

        self.extra_chi2_constraints = []

        self.expBinBoundaries = []

        self.name = 'MiniCombine'

    #____________________________________________________________________
    def SetSM(
            self,
            SM
            ):
        self.SM = SM
        self.SMIsSet = True

    #____________________________________________________________________
    def SetCombinationResult(
            self,
            scanDir,
            excludeLastBin=False,
            ):
        self.scanDir = scanDir
        self.signalStrengthNames = GetPOIsFromScandirHeuristic(scanDir)
        self.expBinBoundaries = PhysicsCommands.FigureOutBinning(self.signalStrengthNames)
        if self.lastBinIsOverflow: self.expBinBoundaries[-1] = 10000.

        if excludeLastBin:
            self.signalStrengthNames = self.signalStrengthNames[:-1]
            self.expBinBoundaries = self.expBinBoundaries[:-1]

        self.scans = []
        rawOutputs = PhysicsCommands.GetScanResults(
            self.signalStrengthNames,
            self.scanDir,
            )
        for rawOutput in rawOutputs:
            scan = Container()
            scan.mus, scan.deltaNLLs = PhysicsCommands.FilterScan(rawOutput)
            scan.unc = PhysicsCommands.FindMinimaAndErrors( scan.mus, scan.deltaNLLs )
            self.scans.append(scan)

        self.mus_data      = [ s.unc['min'] for s in self.scans ]
        self.deltamus_data = [ s.unc['symmError'] for s in self.scans ]

        print '\nSet combination result:'
        print '  expBinBoundaries =', self.expBinBoundaries
        print '  mus_data         =', self.mus_data
        print '  deltamus_data    =', self.deltamus_data


    #____________________________________________________________________
    def SetCorrelationMatrix(
            self,
            corrMatFile,
            forceDiagonalMatrix = False,
            ):

        corrMat = []
        corrMatFp = ROOT.TFile.Open(corrMatFile)

        if forceDiagonalMatrix:
            Commands.Warning( 'Overwriting with square matrix for testing purposes:' )
            corrMat = [ [ 1.0 if i==j else 0.0 for j in xrange(len(self.signalStrengthNames)) ] for i in xrange(len(self.signalStrengthNames)) ]

        else:
            with Commands.OpenRootFile(corrMatFile) as corrMatFp:
                fit = corrMatFp.Get('fit')
                for poi1 in self.signalStrengthNames:
                    corrMatRow = []
                    for poi2 in self.signalStrengthNames:
                        corrMatRow.append( fit.correlation( poi1, poi2 ) )
                    corrMat.append( corrMatRow )

        self.corrMat = corrMat
        print '\nFound corrMat from {0}:'.format(corrMatFile)
        print numpy.array(corrMat)


    #____________________________________________________________________
    def Parametrize(
            self,
            derivedTheoryContainers,
            includeLinearTerms
            ):

        self.includeLinearTerms = includeLinearTerms
        self.derivedTheoryContainers = derivedTheoryContainers

        self.parametrization = Parametrization( verbose=False )
        # self.parametrization.SetSM( self.SM )
        self.parametrization.Parametrize(
            derivedTheoryContainers,
            includeLinearTerms = self.includeLinearTerms,
            couplingsToParametrize = self.couplings
            )


    #____________________________________________________________________
    def Link(
            self
            ):
        """Link together the given SM, scanDir, corrMat and parametrization"""

        self.nBinsExp = len(self.expBinBoundaries)-1

        # Attempt to calculate cross sections
        if not self.SMIsSet:
            Commands.ThrowError( 'Need a SM container supplied' )

        self.SMintegralFn = TheoryCommands.GetIntegral(
            self.SM.binBoundaries,
            self.SM.crosssection
            )
        self.XS_SM_tot = self.SMintegralFn( self.expBinBoundaries[0], self.expBinBoundaries[-1] )

        self.XSs_SM = []
        for i in xrange(self.nBinsExp):
            if self.lastBinIsOverflow and i == self.nBinsExp-1:
                self.XSs_SM.append(
                    self.SMintegralFn( self.expBinBoundaries[i], self.expBinBoundaries[i+1],
                        verbose=False
                        )
                    )
            else:
                self.XSs_SM.append(
                    self.SMintegralFn( self.expBinBoundaries[i], self.expBinBoundaries[i+1] )
                      / ( self.expBinBoundaries[i+1] - self.expBinBoundaries[i] )
                    )

        print '  XSs_SM      =', self.XSs_SM

        # expIntegral = TheoryCommands.GetIntegral( self.expBinBoundaries, self.XSs_SM )
        # print '  int(XSs_SM) =', expIntegral( 0., 500. ), ' (last bin not properly done)'

        # Use SMXSs to compute data XSs
        self.XSs_data      = [ mu * SMXS for mu, SMXS in zip( self.mus_data,      self.XSs_SM ) ]
        self.deltaXSs_data = [ mu * SMXS for mu, SMXS in zip( self.deltamus_data, self.XSs_SM ) ]

        #  --> ONLY IF DIRECTLY USING XSs IN THE CHI2!!! Which is not often!!
        # # Compute covariance using deltaXSs_data
        # self.covMat = [ [ 0. for i in xrange(self.nBinsExp) ] for j in xrange(self.nBinsExp) ]
        # for i in xrange(self.nBinsExp):
        #     for j in xrange(self.nBinsExp):
        #         self.covMat[i][j] = self.corrMat[i][j] * self.deltaXSs_data[i] * self.deltaXSs_data[j]
        # self.covMat_inversed_npArray = numpy.linalg.inv( numpy.array(self.covMat) )

        # Compute covariance using deltamus_data
        self.covMat = [ [ 0. for i in xrange(self.nBinsExp) ] for j in xrange(self.nBinsExp) ]
        for i in xrange(self.nBinsExp):
            for j in xrange(self.nBinsExp):
                self.covMat[i][j] = self.corrMat[i][j] * self.deltamus_data[i] * self.deltamus_data[j]
        self.covMat_inversed_npArray = numpy.linalg.inv( numpy.array(self.covMat) )


    #____________________________________________________________________
    def AddMuConstraint( self, mu, unc, left, right ):

        print '\nPlugging in an uncorrelated constraint mu_{{{0}-{1}}} = {2} +- {3}'.format(
            left, right, mu, unc
            )


        def constr(self):

            XS_param = self.parametrization.EvaluateForBinning(
                self.SM.binBoundaries, [ left, right ],
                returnRatios = False,
                verbose = self.verbose,
                lastBinIsOverflow = False, # Divide by bin width
                )[0]

            XS_SM = self.SMintegralFn( left, right ) / ( right - left )
            mu_param = XS_param / XS_SM

            delta_chi2 = ( mu - mu_param )**2 / (unc**2)

            # print '  Modifying chi2 between {0:.1f} - {1:.1f} using mu = {2:.2f} +- {3:.2f}: mu_param = {4:.2f}, delta_chi2 = {5:.2f}'.format(
            #     left, right,
            #     mu, unc,
            #     mu_param, delta_chi2
            #     )

            return delta_chi2

        constr.suffix = '_{0}to{1}_unc{2}'.format(
            Commands.ConvertFloatToStr(left), Commands.ConvertFloatToStr(right), Commands.ConvertFloatToStr(unc, nDecimals=2 )
            )

        self.extra_chi2_constraints.append( constr )


    #____________________________________________________________________
    def AddAsymMuConstraint( self, mu, unc_up, unc_down, left, right ):
        print '\nPlugging in an uncorrelated asymmetric constraint mu_{{{0}-{1}}} = {2}  {3:+.2f}/{4:+.2f}'.format( left, right, mu, unc_up, unc_down )

        def constr(self):

            XS_param = self.parametrization.EvaluateForBinning(
                self.SM.binBoundaries, [ left, right ],
                returnRatios = False,
                verbose = self.verbose,
                lastBinIsOverflow = False, # Divide by bin width
                )[0]

            XS_SM = self.SMintegralFn( left, right ) / ( right - left )
            mu_param = XS_param / XS_SM

            if mu > mu_param:
                delta_chi2 = ( mu - mu_param )**2 / (unc_up**2)
            else:
                delta_chi2 = ( mu - mu_param )**2 / (unc_down**2)

            return delta_chi2

        constr.suffix = '_{0}to{1}_uncup{2}_uncdown{3}'.format(
            Commands.ConvertFloatToStr(left),
            Commands.ConvertFloatToStr(right),
            Commands.ConvertFloatToStr(unc_up, nDecimals=2 ),
            Commands.ConvertFloatToStr(unc_down, nDecimals=2 )
            )

        self.extra_chi2_constraints.append( constr )


    #____________________________________________________________________
    def AddOverflowMuConstraint( self, mu, unc, left ):
        right = self.SM.binBoundaries[-1]

        print '\nPlugging in an uncorrelated constraint mu_{{{0}-{1}}} = {2} +- {3}'.format(
            left, right, mu, unc
            )

        def constr(self):

            XS_param = self.parametrization.EvaluateForBinning(
                self.SM.binBoundaries, [ left, right ],
                returnRatios = False,
                verbose = self.verbose,
                lastBinIsOverflow = True, # Don't divide by bin width
                )[0]

            XS_data = self.SMintegralFn( left, right )

            mu_param = XS_param / XS_data

            delta_chi2 = ( mu - mu_param )**2 / (unc**2)

            return delta_chi2

        constr.suffix = '_GT{0}_unc{1}'.format(
            Commands.ConvertFloatToStr(left), Commands.ConvertFloatToStr( unc, nDecimals=2 )
            )

        self.extra_chi2_constraints.append( constr )


    #____________________________________________________________________
    def AddAsymMuOverflowConstraint( self, mu, unc_up, unc_down, left ):
        right = self.SM.binBoundaries[-1]
        print '\nPlugging in an uncorrelated asymmetric constraint mu_{{{0}-{1}}} = {2}  {3:+.2f}/{4:+.2f}'.format( left, right, mu, unc_up, unc_down )

        def constr(self):

            XS_param = self.parametrization.EvaluateForBinning(
                self.SM.binBoundaries, [ left, right ],
                returnRatios = False,
                verbose = self.verbose,
                lastBinIsOverflow = True, # Don't divide by bin width
                )[0]

            XS_SM = self.SMintegralFn( left, right )
            mu_param = XS_param / XS_SM

            if mu > mu_param:
                delta_chi2 = ( mu - mu_param )**2 / (unc_up**2)
            else:
                delta_chi2 = ( mu - mu_param )**2 / (unc_down**2)

            return delta_chi2

        constr.suffix = '_{0}to{1}_uncup{2}_uncdown{3}'.format(
            Commands.ConvertFloatToStr(left),
            Commands.ConvertFloatToStr(right),
            Commands.ConvertFloatToStr(unc_up, nDecimals=2 ),
            Commands.ConvertFloatToStr(unc_down, nDecimals=2 )
            )

        self.extra_chi2_constraints.append( constr )


    #____________________________________________________________________
    def ClearConstraints( self ):
        del self.chi2_function
        del self.extra_chi2_constraints
        self.extra_chi2_constraints = []


    #____________________________________________________________________
    def BuildChi2Function(
            self
            ):

        def new_chi2_function( inputTuple, verbose=False ):

            for coupling, couplingVal in zip( self.couplings, inputTuple ):
                if verbose: print '  Setting {0} to {1} in parametrization'.format( coupling, couplingVal )
                setattr( self.parametrization, coupling, couplingVal )

            XSs_parametrization = self.parametrization.EvaluateForBinning(
                self.SM.binBoundaries, self.expBinBoundaries,
                returnRatios = False,
                verbose = verbose,
                lastBinIsOverflow = self.lastBinIsOverflow,
                )

            mus_parametrization = [ XS / SMXS for XS, SMXS in zip( XSs_parametrization, self.XSs_SM ) ]

            if verbose: print '  Found XSs_parametrization =', XSs_parametrization
            if verbose: print '  Found mus_parametrization =', mus_parametrization

            x_column = numpy.array( [ mu_param - mu_data for mu_param, mu_data in zip( mus_parametrization, self.mus_data ) ] ).T
            chi2 = x_column.T.dot(  self.covMat_inversed_npArray.dot( x_column )  )

            for chi2_constraint in self.extra_chi2_constraints:
                chi2 += chi2_constraint(self)

            return chi2

        self.chi2_function = new_chi2_function


    #____________________________________________________________________
    def ScanGrid(
            self,
            nPointsX, xMin, xMax,
            nPointsY, yMin, yMax
            ):

        # First do a bestfit
        SM_tuple = [ float(getattr( self.SM, coupling )) for coupling in self.couplings ]

        print '\nEvaluating chi2 at SM:'
        print '  {0:10} = {1}'.format( self.couplings[0], SM_tuple[0] )
        print '  {0:10} = {1}'.format( self.couplings[1], SM_tuple[1] )
        print '  {0:10} = {1}'.format( 'chi2', self.chi2_function( SM_tuple, verbose=True ) )

        bestfit = minimize( self.chi2_function, SM_tuple, method='Nelder-Mead', tol=1e-6 )
        bestfit_chi2 = bestfit.fun

        print '\nFound bestfit:'
        print '  {0:10} = {1}'.format( self.couplings[0], bestfit.x[0] )
        print '  {0:10} = {1}'.format( self.couplings[1], bestfit.x[1] )
        print '  {0:10} = {1}'.format( 'chi2', bestfit_chi2 )


        xBinBoundaries = [ xMin + i*(xMax-xMin)/float(nPointsX) for i in xrange(nPointsX+1) ]
        yBinBoundaries = [ yMin + i*(yMax-yMin)/float(nPointsY) for i in xrange(nPointsY+1) ]
        xBinCenters = [ 0.5*( xBinBoundaries[i] + xBinBoundaries[i+1] ) for i in xrange(nPointsX) ]
        yBinCenters = [ 0.5*( yBinBoundaries[i] + yBinBoundaries[i+1] ) for i in xrange(nPointsY) ]

        H2 = ROOT.TH2F(
            'H2', '',
            nPointsX, array( 'f', xBinBoundaries ),
            nPointsY, array( 'f', yBinBoundaries ),
            )
        ROOT.SetOwnership( H2, False )

        for ix, x in enumerate(xBinCenters):
            for iy, y in enumerate(yBinCenters):

                # if ix % 10 == 0 and iy % 10 == 0:
                #     print '\nEvaluating chi2 for ix = {0}, iy = {1}'.format( ix, iy )
                #     print '  x = {0}, y = {1}, chi2_function = {2}, delta_chi2 = {3}'.format(
                #         x, y, self.chi2_function([ x, y ]), self.chi2_function([ x, y ]) - bestfit_chi2
                #         )

                delta_chi2 = self.chi2_function(( x, y )) - bestfit_chi2

                if delta_chi2 < 0.:
                    Commands.Warning( 'deltachi2 for x = {0}, y = {1} is negative: {2}; setting to 0.0'.format( x, y, delta_chi2 ) )
                    delta_chi2 = 0.
                if delta_chi2 > 1000.: delta_chi2 = 1000.

                H2.SetBinContent( ix+1, iy+1, delta_chi2 )

        bestfitPoint = ROOT.TGraph( 1, array( 'f', [bestfit.x[0]] ), array( 'f', [bestfit.x[1]] ) )
        bestfitPoint.SetName('bestfitPoint')
        bestfitPoint.SetMarkerStyle(34)
        bestfitPoint.SetMarkerSize(0.9)
        ROOT.SetOwnership( bestfitPoint, False )

        plotContainer = Container()
        plotContainer.name            = '_'.join(self.couplings)
        plotContainer.H2              = H2
        plotContainer.contours_1sigma = TheoryCommands.GetContoursFromTH2( H2, 2.30 )
        plotContainer.contours_2sigma = TheoryCommands.GetContoursFromTH2( H2, 6.18 )
        plotContainer.bestfitPoint    = bestfitPoint

        for constr in self.extra_chi2_constraints:
            plotContainer.name += constr.suffix


        PlotCommands.PlotSingle2DHistogram(
            plotContainer,
            xMin, xMax,
            yMin, yMax,
            xTitle = self.TITLES.get( self.couplings[0], self.couplings[0] ),
            yTitle = self.TITLES.get( self.couplings[1], self.couplings[1] ),
            plotname = self.name,
            # doPNG=True, doROOT=True
            )


#____________________________________________________________________
def GetPOIsFromScandirHeuristic(
        scanDir
        ):

    shFiles = glob( scanDir + '/*.sh' )

    # Strip the _\d.sh ending
    shFiles = [ basename(f).replace('/','') for f in shFiles if f.endswith('_0.sh') ]

    # Look for the string between job_SCAN_..._Mon00...
    POIs = []
    for f in shFiles:
        match = re.match( r'job_SCAN_(.*?)_[A-Z][a-z][a-z]\d\d', f )
        if not match:
            print 'Could not process {0}; continuing'.format( f )
            continue
        POI = match.group(1)
        POIs.append( POI )

    POIs = list(set(POIs))
    POIs.sort( key = Commands.POIsorter )

    return POIs


