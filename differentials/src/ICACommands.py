#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, numpy
from numpy import std

from os.path import *
from scipy import linalg
from array import array
from math import exp, log, sqrt
from copy import deepcopy
from glob import glob
from shutil import copyfile

from time import strftime
datestr = strftime( '%b%d' )

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
ROOT.gStyle.SetOptStat(0)
c = ROOT.TCanvas( 'c', 'c', 1000, 800 )

CLeftMargin   = 0.15
CRightMargin  = 0.03
CBottomMargin = 0.15
CTopMargin    = 0.03
def set_cmargins(
    LeftMargin   = CLeftMargin,
    RightMargin  = CRightMargin,
    BottomMargin = CBottomMargin,
    TopMargin    = CTopMargin,
    ):
    c.SetLeftMargin( LeftMargin )
    c.SetRightMargin( RightMargin )
    c.SetBottomMargin( BottomMargin )
    c.SetTopMargin( TopMargin )

def get_plot_base(
    xMin = 0, xMax = 1,
    yMin = 0, yMax = 1,
    xTitle = 'x', yTitle = 'y',
    SetTitleSizes = True,
    ):

    base = ROOT.TH1F()
    ROOT.SetOwnership( base, False )
    base.SetName( 'genericbase' )
    base.GetXaxis().SetLimits( xMin, xMax )
    base.SetMinimum( yMin )
    base.SetMaximum( yMax )
    base.SetMarkerColor(0)
    base.GetXaxis().SetTitle( xTitle )
    base.GetYaxis().SetTitle( yTitle )

    if SetTitleSizes:
        base.GetXaxis().SetTitleSize( 0.06 )
        base.GetYaxis().SetTitleSize( 0.06 )

    return base

CPlotDir = 'plots_{0}'.format(datestr)
def save_c( outname, asPDF=True, asPNG=False, asROOT=False ):
    if not isdir(CPlotDir): os.makedirs( CPlotDir )
    outname = join( CPlotDir, outname.replace('.pdf','').replace('.png','').replace('.root','') )
    if asPDF: c.SaveAs( outname + '.pdf' )
    if asPNG: c.SaveAs( outname + '.png' )
    if asROOT: c.SaveAs( outname + '.proot' )


########################################
# Main
########################################

def main():
    test_methods()





def test_methods(
        doWhitening = True,
        doICA       = True,
        doFastICA   = True,
        verbose     = True
        ):

    xPoints, yPoints = generate_data()

    data = numpy.array(
        [ [ x, y ] for x,y in zip( xPoints, yPoints ) ]
        )

    w_init = get_winit( 2, data.dtype )


    if doICA:
        chap( 'Running FastICA_byThomas' )
        res = fast_ica_by_thomas( deepcopy(data), n_components=2, w_init=deepcopy(w_init) )
        res.name = 'ICAbyThomas'
        draw_data( res, verbose )

    if doFastICA:
        chap( 'Running FastICA from scikit package' )

        from sklearn.decomposition import FastICA, PCA

        # Compute ICA
        ica = FastICA( n_components=2, w_init=deepcopy(w_init) )
        S_ = ica.fit_transform( deepcopy(data) )  # Reconstruct signals
        A_ = ica.mixing_  # Get estimated mixing matrix

        FastICA_res = Container( name = 'FastICA' )
        FastICA_res.xPoints = xPoints
        FastICA_res.yPoints = yPoints
        FastICA_res.xPointsRotated = list(S_[:,0])
        FastICA_res.yPointsRotated = list(S_[:,1])

        data_mappedBack = numpy.array( numpy.dot(S_, A_.T) + ica.mean_ )
        FastICA_res.xPointsRotatedBack = list(data_mappedBack[:,0])
        FastICA_res.yPointsRotatedBack = list(data_mappedBack[:,1])

        FastICA_res.rotationMatrix = A_

        draw_data( FastICA_res, verbose )

    if doWhitening:
        chap( 'Running a basic whitening procedure' )
        whitening_res = whitening( deepcopy(data), verbose=True )
        draw_data( whitening_res, verbose )


    # ======================================
    # Make another subdir with all the rotated point plots

    rotatedDataFiles = glob( CPlotDir + '/*_rotateddata.pdf' )

    if len(rotatedDataFiles) > 0:
        subdir = join( CPlotDir, 'copiesOfRotatedPoints' )
        if not isdir( subdir ): os.makedirs( subdir )
        for src in rotatedDataFiles:
            dst = join( subdir, basename(src) )
            copyfile( src, dst )


########################################
# Only a whitening operation
########################################

def whitening( data, verbose=False ):

    xPoints = list(data[:,0])
    yPoints = list(data[:,1])

    # ======================================
    # Calculate basic quantities

    covMatMeasured = numpy.cov( xPoints, yPoints )

    if verbose:
        print '\nMeasured covariance matrix:'
        print covMatMeasured
        print ''

    xStdMeasured = sqrt(covMatMeasured[0][0])
    yStdMeasured = sqrt(covMatMeasured[1][1])

    xMeanMeasured = mean( xPoints )
    yMeanMeasured = mean( yPoints )

    corr = covMatMeasured[0][1] / ( xStdMeasured * yStdMeasured )


    # ======================================
    # Do the whitening

    eigenValues, eigenVectors = numpy.linalg.eig(
        numpy.array( covMatMeasured )
        )

    # Convert to ordinary lists
    eigenValues = list( eigenValues )
    eigenVectors = [ list(eigenVectors[:,i]) for i in xrange(eigenVectors.shape[1]) ]

    # Make eigenvectors unit length (should already be the case but to be sure)
    eigenVectors = [ normalize(eigenVector) for eigenVector in eigenVectors ]

    # Sort according to eigenvalues
    eigenValues, eigenVectors = ( list(t) for t in zip(*sorted( zip(eigenValues, eigenVectors), reverse=True )) )

    whiteningMatrix = get_whitening_matrix( eigenValues, eigenVectors )

    if verbose:
        print '\neigenValues:'
        print eigenValues

        print '\neigenVectors:'
        for eigenVector in eigenVectors:
            print numpy_column( eigenVector )

        print '\nwhiteningMatrix:'
        print whiteningMatrix

        print '\nwhiteningMatrix.T:'
        print whiteningMatrix.T

        print '\nnumpy.linalg.inv(whiteningMatrix):'
        print numpy.linalg.inv(whiteningMatrix)


    xPointsWhitened, yPointsWhitened = apply_matrix_on_data(
        whiteningMatrix, xPoints, yPoints
        )

    ret = Container( name = 'onlyWhitened' )
    ret.xPoints = xPoints
    ret.yPoints = yPoints
    ret.xPointsRotated = xPointsWhitened
    ret.yPointsRotated = yPointsWhitened

    # This is still broken
    ret.xPointsRotatedBack, ret.yPointsRotatedBack = apply_matrix_on_data(
        numpy.linalg.inv(whiteningMatrix),
        xPointsWhitened,
        yPointsWhitened,
        whiten=False,
        # verbose=True
        )
    ret.xPointsRotatedBack = [ x + xMeanMeasured for x in ret.xPointsRotatedBack ]
    ret.yPointsRotatedBack = [ y + yMeanMeasured for y in ret.yPointsRotatedBack ]

    ret.rotationMatrix = whiteningMatrix

    return ret



def get_whitening_matrix( eigenValues, eigenVectors ):
    # sqrtD = numpy.diag( [ sqrt(val) for val in eigenValues ] )
    # ET    = numpy.vstack( eigenVectors ) # Works since we need the transpose
    # V = numpy.dot( sqrtD, ET )

    oneOverSqrtD = numpy.diag( [ 1./sqrt(val) for val in eigenValues ] )
    E            = numpy.hstack( [ numpy_column(vec) for vec in eigenVectors ] )
    ET           = E.T

    #   E * D^-0.5 * E^T
    V = E.dot( oneOverSqrtD.dot( ET ) )

    return V


def apply_matrix_on_data( matrix, xPoints, yPoints, whiten=True, verbose=False ):
    
    if whiten:
        xMean = mean(xPoints)
        yMean = mean(yPoints)
        xStd  = std(xPoints)
        yStd  = std(yPoints)

    xPointsWhitened = []
    yPointsWhitened = []

    for x, y in zip( xPoints, yPoints ):
        if whiten:
            x = ( x - xMean ) # / xStd
            y = ( y - yMean ) # / yStd
        whitenedPoint = matrix.dot( numpy_column([x, y]) )
        xPointsWhitened.append( whitenedPoint[0][0] )
        yPointsWhitened.append( whitenedPoint[1][0] )

        if verbose:
            print '-'*70
            print 'matrix = ', matrix
            print 'dot vector = ', numpy_column([x, y])
            print 'outcome = ', whitenedPoint

        # print whitenedPoint

    return xPointsWhitened, yPointsWhitened


########################################
# Minimal FastICA implementation
########################################

def fast_ica_by_thomas( data, n_components, verbose=False, w_init=None ):

    # The hidden transposer
    data = data.T

    # ======================================
    # Whiten data

    whiteningMatrix, X_whitened = fast_ica_whiten( data, n_components )

    if verbose:
        print '\nX_whitened:'
        print X_whitened
        print '\nwhiteningMatrix:'
        print whiteningMatrix


    # ======================================
    # Generate w_init

    if w_init == None:

        random_state = numpy.random.mtrand._rand

        w_init = numpy.asarray(
            random_state.normal(
                size=(n_components, n_components)
                ),
            dtype = X_whitened.dtype
            )

    if verbose:
        print '\nw_init:'
        print w_init


    # ======================================
    # ICA

    if verbose: chap( 'Running FastICA_ica_par' )

    W, nIterations = fast_ica_ica_par( X_whitened, w_init )

    if verbose:
        print '\nFound W in {0} iterations:'.format( nIterations )
        print W


    # ======================================
    # Return

    ret = Container()

    ret.W                = W
    ret.mixingMatrix     = W

    ret.K                = whiteningMatrix
    ret.whiteningMatrix  = whiteningMatrix

    ret.unmixingMatrix   = numpy.dot( W, whiteningMatrix )

    ret.X_whitened       = X_whitened
    ret.data             = data
    ret.X                = data
    ret.n_components     = n_components

    ret.nIterations      = nIterations

    ret.rotatedData      = numpy.dot( numpy.dot( W, whiteningMatrix), data ).T
    ret.rotationMatrix   = linalg.pinv( numpy.dot( W, whiteningMatrix ) )

    ret.xPoints        = list( ret.data[0,:] )
    ret.yPoints        = list( ret.data[1,:] )
    ret.xPointsRotated = list( ret.rotatedData[:,0] )
    ret.yPointsRotated = list( ret.rotatedData[:,1] )

    ret.dataRotatedBack = numpy.dot( ret.rotatedData, ret.rotationMatrix.T ) + numpy.array([ mean(l) for l in [ ret.xPoints, ret.yPoints ] ])
    ret.xPointsRotatedBack = list( ret.dataRotatedBack[:,0] )
    ret.yPointsRotatedBack = list( ret.dataRotatedBack[:,1] )

    return ret


def fast_ica_whiten( X, n_components ):

    # Centering the columns (ie the variables)
    X_mean = X.mean(axis=-1)
    X -= X_mean[:, numpy.newaxis]

    # Whitening and preprocessing by PCA
    u, d, _ = linalg.svd(X, full_matrices=False)
    del _

    K = (u / d).T[:n_components]  # see (6.33) p.140
    del u, d

    X1 = numpy.dot(K, X)

    # see (13.6) p.267 Here X1 is white and data
    # in X has been projected onto a subspace by PCA

    X1 *= numpy.sqrt( n_components ) # T: Originally p from https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/decomposition/fastica_.py#L291

    return K, X1


def fast_ica_ica_par(
        X,
        w_init,
        tol      = 1e-4,
        max_iter = 200,
        ):

    """Parallel FastICA.
    Used internally by FastICA --main loop
    """

    W = fast_ica_sym_decorrelation(w_init)
    del w_init
    p_ = float(X.shape[1])
    for ii in xrange(max_iter):

        gwtx, g_wtx = fast_ica_logcosh( numpy.dot(W, X) )

        W1 = fast_ica_sym_decorrelation(
            numpy.dot(gwtx, X.T) / p_
            - g_wtx[:, numpy.newaxis] * W
            )

        del gwtx, g_wtx

        # builtin max, abs are faster than numpy counter parts.
        lim = max(abs(abs(numpy.diag(numpy.dot(W1, W.T))) - 1))

        W = W1

        if lim < tol:
            break

    else:
        print 'FastICA did not converge. Consider increasing tolerance or the maximum number of iterations.'

    return W, ii + 1


def fast_ica_sym_decorrelation(W):
    """ Symmetric decorrelation
    i.e. W <- (W * W.T) ^{-1/2} * W
    """
    s, u = linalg.eigh(numpy.dot(W, W.T))
    # u (resp. s) contains the eigenvectors (resp. square roots of
    # the eigenvalues) of W * W.T
    return numpy.dot(numpy.dot(u * (1. / numpy.sqrt(s)), u.T), W)


def fast_ica_logcosh(x, fun_args={}):
    alpha = fun_args.get('alpha', 1.0)  # comment it out?

    x *= alpha
    gx = numpy.tanh(x, x)  # apply the tanh inplace
    g_x = numpy.empty(x.shape[0])
    # print 'gx.shape = ', gx.shape
    # print 'g_x.shape = ', g_x.shape
    # XXX compute in chunks to avoid extra allocation
    for i, gx_i in enumerate(gx):  # please don't vectorize.
        g_x[i] = (alpha * (1 - gx_i ** 2)).mean()
    return gx, g_x



########################################
# General functions
########################################

def draw_data( res, verbose=False ):

    xPoints = res.xPoints
    yPoints = res.yPoints

    # ======================================
    # Calculate covariance

    covMatMeasured = numpy.cov( xPoints, yPoints )

    xStdMeasured = sqrt(covMatMeasured[0][0])
    yStdMeasured = sqrt(covMatMeasured[1][1])

    xMeanMeasured = mean( xPoints )
    yMeanMeasured = mean( yPoints )

    corr = covMatMeasured[0][1] / ( xStdMeasured * yStdMeasured )

    if verbose:

        print '\nxMeanMeasured:'
        print xMeanMeasured
        print '\nxStdMeasured:'
        print xStdMeasured

        print '\nyMeanMeasured:'
        print yMeanMeasured
        print '\nyStdMeasured:'
        print yStdMeasured

        print '\nMeasured covariance matrix:'
        print covMatMeasured

        print '\nMeasured correlation matrix:'
        print corr

        print ''


    # ======================================
    # Make basic plot of data

    c.Clear()
    set_cmargins()

    xMin = xMeanMeasured - 3.*xStdMeasured
    xMax = xMeanMeasured + 3.*xStdMeasured
    yMin = yMeanMeasured - 3.*xStdMeasured
    yMax = yMeanMeasured + 3.*xStdMeasured
    base = get_plot_base(
        xMin = xMin, xMax = xMax, yMin = yMin, yMax = yMax,
        xTitle = 'x', yTitle = 'y',
        )
    base.Draw('P')

    Tg = ROOT.TGraph(
        len(xPoints),
        array( 'd', xPoints ),
        array( 'd', yPoints ),
        )
    Tg.SetMarkerStyle(8)
    Tg.Draw('PSAME')


    # Plot correlation line
    corrLineFn = lambda x: yMeanMeasured + corr*yStdMeasured/xStdMeasured * ( x - xMeanMeasured )
    corrLine = ROOT.TLine(
        xMin, corrLineFn(xMin),
        xMax, corrLineFn(xMax),
        )
    ROOT.SetOwnership( corrLine, False )
    corrLine.SetLineWidth(2)
    corrLine.SetLineColor(2)
    corrLine.Draw()

    xMeanLine = ROOT.TLine( xMeanMeasured, yMin, xMeanMeasured, yMax )
    ROOT.SetOwnership( xMeanLine, False )
    xMeanLine.Draw()
    yMeanLine = ROOT.TLine( xMin, yMeanMeasured, xMax, yMeanMeasured )
    ROOT.SetOwnership( yMeanLine, False )
    yMeanLine.Draw()

    rhoText = ROOT.TLatex()
    rhoText.SetNDC()
    rhoText.SetTextAlign(33)
    rhoText.DrawLatex( 1-CRightMargin-0.01, 1-CTopMargin-0.01, '#rho = {0:.4f}'.format( corr ) )

    save_c( '{0}_rawdata'.format(res.name) )


    # ======================================
    # Draw resulted rotated back

    if hasattr( res, 'xPointsRotatedBack' ):

        TgRotatedBack = ROOT.TGraph(
            len(res.xPointsRotatedBack),
            array( 'f', res.xPointsRotatedBack ),
            array( 'f', res.yPointsRotatedBack )
            )

        # if res.name == 'onlyWhitened':
        #     x = ROOT.Double(0.)
        #     y = ROOT.Double(0.)
        #     for iPoint in xrange( len(res.xPointsRotatedBack) ):
        #         TgRotatedBack.GetPoint( iPoint, x, y )
        #         print 'Point {0}: x = {1}, y = {2}'.format( iPoint, x, y )

        ROOT.SetOwnership( TgRotatedBack, False )
        TgRotatedBack.SetMarkerStyle(8)
        TgRotatedBack.SetMarkerSize(0.3)
        TgRotatedBack.SetMarkerColor(2)
        TgRotatedBack.Draw('PSAME')

        save_c( '{0}_rawdata_rotatedback'.format( res.name ) )


    # ======================================
    # Draw rotated result

    if hasattr( res, 'xPointsRotated' ):

        c.Clear()
        set_cmargins()

        xPointsRotated = res.xPointsRotated
        yPointsRotated = res.yPointsRotated

        covMatRotated = numpy.cov( xPointsRotated, yPointsRotated )
        xStdRotated   = sqrt(covMatRotated[0][0])
        yStdRotated   = sqrt(covMatRotated[1][1])
        xMeanRotated  = mean( xPointsRotated )
        yMeanRotated  = mean( yPointsRotated )
        corrRotated   = covMatRotated[0][1] / ( xStdRotated * yStdRotated )

        if verbose:

            print '\n----- Rotated system:'

            print '\ncovMatRotated:'
            print covMatRotated

            print '\nxMeanRotated:'
            print xMeanRotated
            print '\nxStdRotated:'
            print xStdRotated

            print '\nyMeanRotated:'
            print yMeanRotated
            print '\nyStdRotated:'
            print yStdRotated

            print '\ncorrRotated:'
            print corrRotated


        xMinRotated = xMeanRotated - 3.*xStdRotated
        xMaxRotated = xMeanRotated + 3.*xStdRotated
        yMinRotated = yMeanRotated - 3.*xStdRotated
        yMaxRotated = yMeanRotated + 3.*xStdRotated
        baseRotated = get_plot_base(
            xMin = xMinRotated, xMax = xMaxRotated, yMin = yMinRotated, yMax = yMaxRotated,
            xTitle = 'xRotated', yTitle = 'yRotated',
            )
        baseRotated.Draw('P')

        TgRotated = ROOT.TGraph(
            len(xPointsRotated),
            array( 'd', xPointsRotated ),
            array( 'd', yPointsRotated ),
            )
        TgRotated.SetMarkerStyle(8)
        TgRotated.SetMarkerSize(0.8)
        TgRotated.SetMarkerColor(8)

        TgRotated.Draw('PSAME')

        save_c( '{0}_rotateddata'.format(res.name) )


    if hasattr( res, 'rotationMatrix' ):
        print '\nFound rotation matrix:'
        print res.rotationMatrix




def mean( l ):
    return sum(l)/len(l)

def normalize( l ):
    norm = sqrt(sum([ val**2 for val in l ]))
    l = [ val/norm for val in l ]
    return l

def numpy_column( l ):
    return numpy.array(l)[:,numpy.newaxis]

def numpy_row( l ):
    return numpy.array(l)[:,numpy.newaxis]    


def get_winit( n_components, dtype ):
    random_state = numpy.random.mtrand._rand

    w_init = numpy.asarray(
        random_state.normal(
            size=(n_components, n_components)
            ),
        dtype = dtype
        )

    return w_init


def generate_data(
        xMean = 50.,
        xStd  = 10.,
        yMean = 100.,
        yStd  = 5.,
        nGeneratedPoints = 1000.,
        ):

    corr = 0.9
    covMat = [
        [ xStd**2,         corr*xStd*yStd ],
        [ corr*xStd*yStd,  yStd**2        ]
        ]

    numpy.random.seed(1001)
    xPoints, yPoints = numpy.random.multivariate_normal(
        [ xMean, yMean ],
        covMat,
        1000
        ).T

    xPoints = list(xPoints)
    yPoints = list(yPoints)

    return xPoints, yPoints


def chap( text ):
    print '\n' + '='*70
    print text
    print ''


class Container:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)




########################################
# End of Main
########################################
if __name__ == "__main__":
    main()