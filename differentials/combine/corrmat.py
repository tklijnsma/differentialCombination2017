import os, shutil
import logging
from os.path import *
import glob, re, copy
from collections import namedtuple

import differentials
import differentials.core as core

import combine_utils as utils
import combine


class CombineCorrelationMatrix(object):

    def __init__(self, postfit=None):
        super(CombineCorrelationMatrix, self).__init__()
        self.postfit = postfit

    



# #____________________________________________________________________
# def ComputeCorrMatrix(
#         ws,
#         redoPostfit = True,
#         postfitFilename = None,
#         onBatch = True,
#         asimov = False,
#         ):

#     wsTag = basename(ws).replace('/','').replace('.root','')
#     directory = 'corrMat_{0}_{1}'.format( datestr, wsTag )
#     corrmatFilename = join( directory, 'higgsCombine_CORRMAT_{0}.MultiDimFit.mH125.root'.format( wsTag ) )

#     if postfitFilename is None:
#         postfitFilename = join( directory, 'higgsCombine_POSTFIT_{0}.MultiDimFit.mH125.root'.format( wsTag ) )
#     else:
#         if not isfile(postfitFilename):
#             ThrowError( 'Given ws {0} does not exist'.format(postfitFilename) )

#     if redoPostfit:
#         directory = AppendNumberToDirNameUntilItDoesNotExistAnymore(directory)
#         # First regular best fit
#         BasicBestfit(
#             ws,
#             onBatch = onBatch,
#             directory = directory,
#             extraOptions = [
#                 '-m 125',
#                 '--floatOtherPOIs=1',
#                 # '--computeCovarianceMatrix=1',
#                 '--saveWorkspace',
#                 '-n _POSTFIT_{0}'.format( wsTag ),
#                 ( '-t -1' if asimov else '' )
#                 ]
#             )

#     if IsTestMode():
#         pdfIndicesToFreeze = [ 'some', 'pdfs' ]
#     else:
#         pdfIndicesToFreeze = ListOfPDFIndicesToFreeze( postfitFilename, verbose=False )

#     BasicBestfit(
#         postfitFilename,
#         onBatch = onBatch,
#         directory = directory,
#         extraOptions = [
#             '-m 125',
#             '--floatOtherPOIs=1',
#             '--algo none',
#             '--snapshotName MultiDimFit',
#             '--saveWorkspace',
#             '--computeCovarianceMatrix=1',
#             '--freezeNuisances {0}'.format( ','.join(pdfIndicesToFreeze) ),
#             '-n _CORRMAT_{0}'.format( wsTag ),
#             ( '-t -1' if asimov else '' )
#             ]
#         )
