#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, tempfile, shutil, re, subprocess, sys, traceback
from os.path import *
from copy import deepcopy
from glob import glob

from time import strftime
datestr = strftime( '%b%d' )

from Container import Container

import ROOT


import Commands

import PlotCommands

import TheoryCommands
from TheoryCommands import c
from TheoryCommands import SetCMargins
from TheoryCommands import SetPlotDir
from TheoryCommands import SaveAsRoot
from TheoryCommands import SaveC
from TheoryCommands import GetUniqueRootName
from TheoryCommands import GetPlotBase



########################################
# Main
########################################


# class TwoStageCombineScan(CombineScan):
#     pass


class CombineScan(Container):
    """Class that parses a CombineTool.py command for a scan"""

    # ======================================
    # Physics options

    doFastscan = False
    asimov     = False

    # ----------------------

    # # LUMISTUDY = True
    # LUMISTUDY = False

    # # ----------------------

    # UNCORRELATED_THEORY_UNCERTAINTIES = False

    # NO_THEORY_UNCERTAINTIES           = False

    # # ----------------------

    # PROFILE_TOTAL_XS                  = False

    # ONEDIM_SCAN_TOTALXS               = False

    # FIT_ONLY_NORMALIZATION            = False

    # # ----------------------

    # INCLUDE_BR_COUPLING_DEPENDENCY    = False

    # FIX_KAPPAV                        = False

    # MAX_KAPPAV_ONE                    = False

    # # ----------------------

    # DO_BR_UNCERTAINTIES               = False

    # FIT_RATIO_OF_BRS                  = False

    # ONEDIM_SCAN_RATIO_OF_BRS          = False

    # # ----------------------

    # DO_ONLY_ONE_KAPPA                 = False


    # ======================================
    # Combine options

    tags                              = []
    subDirectory                      = None

    SQUARE_DIST_POI_STEP              = True
    FLOAT_OTHER_POIS                  = True
    ALGO                              = 'grid'
    SAVE_WORKSPACE                    = False

    # Very basic, unlikely to ever change
    METHOD                            = 'MultiDimFit'
    DEFAULT_MASS                      = 125.
    SAVE_NLL                          = True
    SAVE_INACTIVE_POI                 = True

    default_minimizer_settings        = True
    cminDefaultMinimizerType          = 'Minuit2'
    cminDefaultMinimizerAlgo          = 'migrad'

    # Default settings for jobs
    fromPostfit                       = True
    onBatch                           = True
    queue                             = 'all.q'
    jobPriority                       = 0
    nPointsPerJob                     = None

    # These should be filled before being to submit a command
    name                        = None
    datacard                    = None
    POIs                        = []
    PhysicsModelParameters      = []
    PhysicsModelParameterRanges = []
    floatNuisances              = []
    freezeNuisances             = []


    # ======================================
    # Variables to help in the class

    deltaNLLCutOff              = 25.
    postfitRootFileBasename     = ''
    taskNamePrefix              = ''
    extraOptions = []

    # Main flag that tells the class to first perform a fastscan, on which the profiled scan will be based
    APPLY_FASTSCAN_FILTER       = False

    # Setter and getter methods for the postfitWS attribute
    _postfitWS = None
    def get_postfitWS( self ):
        return self._postfitWS
    def set_postfitWS( self, postfitWS ):
        print '\nSetting postfitWS to \'{0}\''.format(postfitWS)
        self._postfitWS = postfitWS
    postfitWS = property( get_postfitWS, set_postfitWS )

    # Setter and getter methods for the datacard attribute
    _datacard = None
    def get_datacard( self ):
        return self._datacard
    def set_datacard( self, datacard ):
        print '\nSetting datacard to \'{0}\''.format(datacard)
        self._datacard = abspath(datacard)
    datacard = property( get_datacard, set_datacard )


    # Setter and getter methods for the datacard attribute
    _fastscanRootFile = None
    def get_fastscanRootFile( self ):
        return self._fastscanRootFile
    def set_fastscanRootFile( self, fastscanRootFile ):
        print '\nSetting fastscanRootFile to \'{0}\''.format(fastscanRootFile)
        self._fastscanRootFile = fastscanRootFile
    fastscanRootFile = property( get_fastscanRootFile, set_fastscanRootFile )


    #____________________________________________________________________
    def __init__( self, container=None ):
        super( CombineScan, self ).__init__()

        # ======================================
        # Overwrite with attributes from the input container

        if not container is None:
            for attribute in container.ListAttributes( onlyVariables=True ):
                if attribute.startswith('get_') or attribute.startswith('set_'):
                    continue
                elif callable(getattr( container, attribute )):
                    continue
                elif attribute in [ 'datacard', 'fastscanRootFile', 'postfitWS' ]:
                    continue
                # print 'Copying attribute \'{0}\''.format(attribute)
                setattr( self, attribute, getattr( container, attribute ) )


    #____________________________________________________________________
    def MakeSubdirectory(self):
        if isdir( self.subDirectory ):
            self.subDirectory = Commands.AppendNumberToDirNameUntilItDoesNotExistAnymore(self.subDirectory)
        if not Commands.IsTestMode():
            os.makedirs(self.subDirectory)


    #____________________________________________________________________
    def Run( self ):

        # Make the subDirectory now, avoid race condition while waiting for the postfit
        if self.asimov:
            self.tags.append('asimov')
        if len(self.tags) > 0:
            self.subDirectory += '_' + '_'.join(self.tags)
        self.subDirectory = Commands.AppendNumberToDirNameUntilItDoesNotExistAnymore(self.subDirectory)
        self.MakeSubdirectory()
        
        if self.APPLY_FASTSCAN_FILTER:

            # ======================================
            # Make sure there is a self.postfitWS

            if self.postfitWS is None:
                self.Chapter( 'Postfit not given - Creating postfit' )

                # Create new instance from self
                # print 'Creating copied instance...'
                postfitScan = CombineScan(self)

                # Overwrite some settings
                postfitScan.subDirectory = 'postfitWSs_{0}'.format(datestr)
                postfitScan.onBatch = False

                self.postfitWS = postfitScan.CreatePostfit()


            # ======================================
            # Make sure there is a fastscan and determine accepted points

            if not self.asimov:
                # Now remove any overwriting of PhysicsModelParameters
                Commands.Warning( 'Deleting any previously PhysicsModelParameters (They will be read from the postfit)' )
                self.PhysicsModelParameters = []

            if self.fastscanRootFile is None:
                self.Chapter( 'No fastscanRootFile was given - Creating fastscanRootFile' )
                self.fastscanRootFile = self.DetermineRelevantPointsFromFastScan()

            if not Commands.IsTestMode():
                print '\nMaking basic plot of fast scan result'
                self.PlotFastScan(self.fastscanRootFile)

            print '\nDetermining accepted points from fastscanRootFile:', self.fastscanRootFile
            acceptedPoints = self.GetListOfAcceptedPoints( self.fastscanRootFile )


            # ======================================
            # Submit the scan based on the accepted points

            self.Chapter( 'Submitting profiled scan based on accepted points' )
            self.SubmitScan( acceptedPoints )


        else:
            self.SubmitScan()


    #____________________________________________________________________
    def PlotFastScan( self, fastscanRootFile ):

        print '\nAttempting to make a quick plot of {0}'.format(fastscanRootFile)

        container = TheoryCommands.GetTH2FromListOfRootFiles(
            [ fastscanRootFile ],
            self.POIs[0],
            self.POIs[1],
            verbose = False,
            defaultHValue = 999.,
            )
        container.name = basename(fastscanRootFile).replace('.root','')

        if 'ct' in self.POIs and 'cg' in self.POIs:
            x_ranges = [ -8.5, 8.5 ]
            y_ranges = [ -0.65, 0.65 ]
        else:
            x_ranges = [ -8.5, 8.5 ]
            y_ranges = [ -0.65, 0.65 ]

        PlotCommands.PlotSingle2DHistogram(
            container,
            container.xBinBoundaries[0], container.xBinBoundaries[-1],
            container.yBinBoundaries[0], container.yBinBoundaries[-1],
            self.POIs[0], self.POIs[1],
            'fastscan',
            zMax = 100.,
            palette = 'twocolor',
            getCustomContour = self.deltaNLLCutOff
            )


    #____________________________________________________________________
    def Chapter( self, txt ):
        print '\n' + '='*70
        print txt


    #____________________________________________________________________
    def CreatePostfit( self ):

        finalPostfitWS = abspath(join( self.subDirectory, 'POSTFIT_' + ( 'ASIMOV_' if self.asimov else '' ) + basename(self.datacard).replace('/','') ))

        cmd = self.ParseBestfitCommand()

        with RunInDirectory( subDirectory = self.subDirectory ):
            Commands.executeCommand( cmd )
            Commands.movefile( self.postfitRootFileBasename, finalPostfitWS )

        # Return path to the output root file
        return finalPostfitWS


    #____________________________________________________________________
    def DetermineRelevantPointsFromFastScan( self ):

        try:
            _doFastscan = self.doFastscan
            self.doFastscan = True

            _onBatch = self.onBatch
            self.onBatch = False

            _subDirectory = self.subDirectory
            self.subDirectory = join( self.subDirectory, 'postfit_and_fastscan' )

            srcPostfitWS = self.postfitWS
            self.postfitWS = abspath(join( self.subDirectory, basename(self.postfitWS) ))

            with RunInDirectory( subDirectory = self.subDirectory ):
                Commands.copyfile( srcPostfitWS, self.postfitWS )
                fastscanRootFile, cmd = self.ParseScanCommand()
                Commands.executeCommand( cmd )

        finally:
            self.doFastscan   = _doFastscan
            self.onBatch      = _onBatch
            self.subDirectory = _subDirectory


        return fastscanRootFile

    #____________________________________________________________________
    def SubmitScan( self, acceptedPoints=None ):
        self.PrepareForCommandCompilation( mode = 'scan' )

        try:
            _extraOptions = deepcopy(self.extraOptions)

            with RunInDirectory( subDirectory = self.subDirectory ):
                if acceptedPoints is None:
                    self.extraOptions.append( '--split-points {0}'.format(self.nPointsPerJob) )
                    cmd = self.ParseScanCommand()
                    Commands.executeCommand( cmd )
                else:
                    for iChunk, chunk in enumerate(chunks( [ container.iPoint for container in acceptedPoints ], self.nPointsPerJob )):
                        print '\nJob', iChunk
                        self.extraOptions = _extraOptions + [ '--doPoints ' + ','.join([ str(i) for i in chunk ]) ]
                        self.taskNamePrefix = '_' + str(iChunk)
                        cmd = self.ParseScanCommand()
                        Commands.executeCommand( cmd )

        finally:
            self.extraOptions = _extraOptions

    #____________________________________________________________________
    def PrepareForCommandCompilation( self, mode ):

        if mode == 'postfit':
            if self.datacard is None:
                Commands.ThrowError( 'Set the \'datacard\' attribute to /path/to/datacard.root' )
            elif not isfile(self.datacard):
                Commands.ThrowError( '{0} does not exist'.format(self.datacard) )
            self.datacard = abspath(self.datacard)

            if self.asimov and len(self.PhysicsModelParameters) == 0:
                Commands.ThrowError( 'PhysicsModelParameters HAS to be set when running on asimov, otherwise behavior is unspecfied!!' )

            if self.name is None:
                self.name = basename(self.datacard).replace('/','').replace('.root','')
                Commands.Warning( 'No name given; Will fill \'{0}\''.format(self.name) )

        elif mode == 'scan':

            if self.fromPostfit:
                dc = self.postfitWS
            else:
                dc = self.datacard

            if dc is None:
                Commands.ThrowError( 'Set the \'postfitWS\' attribute to /path/to/postfitWS.root' )
            elif Commands.IsTestMode():
                pass
            elif not isfile(dc):
                Commands.ThrowError( '{0} does not exist'.format(dc) )
            else:
                dc = abspath(dc)

            if self.name is None:
                self.name = basename(dc).replace('/','').replace('.root','')
                Commands.Warning( 'No name given; Will fill \'{0}\''.format(self.name) )


            if self.nPoints is None:
                Commands.ThrowError( 'Mode \'scan\' needs attribute \'nPoints\'' )

            if not self.doFastscan and self.nPointsPerJob is None:
                Commands.ThrowError( 'Mode \'scan\' needs attribute \'nPointsPerJob\' for profiled scans' )


        else:
            Commands.ThrowError( 'Mode \'{0}\' is not implemented'.format(mode) )


        if len( self.POIs ) == 0:
            Commands.Warning( 'No POIs are set; this will take the pre-defined POI set from the datacard' )

        if len( self.PhysicsModelParameters ) == 0:
            Commands.Warning( 'No physics model parameters were overwritten; all the default values are used' )

        if len( self.PhysicsModelParameterRanges ) == 0:
            if not mode == 'scan':
                Commands.Warning( 'No physics model parameter ranges were overwritten; all the default ranges are used' )
            else:
                Commands.ThrowError( 'No physics model parameter ranges were overwritten; not allowed for mode \'scan\'' )


    #____________________________________________________________________
    def ParseScanCommand( self ):
        self.PrepareForCommandCompilation( mode = 'scan' )
        commandType = 'FASTSCAN' if self.doFastscan else 'SCAN'
        taskName = commandType + self.taskNamePrefix + '_' + self.name

        cmd = []

        if self.onBatch:
            cmd.append( 'combineTool.py' )
        else:
            cmd.append( 'combine' )

        if self.fromPostfit:
            cmd.append(self.postfitWS)
        else:
            cmd.append(self.datacard)

        cmd.extend([
            '-n _{0}'.format( taskName ),
            '--algo=grid',
            '--points={0}'.format(self.nPoints)
            ])

        # Set all the common stuff in separate method
        self.CommonCommandSettings(cmd, taskName)

        if self.doFastscan:
            cmd.append( '--fastScan' )
            output = abspath( 'higgsCombine_' + taskName + '.MultiDimFit.mH125.root' )
        elif self.fromPostfit:
            cmd.append( '--snapshotName MultiDimFit' )
            cmd.append( '--skipInitialFit' )
            

        # ======================================
        # Settings for grid

        if self.doFastscan:
            return output, cmd
        else:
            return cmd


    #____________________________________________________________________
    def CommonCommandSettings(self, cmd, taskName):

        if self.asimov:
            cmd.append( '-t -1' )
        if len(self.POIs) > 0:
            cmd.append( '-P ' + ' -P '.join(self.POIs) )
        if len(self.PhysicsModelParameters) > 0:
            cmd.append( '--setPhysicsModelParameters '      + ','.join(self.PhysicsModelParameters) )
        if len(self.PhysicsModelParameterRanges) > 0:
            cmd.append( '--setPhysicsModelParameterRanges ' + ':'.join(self.PhysicsModelParameterRanges) )
        if len(self.floatNuisances) > 0:
            cmd.append( '--floatNuisances ' + ','.join(self.floatNuisances) )
        if len(self.freezeNuisances) > 0:
            cmd.append( '--freezeNuisances ' + ','.join(self.freezeNuisances) )

        cmd.extend([
            '-M {0}'.format( self.METHOD ),
            '-m {0}'.format( self.DEFAULT_MASS ),
            ])

        if self.default_minimizer_settings:
            cmd.extend([
                '--cminDefaultMinimizerType {0}'.format( self.cminDefaultMinimizerType ),
                '--cminDefaultMinimizerAlgo {0}'.format( self.cminDefaultMinimizerAlgo ),            
                ])

        if self.SAVE_NLL:
            cmd.append( '--saveNLL' )

        if self.SAVE_INACTIVE_POI:
            cmd.append( '--saveInactivePOI 1' )
        else:
            cmd.append( '--saveInactivePOI 0' )

        if self.FLOAT_OTHER_POIS:
            cmd.append( '--floatOtherPOIs=1' )
        else:
            cmd.append( '--floatOtherPOIs=0' )

        # For scans: Hard-coded *not* saving the workspace
        # cmd.append( '--saveWorkspace' )

        if len(self.extraOptions) > 0:
            cmd.extend(self.extraOptions)

        if self.onBatch:
            if 't3' in os.environ['HOSTNAME']:

                if not self.queue in [ 'all.q', 'long.q', 'short.q' ]:
                    Commands.ThrowError( 'Queue \'{0}\' is not available on PSI'.format(self.queue) )
                if self.jobPriority != 0:
                    cmd.append(
                        '--job-mode psi --task-name {0} --sub-opts=\'-q {1} -p {2}\' '.format( taskName, self.queue, self.jobPriority ),
                        )
                else:
                    cmd.append(
                        '--job-mode psi --task-name {0} --sub-opts=\'-q {1}\' '.format( taskName, self.queue ),
                        )
            else:
                Commands.ThrowError( 'Only jobs submitted from T3 are implemented now' )


    #____________________________________________________________________
    def GetListOfAcceptedPoints( self, fastscanRootFile, deltaNLLCutOff = None, verbose=False ):

        if deltaNLLCutOff is None:
            deltaNLLCutOff = self.deltaNLLCutOff

        if not isfile(fastscanRootFile):
            if Commands.IsTestMode():
                print '\n[TESTMODE] No file \'{0}\'; Returning some bogus accepted points'.format(fastscanRootFile)
                return [ Container(iPoint=i) for i in xrange(5) ]
            else:
                Commands.ThrowError( 'File \'{0}\' does not exist'.format(fastscanRootFile) )

        fastscanRootFp = ROOT.TFile.Open( fastscanRootFile )

        if not fastscanRootFp.GetListOfKeys().Contains('limit'):
            fastscanRootFp.Close()
            Commands.ThrowError( 'There is no tree \'limit\' in', fastscanRootFile )

        acceptedPoints = []
        rejectedPoints = []
        tree = fastscanRootFp.Get('limit')
        for iEvent, event in enumerate(tree):

            container = Container()
            container.iPoint   = iEvent
            container.deltaNLL = event.deltaNLL
            container.POIvals  = [ getattr( event, POI ) for POI in self.POIs ]
            
            if iEvent == 0:
                self.bestfit = container
                continue

            if container.deltaNLL <= deltaNLLCutOff:
                acceptedPoints.append( container )
            else:
                rejectedPoints.append( container )

        fastscanRootFp.Close()

        if verbose:
            print '\nRejected points:'
            for container in rejectedPoints:
                line = [
                    '{0:7}'.format(container.iPoint),
                    'deltaNLL = {0:+10.2f}'.format( container.deltaNLL )
                    ]
                for POI, POIval in zip( self.POIs, container.POIvals ):
                    line.append( '{0:10} = {1:+7.2f}'.format( POI, POIval ) )
                print ' | '.join(line)

            print '\nAccepted points:'
            for container in acceptedPoints:
                line = [
                    '{0:7}'.format(container.iPoint),
                    'deltaNLL = {0:+10.2f}'.format( container.deltaNLL )
                    ]
                for POI, POIval in zip( self.POIs, container.POIvals ):
                    line.append( '{0:10} = {1:+7.2f}'.format( POI, POIval ) )
                print ' | '.join(line)

        return acceptedPoints


    #____________________________________________________________________
    def ParseBestfitCommand( self ):

        self.PrepareForCommandCompilation( mode = 'postfit' )

        cmd = []

        if self.onBatch:
            cmd.append( 'combineTool.py' )
        else:
            cmd.append( 'combine' )

        taskName = 'POSTFIT_' + ( 'ASIMOV_' if self.asimov else '' ) + self.name

        cmd.extend([
            self.datacard,
            '-n _{0}'.format(taskName),
            # '--algo=grid',
            # '--points={0}'.format(self.nPoints)
            ])

        self.CommonCommandSettings(cmd, taskName)

        # For the postfit: Hard-coded saving the workspace
        cmd.append( '--saveWorkspace' )

        self.postfitRootFileBasename = 'higgsCombine_{0}.MultiDimFit.mH125.root'.format( taskName )

        return cmd


    #____________________________________________________________________
    def ParseBestfitCommand_old( self ):
        self.PrepareForCommandCompilation( mode = 'postfit' )

        cmd = []

        if self.onBatch:
            cmd.append( 'combineTool.py' )
        else:
            cmd.append( 'combine' )

        
        scanName = '_POSTFIT_' + ( 'ASIMOV_' if self.asimov else '' ) + self.name
        cmd.extend([
            # 'combine',
            self.datacard,
            '-n {0}'.format( scanName ),
            ])
        self.postfitRootFileBasename = 'higgsCombine{0}.MultiDimFit.mH125.root'.format( scanName )

        if self.doFastscan:
            cmd.append( '--fastScan' )
        if self.asimov:
            cmd.append( '-t -1' )


        if len(self.POIs) > 0:
            cmd.append( '-P ' + ' -P '.join(self.POIs) )

        if len(self.PhysicsModelParameters) > 0:
            cmd.append( '--setPhysicsModelParameters '      + ','.join(self.PhysicsModelParameters) )

        if len(self.PhysicsModelParameterRanges) > 0:
            cmd.append( '--setPhysicsModelParameterRanges ' + ':'.join(self.PhysicsModelParameterRanges) )


        # ======================================
        # Very basic combine options

        cmd.extend([
            '-M {0}'.format( self.METHOD ),
            '--cminDefaultMinimizerType {0}'.format( self.cminDefaultMinimizerType ),
            '--cminDefaultMinimizerAlgo {0}'.format( self.cminDefaultMinimizerAlgo ),
            '-m {0}'.format( self.DEFAULT_MASS ),
            ])

        if self.SAVE_NLL:
            cmd.append( '--saveNLL' )

        if self.SAVE_INACTIVE_POI:
            cmd.append( '--saveInactivePOI 1' )
        else:
            cmd.append( '--saveInactivePOI 0' )

        if self.FLOAT_OTHER_POIS:
            cmd.append( '--floatOtherPOIs=1' )
        else:
            cmd.append( '--floatOtherPOIs=0' )

        # For the postfit: Hard-coded saving the workspace
        cmd.append( '--saveWorkspace' )


        # ======================================
        # Settings for grid

        if self.onBatch:
            if 't3' in os.environ['HOSTNAME']:
                taskName = 'POSTFIT_' + self.name

                # Hardcode short queue for postfit jobs
                queue = 'short.q'

                if not queue in [ 'all.q', 'long.q', 'short.q' ]:
                    Commands.ThrowError( 'Queue \'{0}\' is not available on PSI'.format(queue) )
                if self.jobPriority != 0:
                    cmd.append(
                        '--job-mode psi --task-name {0} --sub-opts=\'-q {1} -p {2}\' '.format( taskName, queue, self.jobPriority ),
                        )
                else:
                    cmd.append(
                        '--job-mode psi --task-name {0} --sub-opts=\'-q {1}\' '.format( taskName, queue ),
                        )
            else:
                Commands.ThrowError( 'Only jobs submitted from T3 are implemented now' )

        return cmd




#____________________________________________________________________
class RunInDirectory():
    """Context manager to create and go into and out of a directory"""

    def __init__(self, subDirectory=None ):
        self._active = False
        if not subDirectory is None:
            self._active = True
            self.backDir = os.getcwd()
            self.subDirectory = subDirectory

    def __enter__(self):
        if self._active:
            if Commands.IsTestMode():
                print '\n[TESTMODE] Would now create/go into \'{0}\''.format(self.subDirectory)
            else:
                print ''
                if not isdir( self.subDirectory ):
                    print 'Creating \'{0}\''.format(self.subDirectory)
                    os.makedirs( self.subDirectory )
                print 'Entering \'{0}\''.format(self.subDirectory)
                os.chdir( self.subDirectory )
        return self

    def __exit__(self, *args):
        if self._active:
            os.chdir( self.backDir )


#____________________________________________________________________
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]
