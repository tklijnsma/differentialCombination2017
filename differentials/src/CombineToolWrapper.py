#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, tempfile, shutil, re, subprocess, sys, traceback
from os.path import *
from copy import deepcopy, copy
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

#____________________________________________________________________
class CombineConfig(Container):

    SQUARE_DIST_POI_STEP              = True
    FLOAT_OTHER_POIS                  = True
    ALGO                              = 'grid'
    SAVE_WORKSPACE                    = False

    # Very basic, unlikely to ever change
    METHOD                            = 'MultiDimFit'
    DEFAULT_MASS                      = 125.
    SAVE_NLL                          = True
    SAVE_INACTIVE_POI                 = True

    # Default settings for jobs
    fromPostfit                       = False
    onBatch                           = True
    queue                             = 'all.q'
    jobPriority                       = 0

    def __init__(self, args, *otherargs, **kwargs):
        super(CombineConfig, self).__init__(*otherargs, **kwargs)
        self.datacard                    = 'somedatacard.root'
        self.subDirectory                = ''

        self.POIs                        = []
        self.PhysicsModelParameters      = []
        self.hardPhysicsModelParameters = []
        self.PhysicsModelParameterRanges = []
        self.floatNuisances              = []
        self.freezeNuisances             = []

        self.extraOptions                = []

        self.nPointsPerJob               = 8

        self.minimizer_settings = [
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            ]

        # Save the command line args in the config
        self.args = args

        if args.asimov:
            self.asimov = True
        else:
            self.asimov = False


    def make_unique_directory(self):
        self.subDirectory = Commands.make_unique_directory(self.subDirectory)

    def get_name(self):
        return basename(self.datacard).replace('.root','')

    def fix_parameter_at_value(self, name, val):
        self.PhysicsModelParameters.append('{0}={1}'.format(name, val))
        self.hardPhysicsModelParameters.append('{0}={1}'.format(name, val))
        self.freezeNuisances.append(name)


#____________________________________________________________________
class BaseCombineScan(Container):

    debug = True
    verbose = True

    def __init__(self, config):
        self.input = config

        # Attributes that may need to be altered from call to call
        self.onBatch = self.input.onBatch
        self.datacard = abspath(self.input.datacard)
        self.nPoints = self.input.nPoints
        self.nPointsPerJob = self.input.nPointsPerJob
        self.subDirectory = self.input.subDirectory
        self.extraOptions = self.input.extraOptions


    def print_debug(self, txt):
        if self.debug:
            print '[CTW DEBUG]',txt

    def print_info(self, txt):
        if self.verbose:
            print '[CTW INFO]',txt


    def get_output(self):
        with Commands.EnterDirectory(self.subDirectory, verbose=False):
            output = 'higgsCombine{0}.{1}.mH{2}.root'.format(self.get_task_name(), self.input.METHOD, int(self.input.DEFAULT_MASS))
            output = abspath(output)
        return output

    def parse_command(self):
        taskName = self.get_task_name()

        cmd = []

        if self.input.onBatch:
            cmd.append( 'combineTool.py' )
        else:
            cmd.append( 'combine' )

        cmd.extend([
            self.datacard,
            '-n {0}'.format(taskName),
            '-M {0}'.format( self.input.METHOD ),
            '-m {0}'.format( self.input.DEFAULT_MASS ),
            ])

        if self.input.asimov:
            cmd.append( '-t -1' )

        cmd.extend(self.input.minimizer_settings)

        if self.input.SAVE_NLL:
            cmd.append( '--saveNLL' )

        if self.input.SAVE_INACTIVE_POI:
            cmd.append( '--saveInactivePOI 1' )
        else:
            cmd.append( '--saveInactivePOI 0' )

        if self.input.FLOAT_OTHER_POIS:
            cmd.append( '--floatOtherPOIs=1' )
        else:
            cmd.append( '--floatOtherPOIs=0' )

        if len(self.extraOptions) > 0:
            cmd.extend(self.extraOptions)

        if self.onBatch:
            if 't3' in os.environ['HOSTNAME']:
                if not self.input.queue in [ 'all.q', 'long.q', 'short.q' ]:
                    Commands.throw_error( 'Queue \'{0}\' is not available on PSI'.format(self.input.queue) )
                if self.input.jobPriority != 0:
                    cmd.append(
                        '--job-mode psi --task-name {0} --sub-opts=\'-q {1} -p {2}\' '.format( taskName, self.input.queue, self.input.jobPriority ),
                        )
                else:
                    cmd.append(
                        '--job-mode psi --task-name {0} --sub-opts=\'-q {1}\' '.format( taskName, self.input.queue ),
                        )
            else:
                Commands.throw_error( 'Only jobs submitted from T3 are implemented now' )

        return cmd


    def run(self):
        cmd = self.parse_command()
        with Commands.EnterDirectory(self.subDirectory):
            Commands.execute_command(cmd)


#____________________________________________________________________
class CombinePostfit(BaseCombineScan):

    def __init__(self, *args, **kwargs):
        super(CombinePostfit, self).__init__(*args, **kwargs)
        self.subDirectory = join(self.subDirectory, 'postfit_and_fastscan')
        self.onBatch = False

    def get_task_name(self):
        return '_POSTFIT_' + ( 'ASIMOV_' if self.input.asimov else '' ) + self.input.get_name()

    def parse_command(self):
        cmd = super(CombinePostfit, self).parse_command()

        if self.input.asimov and len(self.input.PhysicsModelParameters)==0:
            Commands.throw_error('PhysicsModelParameters *have* to be set when running asimov, or the best fit may make no sense')

        if len(self.input.POIs) > 0:
            cmd.append( '-P ' + ' -P '.join(self.input.POIs) )
        if len(self.input.PhysicsModelParameterRanges) > 0:
            cmd.append( '--setPhysicsModelParameterRanges ' + ':'.join(self.input.PhysicsModelParameterRanges) )
        if len(self.input.floatNuisances) > 0:
            cmd.append( '--floatNuisances ' + ','.join(self.input.floatNuisances) )
        if len(self.input.freezeNuisances) > 0:
            cmd.append( '--freezeNuisances ' + ','.join(self.input.freezeNuisances) )

        allPhysicsModelParameters = self.input.PhysicsModelParameters + self.input.hardPhysicsModelParameters
        if len(allPhysicsModelParameters) > 0:
            cmd.append( '--setPhysicsModelParameters ' + ','.join(allPhysicsModelParameters) )

        # For the postfit: Hard-coded saving the workspace
        cmd.append( '--saveWorkspace' )

        return cmd


#____________________________________________________________________
class CombineScan(BaseCombineScan):
    def __init__(self, *args, **kwargs):
        super(CombineScan, self).__init__(*args, **kwargs)
        self.freezeNuisances = self.input.freezeNuisances[:]

    def get_task_name(self):
        return '_SCAN_' + ( 'ASIMOV_' if self.input.asimov else '' ) + self.input.get_name()

    def parse_command(self):
        cmd = super(CombineScan, self).parse_command()
        cmd.extend([
            '--algo=grid',
            '--points={0}'.format(self.nPoints),
            '--split-points {0}'.format(self.nPointsPerJob)
            ])
        if len(self.input.POIs) > 0:
            cmd.append( '-P ' + ' -P '.join(self.input.POIs) )
        if len(self.input.PhysicsModelParameterRanges) > 0:
            cmd.append( '--setPhysicsModelParameterRanges ' + ':'.join(self.input.PhysicsModelParameterRanges) )
        if len(self.input.floatNuisances) > 0:
            cmd.append( '--floatNuisances ' + ','.join(self.input.floatNuisances) )
        if len(self.freezeNuisances) > 0:
            cmd.append( '--freezeNuisances ' + ','.join(self.freezeNuisances) )        

        allPhysicsModelParameters = self.input.PhysicsModelParameters + self.input.hardPhysicsModelParameters
        if len(allPhysicsModelParameters) > 0:
            cmd.append( '--setPhysicsModelParameters ' + ','.join(allPhysicsModelParameters) )

        return cmd


#____________________________________________________________________
class CombineScanFromPostFit(BaseCombineScan):
    def __init__(self, *args, **kwargs):
        super(CombineScanFromPostFit, self).__init__(*args, **kwargs)
        self.freezeNuisances = self.input.freezeNuisances[:]
    
    def get_task_name(self):
        return '_SCAN_' + ( 'ASIMOV_' if self.input.asimov else '' ) + self.input.get_name()

    def parse_command(self):
        cmd = super(CombineScanFromPostFit, self).parse_command()
        cmd.extend([
            '--algo=grid',
            '--points={0}'.format(self.nPoints),
            '--snapshotName MultiDimFit',
            '--skipInitialFit',
            '--split-points {0}'.format(self.nPointsPerJob)
            ])

        if len(self.input.POIs) > 0:
            cmd.append( '-P ' + ' -P '.join(self.input.POIs) )
        if len(self.input.PhysicsModelParameterRanges) > 0:
            cmd.append( '--setPhysicsModelParameterRanges ' + ':'.join(self.input.PhysicsModelParameterRanges) )
        if len(self.input.floatNuisances) > 0:
            cmd.append( '--floatNuisances ' + ','.join(self.input.floatNuisances) )
        if len(self.freezeNuisances) > 0:
            cmd.append( '--freezeNuisances ' + ','.join(self.freezeNuisances) )
        
        # Not passing all the PhysicsModelParameters; most should be read from the postfit
        # if len(self.input.PhysicsModelParameters) > 0:
        #     cmd.append( '--setPhysicsModelParameters '      + ','.join(self.input.PhysicsModelParameters) )
        # Only some specific 'hard' ones are passed
        if len(self.input.hardPhysicsModelParameters) > 0:
            cmd.append( '--setPhysicsModelParameters ' + ','.join(self.input.hardPhysicsModelParameters) )

        return cmd

    def run(self, postfitWS=None):
        if not(postfitWS is None):
            self.datacard = postfitWS
        super(CombineScanFromPostFit, self).run()


#____________________________________________________________________
class CombineFastScan(BaseCombineScan):

    def __init__(self, *args, **kwargs):
        super(CombineFastScan, self).__init__(*args, **kwargs)
        self.subDirectory = join(self.subDirectory, 'postfit_and_fastscan')
        self.onBatch = False

    def get_task_name(self):
        return '_FASTSCAN_' + ( 'ASIMOV_' if self.input.asimov else '' ) + self.input.get_name()

    def parse_command(self):
        cmd = super(CombineFastScan, self).parse_command()
        cmd.extend([
            '--algo=grid',
            '--points={0}'.format(self.nPoints),
            '--fastScan',
            '--snapshotName MultiDimFit',
            '--skipInitialFit',
            ])

        if len(self.input.POIs) > 0:
            cmd.append( '-P ' + ' -P '.join(self.input.POIs) )
        if len(self.input.PhysicsModelParameterRanges) > 0:
            cmd.append( '--setPhysicsModelParameterRanges ' + ':'.join(self.input.PhysicsModelParameterRanges) )
        if len(self.input.floatNuisances) > 0:
            cmd.append( '--floatNuisances ' + ','.join(self.input.floatNuisances) )
        if len(self.input.freezeNuisances) > 0:
            cmd.append( '--freezeNuisances ' + ','.join(self.input.freezeNuisances) )
        
        # Not passing all the PhysicsModelParameters; most should be read from the postfit
        # if len(self.input.PhysicsModelParameters) > 0:
        #     cmd.append( '--setPhysicsModelParameters '      + ','.join(self.input.PhysicsModelParameters) )
        # Only some specific 'hard' ones are passed
        if len(self.input.hardPhysicsModelParameters) > 0:
            cmd.append( '--setPhysicsModelParameters ' + ','.join(self.input.hardPhysicsModelParameters) )

        return cmd

    def run(self, postfitWS):
        self.datacard = postfitWS
        super(CombineFastScan, self).run()


#____________________________________________________________________
class CombinePointwiseScan(BaseCombineScan):
    """docstring for CombinePointwiseScan"""

    def __init__(self, *args, **kwargs):
        super(CombinePointwiseScan, self).__init__(*args, **kwargs)
        if hasattr(self.input, 'deltaNLLCutOff'):
            self.deltaNLLCutOff = self.input.deltaNLLCutOff
        else:
            self.deltaNLLCutOff = 50.
        self._outputs = []

    def get_task_name_without_number(self):
        return '_SCAN_' + ( 'ASIMOV_' if self.input.asimov else '' ) + self.input.get_name()

    def parse_command(self):
        cmd = super(CombinePointwiseScan, self).parse_command()
        cmd.extend([
            '--algo=grid',
            '--points={0}'.format(self.nPoints),
            '--snapshotName MultiDimFit',
            '--skipInitialFit',
            ])

        if len(self.input.POIs) > 0:
            cmd.append( '-P ' + ' -P '.join(self.input.POIs) )
        if len(self.input.PhysicsModelParameterRanges) > 0:
            cmd.append( '--setPhysicsModelParameterRanges ' + ':'.join(self.input.PhysicsModelParameterRanges) )
        if len(self.input.floatNuisances) > 0:
            cmd.append( '--floatNuisances ' + ','.join(self.input.floatNuisances) )
        if len(self.input.freezeNuisances) > 0:
            cmd.append( '--freezeNuisances ' + ','.join(self.input.freezeNuisances) )

        # Not passing all the PhysicsModelParameters; most should be read from the postfit
        # if len(self.input.PhysicsModelParameters) > 0:
        #     cmd.append( '--setPhysicsModelParameters '      + ','.join(self.input.PhysicsModelParameters) )
        # Only some specific 'hard' ones are passed
        if len(self.input.hardPhysicsModelParameters) > 0:
            cmd.append( '--setPhysicsModelParameters ' + ','.join(self.input.hardPhysicsModelParameters) )

        return cmd

    def run(self, postfitWS, fastscanFile):
        self.datacard = postfitWS
        self.fastscanFile = fastscanFile

        # Base extraOptions
        _extraOptions = self.extraOptions

        accepted_points = self.list_accepted_points(self.fastscanFile)

        with Commands.EnterDirectory(self.subDirectory):
            for iChunk, chunk in enumerate(chunks(accepted_points, self.nPointsPerJob)):
                print '\nJob', iChunk
                self.extraOptions = _extraOptions + [ '--doPoints ' + ','.join([ str(i) for i in chunk ]) ]
                self.get_task_name = lambda: self.get_task_name_without_number() + '_' + str(iChunk)

                cmd = self.parse_command()
                Commands.execute_command( cmd )


    def list_accepted_points(self, fastscanFile):

        if not isfile(fastscanFile):
            if Commands.is_test_mode():
                self.print_info('[TESTMODE] No file \'{0}\'; Returning some bogus accepted points'.format(fastscanFile))
                # return [ Container(iPoint=i) for i in xrange(2*self.nPointsPerJob) ]
                return range(2*self.nPointsPerJob)
            else:
                Commands.throw_error( 'File \'{0}\' does not exist'.format(fastscanFile) )


        with Commands.OpenRootFile(fastscanFile) as fastscanFp:

            if not fastscanFp.GetListOfKeys().Contains('limit'):
                Commands.throw_error( 'There is no tree \'limit\' in', fastscanFile )

            acceptedPoints = []
            rejectedPoints = []
            tree = fastscanFp.Get('limit')
            for iEvent, event in enumerate(tree):
                container = Container()
                container.iPoint   = iEvent
                container.deltaNLL = event.deltaNLL
                container.POIvals  = [ getattr( event, POI ) for POI in self.input.POIs ]
                if iEvent == 0:
                    self.bestfit = container
                    continue
                if container.deltaNLL <= self.deltaNLLCutOff:
                    acceptedPoints.append( container )
                else:
                    rejectedPoints.append( container )

        if self.verbose:
            self.print_info('Rejected points:')
            for container in rejectedPoints:
                line = [
                    '{0:7}'.format(container.iPoint),
                    'deltaNLL = {0:+10.2f}'.format( container.deltaNLL )
                    ]
                for POI, POIval in zip( self.input.POIs, container.POIvals ):
                    line.append( '{0:10} = {1:+7.2f}'.format( POI, POIval ) )
                self.print_info(' | '.join(line))

            self.print_info('Accepted points:')
            for container in acceptedPoints:
                line = [
                    '{0:7}'.format(container.iPoint),
                    'deltaNLL = {0:+10.2f}'.format( container.deltaNLL )
                    ]
                for POI, POIval in zip( self.input.POIs, container.POIvals ):
                    line.append( '{0:10} = {1:+7.2f}'.format( POI, POIval ) )
                self.print_info(' | '.join(line))

        return [ c.iPoint for c in acceptedPoints ]


#____________________________________________________________________

class CombineScan_OLD(Container):
    """Class that parses a CombineTool.py command for a scan"""

    doFastscan = False
    asimov     = False

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
    extraOptions = []


    # ======================================
    # Variables to help in the class

    deltaNLLCutOff              = 25.
    postfitRootFileBasename     = ''
    taskNamePrefix              = ''

    # Main flag that tells the class to first perform a fastscan, on which the profiled scan will be based
    APPLY_FASTSCAN_FILTER       = False

    # Setter and getter methods for the postfitWS attribute
    _postfitWS = None
    def get_postfit_ws( self ):
        return self._postfitWS
    def set_postfit_ws( self, postfitWS ):
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
    def get_fastscan_root_file( self ):
        return self._fastscanRootFile
    def set_fastscan_root_file( self, fastscanRootFile ):
        print '\nSetting fastscanRootFile to \'{0}\''.format(fastscanRootFile)
        self._fastscanRootFile = fastscanRootFile
    fastscanRootFile = property( get_fastscanRootFile, set_fastscanRootFile )


    #____________________________________________________________________
    def __init__( self, container=None ):
        super( CombineScan, self ).__init__()

        # ======================================
        # Overwrite with attributes from the input container

        if not container is None:
            for attribute in container.list_attributes( onlyVariables=True ):
                if attribute.startswith('get_') or attribute.startswith('set_'):
                    continue
                elif callable(getattr( container, attribute )):
                    continue
                elif attribute in [ 'datacard', 'fastscanRootFile', 'postfitWS' ]:
                    continue
                # print 'Copying attribute \'{0}\''.format(attribute)
                setattr( self, attribute, getattr( container, attribute ) )


    #____________________________________________________________________
    def make_subdirectory(self):
        if isdir( self.subDirectory ):
            self.subDirectory = Commands.make_unique_directory(self.subDirectory)
        if not Commands.is_test_mode():
            os.makedirs(self.subDirectory)


    #____________________________________________________________________
    def run( self ):

        # Make the subDirectory now, avoid race condition while waiting for the postfit
        if self.asimov:
            self.tags.append('asimov')
        if len(self.tags) > 0:
            self.subDirectory += '_' + '_'.join(self.tags)
        self.subDirectory = Commands.make_unique_directory(self.subDirectory)
        self.make_subdirectory()
        
        if self.APPLY_FASTSCAN_FILTER:

            # ======================================
            # Make sure there is a self.postfitWS

            if self.postfitWS is None:
                self.chapter( 'Postfit not given - Creating postfit' )

                # Create new instance from self
                # print 'Creating copied instance...'
                postfitScan = CombineScan(self)

                # Overwrite some settings
                postfitScan.subDirectory = 'postfitWSs_{0}'.format(datestr)
                postfitScan.onBatch = False

                self.postfitWS = postfitScan.create_postfit()


            # ======================================
            # Make sure there is a fastscan and determine accepted points

            if not self.asimov:
                # Now remove any overwriting of PhysicsModelParameters
                Commands.warning( 'Deleting any previously PhysicsModelParameters (They will be read from the postfit)' )
                self.PhysicsModelParameters = []

            if self.fastscanRootFile is None:
                self.chapter( 'No fastscanRootFile was given - Creating fastscanRootFile' )
                self.fastscanRootFile = self.determine_relevant_points_from_fast_scan()

            if not Commands.is_test_mode():
                print '\nMaking basic plot of fast scan result'
                self.plot_fast_scan(self.fastscanRootFile)

            print '\nDetermining accepted points from fastscanRootFile:', self.fastscanRootFile
            acceptedPoints = self.get_list_of_accepted_points( self.fastscanRootFile )


            # ======================================
            # Submit the scan based on the accepted points

            self.chapter( 'Submitting profiled scan based on accepted points' )
            self.submit_scan( acceptedPoints )


        else:
            self.submit_scan()


    #____________________________________________________________________
    def plot_fast_scan( self, fastscanRootFile ):

        print '\nAttempting to make a quick plot of {0}'.format(fastscanRootFile)

        container = TheoryCommands.get_TH2_from_list_of_root_files(
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

        PlotCommands.plot_single_TH2(
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
    def chapter( self, txt ):
        print '\n' + '='*70
        print txt


    #____________________________________________________________________
    def create_postfit( self ):

        finalPostfitWS = abspath(join( self.subDirectory, 'POSTFIT_' + ( 'ASIMOV_' if self.asimov else '' ) + basename(self.datacard).replace('/','') ))

        cmd = self.parse_bestfit_command()

        with RunInDirectory( subDirectory = self.subDirectory ):
            Commands.execute_command( cmd )
            Commands.movefile( self.postfitRootFileBasename, finalPostfitWS )

        # Return path to the output root file
        return finalPostfitWS


    #____________________________________________________________________
    def determine_relevant_points_from_fast_scan( self ):

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
                fastscanRootFile, cmd = self.parse_scan_command()
                Commands.execute_command( cmd )

        finally:
            self.doFastscan   = _doFastscan
            self.onBatch      = _onBatch
            self.subDirectory = _subDirectory


        return fastscanRootFile

    #____________________________________________________________________
    def submit_scan( self, acceptedPoints=None ):
        self.prepare_for_command_compilation( mode = 'scan' )

        try:
            _extraOptions = deepcopy(self.extraOptions)

            with RunInDirectory( subDirectory = self.subDirectory ):
                if acceptedPoints is None:
                    self.extraOptions.append( '--split-points {0}'.format(self.nPointsPerJob) )
                    cmd = self.parse_scan_command()
                    Commands.execute_command( cmd )
                else:
                    for iChunk, chunk in enumerate(chunks( [ container.iPoint for container in acceptedPoints ], self.nPointsPerJob )):
                        print '\nJob', iChunk
                        self.extraOptions = _extraOptions + [ '--doPoints ' + ','.join([ str(i) for i in chunk ]) ]
                        self.taskNamePrefix = '_' + str(iChunk)
                        cmd = self.parse_scan_command()
                        Commands.execute_command( cmd )

        finally:
            self.extraOptions = _extraOptions

    #____________________________________________________________________
    def prepare_for_command_compilation( self, mode ):

        if mode == 'postfit':
            if self.datacard is None:
                Commands.throw_error( 'Set the \'datacard\' attribute to /path/to/datacard.root' )
            elif not isfile(self.datacard):
                Commands.throw_error( '{0} does not exist'.format(self.datacard) )
            self.datacard = abspath(self.datacard)

            if self.asimov and len(self.PhysicsModelParameters) == 0:
                Commands.throw_error( 'PhysicsModelParameters HAS to be set when running on asimov, otherwise behavior is unspecfied!!' )

            if self.name is None:
                self.name = basename(self.datacard).replace('/','').replace('.root','')
                Commands.warning( 'No name given; Will fill \'{0}\''.format(self.name) )

        elif mode == 'scan':

            if self.fromPostfit:
                dc = self.postfitWS
            else:
                dc = self.datacard

            if dc is None:
                Commands.throw_error( 'Set the \'postfitWS\' attribute to /path/to/postfitWS.root' )
            elif Commands.is_test_mode():
                pass
            elif not isfile(dc):
                Commands.throw_error( '{0} does not exist'.format(dc) )
            else:
                dc = abspath(dc)

            if self.name is None:
                self.name = basename(dc).replace('/','').replace('.root','')
                Commands.warning( 'No name given; Will fill \'{0}\''.format(self.name) )


            if self.nPoints is None:
                Commands.throw_error( 'Mode \'scan\' needs attribute \'nPoints\'' )

            if not self.doFastscan and self.nPointsPerJob is None:
                Commands.throw_error( 'Mode \'scan\' needs attribute \'nPointsPerJob\' for profiled scans' )


        else:
            Commands.throw_error( 'Mode \'{0}\' is not implemented'.format(mode) )


        if len( self.POIs ) == 0:
            Commands.warning( 'No POIs are set; this will take the pre-defined POI set from the datacard' )

        if len( self.PhysicsModelParameters ) == 0:
            Commands.warning( 'No physics model parameters were overwritten; all the default values are used' )

        if len( self.PhysicsModelParameterRanges ) == 0:
            if not mode == 'scan':
                Commands.warning( 'No physics model parameter ranges were overwritten; all the default ranges are used' )
            else:
                Commands.throw_error( 'No physics model parameter ranges were overwritten; not allowed for mode \'scan\'' )


    #____________________________________________________________________
    def parse_scan_command( self ):
        self.prepare_for_command_compilation( mode = 'scan' )
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
        self.common_command_settings(cmd, taskName)

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
    def common_command_settings(self, cmd, taskName):

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
                    Commands.throw_error( 'Queue \'{0}\' is not available on PSI'.format(self.queue) )
                if self.jobPriority != 0:
                    cmd.append(
                        '--job-mode psi --task-name {0} --sub-opts=\'-q {1} -p {2}\' '.format( taskName, self.queue, self.jobPriority ),
                        )
                else:
                    cmd.append(
                        '--job-mode psi --task-name {0} --sub-opts=\'-q {1}\' '.format( taskName, self.queue ),
                        )
            else:
                Commands.throw_error( 'Only jobs submitted from T3 are implemented now' )


    #____________________________________________________________________
    def get_list_of_accepted_points( self, fastscanRootFile, deltaNLLCutOff = None, verbose=False ):

        if deltaNLLCutOff is None:
            deltaNLLCutOff = self.deltaNLLCutOff

        if not isfile(fastscanRootFile):
            if Commands.is_test_mode():
                print '\n[TESTMODE] No file \'{0}\'; Returning some bogus accepted points'.format(fastscanRootFile)
                return [ Container(iPoint=i) for i in xrange(5) ]
            else:
                Commands.throw_error( 'File \'{0}\' does not exist'.format(fastscanRootFile) )

        fastscanRootFp = ROOT.TFile.Open( fastscanRootFile )

        if not fastscanRootFp.GetListOfKeys().Contains('limit'):
            fastscanRootFp.Close()
            Commands.throw_error( 'There is no tree \'limit\' in', fastscanRootFile )

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
    def parse_bestfit_command( self ):

        self.prepare_for_command_compilation( mode = 'postfit' )

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

        self.common_command_settings(cmd, taskName)

        # For the postfit: Hard-coded saving the workspace
        cmd.append( '--saveWorkspace' )

        self.postfitRootFileBasename = 'higgsCombine_{0}.MultiDimFit.mH125.root'.format( taskName )

        return cmd


    #____________________________________________________________________
    def parse_bestfit_command_old( self ):
        self.prepare_for_command_compilation( mode = 'postfit' )

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
                    Commands.throw_error( 'Queue \'{0}\' is not available on PSI'.format(queue) )
                if self.jobPriority != 0:
                    cmd.append(
                        '--job-mode psi --task-name {0} --sub-opts=\'-q {1} -p {2}\' '.format( taskName, queue, self.jobPriority ),
                        )
                else:
                    cmd.append(
                        '--job-mode psi --task-name {0} --sub-opts=\'-q {1}\' '.format( taskName, queue ),
                        )
            else:
                Commands.throw_error( 'Only jobs submitted from T3 are implemented now' )

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
            if Commands.is_test_mode():
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
