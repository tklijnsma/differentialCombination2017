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

import logging
import differentials
import differentials.core as core


########################################
# Main
########################################

class CombineConfig(object):

    def __init__(self, args, *otherargs, **kwargs):
        super(CombineConfig, self).__init__()

        self.SQUARE_DIST_POI_STEP        = True
        self.FLOAT_OTHER_POIS            = True
        self.ALGO                        = 'grid'
        self.SAVE_WORKSPACE              = False

        # Very basic, unlikely to ever change
        self.METHOD                      = 'MultiDimFit'
        self.DEFAULT_MASS                = 125.
        self.SAVE_NLL                    = True
        self.SAVE_INACTIVE_POI           = True

        # Default settings for jobs
        self.fromPostfit                 = False
        self.onBatch                     = True
        self.queue                       = 'all.q'
        self.jobPriority                 = 0

        self.datacard                    = 'somedatacard.root'
        self.subDirectory                = ''

        self.POIs                        = []
        self.PhysicsModelParameters      = []
        self.hardPhysicsModelParameters  = []
        self.PhysicsModelParameterRanges = []
        self.floatNuisances              = []
        self.freezeNuisances             = []

        self.extraOptions                = []

        self.nPoints                     = 100
        self.nPointsPerJob               = 8

        self.minimizer_settings = [
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            ]

        # Save the command line args in the config
        self.args = args
        self.tags = []

        if args.asimov:
            self.asimov = True
        else:
            self.asimov = False


    def make_unique_directory(self):
        for tag in self.tags:
            self.subDirectory += '_' + tag
        if self.asimov:
            self.subDirectory += '_asimov'
        self.subDirectory = core.make_unique_directory(self.subDirectory)

    def get_name(self):
        return basename(self.datacard).replace('.root','')

    def fix_parameter_at_value(self, name, val):
        self.PhysicsModelParameters.append('{0}={1}'.format(name, val))
        self.hardPhysicsModelParameters.append('{0}={1}'.format(name, val))
        self.freezeNuisances.append(name)


#____________________________________________________________________
class BaseCombineScan(object):

    debug = True
    verbose = True

    def __init__(self, config):
        self.input = deepcopy(config)

        # Attributes that may need to be altered from call to call
        self.onBatch = self.input.onBatch
        self.datacard = abspath(self.input.datacard)
        self.nPoints = self.input.nPoints
        self.nPointsPerJob = self.input.nPointsPerJob
        self.subDirectory = self.input.subDirectory
        self.extraOptions = self.input.extraOptions

        self.set_only_hard_physics_model_parameters = False

        self.freezeNuisances = self.input.freezeNuisances[:]

    def get_task_name(self):
        return '_UNSPECIFIED_' + ( 'ASIMOV_' if self.input.asimov else '' ) + self.input.get_name()

    def print_debug(self, txt):
        if self.debug:
            print '[CTW DEBUG]',txt

    def print_info(self, txt):
        if self.verbose:
            print '[CTW INFO]',txt


    def get_output(self):
        with core.enterdirectory(self.subDirectory, verbose=False):
            if core.is_testmode():
                output = 'higgsCombine_[TESTMODE].root'
            else:
                output = 'higgsCombine{0}.{1}.mH{2}.root'.format(self.get_task_name(), self.input.METHOD, int(self.input.DEFAULT_MASS))
                output = abspath(output)
        return output

    def get_parameter_settings(self):
        cmd = []
        if len(self.input.POIs) > 0:
            cmd.append( '-P ' + ' -P '.join(self.input.POIs) )
        if len(self.input.PhysicsModelParameterRanges) > 0:
            cmd.append( '--setPhysicsModelParameterRanges ' + ':'.join(self.input.PhysicsModelParameterRanges) )
        if len(self.input.floatNuisances) > 0:
            cmd.append( '--floatNuisances ' + ','.join(self.input.floatNuisances) )
        if len(self.freezeNuisances) > 0:
            cmd.append( '--freezeNuisances ' + ','.join(self.freezeNuisances) )
        cmd.extend(self.set_physics_model_parameters())
        return cmd

    def set_physics_model_parameters(self):
        cmd = []
        if not(self.input.asimov) and self.set_only_hard_physics_model_parameters:
            if len(self.input.hardPhysicsModelParameters) > 0:
                cmd.append('--setPhysicsModelParameters ' + ','.join(self.input.hardPhysicsModelParameters))
        else:
            # Do all parameters by default; overwrite with only hard in certain cases
            allPhysicsModelParameters = self.input.PhysicsModelParameters + self.input.hardPhysicsModelParameters
            if len(allPhysicsModelParameters) > 0:
                cmd.append('--setPhysicsModelParameters ' + ','.join(allPhysicsModelParameters))
        return cmd

    def parse_command(self):
        taskName = self.get_task_name()

        cmd = []

        if self.onBatch:
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
                    raise RuntimeError('Queue \'{0}\' is not available on PSI'.format(self.input.queue))
                if self.input.jobPriority != 0:
                    cmd.append(
                        '--job-mode psi --task-name {0} --sub-opts=\'-q {1} -p {2}\' '
                        .format( taskName, self.input.queue, self.input.jobPriority ),
                        )
                else:
                    cmd.append(
                        '--job-mode psi --task-name {0} --sub-opts=\'-q {1}\' '
                        .format( taskName, self.input.queue ),
                        )
            else:
                raise NotImplementedError( 'Only jobs submitted from T3 are implemented now' )

        # Do POIs, ranges, etc.
        cmd.extend(self.get_parameter_settings())

        return cmd


    def execute_command(self, cmd):
        if self.onBatch:
            output = core.execute(cmd, py_capture_output=True)
            logging.info('Output of cmd {0}'.format(cmd))
            logging.info(output)
        else:
            core.execute(cmd)
            output = ''
        if core.is_testmode():
            output = '\nOUTPUT: some output but this is testmode'
        return output

    def run(self):
        cmd = self.parse_command()
        with core.enterdirectory(self.subDirectory):
            output = self.execute_command(cmd)
        self.register_jobids_in_jobmanager(output)

    def register_jobids_in_jobmanager(self, submission_output):
        if core.is_testmode():
            print '[TESTMODE] Not writing any jobmanager files'
            return
        if not self.onBatch:
            print '{0} was not on batch; not registering jobs.'.format(self)
            return

        # Your job 8086766 ("job__SCAN_ASIMOV_hgg_Top_reweighted_nominal_148_0.sh") has been submitted
        jobids = re.findall(r'Your job (\d+)', submission_output)
        if len(jobids) == 0:
            print '\nNo jobids were found in the passed submission output; nothing to register for the jobmanager'

        header = [
            basename(self.subDirectory),
            'datacard: {0}'.format(self.datacard),
            'subDirectory: {0}'.format(self.subDirectory),
            'registration time: {0}'.format(strftime('%y-%m-%d %H:%M:%S')),
            'example cmd:\n\n{0}\n'.format(
                '\n    '.join(self.parse_command())
                )
            ]
        contents = '\n'.join(header) + '\n' + '\n'.join(jobids) + '\n'

        _, jobman_file = tempfile.mkstemp(
            prefix = 'tklijnsm_queuegroup_',
            suffix = '.jobman',
            dir    = '/tmp'
            )
        print 'Dumping following jobmanager contents to {0}:\n'.format(jobman_file)
        print contents + '\n'
        with open(jobman_file, 'w') as jobman_fp:
            jobman_fp.write(contents)


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
            raise ValueError(
                'PhysicsModelParameters *have* to be set when running asimov, '
                'or the best fit may make no sense'
                )

        # For the postfit: Hard-coded saving the workspace
        cmd.append( '--saveWorkspace' )

        return cmd

#____________________________________________________________________
class CombineOnlyPostfit(BaseCombineScan):

    def __init__(self, *args, **kwargs):
        super(CombineOnlyPostfit, self).__init__(*args, **kwargs)
        self.subDirectory = 'out/postfits_{0}'.format(core.datestr())
        self.onBatch = True

    def get_task_name(self):
        return '_POSTFIT_' + ( 'ASIMOV_' if self.input.asimov else '' ) + self.input.get_name()

    def parse_command(self):
        cmd = super(CombineOnlyPostfit, self).parse_command()

        if self.input.asimov and len(self.input.PhysicsModelParameters)==0:
            raise ValueError(
                'PhysicsModelParameters *have* to be set when running asimov, '
                'or the best fit may make no sense'
                )

        # For the postfit: Hard-coded saving the workspace
        cmd.append( '--saveWorkspace' )

        return cmd

#____________________________________________________________________
class CombineCorrMat(BaseCombineScan):

    def __init__(self, *args, **kwargs):
        super(CombineCorrMat, self).__init__(*args, **kwargs)
        self.subDirectory = 'out/corrmats_{0}'.format(core.datestr())

    def get_task_name(self):
        return '_CORRMAT_' + ( 'ASIMOV_' if self.input.asimov else '' ) + self.input.get_name()

    def parse_command(self):
        cmd = super(CombineCorrMat, self).parse_command()
        if self.input.asimov and len(self.input.PhysicsModelParameters)==0:
            raise ValueError(
                'PhysicsModelParameters *have* to be set when running asimov, '
                'or the best fit may make no sense'
                )

        if core.is_testmode():
            pdf_vars_to_freeze = [ 'some', 'pdfs' ]
        else:
            # pdf_vars_to_freeze = ListOfPDFIndicesToFreeze( postfitFilename, verbose=False )
            pdf_vars_to_freeze = differentials.pdffreezer.PDFFreezer(self.datacard).get_vars_to_freeze()

        cmd.extend([
            '--algo none',
            '--snapshotName MultiDimFit',
            '--saveWorkspace',
            # '--skipInitialFit',
            '--computeCovarianceMatrix=1',
            '--freezeNuisances {0}'.format( ','.join(pdf_vars_to_freeze) ),
            ])

        return cmd

#____________________________________________________________________
class CombineScan(BaseCombineScan):
    def __init__(self, *args, **kwargs):
        super(CombineScan, self).__init__(*args, **kwargs)

    def get_task_name(self):
        return '_SCAN_' + ( 'ASIMOV_' if self.input.asimov else '' ) + self.input.get_name()

    def parse_command(self):
        cmd = super(CombineScan, self).parse_command()
        cmd.extend([
            '--algo=grid',
            '--points={0}'.format(self.nPoints),
            '--split-points {0}'.format(self.nPointsPerJob)
            ])
        return cmd

#____________________________________________________________________
class CombineScanFromPostFit(BaseCombineScan):
    def __init__(self, *args, **kwargs):
        super(CombineScanFromPostFit, self).__init__(*args, **kwargs)
        self.set_only_hard_physics_model_parameters = True

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
        self.set_only_hard_physics_model_parameters = True

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
        self.set_only_hard_physics_model_parameters = True
        self.onBatch = True

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
        return cmd

    def run(self, postfitWS, fastscanFile):
        self.datacard = postfitWS
        self.fastscanFile = fastscanFile

        # Base extraOptions
        _extraOptions = self.extraOptions

        accepted_points = self.list_accepted_points(self.fastscanFile)

        submission_outputs = ''
        with core.enterdirectory(self.subDirectory):
            for iChunk, chunk in enumerate(chunks(accepted_points, self.nPointsPerJob)):
                print '\nJob', iChunk
                self.extraOptions = _extraOptions + [ '--doPoints ' + ','.join([ str(i) for i in chunk ]) ]
                self.get_task_name = lambda: self.get_task_name_without_number() + '_' + str(iChunk)
                cmd = self.parse_command()
                output = self.execute_command(cmd)
                submission_outputs += '\n' + output
        self.register_jobids_in_jobmanager(submission_outputs)


    def list_accepted_points(self, fastscanFile):

        if not isfile(fastscanFile):
            if core.is_testmode():
                self.print_info('[TESTMODE] No file \'{0}\'; Returning some bogus accepted points'.format(fastscanFile))
                # return [ Container(iPoint=i) for i in xrange(2*self.nPointsPerJob) ]
                return range(2*self.nPointsPerJob)
            else:
                raise ValueError('File \'{0}\' does not exist'.format(fastscanFile))


        with core.openroot(fastscanFile) as fastscanFp:
            if not fastscanFp.GetListOfKeys().Contains('limit'):
                raise ValueError('There is no tree \'limit\' in', fastscanFile)

            acceptedPoints = []
            rejectedPoints = []
            tree = fastscanFp.Get('limit')
            for iEvent, event in enumerate(tree):
                container = core.AttrDict()
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
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]
