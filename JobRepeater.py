#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, re, sys, argparse
from os.path import *
from glob import glob
from cStringIO import StringIO
from math import ceil

sys.path.append('src')
import Commands
from Container import Container

import ROOT



########################################
# Main
########################################

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument( 'scandir', type=str, default='default' )
    parser.add_argument( '--test', action='store_true', help='boolean')
    args = parser.parse_args()

    if args.test:
        Commands.TestMode()

    scandir = args.scandir
    verbose = True


    # ======================================
    # Find list of jobs that need rescanning

    allRootFiles = glob( scandir + '/higgsCombine*.root' )

    rootFiles = []
    for rootFile in allRootFiles:

        with Capturing() as output:
            rootFp = ROOT.TFile.Open( rootFile )

        if not rootFp.GetListOfKeys().Contains( 'limit' ):
            rootFiles.append( rootFile )
            # if verbose: print '    No tree \'limit\' in \'{0}\''.format( rootFile )
        rootFp.Close()

    # print 'root files to rescan:'
    # print '  ' + '\n  '.join(rootFiles)



    # ======================================
    # Gather list of jobs that are currently running

    qstatOutput = Commands.executeCommand( 'qstat -u \'tklijnsm\'', captureOutput=True, ignoreTestmode=True )
    jobIDs = re.findall( r'\n(\d+)\s', qstatOutput )

    runningShFiles = []
    for jobID in jobIDs:
        qstatPerJob = Commands.executeCommand( 'qstat -j {0}'.format(jobID), captureOutput=True, ignoreTestmode=True )
        scriptFile = re.search( r'script_file:\s+(.*)\n', qstatPerJob ).group(1)
        runningShFiles.append( scriptFile.strip() )

    # print '\n\nSkipping currently running sh files:'
    # print '  ' + '\n  '.join( runningShFiles )
    # sys.exit()



    # ======================================
    # Find relevant .sh file for the failed job (try to loop only once)

    containers = []
    for rootFile in rootFiles:
        container = Container()
        container.rootFile = rootFile
        container.name     = re.search( r'higgsCombine(.*)\.MultiDimFit\.mH\d+\.root', rootFile ).group(1)
        container.shFile   = None
        containers.append(container)

    allShFiles = glob( scandir + '/*.sh' )

    # Pretty inefficient loop but not sure how to do this otherwise
    for shFile in allShFiles:
        with open( shFile, 'r' ) as shFp:
            fulltext = shFp.read()

        for container in containers:
            if container.name in fulltext:
                container.shFile = shFile
                # if verbose:
                #     print 'Matched:'
                #     print '  ', container.rootFile
                #     print '  ', container.shFile
                break

        del fulltext


    # print ''
    # for container in containers:
    #     print 'Found following shFile for ', container.rootFile
    #     print '  ', container.shFile


    # ======================================
    # Resubmit the jobs with number of points split over more files

    splitFactor = 3
    for container in containers:

        if abspath(container.shFile) in runningShFiles:
            print 'Skipping \'{0}\' (already running in job {1})'.format(
                basename(container.shFile),
                jobIDs[ runningShFiles.index(abspath(container.shFile)) ]
                )
            continue

        with open( container.shFile, 'r' ) as shFp:
            shText = shFp.read()

        # Determine first point and last point
        firstPoint = int( search( r'\-\-firstPoint\s(\d+)\s', shText ).group(1) )
        lastPoint  = int( search( r'\-\-lastPoint\s(\d+)\s', shText ).group(1) )

        nPointsOriginal = lastPoint - firstPoint + 1
        nPointsPerJob   = int(ceil( float(nPointsOriginal) / float(splitFactor) ))

        pointBoundaries = range( firstPoint, lastPoint, nPointsPerJob )
        if not pointBoundaries[-1] == lastPoint:
            pointBoundaries.append( lastPoint )


        nJobs = len(pointBoundaries)-1
        for iJob in xrange(nJobs):

            newShFile = container.shFile.replace( '.sh', '_subsplit{0}.sh'.format(iJob) )

            newFirstPoint = pointBoundaries[iJob]

            if iJob == nJobs-1: # Only do the last bound inclusive
                newLastPoint  = pointBoundaries[iJob+1]
            else:
                newLastPoint  = pointBoundaries[iJob+1]-1

            newShText = shText
            newShText = newShText.replace( '--firstPoint {0}'.format( firstPoint ), '--firstPoint {0}'.format( newFirstPoint ) )
            newShText = newShText.replace( '--lastPoint {0}'.format( lastPoint ), '--lastPoint {0}'.format( newLastPoint ) )
            newShText = newShText.replace( 'POINTS.{0}.{1}'.format( firstPoint, lastPoint ), 'POINTS.{0}.{1}'.format( newFirstPoint, newLastPoint ) )

            if Commands.IsTestMode():
                print '{0}Would now open new shFile {1} with following contents:{2}'.format( '\033[94m', newShFile, '\033[0m' )
                print newShText
            else:
                with open( newShFile, 'w' ) as newShFp:
                    newShFp.write( newShText )
                Commands.executeCommand( 'qsub -q all.q {0}'.format(newShFile) )


class Capturing(list): # <- Inherits from list, so isinstance( self, list ) == True
    
    def __init__(self):

        # Open some fd's to /dev/null
        self.null_fd_out = os.open( os.devnull, os.O_RDWR )
        self.null_fd_err = os.open( os.devnull, os.O_RDWR )

        # Save original fd's:
        self.save_fd_out = os.dup(1)
        self.save_fd_err = os.dup(2)



    def __enter__(self):

        # Assign the null pointers to stdout and stderr.
        os.dup2( self.null_fd_out, 1 )
        os.dup2( self.null_fd_err, 2 )


        # Save original file descriptors (this is a copy operation)
        self._stdout   = sys.stdout
        self._stderr   = sys.stderr

        # Open new file-descriptor-like strings for stdout and stderr
        sys.stdout     = StringIO()
        sys.stderr     = StringIO()

        # Save these to the class
        self._stringio_out = sys.stdout
        self._stringio_err = sys.stderr

        return self


    def __exit__(self, *args):

        # Remember self is a list
        self.extend( self._stringio_out.getvalue().splitlines() )
        self.extend( self._stringio_err.getvalue().splitlines() )

        # Be nice and free up the now useless variables
        del self._stringio_out
        del self._stringio_err

        # Restore regular file descriptors
        sys.stdout = self._stdout
        sys.stderr = self._stderr


        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2( self.save_fd_out , 1 )
        os.dup2( self.save_fd_err , 2 )
        # Close the null files
        os.close( self.null_fd_out )
        os.close( self.null_fd_err )


def search( pat, text ):
    match = re.search( pat, text )
    if not match:
        print '[ERROR] Could not match pattern \'{0}\' to text \'{1}\''.format( pat, text )
        sys.exit()
    return match




def old():

    jobToRepeat = 'manual_Scan_coupling_Aug09_1/job_SCAN_kappackappab_Aug09_combinedCard_Jul26_CouplingModel_Yukawa_withTheoryUncertainties_27.sh'


    with open( jobToRepeat, 'r' ) as shFp:
        fulltext = shFp.read()
    lines = [ line for line in fulltext.split('\n') if not len(line) == 0 ]

    print fulltext
    print '\n'
    print lines
    print '\n\n------------------------------'

    # Determine the job number
    jobNumber = search( r'_(\d+).sh', jobToRepeat ).group(1)

    print jobNumber

    # Determine first point and last point
    firstPoint = int( search( r'\-\-firstPoint\s(\d+)\s', fulltext ).group(1) )
    lastPoint  = int( search( r'\-\-lastPoint\s(\d+)\s', fulltext ).group(1) )

    print firstPoint
    print lastPoint



    newNPointsPerJob = 4
    newPointBoundaries = range( firstPoint, lastPoint, newNPointsPerJob ) + [ lastPoint ]
    if newPointBoundaries[-1] == newPointBoundaries[-2]:
        newPointBoundaries = newPointBoundaries[:-1]


    print newPointBoundaries


    nNewJobs = len(newPointBoundaries) - 1
    oldStdoutdir = search( r'#\$ -o ([\/\_\w]+)(\s)', fulltext ).group(1)

    for iJob in xrange(nNewJobs):

        newFilename = jobToRepeat.replace( '{0}.sh'.format(jobNumber), '{0}_{1}.sh'.format( jobNumber, iJob ) )
        stdoutDir = abspath( dirname( newFilename ) )

        newFirstPoint = newPointBoundaries[iJob]
        newLastPoint  = newPointBoundaries[iJob+1]

        newFulltext = fulltext

        newFulltext = newFulltext.replace( oldStdoutdir, stdoutDir )
        if oldStdoutdir.endswith('/'):
            newFulltext = newFulltext.replace( oldStdoutdir[:-1], stdoutDir )

        newFulltext = newFulltext.replace( '--firstPoint {0}'.format( firstPoint ), '--firstPoint {0}'.format( newFirstPoint ) )
        newFulltext = newFulltext.replace( '--lastPoint {0}'.format( lastPoint ), '--lastPoint {0}'.format( newLastPoint ) )

        newFulltext = newFulltext.replace( 'POINTS.{0}.{1}'.format( firstPoint, lastPoint ), 'POINTS.{0}.{1}'.format( newFirstPoint, newLastPoint ) )


        print '\n-------------------------'
        print newFilename
        print newFulltext

        with open( newFilename, 'w' ) as newFp:
            newFp.write( newFulltext )

        Commands.executeCommand( 'qsub -q all.q {0}'.format(newFilename) )







########################################
# End of Main
########################################
if __name__ == "__main__":
    main()