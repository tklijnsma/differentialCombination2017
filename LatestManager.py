#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, sys, argparse
from os.path import *
from glob import glob

sys.path.append('src')
import Commands

import LatestPaths

from time import strftime
datestr = strftime( '%b%d' )


########################################
# Main
########################################

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument( 'pattern', type=str, default='*', nargs='?' )
    parser.add_argument( '--tarball', action='store_true' )
    parser.add_argument( '--test', action='store_true' )
    parser.add_argument( '--shrinkLogFiles', action='store_true' )
    # parser.add_argument( '--list', metavar='N', type=str, nargs='+', help='list of strings' )
    args = parser.parse_args()

    if args.test:
        Commands.test_mode()

    if args.tarball:
        TarballLatestPaths()

    if args.shrinkLogFiles:
        directories = [ d for d in glob( 'out/*' ) if isdir(d) ]
        for directory in directories:
            logfiles = glob( directory + '/*.o*' ) + glob( directory + '/*.e*' )
            for logfile in logfiles:
                ShrinkLogFile( logfile, HRToBytes( 100, 'mb' ) )

    else:
        for obj, inode in ListLatestPaths( args.pattern ):
            print '{0:56} : {1}'.format( obj, inode )


#____________________________________________________________________
def ListLatestPaths( pattern='*' ):

    allObjectsInLatestPaths = dir(LatestPaths)

    ret = []

    for obj in allObjectsInLatestPaths:

        if callable(obj):
            continue

        inode = getattr( LatestPaths, obj )

        if not isinstance( inode, basestring ):
            continue

        # if inode.startswith( 'suppliedInput' ):
        #     continue
        if inode.endswith( '.pyc' ):
            continue
        if not( isfile(inode) or isdir(inode) ):
            continue
        if not pattern == '*' and not pattern in obj + ' : ' + inode:
            continue

        if isdir(inode):
            log_files = glob( inode + '/*.o*' ) + glob( inode + '/*.e*' )
            for log_file in log_files:
                ShrinkLogFile(log_file)

        ret.append( ( obj, inode ) )

    return ret


#____________________________________________________________________
def TarballLatestPaths():

    LatestInodes = [ inode for _, inode in ListLatestPaths() ]

    tarballOut = join( os.getcwd(), 'tarball_out_{0}.tar'.format(datestr) )

    cmd = 'tar cf {0} {1}'.format( tarballOut, ' '.join(LatestInodes) )

    # Commands.test_mode()
    Commands.execute_command( cmd )


#____________________________________________________________________
def ShrinkLogFile( logFile, cutoff=None ):
    if cutoff is None: cutoff = HRToBytes( 10, 'mb' )
    if getsize(logFile) > cutoff:
        if Commands.is_test_mode():
            print 'Would now shrink {0} (is {1})'.format( logFile, bytesToHR(getsize(logFile)) )
        else:
            print 'Shrinking {0} (is {1})'.format( logFile, bytesToHR(getsize(logFile)) )
            beginning = Commands.execute_command( 'head -n 300 {0}'.format(logFile), captureOutput=True )
            end       = Commands.execute_command( 'tail -n 300 {0}'.format(logFile), captureOutput=True )
            centerTag = '\n<<<<<<<<<\nFILE LIMITED\n{0}\n>>>>>>>>>\n'.format( Commands.tag_git_commit_and_module() )
            newText = beginning[:-1] + '\n' + centerTag + '\n' + end
            with open( logFile, 'w' ) as logFp:
                logFp.write( newText )


#____________________________________________________________________
def bytesToHR( num, suffix='B' ):
    num = float(num)
    for unit in [ '', 'K', 'M', 'G' ]:
        if abs(num) < 1024.0:
            return "{0:.1f} {1}{2}".format( num, unit, suffix )
        num /= 1024.0
    return "{0:.1f} {1}{2}".format( num, 'T', suffix )


#____________________________________________________________________
def HRToBytes( num, unit ):
    num = float(num)
    conversionToBytes = {
        'kb' : 1024,
        'mb' : 1024*1024,
        'gb' : 1024*1024*1024,
        'tb' : 1024*1024*1024*2014,
        'k' : 1024,
        'm' : 1024*1024,
        'g' : 1024*1024*1024,
        't' : 1024*1024*1024*2014,
        }
    unit = unit.lower()
    if not unit in conversionToBytes:
        Tools.throw_error( 'Unit \'{0}\' has no known conversion factor to bytes'.format(unit) )
    return int( num*conversionToBytes[unit] )


########################################
# End of Main
########################################
if __name__ == "__main__":
    main()