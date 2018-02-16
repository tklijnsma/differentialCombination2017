#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, shutil, re, subprocess, sys, traceback, itertools
from os.path import *
from glob import glob
from datetime import datetime
import random

from time import strftime
datestr = strftime( '%b%d' )
datestr_detailed = strftime( '%y-%m-%d %H:%M:%S' )

src = os.path.basename(__file__)
from subprocess import check_output
currentcommit = check_output(['git', 'log', '-1', '--oneline' ])

from Container import Container

import ROOT

########################################
# Main
########################################

TESTMODE = False
def TestMode(flag=True):
    global TESTMODE
    TESTMODE = flag
def IsTestMode():
    global TESTMODE
    return TESTMODE

DISABLE_WARNINGS = False
def DisableWarnings(flag=True):
    global DISABLE_WARNINGS
    DISABLE_WARNINGS = flag

TEMPJOBDIR = abspath( 'tmp' )
def SetTempJobDir( newdirname='tmp' ):
    global TEMPJOBDIR
    TEMPJOBDIR = newdirname
def GetTempJobDir():
    global TEMPJOBDIR
    return TEMPJOBDIR

#____________________________________________________________________
def FormatStrToWidth(text, width):
    if len(text) > width:
        text = text[:width-3] + '...'
    ret = '{0:{width}}'.format( text, width=width )
    return ret

def PrintTable(table, minColWidth=1, maxColWidth=20, sep='  ', newline_sep='\n'):
    nRows = len(table)
    nCols = len(table[0])

    maxColWidths = []
    for iCol in xrange(nCols):
        maxWidth = 0
        for iRow in xrange(nRows):
            entry = table[iRow][iCol]
            if not isinstance( entry, basestring ):
                ThrowError( 'Entry of a table is not a string' )
            # entry = escape_ansi( entry )
            if len(entry) > maxWidth:
                maxWidth = len(entry)
        maxColWidths.append( min( maxWidth, maxColWidth ) )

    out = []
    for iRow in xrange(nRows):
        line = []
        for iCol in xrange(nCols):
            entry = table[iRow][iCol]
            line.append( FormatStrToWidth( entry, maxColWidths[iCol] ) )
        out.append( sep.join(line) )
    return newline_sep.join(out)

#____________________________________________________________________
def AppendNumberToDirNameUntilItDoesNotExistAnymore(
        dirName,
        nAttempts = 100,
        ):

    dirName = abspath(dirName)

    if not isdir(dirName):
        return dirName

    dirName += '_{0}'
    for iAttempt in xrange(nAttempts):
        if not isdir( dirName.format(iAttempt) ):
            dirName = dirName.format(iAttempt)
            break
    else:
        ThrowError( 'Could not create a unique directory for {0}'.format(dirName.format('X')) )
        sys.exit()

    print '[info] New directory: {0}'.format( dirName )
    return dirName




#____________________________________________________________________
def ConvertFloatToStr( number, nDecimals=None ):
    number = float(number)
    if not nDecimals is None:
        string = '{:.{nDecimals}f}'.format( number, nDecimals=nDecimals ).replace('-','m').replace('.','p')
        return string
    if number.is_integer():
        number = int(number)
    string = str(number).replace('-','m').replace('.','p')
    return string

#____________________________________________________________________
def ConvertStrToFloat( string ):
    string = str(string)
    number = string.replace('m','-').replace('p','.')
    number = float(number)
    return number

#____________________________________________________________________
def getRangeFromStr(text):
    regular_match = re.search(r'([\dpm\.\-]+)_([\dpm\.\-]+)', text)
    overflow_match = re.search(r'(GE|GT)([\dpm\.\-]+)', text)

    if regular_match:
        left = ConvertStrToFloat(regular_match.group(1))
        right = ConvertStrToFloat(regular_match.group(2))
    elif overflow_match:
        left = ConvertStrToFloat(overflow_match.group(2))
        right = 'INF'
    else:
        left = 'UNDEFINED'
        right = 'UNDEFINED'

    return left, right

def rangeSorter(text):
    left, right = getRangeFromStr(text)
    if left == 'UNDEFINED':
        return 900000
    elif right == 'INF':
        return 800000
    else:
        return left

def POIsorter( POI ):
    _1, _2, Range = InterpretPOI(POI)
    if Range[0] == '-INF':
        return -100000
    elif len(Range) > 1 and Range[1] == 'INF':
        return 900000
    else:
        return Range[0]

def SortPOIs( POIs ):
    POIs.sort( key = POIsorter )


#____________________________________________________________________
def executeCommand( cmd, captureOutput=False, ignoreTestmode=False ):

    if not isinstance( cmd, basestring ):
        cmd = [ l for l in cmd if not len(l.strip()) == 0 ]
        cmdStr = '\n    '.join( cmd )
        cmdExec = ' '.join(cmd)
    else:
        cmdStr = cmd
        cmdExec = cmd

    if TESTMODE and not ignoreTestmode:
        print '\n[TESTMODE] ' + cmdStr + '\n'
    else:
        if not captureOutput:
            print '\n[EXECUTING] ' + cmdStr + '\n'
            os.system( cmdExec )
        else:
            output = subprocess.check_output(
                cmd,
                shell=True,
                )
            return output


#____________________________________________________________________
def copyfile( src, dst, verbose=True ):
    # src = abspath(src)
    # dst = abspath(dst)
    if IsTestMode():
        print '\n[TESTMODE] Would now copy\n  {0}\n  to\n  {1}'.format( src, dst )
    else:
        if verbose: print '\n[EXECUTING] Copying\n  {0}\n  to\n  {1}'.format( src, dst )
        shutil.copyfile( src, dst )

#____________________________________________________________________
def movefile( src, dst, verbose=True ):
    # src = relpath( src, '.' )
    # dst = relpath( dst, '.' )
    if IsTestMode():
        print '\n[TESTMODE] Would now move\n  {0}\n  to\n  {1}'.format( src, dst )
    else:
        if verbose: print '\n[EXECUTING] Moving\n  {0}\n  to\n  {1}'.format( src, dst )
        os.rename( src, dst )

#____________________________________________________________________
def newColorCycle():
    return itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )




class AnalysisError(Exception):
    pass
#____________________________________________________________________
def ThrowError(
    errstr = '',
    throwException = True
    ):

    if throwException:
        raise AnalysisError( errstr )

    else:
        stack = traceback.extract_stack(None, 2)[0]
        linenumber = stack[1]
        funcname = stack[2]

        cwd = abspath( os.getcwd() )
        modulefilename = relpath( stack[0], cwd )

        print 'ERROR in {0}:{1} {2}:\n    {3}'.format( modulefilename, linenumber, funcname, errstr )

#____________________________________________________________________
def Warning(
        warningStr,
        ):

    if DISABLE_WARNINGS:
        return

    stack = traceback.extract_stack(None, 2)[0]
    linenumber = stack[1]
    funcname = stack[2]

    cwd = abspath( os.getcwd() )
    modulefilename = relpath( stack[0], cwd )

    print '\n[WARNING {0}:{1} L{2}] '.format(modulefilename, funcname, linenumber) + warningStr

#____________________________________________________________________
def TagGitCommitAndModule():

    stack = traceback.extract_stack(None, 2)
    stack = stack[0]
    linenumber = stack[1]
    funcname = stack[2]

    cwd = abspath( os.getcwd() )
    modulefilename = relpath( stack[0], cwd )

    ret = 'Generated on {0} by {1}; current git commit: {2}'.format( datestr_detailed, modulefilename, currentcommit.replace('\n','') )

    return ret


#____________________________________________________________________
class EnterDirectory():
    """Context manager to (create and) go into and out of a directory"""

    def __init__(self, subDirectory=None, verbose=True ):
        self.verbose = verbose
        self._active = False
        if not subDirectory is None and not subDirectory == '':
            self._active = True
            self.backDir = os.getcwd()
            self.subDirectory = subDirectory

    def __enter__(self):
        if self._active:
            if IsTestMode():
                if self.verbose: print '\n[TESTMODE] Would now create/go into \'{0}\''.format(self.subDirectory)
            else:
                if self.verbose: print ''
                if not isdir( self.subDirectory ):
                    if self.verbose: print 'Creating \'{0}\''.format( relpath( self.subDirectory, self.backDir ) )
                    os.makedirs( self.subDirectory )
                if self.verbose: print 'Entering \'{0}\''.format( relpath( self.subDirectory, self.backDir ) )
                os.chdir( self.subDirectory )
        return self

    def __exit__(self, *args):
        if self._active:
            os.chdir( self.backDir )

#____________________________________________________________________
class OpenRootFile():
    """Context manager to safely open and close root files"""

    def __init__(self, rootFile ):
        self._rootFile = rootFile

    def __enter__(self):
        if not isfile(self._rootFile):
            ThrowError( 'File {0} does not exist'.format(self._rootFile) )
        self._rootFp = ROOT.TFile.Open(self._rootFile)
        return self._rootFp

    def __exit__(self, *args):
        self._rootFp.Close()


#____________________________________________________________________
def fileno(file_or_fd):
    fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
    if not isinstance(fd, int):
        raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
    return fd


class RedirectStdout():
    """Context manager to capture stdout from ROOT/C++ prints"""

    def __init__( self, verbose=False ):
        self.stdout_fd = fileno(sys.stdout)
        self.enableDebugPrint = verbose
        self.isRedirected = False
        pass
        
    def __enter__( self ):
        self.debugPrint( 'Entering' )

        self.captured_fd_r, self.captured_fd_w = os.pipe()
        self.debugPrint( '  Opened read: {0}, and write: {1}'.format( self.captured_fd_r, self.captured_fd_w ) )

        # Copy stdout
        self.debugPrint( '  Copying stdout' )
        self.copied_stdout = os.fdopen( os.dup(self.stdout_fd), 'wb' )
        sys.stdout.flush() # To flush library buffers that dup2 knows nothing about
        
        # Overwrite stdout_fd with the target
        self.debugPrint( '  Overwriting target ({0}) with stdout_fd ({1})'.format( fileno(self.captured_fd_w), self.stdout_fd ) )        
        os.dup2( fileno(self.captured_fd_w), self.stdout_fd )
        self.isRedirected = True

        os.close( self.captured_fd_w )

        return self


    def __exit__(self, *args):
        sys.stdout.flush()
        os.dup2( self.copied_stdout.fileno(), self.stdout_fd )  # $ exec >&copied


    def read( self ):
        sys.stdout.flush()
        self.debugPrint( '  Draining pipe' )

        # Without this line the reading does not end - is that 'deadlock'?
        os.close(self.stdout_fd)

        captured_str = ''
        while True:
            data = os.read( self.captured_fd_r, 1024)
            if not data:
                break
            captured_str += data
            self.debugPrint( '\n  captured_str: ' + captured_str )

        self.debugPrint( '  Draining completed' )
        

        return captured_str


    def debugPrint( self, text ):
        if self.enableDebugPrint:
            if self.isRedirected:
                os.write( fileno(self.copied_stdout), text + '\n' )
            else:
                os.write( fileno(self.stdout_fd), text + '\n' )




#____________________________________________________________________
def GetCMSLabel(
        text='Preliminary',
        x=None,
        y=None,
        textSize   = 0.06,
        textOffset = None,
        drawImmediately=True
        ):

    if not text in [ 'Preliminary', 'Supplementary' ]:
        Commands.Warning( 'Label \'{0}\' is not a standard label!'.format(text) )
    l = ROOT.TLatex()
    ROOT.SetOwnership( l, False )
    l.SetNDC()
    l.SetTextAlign(11)
    l.SetTextFont(42)
    l.SetTextSize(textSize)
    latexStr = '#bf{{CMS}} #it{{#scale[0.75]{{{0}}}}}'.format(text)

    if textOffset is None:
        textOffset = 0.25 * textSize

    if not drawImmediately:
        return latexStr, l
    else:
        if x is None:
            x = ROOT.gPad.GetLeftMargin()
        if y is None:
            y = 1. - ROOT.gPad.GetTopMargin() + textOffset
        l.DrawLatex( x, y, latexStr )

#____________________________________________________________________
def GetCMSLumi(
        lumi=35.9,
        x=None,
        y=None,
        textSize=0.05,
        textOffset = None,
        drawImmediately=True
        ):

    l = ROOT.TLatex()
    ROOT.SetOwnership( l, False )
    l.SetNDC()
    l.SetTextAlign(31)
    l.SetTextFont(42)
    l.SetTextSize(textSize)
    latexStr = '{0:.1f} fb^{{-1}} (13 TeV)'.format(lumi)

    if textOffset is None:
        textOffset = 0.25 * textSize

    if not drawImmediately:
        return latexStr, l
    else:
        if x is None:
            x = 1. - ROOT.gPad.GetRightMargin()
        if y is None:
            y = 1. - ROOT.gPad.GetTopMargin() + textOffset
        l.DrawLatex( x, y, latexStr )

def __uniqueid__():
    mynow=datetime.now
    sft=datetime.strftime
    # store old datetime each time in order to check if we generate during same microsecond (glucky wallet !)
    # or if daylight savings event occurs (when clocks are adjusted backward) [rarely detected at this level]
    old_time=mynow() # fake init - on very speed machine it could increase your seed to seed + 1... but we have our contingency :)
    # manage seed
    seed_range_bits=14 # max range for seed
    seed_max_value=2**seed_range_bits - 1 # seed could not exceed 2**nbbits - 1
    # get random seed
    seed=random.getrandbits(seed_range_bits)
    current_seed=str(seed)
    # producing new ids
    while True:
        # get current time 
        current_time=mynow()
        if current_time <= old_time:
            # previous id generated in the same microsecond or Daylight saving time event occurs (when clocks are adjusted backward)
            seed = max(1,(seed + 1) % seed_max_value)
            current_seed=str(seed)
        # generate new id (concatenate seed and timestamp as numbers)
        #newid=hex(int(''.join([sft(current_time,'%f%S%M%H%d%m%Y'),current_seed])))[2:-1]
        newid=int(''.join([sft(current_time,'%f%S%M%H%d%m%Y'),current_seed]))
        # save current time
        old_time=current_time
        # return a new id
        yield newid


########################################
# Execute this if called directly
########################################

if __name__ == '__main__':
    print 'Import functions from this file'
    print 'Available functions:'
    functions = [ name for (name, thing) in locals().items() if callable(thing) ]
    print '\n'.join( functions )