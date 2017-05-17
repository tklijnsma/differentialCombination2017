#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import Commands

import os, tempfile, shutil, re, glob, itertools
from os.path import *
from operator import itemgetter
from array import array

from time import strftime
datestr = strftime( '%b%d' )

import ROOT


########################################
# Commands
########################################



def ChangeOutsideAcceptanceToSignalProcess(
    datacard = 'suppliedInput/PTH/hzz4l_comb_13TeV_xs.txt',
    ):

    with open( datacard, 'r' ) as datacardFp:
        lines = datacardFp.readlines()

    # Find line with "process      <number>"


    for iLine, line in enumerate(lines):
        line = line.strip()
        if not line.startswith('process'): continue

        components = [ l.strip() for l in line.split() ]
        try:
            eval( components[1] )
        except:
            continue

        break

    iLineProcess = iLine
    print 'Found line with following components:'
    print components

    print 'Shifted list (white space unfortunately not preserved):'
    components = [ components[0] ] + [ str(int(c)-1) for c in components[1:] ]
    print components

    line = ' '.join(components) + '\n'
    print 'Inserting line:'
    print line


    # Insert in place of old process line
    lines[iLineProcess] = line

    datacardOut = datacard.replace( '.txt', '_processesShifted.txt' )
    with open( datacardOut, 'w' ) as datacardFp:
        datacardFp.write( ''.join(lines) )






def RenameNormalizationsInHzz():

    # Compile list of all faulty normalization names
    signalprocesses, processes, bins = Commands.ListProcesses( 'suppliedInput/PTH/hzz4l_comb_13TeV_xs.txt' )
    oldNames = []
    newNames = []
    for process, bin in itertools.product( signalprocesses, bins ):
        oldNames.append( '{0}_{1}_norm'.format( process, bin.replace( 'PTH_', '' ) ) )
        newNames.append( '{0}_{1}_norm'.format( process, bin ) )


    outdir = 'suppliedInput/PTH_normChanged'
    if not isdir(outdir): os.makedirs(outdir)

    for rootFile in glob( 'suppliedInput/PTH/*Databin*.root' ):
        outRootFile = join( outdir, basename(rootFile) )

        print 'Processing {0}'.format( rootFile )

        ChangeNameOfElementInWs(
            rootFile,
            outRootFile,
            oldNames,
            newNames,
            )

        break



def ChangeNameOfElementInWs(
    wsFile,
    wsFileOut,
    oldNames,
    newNames,
    wsName = 'w',
    ):

    assert( not isinstance( oldNames, basestring ) and not isinstance( newNames, basestring ) )
    if not len(oldNames) == len(newNames):
        print 'ERROR in ChangeNameOfElementInWs(): Supply lists of the same length'
        print '    oldNames = ', oldNames
        print '    newNames = ', newNames


    wsFp = ROOT.TFile.Open( wsFile )
    w = wsFp.Get(wsName)




    for varName, newVarName in zip( oldNames, newNames ):
        # print 'Trying to get \'{0}\''.format( varName )
        variable = GetVarFromWs( w, varName )

        if not variable:
            continue
        else:
            print 'Trying to import \'{0}\' as \'{1}\''.format( varName, newVarName )
            variable.SetName( newVarName )
            getattr( w, 'import' )( variable )


    # if not success:
    #     print 'Failure in ChangeNameOfElementInWs():'
    #     print '    Element \'{0}\' does not exist in \'{1}\''.format(
    #         oldName, wsFile
    #         )

    newWsFp = ROOT.TFile.Open( wsFileOut, 'RECREATE' )
    newWsFp.WriteTObject( w )
    newWsFp.Close()


    wsFp.Close()
    
    


def GetVarFromWs( w, varName ):
    success = False
    for getter in [ 'pdf', 'function', 'var', 'cat', 'genobj' ]:
        try:
            variable = getattr( w, getter )( varName )
            variable.Print()
            success = True
        except ReferenceError:
            continue
        break

    if not success:
        return None
    else:
        return variable


########################################
# Execute this if called directly
########################################

if __name__ == '__main__':
    print 'Import functions from this file'
    print 'Available functions:'
    functions = [ name for (name, thing) in locals().items() if callable(thing) ]
    print '\n'.join( functions )