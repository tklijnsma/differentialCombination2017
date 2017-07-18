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







# Returns a function that can be evaluated for an integral from a to b
def GetIntegralTrapzoidal( xs, ys ):

    nPoints = len(xs)

    linearInterpolate = lambda a, x1, x2, y1, y2: \
        y1 + ( a - x1 ) / ( x2 - x1 ) * ( y2 - y1 )

    # Function to return
    def integralfunction( a, b, verbose=False ):

        if verbose: print '\n\nInterpolation function called with a = {0} and b = {1}'.format( a, b )

        aInterpolated = False
        bInterpolated = False

        # Extrapolating
        if a < xs[0]:
            if verbose: print 'Warning: Have to extrapolate for integral (a={0}, min x={1})'.format( a, xs[0] )
            ya = linearInterpolate( a, xs[0], xs[1], ys[0], ys[1] )
            ia = -1
            aInterpolated = True
        if b > xs[-1]:
            if verbose: print 'Warning: Have to extrapolate for integral (b={0}, max x={1})'.format( b, xs[-1] )
            yb = linearInterpolate( b, xs[-2], xs[-1], ys[-2], ys[-1] )
            ib = nPoints
            bInterpolated = True

        # Checking if bounds are defined exactly in the lists
        if a in xs:
            ia = xs.index(a)
            ya = ys[ia]
            aInterpolated = True
        if b in xs:
            ib = xs.index(b)
            yb = ys[ib]
            ib -= 1 # Otherwise the last value is double later
            bInterpolated = True

        # Interpolate linearly where a and b are
        for iPoint in xrange(nPoints-1):
            if aInterpolated and bInterpolated: break

            if not aInterpolated and a < xs[iPoint+1]:
                ya = linearInterpolate( a, xs[iPoint], xs[iPoint+1], ys[iPoint], ys[iPoint+1] )
                ia = iPoint
                aInterpolated = True
                if verbose: print 'Interpolated a = {0}: ia = {1}, ya = {2}'.format( a, ia, ya )

            if not bInterpolated and b < xs[iPoint+1]:
                yb = linearInterpolate( b, xs[iPoint], xs[iPoint+1], ys[iPoint], ys[iPoint+1] )
                ib = iPoint
                bInterpolated = True
                if verbose: print 'Interpolated b = {0}: ib = {1}, yb = {2}'.format( b, ib, yb )


        if not aInterpolated or not bInterpolated:
            Commands.ThrowError( 'Somehow this case was not interpolated; this is wrong' )
            return 0.

        x_integration = [ a  ] + xs[ia+1:ib+1] + [ b  ]
        y_integration = [ ya ] + ys[ia+1:ib+1] + [ yb ]
        
        if verbose:
            print 'Integrating over points'
            print 'x values:'
            print '    ' + '\n    '.join( [ str(f) for f in  x_integration ] )
            print 'y values:'
            print '    ' + '\n    '.join( [ str(f) for f in  y_integration ] )

        # Do trapezoidal integration
        integral = 0.
        for iPoint in xrange(len(x_integration)-1):
            dx = x_integration[iPoint+1] - x_integration[iPoint]
            y1 = y_integration[iPoint]
            y2 = y_integration[iPoint+1]
            integral += dx * 0.5*( y1 + y2 )

        if verbose: print 'Found integral = {0}'.format( integral )

        return integral

    return integralfunction




def hzz_T2WS(
    datacard,
    extraOptions=None,
    ):

    outputDir = abspath( 'workspaces_{0}'.format(datestr) )
    if not isdir( outputDir ): os.makedirs( outputDir )
    outputWS = join( outputDir, basename(datacard).replace( '.txt', '.root' ) )

    signalprocesses, processes, bins = Commands.ListProcesses( datacard )
    # cats = list(set([ b.split('cat')[1] for b in bins ]))

    cmd = []
    cmd.append( 'text2workspace.py' )
    cmd.append( datacard )
    cmd.append( '-o {0}'.format(outputWS) )
    cmd.append( '-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel' )
    cmd.append( '--PO verbose' )
    cmd.append( '--PO \'higgsMassRange=123,127\'' )

    parRange = [ -1.0, 3.0 ]
    parName = lambda process: 'r_' + process

    for iProcess, process in enumerate(signalprocesses):

        # Dummy to create a RooRealVar for r_$PROCESS
        cmd.append( '--PO \'map=dummy_{parName}:{parName}[1.0,{down},{up}]\''.format(
            parName = parName(process),
            down = parRange[0],
            up = parRange[1],
            ))
    
        # Actual scaling expression that uses r_$PROCESS and the pre-created $PROCESS_$BIN_norm
        for bin in bins:

            cmd.append((
                '--PO \'map={binName}/{processName}:'
                '{scalingParName}=expr::{scalingParName}("@0*@1",'
                '{parName},{normName})\'').format(
                    binName        = bin,
                    processName    = process,
                    parName        = parName(process),
                    scalingParName = parName(process) + '_' + bin,
                    # normName       = '{0}_{1}_norm'.format( process, bin.replace('PTH_','').replace('GE','GT') )
                    normName       = '{0}_norm'.format( process )
                    )
                )


    if not extraOptions:
        pass
    elif isinstance( extraOptions, basestring ):
        cmd.append( extraOptions )
    else:
        cmd.extend( extraOptions )

    Commands.executeCommand( cmd )



def hzz_manualMapping_T2WS(
    datacard,
    extraOptions=None,
    ):

    outputDir = abspath( 'workspaces_{0}'.format(datestr) )
    if not isdir( outputDir ): os.makedirs( outputDir )
    outputWS = join( outputDir, basename(datacard).replace( '.txt', '.root' ) )

    signalprocesses, processes, bins = ListProcesses( datacard )
    cats = list(set([ b.split('cat')[1] for b in bins ]))

    cmd = []
    cmd.append( 'text2workspace.py' )
    cmd.append( datacard )
    cmd.append( '-o {0}'.format(outputWS) )
    cmd.append( '-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel' )
    cmd.append( '--PO verbose' )
    cmd.append( '--PO \'higgsMassRange=123,127\'' )

    parRange = [ -1.0, 3.0 ]
    parName = lambda process: 'r_' + process
    fracParName = lambda process, cat: 'fracVar_{0}_cat{1}'.format( process, cat )


    # Create xs_hzz_smH_PTH_GT350 --> This was simply forgotten in the WS from David (now set equal to 200_350)
    cmd.append( '--PO \'map=dummy_xs_hzz_smH_PTH_GT350:xs_hzz_smH_PTH_GT350[0.0710184,0.0710184,0.0710184]\'' )

    for iProcess, process in enumerate(signalprocesses):

        # cmd += "--PO 'map=dummySqrtTau:sqrtTau[0.1]' "
        # cmd += "--PO 'map=.*/InsideAcceptance_genNjets2p5_m0p5to0p5:r0=expr::r0(\"@0+@1\",r0u,deltar0)' "

        # (SigmaBin0,fracSM4eBin0,K1Bin0,xs_hzz_smH_PTH_0_15) formula="(@0*@3*@1*@2)"

        # Dummy to create a RooRealVar for r_$PROCESS
        cmd.append( '--PO \'map=dummy_{parName}:{parName}[1.0,{down},{up}]\''.format(
            parName = parName(process),
            down = parRange[0],
            up = parRange[1],
            ))

        # Create expression for the fraction
        cmd.append((
            '--PO \'map=dummy_{parName}:'
            '{parName}=expr::{parName}("@0*@1",{fracSM4e},{K1})\'').format(
                parName  = fracParName(process,'4e'),
                fracSM4e = 'fracSM4eBin{0}'.format(iProcess),
                K1       = 'K1Bin{0}'.format(iProcess),
                )
            )
        cmd.append((
            '--PO \'map=dummy_{parName}:'
            '{parName}=expr::{parName}('
            '"((1.0-{fracSM4e}*{K1})*(1.0-{K2}*{fracSM4mu}/(1.0-{fracSM4e})))",'
            '{fracSM4e},{fracSM4mu},{K1},{K2})\'').format(
                parName   = fracParName(process,'2e2mu'),
                fracSM4e  = 'fracSM4eBin{0}'.format(iProcess),
                fracSM4mu = 'fracSM4muBin{0}'.format(iProcess),
                K1        = 'K1Bin{0}'.format(iProcess),
                K2        = 'K2Bin{0}'.format(iProcess),
                )
            )
        cmd.append((
            '--PO \'map=dummy_{parName}:'
            '{parName}=expr::{parName}('
            '"((1.0-{fracSM4e}*{K1})*({K2}*{fracSM4mu}/(1.0-{fracSM4e})))",'
            '{fracSM4e},{fracSM4mu},{K1},{K2})\'').format(
                parName   = fracParName(process,'4mu'),
                fracSM4e  = 'fracSM4eBin{0}'.format(iProcess),
                fracSM4mu = 'fracSM4muBin{0}'.format(iProcess),
                K1        = 'K1Bin{0}'.format(iProcess),
                K2        = 'K2Bin{0}'.format(iProcess),
                )
            )
    
        # Actual scaling expression that uses r_$PROCESS
        for cat in cats:

            cmd.append((
                '--PO \'map=.*{cat}.*/{processName}:'
                '{scalingParName}=expr::{scalingParName}("@0*@1*@2",'
                '{parName},{xsName},{fractionName})\'').format(
                    cat            = cat,
                    processName    = process,
                    parName        = parName(process),
                    scalingParName = parName(process) + '_cat' + cat,
                    xsName         = 'xs_hzz_{0}'.format(process),
                    fractionName   = fracParName( process, cat ),
                    )
                )

        # cmd.append( '--PO \'map=.*/{0}:{1}[1.0,{2},{3}]\''.format(
        #     process,
        #     parName( process ),
        #     parRange[0], parRange[1],
        #     ))

    if not extraOptions:
        pass
    elif isinstance( extraOptions, basestring ):
        cmd.append( extraOptions )
    else:
        cmd.extend( extraOptions )

    executeCommand( cmd )









########################################
# Execute this if called directly
########################################

if __name__ == '__main__':
    print 'Import functions from this file'
    print 'Available functions:'
    functions = [ name for (name, thing) in locals().items() if callable(thing) ]
    print '\n'.join( functions )