#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, tempfile, shutil, re
from os.path import *

from time import strftime
datestr = strftime( '%b%d' )

import ROOT

########################################
# Main
########################################

TESTMODE = False
def TestMode( flag=True ):
    global TESTMODE
    TESTMODE = flag

TEMPJOBDIR = 'tmp'



def BasicT2WS(
    datacard
    ):

    outputDir = abspath( 'workspaces_{0}'.format(datestr) )
    if not isdir( outputDir ): os.makedirs( outputDir )
    outputWS = join( outputDir, basename(datacard).replace( '.txt', '.root' ) )

    processes, bins = ListProcesses( datacard )

    cmd = []
    cmd.append( 'text2workspace.py' )
    cmd.append( datacard )
    cmd.append( '-o {0}'.format(outputWS) )
    cmd.append( '-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel' )
    cmd.append( '--PO verbose' )
    cmd.append( '--PO \'higgsMassRange=123,127\'' )

    parRange = [ -1.0, 3.0 ]
    # parName = lambda process: 'r_' + process.split('_',1)[1]
    parName = lambda process: 'r_' + process

    for process in processes:
        cmd.append( '--PO \'map=.*/{0}:{1}[1.0,{2},{3}]\''.format(
            process,
            parName( process ),
            parRange[0], parRange[1],
            ))

    executeCommand( cmd )



def BasicBestfit(
        datacard,
        setPOIs = True
        ):

    cmd = [
        'combine',
        datacard,
        '-M MultiDimFit',
    
        '-m 125',
        '--floatOtherPOIs=1',
        # '--saveWorkspace',
        '--minimizerStrategy 2',
        ]

    if setPOIs:

        datacardFp = ROOT.TFile.Open( datacard )
        w = datacardFp.Get('w')
        POIlist = ROOT.RooArgList( w.set('POI') )

        parNames = []
        for i in xrange( POIlist.getSize() ):
            parName = POIlist[i].GetName()
            if parName.startswith('r_'):
                parNames.append( parName + '=1.0' )

        cmd.append( '--setPhysicsModelParameters ' + ','.join(parNames) )


    executeCommand( cmd )



def ListPOIs( datacardRootFile ):

    datacardFp = ROOT.TFile.Open( datacardRootFile )
    w = datacardFp.Get('w')
    POIlist = ROOT.RooArgList( w.set('POI') )

    parNames = []
    for i in xrange( POIlist.getSize() ):
        parName = POIlist[i].GetName()
        if parName.startswith('r_'):
            parNames.append( parName )

    datacardFp.Close()
    return parNames




def ListProcesses( datacardFile ):

    # shapes
    # smH_PTH_30p0_45p0 hgg_PTH_0p0_15p0_SigmaMpTTag0
    # CMS-HGG_sigfit_differential_pT_moriond17_reminiaod_extrabin_corrections_newsysts_v5.root
    # wsig_13TeV:hggpdfsmrel_13TeV_InsideAcceptance_genPt_30p0to45p0_SigmaMpTTag_0_recoPt_0p0to15p0

    with open( datacardFile, 'r' ) as datacardFp:
        lines = [ i.strip() for i in datacardFp.readlines() if i.strip().startswith('shapes') ]

    allProcesses = list(set([ l.split()[1] for l in lines ]))
    allProcesses.sort()

    processes = []
    for process in allProcesses:
        components = process.split('_')
        if not ( len(components) == 3 or len(components) == 4 ): continue
        if not components[0].endswith('H'): continue
        processes.append( process )


    bins = list(set([ l.split()[2] for l in lines ]))
    bins.sort()

    return processes, bins




def executeCommand( cmd ):

    if not isinstance( cmd, basestring ):
        cmdStr = '\n    '.join( cmd )
        cmdExec = ' '.join(cmd)
    else:
        cmdStr = cmd
        cmdExec = cmd

    if TESTMODE:
        print '\nTESTMODE: ' + cmdStr
    else:
        print '\nEXECUTING: ' + cmdStr
        os.system( cmdExec )




def executeCommandOnBatch( cmd ):
    print 'Not yet implemented'
    return
    # global TEMPJOBDIR
    # if not isdir(TEMPJOBDIR): os.makedirs(TEMPJOBDIR)



def BasicCombineTool(
        datacard,
        POIpattern    = '*',
        POIRange      = [ -1.0, 3.0 ],
        nPoints       = 100,
        nPointsPerJob = 3,
        queue         = '1nh',
        notOnBatch    = False,
        jobDirectory  = None,
        ):

    datacard = abspath( datacard )

    # cmd = 'combineTool.py '
    # cmd += rootfile
    # cmd += ' -n {0}'.format(scanName) 
    # cmd += ' -M MultiDimFit ' 
    # cmd += ' --cminDefaultMinimizerType Minuit2 ' 
    # cmd += ' --cminDefaultMinimizerAlgo migrad ' 
    # cmd += ' --algo=grid  ' 
    # cmd += ' --floatOtherPOIs=1 ' 
    # cmd += ' -P "{0}" '.format(POI) 
    # cmd += ' --setPhysicsModelParameterRanges "{0}"={1:.3f},{2:.3f} '.format( POI, Range[0], Range[1]) 
    # cmd += ' -m 125.00 ' 
    # cmd += ' --squareDistPoi ' 
    # cmd += ' --saveNLL ' 
    # cmd += ' --saveInactivePOI 1 ' 
    # cmd += ' --saveSpecifiedFunc {0} '.format(saveFunctions) 
    # cmd += ' --points={0} '.format(nPoints) 
    # cmd += ' --split-points {0} '.format(nPointsPerJob) 
    # cmd += ' --job-mode lxbatch --task-name {0} --sub-opts=\'-q 1nh\' '.format(scanName) 

    scanName = basename(datacard).replace('.root','')
    scanName = re.sub( r'\W', '', scanName )
    scanName = 'SCAN_{0}_{1}'.format( datestr, scanName )

    POIs = ListPOIs( datacard )
    if not POIpattern == '*':
        allPOIs = POIs[:]
        POIs = []
        for iterPOI in allPOIs:
            if re.search( POIpattern, iterPOI ):
                POIs.append( iterPOI )
        if len(POIs) == 0:
            print 'ERROR: Pattern \'{0}\' does not match any of the available POIs in the workspace:'.format( POIpattern )
            print '\n'.join( allPOIs )
            return

    currentdir = os.getcwd()
    if not TESTMODE:
        if not jobDirectory:
            jobDirectory = TEMPJOBDIR
        if not isdir( jobDirectory ):
            print 'Creating directory {0}'.format( jobDirectory )
            os.makedirs( jobDirectory )
        print 'Moving to directory {0}'.format( jobDirectory )
        os.chdir( jobDirectory )


    for POI in POIs:

        cmd = [
            'combineTool.py',
            datacard,
            '-n {0}'.format( scanName ),
            '-M MultiDimFit',
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            '--algo=grid',
            '--floatOtherPOIs=1',
            '--setPhysicsModelParameterRanges "{0}"={1:.3f},{2:.3f} '.format( POI, POIRange[0], POIRange[1] ),
            '-m 125.00',
            '--squareDistPoi',
            '--saveNLL',
            '--saveInactivePOI 1',
            '--points={0} '.format(nPoints),
            ]

        if not notOnBatch:
            cmd.append(
                '--split-points {0} '.format(nPointsPerJob)
                )

            if 't3' in os.environ['HOSTNAME']:

                if not queue in [ 'all.q', 'long.q', 'short.q' ]:
                    print 'Queue \'{0}\' is not available on PSI'.format(queue)
                    return

                cmd.append(
                    '--job-mode SGE --task-name {0} --sub-opts=\'-q {1}\' '.format( scanName, queue ),
                    )

            else:

                if not queue in [ '8nm', '1nh', '8nh', '1nd', '2nd' ]:
                    print 'Queue \'{0}\' is not available on lxplus'.format(queue)
                    return

                cmd.append(
                    '--job-mode lxbatch --task-name {0} --sub-opts=\'-q {1}\' '.format( scanName, queue ),
                    )

        executeCommand( cmd )

    os.chdir( currentdir )





########################################
# Execute this if called directly
########################################

if __name__ == '__main__':
    print 'Import functions from this file'
    print 'Available functions:'
    functions = [ name for (name, thing) in locals().items() if callable(thing) ]
    print '\n'.join( functions )