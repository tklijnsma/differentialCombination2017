#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, tempfile, shutil, re, subprocess, sys, traceback
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

def IsTestMode():
    global TESTMODE
    return TESTMODE


TEMPJOBDIR = abspath( 'tmp' )


def BasicCombineCards(
    outputFile,
    *inputList
    ):

    cmd = [ 'combineCards.py' ]
    for datacard in inputList:
        cmd.append( datacard )

    cmd.append( '> {0}'.format( outputFile ) )

    executeCommand( cmd )



def BasicT2WS(
    datacard,
    extraOptions=None,
    ):

    outputDir = abspath( 'workspaces_{0}'.format(datestr) )
    if not isdir( outputDir ): os.makedirs( outputDir )
    outputWS = join( outputDir, basename(datacard).replace( '.txt', '.root' ) )

    signalprocesses, processes, bins = ListProcesses( datacard )

    cmd = []
    cmd.append( 'text2workspace.py' )
    cmd.append( datacard )
    cmd.append( '-o {0}'.format(outputWS) )
    cmd.append( '-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel' )
    cmd.append( '--PO verbose' )
    cmd.append( '--PO \'higgsMassRange=123,127\'' )

    parRange = [ -1.0, 3.0 ]
    parName = lambda process: 'r_' + process

    for process in signalprocesses:
        cmd.append( '--PO \'map=.*/{0}:{1}[1.0,{2},{3}]\''.format(
            process,
            parName( process ),
            parRange[0], parRange[1],
            ))

    if not extraOptions:
        pass
    elif isinstance( extraOptions, basestring ):
        cmd.append( extraOptions )
    else:
        cmd.extend( extraOptions )

    executeCommand( cmd )


def hzz_T2WS(
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
                    normName       = '{0}_{1}_norm'.format( process, bin.replace('PTH_','').replace('GE','GT') )
                    )
                )


    if not extraOptions:
        pass
    elif isinstance( extraOptions, basestring ):
        cmd.append( extraOptions )
    else:
        cmd.extend( extraOptions )

    executeCommand( cmd )



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



def BasicBestfit(
        datacard,
        setPOIs = True,
        onBatch = False,
        extraOptions = None,
        ):
    
    datacard = abspath( datacard )

    cmd = [
        'combine',
        datacard,
        '-M MultiDimFit',
        # '--saveWorkspace',
        # '--minimizerStrategy 2',
        # '-v 2',

        # '-m 125',
        # '--floatOtherPOIs=1',
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


    if not extraOptions:
        pass
    elif isinstance( extraOptions, basestring ):
        cmd.append( extraOptions )
    else:
        cmd.extend( extraOptions )


    if onBatch:

        if not TESTMODE:

            if not isdir(TEMPJOBDIR): os.makedirs(TEMPJOBDIR)
            osHandle, shFile =  tempfile.mkstemp(
                prefix = 'basicbestfitjob_',
                suffix = '.sh',
                dir = TEMPJOBDIR
                )
            shFile = abspath( shFile )

            cmsswPath = join( os.environ['CMSSW_BASE'], 'src' )

            with open( shFile, 'w' ) as shFp:
                shFp.write( '#$ -o {0} \n'.format( TEMPJOBDIR ) )
                shFp.write( '#$ -e {0} \n'.format( TEMPJOBDIR ) )
                shFp.write( 'export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch/ \n' )
                shFp.write( 'source /cvmfs/cms.cern.ch/cmsset_default.sh \n' )
                shFp.write( 'source /swshare/psit3/etc/profile.d/cms_ui_env.sh \n' )
                shFp.write( 'cd {0} \n'.format( cmsswPath ) )
                shFp.write( 'eval `scramv1 runtime -sh` \n')
                shFp.write( 'cd {0} \n'.format( TEMPJOBDIR ) )
                shFp.write( ' '.join(cmd) )

            qsubCmd = 'qsub -q short.q {0}'.format( shFile )
            executeCommand( qsubCmd )

        else:
            print 'TESTMODE + onBatch not implemented'

    else:
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

    # 09-05-2017: shapes lines may have stars - use process and bins lines instead

    # with open( datacardFile, 'r' ) as datacardFp:
    #     lines = [ i.strip() for i in datacardFp.readlines() if i.strip().startswith('shapes') ]

    # allProcesses = list(set([ l.split()[1] for l in lines ]))
    # allProcesses.sort()

    # processes = []
    # for process in allProcesses:
    #     components = process.split('_')
    #     if not ( len(components) == 3 or len(components) == 4 ): continue
    #     if not components[0].endswith('H'): continue
    #     processes.append( process )


    # bins = list(set([ l.split()[2] for l in lines ]))
    # bins.sort()

    # return processes, bins

    # with open( datacardFile, 'r' ) as datacardFp:
    #     lines = [ i.strip() for i in datacardFp.readlines() if i.strip().startswith('shapes') ]


    with open( datacardFile, 'r' ) as datacardFp:
        lines = [ i.strip() for i in datacardFp.readlines() ]

        for line in lines:
            if line.startswith('bin'):
                binline = line
                break

        for line in lines:
            if line.startswith('process'):
                processline = line
                break

    bins = binline.split()[1:]
    bins.sort()

    processes = list(set(processline.split()[1:]))
    processes.sort()

    signalprocesses = [ p for p in processes if p.split('_')[0].endswith('H') and not p == 'nonResH' ]

    return signalprocesses, processes, bins, 



def executeCommand( cmd, captureOutput=False ):

    if not isinstance( cmd, basestring ):
        cmdStr = '\n    '.join( cmd )
        cmdExec = ' '.join(cmd)
    else:
        cmdStr = cmd
        cmdExec = cmd

    if TESTMODE:
        print '\nTESTMODE: ' + cmdStr + '\n'
    else:
        if not captureOutput:
            print '\nEXECUTING: ' + cmdStr + '\n'
            os.system( cmdExec )
        else:
            print 'THIS IS UNTESTED'
            output = subprocess.check_output(
                cmd,
                shell=True,
                )
            print output





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

    scanName = basename(datacard).replace('.root','')
    scanName = re.sub( r'\W', '', scanName )
    scanName = 'SCAN_{{0}}_{0}_{1}'.format( datestr, scanName )

    allPOIs = ListPOIs( datacard )
    if not POIpattern == '*':
        doPOIs = []
        for iterPOI in allPOIs:
            if re.search( POIpattern, iterPOI ):
                doPOIs.append( iterPOI )
        if len(doPOIs) == 0:
            print 'ERROR: Pattern \'{0}\' does not match any of the available POIs in the workspace:'.format( POIpattern )
            print '\n'.join( allPOIs )
            return
    else:
        doPOIs = allPOIs[:]

    currentdir = os.getcwd()
    if not TESTMODE:
        if not jobDirectory:
            jobDirectory = TEMPJOBDIR
        if not isdir( jobDirectory ):
            print 'Creating directory {0}'.format( jobDirectory )
            os.makedirs( jobDirectory )
        print 'Moving to directory {0}'.format( jobDirectory )
        os.chdir( jobDirectory )


    for POI in doPOIs:

        cmd = [
            'combineTool.py',
            datacard,
            '-n {0}'.format( scanName.format(POI) ),
            '-M MultiDimFit',
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            '--algo=grid',
            '--floatOtherPOIs=1',
            '-P "{0}"'.format( POI ),
            '--setPhysicsModelParameterRanges "{0}"={1:.3f},{2:.3f} '.format( POI, POIRange[0], POIRange[1] ),
            '--setPhysicsModelParameters {0}'.format( ','.join([ iterPOI + '=1.0' for iterPOI in allPOIs ]) ),
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
                    '--job-mode psi --task-name {0} --sub-opts=\'-q {1}\' '.format( scanName.format(POI), queue ),
                    )

            else:

                if not queue in [ '8nm', '1nh', '8nh', '1nd', '2nd' ]:
                    print 'Queue \'{0}\' is not available on lxplus'.format(queue)
                    return

                cmd.append(
                    '--job-mode lxbatch --task-name {0} --sub-opts=\'-q {1}\' '.format( scanName.format(POI), queue ),
                    )

        executeCommand( cmd )

    os.chdir( currentdir )


def ConvertTChainToArray(
    treeName,
    rootFileList,
    variablePattern = '*',
    ):

    if len(rootFileList) == 0:
        ThrowError( 'rootFileList has length 0' )
        sys.exit()

    # Open first file to get list of all variables
    rootFp = ROOT.TFile.Open( rootFileList[0] )
    tree = rootFp.Get( treeName )

    allVarObjArray = tree.GetListOfBranches()
    nAllVars = allVarObjArray.GetEntries()

    useVars = []
    for iVar in xrange(nAllVars):
        varName = allVarObjArray[iVar].GetTitle()
        if not variablePattern == '*':
            if re.search( variablePattern, varName ):
                useVars.append( varName.split('/',1)[0] )
        else:
            useVars.append( varName.split('/',1)[0] )

    rootFp.Close()

    # Now read the entries from the chain
    chain = ROOT.TChain( 'limit' )
    for rootFile in rootFileList:
        chain.Add( rootFile )


    res = {}
    for varName in useVars:
        res[varName] = []

    for event in chain:
        for varName in useVars:
            res[varName].append( getattr( event, varName ) )

    return res


def InterpretPOI( POI ):
    # Assume it starts with 'r_'
    if not POI.startswith('r_'):
        print 'ERROR: POI {0} does not start with \'r_\''.format( POI )
        return
    components = POI.split('_')[1:]

    if len(components) <= 2:
        print 'ERROR: POI {0} has too few components'.format( POI )
        return

    productionMode  = components[0]
    if not productionMode.endswith('H'):
        print 'ERROR: production mode {0} should end with \'H\''
        return

    observableName  = components[1]

    observableRange = components[2:]

    if len(observableRange) == 1:
        rangeStr = observableRange[0]
        if 'GT' in rangeStr or 'GE' in rangeStr :
            Range = [ float( rangeStr.replace('GT','').replace('GE','') ), 'INF' ]
        elif 'LT' in rangeStr or 'LE' in rangeStr :
            Range = [ '-INF', float( rangeStr.replace('LT','').replace('LE','') ) ]
        else:
            Range = [ float( rangeStr ) ]
    elif len(observableRange) == 2:
        Range = [ float(i) for i in observableRange ]
    else:
        print 'ERROR: Could not make sense of observable range \'{0}\''.format( '_'.join(observableRange) )
        return

    return productionMode, observableName, Range



def ThrowError(
    errstr = ''
    ):

    stack = traceback.extract_stack(None, 2)[0]
    linenumber = stack[1]
    funcname = stack[2]

    cwd = abspath( os.getcwd() )
    modulefilename = relpath( stack[0], cwd )

    print 'ERROR in {0}:{1} {2}:\n    {3}'.format( modulefilename, linenumber, funcname, errstr )




########################################
# Execute this if called directly
########################################

if __name__ == '__main__':
    print 'Import functions from this file'
    print 'Available functions:'
    functions = [ name for (name, thing) in locals().items() if callable(thing) ]
    print '\n'.join( functions )