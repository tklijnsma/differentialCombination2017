#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, tempfile, shutil, re, subprocess, sys, traceback
from os.path import *
from glob import glob

from time import strftime
datestr = strftime( '%b%d' )

from Container import Container

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
def SetTempJobDir( newdirname='tmp' ):
    global TEMPJOBDIR
    TEMPJOBDIR = newdirname
def GetTempJobDir():
    global TEMPJOBDIR
    return TEMPJOBDIR

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
    outputWS=None,
    outputDir=None,
    autoMaps=False,
    manualMaps=None,
    smartMaps=None,
    verbose=False,
    ):

    if outputDir is None:
        outputDir = abspath( 'workspaces_{0}'.format(datestr) )
    if not isdir( outputDir ): os.makedirs( outputDir )
    
    if not outputWS:
        outputWS = basename(datacard).replace( '.txt', '.root' )
    outputWS = join( outputDir, outputWS )

    signalprocesses, processes, bins = ListProcesses( datacard )

    cmd = []
    cmd.append( 'text2workspace.py' )
    cmd.append( datacard )
    cmd.append( '-o {0}'.format(outputWS) )
    cmd.append( '-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel' )
    cmd.append( '--PO verbose' )
    cmd.append( '--PO \'higgsMassRange=123,127\'' )


    if not manualMaps and not smartMaps: autoMaps = True

    if smartMaps:

        # Example procPat:
        # '--PO \'map=.*/InsideAcceptance_genPt_200p0_350p0:r_xH_PTH_200_350[1.0,-1.0,4.0]\'',
        # Want to replace this by:
        # '--PO \'map=.*/InsideAcceptance_genPt_([\dpm]+)p0_([\dpm]+)p0:r_xH_PTH_\1_\2[1.0,-1.0,4.0]\'',

        # in smartMap form: ( r'.*/InsideAcceptance_genPt_([\dm]+)p0_([\dm]+)p0', r'r_xH_PTH_\1_\2[1.0,-1.0,4.0]' )

        # Manual maps should override smart maps; gather all the patters that are already in a manualMap
        if manualMaps:
            manualMapPats = []
            for manualMap in manualMaps:
                match = re.search( r'map=(.*):', manualMap )
                if not match: continue
                manualMapPats.append( match.group(1) )

        newMaps = []
        for binprocPat, yieldParPat in smartMaps:

            for proc in signalprocesses:
                for bin in bins:

                    binprocStr = '{0}/{1}'.format( bin, proc )

                    manualMapAvailable = False
                    if manualMaps:
                        for manualMapPat in manualMapPats:
                            if re.match( manualMapPat, binprocStr ):
                                manualMapAvailable = True

                    if manualMapAvailable:
                        continue
                    elif not re.match( binprocPat, binprocStr ):
                        continue

                    if verbose: print 'Pattern \"{0}\" matches with \"{1}\"'.format( binprocPat, binprocStr )

                    yieldPar = re.sub( binprocPat, yieldParPat, binprocStr )
                    newMap = '--PO \'map={0}:{1}\''.format( binprocStr, yieldPar )

                    newMaps.append( newMap )


        for newMap in newMaps:
            cmd.append(newMap)

    if manualMaps:
        for manualMap in manualMaps:
            cmd.append( manualMap )

    # Default option
    if autoMaps:
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


def BasicT2WSwithModel(
    datacard,
    pathToModel,
    modelName=None,
    suffix=None,
    extraOptions=None,
    smartMaps=None,
    manualMaps=None,
    verbose=True,
    ):

    # ======================================
    # Copy model to compiled directory
    # (scram b takes unnecessarily long)

    if not TESTMODE:
        dst = join( os.environ['CMSSW_BASE'], 'bin', os.environ['SCRAM_ARCH'], basename(pathToModel) )
        print 'Copying\n    {0}\n    to\n    {1}'.format( pathToModel, dst )
        shutil.copyfile( pathToModel, dst )


    # ======================================
    # Build command

    outputDir = abspath( 'workspaces_{0}'.format(datestr) )
    if not isdir( outputDir ): os.makedirs( outputDir )
    outputWS = join( outputDir, basename(datacard).replace( '.txt', '_{0}.root'.format( basename(pathToModel).replace('.py','') ) ) )
    
    if suffix != None:
        outputWS = outputWS.replace( '.root', '_{0}.root'.format(suffix) )

    signalprocesses, processes, bins = ListProcesses( datacard )

    moduleName = basename(pathToModel).replace('.py','')
    if modelName == None:
        modelName = moduleName[0].lower() + moduleName[1:]

    cmd = []
    cmd.append( 'text2workspace.py' )
    cmd.append( datacard )
    cmd.append( '-o {0}'.format(outputWS) )
    cmd.append( '-P {0}:{1}'.format( moduleName, modelName ) )

    if not extraOptions:
        pass
    elif isinstance( extraOptions, basestring ):
        cmd.append( extraOptions )
    else:
        cmd.extend( extraOptions )


    # ======================================
    # Possibility to use maps here as well

    if smartMaps:

        # Manual maps should override smart maps; gather all the patters that are already in a manualMap
        if manualMaps:
            manualMapPats = []
            for manualMap in manualMaps:
                match = re.search( r'map=(.*):', manualMap )
                if not match: continue
                manualMapPats.append( match.group(1) )

        newMaps = []
        for binprocPat, yieldParPat in smartMaps:

            for proc in signalprocesses:
                for bin in bins:

                    binprocStr = '{0}/{1}'.format( bin, proc )

                    manualMapAvailable = False
                    if manualMaps:
                        for manualMapPat in manualMapPats:
                            if re.match( manualMapPat, binprocStr ):
                                manualMapAvailable = True

                    if manualMapAvailable:
                        continue
                    elif not re.match( binprocPat, binprocStr ):
                        continue

                    if verbose: print 'Pattern \"{0}\" matches with \"{1}\"'.format( binprocPat, binprocStr )

                    yieldPar = re.sub( binprocPat, yieldParPat, binprocStr )
                    newMap = '--PO \'map={0}:{1}\''.format( binprocStr, yieldPar )

                    newMaps.append( newMap )


        for newMap in newMaps:
            cmd.append(newMap)

    if manualMaps:
        for manualMap in manualMaps:
            cmd.append( manualMap )



    executeCommand( cmd )

    

def BasicGenericCombineCommand(
        cmd,
        batchJobSubDir = None,
        onBatch        = False,
        sendEmail      = True,
        ):
    
    if onBatch:

        if not batchJobSubDir:
            batchJobSubDir = 'job_{0}'.format( strftime( '%H%M%S' ) )
        tempdir = abspath( join( TEMPJOBDIR, batchJobSubDir ) )
        tempdir = AppendNumberToDirNameUntilItDoesNotExistAnymore( tempdir )

        if not TESTMODE:

            if not isdir(tempdir): os.makedirs(tempdir)

            osHandle, shFile =  tempfile.mkstemp(
                prefix = 'basicbestfitjob_',
                suffix = '.sh',
                dir = tempdir
                )
            shFile = abspath( shFile )

            cmsswPath = join( os.environ['CMSSW_BASE'], 'src' )

            with open( shFile, 'w' ) as shFp:
                if sendEmail:
                    shFp.write( '#$ -M tklijnsm@gmail.com \n' )
                    shFp.write( '#$ -m eas \n' )
                shFp.write( '#$ -o {0} \n'.format( tempdir ) )
                shFp.write( '#$ -e {0} \n'.format( tempdir ) )
                shFp.write( 'export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch/ \n' )
                shFp.write( 'source /cvmfs/cms.cern.ch/cmsset_default.sh \n' )
                shFp.write( 'source /swshare/psit3/etc/profile.d/cms_ui_env.sh \n' )
                shFp.write( 'cd {0} \n'.format( cmsswPath ) )
                shFp.write( 'eval `scramv1 runtime -sh` \n')
                shFp.write( 'cd {0} \n'.format( tempdir ) )
                shFp.write( ' '.join(cmd) )

            qsubCmd = 'qsub -q short.q {0}'.format( shFile )
            executeCommand( qsubCmd )

        else:
            print 'TESTMODE + onBatch not implemented'

    else:
        executeCommand( cmd )




def BasicBestfit(
        datacard,
        setPOIs = True,
        onBatch = False,
        batchJobSubDir = None,
        sendEmail = True,
        extraOptions = None,
        ):
    
    datacard = abspath( datacard )

    cmd = [
        'combine',
        datacard,
        '-M MultiDimFit',
        '--saveNLL',
        # '--saveWorkspace',
        '--minimizerStrategy 2',
        '-v 2',
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

            tempdir = abspath(TEMPJOBDIR)
            if not batchJobSubDir == None:
                tempdir = abspath( join( TEMPJOBDIR, batchJobSubDir ) )
            tempdir = AppendNumberToDirNameUntilItDoesNotExistAnymore( tempdir )

            if not isdir(tempdir): os.makedirs(tempdir)

            osHandle, shFile =  tempfile.mkstemp(
                prefix = 'basicbestfitjob_',
                suffix = '.sh',
                dir = tempdir
                )
            shFile = abspath( shFile )

            cmsswPath = join( os.environ['CMSSW_BASE'], 'src' )

            with open( shFile, 'w' ) as shFp:
                if sendEmail:
                    shFp.write( '#$ -M tklijnsm@gmail.com \n' )
                    shFp.write( '#$ -m eas \n' )
                shFp.write( '#$ -o {0} \n'.format( tempdir ) )
                shFp.write( '#$ -e {0} \n'.format( tempdir ) )
                shFp.write( 'export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch/ \n' )
                shFp.write( 'source /cvmfs/cms.cern.ch/cmsset_default.sh \n' )
                shFp.write( 'source /swshare/psit3/etc/profile.d/cms_ui_env.sh \n' )
                shFp.write( 'cd {0} \n'.format( cmsswPath ) )
                shFp.write( 'eval `scramv1 runtime -sh` \n')
                shFp.write( 'cd {0} \n'.format( tempdir ) )
                shFp.write( ' '.join(cmd) )

            qsubCmd = 'qsub -q short.q {0}'.format( shFile )
            executeCommand( qsubCmd )

        else:
            print 'TESTMODE + onBatch not implemented'

    else:

        moveBack = False

        jobdir = GetTempJobDir()
        if basename(jobdir).replace('/','') != 'tmp':
            if not isdir(jobdir): os.makedirs(jobdir)
            backDir = os.getcwd()
            os.chdir(jobdir)
            moveBack = True

        try:
            executeCommand( cmd )
        finally:
            if moveBack: os.chdir( backDir )



def ListOfPDFIndicesToFreeze( postfitFile, freezeAllIndices=False, verbose=False, snapshotName='MultiDimFit' ):

    wsFp = ROOT.TFile.Open( postfitFile )
    ws   = wsFp.Get('w')
    loadedSnapshot = ws.loadSnapshot(snapshotName)

    if not loadedSnapshot:
        ThrowError( 'Could not load {0} snapshot - are you passing the post fit workspace?'.format( snapshotName ), throwException=True )

    varsToFreeze = []
    varsToFloat = []
    
    # Find all pdf indexes 
    allCats = ws.allCats()
    catitr = allCats.createIterator()
    cat = catitr.Next()
    while cat:
        if cat.GetName().startswith("pdfindex"):
            varsToFreeze.append( cat.GetName() )
        cat = catitr.Next()
        
    # Find all background pdfs
    allPdfs = ws.allPdfs()

    if verbose:
        nPdfs = allPdfs.getSize()
        print 'Looping over {0} pdfs in the workspace'.format(nPdfs)

    pdfitr = allPdfs.createIterator()
    pdf = pdfitr.Next()
    while pdf:
        if pdf.GetName().startswith("shapeBkg_bkg"):
            # bgks from hzz are RooHistPdfs, not RooMultiPdfs
            if hasattr( pdf, 'getNumPdfs' ):
                print pdf
                pdf.Print()
                # Loop over all shapes in the envelope
                for ishape in xrange(pdf.getNumPdfs()):

                    shape = pdf.getPdf( ishape )
                    observables = ROOT.RooArgList(shape.getObservables( ws.allVars() ))
                    observables = filter(lambda x: not x.startswith('CMS_hgg_mass'), map(lambda x: observables[x].GetName(), xrange(observables.getSize()) ) )

                    # Freeze all pdf parameters except those from the best fit function
                    if ishape == pdf.getCurrentIndex() and not freezeAllIndices:
                        varsToFloat.extend( observables )
                    else:
                        varsToFreeze.extend( observables )
        pdf = pdfitr.Next()

    # For some reason, the first variable is never frozen; simply append it again at end of list
    varsToFreeze.append( varsToFreeze[0] )

    if verbose:
        print '\n\nFreezing the following variables:'
        for i in varsToFreeze:
            print i

        print '\n\nFloating the following variables:'
        for i in varsToFloat:
            print i

    wsFp.Close()

    return varsToFreeze



def ConvertFloatToStr( number ):
    number = float(number)
    if number.is_integer():
        number = int(number)
    string = str(number).replace('-','m').replace('.','p')
    return string

def ConvertStrToFloat( string ):
    string = str(string)
    number = string.replace('m','-').replace('p','.')
    number = float(number)
    return number


def ListPOIs( datacardRootFile, nofilter=False ):

    datacardFp = ROOT.TFile.Open( datacardRootFile )
    w = datacardFp.Get('w')
    POIlist = ROOT.RooArgList( w.set('POI') )

    parNames = []
    for i in xrange( POIlist.getSize() ):
        parName = POIlist[i].GetName()
        if nofilter:
            parNames.append( parName )
        elif parName.startswith('r_'):
            parNames.append( parName )

    datacardFp.Close()
    return parNames


def ListSet(
        datacardRootFile,
        setName='POI',
        pattern='*',
        ):
    
    datacardFp = ROOT.TFile.Open( datacardRootFile )
    w = datacardFp.Get('w')
    argset = w.set(setName)
    if not argset:
        print 'No set \'{0}\' in {1}'.format( setName, datacardRootFile )
        datacardFp.Close()
        return []

    arglist = ROOT.RooArgList( argset )

    varNames = []
    for i in xrange( arglist.getSize() ):
        varName = arglist[i].GetName()
        if pattern != '*':
            if re.search( pattern, varName ):
                varNames.append( varName )
        else:
            varNames.append( varName )

    datacardFp.Close()
    return varNames


def ReadBinBoundariesFromWS( datacardRootFile, theory = True ):

    if theory == True:
        setName = 'theoryBinBoundaries'
    else:
        setName = 'expBinBoundaries'

    datacardFp = ROOT.TFile.Open( datacardRootFile )
    w = datacardFp.Get('w')
    argSet = w.set(setName)
    if not argSet:
        ThrowError( 'No set \'{0}\' in {1}'.format( setName, datacardRootFile ) )
        datacardFp.Close()
        sys.exit()
    binBoundaryList = ROOT.RooArgList(argSet)

    parValues = []
    for i in xrange( binBoundaryList.getSize() ):
        parValue = binBoundaryList[i].getVal()
        parValues.append( parValue )

    parValues = list(set(parValues))
    parValues.sort()

    datacardFp.Close()
    return parValues


def ReadTheoryBinBoundariesFromWS( datacardRootFile ):
    return ReadBinBoundariesFromWS( datacardRootFile, theory = True )    
def ReadExpBinBoundariesFromWS( datacardRootFile ):
    return ReadBinBoundariesFromWS( datacardRootFile, theory = False )


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

    processes = list(set(processline.split()[1:]))
    processes.sort()

    signalprocesses = [ p for p in processes if p.split('_')[0].endswith('H') and not p == 'nonResH' and not 'OutsideAcceptance' in p ]
    signalprocesses.sort( key = lambda i: InterpretPOI(i)[2][0] )

    return signalprocesses, processes, bins



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
            output = subprocess.check_output(
                cmd,
                shell=True,
                )
            return output





def executeCommandOnBatch( cmd ):
    print 'Not yet implemented'
    return
    # global TEMPJOBDIR
    # if not isdir(TEMPJOBDIR): os.makedirs(TEMPJOBDIR)




def MultiDimCombineTool(
        datacard,
        nPoints       = 100,
        nPointsPerJob = 3,
        queue         = '1nh',
        notOnBatch    = False,
        jobDirectory  = None,
        fastscan      = False,
        asimov        = False,
        jobPriority   = 0,
        extraOptions  = [],
    ):

    datacard = abspath( datacard )

    scanName = basename(datacard).replace('.root','')
    scanName = re.sub( r'\W', '', scanName )
    scanName = 'SCAN_{0}_{1}_{2}'.format( ''.join(ListPOIs(datacard,nofilter=True)), datestr, scanName )

    currentdir = os.getcwd()
    if not TESTMODE:
        if not jobDirectory:
            jobDirectory = TEMPJOBDIR
        if not isdir( jobDirectory ):
            print 'Creating directory {0}'.format( jobDirectory )
            os.makedirs( jobDirectory )
        print 'Moving to directory {0}'.format( jobDirectory )
        os.chdir( jobDirectory )


    cmd = [
        'combineTool.py',
        datacard,
        '-n {0}'.format( scanName ),
        '-M MultiDimFit',
        '--cminDefaultMinimizerType Minuit2',
        '--cminDefaultMinimizerAlgo migrad',
        '--algo=grid',
        '--floatOtherPOIs=1',
        # '-P "{0}"'.format( POI ),
        # '--setPhysicsModelParameterRanges "{0}"={1:.3f},{2:.3f} '.format( POI, POIRange[0], POIRange[1] ),
        # '--setPhysicsModelParameters {0}'.format( ','.join([ iterPOI + '=1.0' for iterPOI in allPOIs ]) ),
        '-m 125.00',
        '--saveNLL',
        '--saveInactivePOI 1',
        '--points={0} '.format(nPoints),
        # '--minimizerStrategy 2',
        ]

    if fastscan:
        cmd.append( '--fastScan' )
    if asimov:
        cmd.append( '-t -1' )

    if not extraOptions:
        pass
    elif isinstance( extraOptions, basestring ):
        cmd.append( extraOptions )
    else:
        cmd.extend( extraOptions )


    if not notOnBatch:
        cmd.append(
            '--split-points {0} '.format(nPointsPerJob)
            )

        if 't3' in os.environ['HOSTNAME']:

            if not queue in [ 'all.q', 'long.q', 'short.q' ]:
                print 'Queue \'{0}\' is not available on PSI'.format(queue)
                return

            if jobPriority != 0:
                cmd.append(
                    '--job-mode psi --task-name {0} --sub-opts=\'-q {1} -p {2}\' '.format( scanName, queue, jobPriority ),
                    )
            else:
                cmd.append(
                    '--job-mode psi --task-name {0} --sub-opts=\'-q {1}\' '.format( scanName, queue ),
                    )

        else:

            if not queue in [ '8nm', '1nh', '8nh', '1nd', '2nd' ]:
                print 'Queue \'{0}\' is not available on lxplus'.format(queue)
                if not IsTestMode():
                    return

            cmd.append(
                '--job-mode lxbatch --task-name {0} --sub-opts=\'-q {1}\' '.format( scanName, queue ),
                )

    executeCommand( cmd )

    os.chdir( currentdir )

    



def BasicCombineTool(
        datacard,
        POIpattern    = '*',
        POIRange      = [ -1.0, 4.0 ],
        nPoints       = 100,
        nPointsPerJob = 3,
        queue         = '1nh',
        notOnBatch    = False,
        jobDirectory  = None,
        asimov        = False,
        fastscan      = False,
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
    if not jobDirectory:
        jobDirectory = TEMPJOBDIR
    jobDirectory = AppendNumberToDirNameUntilItDoesNotExistAnymore( jobDirectory )

    if not TESTMODE:
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

        if asimov:
            cmd.append( '-t -1' )
        if fastscan:
            cmd.append( '--fastScan' )

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
    rootFileList,
    treeName = 'limit',
    variablePattern = '*',
    returnStyle = 'dict',
    verbose = False
    ):

    if len(rootFileList) == 0:
        ThrowError( 'rootFileList has length 0' )
        sys.exit()

    if verbose: 'Looking for at least one filled root file to obtain the variable list from...'
    foundTree = False
    for rootFile in rootFileList:

        rootFp = ROOT.TFile.Open( rootFile )
        allKeys = rootFp.GetListOfKeys()

        if not allKeys.Contains( 'limit' ):
            if verbose: print '    No tree \'{0}\' in \'{1}\''.format( treeName, rootFile )
            rootFp.Close()
        else:
            tree = rootFp.Get( treeName )
            allVarObjArray = tree.GetListOfBranches()
            treeLoaded = True
            if verbose: print '    Found tree in {0}'.format( rootFile )
            break
    else:
        ThrowError( 'Not a single root file had a tree called \'{0}\'; Cannot extract any data', throwException=True )


    # Get the variable list from this one file
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
    chain = ROOT.TChain( treeName )
    for rootFile in rootFileList:
        chain.Add( rootFile )


    # Return an object

    if returnStyle == 'dict':

        res = {}
        for varName in useVars:
            res[varName] = []

        for event in chain:
            for varName in useVars:
                res[varName].append( getattr( event, varName ) )


    elif returnStyle == 'dictPerPoint':

        res = []
        for event in chain:
            entry = {}
            for varName in useVars:
                entry[varName] = getattr( event, varName )
            res.append( entry )


    elif returnStyle == 'container':

        res = Container()

        for varName in useVars:
            setattr( res. varName, [] )

        for event in chain:
            for varName in useVars:
                getattr( res, varName ).append( getattr( event, varName ) )


    elif returnStyle == 'containerPerPoint':

        res = []

        for event in chain:
            entry = Container()
            for varName in useVars:
                setattr( entry, varName, getattr( event, varName ) )
            res.append( entry )


    return res


def InterpretPOI( POI ):

    # Assume it starts with 'r_'
    if not POI.startswith('r_'):
        # print 'WARNING: POI {0} does not start with \'r_\'; Will now try to manually add it'.format( POI )
        POI = 'r_' + POI
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



class AnalysisError(Exception):
    pass

def ThrowError(
    errstr = '',
    throwException = False
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




########################################
# Execute this if called directly
########################################

if __name__ == '__main__':
    print 'Import functions from this file'
    print 'Available functions:'
    functions = [ name for (name, thing) in locals().items() if callable(thing) ]
    print '\n'.join( functions )