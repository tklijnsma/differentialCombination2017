#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, tempfile, shutil, re, subprocess, sys, traceback, itertools
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
def GlobRootFiles(path):
    if path.endswith('/'):
        path += '*.root'
    else:
        path += '/*.root'
    return glob(path)

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
def BasicCombineCards(
    outputFile,
    *inputList
    ):

    cmd = [ 'combineCards.py' ]
    for datacard in inputList:
        cmd.append( datacard )

    cmd.append( '> {0}'.format( outputFile ) )

    executeCommand( cmd )


#____________________________________________________________________
def BasicT2WS(
        datacard,
        extraOptions=None,
        outputWS=None,
        outputDir=None,
        autoMaps=False,
        manualMaps=None,
        smartMaps=None,
        verbose=False,
        suffix='',
        ):

    if outputDir is None:
        outputDir = abspath( 'out/workspaces_{0}'.format(datestr) )
    if not isdir( outputDir ): os.makedirs( outputDir )
    
    if not outputWS:
        outputWS = basename(datacard).replace( '.txt', '.root' )
    outputWS = outputWS.replace( '.root', suffix + '.root' )
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

#____________________________________________________________________
def CopyPhysicsModels():
    """
    Copies models to compiled directory
    (scram b takes unnecessarily long)
    """

    physicsModelsDir = 'physicsModels'
    dst = join( os.environ['CMSSW_BASE'], 'bin', os.environ['SCRAM_ARCH'], basename(physicsModelsDir) )

    if IsTestMode():
        print 'Would now copy {0} to {1}'.format( physicsModelsDir, dst )
    else:
        print 'Copying {0} to {1}'.format( physicsModelsDir, dst )
        if isdir(dst):
            shutil.rmtree(dst)
        shutil.copytree( physicsModelsDir, dst )

#____________________________________________________________________
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

    # if not TESTMODE:
    #     dst = join( os.environ['CMSSW_BASE'], 'bin', os.environ['SCRAM_ARCH'], basename(pathToModel) )
    #     print 'Copying\n    {0}\n    to\n    {1}'.format( pathToModel, dst )
    #     shutil.copyfile( pathToModel, dst )
    CopyPhysicsModels()


    # ======================================
    # Build command

    outputDir = abspath( 'out/workspaces_{0}'.format(datestr) )
    if not isdir( outputDir ): os.makedirs( outputDir )
    outputWS = join( outputDir, basename(datacard).replace( '.txt', '_{0}.root'.format( basename(pathToModel).replace('.py','') ) ) )
    
    if suffix != None:
        outputWS = outputWS.replace( '.root', '_{0}.root'.format(suffix) )

    moduleName = basename(pathToModel).replace('.py','')
    if modelName == None:
        modelName = moduleName[0].lower() + moduleName[1:]

    if basename(dirname(pathToModel)) == 'physicsModels':
        moduleName = 'physicsModels.' + moduleName

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

        signalprocesses, processes, bins = ListProcesses( datacard )

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

    
#____________________________________________________________________
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



#____________________________________________________________________
def BasicBestfit(
        datacard,
        setPOIs = True,
        onBatch = False,
        sendEmail = True,
        extraOptions = None,
        directory = None,
        ):
    
    datacard = abspath( datacard )

    cmd = [
        'combine',
        datacard,
        '-M MultiDimFit',
        '--saveNLL',
        '--minimizerStrategy 2',
        '-v 2',
        ]

    if setPOIs:
        if IsTestMode():
            parNames = [ 'POIs', 'in', 'the', 'POI', 'set' ]
        else:
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

        if directory is None:
            directory = abspath( GetTempJobDir() )            
        # directory = AppendNumberToDirNameUntilItDoesNotExistAnymore(directory)

        cmsswPath = join( os.environ['CMSSW_BASE'], 'src' )
        shText = []
        if sendEmail:
            shText.extend([
            '#$ -M tklijnsm@gmail.com',
            '#$ -m eas',
            ])
        shText.extend([
            '#$ -o {0}'.format( directory ),
            '#$ -e {0}'.format( directory ),
            'export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch/',
            'source /cvmfs/cms.cern.ch/cmsset_default.sh',
            'source /swshare/psit3/etc/profile.d/cms_ui_env.sh',
            'cd {0}'.format( cmsswPath ),
            'eval `scramv1 runtime -sh`',
            'cd {0}'.format( directory ),
            ' '.join(cmd),
            ])

        print ''
        if not TESTMODE:
            print '[info] Creating', directory
            os.makedirs(directory)

            osHandle, shFile =  tempfile.mkstemp(
                prefix = 'basicbestfitjob_',
                suffix = '.sh',
                dir = directory
                )
            shFile = abspath( shFile )

            with open( shFile, 'w' ) as shFp:
                shFp.write( '\n'.join(shText) )

            qsubCmd = 'qsub -q short.q {0}'.format( shFile )
            executeCommand( qsubCmd )

        else:
            print '[TESTMODE] Would now create', relpath( directory, os.getcwd() )
            print '[TESTMODE] Would now create {0}/some_script.sh, with contents:'.format( relpath( directory, os.getcwd() ) )
            print '\n' + '\n'.join(shText)
            print '[TESTMODE] Would now execute \'qsub -q short.q {0}/some_script.sh\''.format( relpath( directory, os.getcwd() ) )

    else:
        with EnterDirectory( directory ):
            executeCommand( cmd )


#____________________________________________________________________
def ComputeCorrMatrix(
        ws,
        redoPostfit = True,
        postfitFilename = None,
        onBatch = True,
        asimov = False,
        ):

    wsTag = basename(ws).replace('/','').replace('.root','')
    directory = 'corrMat_{0}_{1}'.format( datestr, wsTag )
    corrmatFilename = join( directory, 'higgsCombine_CORRMAT_{0}.MultiDimFit.mH125.root'.format( wsTag ) )

    if postfitFilename is None:
        postfitFilename = join( directory, 'higgsCombine_POSTFIT_{0}.MultiDimFit.mH125.root'.format( wsTag ) )
    else:
        if not isfile(postfitFilename):
            ThrowError( 'Given ws {0} does not exist'.format(postfitFilename) )

    if redoPostfit:
        directory = AppendNumberToDirNameUntilItDoesNotExistAnymore(directory)
        # First regular best fit
        BasicBestfit(
            ws,
            onBatch = onBatch,
            directory = directory,
            extraOptions = [
                '-m 125',
                '--floatOtherPOIs=1',
                # '--computeCovarianceMatrix=1',
                '--saveWorkspace',
                '-n _POSTFIT_{0}'.format( wsTag ),
                ( '-t -1' if asimov else '' )
                ]
            )

    if IsTestMode():
        pdfIndicesToFreeze = [ 'some', 'pdfs' ]
    else:
        pdfIndicesToFreeze = ListOfPDFIndicesToFreeze( postfitFilename, verbose=False )

    BasicBestfit(
        postfitFilename,
        onBatch = onBatch,
        directory = directory,
        extraOptions = [
            '-m 125',
            '--floatOtherPOIs=1',
            '--algo none',
            '--snapshotName MultiDimFit',
            '--saveWorkspace',
            '--computeCovarianceMatrix=1',
            '--freezeNuisances {0}'.format( ','.join(pdfIndicesToFreeze) ),
            '-n _CORRMAT_{0}'.format( wsTag ),
            ( '-t -1' if asimov else '' )
            ]
        )

#____________________________________________________________________
def ListOfPDFIndicesToFreeze( postfitFile, freezeAllIndices=False, verbose=False, snapshotName='MultiDimFit', returnAlsoFloating=False ):

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

    if returnAlsoFloating:
        return varsToFreeze, varsToFloat
    else:
        return varsToFreeze


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
def ListSet(
        datacardRootFile,
        setName='POI',
        pattern='*',
        ):

    closeFp = False
    if isinstance( datacardRootFile, basestring ):
        datacardFp = ROOT.TFile.Open( datacardRootFile )
        w = datacardFp.Get('w')
        closeFp = True
    elif isinstance( datacardRootFile, ROOT.RooWorkspace ):
        w = datacardRootFile

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

    if closeFp: datacardFp.Close()
    return varNames

#____________________________________________________________________
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

#____________________________________________________________________
def ReadTheoryBinBoundariesFromWS( datacardRootFile ):
    return ReadBinBoundariesFromWS( datacardRootFile, theory = True )    
def ReadExpBinBoundariesFromWS( datacardRootFile ):
    return ReadBinBoundariesFromWS( datacardRootFile, theory = False )

#____________________________________________________________________
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
def executeCommandOnBatch( cmd ):
    print 'Not yet implemented'
    return
    # global TEMPJOBDIR
    # if not isdir(TEMPJOBDIR): os.makedirs(TEMPJOBDIR)

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

    if not jobDirectory:
        jobDirectory = TEMPJOBDIR
    if not TESTMODE:
        if not isdir( jobDirectory ):
            print 'Creating directory {0}'.format( jobDirectory )
            os.makedirs( jobDirectory )
        print 'Moving to directory {0}'.format( jobDirectory )
        os.chdir( jobDirectory )
    else:
        print '\nTESTMODE: Would now create {0}'.format(jobDirectory)


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


#____________________________________________________________________    
def BasicCombineTool(
        datacard,
        POIpattern    = '*',
        POIRange      = None,
        nPoints       = 100,
        nPointsPerJob = 3,
        queue         = '1nh',
        notOnBatch    = False,
        jobDirectory  = None,
        asimov        = False,
        fastscan      = False,
        extraOptions  = None,
        disableFloatOtherPOIs = False,
        physicsModelParameterRanges = [],
        setPhysicsModelParameters = True
        ):

    datacard = abspath( datacard )

    scanName = basename(datacard).replace('.root','')
    scanName = re.sub( r'\W', '', scanName )
    scanName = 'SCAN_{{0}}_{0}_{1}'.format( datestr, scanName )

    if POIRange is None:
        POIRange = [ -1.0, 4.0 ]

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
    if asimov: jobDirectory += '_asimov'
    jobDirectory = AppendNumberToDirNameUntilItDoesNotExistAnymore( jobDirectory )

    if not TESTMODE:
        if not isdir( jobDirectory ):
            print 'Creating directory {0}'.format( jobDirectory )
            os.makedirs( jobDirectory )
        print 'Moving to directory {0}'.format( jobDirectory )
        os.chdir( jobDirectory )
    else:
        print '\nTESTMODE: Would now create {0}'.format(jobDirectory)


    for POI in doPOIs:

        cmd = [
            'combineTool.py',
            datacard,
            '-n {0}'.format( scanName.format(POI) ),
            '-M MultiDimFit',
            '--cminDefaultMinimizerType Minuit2',
            '--cminDefaultMinimizerAlgo migrad',
            '--algo=grid',
            ( '--floatOtherPOIs=1' if not disableFloatOtherPOIs else '' ),
            '-P "{0}"'.format( POI ),
            # '--setPhysicsModelParameterRanges "{0}"={1:.3f},{2:.3f} '.format( POI, POIRange[0], POIRange[1] ),
            '-m 125.00',
            '--squareDistPoi',
            '--saveNLL',
            '--saveInactivePOI 1',
            '--points={0} '.format(nPoints),
            ]

        if setPhysicsModelParameters:
            cmd.append('--setPhysicsModelParameters {0}'.format( ','.join([ iterPOI + '=1.0' for iterPOI in allPOIs ]) ))

        physicsModelParameterRanges.insert( 0, [ POI, POIRange[0], POIRange[1] ] )
        physicsModelParameterRangesStr = ':'.join([ '"{0}"={1:.3f},{2:.3f}'.format(parName, left, right) for parName, left, right in physicsModelParameterRanges ])
        cmd.append( '--setPhysicsModelParameterRanges ' + physicsModelParameterRangesStr )

        if asimov:
            cmd.append( '-t -1' )
        if fastscan:
            cmd.append( '--fastScan' )

        if not extraOptions is None:
            cmd.extend(extraOptions)

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

#____________________________________________________________________
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
        ThrowError(
            'Not a single root file had a tree called \'{0}\'; Cannot extract any data'.format(treeName)
            + '\n  First root file:\n  {0}'.format( rootFileList[0] )
            + '\n  Last root file:\n  {0}'.format( rootFileList[-1] )
            ,
            throwException=True
            )


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


#____________________________________________________________________
def newColorCycle():
    return itertools.cycle( range(2,5) + range(6,10) + range(40,50) + [ 30, 32, 33, 35, 38, 39 ] )


#____________________________________________________________________
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
            Range = [ ConvertStrToFloat( rangeStr.replace('GT','').replace('GE','') ), 'INF' ]
        elif 'LT' in rangeStr or 'LE' in rangeStr :
            Range = [ '-INF', ConvertStrToFloat( rangeStr.replace('LT','').replace('LE','') ) ]
        else:
            Range = [ ConvertStrToFloat( rangeStr ) ]
    elif len(observableRange) == 2:
        Range = [ ConvertStrToFloat(i) for i in observableRange ]
    else:
        print 'ERROR: Could not make sense of observable range \'{0}\''.format( '_'.join(observableRange) )
        return

    return productionMode, observableName, Range



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

    def __init__(self, subDirectory=None ):
        self._active = False
        if not subDirectory is None and not subDirectory == '':
            self._active = True
            self.backDir = os.getcwd()
            self.subDirectory = subDirectory

    def __enter__(self):
        if self._active:
            if IsTestMode():
                print '\n[TESTMODE] Would now create/go into \'{0}\''.format(self.subDirectory)
            else:
                print ''
                if not isdir( self.subDirectory ):
                    print 'Creating \'{0}\''.format( relpath( self.subDirectory, self.backDir ) )
                    os.makedirs( self.subDirectory )
                print 'Entering \'{0}\''.format( relpath( self.subDirectory, self.backDir ) )
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