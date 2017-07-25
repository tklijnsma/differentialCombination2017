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
    manualMaps=None,
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

    if manualMaps:

        for manualMap in manualMaps:
            cmd.append( manualMap )

    else:
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
    extraOptions=None,
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




def MultiDimCombineTool(
        datacard,
        nPoints       = 100,
        nPointsPerJob = 3,
        queue         = '1nh',
        notOnBatch    = False,
        jobDirectory  = None,
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

            cmd.append(
                '--job-mode psi --task-name {0} --sub-opts=\'-q {1}\' '.format( scanName, queue ),
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
    rootFileList,
    treeName = 'limit',
    variablePattern = '*',
    ):

    if len(rootFileList) == 0:
        ThrowError( 'rootFileList has length 0' )
        sys.exit()

    foundTree = False
    for rootFile in rootFileList:

        rootFp = ROOT.TFile.Open( rootFile )
        allKeys = rootFp.GetListOfKeys()

        if not allKeys.Contains( 'limit' ):
            print 'No tree \'{0}\' in \'{1}\''.format( treeName, rootFile )
            rootFp.Close()
        else:
            tree = rootFp.Get( treeName )
            allVarObjArray = tree.GetListOfBranches()
            treeLoaded = True
            print 'Found tree in {0}'.format( rootFile )
            break


    # tryFile = 0
    # treeLoaded = False
    # while tryFile < len(rootFileList):
    #     try:
    #         # Open first file to get list of all variables
    #         rootFp = ROOT.TFile.Open( rootFileList[tryFile] )
    #         allkeys = rootFp.GetListOfKeys()

    #         print rootFileList[tryFile]

    #         allkeys.Print()

    #         print allkeys.Contains( 'limit' )
    #         sys.exit()

    #         tree = rootFp.Get( treeName )
    #         allVarObjArray = tree.GetListOfBranches()
    #         treeLoaded = True
    #         print 'Found tree in {0}'.format( rootFileList[tryFile] )
    #     except AttributeError:
    #         print 'No tree \'{0}\' in \'{1}\''.format( treeName, rootFileList[tryFile] )
    #         rootFp.Close()
    #     finally:
    #         tryFile += 1

    # if not treeLoaded:
    #     print 'ERROR: Could not load tree \{0}\' for any file'.format( treeName )
    #     return


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
        print 'WARNING: POI {0} does not start with \'r_\'; Will now try to manually add it'.format( POI )
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