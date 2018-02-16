def BasicCombineCards(
    outputFile,
    *inputList
    ):

    cmd = [ 'combineCards.py' ]
    for datacard in inputList:
        cmd.append( datacard )

    cmd.append( '> {0}'.format( outputFile ) )

    executeCommand( cmd )


def ProcessMaps(datacard, smartMaps=None, manualMaps=None, verbose=True):

    if smartMaps is None and manualMaps is None:
        # raise ValueError('Supply smartMaps and/or manualMaps')
        return []

    signalprocesses, processes, bins = ListProcesses( datacard )

    maps = []

    if smartMaps:
        # Example procPat:
        # '--PO \'map=.*/InsideAcceptance_genPt_200p0_350p0:r_xH_PTH_200_350[1.0,-1.0,4.0]\'',
        # Want to replace this by:
        # '--PO \'map=.*/InsideAcceptance_genPt_([\dpm]+)p0_([\dpm]+)p0:r_xH_PTH_\1_\2[1.0,-1.0,4.0]\'',
        # in smartMap form: ( r'.*/InsideAcceptance_genPt_([\dm]+)p0_([\dm]+)p0', r'r_xH_PTH_\1_\2[1.0,-1.0,4.0]' )

        # Manual maps should have priority over smart maps;
        # gather all the patterns that are already in a manualMap
        if manualMaps:
            manualMapPats = []
            for manualMap in manualMaps:
                match = re.search( r'map=(.*):', manualMap )
                if not match: continue
                manualMapPats.append( match.group(1) )

        for binprocPat, yieldParPat in smartMaps:
            for proc in signalprocesses:
                for bin in bins:
                    binprocStr = '{0}/{1}'.format( bin, proc )

                    # Check if there are manual patterns that also match this string
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
                    maps.append( '--PO \'map={0}:{1}\''.format( binprocStr, yieldPar ))

    # Add manual maps after (the last matching pattern is applied in MultiSignalModel)
    if manualMaps:
        for manualMap in manualMaps:
            maps.append( manualMap )

    return maps


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

    cmd = []
    cmd.append( 'text2workspace.py' )
    cmd.append( datacard )
    cmd.append( '-o {0}'.format(outputWS) )
    cmd.append( '-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel' )
    cmd.append( '--PO verbose' )
    cmd.append( '--PO \'higgsMassRange=123,127\'' )

    maps = ProcessMaps(datacard, smartMaps, manualMaps, verbose)
    for _map in maps:
        cmd.append( _map )

    if not extraOptions:
        pass
    elif isinstance( extraOptions, basestring ):
        cmd.append( extraOptions )
    else:
        cmd.extend( extraOptions )

    executeCommand( cmd )

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

    # Possibility to use maps here as well
    maps = ProcessMaps(datacard, smartMaps, manualMaps, verbose)
    for _map in maps:
        cmd.append( _map )

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

