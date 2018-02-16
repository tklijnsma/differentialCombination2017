
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
def GlobRootFiles(path):
    if path.endswith('/'):
        path += '*.root'
    else:
        path += '/*.root'
    return glob(path)









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

