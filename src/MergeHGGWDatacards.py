#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, re
from os.path import *
from copy import deepcopy

from time import strftime
datestr = strftime( '%b%d' )

import ROOT
import Commands


class Container:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


########################################
# Main
########################################


TEMPCARDDIR = 'tempcards_{0}'.format( datestr )
def makeTempDir( newTempDir=None ):
    global TEMPCARDDIR
    if newTempDir: TEMPCARDDIR = newTempDir
    if not isdir(TEMPCARDDIR): os.makedirs(TEMPCARDDIR)


def main():

    print 'Don\'t execute this.'

    # xhDatacard  = '../suppliedInput/fromVittorio/differential_pT_moriond17_HxOnly_July19/Datacard_13TeV_differential_pT_moriond17_HxOnly.txt'    
    # gghDatacard = '../suppliedInput/fromVittorio/pT_ggHonly_Jun26/Datacard_13TeV_differential_pT_moriond17_ggHonly_v2.txt'

    # Merge_xH_ggH_hgg_cards(
    #         xhDatacard,
    #         gghDatacard,
    #         )


########################################
# njets combination
########################################

def RenameProcesses_Hgg_differentials(
        indatacard,
        globalReplace = None,
        writeToFile = True,
        outTag = '_renamedProcesses'
        ):

    with open( indatacard, 'r' ) as indatacardFp:
        indatacardTxt = indatacardFp.read()


    # ======================================
    # Renaming

    outdatacardTxt = indatacardTxt

    # Replace multiple times to cover neighboring matches

    # InsideAcceptance_myGenNjets2p5_3p5to100p0
    # InsideAcceptance_myGenNjets2p5_2p5to3p5

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)InsideAcceptance_myGenNjets2p5_3p5to100p0(\W)',
            r'\1smH_NJ_GE4\2',
            outdatacardTxt
            )

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)InsideAcceptance_myGenNjets2p5_([m\d]+)p5to(\d+)p5(\W)',
            r'\1smH_NJ_\3\4',
            outdatacardTxt
            )



    # Slightly updated conventions for the NNLOPS datacards (from Nov 10)

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)InsideAcceptance_genNjets2p5_3p5_100p0(\W)',
            r'\1smH_NJ_GE4\2',
            outdatacardTxt
            )

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)InsideAcceptance_genNjets2p5_([m\d]+)p5_(\d+)p5(\W)',
            r'\1smH_NJ_\3\4',
            outdatacardTxt
            )


    # Rapidity renaming Nov12

    # Has to be done mosty manually... David used "0p30", Vittorio used "0p3"

    # smH_YH_0p0_0p15 
    # smH_YH_0p15_0p30
    # smH_YH_0p30_0p60
    # smH_YH_0p60_0p90
    # smH_YH_0p90_1p20
    # smH_YH_1p20_2p50

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)InsideAcceptance_genAbsRapidity_0p0_0p15(\W)',
            r'\1smH_YH_0p0_0p15\2',
            outdatacardTxt
            )

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)InsideAcceptance_genAbsRapidity_0p15_0p3(\W)',
            r'\1smH_YH_0p15_0p30\2',
            outdatacardTxt
            )

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)InsideAcceptance_genAbsRapidity_0p3_0p6(\W)',
            r'\1smH_YH_0p30_0p60\2',
            outdatacardTxt
            )

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)InsideAcceptance_genAbsRapidity_0p6_0p9(\W)',
            r'\1smH_YH_0p60_0p90\2',
            outdatacardTxt
            )

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)InsideAcceptance_genAbsRapidity_0p9_3p0(\W)',
            r'\1smH_YH_0p90_2p50\2',
            outdatacardTxt
            )


    # ======================================
    # Process simple global replacements

    if not globalReplace == None:
        for string, replacement in globalReplace:
            outdatacardTxt = outdatacardTxt.replace( string, replacement )


    # ======================================
    # Write to file / return

    if writeToFile:

        outdatacard = indatacard.replace( '.txt', '{0}.txt'.format(outTag) )

        if Commands.IsTestMode():
            print '\nContents of the new datacard:\n\n'
            print outdatacardTxt
            print '\n\nWould now open \'{0}\' to dump these contents'.format( outdatacard )
        else:
            with open( outdatacard, 'w' ) as outdatacardFp:
                outdatacardFp.write( outdatacardTxt )
            print '[info] Wrote output to ' + outdatacard

    return outdatacardTxt





########################################
# pt combination
########################################


def RenameProcesses_Aug21(
    indatacard,
    renameOutsideAcceptance = True,
    globalReplace = None,
    writeToFile = True,
    outTag = '_renamedProcesses'
    ):

    with open( indatacard, 'r' ) as indatacardFp:
        indatacardTxt = indatacardFp.read()


    # ======================================
    # Renaming

    outdatacardTxt = indatacardTxt

    # Replace multiple times to cover neighboring matches


    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)gghInsideAcceptance_genPt_350p0_10000p0(\W)',
            r'\1ggH_PTH_GT350\2',
            outdatacardTxt
            )

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)gghInsideAcceptance_genPt_(\d+)p0_(\d+)p0(\W)',
            r'\1ggH_PTH_\2_\3\4',
            outdatacardTxt
            )

    if renameOutsideAcceptance:
        outdatacardTxt = re.sub(
            r'(\W)gghOutsideAcceptance(\W)',
            r'\1ggH_OutsideAcceptance\2',
            outdatacardTxt
            )


    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)hxInsideAcceptance_genPt_350p0_10000p0(\W)',
            r'\1xH_PTH_GT350\2',
            outdatacardTxt
            )

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)hxInsideAcceptance_genPt_(\d+)p0_(\d+)p0(\W)',
            r'\1xH_PTH_\2_\3\4',
            outdatacardTxt
            )

    if renameOutsideAcceptance:
        outdatacardTxt = re.sub(
            r'(\W)hxOutsideAcceptance(\W)',
            r'\1xH_OutsideAcceptance\2',
            outdatacardTxt
            )

    # ----------------
    # For SMH card

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)InsideAcceptance_genPt_350p0_10000p0(\W)',
            r'\1smH_PTH_GT350\2',
            outdatacardTxt
            )

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)InsideAcceptance_genPt_(\d+)p0_(\d+)p0(\W)',
            r'\1smH_PTH_\2_\3\4',
            outdatacardTxt
            )

    if renameOutsideAcceptance:
        outdatacardTxt = re.sub(
            r'(\W)OutsideAcceptance(\W)',
            r'\1OutsideAcceptance\2',
            outdatacardTxt
            )


    # ======================================
    # Process simple global replacements

    if not globalReplace == None:
        for string, replacement in globalReplace:
            outdatacardTxt = outdatacardTxt.replace( string, replacement )


    # ======================================
    # Write to file / return

    if writeToFile:
        outdatacard = indatacard.replace( '.txt', '{0}.txt'.format(outTag) )

        if Commands.IsTestMode():
            print '\nContents of the new datacard:\n\n'
            print outdatacardTxt
            print '\n\nWould now open \'{0}\' to dump these contents'.format( outdatacard )
        else:
            with open( outdatacard, 'w' ) as outdatacardFp:
                outdatacardFp.write( outdatacardTxt )
            print '[info] Wrote output to ' + outdatacard

    return outdatacardTxt



#____________________________________________________________________
def sorter( process ):

    if not isinstance( process, basestring ):
        process = process[0]

    ret = 100000

    if process.startswith( 'ggH' ):
        ret -= 10000
    elif process.startswith( 'xH' ):
        ret -= 20000
    elif process.startswith( 'smH' ):
        ret -= 30000

    if process.startswith( 'ggH' ) or process.startswith( 'xH' ) or process.startswith( 'smH' ):
        components = process.split('_')
        if len(components) >= 3:
            leftBound = int( components[2].replace('GT','') )
            ret += leftBound

    else:
        ret += 1000

    return ret


def RenumberProcessesHZZ_Aug21(
        datacardFile,
        ):

    with open( datacardFile, 'r' ) as datacardFp:
        lines = datacardFp.readlines()

    firstMatch = True
    for iLine, line in enumerate(lines):

        if line.startswith( 'process ' ):
            if firstMatch:
                names = line.split()[1:]
                firstMatch = False
                iNameLine = iLine
            else:
                numbers = line.split()[1:]
                iNumberLine = iLine
                break

    else:
        print 'Reached end of for-loop; Something is wrong'

    signals = []
    bkgs    = []
    for process in list(set(names)):
        if process.startswith('ggH') or process.startswith('xH') or process.startswith('smH'):
            signals.append(process)
        else:
            bkgs.append(process)

    signals.sort( key=sorter )
    bkgs.sort( key=sorter )

    numberDict = {}
    for iSignal, signal in enumerate( signals ):
        numberDict[signal] = -iSignal
    for iBkg, bkg in enumerate( bkgs ):
        numberDict[bkg] = iBkg + 1


    newNumberLine = '{0:39}'.format('process')
    for name in names:
        newNumberLine += '{0:<25}'.format( numberDict[name] )
    lines[iNumberLine] = newNumberLine + '\n'


    out = datacardFile.replace( '.txt', '_processesShifted.txt' )
    with open( out, 'w' ) as outFp:
        outFp.write( ''.join(lines) )
    print '[info] Wrote output to \'{0}\''.format( out )


#____________________________________________________________________
# There is no need for this function, channel renaming is only for convenience
def RenameProcessesHZZchannels(
    process,
    indatacard,
    outdatacard=None,
    globalReplace=None,
    ):

    with open( indatacard, 'r' ) as indatacardFp:
        indatacardTxt = indatacardFp.read()


    # ======================================
    # Renaming

    outdatacardTxt = indatacardTxt

    # Replace multiple times to cover neighboring matches

    # HIER VERDER

    # ch1 --> smH_... enzo doen




    # for i in xrange(3):
    #     outdatacardTxt = re.sub(
    #         r'(\W)InsideAcceptance_genPt_350p0_10000p0(\W)',
    #         r'\1{0}_PTH_GT350\2'.format(process),
    #         outdatacardTxt
    #         )

    # for i in xrange(3):
    #     outdatacardTxt = re.sub(
    #         r'(\W)InsideAcceptance_genPt_(\d+)p0_(\d+)p0(\W)',
    #         r'\1{0}_PTH_\2_\3\4'.format(process),
    #         outdatacardTxt
    #         )


    # ======================================
    # Process simple global replacements

    if not globalReplace == None:
        for string, replacement in globalReplace:
            outdatacardTxt = outdatacardTxt.replace( string, replacement )


    # ======================================
    # Write to file / return

    if outdatacard == 'auto':
        outdatacard = indatacard.replace( '.txt', '_renamedProcesses_{0}.txt'.format(datestr) )
        with open( outdatacard, 'w' ) as outdatacardFp:
            outdatacardFp.write( outdatacardTxt )
        print '[info] Wrote output to ' + outdatacard
    elif not outdatacard == None:
        outdatacard = join( TEMPCARDDIR, basename(indatacard).replace( '.txt', '_renamedProcesses.txt' ) )
        with open( outdatacard, 'w' ) as outdatacardFp:
            outdatacardFp.write( outdatacardTxt )
        print '[info] Wrote output to ' + outdatacard

    return outdatacardTxt



def RenameProductionModeHgg( process, wsFile ):

    outRootFile = wsFile.replace( '.root', '_FullyRenamed_{0}.root'.format(datestr) )
    fout = ROOT.TFile( outRootFile, 'recreate' )

    # ======================================
    # Get input ws

    wsFp   = ROOT.TFile.Open( wsFile )
    w      = ROOT.gDirectory.Get('wsig_13TeV')

    componentArgset   = w.components()
    nComponents       = componentArgset.getSize()

    componentIterator = w.componentIterator()

    for i in xrange(nComponents):
        element = componentIterator.Next()
        name = element.GetName()

        if 'InsideAcceptance' in name or 'OutsideAcceptance' in name:
            newname = name.replace( 'InsideAcceptance', '{0}_InsideAcceptance'.format(process) ).replace( 'OutsideAcceptance', '{0}_OutsideAcceptance'.format(process) )
            print 'Replacing \'{0}\' with \'{1}\''.format( name, newname )
            element.SetName( newname )
        else:
            print 'Not renaming \'{0}\''.format(name)

        # newname = process + '_' + name
        # print 'Replacing \'{0}\' with \'{1}\''.format( name, newname )
        # element.SetName( newname )



    # ======================================
    # Write to file

    fout.WriteTObject(w)
    fout.Close()

    wsFp.Close()





def RenameProductionModeHgg_PerMethod( process, wsFile ):

    # ======================================
    # Arrange IO for output

    if Commands.IsTestMode():
        makeTempDir( 'RenamedTests_{0}'.format(datestr) )
        outRootFile = join( TEMPCARDDIR, basename(wsFile.replace('.root','_feaRenamed.root')) )
    else:
        outRootFile = wsFile.replace('.root','_FullyRenamed.root')

    fout = ROOT.TFile( outRootFile, 'recreate' )


    # ======================================
    # Get input ws

    wsFp   = ROOT.TFile.Open( wsFile )
    w      = ROOT.gDirectory.Get('wsig_13TeV')



    allMethods = [
        'allVars',
        'allCats',
        'allFunctions',
        'allCatFunctions',
        'allPdfs',
        'allResolutionModels',
        'allData',
        'allEmbeddedData',
        'allGenericObjects',
        ]

    for method in allMethods:

        argset   = getattr( w, method )()

        elements = []
        try:
            iterator = argset.createIterator()
            for i in xrange(argset.getSize()):
                element = iterator.Next()
                elements.append(element)
        except AttributeError:
            element = argset.begin()
            for i in xrange(len(argset)):
                elements.append(element)
                try:
                    next(element)
                except StopIteration:
                    break


        for element in elements:

            name = element.GetName()

            if 'const' in name: print 'Next replaced process should be \'{0}\' to \'{1}\''.format(
                name,
                name.replace( 'InsideAcceptance', '{0}_InsideAcceptance'.format(process) ),
                )

            if 'InsideAcceptance' in name:
                newname = name.replace( 'InsideAcceptance', '{0}_InsideAcceptance'.format(process) )
                print 'Replacing \'{0}\' with \'{1}\''.format( name, newname )
                element.SetName( newname )

            if 'OutsideAcceptance' in name:
                newname = name.replace( 'OutsideAcceptance', '{0}_OutsideAcceptance'.format(process) )
                print 'Replacing \'{0}\' with \'{1}\''.format( name, newname )
                element.SetName( newname )


    # ======================================
    # Write to file

    fout.WriteTObject(w)
    fout.Close()

    wsFp.Close()





def Rename_fea( process, wsFile ):

    # ======================================
    # Arrange IO for output

    if Commands.IsTestMode():
        makeTempDir( 'Renamed_fea_tests_{0}'.format(datestr) )
        outRootFile = join( TEMPCARDDIR, basename(wsFile.replace('.root','_feaRenamed.root')) )
    else:
        outRootFile = wsFile.replace('.root','_feaRenamed.root')

    fout = ROOT.TFile( outRootFile, 'recreate' )


    # ======================================
    # Get input ws

    wsFp   = ROOT.TFile.Open( wsFile )
    w      = ROOT.gDirectory.Get('wsig_13TeV')

    allFunctionsArgset   = w.allFunctions()
    allFunctionsIterator = allFunctionsArgset.createIterator()

    for i in xrange(allFunctionsArgset.getSize()):
        element = allFunctionsIterator.Next()

        name = element.GetName()
        if 'fea_' in name:
            newname = name.replace( 'fea_', 'fea_{0}_'.format(process) )
            print 'Renaming \'{0}\' to \'{1}\''.format( name, newname )
            element.SetName( newname )



    # ======================================
    # Write to file

    fout.WriteTObject(w)
    fout.Close()

    wsFp.Close()





def Merge_xH_ggH_hgg_cards(
        xhDatacard,
        gghDatacard,
        ):
    
    # Basic rename
    xhDCrenamed  = RenameProcesses( 'xH',  xhDatacard,  renameOutsideAcceptance=True,
        globalReplace = [(
            'CMS-HGG_sigfit_differential_pT_moriond17_HxOnly.root',
            'CMS-HGG_sigfit_differential_pT_moriond17_HxOnly_feaRenamed.root'
            )]
        )
    gghDCrenamed = RenameProcesses( 'ggH', gghDatacard, renameOutsideAcceptance=True,
        globalReplace = [(
            'CMS-HGG_sigfit_differential_pT_moriond17_ggHonly_v2.root',
            'CMS-HGG_sigfit_differential_pT_moriond17_ggHonly_v2_feaRenamed.root'
            )]
        )


    # Read into container
    xh  = GetDatacardContainer( xhDCrenamed )
    ggh = GetDatacardContainer( gghDCrenamed )
    xh.process = 'xH'
    ggh.process = 'ggH'

    # Include one directory upwards in the path of the root files
    ExtendPathOfRootFiles( xh, dirname(xhDatacard) )
    ExtendPathOfRootFiles( ggh, dirname(gghDatacard) )

    # Merge
    smh = MergeCards( ggh, xh )

    outpath = join( dirname(gghDatacard), '..', 'hggMerged_{0}.txt'.format(datestr) )

    ParseDataContainer( smh, writeToFile = outpath )



def MergeCards( c1_original, c2_original ):
    c1 = deepcopy(c1_original)
    c2 = deepcopy(c2_original)

    merged = deepcopy( c1 )

    # ======================================
    # Add shapes lines (only add lines that start with the process name, the rest are duplicates)

    for shapesLine in c2.shapesLines:
        if shapesLine.startswith( 'shapes {0}'.format(c2.process) ):
            merged.shapesLines.append( shapesLine )


    # ======================================
    # Renumber the processes in the second datacard

    mostNegativeProcessNumber = abs(min([ int(i) for i in c1.processNumbers ]))
    subtractFromStr = lambda inNumber, subNum=mostNegativeProcessNumber: str( int(inNumber) - subNum )

    for iProcess, process in enumerate(c2.processNames):
        if process.startswith( c2.process ):
            c2.processNumbers[iProcess] = subtractFromStr(c2.processNumbers[iProcess])


    # ======================================
    # Add bin/process entries (for number, rate and nuis pars) to the merged datacard

    entries = []
    # nEntries = len(c2.processNames) * len(c2.binNames)

    splitBinLine            = c2.binLine.split()[1:]
    splitProcessNamesLine   = c2.processNamesLine.split()[1:]
    splitProcessNumbersLine = [ subtractFromStr(i) for i in c2.processNumbersLine.split()[1:] ]
    splitRateLine           = c2.rateLine.split()[1:]

    nEntries = len(splitProcessNamesLine)

    if not len(splitBinLine) == nEntries:
        print '[error] List \'splitBinLine\' contains {0} entries; expected {1}'.format( len(splitBinLine), nEntries )
    if not len(splitProcessNamesLine) == nEntries:
        print '[error] List \'splitProcessNamesLine\' contains {0} entries; expected {1}'.format( len(splitProcessNamesLine), nEntries )
    if not len(splitProcessNumbersLine) == nEntries:
        print '[error] List \'splitProcessNumbersLine\' contains {0} entries; expected {1}'.format( len(splitProcessNumbersLine), nEntries )
    if not len(splitRateLine) == nEntries:
        print '[error] List \'splitRateLine\' contains {0} entries; expected {1}'.format( len(splitRateLine), nEntries )


    for iEntry in xrange(nEntries):
        entry = Container()
        entry.bin           = splitBinLine[iEntry]
        entry.processName   = splitProcessNamesLine[iEntry]
        entry.processNumber = splitProcessNumbersLine[iEntry]
        entry.rate          = splitRateLine[iEntry]
        entry.nuisDict      = { nuis.name: nuis.vals[iEntry] for nuis in c2.nuisanceContainers }

        entries.append(entry)


    # ======================================
    # Define selection
    # (Now very simple but could be expanded)

    def acceptEntry( entry ):
        if entry.processName.startswith( c2.process ):
            return True
        else:
            return False


    # ======================================
    # Add them to the merged datacard

    mergedSplitBinLine            = merged.binLine.split()[1:]
    mergedSplitProcessNamesLine   = merged.processNamesLine.split()[1:]
    mergedSplitProcessNumbersLine = merged.processNumbersLine.split()[1:]
    mergedSplitRateLine           = merged.rateLine.split()[1:]
    mergedSplitNuisanceLines      = [ nuisLine.split()[2:] for nuisLine in merged.nuisanceLines ]

    for entry in entries:

        if not acceptEntry(entry):
            continue

        mergedSplitBinLine.append(            entry.bin )
        mergedSplitProcessNamesLine.append(   entry.processName )
        mergedSplitProcessNumbersLine.append( entry.processNumber )
        mergedSplitRateLine.append(           entry.rate )

        for nuis in merged.nuisanceContainers:
            nuis.vals.append(  entry.nuisDict[nuis.name]  )

    merged.splitBinLine            = mergedSplitBinLine
    merged.splitProcessNamesLine   = mergedSplitProcessNamesLine
    merged.splitProcessNumbersLine = mergedSplitProcessNumbersLine
    merged.splitRateLine           = mergedSplitRateLine

    # Also include the parsed ones
    merged.binLine            = 'bin '     + ' '.join(mergedSplitBinLine)
    merged.processNamesLine   = 'process ' + ' '.join(mergedSplitProcessNamesLine)
    merged.processNumbersLine = 'process ' + ' '.join(mergedSplitProcessNumbersLine)
    merged.rateLine           = 'rate '    + ' '.join(mergedSplitRateLine)

    return merged




def ExtendPathOfRootFiles( card, basepath ):
    card.unextendedShapesLines = card.shapesLines
    card.shapesLines = []
    for shapesLine in card.unextendedShapesLines:
        components = [ i.strip() for i in shapesLine.split(' ') if len(i.strip())>0 ]

        if not len(components) == 5:
            print '[error] The following shapesLine makes no sense:'
            print shapesLine
            return

        rootPath = components[3]
        newRootPath = relpath( join( basepath, rootPath ), join( basepath, '..' ) )

        components[3] = newRootPath

        card.shapesLines.append(
            ' '.join(components)
            )


def GetDatacardContainer(
        datacardTxt
        ):

    
    if isfile( datacardTxt ):
        print 'Given argument \'{0}\' was found to be a valid path; Trying to read datacard from it'.format(datacardTxt)
        datacardFile = datacardTxt
        with open( datacardFile, 'r' ) as datacardFp:
            datacardTxt = datacardFp.read()


    ret = Container()
    ret.shapesLines = []
    ret.nuisanceLines = []
    ret.discreteNuisanceLines = []
    ret.paramNuisanceLines = []


    lines = [ i.strip() for i in datacardTxt.split('\n') if len(i.strip())>0 and not i.strip().startswith('#') ]

    binLineCounter = 0
    processLineCounter = 0
    for line in lines:

        if line.startswith('imax'):
            ret.imaxLine = line
            continue

        elif line.startswith('jmax'):
            ret.jmaxLine = line
            continue

        elif line.startswith('kmax'):
            ret.kmaxLine = line
            continue

        elif line.startswith('shapes'):
            ret.shapesLines.append(line)
            continue

        elif line.startswith('bin ') and binLineCounter == 0:
            ret.firstBinLine = line
            ret.binNames = line.split()[1:]
            binLineCounter += 1
            continue

        elif line.startswith('bin ') and binLineCounter == 1:
            ret.binLine = line
            binLineCounter += 1
            continue

        elif line.startswith('observation'):
            ret.observationLine = line
            continue

        elif line.startswith('process') and processLineCounter == 0:
            ret.processNamesLine = line
            processes = line.split()[1:]
            processLineCounter += 1
            continue

        elif line.startswith('process') and processLineCounter == 1:
            ret.processNumbersLine = line
            processNumbers = line.split()[1:]
            processLineCounter += 1
            continue

        elif line.startswith('rate'):
            ret.rateLine = line
            continue

        elif line.startswith('pdfindex_'):
            ret.discreteNuisanceLines.append( line )
            continue

        else:

            # Line has to be splittable
            if not ' ' in line: continue

            # Check for normal nuisances
            components = line.split()

            # Not a whole lot of nuisance parameter types implemented
            if components[1] in [ 'lnN', 'gmN', 'lnU' ]:
                ret.nuisanceLines.append( line )
            elif components[1] == 'param':
                ret.paramNuisanceLines.append( line )

    
    # Create sub-containers per nuisance parameter
    ret.nuisanceContainers = []
    for nuisLine in ret.nuisanceLines:
        components = nuisLine.split()

        nuisContainer = Container()
        nuisContainer.name = components[0]
        nuisContainer.type = components[1]
        nuisContainer.vals = components[2:]
        ret.nuisanceContainers.append( nuisContainer )


    # Process names is still a merged mess
    processNumberPairs = [ ( name, number ) for name, number in zip( processes, processNumbers ) ]
    processNumberPairs = list(set(processNumberPairs))
    processNumberPairs.sort( key=lambda pair: pair[1] )

    # Convenient formats
    # ret.processNumberPairs = processNumberPairs
    # ret.processNumberDict  = { name:number for name, number in processNumberPairs }
    ret.processNames       = [ pair[0] for pair in processNumberPairs ]
    ret.processNumbers     = [ pair[1] for pair in processNumberPairs ]
    # ret.processNumberDict  = { name:number for name, number in zip( ret.processNames, processNumbers ) }

    return ret


def ParseDataContainer(
        c,
        writeToFile = None,
        ):


    # ======================================
    # Function that parses a line

    col1Width   = 40
    col2Width   = 40
    colGE3Width = 40
    def parseLine(
            line,
            col1Width   = col1Width,
            col2Width   = col2Width,
            colGE3Width = colGE3Width,
            ):
        s = ''
        s += '{0:{1}s} '.format( line[0], col1Width )
        s += '{0:{1}s} '.format( line[1], col2Width )
        for component in line[2:]:
            s += '{0:{1}s} '.format( component, colGE3Width )
        return s


    c.cardTxtLines = []
    def w( text='' ):
        c.cardTxtLines.append( text + '\n' )


    w( c.imaxLine )
    w( c.jmaxLine )
    w( c.kmaxLine )
    
    w()

    for shapesLine in c.shapesLines:
        w( parseLine( shapesLine.split(), 7, 30, 45  ) )

    w()

    w( c.firstBinLine )
    w( c.observationLine )

    if all( hasattr( c, listName ) for listName in [ 'splitBinLine', 'splitProcessNamesLine', 'splitProcessNumbersLine', 'splitRateLine' ] ):
        w( parseLine( [ 'bin', '' ]     + c.splitBinLine ))
        w( parseLine( [ 'process', '' ] + c.splitProcessNamesLine ))
        w( parseLine( [ 'process', '' ] + c.splitProcessNumbersLine ))
        w( parseLine( [ 'rate', '' ]    + c.splitRateLine ))
    else:
        w( parseLine( [ 'bin', '' ]     + c.binLine.split()[1:] ))
        w( parseLine( [ 'process', '' ] + c.processNamesLine.split()[1:] ))
        w( parseLine( [ 'process', '' ] + c.processNumbersLine.split()[1:] ))
        w( parseLine( [ 'rate', '' ]    + c.rateLine.split()[1:] ))


    w()

    for discreteNuisanceLine in c.discreteNuisanceLines:
        w( discreteNuisanceLine )

    w()

    for nuis in c.nuisanceContainers:
        w( parseLine( [ nuis.name, nuis.type ] + nuis.vals ) )

    w()

    for paramNuisanceLine in c.paramNuisanceLines:
        w(paramNuisanceLine)

    if not writeToFile == None:
        with open( writeToFile, 'w' ) as outFp:
            outFp.writelines( c.cardTxtLines )
        print '[info] Wrote output to ' + normpath(writeToFile)




def RenameProcesses(
    process,
    indatacard,
    outdatacard=None,
    renameOutsideAcceptance=False,
    globalReplace=None,
    ):

    with open( indatacard, 'r' ) as indatacardFp:
        indatacardTxt = indatacardFp.read()


    # ======================================
    # Renaming

    outdatacardTxt = indatacardTxt

    # Replace multiple times to cover neighboring matches

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)InsideAcceptance_genPt_350p0_10000p0(\W)',
            r'\1{0}_PTH_GT350\2'.format(process),
            outdatacardTxt
            )

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)InsideAcceptance_genPt_(\d+)p0_(\d+)p0(\W)',
            r'\1{0}_PTH_\2_\3\4'.format(process),
            outdatacardTxt
            )


    # ======================================
    # In the old datacards there is still the word "to" instead of an underscore

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)InsideAcceptance_genPt_350p0to10000p0(\W)',
            r'\1{0}_PTH_GT350\2'.format(process),
            outdatacardTxt
            )

    for i in xrange(3):
        outdatacardTxt = re.sub(
            r'(\W)InsideAcceptance_genPt_(\d+)p0to(\d+)p0(\W)',
            r'\1{0}_PTH_\2_\3\4'.format(process),
            outdatacardTxt
            )

    if renameOutsideAcceptance:
        outdatacardTxt = re.sub(
            r'(\W)OutsideAcceptance(\W)',
            r'\1{0}_OutsideAcceptance\2'.format(process),
            outdatacardTxt
            )


    # ======================================
    # Process simple global replacements

    if not globalReplace == None:
        for string, replacement in globalReplace:
            outdatacardTxt = outdatacardTxt.replace( string, replacement )


    # ======================================
    # Write to file / return

    if outdatacard == 'auto':
        outdatacard = indatacard.replace( '.txt', '_renamedProcesses_{0}.txt'.format(datestr) )
        with open( outdatacard, 'w' ) as outdatacardFp:
            outdatacardFp.write( outdatacardTxt )
        print '[info] Wrote output to ' + outdatacard
    elif not outdatacard == None:
        outdatacard = join( TEMPCARDDIR, basename(indatacard).replace( '.txt', '_renamedProcesses.txt' ) )
        with open( outdatacard, 'w' ) as outdatacardFp:
            outdatacardFp.write( outdatacardTxt )
        print '[info] Wrote output to ' + outdatacard

    return outdatacardTxt






########################################
# End of Main
########################################
if __name__ == "__main__":
    main()