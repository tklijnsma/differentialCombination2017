#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os
from os.path import *
import argparse

import ROOT
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
ROOT.gSystem.AddIncludePath(" -I../HiggsAnalysis/CombinedLimit/interface/ ")

ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")

c = ROOT.TCanvas( 'c', 'c', 800, 600 )


########################################
# Main
########################################

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument( 'rootfile', type=str )
    parser.add_argument( 'varsToPrint', type=str, nargs='*' )
    parser.add_argument( '--MultiDimFit', action='store_true', help='boolean')
    parser.add_argument( '--list', action='store_true', help='boolean')
    parser.add_argument( '--printGetter', action='store_true', help='boolean')
    parser.add_argument( '--components', action='store_true', help='boolean')
    args = parser.parse_args()


    rootFp = ROOT.TFile.Open( args.rootfile )
    wNames = [ 'w', 'wsig_8TeV', 'wsig_13TeV', 'multipdf' ]
    foundWorkspace = False
    for wName in wNames:
        w = rootFp.Get( wName )
        if hasattr( w, 'allVars' ):
            foundWorkspace = True
            break
    if not foundWorkspace:
        print 'No easy-named workspace found; looping over objects'
        allObjs = [ key.GetName() for key in ROOT.gDirectory.GetListOfKeys() ]
        for wName in allObjs:
            w = rootFp.Get( wName )
            if isinstance( w, ROOT.RooWorkspace ):
                print 'Found RooWorkspace with name \'{0}\''.format( wName )
                foundWorkspace = True
                break
    if not foundWorkspace:
        print 'No RooWorkspace found; exiting'
        return


    if args.MultiDimFit:
        w.loadSnapshot('MultiDimFit')


    if args.components:


        componentArgset   = w.components()
        nComponents       = componentArgset.getSize()

        componentIterator = w.componentIterator()

        names = []
        for i in xrange(nComponents):
            element = componentIterator.Next()
            names.append( element.GetName() )

        names.sort()

        print 'List of all components in \'{0}\'\n'.format( wName )
        for name in names:
            print name

        return


    if args.list:

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

            print '\n' + '='*70
            print 'All elements in {0}.{1}():'.format( wName, method )

            argset   = getattr( w, method )()

            names = []
            try:
                iterator = argset.createIterator()
                for i in xrange(argset.getSize()):
                    element = iterator.Next()
                    names.append( element.GetName() )
            except AttributeError:
                item = argset.begin()
                for i in xrange(len(argset)):
                    names.append( item.GetName() )
                    try:
                        next(item)
                    except StopIteration:
                        break

            names.sort()
            for name in names:
                print name

        return



    if len(args.varsToPrint) == 0:
        print '\nDoing {0}.Print():'.format( wName )
        w.Print()
            
    else:
        for thingName in args.varsToPrint:

            print '\nDoing {0}.Print():'.format( thingName )
            for getter in [ 'pdf', 'function', 'var', 'cat', 'set', 'genobj' ]:
                try:
                    thing = getattr( w, getter )( thingName )
                    thing.Print()
                    if args.printGetter: print '    (Successful getter was \'{0}\')'.format(getter)
                except ReferenceError:
                    continue
                break

    rootFp.Close()



########################################
# End of Main
########################################
if __name__ == "__main__":
    main()