#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, argparse, re
from os.path import *


########################################
# Main
########################################

class Container:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


def main():

    parser = argparse.ArgumentParser()
    # parser.add_argument( '--string', type=str, default='default', help='default string' )
    # parser.add_argument( '--boolean', action='store_true', help='boolean')
    parser.add_argument( 'txtFiles', metavar='N', type=str, nargs='+', help='list of strings' )
    args = parser.parse_args()
    

    pat = r'\s*(\w+)\s:\s*([0-9\-\+\.]+)'

    containers = []
    for txtFile in args.txtFiles:

        with open( txtFile, 'r' ) as txtFp:
            lines = txtFp.readlines()


        container = Container()

        container.file = txtFile


        namematch = re.search( r'/([a-zA-Z0-9_]+)_recofit', abspath(txtFile) )
        if namematch:
            container.name = namematch.group(1)
        else:
            container.name = join( dirname(txtFile), basename(txtFile) )

        container.parDict = {}
        readingParameterLines = False
        for line in lines:
            
            if line.startswith( 'best fit parameter values:' ):
                readingParameterLines = True
                continue

            if line.startswith( 'Done in ' ):
                readingParameterLines = False
                continue

            if readingParameterLines:
                match = re.search( pat, line )

                if not match:
                    print 'Skipping line \'{0}\''.format(line)
                    continue

                container.parDict[ match.group(1).replace('ggH','smH').replace('xH','smH') ] = match.group(2)

        containers.append( container )


    def containerSorter( container ):
        if 'ggH' in container.name:
            return 1
        elif 'xH' in container.name:
            return 2
        elif 'merged' in container.name:
            return 3
        elif 'unsplit' in container.name:
            return 4
        else:
            return 5
    containers.sort( key=containerSorter )



    # ======================================
    # Parse a table

    def keySorter( key ):
        components = key.split('_')
        if len(components) >= 4:
            try:
                leftBound = float( components[3].replace('GT','') )

                match = re.search( r'cat([0-9]+)', key )
                if match:
                    catNumber = float(match.group(1))
                    leftBound += 0.1*catNumber

                return leftBound
            except:
                return -1
        else:
            return -1


    allKeys = []
    for container in containers:
        allKeys += container.parDict.keys()
    allKeys = list(set(allKeys))
    allKeys.sort( key=keySorter )

    table = []

    headerline = [ '{0:30}'.format('') ]
    for container in containers: headerline.append( '{0:12}'.format( container.name ) )
    table.append( ' | '.join(headerline) )


    for key in allKeys:
        line = []

        line.append( '{0:30}'.format(key) )

        for container in containers:
            line.append(
                '{0:12}'.format( container.parDict.get( key, 'x' ) )
                )

        table.append( ' | '.join(line) )


    print '\n'.join(table)





















########################################
# End of Main
########################################
if __name__ == "__main__":
    main()