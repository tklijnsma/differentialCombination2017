#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, sys, argparse
from os.path import *

sys.path.append('src')
import Commands

import LatestPaths


########################################
# Main
########################################

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument( 'pattern', type=str, default='*', nargs='?' )
    # parser.add_argument( '--boolean', action='store_true', help='boolean')
    # parser.add_argument( '--list', metavar='N', type=str, nargs='+', help='list of strings' )
    args = parser.parse_args()

    ListLatestPaths( args.pattern )


#____________________________________________________________________
def ListLatestPaths( pattern='*' ):

    allObjectsInLatestPaths = dir(LatestPaths)

    for obj in allObjectsInLatestPaths:

        if callable(obj):
            continue

        inode = getattr( LatestPaths, obj )

        if not isinstance( inode, basestring ):
            continue

        # if inode.startswith( 'suppliedInput' ):
        #     continue
        if inode.endswith( '.pyc' ):
            continue
        if not( isfile(inode) or isdir(inode) ):
            continue
        if not pattern == '*' and not pattern in obj + ' : ' + inode:
            continue

        print '{0:56} : {1}'.format( obj, inode )


########################################
# End of Main
########################################
if __name__ == "__main__":
    main()