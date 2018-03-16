#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, pexpect, subprocess
from os.path import *

try:
    import AuthorizationAgent.AuthorizationAgent as AuthorizationAgent
    AUTO_PASSWORDS = True
except ImportError:
    AUTO_PASSWORDS = False


########################################
# Main
########################################

def main():

    import argparse
    parser = argparse.ArgumentParser()
    # parser.add_argument( '--string', type=str, default='default', help='default string' )
    parser.add_argument( '--pngs', action='store_true', help='boolean')
    parser.add_argument( '--typePassword', action='store_true', help='boolean')
    # parser.add_argument( '--list', metavar='N', type=str, nargs='+', help='list of strings' )
    args = parser.parse_args()

    if args.typePassword:
        global AUTO_PASSWORDS
        AUTO_PASSWORDS = False

    www_onlxplus = '/afs/cern.ch/user/t/tklijnsm/www/'
    plotdir_onlxplus = join( www_onlxplus, 'differentials2017/ptCombination_v4' )
    if plotdir_onlxplus.endswith('/'): plotdir_onlxplus = plotdir_onlxplus[:-1]


    if args.pngs:
        print '\nCreating .pngs out of all .pdfs'
        createPngsForAllPdfs( 'www/' )


    print '\nCopying in index.php in all subdirectories'
    backdir = abspath(os.getcwd())
    cmd = 'find www/ -type d -exec cp index.php {} \;'
    executeCommand( cmd )


    print '\nSyncing'
    # cmd = 'rsync -avWe ssh --delete-before www/ tklijnsm@lxplus.cern.ch:/afs/cern.ch/user/t/tklijnsm/www/differentials2017/ptCombination'
    cmd = 'rsync -avW --delete-before -e ssh www/ tklijnsm@lxplus.cern.ch:{0}'.format( plotdir_onlxplus )
    output = executeCommandWithPassword( cmd )
    # output = executeCommand( cmd )
    print output


    print '\nSetting permissions for afs webserver'
    cmd = 'afind %s -t d -e "fs setacl -dir {} -acl webserver:afs read"' % plotdir_onlxplus
    sshcmd = 'ssh tklijnsm@lxplus.cern.ch \'{0}\''.format( cmd )
    executeCommandWithPassword( sshcmd )
    # executeCommand( sshcmd )


def executeCommandWithPassword( cmd, expectString='Password:', verbose=False ):
    if not AUTO_PASSWORDS:
        return executeCommand( cmd, getOutput=True )

    child = pexpect.spawn( cmd )
    child.timeout = 300
    if verbose: print 'Waiting for expected line \'{0}\'...'.format( expectString )
    i = child.expect( [pexpect.TIMEOUT, expectString] )
    if i == 0:
        print "\nGot unexpected output: %s %s" % (child.before, child.after)
        return
    else:

        child.sendline( getpw() )

    output = child.read()
    return output


def executeCommand( cmd, getOutput=False ):
    if getOutput:
        output = subprocess.check_output( cmd, shell=True )
        return output.strip()
    else:
        os.system( cmd )


def createPngsForAllPdfs( directory ):

    allPdfs = []
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in [ f for f in filenames if f.endswith(".pdf") ]:
            allPdfs.append( join(dirpath, filename) )

    for pdf in allPdfs:
        png = pdf.replace( '.pdf', '.png' )
        cmd = 'convert -density 150 {0} -quality 90 -trim {1}'.format( pdf, png )
        executeCommand( cmd )

def getpw():
    return AuthorizationAgent.ReadKeyCipherFilePair( 
        keyFile='/mnt/t3nfs01/data01/shome/tklijnsm/.keys/key_Dec08_0fEuNc.key',
        pwFile='/mnt/t3nfs01/data01/shome/tklijnsm/.keys/pw_Dec08_Ed3fqK.pw'
        )




########################################
# End of Main
########################################
if __name__ == "__main__":
    main()