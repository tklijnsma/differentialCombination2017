#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os
import differentials
import logging
import argparse


########################################
# Main
########################################

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument( 'scandirs', metavar='N', type=str, nargs='+', help='list of strings' )
    parser.add_argument( '--allq', action='store_true', help='boolean')
    parser.add_argument( '--shortq', action='store_true', help='boolean')
    parser.add_argument( '--longq', action='store_true', help='boolean')
    args = parser.parse_args()

    logging.getLogger().setLevel(logging.WARNING)

    for scandir in args.scandirs:
        accountant = differentials.scan_accounting.ScanAccountant(scandir)
        if args.allq:
            accountant.only_allq_recommendation()
        elif args.shortq:
            accountant.only_shortq_recommendation()
        elif args.longq:
            accountant.only_longq_recommendation()
        else:
            accountant.overview()


########################################
# End of Main
########################################
if __name__ == "__main__":
    main()