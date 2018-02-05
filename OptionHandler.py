#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os
import importlib
import argparse


########################################
# Main
########################################

def flag_as_option(function):
    function.is_option = True
    return function

class OptionHandler(object):
    """docstring for OptionHandler"""

    registered_functions = {}

    def __init__(self, modules=[]):
        self.parser = argparse.ArgumentParser()
        self.process_modules(modules)

    def set_parser(self, parser):
        self.parser = parser

    def process_modules(self, modules):
        for module in modules:
            self.process_module(module)

    def process_module(self, module):
        mod = importlib.import_module(module)
        for attr in dir(mod):
            if hasattr( getattr(mod, attr), 'is_option' ):
                self.make_option( getattr(mod, attr) )

    def make_option(self, function):
        self.parser.add_argument('--{0}'.format(function.__name__), action='store_true')
        self.registered_functions[function.__name__] = function
        return function

    def parse(self):
        self.args = self.parser.parse_args()

    def execute_functions(self):
        for name, function in self.registered_functions.iteritems():
            if getattr(self.args, name):
                function(self.args)
