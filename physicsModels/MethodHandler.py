#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os
import importlib

########################################
# Main
########################################

def flag_as_method(function):
    function.is_method = True
    return function

class MethodHandler(object):
    """docstring for MethodHandler"""

    registered_functions = {}

    def __init__(self, modules):
        self.process_modules(modules)

    def process_modules(self, modules):
        for module in modules:
            self.process_module(module)

    def process_module(self, module):
        mod = importlib.import_module(module)
        for attr in dir(mod):
            fn = getattr(mod, attr)
            if hasattr(fn, 'is_method'):
                self.registered_functions[fn.__name__] = fn

    def make_methods(self, cls):
        for fn_name, fn in self.registered_functions.iteritems():
            setattr(cls, fn_name, fn)