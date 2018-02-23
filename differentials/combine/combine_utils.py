import os, shutil
import logging
import os.path
import glob, re, copy
from collections import namedtuple

import differentials.core

GLOBAL_OUTDIR = 'out'
def set_global_outdir(new_dir):
    global GLOBAL_OUTDIR
    GLOBAL_OUTDIR = new_dir
def get_global_outdir():
    return GLOBAL_OUTDIR