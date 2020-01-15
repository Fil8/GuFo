# Import modules

import os
import sys
import string
import numpy as np


import argparse
from  argparse import ArgumentParser
import textwrap as _textwrap


# get radiobs install directory
GUFO_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GUFO_DIR = GUFO_PATH+'/scavengers/'
sys.path.append(os.path.join(GUFO_PATH, 'scavengers'))

