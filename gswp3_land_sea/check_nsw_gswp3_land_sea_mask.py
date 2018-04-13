#!/usr/bin/env python
# coding: utf-8

"""
Check if pixel is land or sea
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (13.04.2018)"
__email__ = "mdekauwe@gmail.com"

import sys
import numpy as np

nrows = 20
ncols = 28

row = int(sys.argv[1])
col = int(sys.argv[2])
fname = str(sys.argv[3])
data = np.fromfile(fname, dtype=np.int16).reshape(nrows, ncols)
print(data[row,col])
