#!/usr/bin/env python
# coding: utf-8

"""
Join up sensitivity experiment files.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (08.09.2019)"
__email__ = "mdekauwe@gmail.com"

import pandas as pd
import sys
import numpy as np
import os
import glob as glob


pfts = ["rf", "wsf", "dsf", "grw", "saw"]

for pft in pfts:
    print(pft)
    df = pd.concat(map(pd.read_csv, glob.glob(os.path.join("outputs",
                                    "%s_trait_sensitivity_*_*.csv" % (pft)))))


    ofname = os.path.join("outputs", "%s_trait_sensitivity_all.csv" % (pft))
    df.to_csv(ofname, index=False)
