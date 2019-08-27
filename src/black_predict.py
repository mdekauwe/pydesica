#!/usr/bin/env python
"""
Get Blackman predictions...
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (08.08.2019)"
__email__ = "mdekauwe@gmail.com"

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from get_params import get_params

params = get_params()
pfts = list(params)

for pft_name in pfts:

    p = params[pft_name]
    t_crit = -10**-3 * (p.Cs + p.Cl) * p.p50 / p.gmin
    print(pft_name, t_crit)
