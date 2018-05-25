#!/usr/bin/env python
# coding: utf-8

"""
Apply desica to a single pixel of GSWP3 forcing

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (06.03.2018)"
__email__ = "mdekauwe@gmail.com"

import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
from generate_met_data import generate_met_data
from canopy import Canopy, FarquharC3
from math import isclose
from calc_pet import calc_net_radiation, calc_pet_energy
import constants as c
from desica import Desica
from desica import plot_time_to_mortality
from desica import plot_transpiration
from desica import plot_cwd
from desica import plot_sw
from desica import plot_transpiration_and_pet


if __name__ == "__main__":

    fname = "outputs//drydown_out.csv"
    df = pd.read_csv(fname)

    for i in range(len(out), 48):

        day.psi_soil = out.psi_soil[i:i+48]
        print(len
        print(out.t[i], out.t[i] / 48., out.psi_soil[i])

        if i > 300:
            sys.exit()
