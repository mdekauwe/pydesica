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

from scipy.optimize import curve_fit


def sigmoid(x, a, b):
    y = 1.0 / (1.0 + np.exp(-b * (x - a)))
    return y

def func(x, a):
    g1 = 4.0
    return g1 * np.exp(a * (x+np.min(x)))

if __name__ == "__main__":

    fname = "outputs/drydown_out.csv"
    df = pd.read_csv(fname)

    psi_pd = []
    gsw_pd = []
    for i in range(0, len(df), 48):

        day_psi_soil = df.psi_soil[i:i+48].values
        day_gsw = df.gsw[i:i+48].values

        if np.isnan(np.sum(day_gsw)) == False:

            if len(day_gsw) > 1:
                day_gsw = np.nanmean(day_gsw)
                if np.isnan(day_gsw) == False and \
                   np.isnan(day_psi_soil[12]) == False:
                    gsw_pd.append( day_gsw )
                    psi_pd.append( day_psi_soil[12] ) # 6am

    psi_pd = np.asarray(psi_pd)
    gsw_pd = np.asarray(gsw_pd)
    print(np.max(gsw_pd))
    #popt, pcov = curve_fit(func, psi_pd, gsw_pd)
    popt, pcov = curve_fit(sigmoid, psi_pd, gsw_pd)
    print(popt[0])
    plt.plot(psi_pd, gsw_pd, "ko")
    #plt.plot(psi_pd, func(psi_pd, popt[0]), 'r-')
    plt.plot(psi_pd, sigmoid(psi_pd, popt[0], popt[1]), 'r-')

    plt.show()
