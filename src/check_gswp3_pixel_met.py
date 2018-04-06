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



if __name__ == "__main__":

    time_step = 30
    keep_dry = False
    fname = "gswp3_met/GSWP3_met_2006.csv"
    df = pd.read_csv(fname)

    #plt.plot(df.tmax)
    #plt.plot(df.tmin)
    plt.plot(df.rain)
    #plt.plot(df.rh)
    plt.show()
    sys.exit()

    met = generate_met_data(Tmin=df.tmin[0], Tmax=df.tmax[0],
                            RH=df.rh[0]*100., PPFDmax=df.par[0],
                            precip=df.rain[0], ndays=1, time_step=time_step,
                            keep_dry=keep_dry)
    for i in range(1, len(df)):
        met_df = generate_met_data(Tmin=df.tmin[i], Tmax=df.tmax[i],
                                   RH=df.rh[i]*100., PPFDmax=df.par[i],
                                   precip=df.rain[i], ndays=1,
                                   time_step=time_step, keep_dry=keep_dry)

        met_df.day = i+1
        met = met.append(met_df)
    met = met.reset_index(drop=True)

    """
    #plt.plot(met.tair)
    #plt.plot(met.par)
    plt.plot(met.vpd)
    #plt.plot(met.precip)
    plt.xlim(0, 48)
    plt.show()
    """
    #plt.plot(met.tair)
    #plt.plot(met.par)
    plt.plot(met.vpd)
    #plt.plot(met.precip)
    plt.show()
