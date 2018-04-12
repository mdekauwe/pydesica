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

    time_step = 30
    keep_dry = False
    fname = "gswp3_met/GSWP3_met_10_17.csv"
    df = pd.read_csv(fname)

    years = np.unique(df.year)
    for year in years:
        dfx = df[df.year == year]

        met = generate_met_data(year=year, Tmin=dfx.tmin.iloc[0], Tmax=dfx.tmax.iloc[0],
                                RH=dfx.rh.iloc[0]*100., PPFDmax=df.par.iloc[0],
                                precip=dfx.rain.iloc[0], lat=dfx.lat.iloc[0],
                                lon=dfx.lon.iloc[0], ndays=1, time_step=time_step,
                                keep_dry=keep_dry)

        for i in range(1, len(dfx)):

            met_df = generate_met_data(year=year, Tmin=dfx.tmin.iloc[i],
                                       Tmax=dfx.tmax.iloc[i], RH=dfx.rh.iloc[i]*100.,
                                       PPFDmax=dfx.par.iloc[i], precip=dfx.rain.iloc[i],
                                       lat=dfx.lat.iloc[0], lon=dfx.lon.iloc[0], ndays=1,
                                       time_step=time_step, keep_dry=keep_dry)

            met_df.day = i+1
            met = met.append(met_df)
        met = met.reset_index(drop=True)

        # first year, initialise
        if year == years[0]:
            psi_leaf0 = -1.
            psi_stem0 = -0.5
            sw0 = 0.4
        # else use the last point from previous year
        else:
            psi_leaf0 = out.psi_leaf.iloc[-1]
            psi_stem0 = out.psi_stem.iloc[-1]
            sw0 = out.sw.iloc[-1]

        print(year, round(psi_leaf0, 2), round(psi_stem0, 2), round(sw0, 2))

        AL = 3.        # plant leaf area, m2
        p50 = -4.      # xylem pressure inducing 50% loss of hydraulic
                       # conductivity due to embolism, MPa
        psi_f = -3.    # reference potential for Tuzet model, MPa
        gmin = 10.     # minimum stomatal conductance, mmol m-2 s-1
        Cl = 10000.    # leaf capacitance, mmol MPa-1 (total plant)
        Cs = 120000.   # stem capacitance, mmol MPa-1
        g1 = 4.0       # sensitivity of stomatal conductance to the assimilation
                       # rate, kPa
        g0 = 0.0
        theta_J = 0.85
        Rd25 = 0.92
        Q10 = 1.92
        Vcmax25 = 50.0
        Jmax25 = 100.
        Eav = 58550.0
        deltaSv = 629.26
        Eaj = 29680.
        deltaSj = 631.88
        soil_depth = 2.0
        FAO = False
        F = Canopy(g1=g1, g0=g0, theta_J=theta_J, Rd25=Rd25, Q10=Q10,
                   Vcmax25=Vcmax25, Jmax25=Jmax25, Eav=Eav, deltaSv=deltaSv,
                   Eaj=Eaj, deltaSj=deltaSj)
        D = Desica(psi_stem0=psi_stem0, AL=AL, p50=p50, psi_f=psi_f, gmin=gmin,
                   Cl=Cl, Cs=Cs, F=F, g1=g1, run_twice=True,
                   soil_depth=soil_depth, sw0=sw0, stop_dead=True, FAO=FAO)
        out, day_of_death = D.run_simulation(met)

        odir = "plots"
        if not os.path.exists(odir):
            os.makedirs(odir)

        plot_time_to_mortality(odir, out, time_step, year)
        plot_transpiration(odir, out, year)
        plot_transpiration_and_pet(odir, out, year)
        plot_cwd(odir, out, time_step, year)
        plot_sw(odir, out, time_step, year)

        odir = "outputs"
        if not os.path.exists(odir):
            os.makedirs(odir)

        ofname = os.path.join(odir, "desica_out_%d.csv" % (year))
        out.to_csv(ofname, index=False)
