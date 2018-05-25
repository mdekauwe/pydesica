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

    # Examine how time to death changes as gmin is increased?
    do_sensitivity = False
    lat = -35.76
    lon = 148.0
    met = generate_met_data(Tmin=15, Tmax=35.0, RH=30, ndays=200,
                            lat=lat, lon=lon, time_step=time_step)

    psi_stem0 = 0. # initial stem water potential, MPa
    AL = 3.        # plant leaf area, m2
    p50 = -4.      # xylem pressure inducing 50% loss of hydraulic conductivity
                   # due to embolism, MPa
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
    FAO = False
    year = 2000
    F = Canopy(g1=g1, g0=g0, theta_J=theta_J, Rd25=Rd25, Q10=Q10,
               Vcmax25=Vcmax25, Jmax25=Jmax25, Eav=Eav, deltaSv=deltaSv,
               Eaj=Eaj, deltaSj=deltaSj)
    D = Desica(psi_stem0=psi_stem0, AL=AL, p50=p50, psi_f=psi_f, gmin=gmin,
               Cl=Cl, Cs=Cs, F=F, g1=g1, run_twice=True, stop_dead=True,
               FAO=FAO)
    out, day_of_death = D.run_simulation(met)

    odir = "outputs"
    if not os.path.exists(odir):
        os.makedirs(odir)

    ofname = os.path.join(odir, "drydown_out.csv")
    out.to_csv(ofname, index=False)
