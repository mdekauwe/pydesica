#!/usr/bin/env python
# coding: utf-8

"""
Run each PFT

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (19.09.2019)"
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
#from old_desica import Desica
from desica import Desica
from desica import plot_time_to_mortality
import itertools
import multiprocessing as mp
import random


params = pd.read_csv("outputs/params.csv", index_col=None)
params.index = params["trait"]

pfts = ["rf", "wsf", "dsf", "grw", "saw"]
#pfts = ["rf"]
for pft in pfts:
    print(pft)
    p = params[pft]

    #
    ## Generate trait space...
    #
    # min, max, mean
    lai = {}
    lai["rf"] = (4.78, 6.94, 5.86)
    lai["wsf"] = (3.46, 6.19, 4.83)
    lai["dsf"] = (1.43, 4.75, 3.09)
    lai["grw"] = (1.27, 3.39, 2.33)
    lai["saw"] = (0.34, 1.67, 1.0)

    lai_low, lai_high, lai_mu = lai[pft]

    lat = -35.76
    lon = 148.0
    Tmax = 35.
    RH = 10.
    time_step = 30

    met = generate_met_data(Tmin=15, Tmax=Tmax, RH=RH, ndays=720,
                            lat=lat, lon=lon, time_step=time_step)
    Dmax = np.max(met.vpd)
    Dmean = np.mean(met.vpd)

    g0 = 0.0
    theta_J = 0.85
    Rd25 = 0.92
    Q10 = 1.92
    Vcmax25 = p.Vcmax
    Jmax25 = p.Jmax
    Eav = 58550.0
    deltaSv = 629.26
    Eaj = 29680.
    deltaSj = 631.88
    FAO = False
    psi_stem0 = -0.5
    psi_f = p.psiv
    kp_sat = p.kpsat
    g1 = p.g1
    s50 = p.s50
    sf = p.sf
    AL = lai_mu
    psi_e = -1.32 * c.KPA_2_MPA # Sandy clay loam, MPa
    b = 6.77
    #psi_e = -3.17 * c.KPA_2_MPA  # Silty clay clay loam, MPa
    #b = 10.39                    # Silty clay, SW retention curve param

    F = Canopy(g1=g1, g0=g0, theta_J=theta_J, Rd25=Rd25, Q10=Q10,
               Vcmax25=Vcmax25, Jmax25=Jmax25, Eav=Eav,
               deltaSv=deltaSv, Eaj=Eaj, deltaSj=deltaSj)

    D = Desica(psi_stem0=psi_stem0, psi_f=psi_f, F=F, g1=g1, stop_dead=True,
               FAO=FAO, kp_sat=kp_sat, s50=s50, sf=sf, AL=AL,
               force_refilling=False)

    out, day_of_death = D.run_simulation(met)

    #odir = "/Users/mdekauwe/Desktop/refilling_plots"
    odir = "/Users/mdekauwe/Desktop/new_plots"
    #odir = "/Users/mdekauwe/Desktop/old_plots"
    if not os.path.exists(odir):
        os.makedirs(odir)

    plot_time_to_mortality(odir, out, time_step, to_screen=False, pft=pft)
