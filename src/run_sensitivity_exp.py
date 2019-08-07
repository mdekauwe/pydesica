#!/usr/bin/env python
# coding: utf-8

"""
Explore how traits affect vulnerability to reaching critical water potential

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (07.08.2019)"
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
from get_params import get_params
import itertools

names = ['Tmax', 'RH', 'gmin', 'AL', 'p50', 'Cl', 'Cs', 'day_of_death']
df = pd.DataFrame(columns=names)

params = get_params()
pfts = list(params)
p = params["rf"]

time_step = 30
lat = -35.76
lon = 148.0
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
Kplant = p.Kplant
kp_sat = p.kpsat
g1 = p.g1

Tmaxx = [30, 35, 40, 45]
RHx = [30, 20, 10, 5]

N = 3
ranges = [
    np.linspace(5, 15, N),        # gmin
    np.linspace(0.5, 4, N),       # AL
    np.linspace(-1, -6, N),       # p50
    np.linspace(200, 800, N),     # Cl
    np.linspace(10000, 50000, N), # Cs
]

for Tmax in Tmaxx:
    for RH in RHx:
        met = generate_met_data(Tmin=15, Tmax=Tmax, RH=RH, ndays=500,
                                lat=lat, lon=lon, time_step=time_step)

        for gmin, AL, p50, Cl, Cs in itertools.product(*ranges):

            F = Canopy(g1=g1, g0=g0, theta_J=theta_J, Rd25=Rd25, Q10=Q10,
                       Vcmax25=Vcmax25, Jmax25=Jmax25, Eav=Eav,
                       deltaSv=deltaSv, Eaj=Eaj, deltaSj=deltaSj)

            D = Desica(psi_stem0=psi_stem0, AL=AL, p50=p50, psi_f=psi_f,
                       gmin=gmin, Cl=Cl, Cs=Cs, F=F, g1=g1, stop_dead=True,
                       FAO=FAO, kp_sat=kp_sat)

            out, day_of_death = D.run_simulation(met)
            result = [Tmax, RH, gmin, AL, p50, Cl, Cs, day_of_death]
            print(result)

            s = pd.Series(result, index=df.columns)
            df = df.append(s, ignore_index=True)


odir = "outputs"
if not os.path.exists(odir):
    os.makedirs(odir)

ofname = os.path.join(odir, "trait_sensitivity.csv")
df.to_csv(ofname, index=False)
