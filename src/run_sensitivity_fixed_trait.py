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



odir = "outputs"
if not os.path.exists(odir):
    os.makedirs(odir)

# min, max, mean
lai = {}
lai["rf"] = (4.78, 6.94, 5.86)
lai["wsf"] = (3.46, 6.19, 4.83)
lai["dsf"] = (1.43, 4.75, 3.09)
lai["grw"] = (1.27, 3.39, 2.33)
lai["saw"] = (0.34, 1.67, 1.0)

depth = {}
depth["rf"] = 0.5
depth["wsf"] = 0.5
depth["dsf"] = 0.5
#depth["grw"] = 0.3
#depth["saw"] = 0.3
depth["grw"] = 0.5
depth["saw"] = 0.5

params = get_params()
#pfts = [pft_name]
pfts = list(params.columns)


time_step = 30
lat = -35.76
lon = 148.0



Tmax = 35.
RH = 10.
met = generate_met_data(Tmin=15, Tmax=Tmax, RH=RH, ndays=720,
                        lat=lat, lon=lon, time_step=time_step)
Dmax = np.max(met.vpd)
Dmean = np.mean(met.vpd)


names = ['pft','Tmax', 'Dmax', 'Dmean', 'gmin', 'lai', 'p50', 'Cl', 'Cs', \
         'depth', 'psi_stem', 'cwd', 'plc', 'day_of_death']
df = pd.DataFrame(columns=names)

for pft_name in pfts:
    print(pft_name)
    p = params[pft_name]
    lai_low, lai_high, lai_mu = lai[pft_name]


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

    F = Canopy(g1=g1, g0=g0, theta_J=theta_J, Rd25=Rd25, Q10=Q10,
               Vcmax25=Vcmax25, Jmax25=Jmax25, Eav=Eav,
               deltaSv=deltaSv, Eaj=Eaj, deltaSj=deltaSj)

    D = Desica(psi_stem0=psi_stem0, psi_f=psi_f, F=F, g1=g1, stop_dead=True,
               FAO=FAO, kp_sat=kp_sat, s50=s50, sf=sf, AL=AL)
    D.gmin = p.gmin
    D.p50 = p.p50
    D.Cl = p.Cl
    D.Cs = p.Cs
    D.soil_depth = depth[pft_name]


    out, day_of_death = D.run_simulation(met)
    psi_stem = out.psi_stem.iloc[-1]
    plc = out.plc.iloc[-1]
    cwd = out.cwd.iloc[-1]

    #plot_time_to_mortality(odir, out, time_step, to_screen=True)

    result = [pft_name, Tmax, Dmax, Dmean, p.gmin, AL, p.p50, p.Cl, p.Cs,  \
              depth[pft_name], psi_stem, cwd, plc, day_of_death]

    s = pd.Series(result, index=df.columns)
    df = df.append(s, ignore_index=True)

ofname = os.path.join(odir, "trait_sens_fixed.csv")
df.to_csv(ofname, index=False)
