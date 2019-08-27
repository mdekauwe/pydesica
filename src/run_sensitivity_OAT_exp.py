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

def main(pft_name, p, ranges):

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
    kp_sat = p.kpsat
    g1 = p.g1
    s50 = p.s50
    sf = p.sf
    AL = 2.0

    F = Canopy(g1=g1, g0=g0, theta_J=theta_J, Rd25=Rd25, Q10=Q10,
               Vcmax25=Vcmax25, Jmax25=Jmax25, Eav=Eav,
               deltaSv=deltaSv, Eaj=Eaj, deltaSj=deltaSj)

    D = Desica(psi_stem0=psi_stem0, psi_f=psi_f, F=F, g1=g1, stop_dead=True,
               FAO=FAO, kp_sat=kp_sat, s50=s50, sf=sf, AL=AL)

    names = ['Tmax', 'Dmax', 'Dmean', 'gmin', 'lai', 'p50', 'Cl', 'Cs', \
             'depth', 'psi_stem', 'cwd', 'plc', 'day_of_death']
    df = pd.DataFrame(columns=names)



    Tmax = 35.
    RH = 10.
    met = generate_met_data(Tmin=15, Tmax=Tmax, RH=RH, ndays=720,
                            lat=lat, lon=lon, time_step=time_step)
    Dmax = np.max(met.vpd)
    Dmean = np.mean(met.vpd)

    N = 5
    chg = 1.3
    total_exp = 15625

    count = 0
    last_progress = 9.0
    for gmin, AL, p50, Cl, Cs, soil_depth, in \
        itertools.product(*ranges):


        #"""
        (D.gmin, D.AL, D.p50, D.Cl, D.Cs,
         D.soil_depth) = gmin, AL, p50, Cl, Cs, soil_depth

        out, day_of_death = D.run_simulation(met)
        psi_stem = out.psi_stem.iloc[-1]
        plc = out.plc.iloc[-1]
        cwd = out.cwd.iloc[-1]

        result = [Tmax, Dmax, Dmean, gmin, AL, p50, Cl, Cs,  \
                  soil_depth, psi_stem, cwd, plc, day_of_death]

        s = pd.Series(result, index=df.columns)
        df = df.append(s, ignore_index=True)
        #"""

        count += 1
        progress = (count / total_exp) * 100.0

        if progress > last_progress:
            print(pft_name, "--", round(progress,3), count, ":", total_exp)
            last_progress += 9.


    return df

if __name__ == "__main__":

    odir = "outputs"
    if not os.path.exists(odir):
        os.makedirs(odir)

    # Expecting PFT to be supplied on cmd line, e.g.
    # $ python src/run_sensitivity_OAT_exp.py "rf"
    if len(sys.argv) < 2:
        raise TypeError("Expecting pft name to be supplied on cmd line!")
    pft_name = sys.argv[1]

    # min, max, mean
    lai = {}
    lai["rf"] = (4.78, 6.94, 5.86)
    lai["wsf"] = (3.46, 6.19, 4.83)
    lai["dsf"] = (1.43, 4.75, 3.09)
    lai["grw"] = (1.27, 3.39, 2.33)
    lai["saw"] = (0.34, 1.67, 1.0)

    params = get_params()
    pfts = [pft_name]
    p = params[pft_name]
    lai_low, lai_high, lai_mu = lai[pft_name]

    chg = 1.5
    N = 1
    NN = 40

    # gmin
    print("gmin")
    ranges = [
        np.linspace(p.gmin/chg, p.gmin*chg, NN), # gmin
        np.linspace(lai_mu, lai_mu,  N),         # AL
        np.linspace(p.p50/chg, p.p50*chg, N),    # p50
        np.linspace(p.Cl/chg, p.Cl*chg, N),      # Cl
        np.linspace(p.Cs/chg, p.Cs*chg, N),      # Cs
        np.linspace(1.0, 1.0, N)                 # soil_depth
    ]
    df1 = main(pft_name, p, ranges)

    # LAI
    print("lai")
    ranges = [
        np.linspace(p.gmin/chg, p.gmin*chg, N),  # gmin
        np.linspace(lai_low, lai_high,  NN),     # AL
        np.linspace(p.p50/chg, p.p50*chg, N),    # p50
        np.linspace(p.Cl/chg, p.Cl*chg, N),      # Cl
        np.linspace(p.Cs/chg, p.Cs*chg, N),      # Cs
        np.linspace(1.0, 1.0, N)                 # soil_depth
    ]
    df2 = main(pft_name, p, ranges)

    # p50
    print("p50")
    ranges = [
        np.linspace(p.gmin/chg, p.gmin*chg, N),  # gmin
        np.linspace(lai_mu, lai_mu,  N),         # AL
        np.linspace(p.p50/chg, p.p50*chg, NN),   # p50
        np.linspace(p.Cl/chg, p.Cl*chg, N),      # Cl
        np.linspace(p.Cs/chg, p.Cs*chg, N),      # Cs
        np.linspace(1.0, 1.0, N)                 # soil_depth
    ]
    df3 = main(pft_name, p, ranges)

    # Cl
    print("Cl")
    ranges = [
        np.linspace(p.gmin/chg, p.gmin*chg, N),  # gmin
        np.linspace(lai_mu, lai_mu,  N),         # AL
        np.linspace(p.p50/chg, p.p50*chg, N),    # p50
        np.linspace(p.Cl/chg, p.Cl*chg, NN),     # Cl
        np.linspace(p.Cs/chg, p.Cs*chg, N),      # Cs
        np.linspace(1.0, 1.0, N)                 # soil_depth
    ]
    df4 = main(pft_name, p, ranges)

    # Cs
    print("Cs")
    ranges = [
        np.linspace(p.gmin/chg, p.gmin*chg, N),  # gmin
        np.linspace(lai_mu, lai_mu,  N),         # AL
        np.linspace(p.p50/chg, p.p50*chg, N),    # p50
        np.linspace(p.Cl/chg, p.Cl*chg, N),      # Cl
        np.linspace(p.Cs/chg, p.Cs*chg, NN),     # Cs
        np.linspace(1.0, 1.0, N)                 # soil_depth
    ]
    df5 = main(pft_name, p, ranges)

    # depth
    print("soil depth")
    ranges = [
        np.linspace(p.gmin/chg, p.gmin*chg, N),  # gmin
        np.linspace(lai_mu, lai_mu,  N),         # AL
        np.linspace(p.p50/chg, p.p50*chg, N),    # p50
        np.linspace(p.Cl/chg, p.Cl*chg, N),      # Cl
        np.linspace(p.Cs/chg, p.Cs*chg, N),      # Cs
        np.linspace(0.1, 1.0, NN)               # soil_depth
    ]
    df6 = main(pft_name, p, ranges)

    df = pd.concat([df1, df2, df3, df4, df5, df6])
    ofname = os.path.join(odir, "%s_trait_sens_OAT.csv" % (pft_name))
    df.to_csv(ofname, index=False)
