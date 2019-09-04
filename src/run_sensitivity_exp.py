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
from get_params import get_params
import itertools
import multiprocessing as mp




def main(ncpus=None):

    params = get_params()
    pfts = list(params)

    pfts = ["rf"]

    if ncpus is None: # use them all!
        ncpus = mp.cpu_count()

    pool = mp.Pool(processes=ncpus)

    processes = []
    for i, pft_name in enumerate(pfts):
        #print(i, pft_name)
        px = params[pft_name]
        p = mp.Process(target=worker, args=(pft_name, px, ))
        processes.append(p)

    # Run processes
    for p in processes:
        p.start()

def worker(pft_name, p):

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
             'b', 'psi_e', 'depth', 'psi_stem', 'plc', 'day_of_death']
    df = pd.DataFrame(columns=names)


    # min, max, mean
    lai = {}
    lai["rf"] = (4.78, 6.94, 5.86)
    lai["wsf"] = (3.46, 6.19, 4.83)
    lai["dsf"] = (1.43, 4.75, 3.09)
    lai["grw"] = (1.27, 3.39, 2.33)
    lai["saw"] = (0.34, 1.67, 1.0)

    lai_low, lai_high, lai_mu = lai[pft_name]

    Tmax = 35.
    RH = 10.
    met = generate_met_data(Tmin=15, Tmax=Tmax, RH=RH, ndays=720,
                            lat=lat, lon=lon, time_step=time_step)
    Dmax = np.max(met.vpd)
    Dmean = np.mean(met.vpd)

    N = 5
    chg = 1.35
    total_exp = N**6 * (3**2)

    ranges = [
        np.linspace(p.gmin/chg, p.gmin*chg, N),  # gmin
        np.linspace(lai_low, lai_high,  N),      # AL
        np.linspace(p.p50/chg, p.p50*chg, N),    # p50
        np.linspace(p.Cl/chg, p.Cl*chg, N),      # Cl
        np.linspace(p.Cs/chg, p.Cs*chg, N),      # Cs
        np.linspace(0.1, 1.0, N),                # soil_depth
        np.array([2.79, 6.77, 10.39]),           # b - retension param: sand, sandy clay loam and silty clay
        np.array([-0.68*1e-03, \
                 -1.32*1e-03, -3.17*1e-03])      # psi_e
    ]

    count = 0
    last_progress = 9.0
    for gmin, AL, p50, Cl, Cs, soil_depth, b, psi_e in \
        itertools.product(*ranges):


        """
        (D.gmin, D.AL, D.p50, D.Cl, D.Cs,
         D.soil_depth, D.b,
         D.psi_e) = gmin, AL, p50, Cl, Cs, soil_depth, b, psi_e

        out, day_of_death = D.run_simulation(met)
        psi_stem = out.psi_stem.iloc[-1]
        plc = out.plc.iloc[-1]

        result = [Tmax, Dmax, Dmean, gmin, AL, p50, Cl, Cs,  \
                  b, psi_e, soil_depth, psi_stem, plc, day_of_death]

        s = pd.Series(result, index=df.columns)
        df = df.append(s, ignore_index=True)
        """

        progress = (count / total_exp) * 100.0
        
        if progress > last_progress:
            print(pft_name, "--", round(progress,3), count, ":", total_exp)
            last_progress += 9.
        count += 1


    odir = "outputs"
    if not os.path.exists(odir):
        os.makedirs(odir)

    ofname = os.path.join(odir, "%s_trait_sensitivity.csv" % (pft_name))
    df.to_csv(ofname, index=False)


if __name__ == "__main__":

    #main(ncpus=3)
    main()
