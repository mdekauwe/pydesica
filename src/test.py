#!/usr/bin/env python

"""
Test script, with constant met drivers we get steady state after some time, in
which case the water potential for the leaf follows:

Eleaf = kplant(psi_soil - psi_leaf)
or psi_leaf = psi_soil - Eleaf / kplant

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (23.11.2017)"
__email__ = "mdekauwe@gmail.com"

import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
from canopy import Canopy, FarquharC3
from desica import Desica
from math import isclose

day = np.repeat([1,2,3,4,5], 96)
N = len(day)
par = np.ones(N) * 1500.0
tair = np.ones(N) * 25.0
vpd = np.ones(N) * 2.0
precip = np.ones(N) * 0.0
Ca = np.ones(N) * 400.0
press = np.ones(N) * 101.0
met = pd.DataFrame({'day':day, 'par':par, 'tair':tair,
                    'vpd':vpd, 'precip':precip, 'press':press,
                    'Ca':Ca})

psi_stem0 = 0. # initial stem water potential, MPa
AL = 6.        # plant leaf area, m2
p50 = -4.      # xylem pressure inducing 50% loss of hydraulic conductivity
               # due to embolism, MPa
psi_f = -3.    # reference potential for Tuzet model, MPa
gmin = 10.     # minimum stomatal conductance, mmol m-2 s-1
Cl = 10000.    # leaf capacitance, mmol MPa-1 (total plant)
Cs = 120000.   # stem capacitance, mmol MPa-1
g1 = 4.0       # sensitivity of stomatal conductance to the assimilation
               # rate, kPa

F = Canopy(g1=g1)
D = Desica(psi_stem0=psi_stem0, AL=AL, p50=p50, psi_f=psi_f, gmin=gmin,
           Cl=Cl, Cs=Cs, F=F, g1=g1, nruns=3, stop_dead=True, keep_wet=True)
out = D.run_simulation(met)

# last row
steady = out.tail(1)


x = (steady.psi_soil - steady.Eleaf / steady.kplant).values
y = steady.psi_leaf.values

if isclose(x, y, rel_tol=0.01):
    print("Match", x[0], y[0])
else:
    print(x[0], y[0])
