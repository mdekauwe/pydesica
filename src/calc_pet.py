#!/usr/bin/env python

"""
Calculate potential evapotranspiration

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (05.04.2018)"
__email__ = "mdekauwe@gmail.com"

import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from generate_met_data import generate_met_data

def calc_pet_energy(rnet, G=0.0):

    # Energy-only PET (mm 30 s-1), based on Milly et al. 2016
    # rne-G is in MJ m-2 s-1
    pet = 0.8 * rnet - G

    return pet

def calc_net_radiation(sw_rad, tair, albedo=0.15):

    # Net loss of long-wave radn, Monteith & Unsworth '90, pg 52, eqn 4.17
    net_lw = 107.0 - 0.3 * tair # W m-2

    # Net radiation recieved by a surf, Monteith & Unsw '90, pg 54 eqn 4.21
    #    - note the minus net_lw is correct as eqn 4.17 is reversed in
    #      eqn 4.21, i.e Lu-Ld vs. Ld-Lu
    #    - NB: this formula only really holds for cloudless skies!
    #    - Bounding to zero, as we can't have negative soil evaporation, but you
    #      can have negative net radiation.
    #    - units: W m-2
    net_rad = np.maximum(0.0, (1.0 - albedo) * sw_rad - net_lw)

    return net_rad


if __name__ == "__main__":

    J_TO_MJ = 1.0E-6
    PAR_2_SW = 1.0 / 2.3
    SEC_2_HLFHR = 1800.
    time_step = 30
    met = generate_met_data(Tmin=10, RH=30, ndays=1, time_step=time_step)

    sw_rad = met.par * PAR_2_SW
    rnet = calc_net_radiation(sw_rad, met.tair, albedo=0.15)
    # W m-2 -> MJ m-2 s-1
    rnet *= J_TO_MJ
    
    pet = calc_pet_energy(rnet)
    print(np.sum(pet * SEC_2_HLFHR))
