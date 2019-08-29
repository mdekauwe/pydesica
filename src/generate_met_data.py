#!/usr/bin/env python
# coding: utf-8

"""
Generate some fake met data for testing. Replication of Remko's R module

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
import constants as c

def generate_met_data(year=2000, PPFDmax=2000, RH=30, Tmax=30, Tmin=10,
                      day_length=12, sunrise=8, time_step=15, lag=0.5,
                      ndays=40, lat=-35.76, lon=148.0, co2_conc=400.,
                      air_press=101.0, precip=None, keep_dry=False):

    day_2_sec = 1.0 / (60. * 60. * 24.)
    hlhr_2_sec = 1.0 / 1800.0

    if time_step == 15:
        nx = 96
    elif time_step == 30:
        nx = 48
    elif time_step == 60:
        nx = 24
    elif time_step == 120:
        nx = 12

    r = np.arange(0, (24 * 60 - time_step) + time_step, time_step)
    p = np.zeros(len(r))
    ta = np.zeros(len(r))
    ii = np.argwhere((r > (sunrise*60)) & (r <= (sunrise*60 + day_length*60)))
    relt = (r[ii] - sunrise * 60) / float(day_length * 60)
    p[ii] = generate_par(relt, PPFDmax)

    timeh = np.arange(0, (day_length - time_step / 60.) + 0.25, time_step / 60.)
    td = par_to_diurnal(timeh, Tmax, Tmin, day_length, lag)
    for i,ival in enumerate(ii):
        ta[ival-1] = td[i]

    nnight = len(r) - len(td)
    tdec = ta[max(ii)-1] + (Tmin - ta[max(ii)-1]) * \
            np.arange(1, nnight+1) / nnight
    after = np.arange((max(ii)), len(r))
    ta[after] = tdec[np.arange(len(after))]
    before = np.arange(min(ii)-1)
    ta[before] = tdec[np.arange(len(after), len(tdec))]

    vpd, ea = rh_to_vpd(RH, ta)

    day = np.repeat(np.arange(1, ndays+1), nx)
    new_size = int(len(day) / nx)
    vpd = np.tile(vpd, new_size)
    ea = np.tile(ea, new_size)
    tair = np.tile(ta, new_size)
    par = np.tile(p, new_size)
    sw_rad = par * c.PAR_2_SW

    if keep_dry:
        precip = np.zeros(len(par))
    else:
        if precip is None:
            precip = np.zeros(len(par))
        else:
            if precip > 0.0:
                #print(precip)
                precipx = np.ones(len(par)) * precip
                precip = precipx / nx * hlhr_2_sec
                #print(precip)
                #print(np.sum(precip*1800.))
            else:
                precip = np.zeros(len(par))

    Ca = np.ones(len(par)) * co2_conc # umol mol-1
    press = np.ones(len(par)) * air_press # kPa
    lat = np.ones(len(par)) * lat
    lon = np.ones(len(par)) * lon
    oyear = np.ones(len(par)) * year
    met = pd.DataFrame({'year':oyear, 'day':day, 'par':par, 'tair':tair,
                        'vpd':vpd, 'precip':precip, 'press':press,
                        'Ca':Ca, 'ea':ea, 'sw_rad':sw_rad, 'lat':lat,
                        'lon':lon})

    return met

def rh_to_vpd(rh, tair, pressure=101.0):
    esat_val = esat(tair, pressure)
    e = (rh / 100.0) * esat_val
    vpd = (esat_val - e) / 1000.0
    return vpd, e

def esat(tair, pressure=101.0):
    a = 611.21
    b = 17.502
    c = 240.97
    f = 1.0007 + 3.46 * 10**-8 * pressure * 1000
    esat_val = f * a * (np.exp(b * tair / (c + tair)))
    return esat_val

def generate_par(relt, ymax):
    return (ymax / 2.0) * (1.0 + np.sin((2.0 * np.pi * relt) + 1.5 * np.pi))

def par_to_diurnal(timeh, Tmax, Tmin, daylen, lag):
    return (Tmax - Tmin) * np.sin((np.pi * timeh) / (daylen + 2 * lag)) + Tmin

if __name__ == "__main__":

    met = generate_met_data(Tmin=10, RH=30, ndays=200, time_step=30)

    plt.plot(met.vpd)
    plt.xlim(0, 48)
    plt.show()
