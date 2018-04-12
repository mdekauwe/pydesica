#!/usr/bin/env python
# coding: utf-8

"""
Extract data from GSWP3 forcing for a pixel that we can diurnalise to run
pydesica
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (06.04.2018)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import numpy as np
import xarray as xr
import constants as c

def main(met_dir):
    sec_2_3hr = 60. * 60. * 3.

    #random location
    row = 10
    col = 17

    f = open("GSWP3_met_%d_%d.csv" % (row,col), 'w')
    print("doy,tmin,tmax,rain,par,rh,lat,lon", file=f)

    for year in range(2000, 2011):
        #print(year)
        fname = os.path.join(met_dir,
                             "Tair/GSWP3.BC.Tair.3hrMap.%d.nc" % (year))
        tair_vals = open_file(fname, "Tair")

        # get lat, lon
        lat = tair_vals[0,row,col].lat.values
        lon = tair_vals[0,row,col].lon.values

        fname = os.path.join(met_dir,
                             "Qair/GSWP3.BC.Qair.3hrMap.%d.nc" % (year))
        qair_vals = open_file(fname, "Qair")

        fname = os.path.join(met_dir,
                             "SWdown/GSWP3.BC.SWdown.3hrMap.%d.nc" % (year))
        sw_vals = open_file(fname, "SWdown")

        fname = os.path.join(met_dir,
                             "Rainf/GSWP3.BC.Rainf.3hrMap.%d.nc" % (year))
        rain_vals = open_file(fname, "Rainf")

        doy = 1
        for i in range(0, tair_vals.time.size, 8):

            tmax = np.max(tair_vals[i:i+8,row,col].values) - c.DEG_2_KELVIN
            tmin = np.min(tair_vals[i:i+8,row,col].values) - c.DEG_2_KELVIN
            rain = np.sum(rain_vals[i:i+8,row,col].values * sec_2_3hr)
            par = np.max(sw_vals[i:i+8,row,col].values) * c.SW_2_PAR
            qair = np.max(qair_vals[i:i+8,row,col].values)

            # Get RH
            qair = qair_vals[i:i+8,row,col].values
            tair = tair_vals[i:i+8,row,col].values - c.DEG_2_KELVIN
            rh = np.mean(qair_to_rh(qair, tair))

            print("%d,%f,%f,%f,%f,%f,%f,%f" %
                    (doy,tmin,tmax,rain,par,rh,lat,lon), file=f)
            doy += 1

    f.close()

def open_file(fname, var):

    ds = xr.open_dataset(fname)
    data = ds[var]

    # Australia
    #lat_st = np.argwhere(data.lat.values == -43.75)[0][0]
    #lat_en = np.argwhere(data.lat.values == -10.25)[0][0]
    #lon_st = np.argwhere(data.lon.values == 112.75)[0][0]
    #lon_en = np.argwhere(data.lon.values == 153.75)[0][0]

    # NSW
    lat_st = np.argwhere(data.lat.values == -37.75)[0][0]
    lat_en = np.argwhere(data.lat.values == -27.75)[0][0]
    lon_st = np.argwhere(data.lon.values == 140.25)[0][0]
    lon_en = np.argwhere(data.lon.values == 154.25)[0][0]
    nsw = data[:,lat_st:lat_en,lon_st:lon_en]

    return nsw

def qair_to_rh(qair, tair, press=1013.25):
    es = 6.112 * np.exp((17.67 * tair) / (tair + 243.5))
    e = qair * press / (0.378 * qair + 0.622)
    rh = e / es
    rh = np.where(rh < 0.0, 0.0, rh)
    rh = np.where(rh > 1.0, 1.0, rh)

    return rh


if __name__ == "__main__":

    met_dir = "/g/data1/wd9/MetForcing/Global/GSWP3_2017"
    main(met_dir)
