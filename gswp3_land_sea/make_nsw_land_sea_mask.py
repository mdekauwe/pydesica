#!/usr/bin/env python
# coding: utf-8

"""
Make a land/sea mask for NSW
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (06.04.2018)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

fname = "gswp3_landmask_nomissing.nc"
ds = xr.open_dataset(fname)
data = ds["landsea"]

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
nsw = data[lat_st:lat_en,lon_st:lon_en].astype(np.int16)
nrows, ncols = nsw.shape
print(nrows, ncols, nrows*ncols)
data = nsw.values.reshape(nrows, ncols)
data.tofile("nsw_gswp3_land_sea_mask.bin")
