#!/usr/bin/env python
# coding: utf-8

"""
Plot the pixel output

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (12.03.2018)"
__email__ = "mdekauwe@gmail.com"

import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import os

import constants as c


def plot_time_to_mortality(odir, out, timestep=15, year=None):

    cb = ['#377eb8', '#ff7f00', '#4daf4a', \
          '#f781bf', '#a65628', '#984ea3',\
          '#999999', '#e41a1c', '#dede00']

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()

    ln1 = ax1.plot(out.psi_leaf, ls="-", color=cb[2], label="Leaf")
    ln2 = ax1.plot(out.psi_stem, ls="-", color=cb[1], label="Stem")
    ln3 = ax1.plot(out.psi_soil, ls="-", color=cb[0], label="Soil")
    ln4 = ax2.plot(out.plc, ls='-', color=cb[6],
                   label="PLC")

    # added these three lines
    lns = ln1 + ln2 + ln3 + ln4
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=(0.05,0.08), ncol=2)
    #ax1.legend(numpoints=1, loc="best")
    ax2.set_ylabel(r'PLC (%)')
    #ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("Water potential (MPa)")
    if year is None:
        fig.savefig("%s/time_to_mortality.pdf" % (odir), bbox_inches='tight',
                    pad_inches=0.1)
    else:
        fig.savefig("%s/time_to_mortality_%d.pdf" % (odir, year),
                    bbox_inches='tight', pad_inches=0.1)
    plt.close('all')


def plot_transpiration(odir, out, year=None):

    conv = c.MMOL_2_MOL * c.MOL_WATER_2_G_WATER * c.G_TO_KG * \
            c.SEC_2_HLFHR

    trans = []
    for i in range(0, len(out), 48):
        vals = out["Eplant"][i:i+48]
        vals = vals[~np.isnan(vals)]
        if len(vals) > 0:
            aet = np.sum(vals * conv)
            trans.append(aet)

    cb = ['#377eb8', '#ff7f00', '#4daf4a', \
          '#f781bf', '#a65628', '#984ea3',\
          '#999999', '#e41a1c', '#dede00']

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    ax1 = fig.add_subplot(111)
    ax1.plot(trans, ls="-", color=cb[1])
    ax1.set_ylabel("Transpiration (mm d$^{-1}$)")
    ax1.set_xlabel("Time (days)")

    if year is None:
        fig.savefig("%s/transpiration.pdf" % (odir), bbox_inches='tight',
                    pad_inches=0.1)
    else:
        fig.savefig("%s/transpiration_%d.pdf" % (odir, year),
                    bbox_inches='tight', pad_inches=0.1)
    plt.close('all')

def plot_transpiration_and_pet(odir, out, year=None):

    conv = c.MMOL_2_MOL * c.MOL_WATER_2_G_WATER * c.G_TO_KG * \
            c.SEC_2_HLFHR

    trans = []
    pet = []
    for i in range(0, len(out), 48):
        vals = out["Eplant"][i:i+48]
        vals = vals[~np.isnan(vals)]
        if len(vals) > 0:
            aet = np.sum(vals * conv)
            trans.append(aet)
        vals = out["pet"][i:i+48]
        vals = vals[~np.isnan(vals)]
        if len(vals) > 0:
            pet.append(np.sum(vals * c.SEC_2_HLFHR))

    cb = ['#377eb8', '#ff7f00', '#4daf4a', \
          '#f781bf', '#a65628', '#984ea3',\
          '#999999', '#e41a1c', '#dede00']

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    ax1 = fig.add_subplot(111)
    ax1.plot(trans, ls="-", color=cb[1], label="AET")
    ax1.plot(pet, ls="-", color=cb[2], label="PET")
    ax1.set_ylabel("(mm d$^{-1}$)")
    ax1.set_xlabel("Time (days)")

    if year is None:
        fig.savefig("%s/transpiration_and_pet.pdf" % (odir), bbox_inches='tight',
                    pad_inches=0.1)
    else:
        fig.savefig("%s/transpiration_and_pet_%d.pdf" % (odir, year),
                    bbox_inches='tight', pad_inches=0.1)
    plt.close('all')

def plot_cwd(odir, out, timestep=15, year=None):

    if timestep == 15:
        ndays = out.t / 96
    elif timestep == 30:
        ndays = out.t / 48.
    elif timestep == 60:
        ndays = out.t / 24.
    """
    cwd = []
    cum_sumx = 0.0
    for i in range(0, len(out), 48):
        pet = np.sum(out.pet[i:i+48] * c.SEC_2_HLFHR)
        aet = np.sum(out["Eplant"][i:i+48] * c.MMOL_2_MOL * \
                     c.MOL_WATER_2_G_WATER * c.G_TO_KG * \
                     c.SEC_2_HLFHR)
        cum_sumx += pet - aet
        cwd.append(cum_sumx)
    """
    cwd = []
    cum_sumx = 0.0

    dx = 0.0
    dy = 0.0
    hod = 0
    for i in range(len(out)):
        pet = out.pet[i] * c.SEC_2_HLFHR
        aet = out["Eplant"][i] * c.MMOL_2_MOL * c.MOL_WATER_2_G_WATER * \
                c.G_TO_KG * c.SEC_2_HLFHR
        cum_sumx += pet - aet
        cwd.append(cum_sumx)
        dx += pet
        dy += aet

        hod += 1
        if hod > 47:
            hod = 0.0
            #print(dx, dy)
            dx = 0.0
            dy = 0.0

    cb = ['#377eb8', '#ff7f00', '#4daf4a', \
          '#f781bf', '#a65628', '#984ea3',\
          '#999999', '#e41a1c', '#dede00']

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    ax1 = fig.add_subplot(111)
    #ax1.set_xlim(48)
    ax1.plot(cwd, ls="-", color=cb[1])

    ax1.set_ylabel("Accumulated CWD (mm)")
    ax1.set_xlabel("Time (days)")

    if year is None:
        fig.savefig("%s/cwd.pdf" % (odir), bbox_inches='tight', pad_inches=0.1)
    else:
        fig.savefig("%s/cwd_%d.pdf" % (odir, year), bbox_inches='tight',
                    pad_inches=0.1)
    plt.close('all')

def plot_sw(odir, out, time_step=30, year=None):

    cb = ['#377eb8', '#ff7f00', '#4daf4a', \
          '#f781bf', '#a65628', '#984ea3',\
          '#999999', '#e41a1c', '#dede00']

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    ax1 = fig.add_subplot(111)
    ax1.plot(out.index, out.sw, ls="-", color=cb[1])
    ax1.set_ylabel("SWC (m$^{3}$ m$^{-3}$)")
    ax1.set_xlabel("Time (days)")

    if year is None:
        fig.savefig("%s/sw.pdf" % (odir), bbox_inches='tight',
                    pad_inches=0.1)
    else:
        fig.savefig("%s/sw_%d.pdf" % (odir, year), bbox_inches='tight',
                    pad_inches=0.1)
    plt.close('all')

def get_dates(out):
    st = pd.to_datetime("%d/%d %d:%d:%d" %
                        (int(out.doy[0])+1, int(out.year[0]),
                         0, 0, 0),
                        format= "%j/%Y %H:%M:%S")
    dates = pd.date_range(st, periods=out.shape[0], freq='30min')

    return dates

if __name__ == "__main__":

    fname = "gswp3_met/GSWP3_met_10_17.csv"
    df = pd.read_csv(fname)
    years = np.unique(df.year)
    fdir = "outputs"
    year = years[0]
    out = pd.read_csv(os.path.join(fdir, "desica_out_%d.csv" % (year)))

    for year in years[1:]:
        df = pd.read_csv(os.path.join(fdir, "desica_out_%d.csv" % (year)))
        out = out.append(df)
    out = out.reset_index(drop=True)
    out["dates"] = get_dates(out)
    out.set_index('dates', inplace=True)

    time_step = 30
    odir = "plots"
    plot_time_to_mortality(odir, out, time_step)
    plot_transpiration(odir, out)
    plot_transpiration_and_pet(odir, out)
    plot_cwd(odir, out, time_step)
    plot_sw(odir, out, time_step)
