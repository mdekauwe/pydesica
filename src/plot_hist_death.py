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
import seaborn as sns


def drop_duplicates(df):

    gmin = df.drop_duplicates(subset='gmin')
    lai = df.drop_duplicates(subset='lai')
    p50 = df.drop_duplicates(subset='p50')
    depth = df.drop_duplicates(subset='depth')
    Cl = df.drop_duplicates(subset='Cl')
    Cs = df.drop_duplicates(subset='Cs')

    return pd.concat([gmin, lai, p50, depth, Cl, Cs])


rf = pd.read_csv("outputs/rf_trait_sensitivity_all.csv")
wsf = pd.read_csv("outputs/wsf_trait_sensitivity_all.csv")
dsf = pd.read_csv("outputs/dsf_trait_sensitivity_all.csv")
grw = pd.read_csv("outputs/grw_trait_sensitivity_all.csv")
saw = pd.read_csv("outputs/saw_trait_sensitivity_all.csv")

#rf = pd.read_csv("outputs/rf_trait_sens_OAT.csv")
#wsf = pd.read_csv("outputs/wsf_trait_sens_OAT.csv")
#dsf = pd.read_csv("outputs/dsf_trait_sens_OAT.csv")
#grw = pd.read_csv("outputs/grw_trait_sens_OAT.csv")
#saw = pd.read_csv("outputs/saw_trait_sens_OAT.csv")

rf = rf[rf.day_of_death > 0]
wsf = wsf[wsf.day_of_death > 0]
dsf = dsf[dsf.day_of_death > 0]
grw = grw[grw.day_of_death > 0]
saw = saw[saw.day_of_death > 0]

# temp until re-run
rf = rf[rf.psi_e != 3.0]
wsf = wsf[wsf.psi_e != 3.0]
dsf = dsf[dsf.psi_e != 3.0]
grw = grw[grw.psi_e != 3.0]
saw = saw[saw.psi_e != 3.0]

#rf = drop_duplicates(rf)
#wsf = drop_duplicates(wsf)
#dsf = drop_duplicates(dsf)
#grw = drop_duplicates(grw)
#saw = drop_duplicates(saw)

#rf = rf.drop_duplicates(subset='day_of_death')
#wsf = wsf.drop_duplicates(subset='day_of_death')
#dsf = dsf.drop_duplicates(subset='day_of_death')
#grw = grw.drop_duplicates(subset='day_of_death')
#saw = saw.drop_duplicates(subset='day_of_death')

width = 9
height = 6
fig = plt.figure(figsize=(width, height))
fig.subplots_adjust(hspace=0.13)
fig.subplots_adjust(wspace=0.13)
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.sans-serif'] = "Helvetica"
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['font.size'] = 16
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16

colours = sns.color_palette("Set2", 8)

ax = fig.add_subplot(111)

sns.distplot(rf.day_of_death, ax=ax, rug=False, norm_hist=True,
             kde_kws={"label": "RF"}, kde=True, color=colours[0])
sns.distplot(wsf.day_of_death, ax=ax, rug=False, norm_hist=True,
             kde_kws={"label": "WSF"}, kde=True, color=colours[1])
sns.distplot(dsf.day_of_death, ax=ax, rug=False, norm_hist=True,
             kde_kws={"label": "DSF"}, kde=True, color=colours[2])


ax.tick_params(direction='in', length=4)
ax.set_xlabel("Day of death")
ax.set_ylabel("Probability density")
ax.legend(numpoints=1, ncol=1, loc="best", frameon=False)

ofdir = "/Users/mdekauwe/Desktop"
ofname = "day_of_death_RF_WSF_DSF.png"
fig.savefig(os.path.join(ofdir, ofname),
            bbox_inches='tight', pad_inches=0.1, dpi=300)


width = 9
height = 6
fig = plt.figure(figsize=(width, height))
fig.subplots_adjust(hspace=0.13)
fig.subplots_adjust(wspace=0.13)
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.sans-serif'] = "Helvetica"
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['font.size'] = 16
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16

colours = sns.color_palette("Set2", 8)

ax = fig.add_subplot(111)

sns.distplot(grw.day_of_death, ax=ax, rug=False, norm_hist=True,
             kde_kws={"label": "GRW"}, kde=True, color=colours[3])
sns.distplot(saw.day_of_death, ax=ax,  rug=False, norm_hist=True,
             kde_kws={"label": "SAW"}, kde=True, color=colours[4])

ax.set_xlim(0, 600)

ax.tick_params(direction='in', length=4)
ax.set_xlabel("Day of death")
ax.set_ylabel("Probability density")
ax.legend(numpoints=1, ncol=1, loc="best", frameon=False)

ofdir = "/Users/mdekauwe/Desktop"
ofname = "day_of_death_GRW_SAW.png"
fig.savefig(os.path.join(ofdir, ofname),
            bbox_inches='tight', pad_inches=0.1, dpi=300)


width = 9
height = 6
fig = plt.figure(figsize=(width, height))
fig.subplots_adjust(hspace=0.13)
fig.subplots_adjust(wspace=0.13)
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.sans-serif'] = "Helvetica"
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['font.size'] = 16
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16

colours = sns.color_palette("Set2", 8)

ax = fig.add_subplot(111)

sns.distplot(rf.day_of_death, ax=ax, rug=False, norm_hist=True,
             kde_kws={"label": "RF"}, kde=True, color=colours[0])
sns.distplot(wsf.day_of_death, ax=ax, rug=False, norm_hist=True,
             kde_kws={"label": "WSF"}, kde=True, color=colours[1])
sns.distplot(dsf.day_of_death, ax=ax, rug=False, norm_hist=True,
             kde_kws={"label": "DSF"}, kde=True, color=colours[2])
sns.distplot(grw.day_of_death, ax=ax, rug=False, norm_hist=True,
             kde_kws={"label": "GRW"}, kde=True, color=colours[3])
sns.distplot(saw.day_of_death, ax=ax,  rug=False, norm_hist=True,
             kde_kws={"label": "SAW"}, kde=True, color=colours[4])

#ax.set_xlim(100, 600)
#ax.set_ylim(0.0, 0.005)

plt.xticks([], [])
plt.setp(ax.get_xticklabels(), visible=False)

ax.tick_params(direction='in', length=4)
ax.set_xlabel("Day of hydraulic failure ($\Psi$$_{\mathrm{crit}}$)", labelpad=10)
ax.set_ylabel("Probability density")
ax.legend(numpoints=1, ncol=1, loc="best", frameon=False)

ofdir = "/Users/mdekauwe/Desktop"
ofname = "day_of_death_all.png"
fig.savefig(os.path.join(ofdir, ofname),
            bbox_inches='tight', pad_inches=0.1, dpi=300)
