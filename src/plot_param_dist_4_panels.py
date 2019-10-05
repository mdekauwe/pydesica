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

df = pd.read_csv("outputs/params.csv")
df.index = df["trait"]

gmin = df.loc["gmin"].values[1:]
p50 = df.loc["p50"].values[1:]
Cs = df.loc["Cs"].values[1:]
kpsat = df.loc["kpsat"].values[1:]

width = 9
height = 6
fig = plt.figure(figsize=(width, height))
fig.subplots_adjust(hspace=0.15)
fig.subplots_adjust(wspace=0.3)
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.sans-serif'] = "Helvetica"
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['font.size'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14


ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)


ax1.scatter(np.arange(len(gmin)), gmin, color="royalblue")
ax2.scatter(np.arange(len(p50)), p50, color="royalblue")
ax3.scatter(np.arange(len(Cs)), Cs, color="royalblue")
ax4.scatter(np.arange(len(kpsat)), kpsat, color="royalblue")


plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)

labels = ["", "RF", "WSF", "DSF", "GRW", "SAW"]
ax3.set_xticklabels(labels)
ax4.set_xticklabels(labels)


ax1.set_ylabel("$g_{\mathrm{min}}$\n(mmol m$^{-2}$ s$^{-1}$)")
ax2.set_ylabel("$p_{\mathrm{50}}$\n(MPa)")
ax3.set_ylabel("$C_{\mathrm{s}}$\n(mmol kg$^{-1}$ MPa$^{-1}$)")
ax4.set_ylabel("$K_{\mathrm{psat}}$\n(mmol m$^{-2}$ s$^{-1}$ MPa$^{-1}$)")

ofdir = "/Users/mdekauwe/Desktop"
ofname = "day_of_death_rel_imp.png"
fig.savefig(os.path.join(ofdir, ofname),
            bbox_inches='tight', pad_inches=0.1, dpi=300)
