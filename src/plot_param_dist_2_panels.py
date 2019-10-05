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

width = 6
height = 9
fig = plt.figure(figsize=(width, height))
fig.subplots_adjust(hspace=0.15)
fig.subplots_adjust(wspace=0.3)
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.sans-serif'] = "Helvetica"
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['font.size'] = 16
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16


ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)


ax1.plot(gmin, ls=" ", marker="o", color="royalblue", ms=15)
ax2.plot(p50, ls=" ", marker="o", color="royalblue", ms=15)


plt.setp(ax1.get_xticklabels(), visible=False)

labels = ["RF", "WSF", "DSF", "GRW", "SAW"]
ax1.set_xticks(np.arange(5))
ax2.set_xticks(np.arange(5))
ax2.set_xticklabels(labels)


ax1.set_ylabel("$g_{\mathrm{min}}$\n(mmol m$^{-2}$ s$^{-1}$)")
ax2.set_ylabel("$p_{\mathrm{50}}$\n(MPa)")

ofdir = "/Users/mdekauwe/Desktop"
ofname = "day_of_death_rel_imp.png"
fig.savefig(os.path.join(ofdir, ofname),
            bbox_inches='tight', pad_inches=0.1, dpi=300)
