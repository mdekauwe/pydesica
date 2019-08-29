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

rf = pd.read_csv("outputs/rf_trait_sens_OAT.csv")
#wsf = pd.read_csv("outputs/wsf_trait_sens_OAT.csv")
#dsf = pd.read_csv("outputs/dsf_trait_sens_OAT.csv")
#grw = pd.read_csv("outputs/grw_trait_sens_OAT.csv")
#saw = pd.read_csv("outputs/saw_trait_sens_OAT.csv")



width = 9
height = 7
fig = plt.figure(figsize=(width, height))
fig.subplots_adjust(hspace=0.4)
fig.subplots_adjust(wspace=0.05)
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.sans-serif'] = "Helvetica"
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['font.size'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14

ax1 = fig.add_subplot(321)
ax2 = fig.add_subplot(322)
ax3 = fig.add_subplot(323)
ax4 = fig.add_subplot(324)
ax5 = fig.add_subplot(325)
ax6 = fig.add_subplot(326)

colours = sns.color_palette("Set2", 8)


tmp = rf.drop_duplicates(subset='gmin')
ax1.plot(tmp.gmin, tmp.plc, color=colours[0], marker="o", ls=" ",
         label="RF", ms=3)
#ax1.plot(tmp.gmin, tmp.day_of_death, "ro")

tmp = rf.drop_duplicates(subset='lai')
ax2.plot(tmp.lai, tmp.plc, color=colours[0], marker="o", ls=" ", ms=3)

tmp = rf.drop_duplicates(subset='p50')
ax3.plot(tmp.p50, tmp.plc, color=colours[0], marker="o", ls=" ", ms=3)

tmp = rf.drop_duplicates(subset='depth')
ax4.plot(tmp.depth, tmp.plc, color=colours[0], marker="o", ls=" ",
         ms=3)

tmp = rf.drop_duplicates(subset='Cl')
ax5.plot(tmp.Cl, tmp.plc, color=colours[0], marker="o", ls=" ", ms=3)

tmp = rf.drop_duplicates(subset='Cs')
ax6.plot(tmp.Cs, tmp.plc, color=colours[0], marker="o", ls=" ", ms=3)

ax1.set_ylim(40, 100)
ax2.set_ylim(40, 100)
ax3.set_ylim(40, 100)
ax4.set_ylim(40, 100)
ax5.set_ylim(40, 100)
ax6.set_ylim(40, 100)

plt.setp(ax2.get_yticklabels(), visible=False)
plt.setp(ax4.get_yticklabels(), visible=False)
plt.setp(ax6.get_yticklabels(), visible=False)

ax3.set_ylabel("PLC (%)")

ax1.set_xlabel("gmin (mmol m$^{-2}$ s$^{-1}$)")
ax2.set_xlabel("LAI (m$^{2}$ m$^{-2}$)")
ax3.set_xlabel("p50 (MPa)")
ax4.set_xlabel("Depth (m)")
ax5.set_xlabel("Cl (mmol MPa$^{-1}$)")
ax6.set_xlabel("Cs (mmol MPa$^{-1}$)")

ax1.legend(numpoints=1, ncol=1, loc="best", frameon=False)

ax1.tick_params(direction='in', length=4)
ax2.tick_params(direction='in', length=4)
ax3.tick_params(direction='in', length=4)
ax4.tick_params(direction='in', length=4)
ax5.tick_params(direction='in', length=4)
ax6.tick_params(direction='in', length=4)

#plt.show()
ofdir = "/Users/mdekauwe/Desktop"
ofname = "plc.pdf"
fig.savefig(os.path.join(ofdir, ofname),
            bbox_inches='tight', pad_inches=0.1)
