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


rf = pd.read_csv("outputs/rf_trait_sens_OAT.csv")
#wsf = pd.read_csv("outputs/wsf_trait_sensitivity.csv")
#dsf = pd.read_csv("outputs/dsf_trait_sensitivity.csv")
#grw = pd.read_csv("outputs/grw_trait_sensitivity.csv")
##saw = pd.read_csv("outputs/saw_trait_sensitivity.csv")

rf = rf[rf.min_plc > 0]
#rf = rf.drop_duplicates(subset='min_plc')
#print(len(rf))
#wsf = wsf[wsf.plc > 0]
#dsf = dsf[dsf.plc > 0]
#grw = grw[grw.plc > 0]
#saw = saw[saw.plc > 0]


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

ax1 = fig.add_subplot(321)
ax2 = fig.add_subplot(322)
ax3 = fig.add_subplot(323)
ax4 = fig.add_subplot(324)
ax5 = fig.add_subplot(325)
ax6 = fig.add_subplot(326)


tmp = rf.drop_duplicates(subset='gmin')
print(np.unique(tmp.gmin))
ax1.plot(tmp.gmin, tmp.day_of_death, "ro")

tmp = rf.drop_duplicates(subset='lai')
ax2.plot(tmp.lai, tmp.day_of_death, "ro")

tmp = rf.drop_duplicates(subset='p50')
ax3.plot(tmp.p50, tmp.day_of_death, "ro")

tmp = rf.drop_duplicates(subset='depth')
ax4.plot(tmp.depth, tmp.day_of_death, "ro")

tmp = rf.drop_duplicates(subset='Cl')
ax5.plot(tmp.Cl, tmp.day_of_death, "ro")

tmp = rf.drop_duplicates(subset='Cs')
ax6.plot(tmp.Cs, tmp.day_of_death, "ro")

ax1.set_ylim(40, 70)
ax2.set_ylim(40, 70)
ax3.set_ylim(40, 70)
ax4.set_ylim(40, 70)
ax5.set_ylim(40, 70)
ax6.set_ylim(40, 70)

#ax.legend(numpoints=1, ncol=1, loc="best", frameon=False)
plt.show()
ofdir = "/Users/mdekauwe/Desktop"
ofname = "plc.pdf"
fig.savefig(os.path.join(ofdir, ofname),
            bbox_inches='tight', pad_inches=0.1)
