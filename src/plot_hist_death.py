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

dsf = pd.read_csv("outputs/dsf_trait_sensitivity.csv")
saw = pd.read_csv("outputs/saw_trait_sensitivity.csv")

dsf = dsf[dsf.day_of_death > 0]
saw = saw[saw.day_of_death > 0]

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

ax = fig.add_subplot(111)
#ax.hist(dsf.day_of_death, alpha=0.5, label='DSF', density=True)
#ax.hist(saw.day_of_death, alpha=0.5, label='SAW', density=True)

sns.distplot(dsf.day_of_death, ax=ax, rug=False, norm_hist=True,
             kde_kws={"label": "DSF"})
sns.distplot(saw.day_of_death, ax=ax, rug=False, norm_hist=True,
             kde_kws={"label": "SAW"})
ax.tick_params(direction='in', length=4)
ax.set_xlabel("Day of death")
ax.set_ylabel("Probability density")
ax.legend(numpoints=1, ncol=1, loc="best", frameon=False)
plt.show()
