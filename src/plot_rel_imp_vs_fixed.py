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

def get_rel_imp(chg, fix_eff):
    mu = np.mean(np.abs(chg - fix_eff) / fix_eff * 100.0)
    sigma = np.std(np.abs(chg - fix_eff) / fix_eff * 100.0)

    return mu, sigma

df_fix = pd.read_csv("outputs/trait_sens_fixed.csv")


df = pd.DataFrame(columns=["pft","gmin","lai","p50","depth","Cl","Cs"])

for pft_name in ["rf", "wsf", "dsf", "grw", "saw"]:

    dfx = pd.read_csv("outputs/%s_trait_sens_OAT.csv" % (pft_name))
    dfx = dfx[dfx.day_of_death > 0]
    gmin = dfx.drop_duplicates(subset='gmin')

    fix_eff = df_fix[df_fix.pft == pft_name].day_of_death.values
    chg = gmin.day_of_death.values
    gmin_mu, gmin_sigma = get_rel_imp(chg, fix_eff)

    dfx = pd.read_csv("outputs/%s_trait_sens_OAT.csv" % (pft_name))
    dfx = dfx[dfx.day_of_death > 0]
    lai = dfx.drop_duplicates(subset='lai')

    fix_eff = df_fix[df_fix.pft == pft_name].day_of_death.values
    chg = lai.day_of_death.values
    lai_mu, lai_sigma = get_rel_imp(chg, fix_eff)
    print(fix_eff)
    print(chg)
    print(lai_mu)
    print("\n")
    lai_mu, lai_sigma = get_rel_imp(chg, fix_eff)

    dfx = pd.read_csv("outputs/%s_trait_sens_OAT.csv" % (pft_name))
    dfx = dfx[dfx.day_of_death > 0]
    p50 = dfx.drop_duplicates(subset='p50')

    fix_eff = df_fix[df_fix.pft == pft_name].day_of_death.values
    chg = p50.day_of_death.values
    p50_mu, p50_sigma = get_rel_imp(chg, fix_eff)

    dfx = pd.read_csv("outputs/%s_trait_sens_OAT.csv" % (pft_name))
    dfx = dfx[dfx.day_of_death > 0]
    depth = dfx.drop_duplicates(subset='depth')

    fix_eff = df_fix[df_fix.pft == pft_name].day_of_death.values
    chg = depth.day_of_death.values
    depth_mu, depth_sigma = get_rel_imp(chg, fix_eff)

    dfx = pd.read_csv("outputs/%s_trait_sens_OAT.csv" % (pft_name))
    dfx = dfx[dfx.day_of_death > 0]
    Cl = dfx.drop_duplicates(subset='Cl')

    fix_eff = df_fix[df_fix.pft == pft_name].day_of_death.values
    chg = Cl.day_of_death.values
    Cl_mu, Cl_sigma = get_rel_imp(chg, fix_eff)

    dfx = pd.read_csv("outputs/%s_trait_sens_OAT.csv" % (pft_name))
    dfx = dfx[dfx.day_of_death > 0]
    Cs = dfx.drop_duplicates(subset='Cs')

    fix_eff = df_fix[df_fix.pft == pft_name].day_of_death.values
    chg = Cs.day_of_death.values
    Cs_mu, Cs_sigma = get_rel_imp(chg, fix_eff)

    result = [pft_name, gmin_mu, lai_mu, p50_mu, depth_mu, Cl_mu, Cs_mu]
    s = pd.Series(result, index=df.columns)
    df = df.append(s, ignore_index=True)


df.pft = df.pft.str.upper()


width = 9
height = 6
fig = plt.figure(figsize=(width, height))
fig.subplots_adjust(hspace=0.15)
fig.subplots_adjust(wspace=0.15)
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

sns.barplot(x="pft", y="gmin", data=df, palette="muted", ax=ax1)
sns.barplot(x="pft", y="lai",  data=df, palette="muted", ax=ax2)
sns.barplot(x="pft", y="p50", data=df, palette="muted", ax=ax3)
sns.barplot(x="pft", y="depth", data=df, palette="muted", ax=ax4)
sns.barplot(x="pft", y="Cl", data=df, palette="muted", ax=ax5)
sns.barplot(x="pft", y="Cs", data=df, palette="muted", ax=ax6)
#ax1.legend_.remove()
#ax2.legend_.remove()
#ax3.legend_.remove()
#ax4.legend_.remove()
#ax5.legend_.remove()
#ax6.legend_.remove()

def change_width(ax, new_value) :
    for patch in ax.patches :
        current_width = patch.get_width()
        diff = current_width - new_value

        # we change the bar width
        patch.set_width(new_value)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)

change_width(ax1, .9)
change_width(ax2, .9)
change_width(ax3, .9)
change_width(ax4, .9)
change_width(ax5, .9)
change_width(ax6, .9)


ax1.set_ylim(0,80)
ax2.set_ylim(0,80)
ax3.set_ylim(0,80)
ax4.set_ylim(0,80)
ax5.set_ylim(0,80)
ax6.set_ylim(0,80)

ax1.set_xlim(-1,5)
ax2.set_xlim(-1,5)
ax3.set_xlim(-1,5)
ax4.set_xlim(-1,5)
ax5.set_xlim(-1,5)
ax6.set_xlim(-1,5)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
plt.setp(ax4.get_xticklabels(), visible=False)

#plt.setp(ax2.get_yticklabels(), visible=False)
#plt.setp(ax4.get_yticklabels(), visible=False)
#plt.setp(ax6.get_yticklabels(), visible=False)

ax1.set_xlabel('')
ax2.set_xlabel('')
ax3.set_xlabel('')
ax4.set_xlabel('')
ax5.set_xlabel('')
ax6.set_xlabel('')

ax1.set_ylabel('')
ax2.set_ylabel('')
ax3.set_ylabel('Percentage change in\nday of death', position=(1.5, 0.5))
ax4.set_ylabel('')
ax5.set_ylabel('')
ax6.set_ylabel('')


props = dict(boxstyle='round', facecolor='white', alpha=1.0,
                     ec="white")
ax1.text(0.03, 0.93, "(a) $g_{\mathrm{min}}$",
        transform=ax1.transAxes, fontsize=14, verticalalignment='top',
        bbox=props)
ax2.text(0.03, 0.93, "(b) LAI",
        transform=ax2.transAxes, fontsize=14, verticalalignment='top',
        bbox=props)
ax3.text(0.03, 0.93, "(c) $p_{\mathrm{50}}$",
        transform=ax3.transAxes, fontsize=14, verticalalignment='top',
        bbox=props)
ax4.text(0.03, 0.93, "(d) Depth",
        transform=ax4.transAxes, fontsize=14, verticalalignment='top',
        bbox=props)
ax5.text(0.03, 0.93, "(e) $C_{\mathrm{l}}$",
        transform=ax5.transAxes, fontsize=14, verticalalignment='top',
        bbox=props)
ax6.text(0.03, 0.93, "(f) $C_{\mathrm{s}}$",
        transform=ax6.transAxes, fontsize=14, verticalalignment='top',
        bbox=props)

ofdir = "/Users/mdekauwe/Desktop"
ofname = "day_of_death_rel_imp.pdf"
fig.savefig(os.path.join(ofdir, ofname),
            bbox_inches='tight', pad_inches=0.1)
