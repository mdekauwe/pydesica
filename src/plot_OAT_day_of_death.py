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
from get_params import get_params


def make_plot(df, pft_name):


    width = 9
    height = 7
    fig = plt.figure(figsize=(width, height))
    fig.subplots_adjust(hspace=0.5)
    fig.subplots_adjust(wspace=0.2)
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


    tmp = df.drop_duplicates(subset='gmin')
    ax1.plot(tmp.gmin, tmp.day_of_death, color=colours[0], marker="o", ls=" ",
             ms=3)

    tmp = df.drop_duplicates(subset='lai')
    ax2.plot(tmp.lai, tmp.day_of_death, color=colours[0], marker="o", ls=" ", ms=3)

    tmp = df.drop_duplicates(subset='p50')
    ax3.plot(tmp.p50, tmp.day_of_death, color=colours[0], marker="o", ls=" ", ms=3)

    tmp = df.drop_duplicates(subset='depth')
    ax4.plot(tmp.depth, tmp.day_of_death, color=colours[0], marker="o", ls=" ",
             ms=3)

    tmp = df.drop_duplicates(subset='Cl')
    ax5.plot(tmp.Cl, tmp.day_of_death, color=colours[0], marker="o", ls=" ", ms=3)

    tmp = df.drop_duplicates(subset='Cs')
    ax6.plot(tmp.Cs, tmp.day_of_death, color=colours[0], marker="o", ls=" ", ms=3)


    #ax1.set_ylim(40, 70)
    #ax2.set_ylim(40, 70)
    #ax3.set_ylim(40, 70)
    #ax4.set_ylim(40, 70)
    #ax5.set_ylim(40, 70)
    #ax6.set_ylim(40, 70)

    #plt.setp(ax2.get_yticklabels(), visible=False)
    #plt.setp(ax4.get_yticklabels(), visible=False)
    #plt.setp(ax6.get_yticklabels(), visible=False)

    ax3.set_ylabel("Day of death")

    ax1.set_xlabel("gmin (mmol m$^{-2}$ s$^{-1}$)")
    ax2.set_xlabel("LAI (m$^{2}$ m$^{-2}$)")
    ax3.set_xlabel("p50 (MPa)")
    ax4.set_xlabel("Depth (m)")
    ax5.set_xlabel("Cl (mmol MPa$^{-1}$)")
    ax6.set_xlabel("Cs (mmol MPa$^{-1}$)")

    #ax1.legend(numpoints=1, ncol=1, loc="best", frameon=False)

    ax1.tick_params(direction='in', length=4)
    ax2.tick_params(direction='in', length=4)
    ax3.tick_params(direction='in', length=4)
    ax4.tick_params(direction='in', length=4)
    ax5.tick_params(direction='in', length=4)
    ax6.tick_params(direction='in', length=4)

    #plt.show()
    ofdir = "/Users/mdekauwe/Desktop"
    ofname = "day_of_death_%s.pdf" % (pft_name)
    fig.savefig(os.path.join(ofdir, ofname),
                bbox_inches='tight', pad_inches=0.1)




if __name__ == "__main__":

    params = get_params()
    pfts = list(params)

    for pft_name in pfts:

        p = params[pft_name]

        df = pd.read_csv("outputs/%s_trait_sens_OAT.csv" % (pft_name))
        df = df[df.day_of_death > 0]
        make_plot(df, pft_name)
