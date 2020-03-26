#!/usr/bin/env python
"""
Summarise the params into PFTs
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (08.08.2019)"
__email__ = "mdekauwe@gmail.com"

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from tabulate import tabulate

def get_params():
    fname = "/Users/mdekauwe/Dropbox/ARC drought Trial/Parameters for Martin/parameters/parameters.csv"
    df = pd.read_csv(fname)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

    rf = df[df.spp == "Asm"].mean(axis=0)
    wsf = df[(df.spp == "Egr") | (df.spp == "Evi")].mean(axis=0)
    dsf = df[(df.spp == "Aco") | (df.spp == "Cgu") | (df.spp == "Esi")].mean(axis=0)
    grw = df[(df.spp == "Ebl") | (df.spp == "Ema") | (df.spp == "Eme")].mean(axis=0)
    saw = df[(df.spp == "Aan") | (df.spp == "Ela") | (df.spp == "Epo")].mean(axis=0)

    sites = pd.concat([rf, wsf, dsf, grw, saw], axis=1)
    sites.columns = ["rf","wsf","dsf","grw","saw"]

    # Till its fixed, gap fill the rf Kplant
    sites["rf"]["Kplant"] = df.Kplant.mean()
    #sites = sites.drop(["Cs", "Cl"])
    sites = sites.transpose()


    #sites.Vcmax = sites.Vcmax/1e6
    #sites['ejmax'] = sites.Vcmax * 1.67
    sites['Jmax'] = sites.Vcmax * 1.67


    #sites = sites.rename(columns={'Cbranch_mmol.kg.MPa':'Cs',
    #                              'Cleaf_preTLP_mmol.m2.MPa1':'Cl'})

    sites = sites.rename(columns={'Cwood_mmol.kg.MPa':'Cs',
                                  'Cleaf_preTLP_mmol.m2.MPa1':'Cl'})

    sites = sites.transpose()

    sites.index.name = "trait"
    sites.to_csv("outputs/params.csv")

    return (sites)

if __name__ == "__main__":


    sites = get_params()
    #print(sites)

    sites = sites.round(1)
    sites = sites.drop(['AL', 'SM', 'SLA', 'LM' ,'kpsat'], axis=0)


    sites.loc['gmin',:] /= 2.0 # original was two sided

    sites = sites.reindex(["g1","gmin","Vcmax", "Jmax","psiv","sf",\
                            "Kplant","s50","p50","Cl","Cs"])

    sites.columns = sites.columns.str.upper()

    sites["Definitions"] = ["Stomatal slope", \
                            "Cuticular conductance",\
                            "Value of Vcmax at 25 °C", \
                            "Value of Jmax at 25 °C",\
                            "Reference water potential", \
                            "Shape of response to $\Psi_{l}$",\
                            "Plant hydraulic conductance", \
                            "Stomatal sensitivity parameter",\
                            "$\Psi$ at 50% loss of hydraulic conductivity", \
                            "Leaf capacitance",\
                            "Stem capacitance"]
    sites["Units"] = ["-", \
                      "mmol m^-2^ s^-1^",\
                      "$\mu$mol m^-2^ s^-1^", \
                      "$\mu$mol m^-2^ s^-1^",\
                      "MPa", \
                      "MPa^−1^",\
                      "mmol m^-2^ leaf s^-1^ MPa^-1^", \
                      "% MPa^-1^",\
                      "MPa", \
                      "mmol m^-2^ s^-1^ MPa^-1^",\
                      "mmol m^-2^ s^-1^ MPa^-1^"]
    cols = sites.columns.tolist()
    cols = cols[-2:] + cols[:-2]
    sites = sites[cols]

    sites.index = ["g~1~", "g~min~", "V~cmax~ ", "J~max~", \
                   "$\Psi$~f~", "S~f~", "k~plant~", "S~50~", \
                   "P~50~", "C~l~", "C~s~"]


    print(tabulate(sites, tablefmt="pipe", headers="keys"))
