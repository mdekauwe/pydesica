#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

def fsig_tuzet(psi_leaf, sf, psi_f):
    """
    An empirical logistic function to describe the sensitivity of stomata
    to leaf water potential.

    Sigmoid function assumes that stomata are insensitive to psi_leaf at
    values close to zero and that stomata rapidly close with decreasing
    psi_leaf.

    Parameters:
    -----------
    psi_leaf : float
        leaf water potential (MPa)

    Returns:
    -------
    fw : float
        sensitivity of stomata to leaf water potential [0-1]

    Reference:
    ----------
    * Tuzet et al. (2003) A coupled model of stomatal conductance,
      photosynthesis and transpiration. Plant, Cell and Environment 26,
      10971116

    """
    num = 1.0 + np.exp(sf * psi_f)
    den = 1.0 + np.exp(sf * (psi_f - psi_leaf))
    fw = num / den

    return fw


if __name__ == "__main__":

    psi_leaf = np.linspace(-0.05, -6, 50)

    psi_f = -2.0
    sf = 8.0
    fw_d = fsig_tuzet(psi_leaf, sf, psi_f)

    psi_f = -2.455474
    sf = 2.000000
    fw19 = fsig_tuzet(psi_leaf, sf, psi_f)

    psi_f = -1.668544
    sf = 2.000000
    fw20 = fsig_tuzet(psi_leaf, sf, psi_f)

    plt.plot(psi_leaf, fw_d, label="desica_test")
    plt.plot(psi_leaf, fw19, label="wsf")
    plt.plot(psi_leaf, fw20, label="dsf")
    plt.legend(numpoints=1, loc="best")
    plt.show()
