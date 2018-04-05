#!/usr/bin/env python

"""
Desica model: simple plant hydraulics model with mortality.

Key features:
- plant hydraulic conductance depends on xylem water potential, defined by a
  PLC curve
- Stomatal conductance is modelled via the modified Tuzet model. Use of Tuzet
  allows parameterisations of plants with differing stomatal control over
  psi_leaf, i.e. isohydric vs anisohydric.
- During severe drought, root water uptake ceases (due to the decline in
  soil-to-root conductance), but water loss continues due to the minimum
  conductance (gmin). This leads to a gradual decline in the stem water pool.
- The water storage pool is split between the stem and leaf water storage pools
- Xylem water potential is calculated from stem water storage via a simple
  constant capacitance term (and likewise, for the leaf water pool).
- Model solves two fluxes: flux of water from the stem to the leaves and the
  flux of water from the soil to the stem numerically in every timestep. To
  simply this process we use psi_leaf and psi_stem from the previous timestep
  to find the flux_to_leaf and the flux_to_stem
- Leaves are assumed to be perfectly coupled, so transpiration = Eleaf * Al,
  where Al is plant area (m2).
- Approach follows Xu to solve the psi_leaf, psi_stem without the need for a
  numerical integrator.  This avoids potential numerical instabilities due to
  various dependancies on water potential. This method works well at short
  timesteps (up to about) 10 to 15 mins. NB we may still get some ocillations.

References:
----------
* Duursma & Choat (2017). Fitplc - an R package to fit hydraulic vulnerability
  curves. Journal of Plant Hydraulics, 4, e002.
* Tuzet et al. (2003) A coupled model of stomatal conductance,
  photosynthesis and transpiration. Plant, Cell and Environment 26,
  1097–1116.
* Xu X, Medvigy D, Powers JS, Becknell JM, Guan K (2016) Diversity in plant
  hydraulic traits explains seasonal and inter-annual variations of vegetation
  dynamics in seasonally dry tropical forests. New Phytologist, 212, 80–95.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (05.03.2018)"
__email__ = "mdekauwe@gmail.com"

import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
from generate_met_data import generate_met_data
from canopy import Canopy, FarquharC3
from math import isclose
from calc_pet import calc_net_radiation, calc_pet_energy
import constants as c

class Desica(object):

    def __init__(self, plc_dead=88., soil_depth=1.0, ground_area=1.0,
                 met_timestep=30., sf=8., g1=4., Cs=100000., b=6.,
                 Cl=10000., kp_sat=3., p50=-4., psi_f=-2., s50=30., gmin=10,
                 psi_leaf0=-1., psi_stem0=-0.5, theta_sat=0.5, sw0=0.5, AL=2.5,
                 psi_e=-0.8*1E-03, Ksat=20., Lv=10000., F=None, keep_wet=False,
                 stop_dead=True, nruns=1, rroot=1E-06):

        self.keep_wet = keep_wet
        self.stop_dead = stop_dead
        self.plc_dead = plc_dead
        self.nruns = nruns
        self.soil_depth = soil_depth # depth of soil bucket, m
        self.ground_area = ground_area # m
        self.soil_volume = self.ground_area * self.soil_depth # m3
        self.met_timestep = met_timestep
        self.sf = sf # sensitivity parameter, MPa-1
        self.g1 = g1 # sensitivity of stomatal conductance to the assimilation
                     # rate, kPa
        self.Cs = Cs # stem capacitance, mmol MPa-1
        self.Cl = Cl # leaf capacitance, mmol MPa-1 (total plant)
        self.kp_sat = kp_sat # plant saturated hydraulic conductance
                             # (mmol m-2 s-1 MPa-1)
        self.p50 = p50 # xylem pressure inducing 50% loss of hydraulic
                       # conductivity due to embolism, MPa
        self.psi_f = psi_f # reference potential for Tuzet model, MPa
        self.s50 = s50 # is slope of the curve at P50 used in weibull model,
                       # % MPa-1
        self.gmin = gmin # minimum stomatal conductance, mmol m-2 s-1
        self.psi_leaf0 = psi_leaf0 # initial leaf water potential, MPa
        self.psi_stem0 = psi_stem0 # initial stem water potential, MPa
        self.theta_sat = theta_sat # soil water capacity at saturation (m3 m-3)
        self.sw0 = sw0 # initial soil volumetric water content (m3 m-3)
        self.AL = AL # plant leaf area, m2
        self.lai = AL / self.ground_area # leaf area index, m2 m-2
        self.b = b # empirical coefficient related to the clay content of the
                   # soil (Cosby et al. 1984).
        self.psi_e = psi_e # air entry point water potential (MPa)
        self.Ksat = Ksat # saturated conductivity, mol m-1 s-1 MPa-1
        self.Lv = Lv # root length density, m m-3
        self.F = F
        self.rroot = rroot # mean radius of water absorbing roots, m
        self.timestep_sec = 60. * self.met_timestep / self.nruns


    def run_simulation(self, met=None):
        """
        Main wrapper to control everything

        Parameters:
        -----------
        met : object
            met forcing variables: day; Ca; par; precip; press; tair; vpd

        Returns:
        -------
        out : object
            output dataframe containing calculations for each timestep

        """
        (n, out) = self.initialise_model(met)

        for i in range(1, n):

            out = self.run_timestep(i, met, out)

            # save solutions, use as input for another run,
            # keeping everything else the same, this is so we can solve
            # psi_leaf and psi_stem without the need for a numerical integrator.
            # This approach works well for short timesteps (10-15 mins)
            for j in range(1, self.nruns):
                out.psi_leaf[i-1] = out.psi_leaf[i]
                out.psi_stem[i-1] = out.psi_stem[i]
                out = self.run_timestep(i, met, out)

                sw_rad = met.par[i] * c.PAR_2_SW
                rnet = calc_net_radiation(sw_rad, met.tair[i], albedo=0.15)
                # W m-2 -> MJ m-2 s-1
                rnet *= c.J_TO_MJ
                out.pet[i] = calc_pet_energy(rnet)

            # Stop the simulation if we've died, i.e. reached P88
            if self.stop_dead:
                plc = self.calc_plc(out.kplant[i])
                if plc > self.plc_dead:
                    break

        out["plc"] = self.calc_plc(out.kplant)

        # mmol s-1
        out["Eplant"] = self.AL * out.Eleaf

        out["t"] = np.arange(1, n+1)

        return (out)

    def initialise_model(self, met):
        """
        Set everything up: set initial values, build an output dataframe to
        save things

        Parameters:
        -----------
        met : object
            met forcing variables: day; Ca; par; precip; press; tair; vpd

        Returns:
        -------
        n : int
            number of timesteps in the met file
        out : object
            output dataframe to store things as we go along

        """
        n = len(met)

        out = self.setup_out_df(met)
        out.psi_leaf[0] = self.psi_leaf0
        out.psi_stem[0] = self.psi_stem0
        out.sw[0] = self.sw0
        out.psi_soil[0] = self.calc_swp(self.sw0)
        out.Eleaf[0] = 0.0
        out.pet[0] = 0.0

        # soil hydraulic conductance (mmol m-2 s-1 MPa-1)
        out.ksoil[0] = self.calc_ksoil(out.psi_soil[0])

        return n, out

    def setup_out_df(self, met):
        """
        Create and output dataframe to save things

        Parameters:
        -----------
        met : object
            met forcing variables: day; Ca; par; precip; press; tair; vpd

        Returns:
        -------
        out : object
            output dataframe to store things as we go along.
        """
        dummy = np.ones(len(met)) * np.nan
        out = pd.DataFrame({'Eleaf':dummy,
                            'psi_leaf':dummy,
                            'psi_stem':dummy,
                            'psi_soil':dummy,
                            'sw':dummy,
                            'ksoil':dummy,
                            'kplant':dummy,
                            'flux_to_leaf':dummy,
                            'flux_to_stem':dummy,
                            'ksoil2stem':dummy,
                            'kstem2leaf':dummy,
                            'pet':dummy,
                            'cwd':dummy})

        return out

    def run_timestep(self, i, met, out):

        self.calc_conductances(out, i)

        # modified Tuzet model of stomatal conductance
        mult = (self.g1 / met.Ca[i]) * self.fsig_tuzet(out.psi_leaf[i-1])

        # Calculate photosynthesis and stomatal conductance
        gsw = self.F.canopy(met.Ca[i], met.tair[i], met.par[i],
                            met.vpd[i], mult)

        # Don't add gmin, instead use it as the lower boundary
        gsw = max(self.gmin, c.mol_2_mmol * gsw)

        # Leaf transpiration assuming perfect coupling, mmol m-2 s-1
        out.Eleaf[i] = gsw * (met.vpd[i] / met.press[i])

        out.psi_leaf[i] = self.calc_lwp(out.kstem2leaf[i], out.psi_stem[i-1],
                                        out.psi_leaf[i-1], out.Eleaf[i])

        # Flux from stem to leaf (mmol s-1) = change in leaf storage,
        # plus transpiration
        out.flux_to_leaf[i] = self.calc_flux_to_leaf(out.psi_leaf[i],
                                                     out.psi_leaf[i-1],
                                                     out.Eleaf[i])

        # Update stem water potential
        out.psi_stem[i] = self.update_stem_wp(out.ksoil2stem[i],
                                              out.psi_soil[i-1],
                                              out.flux_to_leaf[i],
                                              out.psi_stem[i-1])

        # Flux from the root to the stem
        out.flux_to_stem[i] = self.calc_flux_to_stem(out.psi_stem[i],
                                                     out.psi_stem[i-1],
                                                     out.flux_to_leaf[i])

        out.sw[i] = self.update_sw_bucket(met.precip[i], out.flux_to_stem[i],
                                          out.sw[i-1])

        # for debugging / comparison
        if self.keep_wet:
            out.sw[i] = out.sw[i-1]

        # Update soil water potential
        out.psi_soil[i] = self.calc_swp(out.sw[i])

        # Update soil-to-root hydraulic conductance (mmol m-2 s-1 MPa-1)
        out.ksoil[i] = self.calc_ksoil(out.psi_soil[i])

        return out

    def calc_conductances(self, out, i):
        """
        Update the simple bucket soil water balance

        Parameters:
        -----------
        out : object
            output dataframe to access previous calculations and update new
            states
        i : int
            current index
        """

        # Plant hydraulic conductance (mmol m-2 s-1 MPa-1). NB. depends on stem
        # water potential from the previous timestep.
        out.kplant[i] = self.kp_sat * self.fsig_hydr(out.psi_stem[i-1])

        # Conductance from root surface to the stem water pool (assumed to be
        # halfway to the leaves)
        kroot2stem = 2.0 * out.kplant[i]

        # Conductance from soil to stem water store (mmol m-2 s-1 MPa-1)
        # (conductances combined in series)
        out.ksoil2stem[i] = 1.0 / (1.0 / out.ksoil[i-1] + 1.0 / kroot2stem)

        # Conductance from stem water store to the leaves (mmol m-2 s-1 MPa-1)
        # assumning the water pool is halfway up the stem
        out.kstem2leaf[i] = 2.0 * out.kplant[i]

    def fsig_hydr(self, psi_stem_prev):
        """
        Calculate the relative conductance as a function of xylem pressure
        using the Weibull (sigmoidal) model based on values of P50
        and S50 which are obtained by fitting curves to measured data.

        Higher values for s50 indicate a steeper response to xylem pressure.

        Parameters:
        -----------
        psi_stem_prev : object
            stem water potential from previous timestep, MPa

        Returns:
        --------
        relk : float
            relative conductance (K/Kmax) as a funcion of xylem pressure (-)

        References:
        -----------
        * Duursma & Choat (2017). Journal of Plant Hydraulics, 4, e002.
        """
        # xylem pressure
        P = np.abs(psi_stem_prev)

        # the xylem pressure (P) x% of the conductivity is lost
        PX = np.abs(self.p50)
        V = (50.0 - 100.) * np.log(1.0 - 50. / 100.)
        p = (P / PX)**((PX * self.s50) / V)

        # relative conductance (K/Kmax) as a funcion of xylem pressure
        relk = (1. - 50. / 100.)**p

        return (relk)

    def calc_lwp(self, kstem2leaf, psi_stem_prev, psi_leaf_prev, Eleaf):
        """
        Calculate leaf water potential, MPa

        This is a simplified equation based on Xu et al., using the water
        potentials from the previous timestep and the fact that we increase
        the temporal resolution to get around the need to solve the dynamic eqn
        with a numerical approach, i.e., Runge-Kutta.

        NB. kstem2leaf which is the stem conductance in CABLE needs to be a
            function of sapwood, tree height.


        Parameters:
        -----------
        kstem2leaf : float
            conductance from stem to leaf, mmol m-2 s-1 MPa-1
        psi_stem_prev : float
            stem water potential from the previous timestep, MPa
        psi_leaf_prev : float
            leaf water potential from the previous timestep, MPa
        Eleaf : float
            transpiration, mmol m-2 s-1

        Returns:
        --------
        psi_leaf : float
            leaf water potential, MPa

        References:
        -----------
        * Xu et al. (2016) New Phytol, 212: 80–95. doi:10.1111/nph.14009; see
          appendix and code. Can write the dynamic equation as:
          dpsi_leaf_dt = b + a*psi_leaf
        """

        ap = -(self.AL * kstem2leaf / self.Cl)
        bp = (self.AL * kstem2leaf * psi_stem_prev - self.AL * Eleaf) / self.Cl
        psi_leaf = ((ap * psi_leaf_prev + bp) * \
                    np.exp(ap * self.timestep_sec) - bp) / ap

        return psi_leaf

    def calc_swp(self, sw):
        """
        Calculate the soil water potential (MPa). The params The parameters b
        and psi_e are estimated from a typical soil moisture release function.

        Parameters:
        -----------
        sw : object
            volumetric soil water content, m3 m-3

        Returns:
        -----------
        psi_soil : float
            soil water potential, MPa

        References:
        -----------
        * Duursma et al. (2008) Tree Physiology 28, 265–276, eqn 10
        """
        return self.psi_e * (sw / self.theta_sat)**-self.b

    def calc_flux_to_stem(self, psi_stem, psi_stem_prev, flux_to_leaf):
        """
        Calculate the flux from the root to the stem, i.e. the root water
        uptake (mmol s-1) = change in stem storage plus flux_to_leaf

        Parameters:
        -----------
        psi_stem : float
            stem water potential, MPa
        psi_stem_prev : float
            stem water potential from the previous timestep, MPa
        flux_to_leaf : float
            flux of water from the stem to leaf

        Returns:
        -------
        flux_to_stem : float
            flux from soil to the stem, mmol s-1
        """
        return (psi_stem - psi_stem_prev) * \
                self.Cs / self.timestep_sec + flux_to_leaf

    def calc_flux_to_leaf(self, psi_leaf, psi_leaf_prev, Eleaf):
        """
        Calculate the flux from the stem to the leaf = change in leaf storage
        plus transpiration

        Parameters:
        -----------
        psi_leaf : float
            leaf water potential, MPa
        psi_leaf_prev : float
            leaf water potential from the previous timestep, MPa
        Eleaf : float
            transpiration, mmol m-2 s-1

        Returns:
        -------
        flux_to_leaf : float
            flux from stem to the leaf, mmol s-1
        """
        return (psi_leaf - psi_leaf_prev) * \
                self.Cl / self.timestep_sec + self.AL * Eleaf

    def update_stem_wp(self, ksoil2stem, psi_soil_prev, flux_to_leaf,
                       psi_stem_prev):
        """
        Calculate the flux from the stem to the leaf = change in leaf storage
        plus transpiration

        This is a simplified equation based on Xu et al., using the water
        potentials from the previous timestep and the fact that we increase
        the temporal resolution to get around the need to solve the dynamic eqn
        with a numerical approach, i.e., Runge-Kutta.

        Parameters:
        -----------
        ksoil2stem : float
            conductance from soil to stem, mmol m-2 s-1 MPa-1
        psi_soil_prev : float
            soil water potential from the previous timestep, MPa
        flux_to_leaf : float
            flux from stem to the leaf
        psi_stem_prev : float
            stem water potential from the previous timestep, MPa

        Returns:
        -------
        psi_stem : float
            new stem water potential from the previous timestep, MPa

        References:
        -----------
        * Xu et al. (2016) New Phytol, 212: 80–95. doi:10.1111/nph.14009; see
          appendix and code
        """

        ap = -(self.AL * ksoil2stem / self.Cs)
        bp = (self.AL * ksoil2stem * psi_soil_prev - flux_to_leaf) / self.Cs
        psi_stem = ((ap * psi_stem_prev + bp) * \
                    np.exp(ap * self.timestep_sec)-bp) / ap

        return psi_stem

    def update_sw_bucket(self, precip, water_loss, sw_prev):
        """
        Update the simple bucket soil water balance

        Parameters:
        -----------
        precip : float
            precipitation (kg m-2 s-1)
        water_loss : float
            flux of water out of the soil (transpiration (kg m-2 timestep-1))
        sw_prev : float
            volumetric soil water from the previous timestep (m3 m-3)
        soil_volume : float
            volume soil water bucket (m3)

        Returns:
        -------
        sw : float
            new volumetric soil water (m3 m-3)
        """
        loss = water_loss * c.MMOL_2_MOL * c.MOL_WATER_2_G_WATER * \
                c.G_TO_KG * self.timestep_sec

        delta_sw = precip - loss
        sw = min(self.theta_sat, \
                 sw_prev + delta_sw / (self.soil_volume * c.M_2_MM))

        return sw

    def calc_plc(self, kp):
        """
        Calculates the percent loss of conductivity, PLC (-)

        Parameters:
        -----------
        kp : float
            plant hydraulic conductance (mmol m-2 s-1 MPa-1)

        Returns:
        -------
        plc : float
            percent loss of conductivity (-)

        """
        return 100.0 * (1.0 - kp / self.kp_sat)

    def fsig_tuzet(self, psi_leaf):
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
          1097–1116

        """
        num = 1.0 + np.exp(self.sf * self.psi_f)
        den = 1.0 + np.exp(self.sf * (self.psi_f - psi_leaf))
        fw = num / den

        return fw

    def calc_ksoil(self, psi_soil):
        """
        Calculate soil hydraulic conductance (mol m-1 s-1 MPa-1)

        Parameters:
        -----------
        psi_soil : float
            soil water potential, MPa

        Returns:
        --------
        Ksoil : float
            soil hydraulic conductance, mmol m-2 s-1 MPa-1

        References:
        -----------
        * Duursma et al. (2008) Tree Physiology 28, 265–276, eqn 9, 8, 7
        """

        # A simple equation relating Ks to psi_s is given by (Campbell 1974)
        Ks = self.Ksat * (self.psi_e / psi_soil)**(2.0 + 3.0 / self.b)
        if isclose(psi_soil, 0.0):
            Ks = self.Ksat

        # the radius of a cylinder of soil to which the root has access, n
        rcyl = 1.0 / np.sqrt(np.pi * self.Lv)

        # root length index, m root m-3 soil surface
        Rl = self.Lv * self.soil_depth
        Ksoil = (Rl / self.lai) * 2. * np.pi * Ks / np.log(rcyl / self.rroot)

        return Ksoil


def plot_time_to_mortality(out, timestep=15):

    if timestep == 15:
        ndays = out.t / 96
    elif timestep == 30:
        ndays = out.t / 96 * 2
    elif timestep == 60:
        ndays = out.t / 96 * 4

    cb = ['#377eb8', '#ff7f00', '#4daf4a', \
          '#f781bf', '#a65628', '#984ea3',\
          '#999999', '#e41a1c', '#dede00']

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()

    ln1 = ax1.plot(ndays, out.psi_leaf, ls="-", color=cb[2], label="Leaf")
    ln2 = ax1.plot(ndays, out.psi_stem, ls="-", color=cb[1], label="Stem")
    ln3 = ax1.plot(ndays, out.psi_soil, ls="-", color=cb[0], label="Soil")
    ln4 = ax2.plot(ndays, out.plc, ls='-', color=cb[6],
                   label="PLC")

    # added these three lines
    lns = ln1 + ln2 + ln3 + ln4
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=(0.05,0.08), ncol=2)
    #ax1.legend(numpoints=1, loc="best")
    ax2.set_ylabel(r'PLC (%)')
    ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("Water potential (MPa)")
    fig.savefig("plots/time_to_mortality.pdf", bbox_inches='tight', pad_inches=0.1)

def plot_swp_sw(out):

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    ax1 = fig.add_subplot(111)

    ax1.plot(out.sw, out.psi_soil, "b.", label="Soil")

    ax1.set_xlabel("Volumetric soil water content (m$^{3}$ m$^{-3}$)")
    ax1.set_ylabel("Soil water potential (MPa)")
    #ax1.legend(numpoints=1, loc="best")
    fig.savefig("plots/sw_swp.pdf", bbox_inches='tight', pad_inches=0.1)

def plot_transpiration(out):

    conv = c.MMOL_2_MOL * c.MOL_WATER_2_G_WATER * c.G_TO_KG * \
            c.SEC_2_HLFHR

    trans = []
    for i in range(0, len(out), 48):
        vals = out["Eplant"][i:i+48]
        vals = vals[~np.isnan(vals)]
        if len(vals) > 0:
            aet = np.sum(vals * conv)
            trans.append(aet)

    cb = ['#377eb8', '#ff7f00', '#4daf4a', \
          '#f781bf', '#a65628', '#984ea3',\
          '#999999', '#e41a1c', '#dede00']

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    ax1 = fig.add_subplot(111)
    ax1.plot(trans, ls="-", color=cb[1])
    ax1.set_ylabel("Transpiration (mm d$^{-1}$)")
    ax1.set_xlabel("Time (days)")
    fig.savefig("plots/transpiration.pdf", bbox_inches='tight', pad_inches=0.1)

def plot_cwd(out, timestep=15):

    if timestep == 15:
        ndays = out.t / 96
    elif timestep == 30:
        ndays = out.t / 96 * 2
    elif timestep == 60:
        ndays = out.t / 96 * 4
    """
    cwd = []
    cum_sumx = 0.0
    for i in range(0, len(out), 48):
        pet = np.sum(out.pet[i:i+48] * c.SEC_2_HLFHR)
        aet = np.sum(out["Eplant"][i:i+48] * c.MMOL_2_MOL * \
                     c.MOL_WATER_2_G_WATER * c.G_TO_KG * \
                     c.SEC_2_HLFHR)
        cum_sumx += pet - aet
        cwd.append(cum_sumx)
    """
    cwd = []
    cum_sumx = 0.0
    for i in range(len(out)):
        pet = out.pet[i] * c.SEC_2_HLFHR
        aet = out["Eplant"][i] * c.MMOL_2_MOL * c.MOL_WATER_2_G_WATER * \
                c.G_TO_KG * c.SEC_2_HLFHR
        cum_sumx += pet - aet
        cwd.append(cum_sumx)

    cb = ['#377eb8', '#ff7f00', '#4daf4a', \
          '#f781bf', '#a65628', '#984ea3',\
          '#999999', '#e41a1c', '#dede00']

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    ax1 = fig.add_subplot(111)
    #ax1.set_xlim(48)
    ax1.plot(ndays, cwd, ls="-", color=cb[1], label="Eplant")

    ax1.set_ylabel("CWD (mm)")
    ax1.set_xlabel("Time (days)")
    ax1.legend(numpoints=1, loc="best")
    fig.savefig("plots/cwd.pdf", bbox_inches='tight', pad_inches=0.1)


if __name__ == "__main__":

    time_step = 30

    met = generate_met_data(Tmin=10, Tmax=30.0, RH=30, ndays=700,
                            time_step=time_step)

    psi_stem0 = 0. # initial stem water potential, MPa
    AL = 6.        # plant leaf area, m2
    p50 = -3.      # xylem pressure inducing 50% loss of hydraulic conductivity
                   # due to embolism, MPa
    psi_f = -2.    # reference potential for Tuzet model, MPa
    gmin = 18.     # minimum stomatal conductance, mmol m-2 s-1
    #Cl = 10000.    # leaf capacitance, mmol MPa-1 (total plant)
    #Cs = 120000.   # stem capacitance, mmol MPa-1
    Cl = 0.1*150000    # leaf capacitance, mmol MPa-1 (total plant)
    Cs = 0.9*150000   # stem capacitance, mmol MPa-1
    g1 = 4.0       # sensitivity of stomatal conductance to the assimilation
                   # rate, kPa
    theta_sat = 0.5
    sw0 = 0.4
    F = Canopy(g1=g1)
    D = Desica(psi_stem0=psi_stem0, AL=AL, p50=p50, psi_f=psi_f, gmin=gmin,
               Cl=Cl, Cs=Cs, F=F, g1=g1, nruns=2, theta_sat=theta_sat,
               sw0=sw0, stop_dead=True)
    out = D.run_simulation(met)

    plot_time_to_mortality(out, time_step)
    plot_swp_sw(out)
    plot_transpiration(out)
    plot_cwd(out, time_step)
