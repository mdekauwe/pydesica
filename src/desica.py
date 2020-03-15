#!/usr/bin/env python
# coding: utf-8

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
  10971116.
* Xu X, Medvigy D, Powers JS, Becknell JM, Guan K (2016) Diversity in plant
  hydraulic traits explains seasonal and inter-annual variations of vegetation
  dynamics in seasonally dry tropical forests. New Phytologist, 212, 8095.

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
from calc_pet import calc_net_radiation, calc_pet_energy, calc_fao_pet
import constants as c

class Desica(object):

    def __init__(self, plc_dead=88., soil_depth=1.0, ground_area=1.0,
                 met_timestep=30., sf=8., g1=4., Cs=100000., b=6.,
                 Cl=10000., kp_sat=4., p50=-4., psi_f=-2., s50=30., gmin=10,
                 psi_leaf0=-1., psi_stem0=-0.5, theta_sat=0.5, sw0=0.5, AL=2.5,
                 psi_e=-0.8*c.KPA_2_MPA, Ksat=20., Lv=10000., F=None,
                 keep_wet=False, stop_dead=True, rroot=1E-06, FAO=False,
                 force_refilling=False):

        self.keep_wet = keep_wet
        self.stop_dead = stop_dead
        self.plc_dead = plc_dead
        self.soil_depth = soil_depth # depth of soil bucket, m
        self.ground_area = ground_area # m
        self.soil_volume = self.ground_area * self.soil_depth # m3
        self.met_timestep = met_timestep
        self.sf = sf # sensitivity parameter, MPa-1
        self.g1 = g1 # sensitivity of stomatal conductance to the assimilation
                     # rate (-)
        self.Cs = Cs # stem capacitance, mmol kg-1 MPa-1
        self.Cl = Cl # leaf capacitance, mmol m-2 leaf MPa-1
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
        self.AL = AL   # plant leaf area, m2
        self.lai = AL / self.ground_area # leaf area index, m2 m-2
        self.b = b # empirical coefficient related to the clay content of the
                   # soil (Cosby et al. 1984).

        self.psi_e = psi_e # air entry point water potential (MPa)
        self.Ksat = Ksat # saturated conductivity, mol m-1 s-1 MPa-1
        self.Lv = Lv # root length density, m m-3
        self.F = F
        self.rroot = rroot # mean radius of water absorbing roots, m
        self.timestep_sec = 60. * self.met_timestep
        self.FAO = FAO
        self.force_refilling = force_refilling

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
        day_of_death = -999.
        (n, out) = self.initialise_model(met)

        hod = 1
        doy = 0

        for i in range(1, n):

            if self.force_refilling:

                # psi_stem can refill 30% of distance to psi_soil
                distance = 0.3

                # Impose pre-dawn (6 am) xylem refilling
                if hod == 12 and out.psi_stem[i-1] > -2:

                    # assuming psi_stem can refill X% of distance to psi_soil
                    refill = np.abs(out.psi_stem[i-1] - out.psi_soil[i-1]) * \
                                    distance
                    out.psi_stem[i-1] += refill

                    # Ensure we can't refill above psi_soil
                    out.psi_stem[i-1] = min(out.psi_stem[i-1],
                                            out.psi_soil[i-1])

                    # assume complete refilling
                    #out.psi_stem[i-1] = out.psi_soil[i-1]

            out = self.run_timestep(i, met, out)

            rnet = calc_net_radiation(i, hod, met.lat[i], met.lon[i],
                                      met.sw_rad[i], met.tair[i], met.ea[i])

            if self.FAO:
                out.pet[i] = calc_fao_pet(rnet, met.vpd[i], met.tair[i])
            else:
                out.pet[i] = calc_pet_energy(rnet)

            # Stop the simulation if we've died, i.e. reached P88
            if self.stop_dead:
                plc = self.calc_plc(out.kplant[i])
                if plc > self.plc_dead:
                    if self.met_timestep == 15:
                        day_of_death = i / 96.
                    elif self.met_timestep == 30:
                        day_of_death = i / 48.
                    elif self.met_timestep == 60:
                        day_of_death = i / 24.
                    break

            out.hod[i] = hod
            out.doy[i] = doy
            out.year[i] = met.year[i]
            hod += 1
            if hod > 47:
                hod = 0
                doy += 1

                # calculate diel range in psi_leaf
                #days_psi_leaf = out.psi_leaf[i-47:i].values
                #idx = np.argwhere(met.par[i-47:i] > 0.0).flatten()
                #diel_psi_leaf = -1 * np.fabs(np.min(days_psi_leaf[idx]) - \
                #                             np.max(days_psi_leaf[idx]))
                #print(diel_psi_leaf)

        out["plc"] = self.calc_plc(out.kplant)

        # mmol s-1
        out["Eplant"] = self.lai * out.Eleaf

        out["t"] = np.arange(1, n+1)

        out = out[~np.isnan(out.year)]

        out["cwd"] = calc_cwd(out)

        return (out, day_of_death)

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
        out.gsw[0] = 0.0
        out.plc[0] = 0.0
        out.kplant[0] = self.kp_sat * self.fsig_hydr(out.psi_stem[0])
        # soil hydraulic conductance (mmol m-2 s-1 MPa-1)
        out.ksoil[0] = self.calc_ksoil(out.psi_soil[0])

        out.hod[0] = 0
        out.doy[0] = 0
        out.year[0] = met.year.iloc[0]

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
        out = pd.DataFrame({'year':dummy,
                            'doy':dummy,
                            'hod':dummy,
                            'Eleaf':dummy,
                            'plc':dummy,
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
                            'gsw':dummy,
                            'cwd':dummy})

        return out

    def run_timestep(self, i, met, out):

        out = self.calc_conductances(out, i)

        # modified Tuzet model of stomatal conductance
        mult = (self.g1 / met.Ca[i]) * self.fsig_tuzet(out.psi_leaf[i-1])

        # Calculate photosynthesis and stomatal conductance
        gsw = self.F.canopy(met.Ca[i], met.tair[i], met.par[i],
                            met.vpd[i], mult)

        out.gsw[i] = gsw

        # Don't add gmin, instead use it as the lower boundary
        gsw = max(self.gmin, c.mol_2_mmol * gsw)

        # Leaf transpiration assuming perfect coupling, mmol m-2 s-1
        out.Eleaf[i] = gsw * (met.vpd[i] / met.press[i])

        # Calculate the leaf water potential assuming psi_stem is constant
        out.psi_leaf[i] = self.calc_psi_leaf(out.kstem2leaf[i],
                                             out.psi_stem[i-1],
                                             out.psi_leaf[i-1],
                                             out.Eleaf[i])

        # Flux from stem to leaf (mmol s-1) = change in leaf storage,
        # plus transpiration
        out.flux_to_leaf[i] = self.calc_flux_to_leaf(out.psi_leaf[i],
                                                     out.psi_leaf[i-1],
                                                     out.Eleaf[i],
                                                     out.kstem2leaf[i])

        # Update stem water potential assuming soil water potential is constant
        out.psi_stem[i] = self.update_stem_wp(out.ksoil2stem[i],
                                              out.psi_soil[i-1],
                                              out.flux_to_leaf[i],
                                              out.psi_stem[i-1])

        # Flux from the soil to the stem = change in storage + flux_to_leaf
        out.flux_to_stem[i] = self.calc_flux_to_stem(out.psi_stem[i],
                                                     out.psi_stem[i-1],
                                                     out.flux_to_leaf[i],
                                                     out.ksoil2stem[i])

        # Update psi_stem
        out.psi_stem[i-1] = out.psi_stem[i]
        #out.psi_stem[i] = self.update_stem_again(out.ksoil2stem[i],
        #                                         out.psi_soil[i-1],
        #                                         out.Eleaf[i],
        #                                         out.psi_stem[i-1])

        # Flux from the soil to the stem = change in storage + flux_to_leaf
        #out.flux_to_stem[i] = self.calc_flux_to_stem_again(out.psi_stem[i],
        #                                                   out.psi_stem[i-1],
        #                                                   out.Eleaf[i],
        #                                                   out.ksoil2stem[i])

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

        # Plant hydraulic conductance. NB. depends on stem water potential 
        # from the previous timestep.
        # mmol m-2 leaf s-1 MPa-1
        out.kplant[i] = self.kp_sat * self.fsig_hydr(out.psi_stem[i-1])

        # Conductance from root surface to the stem water pool (assumed to be
        # halfway to the leaves)
        # mmol m-2 ground area s-1 MPa-1
        kroot2stem = 2.0 * (out.kplant[i] * self.lai)

        # Conductance from soil to stem water store (mmol m-2 s-1 MPa-1)
        # (conductances combined in series)
        out.ksoil2stem[i] = 1.0 / (1.0 / out.ksoil[i-1] + 1.0 / kroot2stem)

        # Conductance from stem water store to the leaves
        # assumning the water pool is halfway up the stem
        # mmol m-2 ground area s-1 MPa-1
        out.kstem2leaf[i] = 2.0 * (out.kplant[i] * self.lai)

        return out

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

        relk = max(1.0e-9, min(1.0, relk))

        return (relk)

    def calc_psi_leaf(self, kstem2leaf, psi_stem_prev, psi_leaf_prev, Eleaf):
        """
        Calculate leaf water potential, MPa

        This is a simplified equation based on Xu et al., using the water
        potentials from the previous timestep and the fact that we increase
        the temporal resolution to get around the need to solve the dynamic eqn
        with a numerical approach, i.e., Runge-Kutta.

        NB. kstem2leaf which is the stem conductance in CABLE needs to be a
            function of sapwood, tree height.
        NB. transpiration is being scaled up by LAI, but in CABLE this is
            done already.

        Parameters:
        -----------
        kstem2leaf : float
            conductance from stem to leaf, mmol m-2 ground area s-1 MPa-1
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
        * Xu et al. (2016) New Phytol, 212: 8095. doi:10.1111/nph.14009; see
          appendix and code. Can write the dynamic equation as:
          dpsi_leaf_dt = b + a*psi_leaf
        """

        # MPa
        psi_leaf_min = -20.0

        # scale up leaf-specific capacitance (mmol m-2 ground area MPa-1 s-1)
        leaf_capac = self.Cl * self.lai

        # there is conductance in the trunk
        if kstem2leaf > 1E-09:
            ap = - kstem2leaf / leaf_capac # unitless
            # MPa
            bp = (psi_stem_prev * kstem2leaf - self.lai * Eleaf) / leaf_capac
            psi_leaf = ((ap * psi_leaf_prev + bp) * \
                        np.exp(ap * self.timestep_sec) - bp) / ap # MPa

        # No conductance in the trunk, change in psi_leaf is due only to
        # transpiration
        else:
            # MPa
            psi_leaf = (psi_leaf_prev - \
                        self.lai * Eleaf * self.timestep_sec) / leaf_capac

        if psi_leaf < psi_leaf_min:
            psi_leaf = psi_leaf_min # MPa

        return psi_leaf

    def calc_flux_to_leaf(self, psi_leaf, psi_leaf_prev, Eleaf, kstem2leaf):
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
            flux from stem to the leaf, mmol m−2 s−1
        """
        # scale up leaf-specific capacitance (mmol m-2 ground area MPa-1 s-1)
        leaf_capac = self.Cl * self.lai

        # there is conductance in the trunk
        if kstem2leaf > 1E-09:

            # sapflow rate from stem to leaf within the time step
            # mmol m−2 s−1
            J_sl =  (psi_leaf - psi_leaf_prev) * \
                     leaf_capac / self.timestep_sec + (self.lai * Eleaf)
        else:
            # no conductance in the trunk

            # sapflow rate from stem to leaf within the time step
            # mmol m−2 s−1
            J_sl = 0.0

        return J_sl

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
        * Xu et al. (2016) New Phytol, 212: 8095. doi:10.1111/nph.14009; see
          appendix and code
        """
        # MPa
        psi_stem_min = 2.0 * self.p50

        # mmol m-2 MPa-1 s-1
        stem_capac = self.Cs * self.scale_up_stem_capac()

        # plant cannot take up water, change of psi_stem is solely due to
        # flux_to_leaf (J_sl)
        if ksoil2stem == 0.0:
            psi_stem = psi_stem_prev - \
                            flux_to_leaf * self.timestep_sec / stem_capac # MPa
        # plant can take up water
        else:
            ap = - ksoil2stem / stem_capac  # unitless
            bp = (ksoil2stem * psi_soil_prev - flux_to_leaf) / stem_capac # MPa
            psi_stem = ((ap * psi_stem_prev + bp) * \
                        np.exp(ap * self.timestep_sec) - bp) / ap # MPa

        if psi_stem < psi_stem_min:
            psi_stem = psi_stem_min # MPa

        return psi_stem

    def calc_flux_to_stem(self, psi_stem, psi_stem_prev, flux_to_leaf,
                          ksoil2stem):
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
            flux from soil to the stem, mmol m−2 s−1
        """
        # mmol m-2 MPa-1 s-1
        stem_capac = self.Cs * self.scale_up_stem_capac()

        if ksoil2stem == 0.0:
            # mmol m−2 s−1
            J_sr = 0.0
        else:
            # mmol m−2 s−1
            J_sr = (psi_stem - psi_stem_prev) * \
                    stem_capac / self.timestep_sec + flux_to_leaf

        return J_sr

    def update_stem_again(self, ksoil2stem, psi_soil_prev, Eleaf,
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
            conductance from soil to stem, mmol m-2 ground area s-1 MPa-1
        psi_soil_prev : float
            soil water potential from the previous timestep, MPa
        Eleaf : float
            transpiration, mmol m-2 s-1
        psi_stem_prev : float
            stem water potential from the previous timestep, MPa

        Returns:
        -------
        psi_stem : float
            new stem water potential from the previous timestep, MPa

        References:
        -----------
        * Xu et al. (2016) New Phytol, 212: 8095. doi:10.1111/nph.14009; see
          appendix and code
        """

        # MPa
        psi_stem_min = 3.0 * self.p50

        # mmol m-2 MPa-1 s-1
        stem_capac = self.Cs * self.scale_up_stem_capac()

        # scale up leaf-specific capacitance (mmol m-2 ground area MPa-1 s-1)
        leaf_capac = self.Cl * self.lai

        # mmol m-2 MPa-1 s-1
        total_capac = stem_capac + leaf_capac

        # plant cannot take up water, change of psi_stem is solely due to
        # flux_to_leaf (J_sl)
        if ksoil2stem == 0.0:
            # MPa
            psi_stem = psi_stem_prev - \
                        (self.lai * Eleaf) * self.timestep_sec / total_capac
        # plant can take up water
        else:
            ap = - ksoil2stem / total_capac  # unitless
            bp = (ksoil2stem * psi_soil_prev - \
                  (self.lai * Eleaf)) / total_capac  # MPa
            psi_stem = ((ap * psi_stem_prev + bp) * \
                        np.exp(ap * self.timestep_sec) - bp) / ap # MPa


        if psi_stem < psi_stem_min:
            psi_stem = psi_stem_min # MPa

        return psi_stem

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
        * Duursma et al. (2008) Tree Physiology 28, 265276, eqn 10
        """
        return self.psi_e * (sw / self.theta_sat)**-self.b # MPa

    def calc_flux_to_stem_again(self, psi_stem, psi_stem_prev, transpiration,
                                ksoil2stem):
        """
        Calculate the flux from the root to the stem, i.e. the root water
        uptake (mmol s-1) = change in stem storage plus flux_to_leaf

        Parameters:
        -----------
        psi_stem : float
            stem water potential, MPa
        psi_stem_prev : float
            stem water potential from the previous timestep, MPa
        Eleaf : float
            transpiration, mmol m-2 s-1
        ksoil2stem : float
            conductance from soil to stem, mmol m-2 ground area s-1 MPa-1

        Returns:
        -------
        flux_to_stem : float
            flux from soil to the stem, mmol m−2 s−1
        """

        # mmol m-2 MPa-1 s-1
        stem_capac = self.Cs * self.scale_up_stem_capac()

        if ksoil2stem == 0.0:
            # mmol m−2 s−1
            J_sr = 0.0
        else:
            # mmol m−2 s−1
            J_sr = (psi_stem - psi_stem_prev) * \
                    stem_capac / self.timestep_sec + (transpiration * self.lai)

        return J_sr

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
        delta_sw = (precip * self.timestep_sec) - loss

        sw = min(self.theta_sat, \
                 sw_prev + delta_sw / (self.soil_volume * c.M_2_MM))
        #print(sw * (self.soil_volume * c.M_2_MM))
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
        return 100.0 * (1.0 - kp / self.kp_sat) # unitless

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
          10971116

        """
        num = 1.0 + np.exp(self.sf * self.psi_f)
        den = 1.0 + np.exp(self.sf * (self.psi_f - psi_leaf))
        fw = num / den
        fw = max(1.0e-9, min(1.0, fw)) # unitless

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
        * Duursma et al. (2008) Tree Physiology 28, 265276, eqn 9, 8, 7
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

    def scale_up_stem_capac(self):
        """
        Calculate relative plant conductance as a function of xylem pressure

        References:
        ==========
        * Togashi et al. (2015) Morphological and moisture availability controls
          of the leaf area‐to‐sapwood area ratio: analysis of measurements on
          Australian trees. Ecology and Evolution, 5, 1263–1270.
        * Simard et al. ( 2011) Mapping forest canopy height globally with
          spaceborne lidar, J. Geophys. Res., 116, G04021.
        * Bowman et al. (2005) Sap flow rates and sapwood density are critical
          factors in within‐ and between‐tree variation in CO2 efflux from stems
          of mature Dacrydium cupressinum trees. New Phytologist, 167: 815-828.
        """

        # We need to scale up the measurements to infer a total stem capacitance
        # These assumption will *clearly* need to be improved!

        # From Togashi et al. it looks like there is no obvious reln btw LA:SA
        # and moisture so taking a rough average LA:SA
        la_sa = 5000.0 # m2 m-2

        # We need a height and in the longterm we will use POP, but for now from
        # Simard et al. and sampling where Butt et al. say there are Eucs...we
        # can make an extremely horrible assumption that height is 20 m
        height = 20.0 # m

        # From Bowman we are taking a rough sapwood density (fig 5), which seems
        # to agree with Xu et al.
        sapwood_density = 500.0 # kg m-3

        # scalar to *roughly* convert stem capacitance. We are broadly matching
        # the spirit of Xu et al. but dropping the root capacitance they include.
        # We could of course include this but given this is pretty rough, this
        # seems OK
        # (kg m-2)
        capac_conv = self.lai * height / la_sa * sapwood_density

        return capac_conv

def plot_time_to_mortality(odir, out, timestep=15, to_screen=False, pft=None):

    if timestep == 15:
        ndays = out.t / 96.
    elif timestep == 30:
        ndays = out.t / 48.
    elif timestep == 60:
        ndays = out.t / 24.

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

    if to_screen:
        plt.show()
    else:
        if pft is None :
            fig.savefig("%s/time_to_mortality.pdf" % (odir),
                        bbox_inches='tight', pad_inches=0.1)
        else:
            fig.savefig("%s/time_to_mortality_%s.pdf" % (odir, pft),
                        bbox_inches='tight', pad_inches=0.1)
        plt.close('all')

def plot_swp_sw(odir, out, to_screen=False, year=None):

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
    if to_screen:
        plt.show()
    else:
        if year is None:
            fig.savefig("%s/sw_swp.pdf" % (odir), bbox_inches='tight',
                        pad_inches=0.1)
        else:
            fig.savefig("%s/sw_swp_%d.pdf" % (odir, year), bbox_inches='tight',
                        pad_inches=0.1)
        plt.close('all')

def plot_swp_ksoil(odir, out, year=None):

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

    ax1.plot(out.psi_soil, out.ksoil, "b.", label="Ksoil")

    ax1.set_ylabel("Soil hyd. cond. (mmol m$^{-2}$ s$^{-1}$ MPa$^{-1}$)")
    ax1.set_xlabel("Soil water potential (MPa)")
    #ax1.legend(numpoints=1, loc="best")
    ax1.set_ylim(0, 10)
    if year is None:
        fig.savefig("%s/swp_ksoil.pdf" % (odir), bbox_inches='tight',
                    pad_inches=0.1)
    else:
        fig.savefig("%s/swp_ksoil_%d.pdf" % (odir, year), bbox_inches='tight',
                    pad_inches=0.1)
    plt.close('all')

def plot_transpiration(odir, out, to_screen=False, year=None):

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

    if to_screen:
        plt.show()
    if year is None:
        fig.savefig("%s/transpiration.pdf" % (odir), bbox_inches='tight',
                    pad_inches=0.1)
    else:
        fig.savefig("%s/transpiration_%d.pdf" % (odir, year),
                    bbox_inches='tight', pad_inches=0.1)
    plt.close('all')

def plot_transpiration_and_pet(odir, out, year=None):

    conv = c.MMOL_2_MOL * c.MOL_WATER_2_G_WATER * c.G_TO_KG * \
            c.SEC_2_HLFHR

    trans = []
    pet = []
    for i in range(0, len(out), 48):
        vals = out["Eplant"][i:i+48]
        vals = vals[~np.isnan(vals)]
        if len(vals) > 0:
            aet = np.sum(vals * conv)
            trans.append(aet)
        vals = out["pet"][i:i+48]
        vals = vals[~np.isnan(vals)]
        if len(vals) > 0:
            pet.append(np.sum(vals * c.SEC_2_HLFHR))

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
    ax1.plot(trans, ls="-", color=cb[1], label="AET")
    ax1.plot(pet, ls="-", color=cb[2], label="PET")
    ax1.set_ylabel("(mm d$^{-1}$)")
    ax1.set_xlabel("Time (days)")

    if year is None:
        fig.savefig("%s/transpiration_and_pet.pdf" % (odir), bbox_inches='tight',
                    pad_inches=0.1)
    else:
        fig.savefig("%s/transpiration_and_pet_%d.pdf" % (odir, year),
                    bbox_inches='tight', pad_inches=0.1)
    plt.close('all')

def plot_cwd(odir, out, timestep=15, year=None):

    if timestep == 15:
        ndays = out.t / 96
    elif timestep == 30:
        ndays = out.t / 48.
    elif timestep == 60:
        ndays = out.t / 24.
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
    cwd = calc_cwd(out)

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
    ax1.plot(ndays, cwd, ls="-", color=cb[1])

    ax1.set_ylabel("Accumulated CWD (mm)")
    ax1.set_xlabel("Time (days)")

    if year is None:
        fig.savefig("%s/cwd.pdf" % (odir), bbox_inches='tight', pad_inches=0.1)
    else:
        fig.savefig("%s/cwd_%d.pdf" % (odir, year), bbox_inches='tight',
                    pad_inches=0.1)
    plt.close('all')

def plot_gmin_sensitvity(odir, gmin, death, year=None):

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
    ax1.plot(gmin, death, ls="-", color=cb[1])
    ax1.set_ylabel("Time to mortality (days)")
    ax1.set_xlabel("g$_{min}$ (mmol m$^{-2}$ s$^{-1}$)")

    if year is None:
        fig.savefig("%s/gmin_sensitivity.pdf" % (odir), bbox_inches='tight',
                    pad_inches=0.1)
    else:
        fig.savefig("%s/gmin_sensitivity_%d.pdf" % (odir, year),
                    bbox_inches='tight', pad_inches=0.1)
    plt.close('all')

def plot_sw(odir, out, time_step=30, year=None):

    if time_step == 15:
        ndays = out.t / 96
    elif time_step == 30:
        ndays = out.t / 48.
    elif time_step == 60:
        ndays = out.t / 24.

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
    ax1.plot(ndays, out.sw, ls="-", color=cb[1])
    ax1.set_ylabel("SWC (m$^{3}$ m$^{-3}$)")
    ax1.set_xlabel("Time (days)")

    if year is None:
        fig.savefig("%s/sw.pdf" % (odir), bbox_inches='tight',
                    pad_inches=0.1)
    else:
        fig.savefig("%s/sw_%d.pdf" % (odir, year), bbox_inches='tight',
                    pad_inches=0.1)
    plt.close('all')

def calc_cwd(out):

    cwd = []
    cum_sumx = 0.0

    dx = 0.0
    dy = 0.0
    hod = 0
    for i in range(len(out)):
        pet = out.pet[i] * c.SEC_2_HLFHR
        aet = out["Eplant"][i] * c.MMOL_2_MOL * c.MOL_WATER_2_G_WATER * \
                c.G_TO_KG * c.SEC_2_HLFHR
        cum_sumx += pet - aet
        cwd.append(cum_sumx)
        dx += pet
        dy += aet

        hod += 1
        if hod > 47:
            hod = 0.0
            #print(dx, dy)
            dx = 0.0
            dy = 0.0

    return cwd

if __name__ == "__main__":

    time_step = 30

    # Examine how time to death changes as gmin is increased?
    do_sensitivity = False
    lat = -35.76
    lon = 148.0
    met = generate_met_data(Tmin=10, Tmax=30.0, RH=30, ndays=300,
                            lat=lat, lon=lon, time_step=time_step)
    psi_e = -0.8 * \
             c.KPA_2_MPA   # Sand, MPa
    b = 6.
    #psi_e = -0.68 * \
    #         c.KPA_2_MPA   # Sand, MPa
    #b = 2.79               # Sand, SW retention curve param
    #psi_e = -1.32 * \
    #         c.KPA_2_MPA  # Sandy clay loam, MPa
    #b = 6.77              # Sandy clay loam, SW retention curve param
    #psi_e = -3.17 * \
    #         c.KPA_2_MPA  # Silty clay clay loam, MPa
    #b = 10.39             # Silty clay, SW retention curve param
    kp_sat = 1.686987 #4   # Tim Brodribb pers comm
    psi_stem0 = 0.         # initial stem water potential, MPa
    AL = 2.0#6.            # plant leaf area, m2
    p50 = -3.002384#-4.    # xylem pressure inducing 50% loss of hydraulic
                           # conductivity due to embolism, MPa
    psi_f = -2.455474 #-3. # reference potential for Tuzet model, MPa
    gmin = 10.             # minimum stomatal conductance, mmol m-2 s-1
    Cl = 342.904821        # leaf capacitance, mmol m-2 leaf MPa-1
    Cs = 53266.089926      # stem capacitance, mmol kg-1 MPa-1
    s50 = 35.26948
    g1 = 3.154297 #4.0     # sensitivity of stomatal conductance to the
                           # assimilation rate (-)
    g0 = 0.0
    theta_J = 0.85
    Rd25 = 0.92
    Q10 = 1.92
    Vcmax25 = 85.0 #50.0
    Jmax25 = 142.0 #100.
    Eav = 58550.0
    deltaSv = 629.26
    Eaj = 29680.
    deltaSj = 631.88
    FAO = False
    year = 2000

    F = Canopy(g1=g1, g0=g0, theta_J=theta_J, Rd25=Rd25, Q10=Q10,
               Vcmax25=Vcmax25, Jmax25=Jmax25, Eav=Eav, deltaSv=deltaSv,
               Eaj=Eaj, deltaSj=deltaSj)
    D = Desica(psi_stem0=psi_stem0, AL=AL, p50=p50, psi_f=psi_f, gmin=gmin,
               Cl=Cl, Cs=Cs, F=F, g1=g1, stop_dead=True, psi_e=psi_e,
               FAO=FAO, kp_sat=kp_sat, b=b, s50=s50)
    out, day_of_death = D.run_simulation(met)

    odir = "plots"
    if not os.path.exists(odir):
        os.makedirs(odir)

    plot_time_to_mortality(odir, out, time_step, to_screen=True)
    #plot_transpiration(odir, out, to_screen=True)
    #plot_swp_sw(odir, out, to_screen=True)
    """
    plot_swp_sw(odir, out)
    plot_swp_ksoil(odir, out)
    plot_transpiration(odir, out)
    plot_transpiration_and_pet(odir, out)
    plot_cwd(odir, out, time_step)
    plot_sw(odir, out, time_step)

    odir = "outputs"
    if not os.path.exists(odir):
        os.makedirs(odir)

    ofname = os.path.join(odir, "desica_out.csv")
    out.to_csv(ofname, index=False)


    if do_sensitivity:
        death = []
        gminx = np.linspace(10, 50, 10)
        for gmin in gminx:
            D = Desica(psi_stem0=psi_stem0, AL=AL, p50=p50, psi_f=psi_f,
                       gmin=gmin, Cl=Cl, Cs=Cs, F=F, g1=g1,
                       stop_dead=True, FAO=FAO)
            out, day_of_death = D.run_simulation(met)
            death.append(day_of_death)
        plot_gmin_sensitvity(odir, gminx, death)
    """
