# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2012-2020 GEM Foundation
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenQuake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>.
"""
Model imports :class:`BaylessEtAl2020Directivity`
"""

import numpy as np
from scipy.special import ndtr
from openquake.hazardlib import const
from openquake.hazardlib.imt import PGA, SA
from openquake.hazardlib.gsim.base import CoeffsTable
from openquake.hazardlib.gsim.directivity.base import DirectivityModel


class BaylessEtAl2020(DirectivityModel):
    """
    Implements the directivity model of Bayless, Somerville and Skarlatoudis
    (2020) "A Rupture Directivity Adjustment Model Applicable to the NGA-West2
    Ground Motion Models and Complex Fault Geometries", USGS Technical Report

    Note here that the rake is used to define the style of faulting in general
    in OQ, which means that SOF is defined as strike slip for rakes < 45. or
    rake > 135., reverse for 45.0 <= rake <= 135. and normal for -45 >= rake
    >= -135. However, the Bayless et al. (2020) directivity model can in theory
    support dipping ruptures with a strongly strike-slip style of faulting
    (i.e. abs(rake) < 30.0 or > 150.0). In the current implementation this
    will switch the style of faulting, so care should be taken not to use
    this model on highly oblique slipping ruptures.
    """
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = const.TRT.ACTIVE_SHALLOW_CRUST

    #: Reference to a :class:`intensity measure component type
    #: <openquake.hazardlib.const.IMC>` this GSIM can calculate mean
    #: and standard
    #: deviation for.
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = {PGA, SA}

    #: Set of
    #: :class:`standard deviation types <openquake.hazardlib.const.StdDev>`
    #: this GSIM can calculate.
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = {const.StdDev.INTER_EVENT,
                                            const.StdDev.INTRA_EVENT}

    #: Set of site parameters names this GSIM needs. The set should include
    #: strings that match names of the attributes of a :class:`site
    #: <openquake.hazardlib.site.Site>` object.
    #: Those attributes are then available in the
    #: :class:`SitesContext` object with the same names.
    REQUIRES_SITES_PARAMETERS = set()

    #: Set of rupture parameters (excluding distance information) required
    #: by GSIM. Supported parameters are:
    REQUIRES_RUPTURE_PARAMETERS = {"mag", "strike", "dip", "rake", "gc_length",
                                   "width", "hypo_depth", "ztor", "zbor",
                                   "gc2u_updip_hypo", "apply_directivity"}

    #: Set of types of distance measures between rupture and sites. Possible
    #: values are:
    #: All the distances are available from the :class:`DistancesContext`
    #: object attributes with same names. Values are in kilometers.
    REQUIRES_DISTANCES = {"ry0", "gc2u", "gc2t"}

    # 
    def get_directivity_adjustment(self, sites, rup, dists, imt, stddevs):
        """

        """
        if str(imt) == "PGA":
            # No directivity adjustment for PGA
            return np.zeros(stddevs[0].shape), stddevs
        # Get style-of-faulting indicator 
        sof = self._get_sof(rup.rake)
        # Get Rmax
        rmax = self._get_rmax(sof, rup.mag)
        median_adjustment = self.get_median_adjustment_factor(rup, dists,
                                                              sites, imt,
                                                              sof, rmax)
        # Adjusted within event standard deviation
        phi = self.get_within_event_stddev(imt, stddevs[1],
                                            rmax, dists.rrup)
        return median_adjustment, [stddevs[0], phi]

    def get_median_adjustment_factor(self, rup, dists, sites, imt, sof, rmax):
        """
        """

        # Get the a(M, T) and b(M, T) coefficients
        aval, bval = self._get_a_b(rup.mag, imt.period, sof)

        dval = max((rup.hypo_depth - rup.ztor)/ np.sin(np.radians(rup.dip)),
                   3.0)
        # T-value of bottom of rupture is effectively the distance from the
        # surface projection of the lower edge of the rupture to the up-dip
        # projection of the upper edge of the rupture to the surface
        tbot = rup.zbor / np.tan(np.radians(rup.dip))
        if not rup.ztor:
            # Rupture reaches the surface, so GC2 T ordinate is zero-centred
            # on fault top edge
            t_offset = 0.0
        else:
            # Rupture is buried, so T ordinate should be re-centred on up-dip
            # projection of rupture (i.e. trace)
            t_offset = rup.ztor / np.tan(np.radians(rup.dip))
        # GC2uprime translates the original gc2u (in which 0 is centred on
        # the starting point of the rupture) and translates it so that
        # gc2u = 0 for the up-dip projection of the hypocentre
        gc2uprime = dists.gc2u - rup.gc2u_updip_hypo
        # S-value is equivalent to GC2 U for sites within the length of the
        # fault
        sval = self._get_sval(gc2uprime, rup)
        # Get fdist
        fdist = self._get_fdist(rup.mag, imt.period, sof, dists.ry0,
                                dists.gc2t + t_offset, rmax)
        # Get fs2
        fs2 = self._get_fs2(sof, rup.rake, dval, sval)
        # Get ftheta
        ftheta = self._get_ftheta(sof, rup.rake, dists.ry0,
                                  dists.gc2t + t_offset, gc2uprime)
        # Get fphi
        fphi = self._get_fphi(sof, rup.dip, dists.gc2t + t_offset,
                              tbot + t_offset, rup.zbor)
        # Get fG (geometric directivity predictor)
        return (aval + bval * (fs2 * ftheta * fphi)) * fdist 

    def get_within_event_stddev(self, imt, phi, rmax, rrup):
        """
        Returns the modifies within event standard deviation
        """
        if ("PGA" in str(imt)) or (imt.period <= 0.2):
            # No adjustment
            return phi
        phi_red = np.zeros(rrup.shape)
        phi_red[rrup < rmax] = self.COEFFS[imt]["e1"]
        return np.sqrt(phi ** 2. - phi_red ** 2.)

    def _get_sval(self, gc2uprime, rup):
        """
        Returns the S-value, defined as being equivalent to GC2-U within
        the length of the fault and and equal to the GC2-U coordinate of the
        nearest rupture trace end point otherwise
        """
        # S-value is equivalent to GC2 U for sites within the length of the
        # fault
        sval = np.copy(gc2uprime)
        # For sites off the edge of the fault it is equal to the GC2 U ordinate
        # of the nearest endpoint
        sval[sval < -rup.gc2u_updip_hypo] = -rup.gc2u_updip_hypo
        diff_gc2u = rup.gc_length - rup.gc2u_updip_hypo
        sval[sval > diff_gc2u] = diff_gc2u
        return sval

    def _get_a_b(self, mag, period, sof):
        # Returns period-dependent model factors
        if sof == 1:
            d1, d2 = self.CONSTANTS["d1_sof1"], self.CONSTANTS["d2_sof1"]
        else:
            d1, d2 = self.CONSTANTS["d1_sof2"], self.CONSTANTS["d2_sof2"]
        fg0 = d1 + d2 * mag
        bmax = self.CONSTANTS["b1"] + self.CONSTANTS["b2"] * mag
        tpeak = 10.0 ** (self.CONSTANTS["c1"] + self.CONSTANTS["c2"] * mag)
        bmt = bmax * np.exp((np.log10(period / tpeak) ** 2) /
                            (-2.0 * (self.CONSTANTS["sigma_g"] ** 2.)))
        amt = -bmt * fg0
        return amt, bmt

    def _get_fdist(self, mag, period, sof, ry0, gc2t, rmax):
        # Returns the directivity taper from SOF, magnitude and distance
        rtaper = np.sqrt(gc2t ** 2. + ry0 ** 2.)
        fdist = np.zeros(ry0.shape)
        idx = rtaper <= rmax
        if np.any(idx):
            fdist[idx] = 1.0 - np.exp(((-4.0 * rmax) / rtaper[idx]) + 4.0)
        return fdist

    def _get_rmax(self, sof, mag):
        # Returns the Rmax term described by equations 3e and 3f
        if sof == 1:
            rsof = 0.0
        else:
            rsof = 20.
        if mag > 7.0:
            rmax = 80.0 - rsof
        else:
            rmax = 20.0 * mag - 60.0 - rsof
        return rmax

    def _get_sof(self, rake):
        # Returns the style-of-faulting index (1 for strike slip,
        # 2 for dip-slip)
        if (np.fabs(rake) < 30.0) or (np.fabs(rake) > 150.):
            # Strike slip
            return 1
        else:
            return 2

    def _get_ftheta(self, sof, rake, ry0, gc2t, gc2u):
        # Returns the azimuthal predictor ftheta - equations 5a and 5b
        if sof == 1:
            theta = np.fabs(np.arctan(gc2t / gc2u))
            return np.fabs(np.cos(2.0 * theta))
        else:
            tmin = 10.0 * (np.fabs(np.cos(np.radians(rake))) + 1.0)
            theta = np.fabs(gc2t)
            theta[theta < tmin] = tmin
            # On the hanging wall side 
            theta[gc2t > 0] = tmin
            theta = np.arctan(np.fabs(theta / ry0))
            return np.sin(theta)

    def _get_fphi(self, sof, dip, gc2t, tbot, dbot):
        # Returns the azimuthal predictor fphi - equations 6a and 6b
        if sof == 1:
            # phi = 0, so f_phi = cos(2.0 * 0) = 1
            return np.ones(gc2t.shape)
        else:
            phi = np.zeros(gc2t.shape)
            idx = gc2t < 0.0
            if np.any(idx):
                phi[idx] = np.radians(dip - 90.) +\
                    np.arctan((np.fabs(gc2t[idx]) + tbot) / dbot)
            idx = np.logical_and(gc2t >= 0.0, gc2t <= tbot)
            if np.any(idx):
                phi[idx] = np.radians(90.0 - dip) -\
                    np.arctan((tbot - gc2t[idx]) / dbot)
            idx = gc2t > tbot
            if np.any(idx):
                phi[idx] = np.radians(90.0 - dip) +\
                    np.arctan((gc2t[idx] - tbot) / dbot)
            phi[phi > (np.pi / 4.0)] = np.pi / 4.0
            return np.cos(2.0 * phi)

    def _get_fs2(self, sof, rake, dval, sval):
        # Returns the rupture travel distance term fs2 (equations 4a and 4b)
         
        s_cos_rake = sval * np.cos(np.radians(rake))
        s_2 = np.sqrt(dval ** 2. + s_cos_rake ** 2.)
        if sof == 1:
            return np.log(np.clip(s_2, -np.inf, 465.))
        else:
            s_2[s_cos_rake < 0.0] = dval
            return np.log(np.clip(s_2, -np.inf, 188.0))

    CONSTANTS = {"b1": 0.5469, "b2": -0.0336, "c1": -1.2090, "c2": 0.2858,
                 "d1_sof1": -4.8300, "d1_sof2": -1.5415,
                 "d2_sof1": 0.9928, "d2_sof2": 0.3946, "sigma_g": 0.4653}

    COEFFS = CoeffsTable(sa_damping=5, table="""\
    imt       e1
    0.01   0.000
    0.20   0.000
    0.25   0.008
    0.30   0.020
    0.40   0.035
    0.50   0.051
    0.75   0.067
    1.00   0.080
    1.50   0.084
    2.00   0.093
    3.00   0.110
    4.00   0.139
    5.00   0.166
    7.50   0.188
    10.0   0.199
    """)
