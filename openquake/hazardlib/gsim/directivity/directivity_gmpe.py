# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2020, GEM Foundation
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


import copy
import numpy as np
from openquake.hazardlib.gsim.base import (GMPE, registry, CoeffsTable,
ADMITTED_STR_PARAMETERS, ADMITTED_FLOAT_PARAMETERS, ADMITTED_SET_PARAMETERS)
from openquake.hazardlib import const
from openquake.hazardlib.imt import from_string, PGV, PGA, SA
from openquake.hazardlib.gsim.directivity.bayless_2020 import BaylessEtAl2020

directivity_registry = {
"BaylessEtAl2020": BaylessEtAl2020
}



class GMPEDirectivity(GMPE):
    """
    This is a wrapper GMPE that can apply directivity adjustments

    :param string gmpe:
        The name of a GMPE class used for the calculation.
    :param string directivity_model:
        The name of the direcitivty model class to be used for calculation.
    :param params:
        A dictionary where the key defines the required modification and the
        value is a list with the required parameters.
    """
    REQUIRES_SITES_PARAMETERS = set()
    REQUIRES_DISTANCES = set()
    REQUIRES_RUPTURE_PARAMETERS = {"apply_directivity",}
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = {}
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = ''
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = {const.StdDev.TOTAL,
                                            const.StdDev.INTER_EVENT,
                                            const.StdDev.INTRA_EVENT}
    DEFINED_FOR_TECTONIC_REGION_TYPE = ''
    DEFINED_FOR_REFERENCE_VELOCITY = None

    def __init__(self, gmpe, directivity_model, **kwargs):
        super().__init__(**kwargs)

        # Create the original GMPE
        self.params = kwargs  # non-gmpe parameters
        self.gmpe = registry[gmpe](**kwargs)
        self.directivity_model = \
            directivity_registry[directivity_model](**kwargs)
        self.set_parameters()
        # Add on the required rupture, distance and site parameters
        self.REQUIRES_RUPTURE_PARAMETERS = (
            self.REQUIRES_RUPTURE_PARAMETERS |
            self.directivity_model.REQUIRES_RUPTURE_PARAMETERS)

        self.REQUIRES_DISTANCES = (
            self.REQUIRES_DISTANCES |
            self.directivity_model.REQUIRES_DISTANCES)

        # Add on the required rupture, distance and site parameters
        self.REQUIRES_SITES_PARAMETERS = (
            self.REQUIRES_SITES_PARAMETERS |
            self.directivity_model.REQUIRES_SITES_PARAMETERS)

        # Define only for intensity types supported by both the GMPE and the
        # directivity model
        self.DEFINED_FOR_INTENSITY_MEASURE_TYPES = (
            self.DEFINED_FOR_INTENSITY_MEASURE_TYPES &
            self.directivity_model.DEFINED_FOR_INTENSITY_MEASURE_TYPES)

        # Add on the directivity model parameters
        self.mean = None

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """

        """
        if not rup.apply_directivity:
            # For the givem rupture the directivity is turned off - so simply
            # return the GMPE
            return self.gmpe.get_mean_and_stddevs(sites, rup, dists, imt,
                                                  stddev_types)
        mean, stddevs = self.gmpe.get_mean_and_stddevs(sites, rup, dists, imt,
            [const.StdDev.INTER_EVENT, const.StdDev.INTRA_EVENT])
        # Apply directivity adjustment
        median_adjustment, stddevs =\
            self.directivity_model.get_directivity_adjustment(sites, rup,
                                                              dists, imt,
                                                              stddevs)
        # Adjust the mean
        mean += median_adjustment
        # If the total standard deviation is required
        if (len(stddev_types) == 1) and\
                (stddev_types[0] == const.StdDev.TOTAL):
            stddevs = [np.sqrt(stddevs[0] ** 2.0 + stddevs[1] ** 2.0)]
        return mean, stddevs
