# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2012-2016 GEM Foundation
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
Module :mod:`~openquake.hazardlib.calc.filters` contain filter functions for
calculators.

Filters are functions (or other callable objects) that should take generators
and return generators. There are two different kinds of filter functions:

1. Source-site filters. Those functions take a generator of two-item tuples,
   each pair consists of seismic source object (that is, an instance of
   a subclass of :class:`~openquake.hazardlib.source.base.BaseSeismicSource`)
   and a site collection (instance of
   :class:`~openquake.hazardlib.site.SiteCollection`).
2. Rupture-site filters. Those also take a generator of pairs, but in this
   case the first item in the pair is a rupture object (instance of
   :class:`~openquake.hazardlib.source.rupture.Rupture`). The second element in
   generator items is still site collection.

The purpose of both kinds of filters is to limit the amount of calculation
to be done based on some criteria, like the distance between the source
and the site. So common design feature of all the filters is the loop over
pairs of the provided generator, filtering the sites collection, and if
there are no items left in it, skipping the pair and continuing to the next
one. If some sites need to be considered together with that source / rupture,
the pair gets generated out, with a (possibly) :meth:`limited
<openquake.hazardlib.site.SiteCollection.filter>` site collection.

Consistency of filters' input and output stream format allows several filters
(obviously, of the same kind) to be chained together.

Filter functions should not make assumptions about the ordering of items
in the original generator or draw more than one pair at once. Ideally, they
should also perform reasonably fast (filtering stage that takes longer than
the actual calculation on unfiltered collection only decreases performance).

Module :mod:`openquake.hazardlib.calc.filters` exports one distance-based
filter function of each kind (see :func:`source_site_distance_filter` and
:func:`rupture_site_distance_filter`) as well as "no operation" filters
(:func:`source_site_noop_filter` and :func:`rupture_site_noop_filter`).
"""


def filter_sites_by_distance_to_rupture(rupture, integration_distance, sites):
    """
    Filter out sites from the collection that are further from the rupture
    than some arbitrary threshold.

    :param rupture:
        Instance of :class:`~openquake.hazardlib.source.rupture.Rupture`
        that was generated by :meth:
        `openquake.hazardlib.source.base.BaseSeismicSource.iter_ruptures`
        of an instance of this class.
    :param integration_distance:
        Threshold distance in km.
    :param sites:
        Instance of :class:`openquake.hazardlib.site.SiteCollection`
        to filter.
    :returns:
        Filtered :class:`~openquake.hazardlib.site.SiteCollection`.

    This function is similar to :meth:`openquake.hazardlib.source.base.BaseSeismicSource.filter_sites_by_distance_to_source`.
    The same notes about filtering criteria apply. Site
    should not be filtered out if it is not further than the integration
    distance from the rupture's surface projection along the great
    circle arc (this is known as Joyner-Boore distance, :meth:`
    openquake.hazardlib.geo.surface.base.BaseQuadrilateralSurface.get_joyner_boore_distance`).
    """
    jb_dist = rupture.surface.get_joyner_boore_distance(sites.mesh)
    return sites.filter(jb_dist <= integration_distance)

def source_site_distance_filter(integration_distance):
    """
    Source-site filter based on distance.

    :param integration_distance:
        Threshold distance in km, this value gets passed straight to
        :meth:`openquake.hazardlib.source.base.BaseSeismicSource.filter_sites_by_distance_to_source`
        which is what is actually used for filtering.
    """
    def filter_func(sources_sites):
        for source, sites in sources_sites:
            s_sites = source.filter_sites_by_distance_to_source(
                integration_distance, sites
            )
            if s_sites is None:
                continue
            yield source, s_sites
    return filter_func


def rupture_site_distance_filter(integration_distance):
    """
    Rupture-site filter based on distance.

    :param integration_distance:
        Threshold distance in km, this value gets passed straight to
        :func:`openquake.hazardlib.calc.filters.filter_sites_by_distance_to_rupture`
        which is what is actually used for filtering.
    """
    def filter_func(ruptures_sites):
        for rupture, sites in ruptures_sites:
            r_sites = filter_sites_by_distance_to_rupture(
                rupture, integration_distance, sites)
            if r_sites is None:
                continue
            yield rupture, r_sites
    return filter_func


#: Transparent source-site "no-op" filter -- behaves like a real filter
#: but never filters anything out and doesn't have any overhead.
source_site_noop_filter = lambda sources_sites: sources_sites

#: Rupture-site "no-op" filter, same as :func:`source_site_noop_filter`.
rupture_site_noop_filter = lambda ruptures_sites: ruptures_sites
