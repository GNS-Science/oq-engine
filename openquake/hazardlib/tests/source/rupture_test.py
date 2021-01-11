# The Hazard Library
# Copyright (C) 2012-2020 GEM Foundation
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import unittest
import numpy
import os
from openquake.hazardlib import const
from openquake.hazardlib.geo import Point, Line
from openquake.hazardlib.geo.surface.planar import PlanarSurface
from openquake.hazardlib.tom import PoissonTOM
from openquake.hazardlib.source.rupture import BaseRupture, \
    ParametricProbabilisticRupture, NonParametricProbabilisticRupture,\
    RuptureHypocenterDistribution
from openquake.hazardlib.pmf import PMF
from openquake.hazardlib.geo.mesh import Mesh
from openquake.hazardlib.geo.surface.simple_fault import SimpleFaultSurface


def make_rupture(rupture_class, **kwargs):
    default_arguments = {
        'mag': 5.5,
        'rake': 123.45,
        'tectonic_region_type': const.TRT.STABLE_CONTINENTAL,
        'hypocenter': Point(5, 6, 7),
        'surface': PlanarSurface(11, 12, Point(0, 0, 1), Point(1, 0, 1),
                                 Point(1, 0, 2), Point(0, 0, 2)),
    }
    default_arguments.update(kwargs)
    kwargs = default_arguments
    rupture = rupture_class(**kwargs)
    for key in kwargs:
        if key != 'pmf':
            # for pmf .pmf is a numpy array whereas pmf is a PMF instance
            assert getattr(rupture, key) is kwargs[key], (
                getattr(rupture, key), kwargs[key])
    return rupture


class RuptureCreationTestCase(unittest.TestCase):
    def assert_failed_creation(self, rupture_class, exc, msg, **kwargs):
        with self.assertRaises(exc) as ae:
            make_rupture(rupture_class, **kwargs)
        self.assertEqual(str(ae.exception), msg)

    def test_negative_magnitude(self):
        self.assert_failed_creation(
            BaseRupture, ValueError,
            'magnitude must be positive',
            mag=-1
        )

    def test_zero_magnitude(self):
        self.assert_failed_creation(
            BaseRupture, ValueError,
            'magnitude must be positive',
            mag=0
        )

    def test_probabilistic_rupture_negative_occurrence_rate(self):
        self.assert_failed_creation(
            ParametricProbabilisticRupture, ValueError,
            'occurrence rate must be positive',
            occurrence_rate=-1, temporal_occurrence_model=PoissonTOM(10)
        )

    def test_probabilistic_rupture_zero_occurrence_rate(self):
        self.assert_failed_creation(
            ParametricProbabilisticRupture, ValueError,
            'occurrence rate must be positive',
            occurrence_rate=0, temporal_occurrence_model=PoissonTOM(10)
        )

    def test_rupture_topo(self):
        rupture = make_rupture(BaseRupture, hypocenter=Point(5, 6, -2))
        self.assertEqual(rupture.hypocenter.depth, -2)


class ParametricProbabilisticRuptureTestCase(unittest.TestCase):
    def test_get_probability_one_or_more(self):
        rupture = make_rupture(ParametricProbabilisticRupture,
                               occurrence_rate=1e-2,
                               temporal_occurrence_model=PoissonTOM(10))
        self.assertAlmostEqual(
            rupture.get_probability_one_or_more_occurrences(), 0.0951626
        )

    def test_get_probability_one_occurrence(self):
        rupture = make_rupture(ParametricProbabilisticRupture,
                               occurrence_rate=0.4,
                               temporal_occurrence_model=PoissonTOM(10))
        self.assertAlmostEqual(rupture.get_probability_one_occurrence(),
                               0.0732626)

    def test_sample_number_of_occurrences(self):
        time_span = 20
        rate = 0.01
        num_samples = 2000
        tom = PoissonTOM(time_span)
        rupture = make_rupture(ParametricProbabilisticRupture,
                               occurrence_rate=rate,
                               temporal_occurrence_model=tom)
        numpy.random.seed(37)
        mean = sum(rupture.sample_number_of_occurrences()
                   for i in range(num_samples)) / float(num_samples)
        self.assertAlmostEqual(mean, rate * time_span, delta=2e-3)

    def test_get_probability_no_exceedance(self):
        rupture = make_rupture(ParametricProbabilisticRupture,
                               occurrence_rate=0.01,
                               temporal_occurrence_model=PoissonTOM(50))
        poes = numpy.array([[0.9, 0.8, 0.7], [0.6, 0.5, 0.4]])
        pne = rupture.get_probability_no_exceedance(poes)
        numpy.testing.assert_allclose(
            pne,
            numpy.array([[0.6376282, 0.6703200, 0.7046881],
                         [0.7408182, 0.7788008, 0.8187308]])
        )


class Cdppvalue(unittest.TestCase):

    def make_rupture_fordpp(self, rupture_class, **kwargs):
        # Create the rupture surface.
        upper_seismogenic_depth = 0.
        lower_seismogenic_depth = 15.
        dip = 90.
        mesh_spacing = 1.
        fault_trace_start = Point(10., 45.2)
        fault_trace_end = Point(10., 45.919457)
        fault_trace = Line([fault_trace_start, fault_trace_end])
        default_arguments = {
            'mag': 7.2,
            'rake': 0.,
            'tectonic_region_type': const.TRT.STABLE_CONTINENTAL,
            'hypocenter': Point(10.0, 45.334898, 10),
            'surface': SimpleFaultSurface.from_fault_data(
                fault_trace, upper_seismogenic_depth, lower_seismogenic_depth,
                dip=dip, mesh_spacing=mesh_spacing),
            'rupture_slip_direction': 0.
        }
        default_arguments.update(kwargs)
        kwargs = default_arguments
        rupture = rupture_class(**kwargs)
        for key in kwargs:
            assert getattr(rupture, key) is kwargs[key]
        return rupture

    def test_get_dppvalue(self):
        rupture = self.make_rupture_fordpp(
            ParametricProbabilisticRupture, occurrence_rate=0.01,
            temporal_occurrence_model=PoissonTOM(50))
        # Load the testing site.
        data_path = os.path.dirname(__file__)
        filename = os.path.join(
            data_path, "./data/geo_cycs_ss3_testing_site.csv")
        data = numpy.genfromtxt(filename,
                                dtype=float, delimiter=',', names=True,
                                skip_header=6675, skip_footer=6673)

        for loc in range(len(data)):
            lon = data[loc][0]
            lat = data[loc][1]
            ref_dpp = data[loc][2]
            dpp = rupture.get_dppvalue(Point(lon, lat))

            self.assertAlmostEqual(dpp, ref_dpp, delta=0.1)

    @unittest.skipUnless('OQ_RUN_SLOW_TESTS' in os.environ, 'slow')
    def test_get_cdppvalue(self):
        rupture = self.make_rupture_fordpp(
            ParametricProbabilisticRupture, occurrence_rate=0.01,
            temporal_occurrence_model=PoissonTOM(50))
        # Load the testing site.
        data_path = os.path.dirname(__file__)
        filename = os.path.join(
            data_path, "./data/geo_cycs_ss3_testing_site.csv")
        data = numpy.genfromtxt(filename,
                                dtype=float, delimiter=',', names=True,
                                skip_header=6675, skip_footer=6673)
        points = []
        for loc in range(len(data)):
            lon = data[loc][0]
            lat = data[loc][1]
            points.append(Point(lon, lat))

        mesh = Mesh.from_points_list(points)
        cdpp = rupture.get_cdppvalue(mesh)
        self.assertAlmostEqual(cdpp[0], data[0][3], delta=0.1)
        self.assertAlmostEqual(cdpp[1], data[1][3], delta=0.1)


class NonParametricProbabilisticRuptureTestCase(unittest.TestCase):
    def assert_failed_creation(self, rupture_class, exc, msg, **kwargs):
        with self.assertRaises(exc) as ae:
            make_rupture(rupture_class, **kwargs)
        self.assertEqual(str(ae.exception), msg)

    def test_creation(self):
        pmf = PMF([(0.8, 0), (0.2, 1)])
        make_rupture(NonParametricProbabilisticRupture, pmf=pmf)

    def test_minimum_number_of_ruptures_is_not_zero(self):
        pmf = PMF([(0.8, 1), (0.2, 2)])
        self.assert_failed_creation(
            NonParametricProbabilisticRupture,
            ValueError, 'minimum number of ruptures must be zero', pmf=pmf
        )

    def test_numbers_of_ruptures_not_in_increasing_order(self):
        pmf = PMF([(0.8, 0), (0.1, 2), (0.1, 1)])
        self.assert_failed_creation(
            NonParametricProbabilisticRupture,
            ValueError,
            'numbers of ruptures must be defined in increasing order', pmf=pmf
        )

    def test_numbers_of_ruptures_not_defined_with_unit_step(self):
        pmf = PMF([(0.8, 0), (0.2, 2)])
        self.assert_failed_creation(
            NonParametricProbabilisticRupture,
            ValueError,
            'numbers of ruptures must be defined with unit step', pmf=pmf
        )

    def test_get_probability_no_exceedance(self):
        pmf = PMF([(0.7, 0), (0.2, 1), (0.1, 2)])
        poes = numpy.array([[0.9, 0.8, 0.7], [0.6, 0.5, 0.4]])
        rup = make_rupture(NonParametricProbabilisticRupture, pmf=pmf)
        pne = rup.get_probability_no_exceedance(poes)
        numpy.testing.assert_allclose(
            pne,
            numpy.array([[0.721, 0.744, 0.769], [0.796, 0.825, 0.856]])
        )

    def test_sample_number_of_occurrences(self):
        pmf = PMF([(0.7, 0), (0.2, 1), (0.1, 2)])
        rup = make_rupture(NonParametricProbabilisticRupture, pmf=pmf)
        numpy.random.seed(123)

        n_samples = 50000
        n_occs = numpy.array([
            rup.sample_number_of_occurrences() for i in range(n_samples)
        ])

        p_occs_0 = float(len(n_occs[n_occs == 0])) / n_samples
        p_occs_1 = float(len(n_occs[n_occs == 1])) / n_samples
        p_occs_2 = float(len(n_occs[n_occs == 2])) / n_samples

        self.assertAlmostEqual(p_occs_0, 0.7, places=2)
        self.assertAlmostEqual(p_occs_1, 0.2, places=2)
        self.assertAlmostEqual(p_occs_2, 0.1, places=2)

class RuptureHypocenterDistributionTestCase(unittest.TestCase):

    def setUp(self):
        # Use a simple vertical planar rupture
        top_left = Point(0.0, 0.0, 0.0)
        top_right = Point(0.0, 0.1, 0.0)
        bottom_right = Point(0.0, 0.1, 10.)
        bottom_left = Point(0.0, 0.0, 10.0)
        self.plane = PlanarSurface.from_corner_points(top_left, top_right,
                                                      bottom_right,
                                                      bottom_left)

    def test_valid_rupture_hypo_distribution(self):
        hypo_dist_input = [[0.2, 0.25, 0.5 * 0.25],
                           [0.5, 0.25, 0.5 * 0.5],
                           [0.8, 0.25, 0.5 * 0.25],
                           [0.2, 0.75, 0.5 * 0.25],
                           [0.5, 0.75, 0.5 * 0.5],
                           [0.8, 0.75, 0.5 * 0.25]]

        hypo_dist = RuptureHypocenterDistribution(hypo_dist_input)
        # Check hypocentre probabilities
        numpy.testing.assert_array_almost_equal(
            numpy.array([0.125, 0.25, 0.125, 0.125, 0.25, 0.125]),
            hypo_dist.hypo_probs, 7)
        # Check hypocentre positions
        target_hypo_pos = numpy.array([[0.2, 0.25],
                                    [0.5, 0.25],
                                    [0.8, 0.25],
                                    [0.2, 0.75],
                                    [0.5, 0.75],
                                    [0.8, 0.75]])
        numpy.testing.assert_array_almost_equal(target_hypo_pos,
                                                hypo_dist.hypo_pos, 7)
        # Check slip distribution
        self.assertListEqual(hypo_dist.slip_list, [])
        # Try repeat with a different slip distribution
        slip_dist2 = [[0.0, 0.5], [180.0, 0.5]]
        hypo_dist2 = RuptureHypocenterDistribution(hypo_dist_input,
                                                   slip_dist2)
        numpy.testing.assert_array_almost_equal(numpy.array([[0.0, 0.5],
                                                             [180., 0.5]]),
                                                hypo_dist2.slip_list, 7.0)
        # Check the rupture counts
        print(len(hypo_dist), len(hypo_dist2))
        self.assertEqual(len(hypo_dist), 6)
        self.assertEqual(len(hypo_dist2), 12)

    def test_generate_hypocentres_vertical(self):
        hypo_dist_input = [[0.2, 0.25, 0.5 * 0.25],
                           [0.5, 0.25, 0.5 * 0.5],
                           [0.8, 0.25, 0.5 * 0.25],
                           [0.2, 0.75, 0.5 * 0.25],
                           [0.5, 0.75, 0.5 * 0.5],
                           [0.8, 0.75, 0.5 * 0.25]]

        hypo_dist = RuptureHypocenterDistribution(hypo_dist_input)
        hypo_lons = numpy.zeros(6, dtype=float)
        hypo_lats = numpy.array([0.02, 0.05, 0.08, 0.02, 0.05, 0.08])
        hypo_depths = numpy.array([2.5, 2.5, 2.5, 7.5, 7.5, 7.5])
        target_hypo_probs = numpy.array([0.125, 0.25, 0.125,
                                         0.125, 0.25, 0.125])
        for i, (hypocenter, slip, prob) in \
            enumerate(hypo_dist.get_parameters_probabilities(self.plane)):
            # Note here that the length of the planar fault is not the same
            # as the distance between the top two points. This creates a
            # discrepancy between the predicted and observed latitudes
            self.assertAlmostEqual(hypocenter.longitude, hypo_lons[i], 3)
            self.assertAlmostEqual(hypocenter.latitude, hypo_lats[i], 3)
            self.assertAlmostEqual(hypocenter.depth, hypo_depths[i], 5)
            self.assertAlmostEqual(prob, target_hypo_probs[i], 7)
            self.assertAlmostEqual(slip, 0.0, 7)

#    def test_generate_updip_vertical(self):
#        hypo_dist_input = [[0.2, 0.25, 0.5 * 0.25],
#                           [0.5, 0.25, 0.5 * 0.5],
#                           [0.8, 0.25, 0.5 * 0.25],
#                           [0.2, 0.75, 0.5 * 0.25],
#                           [0.5, 0.75, 0.5 * 0.5],
#                           [0.8, 0.75, 0.5 * 0.25]]
#
#        hypo_dist = RuptureHypocenterDistribution(hypo_dist_input)
#        hypocenters, hypo_probs = hypo_dist.get_hypocenters_probabilities(
#            self.plane)
#
#        hypo_lons = numpy.zeros(6, dtype=float)
#        hypo_lats = numpy.array([0.02, 0.05, 0.08, 0.02, 0.05, 0.08])
#        hypo_depths = numpy.zeros(6, dtype=float)
#        for i, hypocenter in enumerate(hypocenters):
#            updip = hypo_dist.project_updip(hypocenter, self.plane)
#            self.assertAlmostEqual(updip.longitude, 0.0, 7)
#            self.assertAlmostEqual(updip.latitude, hypo_lats[i], 5)
#            self.assertAlmostEqual(updip.depth, 0.0, 7)
#
#    def test_hypocenter_updip_gc2_vertical(self):
#        target_gc2t = numpy.zeros(6, dtype=float)
#        target_gc2u = numpy.array([0.2 * self.plane.length,
#                                   0.5 * self.plane.length,
#                                   0.8 * self.plane.length,
#                                   0.2 * self.plane.length,
#                                   0.5 * self.plane.length,
#                                   0.8 * self.plane.length])
#
#        hypo_dist_input = [[0.2, 0.25, 0.5 * 0.25],
#                           [0.5, 0.25, 0.5 * 0.5],
#                           [0.8, 0.25, 0.5 * 0.25],
#                           [0.2, 0.75, 0.5 * 0.25],
#                           [0.5, 0.75, 0.5 * 0.5],
#                           [0.8, 0.75, 0.5 * 0.25]]
#
#        hypo_dist = RuptureHypocenterDistribution(hypo_dist_input)
#        for i, (hypocenter, slip, prob) in \
#            enumerate(hypo_dist.get_parameters_probabilities(self.plane)):
#            updip = hypo_dist.project_updip(hypocenter, self.plane)
#            # GC2-T, GC2-U for hypocentre
#            gc2t, gc2u = hypo_dist.get_gc2_point(hypocenter, self.plane)
#            self.assertAlmostEqual(gc2t[0], target_gc2t[i], 7)
#            self.assertAlmostEqual(gc2u[0], target_gc2u[i], 7)
#            # GC2-T, GC2-U for updip
#            gc2t, gc2u = hypo_dist.get_gc2_point(updip, self.plane)
#            self.assertAlmostEqual(gc2t[0], target_gc2t[i], 7)
#            self.assertAlmostEqual(gc2u[0], target_gc2u[i], 7)
#
    def test_hypocenter_updip_dipping_plane(self):
        top_left = Point(0.0, 0.0, 1.0)
        top_right = Point(0.0, 0.1, 1.0)
        bottom_right = Point(0.05, 0.1, 11.0)
        bottom_left = Point(0.05, 0.0, 11.0)
        plane2 = PlanarSurface.from_corner_points(top_left, top_right,
                                                  bottom_right, bottom_left)

        # Target hypocentres
        target_hypos = numpy.array([[0.012497215, 0.019981162, 3.499442896],
                                    [0.012497219, 0.049952904, 3.499442896],
                                    [0.012497227, 0.079924647, 3.499442896],
                                    [0.037491646, 0.019981158, 8.498328687],
                                    [0.037491658, 0.049952895, 8.498328687],
                                    [0.037491680, 0.079924631, 8.498328687]])
        hypo_dist_input = [[0.2, 0.25, 0.5 * 0.25],
                           [0.5, 0.25, 0.5 * 0.5],
                           [0.8, 0.25, 0.5 * 0.25],
                           [0.2, 0.75, 0.5 * 0.25],
                           [0.5, 0.75, 0.5 * 0.5],
                           [0.8, 0.75, 0.5 * 0.25]]
        hypo_dist = RuptureHypocenterDistribution(hypo_dist_input)
        for i, (hypocenter, slip, prob) in\
            enumerate(hypo_dist.get_parameters_probabilities(plane2)):
            self.assertAlmostEqual(hypocenter.longitude, target_hypos[i, 0], 7)
            self.assertAlmostEqual(hypocenter.latitude, target_hypos[i, 1], 7)
            self.assertAlmostEqual(hypocenter.depth, target_hypos[i, 2], 7)
#        # Target updip values
#        target_updips = numpy.array([[-0.005000000, 0.019981161, 0.000000000],
#                                     [-0.005000000, 0.049952902, 0.000000000],
#                                     [-0.005000000, 0.079924643, 0.000000000],
#                                     [-0.005000000, 0.019981152, 0.000000000],
#                                     [-0.005000000, 0.049952881, 0.000000000],
#                                     [-0.005000000, 0.079924609, 0.000000000]])

#            self.assertAlmostEqual(updip.longitude, target_updips[i, 0], 7)
#            self.assertAlmostEqual(updip.latitude, target_updips[i, 1], 7)
#            self.assertAlmostEqual(updip.depth, target_updips[i, 2], 7)
#
#    def test_hypocenter_updip_gc2_dipping(self):
#        top_left = Point(0.0, 0.0, 1.0)
#        top_right = Point(0.0, 0.1, 1.0)
#        bottom_right = Point(0.05, 0.1, 11.0)
#        bottom_left = Point(0.05, 0.0, 11.0)
#        plane2 = PlanarSurface.from_corner_points(top_left, top_right,
#                                                  bottom_right, bottom_left)
#        # Target GC2 Hypocenters
#        target_hypo_gc2 = numpy.array([[1.389626847, 2.221803803],
#                                       [1.389626847, 5.554509507],
#                                       [1.389626847, 8.887215212],
#                                       [4.168880542, 2.221803380],
#                                       [4.168880542, 5.554508450],
#                                       [4.168880542, 8.887213520]])
#        # Target GC2 Updips
#        target_updip_gc2 = numpy.array([[-0.555974633, 2.221803699],
#                                        [-0.555974633, 5.554509248],
#                                        [-0.555974633, 8.887214797],
#                                        [-0.555974633, 2.221802769],
#                                        [-0.555974633, 5.554506923],
#                                        [-0.555974633, 8.887211077]])
#
#        hypo_dist_input = [[0.2, 0.25, 0.5 * 0.25],
#                           [0.5, 0.25, 0.5 * 0.5],
#                           [0.8, 0.25, 0.5 * 0.25],
#                           [0.2, 0.75, 0.5 * 0.25],
#                           [0.5, 0.75, 0.5 * 0.5],
#                           [0.8, 0.75, 0.5 * 0.25]]
#        hypo_dist = RuptureHypocenterDistribution(hypo_dist_input)
#        for i, (hypocenter, slip, prob) in\
#            enumerate(hypo_dist.get_parameters_probabilities(plane2)):
#            updip = hypo_dist.project_updip(hypocenter, plane2)
#            gc2t, gc2u = hypo_dist.get_gc2_point(hypocenter, plane2)
#            self.assertAlmostEqual(gc2t[0], target_hypo_gc2[i, 0], 7)
#            self.assertAlmostEqual(gc2u[0], target_hypo_gc2[i, 1], 7)
#            gc2t, gc2u = hypo_dist.get_gc2_point(updip, plane2)
#            self.assertAlmostEqual(gc2t[0], target_updip_gc2[i, 0], 7)
#            self.assertAlmostEqual(gc2u[0], target_updip_gc2[i, 1], 7)
#   
