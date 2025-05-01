"""
Tests for the `utils` module: mostly concerned with map/mask handling.
"""

import numpy as np
import unittest

from dpluspy import utils 


class TestMaskManips(unittest.TestCase):

    def test_(self):
        pass


class TestBedFiles(unittest.TestCase):

    def test_(self):
        pass


class TestBedGraphFiles(unittest.TestCase):

    def test_(self):
        pass


class TestSubsetting(unittest.TestCase):

    def test_subset_means(self):
        pop_ids = ['A', 'B', 'C']
        means = [np.arange(6)]
        to_A = utils.subset_means(means, pop_ids, ['A'])
        self.assertEqual(to_A[0], np.array([0]))
        to_C = utils.subset_means(means, pop_ids, ['C'])
        self.assertEqual(to_C[0], np.array([5]))
        to_AB = utils.subset_means(means, pop_ids, ['A', 'B'])
        self.assertTrue(np.all(to_AB[0] == np.array([0, 1, 3])))
        # Flipping the order of `to_pops` should give a different result
        to_AC = utils.subset_means(means, pop_ids, ['A', 'C'])
        self.assertTrue(np.all(to_AC[0] == np.array([0, 2, 5])))     
        to_CA = utils.subset_means(means, pop_ids, ['C', 'A'])
        self.assertTrue(np.all(to_CA[0] == np.array([5, 2, 0])))

    def test_subset_varcovs(self):
        pop_ids = ['A', 'B', 'C']
        varcovs = [np.arange(36).reshape((6, 6))]
        to_A = utils.subset_varcovs(varcovs, pop_ids, ['A'])
        to_A_expected = np.array([[0]])
        self.assertTrue(np.all(to_A == to_A_expected))
        to_C = utils.subset_varcovs(varcovs, pop_ids, ['C'])
        to_C_expected = np.array([[35]])
        self.assertTrue(np.all(to_C == to_C_expected))
        to_AC = utils.subset_varcovs(varcovs, pop_ids, ['A', 'C'])
        to_AC_expected = np.array([[0, 2, 5], [12, 14, 17], [30, 32, 35]])
        self.assertTrue(np.all(to_AC == to_AC_expected))
        to_CA = utils.subset_varcovs(varcovs, pop_ids, ['C', 'A'])
        to_CA_expected = np.array([[35, 32, 30], [17, 14, 12], [5, 2, 0]])
        self.assertTrue(np.all(to_CA == to_CA_expected))

    def test_subset_stats_by_bin(self):
        pass


class TestMapFunctions(unittest.TestCase):

    def test_map_functions(self):
        # Map functions should be inverses
        self.assertAlmostEqual(
            utils._map_function(utils._inverse_map_function(1e-8)), 1e-8)
        self.assertAlmostEqual(
            utils._map_function(utils._inverse_map_function(1e-2)), 1e-2)
        self.assertAlmostEqual(
            utils._map_function(utils._inverse_map_function(1)), 1)
        self.assertAlmostEqual(
            utils._inverse_map_function(utils._map_function(1e-8)), 1e-8)
        self.assertAlmostEqual(
            utils._inverse_map_function(utils._map_function(1e-2)), 1e-2)
        self.assertAlmostEqual(
            utils._inverse_map_function(utils._map_function(0.45)), 0.45)


