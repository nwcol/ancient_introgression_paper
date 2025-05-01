"""
Tests for the `parsing` module.
"""

import numpy as np
import unittest
import os

from dpluspy import utils, parsing 


def _count_pairs(site_map, bins, weights=None):
    """
    Get binned pair counts with a naive loop.
    """
    if weights is None:
        weights = np.ones(len(site_map))
    count = np.zeros(len(bins) - 1)
    for i, (map_l, weight_l) in enumerate(zip(site_map[:-1], weights[:-1])):
        for map_r, weight_r in zip(site_map[i + 1:], weights[i + 1:]):
            distance = map_r - map_l
            if distance >= bins[0] and distance < bins[-1]:
                idx = np.digitize(distance, bins) - 1
                count[idx] += weight_l * weight_r
 
    return count
 

def _count_pairs_between(
    site_map_l, 
    site_map_r, 
    bins, 
    weights_l=None,
    weights_r=None
):
    """
    Get binned pair counts between two discontinuous map segments with a naive 
    loop.
    """
    site_map_l = np.asarray(site_map_l)
    site_map_r = np.asarray(site_map_r)
    assert len(site_map_l) > 1
    assert len(site_map_r) > 1
    if weights_l is not None:
        assert len(site_map_l) == len(weights_l)
        assert len(site_map_r) == len(weights_r)
    assert np.all(site_map_l < site_map_r[0])
    if weights_l is None:
        weights_l = np.ones(len(site_map_l))
    if weights_r is None:
        weights_r = np.ones(len(site_map_r))
    count = np.zeros(len(bins) - 1)
    for map_l, weight_l in zip(site_map_l, weights_l):
        for map_r, weight_r in zip(site_map_r, weights_r):
            distance = map_r - map_l
            if distance >= bins[0] and distance < bins[-1]:
                idx = np.digitize(distance, bins) - 1
                count[idx] += weight_l * weight_r
    
    return count 


def sample_map(r, L):
    """
    Construct a random recombination map with length `L` and average rate `r`.
    """
    return np.cumsum(np.random.uniform(0, r * 2, L))


class TestLocusPairCounting(unittest.TestCase):
    ## Randomly samples maps, weights.

    def test_counting_within(self):
        bins = np.concatenate(([0], np.logspace(-6, -1, 10)))
        site_map = sample_map(1e-8, 100)
        result = parsing._count_locus_pairs(site_map, bins)
        naive = _count_pairs(site_map, bins)
        self.assertTrue(np.all(result == naive))

        bins = np.logspace(-6, -1, 10)
        result = parsing._count_locus_pairs(site_map, bins)
        naive = _count_pairs(site_map, bins)
        self.assertTrue(np.all(result == naive))

        site_map = sample_map(1e-6, 100)
        bins = np.logspace(-6, -1, 10)
        result = parsing._count_locus_pairs(site_map, bins)
        naive = _count_pairs(site_map, bins)
        self.assertTrue(np.all(result == naive))

        site_map = sample_map(1e-2, 200)
        bins = np.logspace(-6, -1, 10)
        result = parsing._count_locus_pairs(site_map, bins)
        naive = _count_pairs(site_map, bins)
        self.assertTrue(np.all(result == naive))

        bins = np.logspace(-6, -0.35, 10)
        result = parsing._count_locus_pairs(site_map, bins)
        naive = _count_pairs(site_map, bins)
        self.assertTrue(np.all(result == naive))

    def test_counting_between(self):
        bins0 = np.concatenate(([0], np.logspace(-6, -1, 10)))
        bins = np.logspace(-6, -1, 10)

        full_map = sample_map(1e-8, 200)
        left_map, right_map = full_map[:100], full_map[100:]
        result = parsing._count_locus_pairs_between(left_map, right_map, bins0)
        naive = _count_pairs_between(left_map, right_map, bins0)
        self.assertTrue(np.all(result == naive))
        result = parsing._count_locus_pairs_between(left_map, right_map, bins)
        naive = _count_pairs_between(left_map, right_map, bins)
        self.assertTrue(np.all(result == naive))
        # Longer map distances
        full_map = sample_map(1e-4, 200)
        left_map, right_map = full_map[:100], full_map[100:]
        result = parsing._count_locus_pairs_between(left_map, right_map, bins0)
        naive = _count_pairs_between(left_map, right_map, bins0)
        self.assertTrue(np.all(result == naive))
        result = parsing._count_locus_pairs_between(left_map, right_map, bins)
        naive = _count_pairs_between(left_map, right_map, bins)
        self.assertTrue(np.all(result == naive))
        # Large seperation between segments
        right_map += 0.03
        result = parsing._count_locus_pairs_between(left_map, right_map, bins0)
        naive = _count_pairs_between(left_map, right_map, bins0)
        self.assertTrue(np.all(result == naive))
        result = parsing._count_locus_pairs_between(left_map, right_map, bins)
        naive = _count_pairs_between(left_map, right_map, bins)
        self.assertTrue(np.all(result == naive))
        # Very large seperation between segments
        right_map += 0.30
        result = parsing._count_locus_pairs_between(left_map, right_map, bins0)
        naive = _count_pairs_between(left_map, right_map, bins0)
        self.assertTrue(np.all(result == naive))
        result = parsing._count_locus_pairs_between(left_map, right_map, bins)
        naive = _count_pairs_between(left_map, right_map, bins)
        self.assertTrue(np.all(result == naive))

    def test_weighted_counting_within(self):
        # Because cumulative sums are involved, some small loss of precision is 
        # entailed; I use `np.isclose` rather than testing for equality
        weights = np.random.uniform(0, 1, 100)
        bins = np.concatenate(([0], np.logspace(-6, -1, 10)))
        site_map = sample_map(1e-8, 100)
        result = parsing._count_locus_pairs(site_map, bins, weights=weights)
        naive = _count_pairs(site_map, bins, weights=weights)
        self.assertTrue(np.all(np.isclose(result, naive)))

        bins = np.logspace(-6, -1, 10)
        site_map = sample_map(1e-8, 100)
        result = parsing._count_locus_pairs(site_map, bins, weights=weights)
        naive = _count_pairs(site_map, bins, weights=weights)
        self.assertTrue(np.all(np.isclose(result, naive)))

        bins = np.concatenate(([0], np.logspace(-6, -1, 10)))
        site_map = sample_map(1e-3, 100)
        result = parsing._count_locus_pairs(site_map, bins, weights=weights)
        naive = _count_pairs(site_map, bins, weights=weights)
        self.assertTrue(np.all(np.isclose(result, naive)))
        
        weights = np.random.choice([0, 1], size=200)
        bins = np.concatenate(([0], np.logspace(-6, -1, 10)))
        site_map = sample_map(1e-8, 200)
        result = parsing._count_locus_pairs(site_map, bins, weights=weights)
        naive = _count_pairs(site_map, bins, weights=weights)
        self.assertTrue(np.all(np.isclose(result, naive)))

        site_map = sample_map(1e-3, 200)
        result = parsing._count_locus_pairs(site_map, bins, weights=weights)
        naive = _count_pairs(site_map, bins, weights=weights)
        self.assertTrue(np.all(np.isclose(result, naive)))

    def test_weighted_counting_between(self):
        bins0 = np.concatenate(([0], np.logspace(-6, -1, 10)))
        bins = np.logspace(-6, -1, 10)

        full_map = sample_map(1e-8, 200)
        left_map, right_map = full_map[:100], full_map[100:]
        weightsl = np.random.uniform(0, 1, 100)
        weightsr = np.random.uniform(0, 1, 100)
        result = parsing._count_locus_pairs_between(
            left_map, right_map, bins0, weights_l=weightsl, weights_r=weightsr)
        naive = _count_pairs_between(
            left_map, right_map, bins0, weights_l=weightsl, weights_r=weightsr)
        self.assertTrue(np.all(np.isclose(result, naive)))
        result = parsing._count_locus_pairs_between(
            left_map, right_map, bins, weights_l=weightsl, weights_r=weightsr)
        naive = _count_pairs_between(
            left_map, right_map, bins, weights_l=weightsl, weights_r=weightsr)
        self.assertTrue(np.all(np.isclose(result, naive)))

        weightsl = np.random.choice([0, 1], size=100)
        weightsr = np.random.choice([0, 1], size=100)
        result = parsing._count_locus_pairs_between(
            left_map, right_map, bins0, weights_l=weightsl, weights_r=weightsr)
        naive = _count_pairs_between(
            left_map, right_map, bins0, weights_l=weightsl, weights_r=weightsr)
        self.assertTrue(np.all(np.isclose(result, naive)))
        result = parsing._count_locus_pairs_between(
            left_map, right_map, bins, weights_l=weightsl, weights_r=weightsr)
        naive = _count_pairs_between(
            left_map, right_map, bins, weights_l=weightsl, weights_r=weightsr)
        self.assertTrue(np.all(np.isclose(result, naive)))

        full_map = sample_map(1e-4, 200)
        left_map, right_map = full_map[:100], full_map[100:]
        result = parsing._count_locus_pairs_between(
            left_map, right_map, bins0, weights_l=weightsl, weights_r=weightsr)
        naive = _count_pairs_between(
            left_map, right_map, bins0, weights_l=weightsl, weights_r=weightsr)
        self.assertTrue(np.all(np.isclose(result, naive)))
        result = parsing._count_locus_pairs_between(
            left_map, right_map, bins, weights_l=weightsl, weights_r=weightsr)
        naive = _count_pairs_between(
            left_map, right_map, bins, weights_l=weightsl, weights_r=weightsr)
        self.assertTrue(np.all(np.isclose(result, naive)))

        full_map = sample_map(1e-3, 200)
        left_map, right_map = full_map[:100], full_map[100:]
        result = parsing._count_locus_pairs_between(
            left_map, right_map, bins0, weights_l=weightsl, weights_r=weightsr)
        naive = _count_pairs_between(
            left_map, right_map, bins0, weights_l=weightsl, weights_r=weightsr)
        self.assertTrue(np.all(np.isclose(result, naive)))
        result = parsing._count_locus_pairs_between(
            left_map, right_map, bins, weights_l=weightsl, weights_r=weightsr)
        naive = _count_pairs_between(
            left_map, right_map, bins, weights_l=weightsl, weights_r=weightsr)
        self.assertTrue(np.all(np.isclose(result, naive)))
        # All pairs are out of bin range
        right_map += 0.30
        result = parsing._count_locus_pairs_between(
            left_map, right_map, bins0, weights_l=weightsl, weights_r=weightsr)
        naive = _count_pairs_between(
            left_map, right_map, bins0, weights_l=weightsl, weights_r=weightsr)
        self.assertTrue(np.all(np.isclose(result, naive)))
        result = parsing._count_locus_pairs_between(
            left_map, right_map, bins, weights_l=weightsl, weights_r=weightsr)
        naive = _count_pairs_between(
            left_map, right_map, bins, weights_l=weightsl, weights_r=weightsr)
        self.assertTrue(np.all(np.isclose(result, naive)))

    def test_count_locus_pairs_exceptions(self):
        # Maps with zero values return zero sums
        bins = np.logspace(-6, -1, 10)
        empty = np.zeros(len(bins) - 1)
        result = parsing._count_locus_pairs(np.array([]), bins)
        self.assertTrue(np.all(empty == result))
        # Map/weight length mismatch raises an error
        with self.assertRaises(ValueError):
            site_map = np.array([0, 1e-6, 1e-4])
            weights = np.array([1, 1])
            parsing._count_locus_pairs(site_map, bins, weights=weights)
        # Maps must be monotonically increasing
        with self.assertRaises(ValueError):
            site_map = np.array([0, 1e-6, 1e-4, 1e-8])
            parsing._count_locus_pairs(site_map, bins)

    def test_count_locus_pairs_between_exceptions(self):
        bins = np.logspace(-6, -1, 10)
        empty = np.zeros(len(bins) - 1)
        result = parsing._count_locus_pairs_between(
            np.array([]), np.array([]), bins)
        self.assertTrue(np.all(empty == result))
        # Giving weights for one window only raises an error
        site_map_l = np.array([0, 1e-7])
        site_map_r = np.array([1e-6, 2e-6])
        weights = np.array([1, 0])
        with self.assertRaises(ValueError):
            parsing._count_locus_pairs_between(
                site_map_l, site_map_r, bins, weights_l=weights)
        # Map/length mismatches raise errors
        mismatch = np.array([1])
        with self.assertRaises(ValueError):
            parsing._count_locus_pairs_between(site_map_l, site_map_r, bins, 
                weights_l=weights, weights_r=mismatch)
        with self.assertRaises(ValueError):
            mismatch = np.array([1])
            parsing._count_locus_pairs_between(site_map_l, site_map_r, bins, 
                weights_l=mismatch, weights_r=mismatch)
        # The right map must have higher coords than the left
        mis_map = np.array([5e-8, 2e-7])
        with self.assertRaises(ValueError):
            parsing._count_locus_pairs_between(site_map_l, mis_map, bins)


class TestMaps(unittest.TestCase):

    def test_loading_uniform_hapmap_map(self):
        filepath = os.path.join(os.path.dirname(__file__),
            'test_files/uniform_recmap_5bp.txt')
        rec_map = parsing._load_recombination_map(filepath)
        self.assertEqual(rec_map(0), 0)
        self.assertEqual(rec_map(1), 0)
        self.assertAlmostEqual(rec_map(5), 5e-8)
        self.assertAlmostEqual(rec_map(6), 5e-8)

    def test_loading_uniform_bedgraph_map(self):
        pass

    def test_loading_heterog_hapmap_map(self):
        filepath = os.path.join(os.path.dirname(__file__),
            'test_files/heterog_recmap_5bp.txt')
        rec_map = parsing._load_recombination_map(filepath)
        self.assertEqual(rec_map(0), 0)
        self.assertEqual(rec_map(1), 0)
        self.assertAlmostEqual(rec_map(2), 1e-8)
        self.assertAlmostEqual(rec_map(3), 3e-8)
        self.assertAlmostEqual(rec_map(4), 6.5e-8)
        self.assertAlmostEqual(rec_map(5), 1e-7)
        self.assertAlmostEqual(rec_map(6), 1e-7)

    def test_uniform_map(self):
        recmap = parsing._get_uniform_recombination_map(1e7, 1.5e-8)
        self.assertAlmostEqual(
            recmap(1e7 + 1), 1e7 * utils._map_function(1.5e-8))
        

class TestReadGenotypes(unittest.TestCase):

    def test_x(self):
        pass
