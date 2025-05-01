"""
Tests for the `bootstrapping` module.
"""

import numpy as np
import unittest

from dpluspy import bootstrapping 


class TestMeansAcrossRegions(unittest.TestCase):

    def test_means_across_regions(self):
        # Deterministic tests
        pass

    def test_means_across_replicates(self):
        pass


class TestBootstrap(unittest.TestCase):

    def test_bootstrap(self):
        pass

    def test_bootstrap_with_sampling(self):
        # Construct a dataset by random sampling
        mu = np.array([1e-6, 3e-6, 2e-6])
        Sigma = np.array(
            [[6e-14, 2e-14, 1.5e-14],
             [2.5e-14, 4e-14, 2e-14],
             [2e-14, 1e-14, 7e-14]]
        )
        reps = {}
        for i in range(500):
            sample = np.random.multivariate_normal(mu, Sigma)
            sums = np.stack((sample, np.ones(3)), axis=0)
            reps[i] = {'sums': sums, 'denoms': np.array([1, 1])}



class TestMeansAcrossRegions(unittest.TestCase):

    def test_weighted_means_across_regions(self):
        pass


class TestWeightedBootstrap(unittest.TestCase):

    def test_weighted_bootstrap(self):
        pass


