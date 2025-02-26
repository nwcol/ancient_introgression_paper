
import numpy as np

from h2py import parsing


"""
UNIT TESTS
"""


def test_one_genotype_h2():

    func = parsing._one_genotype_h2
    b = np.array([0, 1])
    r= np.array([0, 0])
    
    # exhaustive test of biallelic genotypes
    assert func(np.array([[1, 1], [1, 1]]), r, b) == np.array([0])

    assert func(np.array([[1, 1], [0, 1]]), r, b) == np.array([0])
    assert func(np.array([[1, 1], [1, 0]]), r, b) == np.array([0])

    assert func(np.array([[1, 1], [0, 0]]), r, b) == np.array([0])

    assert func(np.array([[0, 1], [1, 1]]), r, b) == np.array([0])
    assert func(np.array([[1, 0], [1, 1]]), r, b) == np.array([0])

    assert func(np.array([[0, 1], [0, 1]]), r, b) == np.array([1])
    assert func(np.array([[0, 1], [1, 0]]), r, b) == np.array([1])
    assert func(np.array([[1, 0], [0, 1]]), r, b) == np.array([1])
    assert func(np.array([[1, 0], [1, 0]]), r, b) == np.array([1])

    assert func(np.array([[0, 1], [0, 0]]), r, b) == np.array([0])
    assert func(np.array([[1, 0], [0, 0]]), r, b) == np.array([0])

    assert func(np.array([[0, 0], [1, 1]]), r, b) == np.array([0])

    assert func(np.array([[0, 0], [0, 1]]), r, b) == np.array([0])
    assert func(np.array([[0, 0], [1, 0]]), r, b) == np.array([0])

    assert func(np.array([[0, 0], [0, 0]]), r, b) == np.array([0])

    # some examples of triallelic genotypes
    assert func(np.array([[0, 2], [0, 0]]), r, b) == np.array([0])
    assert func(np.array([[0, 1], [2, 2]]), r, b) == np.array([0])
    assert func(np.array([[2, 1], [1, 1]]), r, b) == np.array([0])
    assert func(np.array([[1, 1], [1, 2]]), r, b) == np.array([0])

    assert func(np.array([[0, 1], [1, 2]]), r, b) == np.array([1])
    assert func(np.array([[0, 2], [0, 1]]), r, b) == np.array([1])
    assert func(np.array([[1, 2], [1, 2]]), r, b) == np.array([1])
    assert func(np.array([[1, 2], [0, 1]]), r, b) == np.array([1])

    return


def test_two_haplotype_h2():
    # iterate over a number of haplotype configurations and check their H2
    func = parsing._two_haplotype_h2
    b = np.array([0, 1])
    r = np.array([0, 0])

    # exhaustive test of biallelic haplotypes
    assert func(np.array([0, 0]), np.array([0, 0]), r, b) == np.array([0])
    assert func(np.array([0, 0]), np.array([0, 1]), r, b) == np.array([0])
    assert func(np.array([0, 0]), np.array([1, 0]), r, b) == np.array([0])
    assert func(np.array([0, 0]), np.array([1, 1]), r, b) == np.array([1])
    
    assert func(np.array([0, 1]), np.array([0, 0]), r, b) == np.array([0])
    assert func(np.array([0, 1]), np.array([0, 1]), r, b) == np.array([0])
    assert func(np.array([0, 1]), np.array([1, 0]), r, b) == np.array([1])
    assert func(np.array([0, 1]), np.array([1, 1]), r, b) == np.array([0])

    assert func(np.array([1, 0]), np.array([0, 0]), r, b) == np.array([0])
    assert func(np.array([1, 0]), np.array([0, 1]), r, b) == np.array([1])
    assert func(np.array([1, 0]), np.array([1, 0]), r, b) == np.array([0])
    assert func(np.array([1, 0]), np.array([1, 1]), r, b) == np.array([0])

    assert func(np.array([1, 1]), np.array([0, 0]), r, b) == np.array([1])
    assert func(np.array([1, 1]), np.array([0, 1]), r, b) == np.array([0])
    assert func(np.array([1, 1]), np.array([1, 0]), r, b) == np.array([0])
    assert func(np.array([1, 1]), np.array([1, 1]), r, b) == np.array([0])

    # some triallelic haplotypes
    assert func(np.array([1, 2]), np.array([0, 2]), r, b) == np.array([0])
    assert func(np.array([0, 2]), np.array([0, 1]), r, b) == np.array([0])
    assert func(np.array([2, 2]), np.array([1, 2]), r, b) == np.array([0])
    assert func(np.array([1, 2]), np.array([1, 2]), r, b) == np.array([0])

    assert func(np.array([1, 0]), np.array([2, 1]), r, b) == np.array([1])
    assert func(np.array([2, 1]), np.array([0, 0]), r, b) == np.array([1])
    assert func(np.array([1, 2]), np.array([2, 1]), r, b) == np.array([1])
    assert func(np.array([0, 1]), np.array([2, 2]), r, b) == np.array([1])

    # binning tests
    h1 = np.array([1, 1])
    h2 = np.array([0, 0])

    assert np.all(
        func(h1, h2, np.array([0, 0]), np.array([0, 1, 2, 3])) 
        == np.array([1, 0, 0])
    )
    assert np.all(
        func(h1, h2, np.array([0, 0.5]), np.array([1, 2, 3])) 
        == np.array([0, 0])
    )
    assert np.all(
        func(h1, h2, np.array([0, 1.5]), np.array([0, 1, 2, 3])) 
        == np.array([0, 1, 0])
    )
    assert np.all(
        func(h1, h2, np.array([0, 1]), np.array([0, 1, 2, 3])) 
        == np.array([0, 1, 0])
    )
    assert np.all(
        func(h1, h2, np.array([0, 2]), np.array([0, 1, 2, 3])) 
        == np.array([0, 0, 1])
    )
    assert np.all(
        func(h1, h2, np.array([0, 10]), np.array([0, 1, 2, 3])) 
        == np.array([0, 0, 0])
    )

    return


def test_two_genotype_h2():

    func = parsing._two_genotype_h2
    bins = np.array([0, 1])
    site_map = np.array([0, 0])

    # exhaustive test of biallelic genotypes
    assert func(
        np.array([[0, 0], [0, 0]]), np.array([[0, 0], [0, 0]]), site_map, bins
    ) == np.array([0])
    assert func(
        np.array([[0, 1], [0, 1]]), np.array([[0, 1], [0, 1]]), site_map, bins
    ) == np.array([0.25])
    assert func(
        np.array([[1, 1], [0, 1]]), np.array([[0, 1], [1, 1]]), site_map, bins
    ) == np.array([0.25])
    assert func(
        np.array([[1, 1], [1, 1]]), np.array([[0, 0], [0, 0]]), site_map, bins
    ) == np.array([1])
    assert func(
        np.array([[1, 1], [0, 0]]), np.array([[0, 0], [1, 1]]), site_map, bins
    ) == np.array([1])



    return


"""
SYSTEM TESTS. Some of these use example data stored in tests/test_data.
"""


def test_haplotype_h2_by_sampling():


    return


def test_genotype_h2_by_sampling():


    return







def __test_example_vcf():
    # tests `compute_statistcs` 
    vcf_file = 'test_data/example.vcf'

    sites, genotypes, sample_ids = prototype.read_genotypes(vcf_file)

    # very simple configuration
    pop_dict1 = {'p0': ['s0'], 'p1': ['s1'], 'p2': ['s2']}
    site_map = np.zeros(6)
    bins = np.array([0, 1])

    hap_stats = prototype.compute_statistics(   
        genotypes,
        sample_ids,
        site_map,
        bins,
        pop_dict=pop_dict1,
        use_haplotypes=False,
        get_two_pop=True
    )

    gen_stats = prototype.compute_statistics(   
        genotypes,
        sample_ids,
        site_map,
        bins,
        pop_dict=pop_dict1,
        use_haplotypes=True,
        get_two_pop=True
    )

    expected_hap_stats = np.array([[3, 0, 0, 6, 0, 10]])

    expected_gen_stats = np.array(
        [
            [3, 5.5, 4.375, 6, 5.75, 10],
            [3, 3.5, 3.75, 4, 3.25, 5]
        ]
    )
    assert np.all(gen_stats['sums'] == expected_gen_stats)

    return


