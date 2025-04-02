
## sampling functions for numerically computing H2


def _sample_haplotype_h2(haplotypes, num_reps=1000):
    # from an array of haplotypes with shape (2, n), sample haplotypes to 
    # compute expected H2 numerically. for numerical validation
    num_sites, num_haps = haplotypes.shape
    assert num_sites == 2

    indices = np.arange(num_haps)
    sum_H2 = 0
    for i in range(num_reps):
        sample = haplotypes[:, np.random.choice(indices, size=2, replace=False)]
        sum_H2 += (
            (sample[0, 0] != sample[0, 1]) 
            & (sample[1, 0] != sample[1, 1])
        )
    
    mean_H2 = sum_H2 / num_reps
    return mean_H2


def _sample_two_pop_haplotype_h2(haplotypes1, haplotypes2, num_reps=1000):
    # for numerical validation
    num_sites1, num_haps1 = haplotypes1.shape
    num_sites2, num_haps2 = haplotypes2.shape
    assert num_sites1 == num_sites2 == 2

    sum_H2 = 0
    for i in range(num_reps):
        sample1 = haplotypes1[:, np.random.randint(num_haps1)]
        sample2 = haplotypes2[:, np.random.randint(num_haps2)]
        sum_H2 += ((sample1[0] != sample2[0]) & (sample1[1] != sample2[1]))
    
    mean_H2 = sum_H2 / num_reps
    return mean_H2


def _sample_genotype_h2(genotypes, num_reps=1000, between=False):
    #
    num_sites, num_samps, _ = genotypes.shape
    assert num_sites == 2

    # precompute H2 for each sample
    sample_H2s = []
    for i in range(num_samps):
        genotype = genotypes[:, i]
        sample_H2s.append(
            (genotype[0, 0] != genotype[0, 1]) 
            & (genotype[1, 0] != genotype[1, 1]) 
        )
    sample_H2s = np.array(sample_H2s)

    sum_H2 = 0

    # allow haplotypes to be sampled from different genomes
    if between:
        for i in range(num_reps):
            index1 = np.random.randint(num_samps)
            index2 = np.random.randint(num_samps)

            if index1 == index2:
                sum_H2 += sample_H2s[index1]
            else:
                sample1 = genotypes[:, index1]
                sample2 = genotypes[:, index2]
                hap1 = sample1[[0, 1], np.random.randint(2, size=2)]
                hap2 = sample2[[0, 1], np.random.randint(2, size=2)]   
                sum_H2 += ((hap1[0] != hap2[0]) & (hap1[1] != hap2[1]))    

    # average over within-sample H2
    else:
        for i in range(num_reps):
            sum_H2 += sample_H2s[np.random.randint(num_samps)]

    mean_H2 = sum_H2 / num_reps
    return mean_H2


def _sample_two_pop_genotype_h2(genotypes1, genotypes2, num_reps=1000):
    # for numerical validation. two-sample genotype H2 is the simpler case!
    num_sites1, num_samps1, _ = genotypes1.shape
    num_sites2, num_samps2, _ = genotypes2.shape
    assert num_sites1 == num_sites2 == 2

    sum_H2 = 0
    for i in range(num_reps):
        sample1 = genotypes1[:, np.random.randint(num_samps1)]
        sample2 = genotypes2[:, np.random.randint(num_samps2)]
        hap1 = sample1[[0, 1], np.random.randint(2, size=2)]
        hap2 = sample2[[0, 1], np.random.randint(2, size=2)]
        sum_H2 += ((hap1[0] != hap2[0]) & (hap1[1] != hap2[1]))
    
    mean_H2 = sum_H2 / num_reps
    return mean_H2