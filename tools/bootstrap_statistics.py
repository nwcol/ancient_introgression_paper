"""
Perform a bootstrap on a set of input .pkl files, each holding sums of ``D+``
for one or more genomic windows.
"""

import argparse
import pickle
import numpy as np
from dpluspy import bootstrapping


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_files', required=True, nargs='*',
        help='Input filepaths')
    parser.add_argument('--weighted', action='store_true',
        help='If given, perform a weighed bootstrap')
    parser.add_argument('--get_reps', action='store_true',
        help='If given, save bootstrap replicates in output')
    parser.add_argument('--num_reps', type=int, default=None,
        help='Number of bootstrap samplings to perform')
    parser.add_argument('-o', '--out_file', required=True, 
        help='Output filepath')
    return parser.parse_args()


def main():
    args = get_args()
    regions = {}
    for f in args.in_files:
        with open(f, 'rb') as fin:
            chrom_regions = pickle.load(fin)
        for region in chrom_regions:
            regions[region] = chrom_regions[region]
    if args.weighted:
        means, varcovs = bootstrapping.weighted_bootstrap(regions)
    else:
        means, varcovs = bootstrapping.bootstrap(regions)
    template = regions[next(iter(regions))]
    stats = {
        'means': means, 
        'varcovs': varcovs, 
        'bins': template['bins'], 
        'pop_ids': template['pop_ids']
    }
    with open(args.out_file, 'wb') as fout:
        pickle.dump(stats, fout)
    # Print a summary
    num_regions = len(regions)
    num_sites = np.sum([regions[x]['denoms'][-1] for x in regions])
    num_pairs = np.sum([regions[x]['denoms'][:-1].sum() for x in regions])
    zero_sites = np.sum([1 if np.sum(regions[x]['sums']) == 0 else 0 for x in regions])
    print(f'{num_regions} regions')
    print(f'{zero_sites} regions with 0 sites')
    print(f'{int(num_sites)} total loci')
    print(f'{int(num_pairs)} total locus pairs')
    if args.weighted:
        mut_sum = np.sum([regions[x]['mut_facs'][-1] for x in regions])
        mut_prod_sum = np.sum([regions[x]['mut_facs'][:-1].sum() for x in regions])
        print(f'{mut_sum / num_sites:.6} total avg u')
        print(f'{mut_prod_sum / num_pairs:.6} total locus pairs')

    return 


main()