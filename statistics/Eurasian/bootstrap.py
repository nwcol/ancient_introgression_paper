
import argparse
import pickle
import numpy as np
from dpluspy import bootstrapping


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_files', required=True, nargs='*',
        help='Input filepaths')
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
    means, varcovs = bootstrapping.weighted_bootstrap(regions)
    template = regions[next(iter(regions))]
    bins = np.loadtxt('bins25.txt')
    stats = {'means': means, 'varcovs': varcovs, 'bins': bins, 
        'pop_ids': template['pop_ids']}
    with open(args.out_file, 'wb') as fout:
        pickle.dump(stats, fout)

    # print a summary
    num_regions = len(regions)
    num_sites = np.sum([regions[x]['denoms'][-1] for x in regions])
    num_pairs = np.sum([regions[x]['denoms'][:-1].sum() for x in regions])
    zero_sites = np.sum([1 if np.sum(regions[x]['sums']) == 0 else 0 for x in regions])

    print(f'{num_regions} regions')
    print(f'{zero_sites} regions with 0 sites')
    print(f'{int(num_sites)} total loci')
    print(f'{int(num_pairs)} total locus pairs')

    mut_sum = np.sum([regions[x]['mut_facs'][-1] for x in regions])
    mut_prod_sum = np.sum([regions[x]['mut_facs'][:-1].sum() for x in regions])
    print(f'{mut_sum / num_sites:.6} total avg u')
    print(f'{mut_prod_sum / num_pairs:.6} total locus pairs')

    return 


main()