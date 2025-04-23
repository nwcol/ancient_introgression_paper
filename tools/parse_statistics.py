"""
Compute ``D+`` and ``H`` from sequences in a VCF file, and emit them in a .pkl
file with this structure:
{
    'chr0:0': {
        'sums': np.ndarray([[...]]),
        'denoms': np.ndarray([...]),
        'pop_ids': [...],
        'bins': np.ndarray([...])
    }
    'chr0:1': {...}
    ...
}
Where '0' and '1' are indices of genomic regions. When `--merge` is not thrown,
pairs of genomic regions '(0, 0)', '(0, 1)', ... will appear in the output file 
instead. Regions with 0 locus pairs spanning less distance than the highest bin 
edge are discarded. 

Regions should be specified in a whitespace-separated plaintext file. This file
should have as many rows as there are genomic regions, with three columns:
left/right locus start, left locus end, and right locus end. 

When a mutation map is provided, 'denoms' is replaced by 'mut_facs', which are 
used for weighting statistics when performing bootstraps.
"""

import argparse
import numpy as np
import pickle

from dpluspy import utils, parsing


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', '--vcf_file', required=True, 
        help='Input VCF file')
    parser.add_argument('-pop', '--pop_file', required=True,
        help='Input population file')
    parser.add_argument('-bed', '--bed_file', required=True, 
        help='Input BED file')
    parser.add_argument('-rec', '--rec_file', required=True,
        help='Input recombination map file')
    parser.add_argument('-bins', '--bin_file', required=True,
        help='Recombination bin file')
    parser.add_argument('-mut', '--mut_file', help='Input mutation map file')
    parser.add_argument('--regions', required=True,
        help='Input region specification')
    parser.add_argument('--merge', action='store_true',
        help='Merge within- and between-region ``D+``')
    parser.add_argument('--phased', action='store_true',
        help='Use the phased (haplotype) ``D+`` estimator')
    parser.add_argument('-chrom', '--chrom', required=True,
        help='Prefix for region names in the output file')
    parser.add_argument('-o', '--out_file', required=True, 
        help='Output filepath')
    return parser.parse_args()


def parse_merged(bins, regions, args):
    """
    Compute ``D+``, merging sums within each region with sums between that
    region and the regions to its right.
    """
    stats = {}

    for i, region in enumerate(regions):
        start, left_end, right_end = region
        within = (start, left_end)
        if left_end == right_end:
            parse_between = False
        else:
            parse_between = True
            between = ((start, left_end), (left_end, right_end))
        # compute denominators
        if args.mut_file:
            denoms, num_sites = parsing.compute_mutation_factors(
                args.mut_file,
                bed_file=args.bed_file,
                rec_map_file=args.rec_file,
                interp_method='linear',
                interval=within,
                r_bins=bins
            )
            if parse_between:
                _denoms, _num_sites = parsing.compute_mutation_factors(
                    args.mut_file,
                    bed_file=args.bed_file,
                    rec_map_file=args.rec_file,
                    interp_method='linear',
                    interval_between=between,
                    r_bins=bins
                )
                denoms += _denoms 
                num_sites += _num_sites
            denom_key = 'mut_facs'
        else:
            denoms = parsing.compute_denominators(
                bed_file=args.bed_file,
                rec_map_file=args.rec_file,
                interp_method='linear',
                interval=within,
                r_bins=bins
            )
            if parse_between:
                denoms += parsing.compute_denominators(
                    bed_file=args.bed_file,
                    rec_map_file=args.rec_file,
                    interp_method='linear',
                    interval_between=between,
                    r_bins=bins
                )
            denom_key = 'denoms'
            num_sites = denoms[-1]
        if num_sites == 0:
            print(utils._current_time(), f'Skipping empty window {i}')
            continue
        sums, pop_ids = parsing.compute_statistics(
            args.vcf_file,
            bed_file=args.bed_file,
            pop_file=args.pop_file,
            rec_map_file=args.rec_file,
            interp_method='linear',
            interval=within,
            r_bins=bins,
            phased=bool(args.phased)
        )
        if parse_between:
            sums += parsing.compute_statistics(
                args.vcf_file,
                bed_file=args.bed_file,
                pop_file=args.pop_file,
                rec_map_file=args.rec_file,
                interp_method='linear',
                interval_between=between,
                r_bins=bins,
                phased=bool(args.phased)
            )[0]
        key = f'{args.chrom}:{i}:{region}'
        stats[key] = {
            'pop_ids': pop_ids,
            'bins': bins,
            'sums': sums,
            denom_key: denoms,
            'num_sites': num_sites
        }
    
    return stats


def parse(bins, regions, args):
    """
    Compute ``D+``, keeping within- and between-region sums separate.
    """
    stats = {}

    for i, reg_i in enumerate(regions):
        for j, reg_j in enumerate(regions):
            # within
            if i == j:
                start, end = reg_i[:2]
                interval = (start, end)
                if args.mut_file is not None:
                    denoms, num_sites = parsing.compute_mutation_factors(
                        args.mut_file,
                        bed_file=args.bed_file,
                        rec_map_file=args.rec_file,
                        interp_method='linear',
                        interval=interval,
                        r_bins=bins
                    )
                    denom_key = 'mut_facs'
                else:
                    denoms = parsing.compute_denominators(
                        bed_file=args.bed_file,
                        rec_map_file=args.rec_file,
                        interp_method='linear',
                        interval=interval,
                        r_bins=bins
                    )
                    denom_key = 'denoms'
                    num_sites = denoms[-1]
                if num_sites == 0:
                    print(utils._current_time(), f'Skipping empty window {i}')
                    continue
                sums, pop_ids = parsing.compute_statistics(
                    args.vcf_file,
                    bed_file=args.bed_file,
                    pop_file=args.pop_file,
                    rec_map_file=args.rec_file,
                    interp_method='linear',
                    interval=interval,
                    r_bins=bins,
                    phased=bool(args.phased)
                )
                key = f'{args.chrom}:({i},{i}):{reg_i}'
                stats[key] = {
                    'pop_ids': pop_ids,
                    'bins': bins,
                    'sums': sums,
                    denom_key: denoms,
                    'num_sites': num_sites
                }
            # between
            elif i < j:
                left_start, left_end, right_limit = reg_i
                right_start, right_end = reg_j[:2]
                if right_start >= right_limit:
                    continue 
                if right_end >= right_limit:
                    print(utils._current_time(),
                        f'Setting right end to {right_limit} to meet limit')
                interval = ((left_start, left_end), (right_start, right_end))
                if args.mut_file is not None:
                    denoms, num_sites = parsing.compute_mutation_factors(
                        args.mut_file,
                        bed_file=args.bed_file,
                        rec_map_file=args.rec_file,
                        interp_method='linear',
                        interval_between=interval,
                        r_bins=bins
                    )
                    denom_key = 'mut_facs'
                else:
                    denoms = parsing.compute_denominators(
                        bed_file=args.bed_file,
                        rec_map_file=args.rec_file,
                        interp_method='linear',
                        interval_between=interval,
                        r_bins=bins
                    )
                    denom_key = 'denoms'
                    num_sites = denoms[-1]
                if np.sum(denoms) == 0:
                    print(utils._current_time(), 
                        f'Skipping empty windows ({i}, {j})')
                    continue
                sums, pop_ids = parsing.compute_statistics(
                    args.vcf_file,
                    bed_file=args.bed_file,
                    pop_file=args.pop_file,
                    rec_map_file=args.rec_file,
                    interp_method='linear',
                    interval_between=interval,
                    r_bins=bins,
                    phased=bool(args.phased)
                )
                key = f'{args.chrom}:({i},{j}):({reg_i}),({reg_j})'
                stats[key] = {
                    'pop_ids': pop_ids,
                    'bins': bins,
                    'sums': sums,
                    denom_key: denoms,
                    'num_sites': num_sites
                }
            else:
                continue

    return stats


def main():
    """
    Call one of the two high-level parsing functions above, depending on the 
    block architecture specified.
    """
    args = get_args()

    bins = np.loadtxt(args.bin_file)
    regions = np.loadtxt(args.regions).astype(np.int64)

    if args.merge:
        stats = parse_merged(bins, regions, args)
    else:
        stats = parse(bins, regions, args)

    with open(args.out_file, 'wb') as fout:
        pickle.dump(stats, fout)

    return  


if __name__ == '__main__':
    main()
