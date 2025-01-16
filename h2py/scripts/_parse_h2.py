
import argparse
import numpy as np
import pickle

from h2py import util, h2_parsing, prototype


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--vcf_file',
        type=str,
        required=True
    )
    parser.add_argument(
        '-b', '--bed_file',
        type=str
    )
    parser.add_argument(
        '-r', '--rec_map_file', type=str, default=None
    )
    parser.add_argument(
        '-m', '--mut_map_file', type=str, default=None
    )
    parser.add_argument(
        '-p', '--pop_file',
        type=str
    )
    parser.add_argument(
        '-R', '--region_file',
        type=str
    )
    parser.add_argument(
        '-o', '--out_file',
        type=str,
        required=True
    )
    parser.add_argument(
        '--r_bins', type=str, default=None
    )
    parser.add_argument(
        '--bp_bins', type=str, default=None
    )
    parser.add_argument(
        '--min_reg_len', type=int, default=0
    )
    parser.add_argument(
        '--compute_denom', type=int, default=1
    )
    parser.add_argument(
        '--compute_snp_denom', type=int, default=0
    )
    parser.add_argument(
        '--compute_two_sample', type=int, default=1
    )
    parser.add_argument(
        '--chrom', type=str, default=1
    )
    parser.add_argument(
        '--verbose', type=str,  default=1
    )
    parser.add_argument(
        '--uniform_r', type=float, default=None
    )
    parser.add_argument(
        '--missing_to_ref', type=int, default=1
    )
    parser.add_argument(
        '--use_haplotypes', type=int, default=0
    )
    return parser.parse_args()


def main():

    args = get_args()

    if args.region_file is not None:
        arr = np.loadtxt(args.region_file, usecols=(0,1,2))
        if arr.ndim == 1:
            region = arr
            regions = None 
        else:
            regions = arr
            region = None

    if args.r_bins is not None:
        r_bins = np.loadtxt(args.r_bins)
        bp_bins = None
    else:
        r_bins = None
        bp_bins = None
    
    stats = prototype.parse_statistics(
        args.vcf_file,
        region=region,
        regions=regions,
        chrom_num=args.chrom,
        rec_map_file=args.rec_map_file,
        uniform_r=args.uniform_r,
        r_bins=r_bins,
        bp_bins=bp_bins,
        bed_file=args.bed_file,
        pop_file=args.pop_file,
        use_haplotypes=args.use_haplotypes,
        get_two_pop=args.compute_two_sample,
        compute_denom=args.compute_denom,
        mut_map_file=None
    )

    with open(args.out_file, 'wb') as fout:
        pickle.dump(stats, fout)

    return


def _main():
    """
    
    """
    args = get_args()

    if args.r_bins is not None:
        r_bins = np.loadtxt(args.r_bins)
        bp_bins = None
    elif args.bp_bins is not None:
        r_bins = None
        bp_bins = np.loadtxt(args.bp_bins).astype(np.int64)
    else:
        r_bins = None
        bp_bins = None

    if args.region_file is not None:
        regions = np.loadtxt(args.region_file, usecols=(0,1,2))
        if regions.ndim == 1:
            regions = regions[np.newaxis]
    else:
        regions = [None]

    if args.verbose:
        print(util.get_time(), f'computing H2 on chromosome {args.chrom}')

    region_stats = {}
    for i, region in enumerate(regions):
        key = f'{args.chrom}_{i}'
        region_stats[key] = h2_parsing.compute_H2(
            args.vcf_file,
            bed_file=args.bed_file,
            rec_map_file=args.rec_map_file,
            mut_map_file=args.mut_map_file,
            r=args.uniform_r,
            pop_file=args.pop_file,
            region=region,
            r_bins=r_bins,
            bp_bins=bp_bins,
            phased=False,
            compute_denom=args.compute_denom,
            compute_snp_denom=args.compute_snp_denom,
            compute_two_sample=args.compute_two_sample,
            missing_to_ref=args.missing_to_ref,
            verbose=args.verbose
        )

    with open(args.out_file, 'wb') as fout:
        pickle.dump(region_stats, fout)

    return 


if __name__ == '__main__':
    main()

