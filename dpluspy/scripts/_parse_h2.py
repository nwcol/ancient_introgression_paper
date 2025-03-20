import argparse
import numpy as np
import pickle

from h2py import parsing


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--vcf_file', type=str, required=True
    )
    parser.add_argument(
        '-o', '--out_file', type=str, required=True
    )
    parser.add_argument(
        '-b', '--bed_file', type=str, default=None
    )
    parser.add_argument(
        '-r', '--rec_map_file', type=str, default=None
    )
    parser.add_argument(
        '--uniform_r', type=float, default=None
    )
    parser.add_argument(
        '-m', '--mut_map_file', type=str, default=None
    )
    parser.add_argument(
        '-p', '--pop_file', type=str, default=None
    )
    parser.add_argument(
        '-R', '--region_file', type=str, default=None
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
        '--chrom_num', type=str, default=1
    )
    ## boolean options
    parser.add_argument(
        '--compute_denoms', action="store_true",
        help="Compute the H2 denominator using .bed sites"
    )
    parser.add_argument(
        "--snp_denom", action="store_true",
        help="Compute the H2 denominator using .vcf sites"
    )
    parser.add_argument(
        '--compute_two_pop', action="store_true",
        help="Compute two-population H2 for each pair of populations"
    )
    parser.add_argument(
        '--use_haplotypes', action="store_true",
        help="Compute haplotype H2, which assumes phase is known"
    )
    return parser.parse_args()


def main():

    args = get_args()

    if args.region_file is not None:
        arr = np.loadtxt(args.region_file, dtype=np.int64)
        if arr.ndim == 1:
            regions = [arr] 
        else:
            regions = arr
    else:
        regions = [None]

    r_bins = None
    bp_bins = None

    if args.r_bins is not None:
        r_bins = np.loadtxt(args.r_bins)
    elif args.bp_bins is not None:
        bp_bins = np.loadtxt(args.bp_bins)

    stats = {}
    
    for i, region in enumerate(regions):
        stats[f"{args.chrom_num}_{i}"] = parsing.parse_statistics(
            args.vcf_file,
            bed_file=args.bed_file,
            region=region,
            r_bins=r_bins,
            bp_bins=bp_bins,
            rec_map_file=args.rec_map_file,
            mut_map_file=args.mut_map_file,
            uniform_r=args.uniform_r,
            pop_file=args.pop_file,
            use_haplotypes=args.use_haplotypes,
            compute_two_pop=args.compute_two_pop,
            compute_denom=args.compute_denoms,
            snp_denom=args.snp_denom
        )

    with open(args.out_file, 'wb') as fout:
        pickle.dump(stats, fout)

    return


if __name__ == '__main__':
    main()

