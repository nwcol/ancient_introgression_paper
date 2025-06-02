"""
Parse ``D+`` from a VCF file, a BED file, and a recombination map.
"""

import argparse
import pickle
from dpluspy import parsing


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', '--vcf_file', required=True, 
        help='Input VCF file')
    parser.add_argument('-bed', '--bed_file', required=True, 
        help='Input BED file')
    parser.add_argument('-pop', '--pop_file', required=True,
        help='Input population file')
    parser.add_argument('-map', '--map_file', required=True,
        help='Input recombination map file')
    parser.add_argument('--pos_col', default='Position(bp)',
        help='Recombination map file column to use as physical coordinates')
    parser.add_argument('--map_col', default='Map(cM)',
        help='Recombination map file column to use as map coordinates')
    parser.add_argument('--interp_method', default='linear',
        help='Recombination map interpolation method')
    parser.add_argument('-bins', '--bins_file', required=True,
        help='File holding recombination distance bin edges (.txt)')
    parser.add_argument('--phased', action='store_true',
        help='If given, uses the phased (haplotype) ``D+`` estimator')
    parser.add_argument('-mut', '--mut_file', help='Input mutation map file')
    parser.add_argument('--mut_col', default=None,
        help='Mutation map column (if BEDGRAPH or BED format)')
    parser.add_argument('--regions_file', required=True,
        help='Input region specification')
    parser.add_argument('-chrom', '--chrom', required=True,
        help='Prefix for region names in the output file')
    parser.add_argument('-o', '--out_file', required=True, 
        help='Output filepath (.pkl archive)')
    return parser.parse_args()


def main():
    args = get_args()
    stats = parsing.parse_statistics(
        args.vcf_file,
        args.bed_file,
        pop_file=args.pop_file,
        rec_map_file=args.map_file,
        pos_col=args.pos_col,
        map_col=args.map_col,
        interp_method=args.interp_method,
        r_bins=args.bins_file,
        phased=bool(args.phased),
        mut_map_file=args.mut_file,
        mut_map_col=args.mut_col,
        regions_file=args.regions_file,
        chrom=args.chrom
    )
    with open(args.out_file, 'wb') as fout:
        pickle.dump(stats, fout)
    return


if __name__ == '__main__':
    main()
