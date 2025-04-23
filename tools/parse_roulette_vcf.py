"""
Parse site-resolution mutation maps and coverage from a Roulette VCF file.
"""

import argparse
import gzip
import numpy as np
import numpy.ma as ma
import re

from dpluspy import utils


# This the coefficient for converting roulette 'MR' rates into haploid
# per-generation mutation rates. Taken from
# https://github.com/vseplyarskiy/Roulette/tree/main/adding_mutation_rate
coeff = 1.015e-7 / 2
rate_info = 'MR'


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_file', required=True,
        help='Pathname of input Roulette VCF file')
    parser.add_argument('-np', '--np_out', required=True,
        help='Pathname of output .np file')
    parser.add_argument('-bed', '--bed_out', required=True,
        help='Pathname of output .np file')
    return parser.parse_args()


def read_seq_len(file):
    """
    Obtain the sequence length from the input VCF file header. 
    """
    seq_lens = dict()
    with gzip.open(file, 'rb') as fin:
        for lineb in fin:
            line = lineb.decode()
            if line.startswith('#'):
                if line.startswith('##contig=<ID'):
                    _, _chrom, _length, *__ = re.split('<|>|,', line)
                    chrom_num = _chrom.split('=')[1]
                    length = int(_length.split('=')[1])
                    # we only want autosomes;
                    if chrom_num.strip('chr').isnumeric():
                        seq_lens[chrom_num] = length
                else:
                    pass
            else:
                chrom_num = line.split()[0]
                break
    seq_len = seq_lens[chrom_num]
    if chrom_num.isnumeric(): 
        chrom_num = f'chr{chrom_num}'

    return chrom_num, seq_len


def read_vcf_file(vcf_file, seq_len, verbosity):
    """
    Read an array of mutation rates from a VCF file. 
    """
    rates = np.zeros(seq_len, dtype=float)
    mask = np.ones(seq_len, dtype=bool)
    counter = 0

    with gzip.open(vcf_file, 'rb') as fin:
        for lineb in fin:
            line = lineb.decode()
            if line.startswith('#'):
                continue
            split_line = line.split()
            pos1 = int(split_line[1])
            pos0 = pos1 - 1
            info = split_line[7]
            split_info = info.split(';')
            if counter == 0:
                info_names = [x.split('=')[0] for x in split_info]
                rate_idx = info_names.index(rate_info)
            rates[pos0] += float(split_info[rate_idx].split('=')[1])
            mask[pos0] = 0   
            if counter % verbosity == 0 and counter > 0:
                print(utils._current_time(), f'Parsed row {counter}')
            counter += 1

    print(utils._current_time(), f'Parsed {counter} rows')
    rates = rates * coeff
    rates[mask] = np.nan

    return rates, mask


def main():
    
    args = get_args()

    verbosity = 10000000
    chrom_num, seq_len = read_seq_len(args.in_file)
    rates, mask = read_vcf_file(args.in_file, seq_len, verbosity)
    np.save(args.np_out, rates)
    regions = utils.mask_to_regions(mask)
    utils._write_bed_file(args.bed_out, regions, chrom_num)

    return


if __name__ == '__main__':
    main()
