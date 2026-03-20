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
    parser.add_argument('-i', '--infile', required=True,
        help='Pathname of input Roulette VCF file')
    parser.add_argument('-o', '--outfile', required=True,
        help='Pathname of output .npy file')
    return parser.parse_args()


def read_seq_length(file):
    """
    Obtain the sequence length from the input VCF file header. 
    """
    seq_lengths = dict()
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
                        seq_lengths[chrom_num] = length
                else:
                    pass
            else:
                chrom_num = line.split()[0]
                break
    seq_length = seq_lengths[chrom_num]
    if chrom_num.isnumeric(): 
        chrom_num = f'chr{chrom_num}'
    return chrom_num, seq_length


def read_vcf_file(vcf_file, seq_length, verbosity):
    """
    Read an array of mutation rates from a VCF file. 
    """
    rates = np.zeros(seq_length, dtype=float)
    mask = np.ones(seq_length, dtype=bool)
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
    return rates


def parse_roulette_vcf(infile, outfile):
    verbosity = 10000000
    chrom_num, seq_length = read_seq_length(infile)
    rates = read_vcf_file(infile, seq_length, verbosity)
    np.save(outfile, rates)
    return


def main():
    args = get_args()
    parse_roulette_vcf(args.infile, args.outfile)
    return


if __name__ == '__main__':
    main()
