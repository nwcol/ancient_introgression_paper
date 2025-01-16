## Write a .bed file recording the positions covered in a .vcf file.

import argparse
import gzip
import numpy as np
import warnings

from h2py import util


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--in_file',
        type=str, 
        required=True
    )
    parser.add_argument(
        '-o', '--out_file',
        type=str,
        required=True
    )
    parser.add_argument(
        "--verbosity", type=int, default=1000000
    )
    return parser.parse_args()


def read_and_write(vcf_file, out_file, min_GQ=None, verbosity=1e6):
    """
    Read and return the regions covered in a .vcf file
    """
    out_func = gzip.open if out_file.endswith('gz') else open
    open_func = gzip.open if vcf_file.endswith('.gz') else open

    out_file = out_func(out_file, "wb")
    i = 0

    with open_func(vcf_file, "rb") as fin:
        for lineb in fin:
            line = lineb.decode()
            if line.startswith('#'):
                continue
            
            split_line = line.split()
            chrom = split_line[0]
            position = int(split_line[1])
            info = split_line[7]
            ref = split_line[3]
            
            if i == 0:
                start = position
            else:
                if position - last_position > 1 or chrom != last_chrom:
                    end = last_position 
                    # convert start position to 0 index
                    start0 = start - 1
                    out_file.write(f"{last_chrom}\t{start0}\t{end}\n".encode())
                    start = position

            if min_GQ is not None:


            if "END" in info:
                infos = [field.split("=") for field in info.split(":")]
                info_dict = {key: value for (key, value) in infos}
                last_position = int(info_dict["END"])
            elif len(ref) > 1:
                last_position = position + len(ref) - 1
            else:
                last_position = position

            last_chrom = chrom
            i += 1

            if i % verbosity == 0:
                print(util.get_time(), f"parsed coverage for {i} sites")

    end = last_position 
    start0 = start - 1
    out_file.write(f"{chrom}\t{start0}\t{end}".encode())
    out_file.close()

    return 


def main():
    # get positions and convert to regions; save
    args = get_args() 

    print(util.get_time(), f'parsing coverage from {args.in_file}')
    read_and_write(args.in_file, args.out_file, verbosity=args.verbosity)
    print(util.get_time(), f'coverage written at {args.out_file}')
    return 0


if __name__ == "__main__":
    main()

