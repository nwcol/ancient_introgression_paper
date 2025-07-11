## Write a .bed file recording the positions covered in a .vcf file.

import argparse
import gzip
import numpy as np

from h2py import utils


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
    parser.add_argument(
        "--min_GQ", type=float, default=30
    )
    parser.add_argument(
        "--skip_non_snp", type=int, default=1
    )
    parser.add_argument(
        "--pass_filters", type=str, nargs="*", default=[".", "PASS"]
    )
    return parser.parse_args()


_GQ_index_cache = {}


def read_and_write(
    vcf_file, 
    out_file, 
    min_GQ=None, 
    snp_only=True,
    verbosity=1e6,
    passing=None
):
    """
    Read and return the regions covered in a .vcf file
    """
    if not passing:
        passing = ["PASS"]

    out_func = gzip.open if out_file.endswith('gz') else open
    open_func = gzip.open if vcf_file.endswith('.gz') else open

    out_file = out_func(out_file, "wb")

    reg_start = None
    i = 0

    with open_func(vcf_file, "rb") as fin:
        for lineb in fin:
            line = lineb.decode()
            if line.startswith('#'):
                continue
                
            if i % verbosity == 0:
                print(utils.get_time(), f"parsing line {i}")
            i += 1
    
            split_line = line.split()
            ref = split_line[3]
            alts = split_line[4]
            filt = split_line[6]
            fmt = split_line[8]
            sample = split_line[9]

            if fmt in _GQ_index_cache:
                GQ_index = _GQ_index_cache[fmt]
            else:
                GQ_index = fmt.split(":").index("GQ")
                _GQ_index_cache[fmt] = GQ_index

            fails = False
            GQ_str = sample.split(":")[GQ_index]
            GQ = float(GQ_str) if GQ_str != "." else 0
            if filt not in passing:
                fails = True
            elif GQ < min_GQ:
                fails = True
            elif snp_only:
                if len(ref) > 1:
                    fails = True
                for alt in alts.split(","):
                    if len(alt) > 1 and alt != "<NON_REF>":
                        fails = True

            if fails:
                continue

            chrom = split_line[0]
            line_start = int(split_line[1]) - 1

            if reg_start is None:
                reg_start = line_start
            else:
                if line_start - last_line_end > 0 or chrom != last_chrom:
                    reg_end = last_line_end 
                    out_file.write(f"{last_chrom}\t{reg_start}\t{reg_end}\n".encode())
                    reg_start = line_start

            info = split_line[7]

            if "END" in info:
                info_pairs = [field.split("=") for field in info.split(":")]
                info_dict = {key: value for (key, value) in info_pairs}
                last_line_end = int(info_dict["END"])
            else:
                last_line_end = line_start + len(ref)

            last_chrom = chrom

    reg_start = reg_start
    reg_end = last_line_end 
    out_file.write(f"{last_chrom}\t{reg_start}\t{reg_end}".encode())
    out_file.close()

    return 


def main():
    # get positions and convert to regions; save
    args = get_args() 

    print(utils.get_time(), f'parsing coverage from {args.in_file}')
    read_and_write(
        args.in_file, 
        args.out_file, 
        verbosity=args.verbosity,
        snp_only=bool(args.skip_non_snp),
        passing=args.pass_filters,
        min_GQ=args.min_GQ
    )
    print(utils.get_time(), f'coverage written at {args.out_file}')
    return 0


if __name__ == "__main__":
    main()

