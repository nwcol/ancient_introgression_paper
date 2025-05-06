"""

"""

import argparse
import numpy as np
import sys
import pandas


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-cyto', '--cyto_band_file', required=True, 
        help='Input band file')
    parser.add_argument('-chrom', '--chrom_num', required=True, 
        help='Chromosome number, prefixed by `chr`')
    parser.add_argument('-size', '--window_size', required=True, type=int,
        help='Desired window size in Mb')
    parser.add_argument('-flank', '--flank', required=True, type=int,
        help='Length of flank around centromeres and chromosome ends in Mb')

    return parser.parse_args()


def _main():

    args = get_args()
    full_df = pandas.read_csv(args.cyto_band_file, sep='\t', 
        names=['chrom', 'start', 'end', 'name', 'stain'])
    df = full_df[full_df['chrom'] == args.chrom_num]
    

    return


def main():
    ## expects chrXX
    size = 5000000
    num = str(sys.argv[1])
    bands = []
    cens = []
    with open("cytoBand.txt", "r") as fin:
        for line in fin:
            splitline = line.split()
            if splitline[0] == num:
                band = [int(splitline[1]), int(splitline[2])]
                bands.append(band)
                if splitline[4] == "acen":
                    cens.append(band)
    start = np.min(bands)
    end = np.max(bands) 
    if len(cens) > 0:
        censtart = np.min(cens)
        cenend = np.max(cens)
        starts0 = np.arange(start, censtart, size)
        starts1 = np.arange(cenend, end, size)
        starts = np.concatenate((starts0, starts1))
        ends0 = starts0 + size
        ends1 = starts1 + size
        ends0[-1] = censtart
        ends1[-1] = end
        ends = np.concatenate((ends0, ends1))
        r_ends = np.array([censtart] * len(starts0) + [end] * len(starts1))
    else:
        starts = np.arange(start, end, size)
        ends = starts + size
        r_ends = [ends[-1]] * len(ends)

    with open(f"windows/windows_{num}.txt", "w") as fout:
        for start, end, r_end in zip(starts, ends, r_ends):
            fout.write(f"{start}\t{end}\t{r_end}\n")

    return 


main()