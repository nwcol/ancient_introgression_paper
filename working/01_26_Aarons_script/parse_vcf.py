import sys
import numpy as np
import argparse
import gzip
import pickle
from datetime import datetime


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    sys.stderr.flush()


def current_time():
    return " [" + datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S") + "]"


def make_parser():
    ADHF = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser("parse_vcf.py", formatter_class=ADHF)
    parser.add_argument("-r", "--recombination_map", required=True, type=str)
    parser.add_argument("-f", "--vcf_file", required=True, type=str)
    parser.add_argument("-m", "--mask_file", required=True, type=str)
    parser.add_argument("-o", "--outfile", required=True, type=str)
    parser.add_argument("--start", required=True, type=int)
    parser.add_argument("--end", required=True, type=int)
    parser.add_argument("--name", required=True, type=str)
    return parser


def get_positions(mask, start, end):
    # bed file with half open interval [start, end)
    _L = 0
    with gzip.open(mask, "rb") as fin:
        for line in fin:
            _, left, right = line.decode().split()
            _L += int(right) - int(left)

    positions = np.empty(_L, dtype=int)
    idx = 0
    with gzip.open(mask, "rb") as fin:
        for line in fin:
            _, left, right = line.decode().split()
            left = int(left)
            right = int(right)
            positions[idx : idx + right - left] = np.arange(left, right)
            idx += right - left
    assert not np.any(positions == 0)
    assert positions[-1] != 0
    return positions[(positions >= start) & (positions < end)]


def get_recombination_values(positions, rec_map, col="Map(cM)", cM=True):
    # load the recombination map
    try:
        fin = gzip.open(rec_map, "rb")
        header = fin.readline()
    except:
        fin = open(rec_map, "rb")
        header = fin.readline()
    fin.close()

    # get the columns of data that we need (positions and map)
    cols = header.decode().split()
    if "Position(bp)" not in cols:
        raise ValueError("Need a column named 'Position(bp)'")
    pos_col = cols.index("Position(bp)")
    if col in cols:
        map_col = cols.index(col)
    else:
        raise ValueError(f"Column '{col}' not in recombination map file")

    data = np.loadtxt(rec_map, skiprows=1, usecols=(pos_col, map_col))

    # interpolate to get r-values at each position
    if positions[0] < data[0, 0] or positions[-1] > data[-1, 0]:
        eprint(current_time(), "Warning: positions are outside of map")
    rvals = np.interp(
        positions, data[:, 0], data[:, 1], left=data[0, 1], right=data[-1, 1]
    )

    # if our map is in centiMorgans, divide by 100
    #if cM:
    #    rvals /= 100
    return rvals


def count_pairs(rvals, bins):
    num_pairs = np.zeros(len(bins) - 1, dtype=int)
    tot_rvals = len(rvals)
    for i, r in enumerate(rvals):
        cum_counts = np.searchsorted(rvals[i + 1 :], bins[1:] + rvals[i])
        # very clever, Nick!
        num_pairs[0] += cum_counts[0]
        num_pairs[1:] += np.diff(cum_counts)
        if i % 1000000 == 0:
            eprint(current_time(), f"at position {i} of {tot_rvals}")
    return num_pairs


def get_vcf_positions_genotypes(vcf, mask, start, end):
    # just from first data column
    data = np.loadtxt(mask, usecols=(1,2))
    left = data[:, 0]
    right = data[:, 1]
    vcf_positions = []
    gts = []
    with gzip.open(vcf, "rb") as fin:
        for line in fin:
            l = line.decode()
            if l.startswith("#"):
                continue
            pos = int(l.split()[1])
            if np.any(np.logical_and(pos >= left, pos < right)):
                # vcf_positions.append(pos)
                gt = l.split()[9].split(":")[0]
                if gt[0] == gt[2]:
                    # hom
                    continue
                    # gts.append(0)
                else:
                    # het
                    vcf_positions.append(pos)
                    gts.append(1)
    vcf_positions, gts = np.array(vcf_positions), np.array(gts)
    gts = gts[(vcf_positions >= start) & (vcf_positions < end)]
    vcf_positions = vcf_positions[(vcf_positions >= start) & (vcf_positions < end)]
    return vcf_positions, gts


def map_function(r):
    """
    Haldane's map function; transforms distance in r to cM.
    """
    return -50 * np.log(1 - 2 * r)


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    mask = args.mask_file
    vcf = args.vcf_file
    rec_map = args.recombination_map
    fname = args.outfile
    start = args.start 
    end = args.end
    name = args.name

    eprint(
        current_time(), f"Parsing VCF: {args.vcf_file}, using {args.recombination_map}"
    )

    positions = get_positions(mask, start, end)
    eprint(current_time(), "Got positions from mask")

    rvals = get_recombination_values(positions, rec_map)
    eprint(current_time(), "Got recombination map values for all positions")

    vcf_positions, gts = get_vcf_positions_genotypes(vcf, mask, start, end)
    eprint(current_time(), "Retrived VCF positions and genotypes")

    _bins = np.concatenate(([0], np.logspace(-6, -1, 16), [2.1544e-1, 4.9999e-1]))
    bins = map_function(_bins)

    denoms = count_pairs(rvals, bins)
    denoms = np.append(denoms, 1)
    eprint(current_time(), "Counted pairs")


    vcf_rvals = get_recombination_values(vcf_positions, rec_map)
    sums = count_pairs(vcf_rvals, bins)
    sums = np.append(sums, 1)
    eprint(current_time(), "Got H2 from jointly heterozygous sites")

    data = {name: 
        {
            "pops": ["X"],
            "bins": _bins, 
            "denoms": denoms, 
            "sums": sums[:, None]
        }
    }
    print(data)

    with open(fname, "wb+") as fout:
        pickle.dump(data, fout)
    eprint(current_time(), "Saved data")
