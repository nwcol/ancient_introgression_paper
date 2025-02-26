## create bedfiles containing 1kb windows that fall in each B bin
import numpy as np
from h2py import util


def main():

    file = "B_map_func_YRI_chr1_iter5_3682852.csv.gz"
    windows, data = util.read_bedgraph(file)

    b = data['B']
    b_bins = np.array([0.50, 0.70, 0.80, 0.90, 0.95, 1.0])

    for i in range(len(b_bins) - 1):
        lower, upper = b_bins[i], b_bins[i + 1]
        where = np.where((b >= lower) & (b < upper))[0]
        where_windows = windows[where]
        merged_windows = util.collapse_regions(where_windows)
        util.write_bedfile(f"b_bin{i}_chr1.bed.gz", merged_windows, chrom_num="chr1")

    return


if __name__ == "__main__":
    main()

