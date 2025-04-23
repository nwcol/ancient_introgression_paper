
import numpy as np
import sys

from dpluspy import utils


def main():
    files = sys.argv[1:]
    tot = 0
    for file in files:
        regions, _ = utils._read_bed_file(file)
        # no_overlap = util.mask_to_regions(util.regions_to_mask(regions))
        num_sites = int(np.sum(np.diff(regions, axis=1)))
        tot += num_sites
        print(f"{file}\t{num_sites}")
    print(f"total\t{tot}")

    return


if __name__ == "__main__":
    main()
