
import numpy as np
import sys

from dpluspy import utils



def main():
    filenames = sys.argv[1:]
    tot = 0
    for filename in filenames:
        regions, _ = utils._read_bed_file(filename)
        # regions = utils._mask_to_regions(utils._regions_to_mask(regions))
        num_sites = int(np.sum(np.diff(regions, axis=1)))
        tot += num_sites
        print(f"{filename}\t{num_sites}")
    print(f"total\t{tot}")

    return


if __name__ == "__main__":
    main()
