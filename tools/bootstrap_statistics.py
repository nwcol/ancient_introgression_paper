"""

"""

import argparse
import numpy as np
import pickle
import warnings

from dpluspy import utils, bootstrapping


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_files', required=True, nargs='*',
        help='Pathnames of input .pkl files')
    parser.add_argument('--weighted', action='store_true',
        help='If given, use the mutation-weighted bootstrap method')
    parser.add_argument('--save_reps', action='store_true',
        help='If given, save bootstrap replicates in the output file')
    parser.add_argument('-o', '--out_file', required=True, 
        help='Output filepath')
    return parser.parse_args()


def main():
    """
    
    """
    args = get_args()
    data = {}
    for filename in args.in_files:
        regions = pickle.load(open(filename, 'rb'))
        for region in regions:
            if region in data:
                warnings.warn(f'Region {region} occurs >1 times in input')
            data[region] = regions[region]
    if args.weighted:
        ret = bootstrapping.weighted_bootstrap(data, get_reps=args.save_reps)
    else:
        ret = bootstrapping.bootstrap(data, get_reps=args.save_reps)
    if args.save_reps:
        output = {
            'means': ret[0],
            'varcovs': ret[1],
            'replicates': ret[2]
        }
    else:
        output = {
            'means': ret[0],
            'varcovs': ret[1]
        }
    example = data[next(iter(data))]
    pop_ids = example['pop_ids']
    bins = example['bins']
    output['bins'] = bins
    output['pop_ids'] = pop_ids
    with open(args.out_file, 'wb') as fout:
        pickle.dump(output, fout)

    return
    

if __name__ == '__main__':
    main()
