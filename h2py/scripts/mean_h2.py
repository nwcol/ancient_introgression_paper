
import argparse
import numpy as np
import pickle

from h2py import h2_parsing


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--in_files',
        nargs='*',
        type=str,
        required=True
    )
    parser.add_argument(
        '-o', '--out_file',
        type=str,
        required=True
    )
    parser.add_argument(
        '-d', '--denom_files',
        nargs='*',
        default=None,
        type=str
    )
    parser.add_argument(
        '--replicates', action="store_true"
    )
    parser.add_argument(
        '--weighted', action="store_true"
    )   
    parser.add_argument(
        "--normalize_to",
        type=float,
        default=None
    )
    return parser.parse_args()


def main():
    """
    
    """
    args = get_args()

    if args.replicates:
        data = []
        for in_file in args.in_files:
            with open(in_file, 'rb') as fin:
                file_data = pickle.load(fin)
            data.append(file_data)
    else:
        data = {}
        for in_file in args.in_files:
            with open(in_file, 'rb') as fin:
                file_data = pickle.load(fin)
            for key in file_data:
                if key in data:
                    raise ValueError(f'{key} occurs twice in input')
                else:
                    data[key] = file_data[key]

    if args.denom_files is not None:
        denoms = {}
        for in_file in args.denom_files:
            with open(in_file, 'rb') as fin:
                file_data = pickle.load(fin)
            for key in file_data:
                denoms[key] = file_data[key] 
    else:
        denoms = None

    mean_stats = h2_parsing.compute_mean_H2(
        data, 
        denominators=denoms,
        weighted=bool(args.weighted),
        normalize_to=args.normalize_to
    )
    out = {'mean': mean_stats}

    with open(args.out_file, 'wb') as fout:
        pickle.dump(out, fout)

    return 


if __name__ == '__main__':
    main()

