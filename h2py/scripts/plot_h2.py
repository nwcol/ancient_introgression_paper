"""
Plot an arbitrary number of H, H2 expectations alongside 0 or 1 empirical vals.
"""
import argparse
import demes
import matplotlib.pyplot as plt
import numpy as np
import pickle
import os

from h2py import h2_parsing, inference, plotting


def get_args():
    # get args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d', '--data_files', 
        nargs='*',
        default=[],
        type=str
    )
    parser.add_argument(
        '-g', '--graph_files', 
        nargs='*',
        default=[],
        type=str
    )
    parser.add_argument(
        '--pop_ids',
        nargs='*',
        default=None,
        type=str,
    )
    parser.add_argument(
        '--two_sample', type=int, default=0
    )
    parser.add_argument(
        '--compute_ll', type=int, default=1
    )
    parser.add_argument(
        '--confidence', type=float, default=0.95
    )
    parser.add_argument(
        '--ylim', type=str, default=None
    )
    parser.add_argument(
        '--share_y', type=int, default=0
    )
    parser.add_argument(
        '--ratios', type=int, default=0
    )
    parser.add_argument(
        '-H', '--plot_H', type=int, default=1
    )
    parser.add_argument(
        '-o', '--out_file', required=True
    )
    parser.add_argument(
        '-u', '--u', type=float, default=None
    )
    parser.add_argument(
        '--labels', nargs='*', default=None
    )
    parser.add_argument(
        '--phased', type=int, default=0
    )
    parser.add_argument(
        "--fill_between", action="store_true"
    )
    parser.add_argument(
        "--title", default=None
    )
    parser.add_argument(
        "--hlines", default=None, type=float, nargs="*"
    )
    parser.add_argument(
        "--n_cols", default=5, type=int
    )
    return parser.parse_args()


def main():
    """
    """
    args = get_args()

    if len(args.graph_files) > 0:
        if args.u is not None:
            u = args.u
        else:
            u = None
            for file in args.graph_files:
                g = demes.load(file)
                if 'opt_info' in g.metadata:
                    u = g.metadata['opt_info']['u']
                elif 'u' in g.metadata:
                    u = g.metadata['u']
                break
            if u is None:
                raise ValueError('please provide a u parameter')
    
    # load statistics
    g = args.graph_files[0] if len(args.graph_files) > 0 else None
    datas = []
    labels = []

    for file in args.data_files:
        with open(file, 'rb') as fin:
            dic = pickle.load(fin)
        for key in dic:
            if args.pop_ids is None:
                data = h2_parsing.subset_statistics(dic[key], graph=g)
            else:
                data = h2_parsing.subset_statistics(dic[key], to_pops=args.pop_ids)
            label = os.path.basename(file) + '-' + key

            if args.ratios:
                data["means"][:-1] /= (data["means"][-1] ** 2)
                data["covs"][:-1] /= (data["means"][-1] ** 4)

            datas.append(data)
            labels.append(label)

    # load graphs
    if len(datas) > 0:
        if args.pop_ids is None:
            data = datas[0]
            models = [inference.moments_H2(g, u=u, data=data, phased=args.phased) 
                      for g in args.graph_files]
        else:
            pop_ids = args.pop_ids
            bins = datas[0]['bins']
            models = [inference.moments_H2(g, u=u, bins=bins, sampled_demes=pop_ids, phased=args.phased) 
                      for g in args.graph_files]
        if args.compute_ll:
            lls = [np.round(inference.compute_ll(m, data, include_H=False), 2)
                   for m in models] 
            labels += [os.path.basename(args.graph_files[i]) + f', ll={lls[i]}'
                       for i in range(len(args.graph_files))]
        else:
            labels += [os.path.basename(args.graph_files[i]) 
                       for i in range(len(args.graph_files))]

    else:
        pop_ids = args.pop_ids
        models = [inference.moments_H2(g, u=u, sampled_demes=pop_ids) 
                  for g in args.graph_files]
        labels += [os.path.basename(g) for g in args.graph_files]

    if args.labels is not None:
        labels = args.labels

    if args.ylim:
        if args.ylim == "0,":
            ylim = (0, None)
        else:
            ylim = [float(x) for x in args.ylim.split(",")]
    else:
        ylim = None
        
    plotting.plot_H2s(
        models=models, 
        datas=datas,
        labels=labels,
        plot_H=args.plot_H,
        ylim=ylim,
        share_y=args.share_y,
        conf=args.confidence,
        fill=bool(args.fill_between),
        title=args.title,
        hlines=args.hlines,
        n_cols=args.n_cols
    )

    plt.savefig(args.out_file, dpi=244, bbox_inches='tight')
    return 


if __name__ == "__main__":    
    main()
