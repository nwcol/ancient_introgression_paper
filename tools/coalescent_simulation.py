## run msprime ancestry, mutation simulations and save resulting .vcf file

import argparse
import demes
import msprime
import numpy as np

from h2py import util, simulation


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-g", "--graph_file", 
        required=True
    )
    parser.add_argument(
        "-o", "--out_file", required=True
    )
    parser.add_argument(
        "-p", "--pop_ids", 
        nargs='*', 
        default=None
    )
    parser.add_argument(
        '-n', '--num_samples',
        type=int,
        default=1
    )
    parser.add_argument(
        "-u", "--mut_rate", 
        type=float,
        default=1.5e-8
    )
    parser.add_argument(
        '--mut_map', type=str, default=None
    )
    parser.add_argument(
        '--bed_file', type=str, default=None
    )
    parser.add_argument(
        "-r", "--rec_rate", 
        type=float,
        default=1e-8
    )
    parser.add_argument(
        '--rec_map', type=str, default=None
    )
    parser.add_argument(
        '-L', '--seq_length', type=int, default=None
    )
    parser.add_argument(
        '--dtwf_time',
        type=int,
        default=1000
    )
    return parser.parse_args()


def load_rec_map(rec_map_file, seq_length):
    ##
    edges, map_vals = util.read_hapmap_rec_map(rec_map_file)
    map_rates = np.diff(map_vals) / np.diff(edges)  # cM/bp
    map_rates /= 100
    edges[0] = 0

    if edges[-1] < seq_length:
        edges[-1] = seq_length
    else:
        map_rates = map_rates[edges < seq_length]
        edges = np.append(edges[edges < seq_length], seq_length)

    discret = 50000
    edges = np.concatenate((np.arange(0, seq_length, discret), [seq_length]))
    map_coords = util.read_recombination_map(rec_map_file, edges)
    map_rates = np.diff(map_coords) / np.diff(edges) / 100
    print(np.mean(map_rates), len(map_rates))
    rec_map = msprime.RateMap(position=edges, rate=map_rates)
    print(util.get_time(), 'loaded rec-map')

    return rec_map


def load_mut_map(mut_map_file, seq_length, mask_file=None):
    ## load a mutation .bedgraph file.

    if mut_map_file.endswith(".npy"):
        print(util.get_time(), "discretizing mutation map")
        raw_muts = np.load(mut_map_file)
        discret = 50000
        edges = np.concatenate((np.arange(0, seq_length, discret), [seq_length]))
        windows = np.stack((edges[:-1], edges[1:]), axis=1)

        if mask_file is not None:
            mask = util.regions_to_mask(util.read_bedfile(mask_file))
        else:
            mask = np.zeros(len(raw_muts))

        muts = simulation.discretize_mut_map(raw_muts, mask, windows)

    else:
        regions, data = util.read_bedgraph(mut_map_file)
        muts = data['mut_rate']

        coords = regions[:, 0]
        coords[0] = 0

        if coords[-1] < seq_length:
            edges = np.append(coords, seq_length)
        else:
            muts = muts[coords < seq_length]
            edges = np.append(coords[coords < seq_length], seq_length)

    mut_map = msprime.RateMap(position=edges, rate=muts)
    print(util.get_time(), 'loaded mut-map')

    return mut_map


def increment1(x):
    """
    Increment 1 to every .vcf site to make it 1-indexed rather than 0-indexed.
    """
    return [_ + 1 for _ in x]


def main():
    """
    Run a coalescent simulation.
    """
    args = get_args()

    if args.seq_length is not None:
        seq_length = args.seq_length
    elif args.bed_file is not None:
        seq_length = util.read_bedfile(args.bed_file)[-1, 1] + 1
    else:
        raise ValueError("sequence length argument required")

    graph = demes.load(args.graph_file)
    demography = msprime.Demography.from_demes(graph)

    if args.pop_ids is None:
        sampled_demes = [d.name for d in graph.demes if d.end_time == 0]
    else:
        sampled_demes = args.pop_ids

    config = {s: args.num_samples for s in sampled_demes}

    if args.rec_map is None:
        rec_rate = args.rec_rate
    else:
        # rec_rate = load_rec_map(args.rec_map, seq_length)
        rec_rate = msprime.RateMap.read_hapmap(
            args.rec_map, 
            sequence_length=seq_length,
            position_col=0,
            map_col=2   
        )

    if args.mut_map is None:
        mut_rate = args.mut_rate
    else:
        mut_rate = load_mut_map(args.mut_map, seq_length, args.bed_file)

    if args.dtwf_time > 0:
        model = [
            msprime.DiscreteTimeWrightFisher(duration=args.dtwf_time),
            msprime.StandardCoalescent(),
        ]
    else:
        model=msprime.StandardCoalescent()

    print(util.get_time(), "simulating ancestry")
    ts = msprime.sim_ancestry(
        samples=config,
        ploidy=2,
        demography=demography,
        recombination_rate=rec_rate,
        discrete_genome=True,
        record_provenance=False,
        sequence_length=seq_length,
        model=model
    )
    print(util.get_time(), "completed ancestry simulation")
    mts = msprime.sim_mutations(
        ts,
        rate=mut_rate,
        record_provenance=False
    )
    print(util.get_time(), "completed mutation simulation")

    sample_names = []
    for deme in sampled_demes:
        if args.num_samples == 1:
            sample_names.append(deme)
        else:
            for n in range(args.num_samples):
                sample_names.append(f'{deme}_{n}')

    with open(args.out_file, 'w') as file:
        mts.write_vcf(
            file,
            individual_names=sample_names,
            contig_id='1',
            position_transform=increment1
        )
    return 0


if __name__ == '__main__':
    main()

