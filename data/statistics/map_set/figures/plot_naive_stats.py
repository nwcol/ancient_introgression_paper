
import pickle
import dpluspy


map_types = [
    'Behrer',
    'omniYRI',
    'ZhouFHS',
    'ZhouJHS'
]
stats = {}
for map_type in map_types:
    if map_type == "hapmap":
        continue
    fname = f"../stats/{map_type}_stats.pkl"
    with open(fname, "rb") as fin:
        _stats = pickle.load(fin)
        to_pops = _stats["pop_ids"][:-2]
        stats[map_type] = dpluspy.bootstrapping.subset_stats(
            _stats, to_pops=to_pops, return_dict=True)

naive_stats = {}
for map_type in map_types:
    if map_type == "hapmap":
        continue
    fname = f"../stats/{map_type}_naive_stats.pkl"
    with open(fname, "rb") as fin:
        _stats = pickle.load(fin)
        to_pops = _stats["pop_ids"][:-2]
        naive_stats[map_type] = dpluspy.bootstrapping.subset_stats(
            _stats, to_pops=to_pops, return_dict=True)


labels = ["weighted", "naive"]


for i, map_type in enumerate(map_types):
    dpluspy.plotting.plot_d_plus_curves(
        means=[stats[map_type]["means"], naive_stats[map_type]["means"]],
        varcovs=[stats[map_type]["varcovs"], naive_stats[map_type]["varcovs"]],
        ax_size=2,
        cols=9,
        pop_ids=stats["Behrer"]["pop_ids"],
        bins=stats["Behrer"]["bins"],
        labels=labels,
        out=f"figure_{6 + i}_{map_type}_weighting",
        show=False,
        title=f"Weighting: {map_type}"
    )
