
import pickle
import dpluspy


map_types = [
    'Behrer',
    'Hinch',
    'omniFIN',
    'omniLWK',
    'omniYRI',
    'pyrhoFIN',
    'pyrhoLWK',
    'pyrhoYRI',
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


labels = ["Behrer", "Hinch"]
dpluspy.plotting.plot_d_plus_curves(
    means=[stats[x]["means"] for x in labels],
    varcovs=[stats[x]["varcovs"] for x in labels],
    ax_size=2,
    cols=9,
    pop_ids=stats["Behrer"]["pop_ids"],
    bins=stats["Behrer"]["bins"],
    labels=labels,
    out="figure_1_hinch_comp",
    show=False,
    title="Behrer and Hinch maps"
)

labels = ["Behrer", "pyrhoYRI", "pyrhoLWK", "pyrhoFIN"]
dpluspy.plotting.plot_d_plus_curves(
    means=[stats[x]["means"] for x in labels],
    varcovs=[stats[x]["varcovs"] for x in labels],
    ax_size=2,
    cols=9,
    pop_ids=stats["Behrer"]["pop_ids"],
    bins=stats["Behrer"]["bins"],
    labels=labels,
    out="figure_2_pyrho_comp",
    show=False,
    title="Behrer and pyrho maps"
)

labels = ["Behrer", "omniYRI", "omniLWK", "omniFIN"]
dpluspy.plotting.plot_d_plus_curves(
    means=[stats[x]["means"] for x in labels],
    varcovs=[stats[x]["varcovs"] for x in labels],
    ax_size=2,
    cols=9,
    pop_ids=stats["Behrer"]["pop_ids"],
    bins=stats["Behrer"]["bins"],
    labels=labels,
    out="figure_3_omni_comp",
    show=False,
    title="Behrer and omni maps"
)

labels = ["Behrer", "ZhouFHS", "ZhouJHS"]
dpluspy.plotting.plot_d_plus_curves(
    means=[stats[x]["means"] for x in labels],
    varcovs=[stats[x]["varcovs"] for x in labels],
    ax_size=2,
    cols=9,
    pop_ids=stats["Behrer"]["pop_ids"],
    bins=stats["Behrer"]["bins"],
    labels=labels,
    out="figure_4_zhou_comp",
    show=False,
    title="Behrer and Zhou maps"
)


map_types = ['Behrer', 'hapmap']
stats = {}
for map_type in map_types:
    fname = f"../stats/{map_type}_skip6_7_stats.pkl"
    with open(fname, "rb") as fin:
        _stats = pickle.load(fin)
        to_pops = _stats["pop_ids"]
        to_pops.pop(-1)
        to_pops.pop(-1)
        stats[map_type] = dpluspy.bootstrapping.subset_stats(
            _stats, to_pops=to_pops, return_dict=True)
        
labels = ["Behrer", "hapmap"]
dpluspy.plotting.plot_d_plus_curves(
    means=[stats[x]["means"] for x in labels],
    varcovs=[stats[x]["varcovs"] for x in labels],
    ax_size=2,
    cols=9,
    pop_ids=stats["Behrer"]["pop_ids"],
    bins=stats["Behrer"]["bins"],
    labels=labels,
    out="figure_5_hapmap_comp",
    show=False,
    title="Behrer and Hapmap maps: excluding chromosomes 6, 7"
)