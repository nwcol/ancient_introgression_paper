## Trim by bin

import numpy as np
import pickle

from dpluspy import bootstrapping


stats = pickle.load(open('Eurasian_stats_weighted.pkl', 'rb'))
stats_trimmed = bootstrapping.subset_statistics(
    stats, min_r=1e-6, max_r=1e-2, return_dict=True)
with open('Eurasian_weighted_trimmed.pkl', 'wb') as fout:
    pickle.dump(stats_trimmed, fout)

