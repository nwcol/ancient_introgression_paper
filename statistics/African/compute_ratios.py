## Compute and print out the ratios D+ / H^2

import numpy as np
import pickle

from dpluspy import utils


stats = pickle.load(open('African_stats_weighted.pkl', 'rb'))
pop_ids = stats['pop_ids']
Dp_names = utils._get_Dplus_names(pop_ids)
H_names = utils._get_H_names(pop_ids)
ratio_names = [f'{x}/{y}^2' for x, y in zip(Dp_names, H_names)]
Dps = stats['means'][-2]
Hs = stats['means'][-1]
ratios = Dps / Hs ** 2
for name, val in zip(ratio_names, ratios):
    print(name, val)

