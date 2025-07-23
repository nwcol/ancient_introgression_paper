# Write ghost-Denisova introgression grid models.

import demes 


T_splits = [9e5, 1e6, 1.2e6, 1.4e6, 1.6e6, 1.8e6, 2e6, 2.5e6]
T_pulses = [1e5, 1.5e5, 2e5, 3e5, 4e5, 5e5]
p_pulses = [0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10]
has_h_n_pulse = [True, False]

num_models = len(T_splits) * len(T_pulses) * len(p_pulses) * len(has_h_n_pulse)
print(num_models)


