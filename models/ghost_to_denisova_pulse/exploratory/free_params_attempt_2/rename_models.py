
import demes


for ii in range(694, 1387):
    in_fname = f"../../../ghost_to_denisova_pulse_grid/results/grid_model_{ii}_fit.yaml"
    out_fname = f"models/ghost_den_start_{ii - 694}.yaml"
    g = demes.load(in_fname)
    demes.dump(g, out_fname)

in_fname = f"../../../ghost_to_denisova_pulse_grid/results/model_0_fit3.yaml"
out_fname = f"models/ghost_den_null_model.yaml"
g = demes.load(in_fname)
demes.dump(g, out_fname)
