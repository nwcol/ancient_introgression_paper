
import demes


for ii in range(694, 1387):
    in_fname = f"../../../ghost_to_denisova_pulse_grid/results/grid_model_{ii}_fit.yaml"
    out_fname = f"models/ghost_den_start_{ii - 694}.yaml"
    try:
        g = demes.load(in_fname)
        g._add_symmetric_migration(demes=["Altai", "Denisova"], rate=1e-4, start_time=2e5)
        demes.dump(g, out_fname)
    except:
        print(f"Missing {in_fname}")

in_fname = f"../../../ghost_to_denisova_pulse_grid/results/model_0_fit3.yaml"
out_fname = f"models/ghost_den_null_model.yaml"
g = demes.load(in_fname)
g._add_symmetric_migration(demes=["Altai", "Denisova"], rate=1e-4, start_time=2e5)
demes.dump(g, out_fname)
