

graph = "full_model_guess.yaml"
params = "full_model_params.yaml"

with open("variables.txt", "w") as fout:
    for u in ("1.1e-8", "1.3e-8", "1.5e-8"):
        for rmap in ["bherer", "zhou-JHS"]:
            for i in range(30):
                vars = ",".join([
                    u,
                    "-1",
                    graph,
                    params,
                    f"accessible_1e-4M_ccds_flank_{rmap}_adj.pkl",
                    f"main_fit_{rmap}_{u}_try_{i}.yaml"
                ]) + "\n"
                fout.write(vars)
            for i in range(500):
                vars = ",".join([
                    u,
                    str(i),
                    graph,
                    params,
                    f"accessible_1e-4M_ccds_flank_{rmap}_adj_boot_reps.pkl",
                    f"boot_fit_{rmap}_{u}_rep_{i}.yaml"
                ]) + "\n"
                fout.write(vars)





