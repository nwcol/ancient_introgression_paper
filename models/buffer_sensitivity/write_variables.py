

graph = "full_model_guess.yaml"
params = "full_model_params.yaml"

with open("variables.txt", "w") as fout:
    for rmap in ["bherer", "zhou-JHS"]:
        for buffer in ["intergenic", "1e-4M_ccds_flank", "1e-3M_ccds_flank"]:
            for i in range(30):
                vars = ",".join([
                    "-1",
                    graph,
                    params,
                    f"accessible_{buffer}_{rmap}_adj.pkl",
                    f"main_fit_{buffer}_{rmap}_try_{i}.yaml"
                ]) + "\n"
                fout.write(vars)
            for i in range(500):
                vars = ",".join([
                    str(i),
                    graph,
                    params,
                    f"accessible_{buffer}_{rmap}_adj_boot_reps.pkl",
                    f"boot_fit_{buffer}_{rmap}_rep_{i}.yaml"
                ]) + "\n"
                fout.write(vars)





