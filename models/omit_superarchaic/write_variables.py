

g = "full_no_sup_guess.yaml"
p = "full_no_sup_params.yaml"

with open("parameters.txt", "w") as fout:
    for rmap in ["bherer", "zhou-JHS"]:
            label = g.replace("_guess.yaml", "")
            for i in range(30):
                vars = ",".join([
                    "-1",
                    g,
                    p,
                    f"accessible_1e-4M_ccds_flank_{rmap}_adj.pkl",
                    f"main_fit_{label}_{rmap}_try_{i}.yaml"
                ]) + "\n"
                fout.write(vars)
            for i in range(500):
                vars = ",".join([
                    str(i),
                    g,
                    p,
                    f"accessible_1e-4M_ccds_flank_{rmap}_adj_boot_reps.pkl",
                    f"boot_fit_{label}_{rmap}_rep_{i}.yaml"
                ]) + "\n"
                fout.write(vars)





