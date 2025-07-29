
dataset = "Bherer"

with open("variables.txt", "w") as fout:
    # We will fit the null model in 10 replicates to be safe about convergence
    for i in range(10):
        fgraph = "model0.yaml"
        model_set = 0
        fdata = f"subset_{dataset}.pkl"
        foutstem = f"model_0_fit{i}"
        s = f"{fgraph},{model_set},{fdata},{foutstem}\n"
        fout.write(s)

    for model in range(1, 694):
        fgraph = f"grid_model_{model}.yaml"
        model_set = 1
        fdata = f"subset_{dataset}.pkl"
        foutstem = f"grid_model_{model}_fit"
        s = f"{fgraph},{model_set},{fdata},{foutstem}\n"
        fout.write(s)

    for model in range(694, 1387):
        fgraph = f"grid_model_{model}.yaml"
        model_set = 2
        fdata = f"subset_{dataset}.pkl"
        foutstem = f"grid_model_{model}_fit"
        s = f"{fgraph},{model_set},{fdata},{foutstem}\n"
        fout.write(s)

