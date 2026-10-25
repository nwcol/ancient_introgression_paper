models = [
    "model_000",
    "model_001",
    "model_010",
    "model_100",
    "model_110",
    "model_101",
    "model_011",
    "model_111"
]

with open("fit_variables.txt", "w") as fout:
    for rep in range(500):
        for model in models:
            graph = model + ".yaml"
            params = model + "_params.yaml"
            data = f"rep_{rep}.pkl"
            output = f"rep_{rep}_{model}"
            params = f"{graph},{params},{data},{output}\n"
            fout.write(params)
