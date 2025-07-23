datasets = ["Bherer", "ZhouJHS"]
models = [
    "alt_los",
    "alt_stu",
    "alt_ust",
    "cha_los",
    "cha_stu",
    "cha_ust",
    "vin_los",
    "vin_stu",
    "vin_ust",
]
with open("variables.txt", "w") as fout:
    for model in models:
        for dataset in datasets:
            fgraph = f"{model}.yaml"
            fparams = f"{model}_params.yaml"
            fdata = f"subset_{dataset}.pkl"
            foutstem = f"{model}_{dataset}"
            s = f"{fgraph},{fparams},{fdata},{foutstem}\n"
            fout.write(s)