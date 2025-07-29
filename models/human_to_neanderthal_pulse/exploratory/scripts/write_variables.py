datasets = ["Bherer", "ZhouJHS"]
models = [
    "alt_kho",
    "alt_mbu",
    "alt_yor",
    "cha_kho",
    "cha_mbu",
    "cha_yor",
    "vin_kho",
    "vin_mbu",
    "vin_yor",
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