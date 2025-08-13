
import dpluspy 


data_file = "../../../../data/statistics/main/subset_stats.pkl"

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

u = 1.3e-8

for model in models:
    print(f"Testing model {model} setup ", end="", flush=True)
    graph_file =  f"{model}.yaml"
    param_file = f"{model}_params.yaml"
    pop_ids, bins, means, varcovs = dpluspy.inference.load_stats(
        data_file, graph=graph_file)
    for i in range(5):
        dpluspy.inference.optimize(
            graph_file,
            param_file,
            means,
            varcovs,
            pop_ids=pop_ids,
            bins=bins,
            u=u,
            perturb=1,
            method="fmin",
            max_iter=1,
            verbose=False,
        )
        print(".", end="", flush=True)
    print("")