reps = list(range(100))
datas = ["Bherer", "ZhouJHS"]
models = ["den_null", "den_contraction"]
with open("bootstrap_variables.txt", "w") as fout:
    for model in models:
        for data in datas:
            for rep in reps:
                fgraph = f"{model}_{data}_MLE.yaml"
                fparams = f"{model}_params.yaml"
                fdata = f"subset_{data}.pkl"
                ftbl = f"{model}_{data}_rep{rep}.csv"
                s = f"{fgraph},{fparams},{fdata},{rep},{ftbl}\n"
                fout.write(s)