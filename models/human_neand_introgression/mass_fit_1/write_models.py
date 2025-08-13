# Exchange Yoruba1 with other samples

models = [
    "human_neand_model_0",
    "human_neand_model_1",
    "human_neand_model_2",
    "human_neand_model_3",
    "human_neand_model_4"
]
samples = {
    "mbu1": "Mbuti1", 
    "kho1": "KhomaniSan1"
}

for key in samples:
    for model in models:
        graph_fname = f"models/{model}_yor1.yaml"
        param_fname = f"models/{model}_yor1_params.yaml"

        with open(graph_fname, "r") as fin:
            graph_str = fin.read()
        graph_str = graph_str.replace("Yoruba1", samples[key])
        out_graph_fname = graph_fname.replace("yor1", key)
        with open(out_graph_fname, "w") as fout:
            fout.write(graph_str)

        with open(param_fname, "r") as fin:
            param_str = fin.read()
        param_str = param_str.replace("Yoruba1", samples[key])
        out_param_fname = param_fname.replace("yor1", key)
        with open(out_param_fname, "w") as fout:
            fout.write(param_str)

