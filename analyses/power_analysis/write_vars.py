 
graph0 = "model_0_0.yaml"
params0 = "params_0_0.yaml"
graph1 = "model_0_1.yaml"
params1 = "params_0_1.yaml"

for data in ("sfs", "h2"):
    with open(f"{data}_variables_model_0.txt", "w") as fout:
        # BLRT fits [to construct null distr]
        for rep in range(500):
            true_model = f"model_0_0.yaml"
            prefix = f"{data}_blrt_model_0_rep_{rep}"
            args = (
                true_model,
                graph0,
                params0,
                graph1,
                params1,
                prefix
            )
            line = ",".join(args) + "\n"
            fout.write(line)

        # Focal fits
        for truth in (0, 1, 2, 5, 10):
            for rep in range(100):
                true_model = f"model_0_{truth}.yaml"
                prefix = f"{data}_focal_model_0_{truth}_rep_{rep}"
                args = (
                    true_model,
                    graph0,
                    params0,
                    graph1,
                    params1,
                    prefix
                )
                line = ",".join(args) + "\n"
                fout.write(line)

# Model 1; Model 0 with "known" (fixed) pulse time

# Same null model (no pulse)
graph0 = "model_0_0.yaml"
params0 = "params_0_0.yaml"

graph1 = "model_1_1.yaml"
params1 = "params_1_1.yaml"

for data in ("sfs", "h2"):
    with open(f"{data}_variables_model_1.txt", "w") as fout:
        # BLRT fits [to construct null distr]
        for rep in range(500):
            true_model = f"model_0_0.yaml"
            prefix = f"{data}_blrt_model_1_rep_{rep}"
            args = (
                true_model,
                graph0,
                params0,
                graph1,
                params1,
                prefix
            )
            line = ",".join(args) + "\n"
            fout.write(line)

        # Focal fits
        for truth in (0, 1, 2, 5, 10):
            for rep in range(100):
                true_model = f"model_0_{truth}.yaml"
                prefix = f"{data}_focal_model_1_{truth}_rep_{rep}"
                args = (
                    true_model,
                    graph0,
                    params0,
                    graph1,
                    params1,
                    prefix
                )
                line = ",".join(args) + "\n"
                fout.write(line)


graph0 = "model_2_0.yaml"
params0 = "params_2_0.yaml"
graph1 = "model_2_1.yaml"
params1 = "params_2_1.yaml"

for data in ("sfs", "h2"):
    with open(f"{data}_variables_model_2.txt", "w") as fout:
        # BLRT fits [to construct null distr]
        for rep in range(500):
            true_model = f"model_2_0.yaml"
            prefix = f"{data}_blrt_model_2_rep_{rep}"
            args = (
                true_model,
                graph0,
                params0,
                graph1,
                params1,
                prefix
            )
            line = ",".join(args) + "\n"
            fout.write(line)

        # Focal fits
        for truth in (0, 1, 2, 5, 10):
            for rep in range(100):
                true_model = f"model_2_{truth}.yaml"
                prefix = f"{data}_focal_model_2_{truth}_rep_{rep}"
                args = (
                    true_model,
                    graph0,
                    params0,
                    graph1,
                    params1,
                    prefix
                )
                line = ",".join(args) + "\n"
                fout.write(line)


