#!/bin/bash
set -e


model=$1
rep=$2


# INPUTS
true_graph="models/model${model}_graph.yaml"
graph0="models/model0_graph.yaml"
params0="models/model0_params.yaml"
graph1="models/model1_graph.yaml"
params1="models/model1_params.yaml"


# OUTPUTS
data_file="model${model}_${rep}_stats.pkl.gz"
fit0="model${model}_${rep}_model0_fit.yaml"
fit1="model${model}_${rep}_model1_fit.yaml"
lrt_tbl="model${model}_${rep}_LRT.csv"
boot_lrt_tbl="model${model}_${rep}_bootstrapLRT.csv"


python scripts/simulate_data.py $true_graph $data_file

python scripts/fit_model.py $data_file $graph0 $params0 $fit0

python scripts/fit_model.py $data_file $graph1 $params1 $fit1

python scripts/adjusted_LRT.py $data_file $fit0 $params0 $fit1 $params1 $lrt_tbl

python scripts/bootstrap_LRT.py $data_file $fit0 $params0 $fit1 $params1 $boot_lrt_tbl

