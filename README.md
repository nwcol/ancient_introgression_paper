## Inferring hominin history with recurrent gene flow from single unphased genomes and a two-locus statistic

This repository holds data and analyses for the paper preprinted at
https://www.biorxiv.org/content/10.64898/2026.04.11.717825v1.

### Requirements

To run the scripts and python notebooks in this repository, you will need to install
`h2py`, a Python package available at [https://github.com/nwcol/h2py](https://github.com/nwcol/h2py)
and its dependencies.
You can find documentation and examples in the linked repository.

```sh
pip install git+https://github.com/nwcol/h2py.git
```

Or, you can clone the `h2py` repository and install it locally with

```sh
git clone https://github.com/nwcol/h2py.git
cd h2py
pip install -r requirements.txt
pip install .
```

### Contents

Small data files and Snakefiles/directions for reproducing sequence data,
recombination and mutation maps, and filters can be found in
[data/](data/). To recompute H<sub>2</sub> statistics, run the Snakefile in
[workflow/](workflow/).

Focal best-fit models can be found as loose Demes-format YAML files in
[models/main/](models/main/) and in [models/main/submodels/](models/main/submodels/).
