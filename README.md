# Inferring hominin history with recurrent gene flow from single unphased genomes and a two-locus statistic

This repository holds data and analyses for the paper preprinted at
https://www.biorxiv.org/content/10.64898/2026.04.11.717825v1.
Fitted demographic models can be found in the `models/main_models` directory.

## Requirements

To run the scripts and python notebooks in this repository, you need to install
`dpluspy`, a Python package available at [https://github.com/nwcol/dpluspy](https://github.com/nwcol/dpluspy),
and its dependencies (see [https://github.com/nwcol/dpluspy/blob/main/requirements.txt](https://github.com/nwcol/dpluspy/blob/main/requirements.txt)).

```sh
pip install git+https://github.com/nwcol/dpluspy.git
```

Or, you can clone the `dpluspy` repository and install it locally with

```sh
git clone https://github.com/nwcol/dpluspy.git
cd dpluspy
pip install -r requirements.txt
pip install .
```
