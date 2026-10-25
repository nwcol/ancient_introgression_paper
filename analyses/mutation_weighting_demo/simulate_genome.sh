#!/bin/bash
unzip data.zip
python simulate_genome.py $@
rm -r data