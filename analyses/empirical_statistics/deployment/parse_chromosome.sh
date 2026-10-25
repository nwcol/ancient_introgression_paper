#!/bin/bash
echo $@
set -e
python parse_chromosome.py $@
