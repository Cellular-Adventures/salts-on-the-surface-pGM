#!/usr/bin/env bash
set -ex

# This is the master script for the capsule. When you click "Reproducible Run",
#the code in this file will execute.

# The below script is technically the first to run, but it is dependent on
# user-interpretation of comments in the experimental logs.
# That does not work in a scripting setup, so an output file
# (/data/processed/consolidated-experiment-logs.csv) is directly included.

# python -u "02_data-processing/01_consolidate-cleanup/consolidate_experiments.py" "$@"

python -u "02_data-processing/01_consolidate-cleanup/extract_fluid_properties.py" "$@"
python -u "02_data-processing/01_consolidate-cleanup/calibrate_pressure_probe.py" "$@"
python -u "02_data-processing/01_consolidate-cleanup/process_literature_data.py" "$@"
python -u "02_data-processing/02_create-dataframe/read_experiment_data.py" "$@"
python -u "02_data-processing/02_create-dataframe/combine_experiment_data.py" "$@"

python -u "03_create-figures/figure1.py" "$@"
python -u "03_create-figures/figure2.py" "$@"
python -u "03_create-figures/figure3.py" "$@"
python -u "03_create-figures/figure5.py" "$@"
python -u "03_create-figures/figureB1.py" "$@"
