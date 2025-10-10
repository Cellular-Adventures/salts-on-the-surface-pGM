# Salts on the Surface - $p_{GM}$
This repository contains the files supporting the paper 'Salts on the Surface' -
(url to be added).
Here, you will find all the Python-based data processing and cleanup, and the 
code to create all figures. This code is supported by the $p_{GM}$ model as set 
up by [Duignan, 2021](https://doi.org/10.1016/j.jcis.2021.04.144), which is here 
implemented in Julia.


# Repository structure

| Folder                    | Content                                       |
| ------                    | -------                                       |
| [**01_data**](#01_data)   |                                               |
| /01_raw                   | Raw measurement data                          |
| /02_processed             | Processed measurement data                    |
| [**02_data-processing**](#02_data-processing)|                            |
| /01_consolidate-cleanup   | Python files used to consolidate separate measurements into a full dataset |
| /02_create-dataframe      | Python files used to summarize data in a `pandas` `DataFrame` |
| [**03_create-figures**](#03_create-figures) | Python files to create individual figures     |
| [**04_pGM-model**](#04_pgm-model)          | Git submodule: [GM_pressure](https://github.com/Cellular-Adventures/GM_pressure) |
| [**05_usage-examples**](#05_usage-examples) |                             |
| ...to be filled           |                                               |
| [**06_figures**](#06_figures) | Folders with the `png`, `svg` and `eps` versions of the figures |

# Usage
## Cloning the repo
<!-- [ ] write instructions on how to clone the repo, taking submodule into account -->
...

## 01_data
### Unpacking data
Find the raw data (zipped) (either in `01_data/01_raw/01_zipped` or somewhere online. TBD). Unpack the data into `01_data/01_raw/02_unpacked`.

## 02_data-processing
### Processing metadata
Run `consolidate_experiments.py` (`expand_input` in old setup) to run through the raw data and make an overview of all experiments, their concentrations and timestamps. This script outputs `consolidated-experiment-logs.csv`

### Lookup and calculate fluid properties
Run `extract_fluid_properties.py` to expand the output of `consolidate_experiments.py` with the fluid properties for each experiment. This script outputs `exp-with-fluid-properties.csv`.

### Pressure probe calibration
Run `calibrate_pressure_probe.py` to get the pressure probe calibration. This is needed to convert the raw signals to a hydrostatic pressure, and thus to a gas fraction. This script will generate figures of the regression and `probe-characterization.yaml` containing the regression parameters.

### Read in the **real** data
Run `read_experiment_data.py` to extract FP and PP data for the experiments. The output of this operation cannot be properly captured in a csv file anymore, so the output is an h5-formatted pandas dataframe `consolidated-measurement-data.h5` with the actual data in the `['data']` field.

### Group data per experiment
Run `process_experiments.py` to create a dataframe `per-experiment-data.h5` that combines all data for a single concentration in a single row. Expands data with mean holdups from PP and FP, standard deviation of those, bubble size distributions, characteristic bubble sizes of that distribution, interfacial area, holdup profile equation parameters

## 03_create-figures
Use the files in `03_create-figures` to reproduce the individual figures.

## 04_pGM-model
The files in `04_pGM-model` are used in `extract_fluid_properties.py` to calculate the GM pressure for each concentration. This folder is a git submodule, with the original repo at https://github.com/Cellular-Adventures/GM_pressure. 

## 05_usage-examples
`05_usage-examples` contains some further examples on how to use the pGM model implementation, both in Python and in Julia.

## 06_figures
All final figures are saved in `png`, `svg` and `eps` format in the folder `06_figures`.

# salts-on-the-surface-pGM
Repository containing the files supporting the Salts on the Surface paper. 
Salt-based pGM model, data analysis and the code to create all figures. The data 
itself probably won't fit in here.
