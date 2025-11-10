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
| [**data**](#data)   |                                               |
| /01_raw                   | Raw measurement data                          |
| /02_processed             | Processed measurement data                    |
| [**code/01_data-processing**](#code/01_data-processing)|                            |
| /01_consolidate-cleanup   | Python files used to consolidate separate measurements into a full dataset |
| /02_create-dataframe      | Python files used to summarize data in a `pandas` `DataFrame` |
| [**code/02_create-figures**](#code/02_create-figures) | Python files to create individual figures     |
| [**code/03_pGM-model**](#code/03_pGM-model)          | Git submodule: [GM_pressure](https://github.com/Cellular-Adventures/GM_pressure) |
| [**results/figures**](#results/figures) | Folders with the `png`, `svg` and `eps` versions of the figures |

# Usage
## Cloning the repo
<!-- [ ] write instructions on how to clone the repo, taking submodule into account -->
Clone the repo in however way you want to. The repo contains the zipped raw data and a submodule for the pGM model.

## Basic requirements
You will need to have `python` and `julia` installed. We rely on `conda` for `python` package management, so also make sure you have that available. The setup was tested with Anaconda, but it might also work with miniconda or other setups.

`julia` is best installed through juliaup. This setup was tested with Windows, where the `julia`install through the Microsoft Store uses juliaup. If you're using another operating system, check the instructions given by the lovely `julia` people.

## Python environment
The `python` environment can be generated through conda:

```shell
conda env create -f environment.yml
conda activate salts-on-the-surface-pGM
```

## Julia environment for pGM
The `julia` environment is a bit more tricky. We rely on [JuliaCall](https://juliapy.github.io/PythonCall.jl/stable/juliacall/) to make the link between Python and Julia, and on [PyJuliaPkg](https://github.com/JuliaPy/PyJuliaPkg) for the package management. These modules also dictate the install instructions.

The packages required for the pGM model are summarized in `juliapkg.json`, but we need to make sure that PyJuliaPkg actually finds it. Because we are working in an anaconda environment (at least in my case), PyJuliaPkg will start looking for this file in there (for me at `C:\Users\rikvolger\anaconda3\envs\salts-on-the-surface-pGM`). Specifically, it will look in the folder `julia_env\pyjuliapkg`, which does not exist yet (if you've not ran any scripts using juliacall). Create this folder, and in there, paste the `juliapkg.json` file.

If you're using `bash` (e.g. Git Bash):
```bash
mkdir -p ~/anaconda3/envs/salts-on-the-surface-pGM/julia_env/pyjuliapkg
cp juliapkg.json ~/anaconda3/envs/salts-on-the-surface-pGM/julia_env/pyjuliapkg/juliapkg.json
```
If you're using cmd or powershell, find the equivalent commands in there and do those.

Now, if you run for example `process_literature_data.py`, it should output a bunch of bubbles in the vague shape of the letters PGM (the bubbles might be boxes with question marks in there if your terminal does not support unicode)

## Data
### Unpacking data
Find the raw data (zipped) in `data/01_raw/01_zipped`. With your unzipper of choice (mine is 7zip), unpack the data into `data/01_raw/02_unpacked`. The literature data is expected in `data/01_raw/02_unpacked/01_Literature-Data`, the other data is expected directly in `data/01_raw/02_unpacked`.

## Data processing
### Processing metadata
Run `consolidate_experiments.py` (`expand_input` in old setup) to run through the raw data and make an overview of all experiments, their concentrations and timestamps. This script outputs `consolidated-experiment-logs.csv`

### Lookup and calculate fluid properties
Run `extract_fluid_properties.py` to expand the output of `consolidated-experiment-logs.csv` with the fluid properties for each experiment. This script outputs `exp-with-fluid-properties.csv`.

### Pressure probe calibration
Run `calibrate_pressure_probe.py` to get the pressure probe calibration. This is needed to convert the raw signals to a hydrostatic pressure, and thus to a gas fraction. This script will generate figures of the regression and `probe-characterization.yaml` containing the regression parameters.

### Read in the **real** data
Run `read_experiment_data.py` to extract FP and PP data for the experiments. The output of this operation cannot be properly captured in a csv file anymore, so the output is an h5-formatted pandas dataframe `consolidated-measurement-data.h5` with the actual data in the `['data']` field.

### Group data per experiment
Run `process_experiments.py` to create a dataframe `per-experiment-data.h5` that combines all data for a single concentration in a single row. Expands data with mean holdups from PP and FP, standard deviation of those, bubble size distributions, characteristic bubble sizes of that distribution, interfacial area, holdup profile equation parameters

## Creating figures
Use the files in `code/02_create-figures` to reproduce the individual figures.

## pGM Model
The files in `code/03_pGM-model` are used in `extract_fluid_properties.py` to calculate the GM pressure for each concentration. This folder is a git submodule, with the original repo at https://github.com/Cellular-Adventures/GM_pressure. 

# Results
The results will appear in the `data/02_processed` folder (for the `.csv` and `.h5` files), and in the `results/figures/[eps|png|svg]` folders for the appropriate figure types. Each figure is saved in all three extensions. `png` for quick inspection, `svg` for post-processing in InkScape and `eps` for direct inclusion in LaTeX papers.