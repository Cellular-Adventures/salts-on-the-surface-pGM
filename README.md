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
| **01_data**               |                                               |
| /01_raw                   | Raw measurement data                          |
| /02_processed             | Processed measurement data                    |
| /03_figures               | Processed data supporting individual figures  |
| **02_data-processing**    |                                               |
| /01_consolidate-cleanup   | Python files used to consolidate separate measurements into a full dataset |
| /02_create-dataframe      | Python files used to summarize data in a `pandas` `DataFrame` |
| /03_export-figure-data     | Python files to export the data supporting individual figures |
| **03_create-figures**     | Python files to create individual figures & notebook to create all |
| **04_pGM-model**          |                                               |
| ...to be filled           |                                               |
| **05_usage-examples**     |                                               |
| ...to be filled           |                                               |
| **06_figures**            | Folders with the png, svg and eps versions of the figures |

# Usage
Find the raw data (zipped) (either in `01_data/01_raw/01_zipped` or somewhere online. TBD). Unpack the data into `01_data/01_raw/02_unpacked`.

Run `consolidate_experiments.py` (`expand_input` in old setup) to run through the raw data and make an overview of all experiments, their concentrations and timestamps.

Run `extract_fluid_properties.py` to expand the output of `consolidate_experiments.py` with the fluid properties for each experiment.
TODO make this output a csv file instead of h5. More insightful for the casual user.

Run `pre_process.py` to extract FP and PP data for the experiments. The output of this operation cannot be properly captured in a csv file anymore, so the output is an h5-formatted pandas dataframe.

Run `process_experiments.py` to create a dataframe that combines all data for a single concentration in a single row. Expands data with mean holdups from PP and FP, standard deviation of those, bubble size distributions, characteristic bubble sizes of that distribution, interfacial area, holdup profile equation parameters

Run `extract_figure_data.py` to create dedicated csv or h5 files of the data underlying each final figure.

Use the files in `03_create-figures` to reproduce the individual figures.

The files in `04_pGM-model` are used in `extract_fluid_properties.py` to calculate the GM pressure for each concentration.

`05_usage-examples` contains some further examples on how to use the pGM model implementation, both in Python and in Julia.

All final figures are saved in `png`, `svg` and `eps` format in the folder `06_figures`.

# salts-on-the-surface-pGM
Repository containing the files supporting the Salts on the Surface paper. 
Salt-based pGM model, data analysis and the code to create all figures. The data 
itself probably won't fit in here.
