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

# salts-on-the-surface-pGM
Repository containing the files supporting the Salts on the Surface paper. 
Salt-based pGM model, data analysis and the code to create all figures. The data 
itself probably won't fit in here.
