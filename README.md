# Spatial insurance of distinct ecological functions 

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Research compendium to reproduce analyses and figures of the following article:

> Spatial insurance of distinct ecological functions, by Mouquet N., Mahaut, L., Thuiller, T., Auber, A., 
> Casajus N., Enquist B.J., Gaüzère P., Loiseau N., Mouillot D., Munoz F., Villéger S. & Violle C.
> Submited to Ecology Letters

## Content

This repository is structured as follow:

- [`data/`](https://github.com/nmouquet/SAFE_insu/tree/main/data):
contains data required to reproduce figures and tables

- [`analysis/`](https://github.com/nmouquet/SAFE_insu/tree/main/analysis/):
contains subfolders organized by theme. Each folder contains R scripts to run 
specific analysis

- [`results/`](https://github.com/nmouquet/SAFE_insu/tree/main/results):
follows the structure of analyses. Contains intermediate results and the 
numeric results used to produce the figures

- [`outputs/`](https://github.com/nmouquet/SAFE_insu/tree/main/ouputs):
contains the figures produced for the article

- [`R/`](https://github.com/nmouquet/SAFE_insu/tree/main/R):
contains R functions developed for this project

- [`DESCRIPTION`](https://github.com/nmouquet/SAFE_insu/tree/main/DESCRIPTION):
contains project metadata (author, date, dependencies, etc.)



## Workflow
    
The script [`analysis/01_Simulations.R`](https://github.com/nmouquet/SAFE_insu/blob/main/analysis/01_Simulations.R) compute functional insurance with simulated metacommunities.

The script [`analysis/02_Divgrass_alps.R`](https://github.com/nmouquet/SAFE_insu/blob/main/analysis/02_Divgrass_alps.R) compute functional insurance for alpine herbaceous communities.

The script [`analysis/03_Birds_Global.R`](https://github.com/nmouquet/SAFE_insu/blob/main/analysis/03_Birds_Global.R) compute functional insurance for birds worldwide.


## Big files 

  Some files were not uploaded to GitHub due to their large size. They are available upon request:

    - `results/wdpa_I_union.RData`
    - `results/wdpa_I_II_union.RData`
    - `results/wdpa_I_IV_union.RData`
    - `data/BIG_FILES/`

## Figures 

Figures are stored in `outputs/`.

The following Figures and Tables can be reproduced with the script indicated in brackets (all in [`analysis/`](https://github.com/nmouquet/RLS_HUM_INT/blob/main/analysis/)):
    
- Figure 1 has been produced with other means.
- Figure 2 has been produced with other means.
- [Fig_3](https://github.com/nmouquet/SAFE_insu/tree/main/outputs), was produced by [`01_Simulations.R`](https://github.com/nmouquet/SAFE_insu/blob/main/analysis/01_Simulations.R)
- [Fig_4b](https://github.com/nmouquet/SAFE_insu/tree/main/outputs), [Fig_4c](https://github.com/nmouquet/SAFE_insu/tree/main/outputs), [Fig_4d](https://github.com/nmouquet/SAFE_insu/tree/main/outputs) were produced by [`02_Divgrass_alps.R`](https://github.com/nmouquet/SAFE_insu/blob/main/analysis/02_Divgrass_alps.R). Note that the Fig_4a is produced with tmap and pasted into the final figure. 
- [Fig_5a_1](https://github.com/nmouquet/SAFE_insu/tree/main/outputs), [Fig_5a_2](https://github.com/nmouquet/SAFE_insu/tree/main/outputs), [Fig_5a_3](https://github.com/nmouquet/SAFE_insu/tree/main/outputs), [Fig_5b](https://github.com/nmouquet/SAFE_insu/tree/main/outputs), [Fig_5c](https://github.com/nmouquet/SAFE_insu/tree/main/outputs) were produced by [`03_Birds_Global.R`](https://github.com/nmouquet/SAFE_insu/blob/main/analysis/03_Birds_Global.R).
- [Fig_6a](https://github.com/nmouquet/SAFE_insu/tree/main/outputs), [Fig_6b](https://github.com/nmouquet/SAFE_insu/tree/main/outputs), were produced by [`03_Birds_Global.R`](https://github.com/nmouquet/SAFE_insu/blob/main/analysis/03_Birds_Global.R)

- [S1_Fig_1](https://github.com/nmouquet/SAFE_insu/tree/main/outputs), was produced by [`01_Simulations.R`](https://github.com/nmouquet/SAFE_insu/blob/main/analysis/01_Simulations.R)
- [S2_Fig_1](https://github.com/nmouquet/SAFE_insu/tree/main/outputs), [S2_Fig_2](https://github.com/nmouquet/SAFE_insu/tree/main/outputs) were produced by [`02_Divgrass_alps`](https://github.com/nmouquet/SAFE_insu/blob/main/analysis/02_Divgrass_alps)
- [S3_Fig_1](https://github.com/nmouquet/SAFE_insu/tree/main/outputs), [S3_Fig_2](https://github.com/nmouquet/SAFE_insu/tree/main/outputs),[S3_Fig_3](https://github.com/nmouquet/SAFE_insu/tree/main/outputs),[S3_Fig_4](https://github.com/nmouquet/SAFE_insu/tree/main/outputs), [S3_Fig_5](https://github.com/nmouquet/SAFE_insu/tree/main/outputs) were produced by [`03_Birds_Global.R`](https://github.com/nmouquet/SAFE_insu/blob/main/analysis/03_Birds_Global.R)

