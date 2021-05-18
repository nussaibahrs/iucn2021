Morphological traits of reef corals predict extinction risk but not
conservation status
================

[![](https://img.shields.io/badge/doi-10.17605/OSF.IO/GB5TY-orange.svg)](https://doi.org/10.17605/OSF.IO/GB5TY)
[![](https://img.shields.io/github/languages/code-size/nussaibahrs/iucn2021.svg)](https://github.com/nussaibahrs/iucn2021)

*Nussaïbah B. Raja, Andreas Lauchstedt, John M. Pandolfi, Sun W. Kim,
Ann F. Budd & Wolfgang Kiessling*

  - [Description](#description)
  - [Requirements](#requirements)
  - [Setup](#setup)
  - [Scripts](#scripts)
  - [Troubleshooting](#troubleshooting)

## Description

This repository contains all the **data** (in `/data`) and **R scripts**
(in `/scripts`) necessary to develop and evaluate a model to predict
extinction risk of corals in modern ocean using fossil data. The outputs
of the scripts are provided in the `/output` and figures in the `/figs`
folder. Please cite the study as:

Raja, NB, Lauchstedt, A, Pandolfi, JM, Kim, SW, Budd, AF, Kiessling, W. Morphological traits of reef corals predict extinction risk but not conservation status. Global Ecol Biogeogr. 2021; 00: 1– 12. https://doi.org/10.1111/geb.13321 
## Requirements

This code was developed in `R 4.0.0`. It is therefore recommended to use
the same or any more up-to-date version of R for reproducing the
analyses in this study.

## Setup

You will need to either use the Rstudio project environment or set your
working directory to the root of this folder.

To install all required depdendencies (packages), run:

``` r
source(file.path("inst","dependencies"))
```

## Scripts

The `scripts/` folder contains all the code generated for the above
mentioned study. The folder contains **3** folders for each of the
models trained in the study and one additional script `functions.R`
containing custom functions used in this study.

  - **01-model\_fmm:** This folder contains all required scripts to
    train and evaluate the fossil morphology model (trained using fossil
    data).

  - **02-model\_mmm:** This folder contains all required scripts to
    train and evaluate the modern morphology model (trained using modern
    morphology data) and the predict the extinction risk of modern
    corals (using both the fmm and mmm).

  - **03-model\_mdmm:** This folder contains all required scripts to
    train and evaluate the modern morphology and distributional model
    (trained using modern data).

## Troubleshooting

The issue tracker is the preferred channel for bug reports. You may also
contact me [by email](mailto:nussaibah.raja.schoob@fau.de).
