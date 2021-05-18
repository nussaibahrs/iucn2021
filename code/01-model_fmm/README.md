# Morphological traits of reef corals predict extinction risk but not conservation status

[![](https://img.shields.io/badge/doi-10.17605/OSF.IO/GB5TY-orange.svg)](https://doi.org/10.17605/OSF.IO/GB5TY)
[![](https://img.shields.io/github/languages/code-size/nussaibahrs/iucn2021.svg)](https://github.com/nussaibahrs/iucn2021)

*Nussaïbah B. Raja, Andreas Lauchstedt, John M. Pandolfi, Sun W. Kim,
Ann F. Budd & Wolfgang Kiessling*

This folder contains all required scripts to train and evaluate the fossil morphology model (trained using fossil data).

* **00-load_data.R:** Loads required data

* **01-analysis_h20_automl_fossil.R:** Runs the AutoML algorithm using fossil data, saves the generated model and calculates the optimised decision threhold for each of the models. 

* **03-results_binary_fossil.R:** Computes the metrics for model evaluation.

## Troubleshooting

The issue tracker is the preferred channel for bug reports. You may also
contact me [by email](mailto:nussaibah.raja.schoob@fau.de).
