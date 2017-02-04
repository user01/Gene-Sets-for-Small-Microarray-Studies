# Neuroimmunology Capstone

## Part I

Part I of the project consists of classifying cellular samples into Cell Type and General Cell Types.

The pipeline operates by:

1. Pre-processing raw data
2. Dimensional Reduction of data (from 20k features to <10)
3. Classification of data (clustering or other methods)
4. Performance Analysis

The pipeline is controlled by `make.js`. Install proper make dependencies with `npm install` and run the process with `npm start` or `node make.js`.

Add components to the pipeline by amending the start of `make.js` and adding the appropriate script.

The process depends on a modern R (>3.3.*) and python (>3.5) environment. The following dependencies are required:

#### R

 * EMCluster
 * Rtsne
 * readr
 * dplyr
 * purrr
 * stringr
 * argparse

#### Python

Currently, no extra packages are required

### Dimensional Reduction format

Scripts must match the form of `dimreduction_{method}.{extension}`.

Scripts must accept parameters of the form `--parameter_name value` and must always be key/value pairs (ie no flags).

Scripts must write to the `results/` directory with tsv data of the form `dimreduction_{method}_[{key}_{value}].tsv`, with every key/value. The data must have the truth labels for each observation (Cell_Type and General_Cell_Type) and >2 other columns of floating values.

### Classification format

Scripts must match the form of `cluster_{classification_method}_{source_name}.{extension}`.

Here `source_name` is of the form `{dimreduction_method}_[{key}_{value}]`. This perfectly matches the `{method}_[{key}_{value}]` set from dimension reduction.

Scripts must accept parameters of the form `--parameter_name value` and must always be key/value pairs (ie no flags).

Scripts must write to the `results/` directory with tsv data of the form `cluster_{classification_method}_{source_name}.tsv`, with every key/value. For each observation index number (idx) truth labels (Cell_Type and General_Cell_Type) and predicted labels (Cell_Type_Predicted General_Cell_Type_Predicted).
