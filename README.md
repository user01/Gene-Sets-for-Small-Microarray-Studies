# Neuroimmunology Capstone

## Dependencies

### R

The process depends on a modern R (v3.3) and several R packages. A convenience script in `dependencies.R` is included to load the packages from CRAN.

* `MASS :: v7.3.45`
* `readr :: v1.0.0`
* `stringr :: v1.1.0`
* `dplyr :: v0.5.0`
* `purrr :: v0.2.2`
* `lubridate :: v1.6.0.9000`
* `argparse :: v1.0.4`
* `Rtsne :: v0.11`
* `randomForest :: v4.6.12`
* `EMCluster :: v0.2.6`
* `sparsediscrim :: v0.2.3`'

### Python

While requiring a python (v3.5) environment, there are no explicit python dependencies.

### Node

Install proper make dependencies with `npm install` with node (>5.0).

## Part I

Part I of the project consists of classifying cellular samples into Cell Type and General Cell Types.

The pipeline operates by:

1. Pre-processing raw data
2. Dimensional Reduction of data (from 20k features to <10)
3. Classification of data (clustering or other methods)
4. Performance Analysis

The pipeline is controlled by `node make_classify.js`.

## Part II

Part II of the project generates gene sets from the data.

Run this code with `node make_pairwise.js`.
