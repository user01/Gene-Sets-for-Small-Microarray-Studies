# Neuroimmunology Capstone

## Dependencies

### R

The process depends on a modern R (>3.3.*) and several R packages. A convenience script in `dependencies.R` is included to load the packages from CRAN.

### Python

While requiring a python (>3.5) environment, there are no explicit python dependencies.

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
