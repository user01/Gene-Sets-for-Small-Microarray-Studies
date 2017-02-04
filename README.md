# Neuroimmunology Capstone

## Part I

Part I of the project consists of classifying cellular samples into Cell Type and General Cell Types.

The pipeline operates by:

1. Pre-processing raw data
2. Dimensional Reduction of data (from 20k features to <10)
3. Classification of data (clustering or other methods)
4. Performance Analysis

The pipeline is controlled by `main.js`. Install proper make dependencies with `npm install` and run the process with `npm start` or `node main.js`.

Add components to the pipeline by amending the start of `main.js` and adding the appropriate script.

### Dimensional Reduction format

Scripts must match the form of `dimreduction_{method}.{extension}`.

Scripts must accept parameters of the form `--parameter_name value` and must always be key/value pairs (ie no flags).

Scripts must write to the `results/` directory with tsv data of the form `dimreduction_{method}_[{key}_{value}].tsv`, with every key/value. The data must have the truth labels for each observation (Cell_Type and General_Cell_Type) and >2 other columns of floating values.
