const concurrency = 7;

// Import Libraries
const R = require('ramda'),
  chalk = require('chalk'),
  pad = require('pad'),
  moment = require('moment'),
  Promise = require('bluebird');

const {
  readJson,
  taskToName,
  taskify,
  make,
  z,
  res,
  pairwise,
  info
} = require('./tools.js');


// *****************************************************************************
// Main task functions
// *****************************************************************************

// Pre-process main data
const load_data = () => {
  return make(res('gene_data_vs_cell_type.tsv'), 'Rscript', ['load_data.R']);
}

const classify_components = (idx, path_components) => {
  const path_pred_rf = pairwise(`predict_${z(idx)}_randomforest.tsv`);
  const args_rf = [
    'classify_randomforest.R',
    '--seed',
    idx,
    '--input',
    path_components,
    '--output',
    path_pred_rf
  ];
  return make(path_pred_rf, 'Rscript', args_rf)
    .then(x => [path_pred_rf]);
}

const score_predictions = (idx, pred_paths, pairs_path) => {
  const path_general = pairwise(`score_general_${z(idx)}.tsv`);
  const path_specific = pairwise(`score_specific_${z(idx)}.tsv`);
  const args = [
    'pairwise_score.R',
    '--inputpairs',
    pairs_path,
    '--inputpredictions',
    pred_paths.join(','),
    '--outputgeneral',
    path_general,
    '--outputspecific',
    path_specific
  ];
  return make(path_general, 'Rscript', args);
}

const bootstrap = (idx) => {
  const path_outputpairs = pairwise(`pairs_${z(idx)}.tsv`);
  const path_outputcomponents = pairwise(`components_${z(idx)}.tsv`);
  const args = [
    'pairwise_bootstrap.R',
    '--input',
    res('gene_data_vs_cell_type.tsv'),
    '--outputpairs',
    path_outputpairs,
    '--outputcomponents',
    path_outputcomponents,
    '--seed',
    idx
  ];

  return make(path_outputpairs, 'Rscript', args)
    .then(x => info(`Boostrap Started for ${idx}`))
    .then(x => classify_components(idx, path_outputcomponents))
    .then(path_preds => score_predictions(idx,path_preds,path_outputpairs))
    .then(x => info(`Boostrap Finished for ${idx}`))
    // .then(path_preds => info(`Prediction files ${path_preds.join('--')}`));
}

const bootstrapAll = () => Promise.map(R.range(0,2), bootstrap, {concurrency});

load_data()
  .then(x => info('Data Preprocessing Finished'))
  .then(bootstrapAll)
