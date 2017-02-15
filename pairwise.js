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

const bootstrap = (idx) => {
  const path_outputpairs = pairwise(`pairs_${idx}.tsv`);
  const path_outputcomponents = pairwise(`components_${idx}.tsv`);
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
    .then(x => info(`Data Preprocessing Finished ${x}`));
}

const bootstrapAll = () => Promise.map(R.range(0,2), bootstrap);

load_data()
  .then(x => info('Data Preprocessing Finished'))
  .then(bootstrapAll)
