const concurrency = 7;
const bootstraps = 5;

// Import Libraries
const R = require('ramda'),
  chalk = require('chalk'),
  pad = require('pad'),
  moment = require('moment'),
  Promise = require('bluebird');


const {
  make,
  z,
  res,
  data,
  readTypes,
  info
} = require('./tools.js');


// *****************************************************************************
// Main task functions
// *****************************************************************************

// Pre-process main data
const load_data = () => {
  return Promise.all([
    make(res('gene_data_vs_cell_type.tsv'), 'Rscript', ['load_data.R']),
    readTypes(data('Gautier_Immgen_Sample_Metadata.tsv'))
  ], {concurrency})
  .then(R.nth(1));
}

const celltypes_to_tasks = (celltypes) => {
  const types_to_tasks = (type) => {
    return R.pipe(
      R.prop(type),
      R.xprod(R.range(0, bootstraps)),
      R.map(([bootstrap, name]) => { return {type, name, bootstrap}; })
    )(celltypes);
  };
  return R.concat(
    types_to_tasks('general'),
    types_to_tasks('specific')
  );
};



const bootstrap = (task) => {
  const path_output = res(`score_${
    z(task.bootstrap)
  }_${
    task.type
  }_${
    task.name.replace(/\s+/g, '.')
  }.tsv`);
  const args_lda = [
    'score_lda.R',
    '--input',
    res('gene_data_vs_cell_type.tsv'),
    '--seed',
    task.bootstrap,
    '--type',
    task.type,
    '--name',
    `"${task.name}"`,
    '--output',
    path_output
  ];
  return make(path_output, 'Rscript', args_lda)
    .then(x => [path_output]);
}

const bootstrapAll = (tasks) => Promise.map(tasks, bootstrap, {concurrency});

load_data()
  .then(celltypes_to_tasks)
  .then(info('Data Loaded'))
  .then(bootstrapAll)
  .then(info('Bootstrap Complete'))
  // .then(score)
  .then(info('All Finished'));
