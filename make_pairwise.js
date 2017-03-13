const concurrency = 2;
const bootstraps = 100;

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
var task_count = 0;
var task_done = 0;

var anyMakeRun = false;
const noteRun = () => {
  anyMakeRun = true;
};

var allMs = 0;
const start = moment();
const noteMs = (ms) => {
  allMs += ms;
};

// Pre-process main data
const load_data = () => {
  return Promise.all([
    make(res('gene_data_vs_cell_type.tsv'), 'Rscript', ['load_data.R'], noteRun, noteMs),
    readTypes(data('Gautier_Immgen_Sample_Metadata.tsv'))
  ], {concurrency})
  .then(R.nth(1));
};

const celltypes_to_tasks = (celltypes) => {
  // console.log(celltypes);
  const types_to_tasks = (type) => {
    return R.pipe(
      R.prop(type),
      R.xprod(R.range(0, bootstraps)),
      R.map(([bootstrap, name]) => { return {type, name, bootstrap}; })
    )(celltypes);
  };
  const tasks = R.concat(
    types_to_tasks('General_Cell_Type'),
    // types_to_tasks('Cell_Type')
    []
  );
  return Promise.resolve(tasks);
};



const bootstrap = (task) => {
  const filename = `.${
    z(task.bootstrap)
  }.${
    task.type
  }.${
    task.name.replace(/\s+/g, '_')
  }.tsv`;
  const path_output = res(`score${filename}`)
  const path_feedback = res(`feedback${filename}`)

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
    path_output,
    '--feedback',
    path_feedback
  ];
  return make(
      path_feedback, // Using feedback allows passes to be skipped if unable to run
      'Rscript',
      args_lda,
      noteRun,
      noteMs,
      `Task: ${task_done}/${task_count} ${Math.floor(100 * task_done/task_count)}%`)
    .then(x => { task_done++ })
    .then(x => task);
}

const bootstrapAll = (tasks) => Promise.map(tasks, bootstrap, {concurrency});


const buildSets = (type) => {
  const path_output = pairwise(`sets_${type}.gmt`);
  const pattern = `score_${type}_\\\\d+.tsv`;
  const args = [
    'pairwise_buildsets.R',
    '--scorespath',
    'results',
    '--scorespattern',
    pattern,
    '--output',
    path_output
  ];
  return make(path_output, 'Rscript', args);
}
const buildSetsAll = () => Promise.map(['General_Cell_Type','Cell_Type'], buildSets, {concurrency});



load_data()
  .then(info('Data Loaded'))
  .then(cell_types => {
    return celltypes_to_tasks(cell_types)
      .then(tasks => {
        task_count = tasks.length;
        return tasks;
      })
      .then(bootstrapAll)
      .then(info('Bootstrap Complete'))
      .then(x => scoreAll(cell_types))
  })
  .then(x => {
    const wallTime = moment.duration(moment().diff(start))
    const procTime = moment.duration(allMs);
    info(`All Tasks Finished`);
    console.log(`${pad(71, chalk.blue("Processing Time:"))}${pad(20,procTime.asSeconds() + ' seconds')}${pad(20,procTime.humanize())}`);
    console.log(`${pad(71, chalk.blue("Wall Time:"))}${pad(20,wallTime.asSeconds() + ' seconds')}${pad(20,wallTime.humanize())}`);
  })
  .then(info('All Finished'));
