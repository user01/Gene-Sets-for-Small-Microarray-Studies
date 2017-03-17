// Import Libraries
const R = require('ramda'),
  chalk = require('chalk'),
  pad = require('pad'),
  moment = require('moment'),
  Promise = require('bluebird');

const {
  readJson,
  taskToName,
  taskArgsToNew,
  taskify,
  make,
  res,
  info
} = require('./tools.js');

const dimensionReductionData = readJson('data_dimensionreduction.json');
const clusterData = readJson('data_cluster.json');
// Change this parameter to adjust the maximum number of cores
const concurrency = 7;


const dimensionReductionTasks = taskify('dimreduction')(dimensionReductionData);
const clusterTasks = R.pipe(
  taskify('cluster'),
  R.xprod(R.__, dimensionReductionTasks),
  R.map(task => {

    const newArgs = taskArgsToNew(task[1], task[0][4]);
    const newName = taskToName(task[1], task[0][4]);
    const output_path = task[0][3].replace('||', newName);

    const fixedArgs = R.pipe(
      R.concat(R.__, ['--output', res(output_path)])
    )(newArgs);
    return R.pipe(
      R.head,
      R.remove(2, 3), // Remove the name, old target
      R.concat(R.__, [output_path, fixedArgs])
    )(task);
  })
)(clusterData);


// Performs action, if file does not exist
var anyMakeRun = false;
const noteRun = () => {
  anyMakeRun = true;
};

var allMs = 0;
const start = moment();
const noteMs = (ms) => {
  allMs += ms;
}

// *****************************************************************************
// Main task functions
// *****************************************************************************

// Pre-process main data
const load_data = () => {
  return make(res('gene_data_vs_cell_type.tsv'), 'Rscript', ['load_data.R'], noteRun, noteMs);
}

const args_to_outfile = (script_name, args) => {
  switch (script_name) {
    case 'dimreduction_pca.R':
      return `dimreduction_pca_dimensions_${args[1]}.tsv`;
    case 'dimreduction_tsne.R':
      return `dimreduction_tsne_perplexity_${args[1]}_${args[3]}.tsv`;
    default:
      console.error(`Unknown script name: ${script_name}`);
  }
  return 'unknown_data.tsv';
}

// Perform all dimension reduction
const dimreduction = () => {
  const correctedTasks = R.map(task => {
    const newName = task[3].replace('||_', '');
    return R.pipe(
      R.remove(3, 1),
      R.insert(2, newName),
      R.remove(3, 1)
    )(task);
  }, dimensionReductionTasks);
  return Promise.map(correctedTasks, (task) => {
    const outputPath = res(task[2]);
    const inputPath = res('gene_data_vs_cell_type.tsv');
    const args = R.pipe(
      R.prepend(task[1]),
      R.concat(R.__, ['--input', inputPath, '--output', outputPath])
    )(task[3]);
    return make(outputPath,
      task[0], args, noteRun, noteMs);
  }, {
      concurrency
    });
};

// Perform all clustering
const cluster = () => {
  return Promise.map(clusterTasks, (task) => {
    return make(res(task[2]),
      task[0], R.prepend(task[1], task[3]), noteRun, noteMs);
  }, {
      concurrency
    });
}

// Read final results
const readResults = () => {
  return make(anyMakeRun ? false : res('results_all.csv'),
    'Rscript', ['read_results.R'], noteRun, noteMs);
}

// *****************************************************************************
// Full Pipeline
// *****************************************************************************
load_data()
  .then(info('Data Preprocessing Finished'))
  .then(dimreduction)
  .then(info('Dimensional Reduction Finished'))
  .then(cluster)
  .then(info('Clustering Finished'))
  .then(readResults)
  .then(x => {
    const wallTime = moment.duration(moment().diff(start))
    const procTime = moment.duration(allMs);
    info(`All Tasks Finished`);
    console.log(`${pad(71, chalk.blue("Processing Time:"))}${pad(20, procTime.asSeconds() + ' seconds')}${pad(20, procTime.humanize())}`);
    console.log(`${pad(71, chalk.blue("Wall Time:"))}${pad(20, wallTime.asSeconds() + ' seconds')}${pad(20, wallTime.humanize())}`);
  });

info(`Starting ${clusterTasks.length + dimensionReductionTasks.length} tasks`);
