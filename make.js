// Dimension reduction methods to run
const dimensionReductionPasses = [
  ['pca', 2],
  ['pca', 8]
];

// Clustering methods
const clusterMethods = [
  'em'
];
// Passes to try with each cluster
const clusterPasses = [
  ['pca_2', 8],
  ['pca_2', 12],
  ['pca_8', 8]
];

// *****************************************************************************
// Do not edit below this line
// *****************************************************************************
// Import Libraries
const path = require('path'),
  R = require('ramda'),
  chalk = require('chalk'),
  pad = require('pad');
moment = require('moment'),
  fs = require('fs'),
  Promise = require('bluebird'),
  spawn = require('child_process').spawn;

var allMs = 0;
const start = moment();

// Helper functions
const res = (filename) => path.join('results', filename);
const data = (filename) => path.join('data', filename);
const info = i => console.log(pad(80, chalk.blue.bold(i), ' '));

// Utility Functions
// Checks if file exists
const fsAccess = (path) => {
  return new Promise(function(resolve, reject) {
    fs.access(path, fs.constants.R_OK, (err) => {
      return (err ? reject : resolve)(path);
    });
  });
};
// Runs command with arguments
const cmd = (bin, args) => {
  return new Promise(function(resolve, reject) {
    const proc = spawn(bin + '', args.map(s => s + ''));
    var stdout = '';
    var stderror = '';
    proc.stdout.on('data', data => stdout += data);
    proc.stderr.on('data', data => stderror += data);
    proc.on('close', code => {
      // if (code == 0) {
      //   resolve(`${bin} ${args.join(' ')}`)
      // } else {
      //   console.log(stdout);
      //   console.log(stderror);
      //   reject(stdout);
      // }
      return (code == 0 ? resolve : reject)(`${bin} ${args.join(' ')}`);
    });
  });
};



// Performs action, if file does not exist
const make = (target, bin, args) => {
  const start = moment();
  return fsAccess(target)
    .catch(x => {
      console.log(` ${pad(19, chalk.blue('WORKING'), ' ')}:${pad(60, chalk.yellow(target), ' ')} : ${args.join(' ')}`);
      return cmd(bin, args);
    })
    .then(x => {
      const ms = moment().diff(start);
      allMs += ms;
      if (ms < 10) return;
      const seconds = Math.round(ms / 1000);
      console.log(` ${pad(19, chalk.green('COMPLETED'), ' ')}:${pad(10,`${seconds} sec`,' ')}${pad(50, chalk.yellow(target), ' ')} : ${args.join(' ')}`);
    })
};

// Main task functions
// Pre-process main data
const load_data = () => {
  return make(res('gene_data_vs_cell_type.tsv'), 'Rscript', ['load_data.R']);
}

// Perform all dimension reduction
const dimreduction = () => {
  return Promise.map(dimensionReductionPasses, (pass) => {
    return make(res(`dimreduced_${pass[0]}_${pass[1]}.tsv`),
      'Rscript', [`dimreduction_${pass[0]}.R`, '--dimensions', pass[1]]);
  }, {
    concurrency: 6
  });
};

// Perform all clustering
const cluster = () => {
  const combos = R.pipe(
    R.xprod(clusterMethods),
    R.map(R.flatten)
  )(clusterPasses);
  return Promise.map(combos, (combo) => {
    return make(res(`cluster_${combo[0]}_${combo[1]}_${combo[2]}c.tsv`),
      'Rscript', [`cluster_${combo[0]}.R`, '--name', combo[1], '--clusters',
        combo[2], '--plot', '--results'
      ]);
  }, {
    concurrency: 6
  });
}

// Read final results
const readResults = () => {
  return make(res('results_all.csv'), 'Rscript', ['read_results.R']);
}

// Full Pipeline
load_data()
  .then(x => info('Data Preprocessing Finished'))
  .then(dimreduction)
  .then(x => info('Dimensional Reduction Finished'))
  .then(cluster)
  .then(x => info('Clustering Finished'))
  .then(readResults)
  .then(x => {
    const wallTime = Math.round(moment().diff(start) / 1000)
    const procTime = Math.round(allMs / 1000);
    info(`All Tasks Finished - Proc: ${procTime}sec, Wall: ${wallTime}sec`);
  });

info('Starting tasks');