// Dimension reduction methods to run
const dimensionReductionData = {
  Rscript: {
    pca: {
      dimensions: [2, 3, 4, 5, 6, 7, 8]
    },
    tsne: {
      perplexity: [20, 25, 30, 35, 40],
      pca: [true, false]
    }
  }
};

// Clustering methods
const clusterData = {
  Rscript: {
    em: {
      clusters: [8, 10, 11, 12, 13, 14, 16]
    }
  }
};


// *****************************************************************************
// Do not edit below this line
// *****************************************************************************
// Import Libraries
const path = require('path'),
  R = require('ramda'),
  chalk = require('chalk'),
  pad = require('pad'),
  moment = require('moment'),
  fs = require('fs'),
  Promise = require('bluebird'),
  spawn = require('child_process').spawn;

var allMs = 0;
const start = moment();

// Helper functions
const res = (filename) => path.join('results', filename);
const data = (filename) => path.join('data', filename);
const info = i => console.log(pad(120, chalk.blue.bold(i), ' '));


var binToExtension = (binary) => {
  switch (binary) {
    case 'Rscript':
      return 'R';
    case 'python':
      return 'py';
    default:
      return 'unknown';
  }
};
var handleParameters = R.pipe(
  R.mapObjIndexed((val, key) => {
    const keys = R.repeat(`--${key}`, R.length(val));
    const values = R.map(R.toString, val);
    return R.zip(keys, values);
  }),
  R.values,
  (values) => {
    return R.reduce(R.xprod, R.head(values), R.tail(values));
  },
  R.map(R.flatten)
);

var filenameObjTo = R.curry((leader, extension, parameters, method_name) => {
  const filename = `${leader}_${method_name}.${extension}`;
  const paramChains = handleParameters(parameters);
  const filenames = R.repeat(filename, R.length(paramChains));
  const methodNames = R.repeat(method_name, R.length(paramChains));
  const filename_methods = R.zip(filenames, methodNames);
  const paramChainsNested = R.map(x => [x], paramChains);
  // `dimreduced_${task[0]}_${task[1]}.tsv`
  const resultFilenames = R.map(chain =>
    `${leader}_${method_name}_||_${chain.map(R.replace('--','')).join('_')}.tsv`,
    paramChains);
  const all = R.pipe(
    R.zip(filename_methods),
    R.map(R.unnest)
  )(resultFilenames)

  return R.pipe(
    R.zip(all),
    // R.zip(filename_methods),
    R.map(R.unnest)
  )(paramChainsNested);
});
var taskify = (leader) => {
  return R.pipe(
    R.mapObjIndexed((filedata, binary) => {
      const fileDatas = R.pipe(
        R.mapObjIndexed(filenameObjTo(leader, binToExtension(binary))),
        R.values,
        R.unnest
      )(filedata);
      return R.pipe(
        R.repeat(R.__, R.length(fileDatas)),
        R.zip(R.__, fileDatas),
        R.map(R.unnest)
      )(binary);
    }),
    R.values,
    R.unnest
  );
};

var taskToName = (task, existingArgs) => {
  const name = R.pipe(
    R.nth(4),
    R.map(R.replace('--', '')),
    R.join('_')
  )(task);
  return R.pipe(
    R.append('--name'),
    // R.append(task)
    R.append(`${task[2]}_${name}`)
  )(existingArgs);
};

var dimensionReductionTasks = taskify('dimreduction')(dimensionReductionData);
var clusterTasks = R.pipe(
  taskify('cluster'),
  R.xprod(R.__, dimensionReductionTasks),
  R.map(task => {
    const newArgs = taskToName(task[1], task[0][4]);
    const newTarget = task[0][3].replace('||', R.last(newArgs));
    return R.pipe(
      R.head,
      R.remove(2, 3), // Remove the name, old target
      R.concat(R.__, [newTarget, newArgs])
    )(task)
  })
)(clusterData);
clusterTasks

// Utility Functions
// Checks if file exists
const fsAccess = (path) => {
  if (path == false) {
    return Promise.reject(path);
  }
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
      if (code == 0) {
        resolve(`${bin} ${args.join(' ')}`)
      } else {
        console.log(`${chalk.red.bold(' !!! Error !!!')} in ${args.join(' ')}`);
        console.log(stdout);
        console.log(stderror);
        reject(stdout);
      }
    });
  });
};



// Performs action, if file does not exist
var anyMakeRun = false;
const make = (target, bin, args) => {
  const start = moment();
  return fsAccess(target)
    .catch(x => {
      anyMakeRun = true;
      console.log(` ${pad(19, chalk.blue('WORKING'), ' ')}:${pad(100, chalk.yellow(target), ' ')} : ${args.join(' ')}`);
      return cmd(bin, args);
    })
    .then(x => {
      const ms = moment().diff(start);
      allMs += ms;
      if (ms < 50) return;
      const seconds = Math.round(ms / 1000);
      console.log(` ${pad(19, chalk.green('COMPLETED'), ' ')}:${pad(10,`${seconds} sec`,' ')}${pad(90, chalk.yellow(target), ' ')} : ${args.join(' ')}`);
    })
    .catch(x => {
      console.log(` ${pad(19, chalk.red('FAILED'), ' ')}:${pad(10,` `,' ')}${pad(90, chalk.yellow(target), ' ')} : ${args.join(' ')}`);
    });
};

// Main task functions
// Pre-process main data
const load_data = () => {
  return make(res('gene_data_vs_cell_type.tsv'), 'Rscript', ['load_data.R']);
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
    return make(res(task[2]),
      task[0], R.prepend(task[1], task[3]));
  }, {
    concurrency: 6
  });
};

// Perform all clustering
const cluster = () => {
  return Promise.map(clusterTasks, (task) => {
    return make(res(task[2]),
      task[0], R.prepend(task[1], task[3]));
  }, {
    concurrency: 6
  });
}

// Read final results
const readResults = () => {
  return make(anyMakeRun ? false : res('results_all.csv'),
    'Rscript', ['read_results.R']);
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
    const wallTime = moment.duration(moment().diff(start))
    const procTime = moment.duration(allMs);
    info(`All Tasks Finished`);
    console.log(`${pad(71, chalk.blue("Processing Time:"))}${pad(20,procTime.asSeconds() + ' seconds')}${pad(20,procTime.humanize())}`);
    console.log(`${pad(71, chalk.blue("Wall Time:"))}${pad(20,wallTime.asSeconds() + ' seconds')}${pad(20,wallTime.humanize())}`);
  });

info(`Starting ${clusterTasks.length + dimensionReductionTasks.length} tasks`);
