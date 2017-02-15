const path = require('path'),
  R = require('ramda'),
  chalk = require('chalk'),
  pad = require('pad'),
  moment = require('moment'),
  fs = require('fs'),
  Promise = require('bluebird'),
  spawn = require('child_process').spawn;


const readJson = R.pipe(
  fs.readFileSync,
  (buf) => buf.toString(),
  JSON.parse
);

const binToExtension = (binary) => {
  switch (binary) {
    case 'Rscript':
      return 'R';
    case 'python':
      return 'py';
    default:
      return 'unknown';
  }
};


const handleParameters = R.pipe(
  R.mapObjIndexed((val, key) => {
    const keys = R.repeat(`--${key}`, R.length(val));
    const values = R.map(i => '' + i, val);
    return R.zip(keys, values);
  }),
  R.values,
  (values) => {
    return R.reduce(R.xprod, R.head(values), R.tail(values));
  },
  R.map(R.flatten)
);

const filenameObjTo = R.curry((leader, extension, parameters, method_name) => {
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
    R.map(R.unnest)
  )(paramChainsNested);
});


const taskify = (leader) => {
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

const taskToName = (task, existingArgs) => {
  const name = R.pipe(
    R.nth(4),
    R.map(R.replace('--', '')),
    R.join('_')
  )(task);
  return R.pipe(
    R.append('--name'),
    R.append(`${task[2]}_${name}`)
  )(existingArgs);
};



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


const make = (target, bin, args, noteRun = () => {}, noteMs = () => {}) => {
  const start = moment();
  const logTarget = () => {
    console.log(` ${pad(140, chalk.yellow(target), ' ')}`);
    console.log(` ${pad(140, chalk.white(args.join(' ')), ' ')}`);
  }
  return fsAccess(target)
    .catch(x => {
      noteRun();
      console.log(` ${pad(80, chalk.blue('WORKING'), ' ')}`);
      logTarget();
      return cmd(bin, args);
    })
    .then(x => {
      const ms = moment().diff(start);
      noteMs(ms);
      if (ms < 50) return;
      const seconds = Math.round(ms / 1000);
      console.log(` ${pad(80, `${seconds} seconds  ` + chalk.green('COMPLETED'), ' ')}`);
      logTarget();
    })
    .catch(err => {
      console.log(` ${pad(19, chalk.red('FAILED'), ' ')}`);
      console.error(err);
      logTarget();
    });
};

module.exports = {
  readJson,
  binToExtension,
  handleParameters,
  filenameObjTo,
  taskify,
  taskToName,
  fsAccess,
  cmd,
  make,
  res: (filename) => path.join('results', filename),
  pairwise: (filename) => path.join('pairwise', filename),
  data: (filename) => path.join('data', filename),
  info: i => console.log(pad(140, chalk.blue.bold(i), ' '))
};
