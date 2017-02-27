const path = require('path'),
  R = require('ramda'),
  chalk = require('chalk'),
  pad = require('pad'),
  csvparse = require('csv-parse'),
  moment = require('moment'),
  fs = require('fs'),
  Promise = require('bluebird'),
  spawn = require('child_process').spawn;

const lineWidth = 120;

const readJson = R.pipe(
  fs.readFileSync,
  (buf) => buf.toString(),
  JSON.parse
);


const readTsv = (tsvFile) => {
  return new Promise(function(resolve, reject) {
    fs.readFile(tsvFile, (err, tsv) => {
      if (err) {
        reject(err);
        return;
      }
      csvparse(tsv.toString(), {
          delimiter: '\t'
        },
        function(err, output) {
          if (err) {
            reject(err);
          } else {
            resolve(output);
          }
        });
    });
  });
};

const readTypes = (typesFile) => {
  return readTsv(typesFile)
    .then(types => {
      const elms = R.tail(types);
      const general = R.pipe(
        R.map(R.nth(2)),
        R.uniq
      )(elms);
      const specific = R.pipe(
        R.map(R.nth(1)),
        R.uniq
      )(elms);
      return { general, specific };
    });
};


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
        reject(`${stdout}\n\n${stderror}`);
      }
    });
  });
};


const make = (target, bin, args, noteRun = () => {}, noteMs = () => {}, message = '') => {
  const start = moment();
  const elapsed = () => {
    const ms = moment().diff(start);
    noteMs(ms);
    if (ms < 50) return;
    return Math.round(ms / 1000);
  }
  const logTarget = (color, response='', info='', err = false) => {
    const header = R.pipe(
      R.filter(s => s.length > 0),
      R.join(' ')
    )([color(response), chalk.gray(message), info]);
    console.log(` ${color('┌-----')}${pad(` ${header} `, lineWidth, color('-'))}`);
    console.log(` ${color('|')} ${chalk.yellow(target)}`);
    console.log(` ${color('|')} ${chalk.white(args.join(' '))}`);
    if (err) {
      console.log('ERROR')
      console.log(err);
    }
    console.log(` ${color(`└${pad('', lineWidth - 5, '-')}`)}`);
  }
  return fsAccess(target)
    .catch(x => {
      noteRun();
      logTarget(chalk.blue, 'WORKING', '');
      return cmd(bin, args);
    })
    .then(x => {
      logTarget(chalk.green, 'COMPLETED', `${elapsed()} seconds`);
    })
    .catch(err => {
      console.log(err);
      logTarget(chalk.red, 'ERROR', `${elapsed()} seconds`, err);
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
  readTypes,
  cmd,
  make,
  z: (s) => pad(5, s+'', '0'),
  res: (filename) => path.join('results', filename),
  data: (filename) => path.join('data', filename),
  info: i => {
    return (innerValue) => {
      console.log(pad(lineWidth, chalk.blue.bold(i), ' '))
      return innerValue;
    }
  }
};
