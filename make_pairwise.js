const concurrency = 6;
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
  readFeedback,
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
    ], {
      concurrency
    })
    .then(R.nth(1));
};

const celltypes_to_tasks = (celltypes) => {
  // console.log(celltypes);
  const types_to_tasks = (type) => {
    return R.pipe(
      R.prop(type),
      R.xprod(R.range(0, bootstraps)),
      R.map(([bootstrap, name]) => {
        return {
          type,
          name,
          bootstrap
        };
      })
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

  const args_lda = [
    'score_lda.py',
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
  return make(
      path_output,
      'python',
      args_lda,
      noteRun,
      noteMs,
      `Task: ${task_done}/${task_count} ${Math.floor(100 * task_done/task_count)}%`)
    .then(x => {
      task_done++
    })
    .then(x => task);
}

const bootstrapAll = (tasks) => Promise.map(tasks, bootstrap, {
  concurrency
});


const buildSets = (task) => {
  const path_output = res('');
  const path_feedback = res(`set.feedback.${task.type}.${task.name}.tsv`);

  const args = [
    'pairwise_buildsets.py',
    '--low',
    '15',
    '--high',
    '200',
    '--count',
    '10',
    '--type',
    task.type,
    '--name',
    `"${task.name}"`,
    '--raw',
    res('gene_data_vs_cell_type.tsv'),
    '--input',
    path_output,
    '--output',
    path_output,
    '--feedback',
    path_feedback
  ];
  return make(path_feedback, 'python', args)
    .then(x => path_feedback);
}
const buildSetsAll = (tasks) => {
  const tasks_without_bootstrap = R.uniqBy(o => `${o.type}-${o.name}`, tasks);
  return Promise.map(R.take(5, tasks_without_bootstrap), buildSets, {
    concurrency
  });
}

const evaluateSet = (feedback) => {
  return readFeedback(feedback);
};

const evaluateSets = (feedbacks) => {
  return Promise.map(feedbacks, evaluateSet, {
    concurrency
  });
};


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
      .then(buildSetsAll)
      .then(info('Sets Built'))
      .then(evaluateSets)
  })
  .then(x => {
    const wallTime = moment.duration(moment().diff(start))
    const procTime = moment.duration(allMs);
    info(`All Tasks Finished`);
    console.log(`${pad(71, chalk.blue("Processing Time:"))}${pad(20,procTime.asSeconds() + ' seconds')}${pad(20,procTime.humanize())}`);
    console.log(`${pad(71, chalk.blue("Wall Time:"))}${pad(20,wallTime.asSeconds() + ' seconds')}${pad(20,wallTime.humanize())}`);
  })
  .then(info('All Finished'));
