// Import Libraries
const R = require('ramda'),
  chalk = require('chalk'),
  path = require('path'),
  fs = require('fs'),
  pad = require('pad'),
  moment = require('moment'),
  Promise = require('bluebird');


const {
  make,
  z,
  data,
  readTsv,
  readTypes,
  info
} = require('./tools.js');

const run_lda_bootstrap = (title,
  bootstraps,
  fraction = 0.66,
  concurrency = 1,
  python_binary = 'python',
  lda_score_measure = 'fmeasure',
  test_specific_cells = false) => {

  const results_directory = path.join('results', title);

  if (!fs.existsSync(results_directory)) {
    fs.mkdirSync(results_directory);
  }

  //override local res to match set
  const res = (filename) => path.join(results_directory, filename);

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

  const localMake = (dependencyPath, binary, args, feedback = '') => {
    return make(dependencyPath, binary, args, noteRun, noteMs, feedback);
  }

  // Pre-process main data
  const load_data = () => {
    const output_path = res('gene_data_vs_cell_type.tsv');
    const output_validation_path = res('gene_data_vs_cell_type.validation.tsv');
    return Promise.all([
        localMake(output_path, 'Rscript', ['load_data.R', '--output', output_path]),
        localMake(output_validation_path, 'Rscript', ['load_data.R',
          '--output', output_validation_path,
          '--inputnormal', data('Amit_Norm_Counts.tsv'),
          '--inputmeta', data('Amit_Sample_Metadata.tsv')
        ]),
        readTypes(data('Gautier_Immgen_Sample_Metadata.tsv'))
      ], {
        concurrency
      })
      .then(R.nth(2));
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
      test_specific_cells ? types_to_tasks('Cell_Type') : []
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
    const path_output = res(`${title}.score${filename}`)

    const args_lda = [
      'score_lda.py',
      '--input',
      res('gene_data_vs_cell_type.tsv'),
      '--seed',
      task.bootstrap,
      '--bootstrap',
      `${fraction}`,
      '--type',
      task.type,
      '--score',
      lda_score_measure,
      '--name',
      `"${task.name}"`,
      '--output',
      path_output
    ];
    return localMake(
        path_output,
        python_binary,
        args_lda,
        `Task: ${task_done}/${task_count} ${Math.floor(100 * task_done/task_count)}%`)
      .then(x => {
        task_done++
      })
      .then(x => task);
  };

  const bootstrapAll = (tasks) => Promise.map(tasks, bootstrap, {
    concurrency
  });


  const buildSets = (task) => {
    const path_output = res('');
    const path_feedback = res(`${title}.set.feedback.${task.type}.${task.name.replace(/\s/g,'_')}.tsv`);

    const args = [
      'pairwise_buildsets.py',
      '--low',
      '15',
      '--high',
      '500',
      '--count',
      '20',
      '--type',
      task.type,
      '--name',
      `"${task.name}"`,
      '--title',
      title,
      '--raw',
      res('gene_data_vs_cell_type.tsv'),
      '--input',
      path_output,
      '--output',
      path_output, //TODO: modify code to accept name parameter
      '--feedback',
      path_feedback
    ];
    return localMake(path_feedback, python_binary, args)
      .then(x => R.merge(task, {
        path_feedback
      }));
  };
  const buildSetsAll = (tasks) => {
    const tasks_without_bootstrap = R.uniqBy(o => `${o.type}-${o.name}`, tasks);
    return Promise.map(tasks_without_bootstrap, buildSets, {
      concurrency
    });
  };

  const scoreSets = (task) => {
    console.log(task);

    const data_path = res('gene_data_vs_cell_type.tsv');
    const data_validation_path = res('gene_data_vs_cell_type.validation.tsv');
    const path_source = res('');
    const path_output = res(`${title}.set.results.${task.type}.${task.name.replace(/\s/g,'_')}.tsv`);

    const args = [
      'score_gsva.R',
      '--input',
      task.path_feedback,
      '--basedata',
      data_path,
      '--validationdata',
      data_validation_path,
      '--source',
      path_source,
      '--output',
      path_output
    ];
    return localMake(path_output, 'Rscript', args)
      .then(x => path_output);
    // return readTsv(feedbackPath)
    //   .then(feedback => {
    //     // console.log(feedback);
    //   })
  }

  const evaluateSets = (feedbackPaths) => {
    return Promise.map(feedbackPaths, scoreSets, {
      concurrency
    });
  };


  const readSets = () => {
    const path_sets = res(`${title}.sets.gmt`);
    const args = [
      'read_pairwise.py',
      '--raw',
      res('gene_data_vs_cell_type.tsv'),
      '--title',
      title,
      '--input',
      res(''),
      '--outputfull',
      res(`${title}.sets.full.tsv`),
      '--outputleader',
      res(`${title}.sets.leaders.tsv`),
      '--outputsets',
      path_sets
    ];
    return make(false, python_binary, args)
      .then(x => path_sets);
  };

  return load_data()
    .then(info(`Data Loaded`))
    .then(cell_types => {
      return celltypes_to_tasks(cell_types)
        .then(tasks => {
          task_count = tasks.length;
          return tasks;
        })
        .then(bootstrapAll)
        .then(info(`Bootstrap ${title} Complete`))
        .then(buildSetsAll)
        .then(info(`Sets ${title} Built`))
        .then(evaluateSets)
        .then(info(`Sets ${title} Evaluated`))
        .then(readSets)
        .then(info(`${title} Complete`))
    })
    .then(x => {
      const wallTime = moment.duration(moment().diff(start))
      const procTime = moment.duration(allMs);
      info(`${title} All Tasks Finished`);
      console.log(`${pad(71, chalk.blue("Processing Time:"))}${pad(20,procTime.asSeconds() + ' seconds')}${pad(20,procTime.humanize())}`);
      console.log(`${pad(71, chalk.blue("Wall Time:"))}${pad(20,wallTime.asSeconds() + ' seconds')}${pad(20,wallTime.humanize())}`);
    })
    .then(info(`${title} All Finished`));

}; // end of LDA bootstrap code


module.exports = {
  run_lda_bootstrap
};
