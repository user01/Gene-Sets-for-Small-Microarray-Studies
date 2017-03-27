const {
  run_lda_bootstrap,
  read_pairwise_union
} = require('./utilities/bootstrap.lda.js');

const runBase = () =>
  run_lda_bootstrap('base',   // Title
                    1,        // Bootstraps
                    1.0,      // Fraction to sample
                    4,        // concurrency
                    'python3',
                    'fmeasure',
                    true);

const runMain = () =>
  run_lda_bootstrap('full',   // Title
                    100,      // Bootstraps
                    0.66,     // Fraction to sample
                    4,        // concurrency
                    'python3',
                    'fmeasure',
                    true);

const runAccuracy = () =>
  run_lda_bootstrap('accuracy',   // Title
                    100,      // Bootstraps
                    0.66,     // Fraction to sample
                    4,        // concurrency
                    'python3',
                    'accuracy',
                    true);

const runSmall = () =>
  run_lda_bootstrap('small',   // Title
                    200,      // Bootstraps
                    0.33,     // Fraction to sample
                    4,        // concurrency
                    'python3',
                    'fmeasure',
                    true);

const runTiny = () =>
  run_lda_bootstrap('tiny',   // Title
                    300,      // Bootstraps
                    0.2,     // Fraction to sample
                    4,        // concurrency
                    'python3',
                    'fmeasure',
                    true);



// runBase()
runBase()
  .then(runMain)
  .then(runAccuracy)
  .then(runSmall)
  .then(runTiny)
  .then(read_pairwise_union);
