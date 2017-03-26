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



// runBase()
runBase()
  .then(runMain)
  .then(read_pairwise_union);
