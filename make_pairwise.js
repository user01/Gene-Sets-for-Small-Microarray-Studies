const {
  run_lda_bootstrap,
  read_pairwise_union
} = require('./utilities/bootstrap.lda.js');

const runBase100_200 = () =>
  run_lda_bootstrap('100.200.base',   // Title
                    1,        // Bootstraps
                    1.0,      // Fraction to sample
                    4,        // concurrency
                    'python3', // python_binary
                    'fmeasure', // scoring mechanic
                    100,        // smallest gene set
                    200,        // largest gene set size
                    300,        // number of sets to create
                    false);

const runBase200_300 = () =>
  run_lda_bootstrap('200.300.base',   // Title
                    1,        // Bootstraps
                    1.0,      // Fraction to sample
                    4,        // concurrency
                    'python3', // python_binary
                    'fmeasure', // scoring mechanic
                    200,        // smallest gene set
                    300,        // largest gene set size
                    300,        // number of sets to create
                    false);


const runBase300_400 = () =>
  run_lda_bootstrap('300.400.base',   // Title
                    1,        // Bootstraps
                    1.0,      // Fraction to sample
                    4,        // concurrency
                    'python3', // python_binary
                    'fmeasure', // scoring mechanic
                    300,        // smallest gene set
                    400,        // largest gene set size
                    300,        // number of sets to create
                    false);

const runMain100_200 = () =>
  run_lda_bootstrap('100.200.main',   // Title
                    100,        // Bootstraps
                    0.66,      // Fraction to sample
                    4,        // concurrency
                    'python3', // python_binary
                    'fmeasure', // scoring mechanic
                    100,        // smallest gene set
                    200,        // largest gene set size
                    300,        // number of sets to create
                    false);

const runMain200_300 = () =>
  run_lda_bootstrap('200.300.main',   // Title
                    100,        // Bootstraps
                    0.66,      // Fraction to sample
                    4,        // concurrency
                    'python3', // python_binary
                    'fmeasure', // scoring mechanic
                    200,        // smallest gene set
                    300,        // largest gene set size
                    300,        // number of sets to create
                    false);


const runMain300_400 = () =>
  run_lda_bootstrap('300.400.main',   // Title
                    100,        // Bootstraps
                    0.66,      // Fraction to sample
                    4,        // concurrency
                    'python3', // python_binary
                    'fmeasure', // scoring mechanic
                    300,        // smallest gene set
                    400,        // largest gene set size
                    300,        // number of sets to create
                    false);


runBase100_200()
  .then(runBase200_300)
  .then(runBase300_400)
  .then(runMain100_200)
  .then(runMain200_300)
  .then(runMain300_400);
