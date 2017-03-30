const {
  run_lda_bootstrap,
  read_pairwise_union
} = require('./utilities/bootstrap.lda.js');


const runBaseline = () =>
  run_lda_bootstrap('baseline',   // Title
                    1,        // Bootstraps
                    1.0,      // Fraction to sample
                    4,        // concurrency
                    'python3', // python_binary
                    'fmeasure', // scoring mechanic
                    100,        // smallest gene set
                    600,        // largest gene set size
                    50,        // number of sets to create
                    8000,       // number of gene pairs
                    false);


const runMain100_200_big = () =>
  run_lda_bootstrap('100_200_main_big',   // Title
                    100,        // Bootstraps
                    0.66,      // Fraction to sample
                    4,        // concurrency
                    'python3', // python_binary
                    'fmeasure', // scoring mechanic
                    100,        // smallest gene set
                    200,        // largest gene set size
                    300,        // number of sets to create
                    8000,       // number of gene pairs
                    false);

const runMain200_300_big = () =>
  run_lda_bootstrap('200_300_main_big',   // Title
                    100,        // Bootstraps
                    0.66,      // Fraction to sample
                    4,        // concurrency
                    'python3', // python_binary
                    'fmeasure', // scoring mechanic
                    200,        // smallest gene set
                    300,        // largest gene set size
                    300,        // number of sets to create
                    8000,       // number of gene pairs
                    false);


const runMain300_400_big = () =>
  run_lda_bootstrap('300_400_main_big',   // Title
                    100,        // Bootstraps
                    0.66,      // Fraction to sample
                    4,        // concurrency
                    'python3', // python_binary
                    'fmeasure', // scoring mechanic
                    300,        // smallest gene set
                    400,        // largest gene set size
                    300,        // number of sets to create
                    8000,       // number of gene pairs
                    false);


Promise.resolve()
        .then(runBaseline)
        .then(runMain200_300_big)
        .then(runMain300_400_big)
        .then(runMain100_200_big)
        .then(read_pairwise_union)
