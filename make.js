const {
  run_lda_bootstrap,
  read_pairwise_union
} = require('./utilities/bootstrap.lda.js');
const R = require('ramda');

const runBaseline = () =>
  run_lda_bootstrap('baseline',   // Title
                    1,        // Bootstraps
                    1.0,      // Fraction to sample
                    4,        // concurrency
                    'python3', // python_binary
                    'fmeasure', // scoring mechanic
                    25,        // smallest gene set
                    200,        // largest gene set size
                    50,        // number of sets to create
                    500,       // number of gene pairs
                    false);


const runTask = (options={}) => {
  const opts = R.merge({
    bootstraps: 100,
    fraction: 0.66,
    scoring: 'fmeasure',
    low: 25,
    high: 200,
    count: 50,
    pairs: 500
  }, options);
  const title_str = R.pipe(
    R.values,
    R.join('.')
  )(opts);
  return () => run_lda_bootstrap(`main_${title_str}`,   // Title
                    opts.bootstraps,        // Bootstraps
                    opts.fraction,      // Fraction to sample
                    4,        // concurrency
                    'python3', // python_binary
                    opts.scoring, // scoring mechanic
                    opts.low,        // smallest gene set
                    opts.high,        // largest gene set size
                    opts.count,        // number of sets to create
                    opts.pairs,       // number of gene pairs
                    false);
};



Promise.resolve()
        .then(runBaseline)
        .then(runTask())
        .then(runTask({low:100,high:200}))
        .then(runTask({low:200,high:300}))
        .then(runTask({low:300,high:400}))
        .then(runTask({low:500,high:600}))
        .then(runTask({bootstraps:200,fraction:0.5}))
        .then(runTask({bootstraps:200,fraction:0.5,low:200,high:300}))
        .then(read_pairwise_union)
