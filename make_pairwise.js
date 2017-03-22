const {
  run_lda_bootstrap
} = require('./utilities/bootstrap.lda.js');

run_lda_bootstrap('base',   // Title
                  1,        // Bootstraps
                  1.0,      // Fraction to sample
                  6)        // concurrency

run_lda_bootstrap('full',   // Title
                  100,      // Bootstraps
                  0.66,     // Fraction to sample
                  6)        // concurrency
