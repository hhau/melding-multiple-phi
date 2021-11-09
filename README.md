# Code to accompany _Combining chains of Bayesian models with Markov melding_

## Re-running this code

- Re-running the code requires [MIMIC-III](https://physionet.org/content/mimiciii/1.4/) with all the standard [derived concept tables](https://github.com/MIT-LCP/mimic-code/tree/main/mimic-iii/concepts) on a Postgres database, running and available at the connection give in [these](scripts/mimic-example/get-baseline-data.R) [three](scripts/mimic-example/get-blood-gasses-and-define-pf-cohort.R) [files](scripts/mimic-example/get-raw-fluids.R).
  - The first time you rerun the analysis you may need to uncomment and run the `CREATE OR REPLACE` commands at the top of [`blood-gasses.sql`](scripts/mimic-example/queries/blood-gasses.sql).
- The R package dependencies are managed by `renv`. Some of the packages are not available on CRAN and are my own forks (the BA template is a fork of `rticles` that I haven't had the time to submit a PR for), but `renv` should install the correct versions.
- The analysis is handled by `GNU Make`, and you will need the latest version (>= 4.3) to run `make`.
  - Running Make with `-j 2` or higher will parallelise much of the computation. The MCMC samplers run 5 parallel chains, so ensure you have at least 5 times the number of CPU cores as the requested number of `Make` jobs.
- I think some parts of the analysis assume certain folders exist but don't necessarily check that they do indeed exist.

## Diagnostics

The MCMC diagnostic plots and tables referred to in the text can be found [here for the little owls example](rmd-reports/2021-06-30_owls-diagnostics.html) and [here for the MIMIC respiratory failure example](rmd-reports/2021-06-15_diagnosis-issues.html).
You will need to download the `.html` files and open them in a browser to view the diagnostics.

## Other things

- The logistic regression model discussed in the discussion is [here](scripts/mimic-example/queries/blood-gasses.sql), and the correspondence from the MIMIC team is [here](https://github.com/MIT-LCP/mimic-code/issues/1033).
