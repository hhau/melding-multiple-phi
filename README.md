# Code to accompany `TITLE TBA`

# Re-running this code

- requires MIMIC, with all the standard dependent tables, and some extra functions, on a Postgres database, running and available at the connection give in [these]() [three]() [files]().
- The R package dependencies are managed by `renv`. Some of the packages are not available on CRAN / are my own forks, but `renv` should install the correct versions.
- The analysis is handled by `GNU Make`, and you will need the latest version (>= 4.3) to run `make`.
- I think some parts of the analysis assume certain folders exist but don't necessarily check that they do indeed exist.

# Diagnostics

The MCMC diagnostic plots and tables referred to in the text can be found [here for the little owls example]() and [here for the MIMIC respiratory failure example]().
You will need to download the `.html` files and open them in a browser to view the diagnostics.