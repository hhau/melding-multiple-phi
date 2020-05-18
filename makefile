RSCRIPT = Rscript
PLOT_SETTINGS = scripts/common/plot-settings.R
TEX_FILES = $(wildcard tex-input/*.tex) \
	$(wildcard tex-input/*/*.tex) \
	$(wildcard tex-input/*/*/*.tex)
MCMC_UTIL = scripts/common/mcmc-util.R
BIBLIOGRAPHY = bibliography/multi-phi-bib.bib

# useful compound make components
PLOTS = plots
RDS = rds
SCRIPTS = scripts

# ALL figures
ALL_PLOTS =

# if you wildcard the all-target, then nothing will happen if the target doesn't
# exist (no target). hard code the target.
# CHANGE THIS:
WRITEUP = multiple-phi.pdf

all : $(WRITEUP)

clean : 
	trash *.aux *.out

## Three Gaussian example basenames
EX_THREE_GAUSSIANS = ex-three-gaussians

### data
MODEL_1_DATA = $(RDS)/$(EX_THREE_GAUSSIANS)/model-1-data.rds
MODEL_1_DATA : scripts/$(EX_THREE_GAUSSIANS)/data-generation-1.R
	$(RSCRIPT) $<

MODEL_3_DATA = $(RDS)/$(EX_THREE_GAUSSIANS)/model-3-data.rds
MODEL_3_DATA : scripts/$(EX_THREE_GAUSSIANS)/data-generation-3.R
	$(RSCRIPT) $<

MODEL_2_DATA = $(RDS)/$(EX_THREE_GAUSSIANS)/model-2-data.rds
MODEL_2_DATA : scripts/$(EX_THREE_GAUSSIANS)/data-generation-2.R $(MODEL_1_DATA) $(MODEL_3_DATA)
	$(RSCRIPT) $<


## Pooling visualisation tests
POOLING_TESTS = pooling-tests
POOLING_SCRIPTS = scripts/$(POOLING_TESTS)

POOLED_PLOT_2D = plots/pooling-tests/pooled-densities-2d.pdf
$(POOLED_PLOT_2D) : $(POOLING_SCRIPTS)/visualisation.R $(POOLING_SCRIPTS)/density-functions.R $(PLOT_SETTINGS)
	$(RSCRIPT) $<

ALL_PLOTS += $(POOLED_PLOT_2D)

################################################################################
## Owls example
OWLS_BASENAME = owls-example
OWLS_DATA = $(wildcard rds/owls-example/*-data.rds)
OWLS_SCRIPTS = $(SCRIPTS)/$(OWLS_BASENAME)
OWLS_RDS = $(RDS)/$(OWLS_BASENAME)
OWLS_PLOTS = $(PLOTS)/$(OWLS_BASENAME)
OWLS_POSTERIOR_SAMPLES = $(wildcard rds/owls-example/*-samples.rds)

$(OWLS_DATA) : $(OWLS_SCRIPTS)/load-and-write-data.R
	$(RSCRIPT) $<

ORIG_IPM_SAMPLES = $(OWLS_RDS)/original-ipm-samples.rds
$(ORIG_IPM_SAMPLES) : $(OWLS_SCRIPTS)/fit-original-ipm.R $(OWLS_SCRIPTS)/models/original-ipm.bug $(OWLS_DATA) $(MCMC_UTIL)
	$(RSCRIPT) $<

FECUNDITY_SUBPOSTERIOR = $(OWLS_RDS)/fecundity-subposterior-samples.rds 
$(FECUNDITY_SUBPOSTERIOR) : $(OWLS_SCRIPTS)/fit-fecundity.R $(OWLS_SCRIPTS)/models/fecundity-model.stan $(FECUNDITY_DATA)
	$(RSCRIPT) $< 

FECUNDITY_DIAGNOSTIC = $(OWLS_PLOTS)/stage-one-diagnostics-fecundity.png
$(FECUNDITY_DIAGNOSTIC) : $(OWLS_SCRIPTS)/diagnostics-stage-one.R
	$(RSCRIPT) $<

CAPTURE_RECAPTURE_DIAGNOSTIC = $(OWLS_PLOTS)/stage-one-diagnostics-capture-recapture.png
$(CAPTURE_RECAPTURE_DIAGNOSTIC) : $(FECUNDITY_DIAGNOSTIC)

ALL_PLOTS += $(FECUNDITY_DIAGNOSTIC) $(CAPTURE_RECAPTURE_DIAGNOSTIC)

CAPTURE_RECAPTURE_SUBPOSTERIOR = $(OWLS_RDS)/capture-recapture-subposterior-samples.rds
$(CAPTURE_RECAPTURE_SUBPOSTERIOR) : $(OWLS_SCRIPTS)/fit-capture-recapture.R $(OWLS_SCRIPTS)/models/capture-recapture.bug $(OWLS_DATA) $(MCMC_UTIL)
	$(RSCRIPT) $<

COUNT_DATA_SUBPOSTERIOR = $(OWLS_RDS)/count-data-subposterior-samples.rds
$(COUNT_DATA_SUBPOSTERIOR) : $(OWLS_SCRIPTS)/fit-count-data.R $(OWLS_SCRIPTS)/models/count-data.bug $(OWLS_DATA) $(MCMC_UTIL)
	$(RSCRIPT) $<

SUBPOSTERIOR_PLOT = $(OWLS_PLOTS)/subposteriors.pdf
$(SUBPOSTERIOR_PLOT) : $(OWLS_SCRIPTS)/plot-subposteriors.R $(PLOT_SETTINGS) $(MCMC_UTIL) $(ORIG_IPM_SAMPLES) $(FECUNDITY_SUBPOSTERIOR) $(CAPTURE_RECAPTURE_SUBPOSTERIOR) $(COUNT_DATA_SUBPOSTERIOR)
	$(RSCRIPT) $<

ALL_PLOTS += $(SUBPOSTERIOR_PLOT)

MELDED_POSTERIOR = $(OWLS_RDS)/melded-posterior-samples.rds
$(MELDED_POSTERIOR) : $(OWLS_SCRIPTS)/mcmc-main-stage-two.R $(OWLS_SCRIPTS)/mcmc-nimble-functions.R $(MCMC_UTIL) $(FECUNDITY_SUBPOSTERIOR) $(CAPTURE_RECAPTURE_SUBPOSTERIOR) $(OWLS_DATA)
	$(RSCRIPT) $<

MELDED_DIAGNOSTIC_PLOT = $(OWLS_PLOTS)/stage-two-diagnostics.png
$(MELDED_DIAGNOSTIC_PLOT) : $(OWLS_SCRIPTS)/diagnostics-stage-two.R $(MELDED_POSTERIOR) $(PLOT_SETTINGS)
	$(RSCRIPT) $<

MELDED_QQ_PLOT = $(OWLS_PLOTS)/orig-meld-qq-compare.pdf
$(MELDED_QQ_PLOT) : $(OWLS_SCRIPTS)/plot-qq-comparison.R $(MELDED_POSTERIOR) $(ORIG_IPM_SAMPLES) $(PLOT_SETTINGS)
	$(RSCRIPT) $<

ALL_PLOTS += $(MELDED_DIAGNOSTIC_PLOT) $(MELDED_QQ_PLOT)

################################################################################

# knitr is becoming more picky about encoding, specify UTF-8 input
$(WRITEUP) : $(wildcard *.rmd) $(TEX_FILES) $(ALL_PLOTS) $(OWLS_DATA) $(BIBLIOGRAPHY)
	$(RSCRIPT) -e "rmarkdown::render(input = Sys.glob('*.rmd'), encoding = 'UTF-8')"