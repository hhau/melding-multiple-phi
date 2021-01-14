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
BASENAME = multiple-phi
WRITEUP = $(BASENAME).pdf

all : $(WRITEUP)

clean : 
	trash \
		$(BASENAME).aux \
		$(BASENAME).out \
		$(BASENAME).pdf \
		$(BASENAME).tex \
		Rplots.pdf

################################################################################
## Pooling visualisation tests
POOLING_TESTS = pooling-tests
POOLING_SCRIPTS = scripts/$(POOLING_TESTS)
POOLING_OUTPUTS = rds/$(POOLING_TESTS)
POOLING_PLOTS = plots/pooling-tests

POOLED_PLOT_2D = $(POOLING_PLOTS)/version-two.pdf
$(POOLING_SCRIPTS)/sub-plot-maker.R : $(POOLING_SCRIPTS)/density-functions.R

$(POOLED_PLOT_2D) : $(POOLING_SCRIPTS)/plot-pooled-priors.R  $(PLOT_SETTINGS) $(POOLING_SCRIPTS)/sub-plot-maker.R
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

CAPTURE_RECAPTURE_SUBPOSTERIOR = $(OWLS_RDS)/capture-recapture-subposterior-samples.rds
$(CAPTURE_RECAPTURE_SUBPOSTERIOR) : $(OWLS_SCRIPTS)/fit-capture-recapture.R $(OWLS_SCRIPTS)/models/capture-recapture.bug $(OWLS_DATA) $(MCMC_UTIL)
	$(RSCRIPT) $<

FECUNDITY_DIAGNOSTIC_PLOT = $(OWLS_PLOTS)/stage-one-diagnostics-fecundity.png
$(FECUNDITY_DIAGNOSTIC_PLOT) : $(OWLS_SCRIPTS)/diagnostics-stage-one-fecundity-plot.R $(PLOT_SETTINGS) $(FECUNDITY_SUBPOSTERIOR)
	$(RSCRIPT) $<

CAPTURE_RECAPTURE_DIAGNOSTIC_PLOT = $(OWLS_PLOTS)/stage-one-diagnostics-capture-recapture.png
$(CAPTURE_RECAPTURE_DIAGNOSTIC_PLOT) : $(OWLS_SCRIPTS)/diagnostics-stage-one-capture-recapture-plot.R $(PLOT_SETTINGS) $(CAPTURE_RECAPTURE_SUBPOSTERIOR)
	$(RSCRIPT) $<

ALL_PLOTS += $(FECUNDITY_DIAGNOSTIC_PLOT) $(CAPTURE_RECAPTURE_DIAGNOSTIC_PLOT)

STAGE_ONE_DIAGNOSTIC_TABLE = tex-input/owls-example/appendix-info/0010-stage-one-diagnostics.tex
$(STAGE_ONE_DIAGNOSTIC_TABLE) : $(OWLS_SCRIPTS)/diagnostics-stage-one-table.R $(FECUNDITY_SUBPOSTERIOR) $(CAPTURE_RECAPTURE_SUBPOSTERIOR)
	$(RSCRIPT) $<

COUNT_DATA_SUBPOSTERIOR = $(OWLS_RDS)/count-data-subposterior-samples.rds
$(COUNT_DATA_SUBPOSTERIOR) : $(OWLS_SCRIPTS)/fit-count-data.R $(OWLS_SCRIPTS)/models/count-data.bug $(OWLS_DATA) $(MCMC_UTIL)
	$(RSCRIPT) $<

SUBPOSTERIOR_PLOT = $(OWLS_PLOTS)/subposteriors.pdf
$(SUBPOSTERIOR_PLOT) : $(OWLS_SCRIPTS)/plot-subposteriors.R $(PLOT_SETTINGS) $(MCMC_UTIL) $(ORIG_IPM_SAMPLES) $(FECUNDITY_SUBPOSTERIOR) $(CAPTURE_RECAPTURE_SUBPOSTERIOR) $(COUNT_DATA_SUBPOSTERIOR)
	$(RSCRIPT) $<

ALL_PLOTS += $(SUBPOSTERIOR_PLOT)

NORMAL_APPROX_MELDED_POSTERIOR = $(OWLS_RDS)/melded-posterior-normal-approx-samples.rds
$(NORMAL_APPROX_MELDED_POSTERIOR) : $(OWLS_SCRIPTS)/fit-normal-approx.R $(MCMC_UTIL) $(CAPTURE_RECAPTURE_SUBPOSTERIOR) $(FECUNDITY_SUBPOSTERIOR) $(OWLS_DATA) $(OWLS_SCRIPTS)/models/count-data-normal-approx.bug
	$(RSCRIPT) $<

MELDED_POSTERIOR = $(OWLS_RDS)/melded-posterior-samples.rds
$(MELDED_POSTERIOR) : $(OWLS_SCRIPTS)/mcmc-main-stage-two.R $(OWLS_SCRIPTS)/mcmc-nimble-functions.R $(MCMC_UTIL) $(FECUNDITY_SUBPOSTERIOR) $(CAPTURE_RECAPTURE_SUBPOSTERIOR) $(OWLS_DATA)
	$(RSCRIPT) $<

MELDED_DIAGNOSTIC_TABLE = tex-input/owls-example/appendix-info/0020-stage-two-diagnostics.tex
$(MELDED_DIAGNOSTIC_TABLE) : $(OWLS_SCRIPTS)/diagnostics-stage-two-table.R $(MELDED_POSTERIOR) 
	$(RSCRIPT) $<

MELDED_DIAGNOSTIC_PLOT = $(OWLS_PLOTS)/stage-two-diagnostics.png
$(MELDED_DIAGNOSTIC_PLOT) : $(OWLS_SCRIPTS)/diagnostics-stage-two-plot.R $(MELDED_POSTERIOR) $(PLOT_SETTINGS)
	$(RSCRIPT) $<

MELDED_QQ_PLOT = $(OWLS_PLOTS)/orig-meld-qq-compare.pdf
$(MELDED_QQ_PLOT) : $(OWLS_SCRIPTS)/plot-qq-comparison.R $(MELDED_POSTERIOR) $(ORIG_IPM_SAMPLES) $(PLOT_SETTINGS)
	$(RSCRIPT) $<

ALL_PLOTS += $(MELDED_DIAGNOSTIC_PLOT) $(MELDED_QQ_PLOT)

################################################################################
## surv-example
SURV_BASENAME = surv-example
SURV_SCRIPTS = $(SCRIPTS)/$(SURV_BASENAME)
SURV_RDS = $(RDS)/$(SURV_BASENAME)
SURV_PLOTS = $(PLOTS)/$(SURV_BASENAME)
SURV_MODELS = $(SURV_SCRIPTS)/models

SURV_GLOBAL_SETTINGS = $(SURV_SCRIPTS)/GLOBALS.R

### Submodel 1
SURV_SIMULATION_SETTINGS_ONE = $(SURV_RDS)/submodel-one-simulation-settings.rds
$(SURV_SIMULATION_SETTINGS_ONE) : $(SURV_SCRIPTS)/simulation-settings-submodel-one.R $(SURV_GLOBAL_SETTINGS)
	$(RSCRIPT) $<

SURV_SUBMODEL_ONE_SIMULATED_DATA = $(SURV_RDS)/submodel-one-simulated-data.rds
$(SURV_SUBMODEL_ONE_SIMULATED_DATA) : $(SURV_SCRIPTS)/simulate-data-submodel-one.R $(SURV_SIMULATION_SETTINGS_ONE) 
	$(RSCRIPT) $<

SURV_SUMODEL_ONE = $(SURV_MODELS)/submodel-one.stan

SURV_SUBMODEL_ONE_SIMULATED_OUTPUT = $(SURV_RDS)/submodel-one-output.rds
$(SURV_SUBMODEL_ONE_SIMULATED_OUTPUT) : $(SURV_SCRIPTS)/fit-submodel-one.R $(SURV_SUBMODEL_ONE_SIMULATED_DATA) $(SURV_SUMODEL_ONE)
	$(RSCRIPT) $<

SURV_ALL_SUBMODEL_ONE_INPUTS = $(SURV_SIMULATION_SETTINGS_ONE) $(SURV_SUBMODEL_ONE_SIMULATED_DATA) $(SURV_SUBMODEL_ONE_SIMULATED_OUTPUT)

SURV_SUBMODEL_ONE_POSTERIOR_PLOT = $(SURV_PLOTS)/submodel-one-posterior.pdf
$(SURV_SUBMODEL_ONE_POSTERIOR_PLOT) : $(SURV_SCRIPTS)/plot-submodel-one.R $(PLOT_SETTINGS) $(MCMC_UTIL) $(SURV_ALL_SUBMODEL_ONE_INPUTS)
	$(RSCRIPT) $<

ALL_PLOTS += $(SURV_SUBMODEL_ONE_POSTERIOR_PLOT)

### Submodel 2

SURV_SUBMODEL_TWO_SIMULATED_DATA = $(SURV_RDS)/submodel-two-simulated-data.rds
$(SURV_SUBMODEL_TWO_SIMULATED_DATA) : $(SURV_SCRIPTS)/simulate-data-submodel-two.R $(SURV_GLOBAL_SETTINGS) $(SURV_SIMULATION_SETTINGS_ONE)
	$(RSCRIPT) $<

### Submodel 3
SURV_SIMULATION_SETTINGS_THREE = $(SURV_RDS)/submodel-three-simulation-settings.rds
$(SURV_SIMULATION_SETTINGS_THREE) : $(SURV_SCRIPTS)/simulation-settings-submodel-three.R $(SURV_GLOBAL_SETTINGS)
	$(RSCRIPT) $<

SURV_SUBMODEL_THREE_SIMULATED_DATA = $(SURV_RDS)/submodel-three-simulated-data.rds
$(SURV_SUBMODEL_THREE_SIMULATED_DATA) : $(SURV_SCRIPTS)/simulate-data-submodel-three.R $(SURV_SIMULATION_SETTINGS_THREE) 
	$(RSCRIPT) $<

SURV_SUMODEL_THREE = $(SURV_MODELS)/submodel-three.stan
SURV_SUBMODEL_THREE_SIMULATED_OUTPUT = $(SURV_RDS)/submodel-three-output.rds
$(SURV_SUBMODEL_THREE_SIMULATED_OUTPUT) : $(SURV_SCRIPTS)/fit-submodel-three.R $(SURV_SUBMODEL_THREE_SIMULATED_DATA) $(SURV_SUMODEL_THREE)
	$(RSCRIPT) $<

SURV_ALL_SUBMODEL_THREE_INPUTS = $(SURV_SIMULATION_SETTINGS_THREE) $(SURV_SUBMODEL_THREE_SIMULATED_DATA) $(SURV_SUBMODEL_THREE_SIMULATED_OUTPUT)

SURV_SUBMODEL_THREE_POSTERIOR_PLOT = $(SURV_PLOTS)/submodel-three-posterior.pdf
$(SURV_SUBMODEL_THREE_POSTERIOR_PLOT) : $(SURV_SCRIPTS)/plot-submodel-three.R $(PLOT_SETTINGS) $(MCMC_UTIL) $(SURV_ALL_SUBMODEL_THREE_INPUTS)
	$(RSCRIPT) $<

### Compositional step / plots

################################################################################

# knitr is becoming more picky about encoding, specify UTF-8 input
$(WRITEUP) : $(wildcard *.rmd) $(TEX_FILES) $(ALL_PLOTS) $(OWLS_DATA) $(BIBLIOGRAPHY)
	$(RSCRIPT) -e "rmarkdown::render(input = Sys.glob('*.rmd'), encoding = 'UTF-8')"