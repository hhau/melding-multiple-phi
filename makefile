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
$(ORIG_IPM_SAMPLES) : $(OWLS_SCRIPTS)/fit-original-ipm.R $(OWLS_SCRIPTS)/models/original-ipm.bug $(OWLS_DATA)
	$(RSCRIPT) $<

FECUNDITY_SUBPOSTERIOR = $(OWLS_RDS)/fecundity-subposterior-samples.rds 
$(FECUNDITY_SUBPOSTERIOR) : $(OWLS_SCRIPTS)/fit-fecundity.R $(OWLS_SCRIPTS)/models/fecundity-model.stan $(FECUNDITY_DATA)
	$(RSCRIPT) $< 

CAPTURE_RECAPTURE_SUBPOSTERIOR = $(OWLS_RDS)/capture-recapture-subposterior-samples.rds
$(CAPTURE_RECAPTURE_SUBPOSTERIOR) : $(OWLS_SCRIPTS)/fit-capture-recapture.R $(OWLS_SCRIPTS)/models/capture-recapture.bug $(OWLS_DATA)
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
$(COUNT_DATA_SUBPOSTERIOR) : $(OWLS_SCRIPTS)/fit-count-data.R $(OWLS_SCRIPTS)/models/count-data.bug $(OWLS_DATA)
	$(RSCRIPT) $<

SUBPOSTERIOR_PLOT = $(OWLS_PLOTS)/subposteriors.pdf
$(SUBPOSTERIOR_PLOT) : $(OWLS_SCRIPTS)/plot-subposteriors.R $(PLOT_SETTINGS) $(MCMC_UTIL) $(ORIG_IPM_SAMPLES) $(FECUNDITY_SUBPOSTERIOR) $(CAPTURE_RECAPTURE_SUBPOSTERIOR) $(COUNT_DATA_SUBPOSTERIOR)
	$(RSCRIPT) $<

ALL_PLOTS += $(SUBPOSTERIOR_PLOT)

NORMAL_APPROX_MELDED_POSTERIOR = $(OWLS_RDS)/melded-posterior-normal-approx-samples.rds
$(NORMAL_APPROX_MELDED_POSTERIOR) : $(OWLS_SCRIPTS)/fit-normal-approx.R $(CAPTURE_RECAPTURE_SUBPOSTERIOR) $(FECUNDITY_SUBPOSTERIOR) $(OWLS_DATA) $(OWLS_SCRIPTS)/models/count-data-normal-approx.bug
	$(RSCRIPT) $<

MELDED_POSTERIOR = $(OWLS_RDS)/melded-posterior-samples.rds
$(MELDED_POSTERIOR) : $(OWLS_SCRIPTS)/mcmc-main-stage-two.R $(OWLS_SCRIPTS)/mcmc-nimble-functions.R $(FECUNDITY_SUBPOSTERIOR) $(CAPTURE_RECAPTURE_SUBPOSTERIOR) $(OWLS_DATA)
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
SURV_TEX = tex-input/$(SURV_BASENAME)

SURV_GLOBAL_SETTINGS = $(SURV_SCRIPTS)/GLOBALS.R

### Submodel 1
SURV_SIMULATION_SETTINGS_ONE = $(SURV_RDS)/submodel-one-simulation-settings.rds
$(SURV_SIMULATION_SETTINGS_ONE) : $(SURV_SCRIPTS)/simulation-settings-submodel-one.R $(SURV_GLOBAL_SETTINGS)
	$(RSCRIPT) $<

SURV_SUBMODEL_ONE_SIMULATED_DATA = $(SURV_RDS)/submodel-one-simulated-data.rds
$(SURV_SUBMODEL_ONE_SIMULATED_DATA) : $(SURV_SCRIPTS)/simulate-data-submodel-one.R $(SURV_SIMULATION_SETTINGS_ONE) 
	$(RSCRIPT) $<

SURV_SUMODEL_ONE = $(SURV_MODELS)/submodel-one.stan

SURV_SUBMODEL_ONE_OUTPUT = $(SURV_RDS)/submodel-one-output.rds
$(SURV_SUBMODEL_ONE_OUTPUT) : $(SURV_SCRIPTS)/fit-submodel-one.R $(SURV_SUBMODEL_ONE_SIMULATED_DATA) $(SURV_SUMODEL_ONE)
	$(RSCRIPT) $<

SURV_ALL_SUBMODEL_ONE_INPUTS = $(SURV_SIMULATION_SETTINGS_ONE) $(SURV_SUBMODEL_ONE_SIMULATED_DATA) $(SURV_SUBMODEL_ONE_OUTPUT)

SURV_SUBMODEL_ONE_POSTERIOR_PLOT = $(SURV_PLOTS)/submodel-one-posterior.pdf
$(SURV_SUBMODEL_ONE_POSTERIOR_PLOT) : $(SURV_SCRIPTS)/plot-submodel-one.R $(PLOT_SETTINGS) $(MCMC_UTIL) $(SURV_ALL_SUBMODEL_ONE_INPUTS)
	$(RSCRIPT) $<

ALL_PLOTS += $(SURV_SUBMODEL_ONE_POSTERIOR_PLOT)

SURV_SUBMODEL_ONE_DIAG_PLOT = $(SURV_PLOTS)/stage-one-submodel-one-diags.png
$(SURV_SUBMODEL_ONE_DIAG_PLOT) : $(SURV_SCRIPTS)/diagnose-submodel-one-plots.R $(PLOT_SETTINGS) $(MCMC_UTIL) $(SURV_SUBMODEL_ONE_OUTPUT)
	$(RSCRIPT) $<

ALL_PLOTS += $(SURV_SUBMODEL_ONE_DIAG_PLOT)

SURV_SUBMODEL_ONE_DIAG_TABLE = $(SURV_TEX)/0080-submodel-one-numeric-diags.tex
$(SURV_SUBMODEL_ONE_DIAG_TABLE) : $(SURV_SCRIPTS)/diagnose-submodel-one-tables.R $(MCMC_UTIL) $(SURV_SUBMODEL_ONE_OUTPUT)
	$(RSCRIPT) $<

TEX_FILES += $(SURV_SUBMODEL_ONE_DIAG_TABLE)

### Submodel 3
SURV_SIMULATION_SETTINGS_THREE = $(SURV_RDS)/submodel-three-simulation-settings.rds
$(SURV_SIMULATION_SETTINGS_THREE) : $(SURV_SCRIPTS)/simulation-settings-submodel-three.R $(SURV_GLOBAL_SETTINGS)
	$(RSCRIPT) $<

SURV_SUBMODEL_THREE_SIMULATED_DATA = $(SURV_RDS)/submodel-three-simulated-data.rds
$(SURV_SUBMODEL_THREE_SIMULATED_DATA) : $(SURV_SCRIPTS)/simulate-data-submodel-three.R $(SURV_SIMULATION_SETTINGS_THREE) $(SURV_SIMULATION_SETTINGS_ONE)
	$(RSCRIPT) $<

SURV_SUMODEL_THREE = $(SURV_MODELS)/submodel-three.stan
SURV_SUBMODEL_THREE_OUTPUT = $(SURV_RDS)/submodel-three-output.rds
$(SURV_SUBMODEL_THREE_OUTPUT) : $(SURV_SCRIPTS)/fit-submodel-three.R $(SURV_SUBMODEL_THREE_SIMULATED_DATA) $(SURV_SUMODEL_THREE)
	$(RSCRIPT) $<

SURV_ALL_SUBMODEL_THREE_INPUTS = $(SURV_SIMULATION_SETTINGS_THREE) $(SURV_SUBMODEL_THREE_SIMULATED_DATA) $(SURV_SUBMODEL_THREE_OUTPUT)

SURV_SUBMODEL_THREE_POSTERIOR_PLOT = $(SURV_PLOTS)/submodel-three-posterior.pdf
$(SURV_SUBMODEL_THREE_POSTERIOR_PLOT) : $(SURV_SCRIPTS)/plot-submodel-three.R $(PLOT_SETTINGS) $(MCMC_UTIL) $(SURV_ALL_SUBMODEL_THREE_INPUTS)
	$(RSCRIPT) $<

ALL_PLOTS += $(SURV_SUBMODEL_THREE_POSTERIOR_PLOT)

SURV_SUBMODEL_THREE_DIAG_PLOT = $(SURV_PLOTS)/stage-one-submodel-three-diags.png
$(SURV_SUBMODEL_THREE_DIAG_PLOT) : $(SURV_SCRIPTS)/diagnose-submodel-three-plots.R $(PLOT_SETTINGS) $(MCMC_UTIL) $(SURV_SUBMODEL_THREE_OUTPUT)
	$(RSCRIPT) $<

ALL_PLOTS += $(SURV_SUBMODEL_THREE_DIAG_PLOT)

SURV_SUBMODEL_THREE_DIAG_TABLE = $(SURV_TEX)/0081-submodel-three-numeric-diags.tex
$(SURV_SUBMODEL_THREE_DIAG_TABLE) : $(SURV_SCRIPTS)/diagnose-submodel-three-tables.R $(MCMC_UTIL) $(SURV_SUBMODEL_THREE_OUTPUT)
	$(RSCRIPT) $<

TEX_FILES += $(SURV_SUBMODEL_THREE_DIAG_TABLE)

### Stage 2
SURV_SUBMODEL_TWO_SIMULATED_DATA = $(SURV_RDS)/submodel-two-simulated-data.rds
$(SURV_SUBMODEL_TWO_SIMULATED_DATA) : $(SURV_SCRIPTS)/simulate-data-submodel-two.R $(SURV_GLOBAL_SETTINGS) $(SURV_SIMULATION_SETTINGS_ONE)
	$(RSCRIPT) $<

### Stage 2 - sampling the melded posterior
SURV_SUBMODEL_TWO_PHI_STEP = $(SURV_MODELS)/submodel-two-phi-step.stan
SURV_SUBMODEL_TWO_PSI_STEP = $(SURV_MODELS)/submodel-two-psi-step.stan
SURV_STAGE_TWO_STAN_MODELS = $(SURV_SUBMODEL_TWO_PHI_STEP) $(SURV_SUBMODEL_TWO_PSI_STEP)
SURV_STAGE_TWO_PHI_12_SAMPLES = $(SURV_RDS)/stage-two-phi-12-samples.rds
SURV_STAGE_TWO_PHI_23_SAMPLES = $(SURV_RDS)/stage-two-phi-23-samples.rds
SURV_STAGE_TWO_PSI_2_SAMPLES = $(SURV_RDS)/stage-two-psi-2-samples.rds
SURV_STAGE_TWO_PSI_1_INDICES = $(SURV_RDS)/stage-two-psi-1-indices.rds
SURV_STAGE_TWO_PSI_3_INDICES = $(SURV_RDS)/stage-two-psi-3-indices.rds

$(SURV_STAGE_TWO_PHI_12_SAMPLES) : $(SURV_SCRIPTS)/fit-stage-two.R $(SURV_SUBMODEL_ONE_OUTPUT) $(SURV_SUBMODEL_THREE_OUTPUT) $(SURV_SUBMODEL_TWO_SIMULATED_DATA) $(SURV_STAGE_TWO_STAN_MODELS) $(SURV_GLOBAL_SETTINGS)
	$(RSCRIPT) $<

$(SURV_STAGE_TWO_PHI_23_SAMPLES) : $(SURV_STAGE_TWO_PHI_12_SAMPLES)

$(SURV_STAGE_TWO_PSI_2_SAMPLES) : $(SURV_STAGE_TWO_PHI_12_SAMPLES)

$(SURV_STAGE_TWO_PSI_1_INDICES)	: $(SURV_STAGE_TWO_PHI_12_SAMPLES)

$(SURV_STAGE_TWO_PSI_3_INDICES)	: $(SURV_STAGE_TWO_PHI_12_SAMPLES)

### Stage 2 - diagnostics
#### This form of target definition + rule will make --dry-run report 
#### incorrectly, but make will only run the relevant file once.
SURV_STAGE_TWO_DIAG_PLOTS = $(wildcard plots/surv-example/stage-two*-diags.png)
$(SURV_STAGE_TWO_DIAG_PLOTS) : $(SURV_SCRIPTS)/diagnose-stage-two-plots.R $(PLOT_SETTINGS) $(SURV_STAGE_TWO_PHI_12_SAMPLES) $(SURV_STAGE_TWO_PHI_23_SAMPLES) $(SURV_STAGE_TWO_PSI_2_SAMPLES)
	$(RSCRIPT) $<

ALL_PLOTS += $(SURV_STAGE_TWO_DIAG_PLOTS)

#### The tables will be automatically picked up by TEX_FILES, under the 
#### assumption that they exist first.
SURV_STAGE_TWO_DIAG_TABLES = $(wildcard tex-input/surv-example/009*-stage-two-**.tex) 
$(SURV_STAGE_TWO_DIAG_TABLES) : $(SURV_SCRIPTS)/diagnose-stage-two-tables.R $(SURV_STAGE_TWO_PHI_12_SAMPLES) $(SURV_STAGE_TWO_PHI_23_SAMPLES) $(SURV_STAGE_TWO_PSI_2_SAMPLES) $(MCMC_UTIL)
	$(RSCRIPT) $<

TEX_FILES += $(SURV_STAGE_TWO_DIAG_TABLES)

### Compositional step / plots
SURV_EXAMPLE_PHI_12_SRINKAGE_PLOT = $(SURV_PLOTS)/phi-12-inter-stage-comparison.pdf
$(SURV_EXAMPLE_PHI_12_SRINKAGE_PLOT) : $(SURV_SCRIPTS)/plot-event-time-shrinkage.R $(SURV_STAGE_TWO_PHI_12_SAMPLES) $(SURV_SUBMODEL_ONE_OUTPUT) $(SURV_SIMULATION_SETTINGS_ONE)
	$(RSCRIPT) $<

ALL_PLOTS += $(SURV_EXAMPLE_PHI_12_SRINKAGE_PLOT)

################################################################################
# knitr is becoming more picky about encoding, specify UTF-8 input
$(WRITEUP) : $(wildcard *.rmd) $(TEX_FILES) $(ALL_PLOTS) $(OWLS_DATA) $(BIBLIOGRAPHY)
	$(RSCRIPT) -e "rmarkdown::render(input = Sys.glob('*.rmd'), encoding = 'UTF-8')"