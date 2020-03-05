RSCRIPT = Rscript
PLOT_SETTINGS = scripts/common/plot-settings.R
TEX_FILES = $(wildcard tex-input/*.tex) \
	$(wildcard tex-input/*/*.tex) \
	$(wildcard tex-input/*/*/*.tex)

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

## Owls example
OWLS_BASENAME = owls-example
OWLS_DATA = $(wildcard rds/owls-example/*-data.rds)

$(OWLS_DATA) : $(SCRIPTS)/$(OWLS_BASENAME)/load-and-write-data.R
	$(RSCRIPT) $<

FECUNDITY_SUBPOSTERIOR = $(RDS)/$(OWLS_BASENAME)/fecundity-subposterior-samples.rds 
$(FECUNDITY_SUBPOSTERIOR) : $(SCRIPTS)/$(OWLS_BASENAME)/fit-fecundity.R $(RDS)/$(OWLS_BASENAME)/models/fecundity-model.stan $(FECUNDITY_DATA)
	$(RSCRIPT) $< 

ORIG_IPM_RESULTS = $(RDS)/$(OWLS_BASENAME)/original-ipm-results.rds
$(ORIG_IPM_RESULTS) : $(SCRIPTS)/$(OWLS_BASENAME)/fit-original-ipm.R $(SCRIPTS)/$(OWLS_BASENAME)/models/original-ipm.bug $(OWLS_DATA)
	$(RSCRIPT) $<

################################################################################

# knitr is becoming more picky about encoding, specify UTF-8 input
$(WRITEUP) : $(wildcard *.rmd) $(TEX_FILES) $(ALL_PLOTS) $(OWLS_DATA)
	$(RSCRIPT) -e "rmarkdown::render(input = Sys.glob('*.rmd'), encoding = 'UTF-8')"