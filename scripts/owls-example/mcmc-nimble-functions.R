library(nimble)

# I need to load the samples in here in this file because they
# now form part of the proposal

# Read in previous stage samples
fecunditiy_submodel_samples <- readRDS("rds/owls-example/fecundity-subposterior-samples.rds")
capture_recapture_submodel_samples <- readRDS("rds/owls-example/capture-recapture-subposterior-samples.rds")

# reflow into iterations x parameters matrix
fecunditiy_proposal <- as.vector(fecunditiy_submodel_samples)
capture_recapture_proposal <- array(
  as.vector(capture_recapture_submodel_samples),
  dim = c(2e4 * 6, 85), 
  dimnames = list(
    NULL,
    dimnames(capture_recapture_submodel_samples)[[3]]
  )
)
# thin this to just the pars of interest
# find out how to save the proposed indices and store them later
capture_recapture_proposal <- capture_recapture_proposal[, c("v[1]", "v[2]")]

# phi_12 = 'v[1]', 'v[2]'
phi_one_two_proposal <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    internal_capture_recapture_proposal <- control$capture_recapture_proposal
  },
  run = function() {
    model_lp_initial <- getLogProb(model, calcNodes)
    proposal_index <- rcat(
      n = 1,
      prob = nimRep(1, each = 1.2e5)
    )
    proposal <- nimC(
      internal_capture_recapture_proposal[proposal_index, 1],
      internal_capture_recapture_proposal[proposal_index, 2]
    )

    model[[target]] <<- proposal 
    model_lp_proposed <- calculate(model, calcNodes)
    log_MH_ratio <- model_lp_proposed - model_lp_initial
    u <- runif(1, 0, 1)
    if (u < exp(log_MH_ratio)) {
      jump <- TRUE
    } else { 
      jump <- FALSE
    }
    # keep the model and mvSaved objects consistent
    if (jump) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else {
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list(reset = function() {})
)

# phi_23 = 'fec'
phi_two_three_proposal <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    internal_fecunditiy_proposal <- control$fecunditiy_proposal
  },
  run = function() {
    model_lp_initial <- getLogProb(model, calcNodes)
    proposal <- internal_fecunditiy_proposal[
      rcat(
        n = 1,
        prob = nimRep(1, each = 3e4)
      )
    ]
    model[[target]] <<- proposal
    model_lp_proposed <- calculate(model, calcNodes)
    log_MH_ratio <- model_lp_proposed - model_lp_initial
    u <- runif(1, 0, 1)
    if (u < exp(log_MH_ratio)) {
      jump <- TRUE
    } else { 
      jump <- FALSE
    }
    # keep the model and mvSaved objects consistent
    if (jump) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else {
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list(reset = function() {})
)

# set some initial values
count_model_inits_function  <- function() {
  list(
    v = c(
      capture_recapture_proposal[sample(1 : 1.2e5, size = 1),],
      rnorm(n = 4, mean = 0, sd = 0.1)
    ),
    fec = fecunditiy_proposal[
      sample(1 : 3e4, size = 1)
    ],
    N1 = array(
      data = sample(1 : 10, size = 26, replace = TRUE),
      dim = c(26),
      dimnames = list(
        sprintf("N1[%d]", 1 : 26)
      )
    ),
    NadSurv = array(
      data = sample(1 : 10, size = 26, replace = TRUE),
      dim = c(26),
      dimnames = list(
        sprintf("NadSurv[%d]", 1 : 26)
      )
    ),
    Nadimm = array(
      data = sample(1 : 10, size = 26, replace = TRUE),
      dim = c(26),
      dimnames = list(
        sprintf("Nadimm[%d]", 1 : 26)
      )
    )
  )
}

# the default sampler doesn't obey support constraints,
# so we need to constrain the lower bound of the slice
# sampler to be positive. 
# In line with the above, I have put the adaptive functions in the 
# 'methods' list, but this may be incorrect, if the internal samplers
# are written with a different structure
positive_discrete_slice_sampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    adaptive <- if (!is.null(control$adaptive)) 
        control$adaptive
    else TRUE
    adaptInterval <- if (!is.null(control$adaptInterval)) 
        control$adaptInterval
    else 200
    width <- if (!is.null(control$sliceWidth)) 
        control$sliceWidth
    else 1
    maxSteps <- if (!is.null(control$sliceMaxSteps)) 
        control$sliceMaxSteps
    else 100
    maxContractions <- if (!is.null(control$maxContractions)) 
        control$maxContractions
    else 1000
    maxContractionsWarning <- if (!is.null(control$maxContractionsWarning)) 
        control$maxContractionsWarning
    else TRUE
    eps <- 1e-15
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    widthOriginal <- width
    timesRan <- 0
    timesAdapted <- 0
    sumJumps <- 0
    discrete <- model$isDiscrete(target)
    if (length(targetAsScalar) > 1) 
        stop("cannot use slice sampler on more than one target node")
  },
  run = function() {
    u <- getLogProb(model, calcNodes) - rexp(1, 1)
    x0 <- model[[target]]
    L <- x0 - runif(1, 0, 1) * width
    if (L < 0) {
        L <- 0
    }
    R <- L + width
    maxStepsL <- floor(runif(1, 0, 1) * maxSteps)
    maxStepsR <- maxSteps - 1 - maxStepsL
    lp <- setAndCalculateTarget(L)
    while (maxStepsL > 0 & !is.nan(lp) & lp >= u) {
        L <- L - width
        if (L < 0) {
            L <- 0
        }
        lp <- setAndCalculateTarget(L)
        maxStepsL <- maxStepsL - 1
    }
    lp <- setAndCalculateTarget(R)
    while (maxStepsR > 0 & !is.nan(lp) & lp >= u) {
        R <- R + width
        lp <- setAndCalculateTarget(R)
        maxStepsR <- maxStepsR - 1
    }
    x1 <- L + runif(1, 0, 1) * (R - L)
    lp <- setAndCalculateTarget(x1)
    numContractions <- 0
    while ((is.nan(lp) | lp < u) & (R - L)/(abs(R) + abs(L) + 
        eps) > eps & numContractions < maxContractions) {
        if (x1 < x0) {
            L <- x1
        }
        else {
            R <- x1
        }
        x1 <- L + runif(1, 0, 1) * (R - L)
        lp <- setAndCalculateTarget(x1)
        numContractions <- numContractions + 1
    }
    if ((R - L)/(abs(R) + abs(L) + eps) <= eps | numContractions == 
        maxContractions) {
        if (maxContractionsWarning) 
            cat("Warning: slice sampler reached maximum number of contractions for '", 
                target, "'. Current parameter value is ", x0, 
                ".\n")
        nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, 
            logProb = TRUE)
    }
    else {
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, 
            logProb = TRUE)
        jumpDist <- abs(x1 - x0)
        if (adaptive) 
            adaptiveProcedure(jumpDist)
    }
  },
  methods = list(
    setAndCalculateTarget = function(value = double()) {
      if (discrete) 
          value <- floor(value)
      model[[target]] <<- value
      lp <- calculate(model, target)
      if (lp == -Inf) 
          return(-Inf)
      lp <- lp + calculate(model, calcNodesNoSelf)
      returnType(double())
      return(lp)
    } ,
    adaptiveProcedure = function(jumpDist = double()) {
      timesRan <<- timesRan + 1
      sumJumps <<- sumJumps + jumpDist
      if (timesRan%%adaptInterval == 0) {
          adaptFactor <- (3/4)^timesAdapted
          meanJump <- sumJumps/timesRan
          width <<- width + (2 * meanJump - width) * adaptFactor
          timesAdapted <<- timesAdapted + 1
          timesRan <<- 0
          sumJumps <<- 0
      }
    },
    reset = function() {
        width <<- widthOriginal
        timesRan <<- 0
        timesAdapted <<- 0
        sumJumps <<- 0
    }
  ),
  where = getLoadingNamespace()
)