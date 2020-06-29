# Time: 1330 - 15:45

## Multiple phi

- coherent needs to not be used, doesn't make sense.


### Pooling

- Can we get approximately the same results, with a lot more intuition, by specifying a conditional extension, such as:
  \pd_{1}(\phi_{2 \cap 3} \mid \phi_{1 \cap 2}) = \text{N}(0, 100^2), which makes some of the maths easier

- frame this as "this is a subjective choice, and difficult to provide consistent guidance." "Hope to provide some intuition for performing this operation", but, for the following non-exhaustive list of reasons, it is almost impossible to provide all purpose guidance for choosing these parameters. 

- for Plot
  - aim is to give readers a little intuition for what is happening for various values of $\lambda$
  - Two one-page plots with different values of $\lambda$ (one page per distribution type?)

- set Eqs (10), (11) back to just using parameters 

### Things for second example

- think of Bayesian models that are do-able, but nontrivial.
- joint + longitudinal models (again?)
- Spatial models?
- Quantile regression.
- Look at Aki's posterior DB for models that could be extended / combined with something?
  - We may have to build the second/third models.

- look at biological part of textbook Rob emailed (Applied Bayesian case studies) 


### possible PF ratio question

- For respiratory failure, this is a much more biologically relevant endpoint (ICU admission as other problems as well)
- The key measures are 
  - pAO_{2} - partial pressure of oxygen?
  - FiO_{2} - something else of oxygen
  - sO_{2} - 
  - What are these measured in (units / percentage), when and where are they measured (ICU / general ward, sometimes in general ward if clinician is worried about respiratory status of patient)
  - Some of these measures are easier to obtain than others. hence, it may be worth modelling them (trying to model the harder to obtain measures as a function of the easier ones, as well as introducing appropriate uncertainty.)
  - regression model a la: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3776410/pdf/nihms505475.pdf

- Ideally, we would have these measurements over time. Then, we could ask, what is the distribution of the time to first (probable) crossing of some threshold. 
- This distribution could then be fed into a survival model, as the survival status.
  - read: https://data.princeton.edu/pop509/ParametricSurvival.pdf
  - https://arxiv.org/pdf/2002.09633.pdf
- For the third model, we could have covariates in the survival model, and have some of them missing (Bayesian imputation)

- This will involve learning about survival models + bayesian missing data / multiple imputation (which is fine, might even be quite interesting)


- make notes on: file:///Users/amanderson/Downloads/env.571.pdf
  - multiple scale

## meeting links:

- build regression model for https://www.intensive.org/epic2/Documents/Estimation%20of%20PO2%20and%20FiO2.pdf
- missing data model from http://www.bias-project.org.uk/Missing2012/Lectures.pdf


- This seems important in ARDS (acute respiratory distress syndrome, and they may have relevant, available data?)
  - paper: https://www.nejm.org/doi/full/10.1056/NEJM200005043421801
  - paper: https://bmjopenrespres.bmj.com/content/6/1/e000420
  - data: http://ardsnet.org/ via https://biolincc.nhlbi.nih.gov/home/