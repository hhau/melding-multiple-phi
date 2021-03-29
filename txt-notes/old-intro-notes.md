# sec 1.1.2

- what are some good examples of methodology papers that use simulated data (ideally data inspired by a real problem?)

- In our second example we consider modelling the time to respiratory failure amongst intensive care unit patients as a function of time-invariant (_baseline_) and time-varying (_longitudinal_) covariates.
  - what is the longitudinal covariate.


- Integrating longitudinal and standard time-to-event data that are available concurrently is popular, and is often referred to as a _joint model_ [@rizopoulos_joint_2012].


In our second example we wish to synthesise data about a longitudinal process with data relating to a time-to-event process, where the event time is determined by a second longitudinal process crossing a pre-specified threshold.


Complicated longitudinal processes require flexible regression models, often including many subject specific parameters, and can be computationally challenging to fit jointly with time to event data.
This motivates a two-stage approach, where the longitudinal submodel is fit first and the fitted values are then used in the survival submodel [@tsiatis_joint_2004; @mauff_joint_2020].

In some settings the time at which the relevant event occurs is inherently uncertain.
This can be due to the definition of the event, for example the time when a measurable but noisy quantity crosses a prescribed threshold, or because the event of interest is inherently unobservable and must be inferred from other information.
_First hitting time_, or _threshold regression_ models [@lee_threshold_2006] are used in threshold crossing settings, and applications include kidney failure [@diggle_real-time_2015], dementia diagnosis [@hashemi_latent_2003], and disease transmission [@yu_semi-parametric_2009].
These methods assume that relevant measurements are realisations from a specific, possibly latent, stochastic process, where the distribution of the threshold crossing time is known.
A conceptually similar idea is presented in @lu_using_1993 and @sweeting_estimating_2010; both use hierarchical regression models instead of stochastic processes, and define the event time as when the estimated regression curve crosses the threshold.
Another setting in which inherently uncertain event times arise is record linkage [@harron_methodological_2015].
Such settings are unique in that the true event time is assumed to be one of a finite set of possible times, with one possible event time arising from each matching record.
Integrating this uncertainty into standard survival models is challenging but possible [@wang_integrative_2020], however these methods do not consider event times that are continuous[^wang].

[^wang]: @wang_integrative_2020 assume that the response is a finite mixture, where each mixture component is a point mass at the possible event times arising due to record linkage. Uncertainty is then incorporated by estimating the mixture weights.

It seems desirable to integrate both uncertain event times and longitudinal data in a survival analysis.
Doing so in a standard joint model is challenging, as incorporating uncertainty in the response requires complex methodology [@giganti_accounting_2020; @oh_considerations_2018; @oh_raking_2021] which do not explicitly include time-varying covariates.
Similarly, incorporating time-varying covariates into a threshold regression model is nontrivial, and assumes the relevant covariates satisfy the Markov property [@lee_threshold_2010].
 <!-- - Another is that integrating time dependent covariates into the FHT/TR models is non straightforward -- as far as I can tell the Markov assumption in Lee and Whitmore is equivalent to restricting the form of the longitudinal submodel to linear forms (non linear longitudinal regression models are non-Markovian / joint model that depend on the integral of the longitudinal model are definitely non-Markovian. This is discussed a bit in @lee_threshold_2010)
    - it would help if our longitudinal model was non-linear so that we could have this benefit -- maybe making the regression quadratic.
  - One motivation comes from the respiratory failure idea and its similarities to the degradation models. (what would the other, longitudinal covariate be?)
    - This would be strengthened by using non-linear regression / splines for the first submodel. -->
Chained Markov melding offers a conceptually straightforward approach to incorporating uncertain event times.
Specifically, we consider the event time as a submodel-derived quantity from a hierarchical regression model akin to @lu_using_1993. 
We call this submodel the _uncertain event time_ submodel and denote it $\pd_{1}(\phi_{1 \cap 2}, \psi_{1}, Y_{1})$, with $\phi_{1 \cap 2}$ the event time.
The survival submodel $\pd_{2}(\phi_{1 \cap 2}, \phi_{2 \cap 3}, \psi_{2}, Y_{2})$ uses $\phi_{1 \cap 2}$, the common quantity, as the response. 
We treat the longitudinal submodel, $\pd_{3}(\phi_{2 \cap 3}, \psi_{3}, Y_{3})$, separately from the survival submodel, as is common in two-stage survival modelling, and denote the subject-specific parameters that also appear in the survival model as $\phi_{2 \cap 3}$.
The high level submodel relationships are displayed as a DAG in Figure \ref{fig:surv-simple-dag}.

# sec 5.3

Standard survival analyses consider the event times and censoring indicators as known, fixed quantities.
However, in some settings the survival times are contaminated with nontrivial measurement error [@gu_semiparametric_2015; @meier_discrete_2003; @oh_considerations_2018; @snapinn_survival_1998;@wang_integrative_2020], or are derived from a, possibly highly complex, model for some other measurable quantity of the individual [see @lee_threshold_2006 for a review; and @szczesniak_dynamic_2020 for recent application].
Integrating the output of such models into standard parametric or Cox proportional hazards survival models is challenging, but we believe it is made conceptually and computational easier by our chained melding framework.

The integration of complex regression models into survival models is routine in the Joint modelling framework [@rizopoulos_joint_2012].
Such model complexity is often warranted due to the wide variety of behaviours of the longitudinal process. 
The associated computational cost motivates two stage estimation procedures [@ye_semiparametric_2008; @mauff_joint_2020] which reduce the computation required.
Such two stage processes are, in effect, considering a separate submodel for the longitudinal data, and we do the same in our example.
We argue that considering the models for the event times and longitudinal process separately is not an unreasonable modelling choice, and choosing to unite them in the survival model is a natural application of the chained melding method.

To highlight the chained melding process, and not obscure it behind individual submodel complexity, we use simpler versions of the threshold crossing, survival, and longitudinal submodels.
We use a subject specific linear model with constant threshold, a parametric Weibull model, and a distinct linear subject specific linear model for the respective submodels.
Data are simulated in a manner that is realistic, but avoids multiple types of censoring in the survival model and missing data.