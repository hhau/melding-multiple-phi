## Combining chains of Bayesian models with Markov melding

A challenge for practitioners of Bayesian inference is specifying a single model that incorporates multiple heterogeneous data sources or types.
It may be easier instead to specify distinct submodels for each source of data, then "join" the submodels together.
When all submodels are related via the same common quantity, which could be a parameter or deterministic function thereof, the submodels can be joined using _Markov melding_ (Goudie et al., 2019).
However, it is unclear how to join submodels when they are linked in a chain structure, i.e. where neighbouring submodels contain common quantities.

We propose _chained Markov melding_ as a generic method for combining chains of submodels into a joint model.
Part of the submodel joining process involves combining, or pooling, the submodel-specific prior distributions for the common quantities.
Standard techniques for pooling distributions, such as linear or logarithmic pooling, are not applicable.
We propose extensions to these standard pooling methods that can be applied to chains of submodels.
We also propose a posterior sampling algorithm that uses the chain structure to incorporate information contained in the submodels in multiple stages, potentially in parallel.

We demonstrate our methodology using two examples.
The first example is an integrated population model used in conservation ecology, where multiple data are required to accurately estimate population immigration and reproduction rates.
We also consider a joint longitudinal and time-to-event model with the added complication that event times are derived from a separate submodel, and hence are uncertain.
We show that failing to incorporate uncertainty in the event times results in overconfident and biased posterior inference.
Chained Markov melding is a conceptually appealing approach to integrating submodels in these settings.