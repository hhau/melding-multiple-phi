Rob meetings notes dump

1. wsre

  - Make changes that are suggested by Rob in response-modified.md
  - rework introduction, argument, framing
    - Setting off point is evidence synthesis, with complicated, deterministic link functions
      - Emphasise the need for both prior marginals to be non-flat and unknown (strange)
    - emphasise that sample re-use / two stage sampling is common, and our choice to use it is slightly synthetic (our examples are simple and we could just sample the joint), but we choose not to, because ?
    - when is it difficult to sample the joint model
    - 'toy' nature of other example (HIV)
    - simulation study to show it works?
      - Unclear how much this is beneficial
    - Need another, much more real, example where this is a problem. (where we have a complex, deterministic, link function, and two strange prior marginals)
      - if i knew how to locate or invent one of these, I think I would have by now.
      - I want a model joining example, where we're joining on some derivable quantity in both models?

2. multiple phi

  - where are we heading in general?
    - Read JRSSC, not there
    - Probably towards Bayesian Analysis.
      - Figure out specific plan for getting the chapter to a state where we could submit to BA
      - Figure out how general the methodology should be written up in the text (M = M, or M = 3)
        - make clear that we are only talking about chains of models, not graphs
    - Where are we going to get some real, novel data to do the survival analysis example on
      - First need to figure out if this works with simulated data
      - Doing something like SBC on the whole melding process is one way to show that it _could_ work?
    - Why is this useful
    - Why should people care?
  - surv example
    - read other literature on uncertain event times + summarise
  - 'not readily available' -> not easily analytically tractable
  - Priors:
    - the aim is to help an unfamiliar, semi-applied researcher figure out what values of lambda are reasonable
    - Properties of approaches, are they externally Bayesian? Figure out
    - Rewrite Section 3.2.1 (linear pooling) to not mention bad idea first
    - Section 3.2.2 Include a more general introduction to Log pooling
    - 3.2.3 Overview of area that makes people feel confident that they understand the choices of lambda?
    - Figure 2:
      - Lambda still sum to 1
      - linear, evenly spaced values of $\lambda_{1}$.
      - Gaussian only
      - 1-D marginals in axes, 2-d marginal on the plot, pooled prior on the plot
    - Cut all of Figure 3 + discussion
  - simulated sruv example
    - censor at last observation time
  - Add signal to ensure that coefficients in survival model are not zero?
    - Maybe not, I think doing something like SBC might be more convincing anyway?
  - What is it we are trying to show? What is improved by adding the uncertain event times?
    - 'Better' quantification of uncertainty? do we even have a way of knowing what it _should_ be? 
  - Real example is to show?
    - Look up lumley third author paper about uncertain event times
    - HIV => AIDS transition event, modelling this from cell counts?
      - Peter?
    - Other examples?
      - Talk to Ed?

3. phd general

  - how do I make WSRE something I can confidently defend in the viva?
  - How do I ensure the novel contribution of the second paper is enough to get a thesis?