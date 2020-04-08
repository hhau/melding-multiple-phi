# Effective Sample Size per Second is misleading (at best)

## Structure

- Thesis
- Example: My recent experience with Nimble
-

- Bayesian statistics has an obsession with speed.
- in some settings this is sensible

- I would like to argue that comparing posterior estimation algorithms in terms of ESS/second performance is not only misleading, but completely misses the point of algorithmic comparison

- Namely, algorithms have audiences, sometimes they are explicitly stated, but many times not.

- Developing new statistical algorithms requires demonstrating that your proposed algorithm achieves the same result as an accepted algorithm in empirical examples.
- We conveniently ignore the parts of this process that _actually_ take the most time.
  - implementation / learning a new framework because one needs to compare your implementation to it
  - debugging / verifying said new implementation
    - Can we ever be sure that stochastic algorithms are free of subtle biases?
    - Inscrutable error messages in new languages.
      - Nimble error message.
  - compilation. (My current experience with Stan/Nimble is that in cases where the examples are simple enough to be compared across implementations, the compilation time dominates)

1. The time taken to learn and implement the model of interest in the programming language / modelling framework of choice
  - No one includes how comprehensible the documentation is in their speed tests. How is it far to compare the tensorflow doc with the Stan doc with the BUGS manual?
  The first is vast, and vastly varies in quality and comprehensibility (parts of the doc in `tf.contrib` are just flat out wrong/outdated, and only occasionally bare some relation to the academic work on which they are based).
  Stan's doc is also vast, but perhaps less so than TF's. It was incredibly difficult to use when it was in one big pdf, and is still not searchable in any kind of usable way (bookdown's search is terrible), and I'm never sure if I should be looking which of the three sections I should be looking in.

2. researcher DFs mean I can always make the speed test say whatever I want it to.