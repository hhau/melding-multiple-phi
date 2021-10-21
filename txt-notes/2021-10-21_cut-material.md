## Coherency of the chained melded model

- _I think a few of these results are wrong. Maybe this should just be a paragraph?_
    - 3 ii) is definitely wrong.
    - 3 i) doesn't required independence
- _Some discussion of the M = 3 case in text, with reference to the general case (if I can figure out what that even means) in the appendix_
- _I think I can also put together a discussion of what it means to be externally Bayesian in this context now?_.
    - _no_

A reasonable requirement for a modular inference method is that the final posterior distribution should not, at least theoretically, depend on the order in which data are observed or integrated into the model.
In the context of belief distributions, @bissiri_general_2016 call this property 'coherence', which we will co-opt.
Specifically, in the $\Nm = 3$ case, it seems desirable that the chained melded model be the same if either $\pd_{1}$ or $\pd_{3}$ are integrated with submodel $\pd_{2}$ first.
We will show that the form proposed in Equation \eqref{eqn:melded-model-full} satisfies this property, whilst the model produced by applying the melding method of @goudie_joining_2019 is, in general, sensitive to the order of integration.
    <!--
        - should we also check p(13)2?
            - (13)2 is the same as original melding with \phi = (\phi_{1 \cap 2}, \phi_{2 \cap 3}),
            - but the strategy doesn't generalise to higher $M$ (M = 4 immediately disproves)

        - say we had 9 models, and the first three formed \pd_{1}, the second three formed \pd_{2} and the remaining $\pd_{3}$. Would we apply melding here twice? No, we would think about this as chained melding with M = 9 models.
    -->

- do we need to know all the submodels ahead of time a la Section 4.4.3 of @jacob_better_2017-1?

### Applying Markov melding twice

The general thing I want to be able to say is


The idea is that they can only possibly be equal if there is prior independence between the common quantities in each submodel.
If this is missing, then for any, non product-of-experts, choice of weight functions

Denote the original melding operator with $\circledast$.
Its output is
\input{tex-input/noncommutativity/0005-def-usual-melded-model.tex}
where $\pd_{\text{pool}}^{12}(\phi_{1 \cap 2}) = g^{12}(\pd_{1}(\phi_{1 \cap 2}), \pd_{2}(\phi_{1 \cap 2}))$ for some pooling function $g^{12}$.
We denote the parameter space of the output as $\boldsymbol{\Theta}_{12} = (\phi_{1 \cap 2}, \phi_{2 \cap 3}, \psi_{1}, \psi_{2}, Y_{1}, Y_{2})$, so that any prior marginal distribution of $\pd_{\text{meld}}^{12}$ can be derived by integrating out the irrelevant components of $\boldsymbol{\Theta}_{12}$.
For example,
\input{tex-input/noncommutativity/0006-example-melded-marginal-definition.tex}
where $\boldsymbol{\Theta}_{12} \setminus \phi_{2 \cap 3}$ is $\boldsymbol{\Theta}_{12}$ without $\phi_{2 \cap 3}$.

To integrate third submodel, we apply the original operator to $\pd_{\text{meld}}^{12}$ and $\pd_{3}$
\input{tex-input/noncommutativity/0007-iterated-application-melding.tex}
so that the parentheses in the superscript of $\pd_{\text{meld}}^{(12)3}$ indicate the order in which the submodels are melded.
As before we define $\pd_{\text{pool}}^{(12)3}(\phi_{2 \cap 3}) = g^{(12)3}(\pd_{\text{meld}}^{12}(\phi_{2 \cap 3}), \pd_{3}(\phi_{2 \cap 3}))$, but for a potentially different choice of pooling function $g^{(12)3}$.

It will be convenient to consider the following, expanded form of Equation \eqref{eqn:iterated-application-melding-two}
\input{tex-input/noncommutativity/0013-expanded-double-melded-model.tex}
and the equivalent expression for $\pd_{\text{meld}}^{1(23)}$, which we derive by careful inspection of the superscripts
\input{tex-input/noncommutativity/0018-symmetric-expanded-double-melded-model.tex}

### Does melding twice produce the same model as the chained melded model?

For Equations \eqref{eqn:expanded-double-melded-model} and \eqref{eqn:symmetric-expanded-double-melded-model} to be equal to the model defined in Equation \eqref{eqn:melded-model-full}, the following equalities must hold:
\input{tex-input/noncommutativity/0040-coinciding-equalities.tex}
It is clear from Equation \eqref{eqn:coinciding-equalities-one} that one necessary condition is for $\phi_{1 \cap 2}$ and $\phi_{2 \cap 3}$ to be _a priori_ independent in $\pd_{2}$.

#### Melded marginal equality

By inspecting Equation \eqref{eqn:coinciding-equalities-one} we note the additional necessary condition that $\pd_{\text{meld}}^{12}(\phi_{2 \cap 3}) = \pd_{2}(\phi_{2 \cap 3})$.
To see when this is true, consider the following derivation
\input{tex-input/noncommutativity/0015-verify-dictatorial-pooling.tex}
Hence, for $\pd_{\text{meld}}^{12}(\phi_{2 \cap 3}) = \pd_{2}(\phi_{2 \cap 3})$ to hold we require $\pd_{\text{pool}}^{12}(\phi_{1 \cap 2}) = \pd_{2}(\phi_{1 \cap 2})$.
This is _dictatorial pooling_, where one submodel's prior marginal is used as a prior in the melded model.
An identical argument can be used to show that $\pd_{\text{meld}}^{23}(\phi_{1 \cap 2}) = \pd_{2}(\phi_{1 \cap 2})$ requires the equivalent choice of dictatorial pooling for $\phi_{2 \cap 3}$, i.e. $\pd_{\text{pool}}^{23}(\phi_{2 \cap 3}) = \pd_{2}(\phi_{2 \cap 3})$.

#### Pooling equality

Without any loss of generality, we examine only the $\pd_{\text{pool}}^{12}(\phi_{1 \cap 2}) \pd_{\text{pool}}^{(12)3}(\phi_{2 \cap 3}) = \pd_{\text{pool}}(\phi_{1 \cap 2}, \phi_{2 \cap 3})$ equality in Equation \eqref{eqn:coinciding-equalities-two} for further necessary conditions.
There are three decisions to be made about the method of pooling used on the right hand side (RHS) of this equality:

1. Form $\pd_{\text{pool}}(\phi_{1 \cap 2}, \phi_{2 \cap 3})$ via the logarithmic pooling method of Equation \eqref{eqn:pooled-prior-overall}.

    Logarithmic pooling results in a RHS proportional to $\pd_{1}(\phi_{1 \cap 2})^{\lambda_{1}} \pd_{2}(\phi_{1 \cap 2}, \phi_{2 \cap 3})^{\lambda_{2}} \pd_{3}(\phi_{2 \cap 3})^{\lambda_{3}}$.
    Hence both $g^{12}$ and $g^{(12)3}$ are logarithmic pooling functions; if either were a linear pooling function we would get more than one term, and if either were dictatorial the left hand side (LHS) would exclude either $\pd_{1}(\phi_{1 \cap 2})$ or $\pd_{3}(\phi_{2 \cap 3})$.
    Thus,
    \input{tex-input/noncommutativity/0041-pooling-equality-rhs-log.tex}
    for arbitrary positives weights $\lambda_{\cdot}$.
    Equation \eqref{eqn:pooling-equality-rhs-log} holds iff $\pd_{2}(\phi_{1 \cap 2}, \phi_{2 \cap 3}) = \pd_{2}(\phi_{1 \cap 2})\pd_{2}(\phi_{2 \cap 3})$ and $\pd_{\text{meld}}^{12}(\phi_{2 \cap 3}) = \pd_{2}(\phi_{2 \cap 3})$.
    We have shown earlier that is only true if dictatorial pooling is used for $g^{12}$.
    But if $g^{12}$ is dictatorial, then the LHS of Equation \eqref{eqn:pooling-equality-rhs-log} would not contain a $\pd_{1}(\phi_{1 \cap 2})$ term, which is a contradiction.
    Hence, Equation \eqref{eqn:pooling-equality-rhs-log} is not true in general.

    A special case of logarithmic pooling is that of product of experts, where $g^{12}(\pd_{1}(\phi_{1 \cap 2}), \pd_{2}(\phi_{1 \cap 2})) = \pd_{1}(\phi_{1 \cap 2})\pd_{2}(\phi_{1 \cap 2})$ and $g^{(12)3}(\pd_{\text{meld}}^{12}(\phi_{2 \cap 3}), \pd_{3}(\phi_{2 \cap 3})) = \pd_{\text{meld}}^{12}(\phi_{2 \cap 3})\pd_{3}(\phi_{2 \cap 3})$.
    In this specific instance all prior terms cancel, and only terms containing data remain in the melded model, so Equation \eqref{eqn:pooling-equality-rhs-log} is trivially true (both sides are equal to 1).

2. From $\pd_{\text{pool}}(\phi_{1 \cap 2}, \phi_{2 \cap 3})$ via the linear pooling method of Equations \eqref{eqn:M-model-linear-pooling-one} -- \eqref{eqn:silly-linear-overall}.

    If linear pooling is used then the RHS contains 4 terms.
    Thus both $g^{12}$ and $g^{1(23)}$ are linear -- all other combinations of pooling functions produce fewer than four terms -- and results in
    \input{tex-input/noncommutativity/0042-pooling-equality-rhs-lin.tex}
    Equation \eqref{eqn:pooling-equality-rhs-lin} again requires $\pd_{\text{meld}}^{12}(\phi_{2 \cap 3}) = \pd_{2}(\phi_{2 \cap 3})$, which we have shown to only be possible under dictatorial pooling.
    Thus we arrive at the same contradiction as before.

3. Use dictatorial pooling.

    There are two valid dictatorial pooling choices for the RHS:

    i. Set $\pd_{\text{pool}}(\phi_{1 \cap 2}, \phi_{2 \cap 3}) = \pd_{2}(\phi_{1 \cap 2}, \phi_{2 \cap 3})$.

        If one chooses $g^{12}(\pd_{1}(\phi_{1 \cap 2}), \pd_{2}(\phi_{1 \cap 2})) = \pd_{2}(\phi_{1 \cap 2})$ and $g^{(12)3}(\pd_{\text{meld}}^{12}(\phi_{2 \cap 3}), \pd_{3}(\phi_{2 \cap 3})) = \pd_{\text{meld}}^{12}(\phi_{2 \cap 3})$, then Equation \eqref{eqn:coinciding-equalities-two} simplifies to $\pd_{2}(\phi_{1 \cap 2}) \pd_{\text{meld}}^{12}(\phi_{2 \cap 3}) = \pd_{2}(\phi_{1 \cap 2}, \phi_{2 \cap 3})$, which is true iff the prior independence assumption is satisfied.
        Additionally, we require $\pd_{\text{meld}}^{12}(\phi_{2 \cap 3}) = \pd_{2}(\phi_{2 \cap 3})$, which as we have noted, is only true if $g^{12}$ is dictatorial.
        In this case, unlike previous cases, $g^{12}$ is exactly the dictatorial pooling function we require.

    i. Set $\pd_{\text{pool}}(\phi_{1 \cap 2}, \phi_{2 \cap 3}) = \pd_{1}(\phi_{1 \cap 2}) \pd_{3}(\phi_{2 \cap 3})$.

        Satisfying Equation \eqref{eqn:coinciding-equalities-two} is possible in this case if $g^{12}(\pd_{1}(\phi_{1 \cap 2}), \pd_{2}(\phi_{1 \cap 2})) = \pd_{1}(\phi_{1 \cap 2})$ and $g^{(12)3}(\pd_{\text{meld}}^{12}(\phi_{2 \cap 3}), \pd_{3}(\phi_{2 \cap 3})) = \pd_{3}(\phi_{2 \cap 3})$.

So in general, applying the original melding operator twice does not result in the same model as \eqref{eqn:melded-model-full}, except in cases where $\phi_{1 \cap 2}$ and $\phi_{2 \cap 3}$ are a priori independent in $\pd_{2}$ and specific forms of dictatorial pooling are used.

### Is the original operator commutative?

Commutativity of the original operator would imply that $(\pd_{1} \circledast \pd_{2}) \circledast \pd_{3} = \pd_{1} \circledast (\pd_{2} \circledast \pd_{3})$.
By carefully considering the indices in Equation \eqref{eqn:expanded-double-melded-model} and \eqref{eqn:symmetric-expanded-double-melded-model}, we find that the original melding operator is only commutative if
\input{tex-input/noncommutativity/0016-commutativity-condition.tex}
which implies the following equalities
\input{tex-input/noncommutativity/0014-orig-melding-commutative-equalities.tex}
Showing one of the equalities in Equation \eqref{eqn:orig-melding-commutative-equalities-1} and \eqref{eqn:orig-melding-commutative-equalities-2} implies its partner equality is also true.
Consider the first equality
\input{tex-input/noncommutativity/0017-pooling-equality.tex}
Assume that $g^{12}$ and $g^{1(23)}$ are both linear or logarithmic pooling functions.
For Equation \eqref{eqn:pooling-equality} to be true, $\pd_{2}(\phi_{1 \cap 2}) = \pd_{\text{meld}}^{23}(\phi_{1 \cap 2})$, which is the same result we require in Equation \eqref{eqn:orig-melding-commutative-equalities-2}.
We have already shown that this is only true when using certain forms of dictatorial pooling.
Hence, the original operator is commutative under in the same settings in which applying the original melding operator twice results in the chained melded model.
