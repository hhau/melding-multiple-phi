> For figure 3 -  it will never tell the full story, but I would explore more whether you can show the gradual shift in pooled density as you change the lambda. You're just looking to give people a feeling of what is going on.

- More rows / values of lambda? Maybe a two, 1 page grids, one for each type of density

> I'm not sure I understand the "Pooling density and dimension" section. I'm not sure what you are trying to do with this (I realise you're also unsure and trying to work it out..). I particularly don't understand the solid red line in Figure 4 - why does it go up? We can discuss.

The idea was:

Say we want each marginal to equally contribute to \pd_{\text{pool}} (\phi_{1 \cap 2}, \phi_{2 \cap 3}) in the _literal_ sense, i.e that for some value(s?) of phi (maybe the pooled mean?), the actual value of the pooled density function is 1/3 of each marginal. The thinking was that higher dimensional density functions have smaller values, which we can see from the blue solid line (the Gaussian density function at 0 gets smaller as the dimension increases). Hence, we could set the weights to be both a function of dimension, so that the actual value of the pooled density was comprised equally of each (different dimensional) marginal. However, the Student-t distribution (As an example of a non-Gaussian density) does not uniformly get smaller as its dimension increases; in a shrinking neighbourhood around 0 it actually increases. This makes it hard to give general advice on specifying the weights as a function of dimension - it depends on where you want to evaluate the pooled prior.

> Did we ever think about pooling conditionals? I seem to recall you worked out it doesn't make sense but can't remember why. Can you remind me?

- I need to check exactly what was done here.

Say there exists a \theta_2 \subset \psi_{2} such that \phi_{1 \cap 2} \indep \phi_{2 \cap 3} \mid \theta_2 (i.e. \theta_2 makes \phi_{1 \cap 2} and \phi_{2 \cap 3} conditionally independent).
We showed that \pd_{\text{meld}}(\cdot \mid \theta_2) was order-independent under the original melding operator. (i.e. you could meld 2 & 3 first, then 1 with 2 & 3, and vice-versa.)  
This was kind of nonsensical as it's a very strange joint/conditional to be interested in (and the corresponding posterior is unlikely to be of interest either). 


Shall we discuss on Monday? 1.30pm say?
