

# Suggestions to Likelihood and Bayesian Inference, 2nd edition (2020)

## Clarification in Example 7.4 (Hardy–Weinberg equilibrium), p. 230

The distinction between the model with one parameter and the one with zero
parameters would be clearer if the text on page 230 read:

> ..., substantially larger than for the Hardy–Weinberg equilibrium **with unknown $\upsilon$** and for the saturated trinomial model. ...

## Addition to Section 7.2.1 (Marginal likelihood and the Bayes factor)

Add the Savage–Dickey density ratio representation of the Bayes factor:

$$
\mathrm{BF}_{12}
= \frac{p(x \mid M_1)}{p(x \mid M_2)}
= \frac{p(\theta_0 \mid x, M_2)}{p(\theta_0 \mid M_2)}.
$$

This could, for example, be included before Example 7.6 (p.233).

## Addition to Appendix A.3.5, p. 352

For two $p$-dimensional random variables $\boldsymbol{X}$ and $\boldsymbol{Y}$,
we have

$$
\mathrm{Cov}(\boldsymbol{Y})
= \mathrm{E}\{\mathrm{Cov}(\boldsymbol{Y} \mid \boldsymbol{X})\} + \mathrm{Cov}\{\mathrm{E}(\boldsymbol{Y} \mid \boldsymbol{X})\}.
$$

This is the law of total covariance, analogous to the law of total variance in (A.20).
