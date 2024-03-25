logisticGP package vignette
================

## Introduction

This is an introduction to `logisticGP` package.

Given a response-covariate pair $\left(y, x^T\right) \in R \times R^p$,
Conditional density estimation (CDE) aims to find the density function
of $y$ conditioned on the $x^T$ slice, denoted by
$f\left(y \mid x^T\right)$. In Bayesian paradigm, we assume prior
knowledge on $f$ and subsequently find posterior estimates of $f$ given
the data.

In [\[our article\]](), we focus on developing a computationally
tractable technique to model the conditional density using Logistic
Gaussian Process (LGP) prior. The main idea behind this formulation is
to use the triangular basis to approximate the Gaussian process (GP)
using a pre-fixed regular grid, as discussed in [Maatouk and Bay
(2017)](https://doi.org/10.1007/s11004-017-9673-2). The logistic density
transform was introduced by [Leonard,
1978](https://doi.org/10.1111/j.2517-6161.1978.tb01655.x) and we utilize
the foundation to model the conditional density of spatially varying
response in the presence of a high dimensional covariate space.

## Vignettes

The vignettes of `logisticGP` package can be found
[here](https://htmlpreview.github.io/?https://github.com/Zmup2eSNntdt/logisticGP/blob/main/logisticGP.html).
