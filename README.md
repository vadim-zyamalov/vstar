
<!-- README.md is generated from README.Rmd. Please edit that file -->

# vstar

<!-- badges: start -->
<!-- badges: end -->

The goal of vstar is to test and estimate Vector Smooth Transition
Autoregression models (VSTAR). The code is based on the specification
from (Teräsvirta and Yang 2014b, 2014a).

(Teräsvirta and Yang 2014b) describes the most general model, allowing
for different transition variables and threshold values for every
equation of the system and every regime. However, in the package
implemented is much simpler, yet suitable specification, which uses a
common transition variable for all regimes and a common threshold value
in each regime.

## Installation

You can install the development version of vstar from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("vadim-zyamalov/vstar")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(vstar)
```

Let’s generate an artifical data with two regimes using the following
matrices: $$A_1' = \left(\matrix{.5 & 0 \cr 0 & .5}\right)$$
$$A_2' = \left(\matrix{.8 & .1 \cr .1 & .8}\right)$$ The transition
variable is a normal random variable with zero mean and unit variance.
The threshold value is a median value of the transition variable. The
smoothness parameter $\gamma$ equals 4. We’ll use a logistic transition
function.

<img src="man/figures/README-data-generation-1.png" width="100%" />

First, we prepare data for the LSTAR-model estimation:

``` r
model.2reg <- vstar.prepare(endo = c("y1", "y2"),
                            p = 1,
                            trans = s,
                            dataset = y.2reg)
```

The linearity test statistics reject the null of linearity:

``` r
round(linearity.test(model.2reg, J = 1), 4)
#>                     stat crit.10 crit.5  crit.1 p.value
#> LM              764.1390  7.7794 9.4877 13.2767       0
#> r-LM            190.2295  1.9478 2.3766  3.3291       0
#> Rao             366.1712  1.9478 2.3766  3.3291       0
#> Bartlett-Wilks 1084.7458  7.7794 9.4877 13.2767       0
```

Next, we obtain the starting values for Nonlinear Least Squares by the
grid search.

``` r
result.2reg <- vstar.grid(dataset = model.2reg,
                          m = 2,
                          gamma.limits = c(3, 5),
                          points = 200,
                          cores = 10)
```

Then we estimate the model with NLS.

``` r
result.nls.2reg <- vstar.nls(result.2reg, tol = 1e-6, verbose = TRUE)
#> Step 1: delta = 0.002780002
#> Step 2: delta = 0.0001072831
#> Step 3: delta = 0
#> Tolerance level 1e-06 achieved in 3 steps.
```

Then we print out the summary of the estimated model.

``` r
summary(result.nls.2reg)
#> Results of L-VSTAR model estimation
#>
#> Number of regimes:      2
#> Number of lags:         1
#> Number of observations: 949
#> Endogenous variables:
#>  y1 y2
#> Transition variable:    custom
#>
#> Additional model components:
#> Constant included:               FALSE
#> Trend included:                  FALSE
#>
#> Transition function: logistic
#> gamma values:
#> [1] 3.341556
#> threshold values:
#> [1] -0.07474272
#>
#> Matrix of coefficients:
#>            R1.y1      R1.y2      R2.y1     R2.y2
#> l1.y1 0.43108671 -0.0373987 0.84092522 0.1328473
#> l1.y2 0.03924485  0.4748471 0.05479257 0.8376930
#>
#> Matrix of coefficients SD:
#>            R1.y1      R1.y2      R2.y1     R2.y2
#> l1.y1 0.02323907 0.02323704 0.03767280 0.0376695
#> l1.y2 0.01924053 0.01923884 0.03049446 0.0304918
#>
#> Residuals covariance matrix:
#>             resid.y1    resid.y2
#> resid.y1 1.041408739 0.006190959
#> resid.y2 0.006190959 1.041226692
```

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-terasvirta2014linearity" class="csl-entry">

Teräsvirta, Timo, and Yukai Yang. 2014a. “Linearity and Misspecification
Tests for Vector Smooth Transition Regression Models.” 2014-04.
Department of Economics and Business Economics, Aarhus University.
<https://ideas.repec.org/p/aah/create/2014-04.html>.

</div>

<div id="ref-terasvirta2014specification" class="csl-entry">

———. 2014b. “Specification, Estimation and Evaluation of Vector Smooth
Transition Autoregressive Models with Applications.” 2014-08. Department
of Economics and Business Economics, Aarhus University.
<https://ideas.repec.org/p/aah/create/2014-08.html>.

</div>

</div>
