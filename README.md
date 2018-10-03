netprobitFE
===========

R package that implements the methods developed in the paper

Dzemski, A: "An empirical model of dyadic link formation in a network with unobserved heterogeneity", Review of Economics and Statistics, forthcoming.

Installation
------------

The easiest way to install the package is to use the `devtools` package. Make sure that you have the `devtools` package installed and type:

``` r
devtools::install_github("adzemski/netprobitFE")
```

Example
-------

Let's use the package to simulate an example network with *N* = 50 agents.

``` r
library(netprobitFE)
# first define parameters of simulation
sim_design <- define_sim_design_jochmans(N = 50, theta = 1, rho = 0.5, type_C_N = "loglog")
# generate a network with these parameters
set.seed(123)
dyadic_data <- as.data.frame(draw_network_jochmans_2018(sim_design)$links)[, c("i", "j", "X1", "Y")]
head(dyadic_data, n = 10)
```

    ##    i  j X1 Y
    ## 1  1  2 -1 0
    ## 2  1  3  1 0
    ## 3  1  4 -1 0
    ## 4  1  5  1 0
    ## 5  1  6 -1 0
    ## 6  1  7  1 0
    ## 7  1  8 -1 0
    ## 8  1  9  1 0
    ## 9  1 10 -1 0
    ## 10 1 11  1 0

Next, we tell the package about what data we are going to use and what model we are going to estimate.

``` r
probit_model <- define_model(dyadic_data, X_names = "X1", col_sender = "i", col_receiver = "j",
                                          col_outcome = "Y")
```

We now estimate the homophily parameter *θ* and the reciprocity parameter *ρ* using the two-stage procedure described in Section 2.2 of the paper.

``` r
# stage 1
fit_ML <- MLE_stage1(probit_model)
fit_ML$theta_hat_MLE
```

    ##       X1 
    ## 1.041016

``` r
# stage 2
fit_ML <- MLE_stage2(fit_ML)
fit_ML$rho_hat_MLE
```

    ##       rho 
    ## 0.5474438

The procedure `ttest_theta` computes a bias-corrected estimate of *θ* and conducts *t*-tests that are robust to the incidental parameter problem. For example we can test significance, i.e., test *θ* against zero:

``` r
fit_ttest_theta <- ttest_theta(fit_ML, theta0 = 0)
fit_ttest_theta[c("theta_hat_bias_corr", "tstat_theta", "pval_theta")]
```

    ## $theta_hat_bias_corr
    ##        X1 
    ## 0.9779992 
    ## 
    ## $tstat_theta
    ##       X1 
    ## 17.15144 
    ## 
    ## $pval_theta
    ##           X1 
    ## 6.131079e-66

The test rejects at the 5% level. This means that `X1` is significant at the 5% level. Note that the bias-corrected estimate is indeed closer to the true value (*θ* = 1) than the uncorrected ML estimate. We can also test against the true value of *θ*:

``` r
ttest_theta(fit_ML, theta0 = 1)[c("theta_hat_bias_corr", "tstat_theta", "pval_theta")]
```

    ## $theta_hat_bias_corr
    ##        X1 
    ## 0.9779992 
    ## 
    ## $tstat_theta
    ##         X1 
    ## -0.3858337 
    ## 
    ## $pval_theta
    ##        X1 
    ## 0.6996198

Unsurprisingly, this test does not reject.

The procedure `ttest_rho` computes a bias-corrected estimate of *ρ* and conducts *t*-tests that are robust to the incidental parameter problem. For example, we can test significance:

``` r
ttest_rho(fit_ML, fit_ttest_theta, rho0 = 0)[c("rho_hat_bias_corr", "pval_rho")]
```

    ## $rho_hat_bias_corr
    ##       rho 
    ## 0.5179216 
    ## 
    ## $pval_rho
    ##          rho 
    ## 5.230526e-17

The test rejects at the 5% level. This means that the within dyad shocks are correlated (95% confidence). Note that the bias-corrected estimate is indeed closer to the true value (*ρ* = 0.5) than the uncorrected ML estimate. We can also test against the true value of *θ*:

``` r
ttest_rho(fit_ML, fit_ttest_theta, rho0 = 0.5)[c("rho_hat_bias_corr", "pval_rho")]
```

    ## $rho_hat_bias_corr
    ##       rho 
    ## 0.5179216 
    ## 
    ## $pval_rho
    ##       rho 
    ## 0.7718001

This test rejects at the 5% level.

Finally, we test model specification using the transitivity test described in Section 5 of the paper. We set the nominal size of the test to *α* = 0.5 and turn bootstrap inference on (as recommended in the paper).

``` r
triangle_test(fit_ML, fit_ttest_theta, bs_inference = "double", 
              bs_iterations = 200, alpha = 0.05)[c("excess_trans", "excess_trans_bias_corr", "teststat",
                                                   "pval_analytic", "rej_bs_perc_t")]
```

    ## $excess_trans
    ## [1] -0.001561924
    ## 
    ## $excess_trans_bias_corr
    ## [1] -7.116839e-05
    ## 
    ## $teststat
    ## [1] -0.3406273
    ## 
    ## $pval_analytic
    ## [1] 0.7333842
    ## 
    ## $rej_bs_perc_t
    ## [1] FALSE

The statistic `excess_trans_bias` is the naive excess transitivity measure *T*<sub>*N*</sub> from Section 5.2, the statistic `excess_trans_bias_corr` is the same statistic with analytic bias correction based on Theorem 3. Finally, `teststat` is the studentized robust test statistic given in equation (5.3). The percentile bootstrap does not reject (`rej_bs_perc_t = FALSE`). The analytical *p*-value is also greater than 5%. For this example, bootstrap and normal approximation both do not reject.

License
-------

This software is open source software and is licensed under the [GNU General Public License 3.0](https://www.gnu.org/licenses/gpl-3.0.en.html).
