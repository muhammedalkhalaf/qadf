## qadf 1.0.2

This release fixes the test statistic itself; the 1.0.1 submission (DOI correction only) should be discarded in favour of this one.

* Bug fix: the test regression used the first difference as the dependent variable, so the coefficient on y(t-1) was rho - 1, not rho; the statistic then subtracted 1 a second time. The quantile autoregression is now estimated in levels, as in Koenker and Xiao (2004), and `rho_tau`, `rho_ols`, `coef_stat` and `half_life` refer to rho.
* Bug fix: the statistic now follows equation (9) of Koenker and Xiao (2004): the density at the quantile is the difference quotient of the fitted conditional quantile at tau +/- h (Hall-Sheather bandwidth), and the regressor y(t-1) is projected off the constant, the lagged differences and, for `model = "ct"`, the trend. The previous version used a kernel estimate on residuals and the OLS moment matrix.
* Bug fix: the critical values now depend on the estimated nuisance parameter delta^2, as in Hansen (1995), interpolated on the grid 0.1, ..., 1. The previous table was indexed by tau, which the limiting distribution does not depend on.
* Bug fix: `delta2` was always sigma^2 because the sum of the lag coefficients matched no column name; it is now the squared correlation between the differenced series and psi_tau of the quantile residuals.
* Lag selection (AIC, BIC and sequential t) now uses the ADF regression on a common sample, with a trend for `model = "ct"`.
* Results agree with the Stata command qadf (SSC) for `model = "c"` on the same simulated series and lag order (t = -1.166 in both).

## Test environments

* Ubuntu 24.04, R 4.3.3 and R-devel, R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes
