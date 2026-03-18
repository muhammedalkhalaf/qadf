\name{qadf}
\alias{qadf}
\alias{print.qadf}
\alias{summary.qadf}
\alias{plot.qadf}
\title{Quantile Augmented Dickey-Fuller Unit Root Test}
\description{
Performs the Quantile ADF unit root test following Koenker & Xiao (2004).
Tests for unit root across the conditional distribution using quantile regression.
}
\usage{
qadf(y, model = c("c", "ct", "nc"), pmax = 4, 
     ic = c("AIC", "BIC", "t-stat"), 
     tau = seq(0.1, 0.9, 0.1))

\method{print}{qadf}(x, ...)
\method{summary}{qadf}(object, ...)
\method{plot}{qadf}(x, ...)
}
\arguments{
  \item{y}{A numeric vector or time series.}
  \item{model}{Character. Deterministic component: "c" for constant (default),
    "ct" for constant and trend, "nc" for no constant.}
  \item{pmax}{Integer. Maximum number of lags for augmentation. Default is 4.}
  \item{ic}{Character. Information criterion for lag selection: "AIC", "BIC", or "t-stat".}
  \item{tau}{Numeric vector. Quantiles to test. Default is seq(0.1, 0.9, 0.1).}
  \item{x, object}{An object of class "qadf".}
  \item{...}{Additional arguments (currently ignored).}
}
\value{
An object of class "qadf" containing:
\item{results}{Data frame with test results for each quantile}
\item{qks}{QKS statistic (supremum of absolute t-statistics)}
\item{qks_cv}{Critical values for QKS statistic}
\item{model}{Model specification}
\item{lags}{Selected lag order}
\item{n}{Sample size}
\item{y}{Original series}
}
\details{
The Quantile ADF test extends the standard ADF test to different quantiles
of the conditional distribution. This allows testing for unit root behavior
that may vary across the distribution (e.g., different persistence in 
expansions vs recessions).

The test statistic at quantile tau is:
\deqn{t_n(\tau) = f_z \sqrt{y'_{-1} P_X y_{-1}} (\hat{\rho}(\tau) - 1)}

Critical values are interpolated from Hansen (1995) response surfaces.

The QKS (Quantile Kolmogorov-Smirnov) statistic is:
\deqn{QKS = \sup_{\tau} |t_n(\tau)|}
}
\references{
Koenker, R., & Xiao, Z. (2004). Unit root quantile autoregression inference.
\emph{Journal of the American Statistical Association}, 99(467), 775-787.
\doi{10.1198/016214504000000296}

Hansen, B. E. (1995). Rethinking the univariate approach to unit root testing.
\emph{Econometric Theory}, 11(5), 1148-1171.
}
\examples{
# Generate a random walk
set.seed(123)
y <- cumsum(rnorm(200))

# Basic QADF test
result <- qadf(y)
print(result)

# Test specific quantiles
result <- qadf(y, tau = c(0.1, 0.5, 0.9))

# With trend
result <- qadf(y, model = "ct")
summary(result)
}
\keyword{htest}
\keyword{ts}
