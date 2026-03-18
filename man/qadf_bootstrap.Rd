\name{qadf_bootstrap}
\alias{qadf_bootstrap}
\title{Bootstrap QADF Test}
\description{
Computes bootstrap p-value for the QKS statistic from a QADF test.
}
\usage{
qadf_bootstrap(object, nboot = 199, seed = NULL)
}
\arguments{
  \item{object}{A "qadf" object from \code{\link{qadf}}.}
  \item{nboot}{Number of bootstrap replications. Default 199.}
  \item{seed}{Random seed for reproducibility.}
}
\value{
Updated "qadf" object with additional elements:
\item{qks_pvalue}{Bootstrap p-value for the QKS statistic}
\item{qks_boot}{Vector of bootstrap QKS statistics}
\item{nboot}{Number of successful bootstrap replications}
}
\details{
The bootstrap procedure follows Koenker and Xiao (2004):
\enumerate{
  \item Fit AR(p) model to first differences
  \item Resample centered residuals
  \item Generate bootstrap series under unit root null
  \item Calculate QKS statistic for each bootstrap sample
  \item P-value is proportion of bootstrap QKS >= observed QKS
}
}
\references{
Koenker, R., & Xiao, Z. (2004). Unit root quantile autoregression inference.
\emph{Journal of the American Statistical Association}, 99(467), 775-787.
}
\examples{
\donttest{
set.seed(123)
y <- cumsum(rnorm(200))
result <- qadf(y)
result <- qadf_bootstrap(result, nboot = 99)
print(result)
}
}
\seealso{
\code{\link{qadf}}
}
\keyword{htest}
