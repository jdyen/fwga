\name{fit_rand}
\alias{fit_rand}
\title{Fitness function returning a constant value}
\description{
Fitness function that takes a quantitative or binary food web and returns a constant value.
}
\usage{
fit_rand(x, S)
}
\arguments{
  \item{x}{A quantitative interaction matrix (food-web; square matrix)}
  \item{S}{Number of species (columns/rows) in the food web}
}
\value{
Returns a constant value regardless of the input matrix. This function is designed as a null model for use with \code{\link{fwga}}, and will lead to random food webs being generated.
}
\references{
  Yen JDL et al. (in review) Highly connected food webs maximize robustness but not ecosystem throughput.
}
\author{
  Jian Yen and Mauricio Cantor
}
\seealso{
  \code{\link{fit_se}}, \code{\link{fit_robust}}, \code{\link{fit_wtst}},
  \code{\link{fit_tst}}
}
\examples{
  # calculate constant-value fitness for a randomly
  #  generated food web
  fit_rand(matrix(runif(100), ncol=10), S=10)
}