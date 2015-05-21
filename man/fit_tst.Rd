\name{fit_tst}
\alias{fit_tst}
\title{Fitness function based on total-system throughput}
\description{
Fitness function that takes a quantitative food web and returns the negative of the total system throughput (the sum of all link weights).
}
\usage{
fit_tst(x, S)
}
\arguments{
  \item{x}{A quantitative interaction matrix (food-web; square matrix)}
  \item{S}{Number of species (columns/rows) in the food web}
}
\value{
Returns the negative of the total system throughput for a given input food web. The negative output is designed for use with \code{\link{fwga}}, which minimizes a chosen fitness function.
}
\references{
  Yen JDL et al. (in review) Highly connected food webs maximize robustness but not ecosystem throughput.
}
\author{
  Jian Yen and Mauricio Cantor
}
\seealso{
  \code{\link{fit_se}}, \code{\link{fit_robust}}, \code{\link{fit_wtst}},
  \code{\link{fit_rand}}
}
\examples{
  # calculate (negative) total throughput for a randomly
  #  generated food web
  fit_tst(matrix(runif(100), ncol=10), S=10)
}