\name{fit_robust}
\alias{fit_robust}
\title{Fitness function based on robustness}
\description{
Fitness function that takes a quantitative or binary food web and returns the negative of the average number of link removals required to cause loss of half of all species.
}
\usage{
fit_robust(x, S)
}
\arguments{
  \item{x}{A quantitative interaction matrix (food-web; square matrix)}
  \item{S}{Number of species (columns/rows) in the food web}
}
\value{
Returns the negative of the average number of link removals required to cause loss of half of all species. A negative value is returned because this function is designed for use with \code{\link{fwga}}, which minimizes a chosen fitness function.
}
\references{
  Yen JDL et al. (in review) Highly connected food webs maximize robustness but not ecosystem throughput.
}
\author{
  Jian Yen and Mauricio Cantor
}
\seealso{
  \code{\link{fit_se}}, \code{\link{fit_wtst}},
  \code{\link{fit_tst}}, \code{\link{fit_rand}}
}
\examples{
  # calculate (negative) robustness for a randomly
  #  generated food web
  fit_robust(matrix(runif(100), ncol=10), S=10)
}