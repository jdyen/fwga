\name{fit_se}
\alias{fit_se}
\title{Fitness function based on secondary extinctions}
\description{
Fitness function that takes a quantitative food web and returns the average number of secondary extinctions following the random removal of a food-web node.
}
\usage{
fit_se(x, S)
}
\arguments{
  \item{x}{A quantitative interaction matrix (food-web; square matrix)}
  \item{S}{Number of species (columns/rows) in the food web}
}
\value{
Returns the average number of secondary extinctions for a given input food web. Designed for use with \code{\link{fwga}}, which minimizes a chosen fitness function.
}
\references{
  Yen JDL et al. (in review) Highly connected food webs maximize robustness but not ecosystem throughput.
}
\author{
  Jian Yen
  Mauricio Cantor
}

\examples{
  fit_se(matrix(runif(100), ncol=10), S=10)
}