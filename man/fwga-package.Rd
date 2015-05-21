\name{fwga-package}
\alias{fwga-package}
\docType{package}
\title{
Genetic algorithms and null models for food webs or other interaction networks
}
\description{
Genetic algorithm for finding food-web structures that optimize one or more metrics. Includes a collection of null models for binary and quantitative matrices, representing food webs or other interaction networks

}
\details{
\tabular{ll}{
Package: \tab fwga\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-05-21\cr
License: \tab GPL-3\cr
}
Email <jdl.yen@gmail.com>}
\author{
Jian Yen
Mauricio Cantor

Maintainer: Jian Yen <jdl.yen@gmail.com>
}
\references{
  Yen JDL et al. (in review) Highly connected food webs maximize robustness but not ecosystem throughput.
}
\keyword{ package }
\keyword{ food web }
\examples{
  # TRIVIAL EXAMPLE - NOT ENOUGH ITERATIONS
  # run the genetic algorithm with default fitness functions and settings
  sim.fws <- fwga(n.iter=5, suggest=matrix(ifelse(runif(100) > 0.7, 1, 0), ncol=10))
  # convert outputs to binary food webs
  fw.out <- lapply(sim.fws$pops, function(x) ifelse(x > 0, 1, 0))
  # remove cannibalism from the food web
  fw.out <- lapply(fw.out, function(x) {diag(x) <- 0; return(x)})
  # calculate a set of summary metrics for the generated food webs
  out.met <- lapply(fw.out, allMetU)
}
