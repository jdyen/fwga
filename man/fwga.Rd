\name{fwga}
\alias{fwga}
\title{Genetic algorithm for optimal food-web structures}
\description{
Simulates optimal food-web structures according to one or more fitness functions
}
\usage{
fwga(n.pop = 200, n.iter = 100,
     fits = list(fit_tst, fit_robust),
     suggest = NULL, k = 5, kp = 0.5,
     n.keep = round(n.pop/5), damp.fact=1)
}
\arguments{
  \item{n.pop}{ The number of populations (food webs) to be used in the genetic algorithm.}
  \item{n.iter}{ The number of iterations to run the genetic algorithm}
  \item{fits}{A list containing one or more fitness functions, which must return a value
  	          to be minimized.}
  \item{suggest}{*Required* A suggested starting population.}
  \item{k}{ Number of populations to consider in the tournament-selection step (default is k=5 and should be adequate for most cases).}
  \item{kp}{ Probability of selecting a sub-optimal population in the tournament-selection process.}
  \item{n.keep}{ Number of populations to retain at each iteration.}
  \item{damp.fact}{ Amount by which to dampen the mutation rates in each iteration}
}
\details{
  This is the main genetic-algorithm function. Default settings should work for many situations but
  damp.fact and n.iter are the two most useful parameters to change.
}
\value{
List of outputs from the genetic algorithm
  \item{pops}{List containing the final n.keep populations; one element for each population.}
  \item{fitness}{List containing matrices of fitness values for each fitness function; there
  	is one matrix for each fitness function, with iterations in rows and populations in columns.}
  \item{doms}{Dominance values from the final step of the genetic algorithm; included for
  	diagnostic purposes only.}
}
\references{
  Yen JDL et al. (in review) Highly connected food webs maximize robustness but not ecosystem throughput.
}
\author{
  Jian Yen
  Mauricio Cantor
}

\examples{
  sim.fws <- fwga(n.iter=5, suggest=matrix(ifelse(runif(100) > 0.7, 1, 0), ncol=10))
}