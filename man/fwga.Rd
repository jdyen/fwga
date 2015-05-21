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
  This is the main genetic-algorithm function and is designed to minimize the chosen fitness functions. If maximization is required, the fitness functions must be multiplied by minus one; it is easiest to do this in the function definition (see \code{\link{fit_robust}}).
  
  Default settings should work for many situations but
  damp.fact and n.iter are the two most useful parameters to change.
  
  Convergence is not guaranteed using this genetic algorithm. A simple method for determining convergence is to plot fitness against iterations. Obvious tapering-off of the fitness functions suggest, but do not guarantee, convergence. See below for an example of this plot.
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
  # run the genetic algorithm with the default fitness functions
  sim.fws <- fwga(n.iter=100, suggest=matrix(ifelse(runif(100) > 0.7, 1, 0), ncol=10))
  # plot fitness against iterations for the first (of two) fitness functions to assess convergence
  plot(sim.fws$fitness[[1]][, 1], type='l', bty='l', xlab="Iteration", ylab="Fitness")
}