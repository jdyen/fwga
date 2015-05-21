\name{extract.metrics}
\alias{extract.metrics}
\title{Extract a set of food-web metrics from the outputs of \code{\link{allMetU}}}
\description{
Helper function to convert list of lists to list of numeric vectors
}
\usage{
extract.metrics(obj)
}
\arguments{
  \item{obj}{ A list of food-web summaries outputted from an
  	          \code{lapply} call to \code{\link{allMetU}}}
}
\details{
  Purely for convenience: converts a list of lists to a list of numeric vectors.
}
\value{
  List of numeric vectors; one element for each food-web metric.
  \item{connect}{Connectance values for each input food web}
  \item{basal}{Number of basal species for each input food web}
  \item{int}{Number of intermediate species for each input food web}
  \item{top}{Number of top species for each input food web}
  \item{gen}{Generality values for each input food web}
  \item{vul}{Vulnerability values for each input food web}
  \item{sdGen}{Standard deviation of generality values for each input food web}
  \item{sdVul}{Standard deviation of vulnerability values for each input food web}
  \item{mean.pl}{Mean path length values for each input food web}
  \item{mean.tl}{Mean trophic level values for each input food web}
  \item{mean.oi}{Mean omnivory index values for each input food web}
  \item{fg.mod}{Modularity values calculated using a fast-greedy algorithm for each input food web}
  \item{sg.mod}{Modularity values calculated using a spinglass algorithm for each input food web}
  \item{le.mod}{Modularity values calculated using a leading-eigenvector algorithm for each input food web}
  \item{rw.mod}{Modularity values calculated using a random-walk algorithm for each input food web}
  \item{clust}{Clutering coefficient values for each input food web}
  \item{trans}{Transitivity values for each input food web}
}
\references{
  Yen JDL et al. (in review) Highly connected food webs maximize robustness but not ecosystem throughput.
}
\author{
  Jian Yen and Mauricio Cantor
}

\examples{
  # simulate food webs using the genetic algorithm
  sim.fws <- fwga(n.iter=5, suggest=matrix(ifelse(runif(100) > 0.7, 1, 0),
                  ncol=10))
  # convert outputs to binary food webs
  fw.out <- lapply(sim.fws$pops, function(x) ifelse(x > 0, 1, 0))
  # remove cannibalism from the food web
  fw.out <- lapply(fw.out, function(x) {diag(x) <- 0; return(x)})
  # calculate a set of summary metrics for the generated food webs
  out.met <- lapply(fw.out, allMetU)
  # convert the summary metrics from a list of lists to a list
  #  of numeric vectors
  out.met <- extract.metrics(out.met)
}