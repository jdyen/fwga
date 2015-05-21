\name{allMetU}
\alias{allMetU}
\title{
Calculate a set of food-web metrics from a given food web
}
\description{
Calculate a set of food-web metrics from a given food web.
}
\usage{
allMetU(A)
}
\arguments{
  \item{A}{ A food web matrix}
}
\details{
  The input food web can be quantitative or binary but some metrics cannot be calculated if loops are present.
}
\value{
List of metric values.
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
  \item{sg.mod}{*currently not implemented* Modularity values calculated using a spinglass algorithm for each input food web}
  \item{le.mod}{Modularity values calculated using a leading-eigenvector algorithm for each input food web}
  \item{rw.mod}{Modularity values calculated using a random-walk algorithm for each input food web}
  \item{clust}{Clutering coefficient values for each input food web}
  \item{trans}{Transitivity values for each input food web}
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
  fw.out <- lapply(sim.fws$pops, function(x) ifelse(x > 0, 1, 0))
  fw.out <- lapply(fw.out, function(x) {diag(x) <- 0; return(x)})
  out.met <- lapply(fw.out, allMetU)
}