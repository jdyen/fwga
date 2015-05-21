\name{null1}
\alias{null1}
\title{Null model 1: permuting link weights among nodes}
\usage{
  null1(edgelist, iter = 1000)
}
\arguments{
  \item{edgelist}{3-column matrix object}

  \item{iter}{Number of random edgelists to be created}
}
\value{
  list of randomized edgelists
}
\description{
  Permutation of link weights among nodes in a edgelist
}
\details{
  Shuffles the link weights among nodes in a network
}
\author{
  Mauricio Cantor
}
\references{
  "Opsahl et al 2008 Physical review letters 101, 168702";
  "Bascompte et al. 2005. PNAS"
}
\seealso{
  \code{\link{null2}}, \code{\link{null3}}, \code{\link{null4_7}},
  \code{\link{mgen}}
}