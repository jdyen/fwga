\name{mat_el}
\alias{mat_el}
\title{Transforms matrix into edgelist}
\usage{
  mat_el(mat, type = "weighted")
}
\arguments{
  \item{mat}{matrix of interactions}

  \item{type}{\code{binary} returns a binary edgelist and
  \code{weighted} returns a weighted edgelist}
}
\description{
  Transforms matrix into edgelist
}
\author{
  Mauricio Cantor
}
\seealso{
  \code{\link{el_mat}}
}