\name{null3}
\alias{null3}
\title{Null model 3: Version of \code{null 2} for quantitative matrices and weighted networks}
\usage{
  null3(mat, iter, ...)
}
\arguments{
  \item{mat}{A binary or quantitative matrix}

  \item{iter}{Number of iteractions, i.e. number of random
  matrices to be created}

  \item{...}{Further arguments from \code{mgen}}
}
\value{
  a list with iter random matrices
}
\description{
  This null model followes the same idea of null2 by
  Bascompte et al. 2003, but can be applied bot to binary
  and quantitative matrices. It uses a slightly different
  way of producing a probabilistic matrix for reference,
  calling the algorithm \code{mgen} (see help(mgen)).
}
\author{
  Mauricio Cantor with modifications by Diego Barneche
}
\seealso{
  \code{\link{null1}}, \code{\link{null2}}, \code{\link{null4_7}},
  \code{\link{mgen}}
}