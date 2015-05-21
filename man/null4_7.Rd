\name{null4_7}
\alias{null4_7}
\title{Null models 4 to 7:}
\usage{
  null4_7(mat, iter, model, ...)
}
\arguments{
  \item{mat}{A quantitative matrix}

  \item{iter}{Number of random matrices to be created}

  \item{model}{Function to be chosen. See \code{Details}.}

  \item{...}{Further arguments from \code{permatswap} or
  \code{permatfull}}
}
\value{
  a list with \code{iter} random matrices
}
\description{
  shuffles XXX among XXX. It calls algorithms from package
  \code{vegan} to generate random matrices.
}
\details{
  If \code{model} = 4, an unrestricted null model is
  called. All cell values are randomly permuted. Only total
  total sum is constrained; row sums, column sums, and
  matrix fill are allow to change.

  If \code{model} = 5, cell values are permuted restricting
  only column sums and total sum.

  If \code{model} = 6, cell values are permuted restricting
  only row sums and total sum.

  If \code{model} = 7, totally restricted null model is
  called. Cell values are permuted restricting all features
  of the original matrix: column sums, row sums, matrix
  fill and total sum.
}
\author{
  Mauricio Cantor and Diego Barneche, using functions
  available in \code{vegan}
}
\references{
  \code{citation("vegan")}
}
\seealso{
  \code{library(vegan)}
}

