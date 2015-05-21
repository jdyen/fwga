\name{null2}
\alias{null2}
\title{Null model 2: Randomize binary matrices based on rows and columns totals}
\usage{
  null2(mat, iter = 1000)
}
\arguments{
  \item{mat}{Binary matrix}

  \item{iter}{Number of iterations}
}
\value{
  mat.t array of randomized matrices
}
\description{
  This null model created random binary matrices by
  resorting the 1's among the matrix cells according to
  marginal totals of rows and columns.
}
\details{
  Each cell has a probability of being filled that is
  proportional to the number of occurrences of individuals
  in sites: cij = 1/2*(Pi/C + Pj/R) where Pi= row sums; Pj
  = column sums; C = number of columns; and R = number of
  rows.
}
\author{
  Mathias Pires with modifications by Diego Barneche and
  Mauricio Cantor
}
\references{
  "Bascompte J, Jordano P, Melian CJ and Olesen JM (2003)
  The nested assembly of plant-animal mutualistic networks.
  PNAS, 100(16):9383-9387"
}
\seealso{
  \code{\link{null1}}, \code{\link{null3}}, \code{\link{null4_7}},
  \code{\link{mgen}}
}