\name{dim.Q}
\alias{dim.Q}
\alias{dim<-.Q}
\alias{dimnames.Q}
\alias{dimnames<-.Q}
\title{Output of quaternions}
\description{
  Sets the dimension of a quaternion.
}
\usage{
\method{dim}{Q}(x)
\method{dim}{Q}(x) <- value
\method{dimnames}{Q}(x)
\method{dimnames}{Q}(x) <- value
}
\arguments{
  \item{x}{a quaternion vector/tensor}
  \item{value}{a valid dimension}
  }
\value{
  the dimension or the x with dimension set appropriately.
}
\details{
Quaternions are stored as a 4 list of containing vectors giving each one
of the for entries. In case a dimension attribute is applied, it is
applied to all four elements of the list.  
}
\seealso{
  \code{\link{Q}}
}
\examples{
x <- Q(1:10,1:10)
dim(x) <- c(2,5)
dim(x)
x
}
\keyword{arith}