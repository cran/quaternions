\name{matrix}
\alias{matrix.default}
\alias{matrix.Q}
\alias{matrix}
\title{A matrix of quaternions}
\description{
  Constructs a matrix of quaternions.
}
\usage{
 matrix(data,...)
 \method{matrix}{default} (data = NA, nrow = 1, ncol = 1, byrow = FALSE,dimnames = NULL,...)
 \method{matrix}{Q}(data,...)
}
\arguments{
  \item{data}{A vector giving the entries of the matrix to be
    constructed}
  \item{nrow}{integer giving the number of rows}
  \item{ncol}{integer giving the number of columns}
  \item{byrow}{boolean, wether to fill the matrix by rows}
  \item{dimnames}{The dimnames of the matrix}
  \item{\dots}{further arguments or again \code{nrow}, \code{ncol},
    \code{byrow} or \code{dimnames}}
}
\value{
  A quaternion with dimension for a matrix in each component.
}
\details{
 The matrix dimension is stored in every component individually.  
}
\seealso{
  \code{\link{Q}}
}
\examples{
A <- matrix(Q(1:4,2:5),nrow=2)
B <- matrix(c(Q(i=c(1,0,0,1))),nrow=2)
A \%*\% B
}
\author{K.Gerald van den Boogaart}
\keyword{arith}