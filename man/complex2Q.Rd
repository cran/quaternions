\name{complex2Q}
\alias{complex2Q}
\title{Embeds complex numbers in real quaternions}
\description{
  The real and imaginary part of the complex vector/tensor is used as
  real and first imaginary part of a quaternions.
}
\usage{
complex2Q(x)
}
\arguments{
  \item{x}{a complex vector or tensor}
  }
\value{
  an quaternion object of class \code{"Q"} containing the given values.
  
}
\details{
  Complex numbers are contain in the quaternions in 4 different
  ways:
  \code{Q(Re(x),Im(x))}, \code{Q(Re(x),j=Im(x))},
  \code{Q(Re(x),k=Im(x))},
  \code{Q(r=x)}. The first of these is modeled by complex2Q.
}
\author{K.Gerald van den Boogaart}
\examples{
complex2Q(1+c(1i,2i))
}
\keyword{arith}