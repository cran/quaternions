\name{Im}
\alias{Im<-}
\alias{Re<-}
\alias{Im<-.complex}
\alias{Re<-.complex}
\title{Real and Imaginary part of complex numbers}
\description{
  Assigning the real and imaginary parts seperately
}
\usage{
Re(z,...) <- value
\method{Re}{complex}(z,...) <- value
Im(z,...) <- value
\method{Im}{complex}(z,...) <- value
}
\arguments{
  \item{z}{a complex number/vector}
  \item{value}{a real number/vector}
  \item{\dots}{further arguments for generic use only}
}
\value{
  This are assignment functions modifying the real or imaginary part of x.
}
\details{
  This functions allow the assigment of the real and imaginariy part of
  a complex number seperatly.
}
\seealso{
  \code{\link{base}{Im}},   \code{\link{base}{Re}}
}
\examples{
x <- 1+3i
x
Im(x)<- 2
x
Re(x)<- -1
x
}
\author{K.Gerald van den Boogaart}
\keyword{arith}