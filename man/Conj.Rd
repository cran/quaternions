\name{Conj}
\alias{Conj.Q}
\title{Conjugation of quaternions}
\description{
  Computes the conjugated real or complex quaternion.
}
\usage{
\method{Conj}{Q}(z)
}
\arguments{
  \item{z}{a quaternion vector to be conjugated}
}
\value{
  A quaternion object representing the conjugated quaternions.
}
\details{
  The conjugated quaternions has have inverse entries in the imaginary
  dimensions. 
}
\seealso{
  \code{\link{adj}}
}
\examples{
q <- Q(1:2,3:4,5:6,7:8)
q
Conj(q)

}
\author{K.Gerald van den Boogaart}
\keyword{arith}