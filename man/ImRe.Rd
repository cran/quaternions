\name{Imaginary}
\alias{Im.Q}
\alias{Re.Q}
\alias{Im<-.Q}
\alias{Re<-.Q}
\title{Real and Imaginary part of the quaternions}
\description{
  The first element of quaterion is its real part and others its
  imaginary part. 
}
\usage{
\method{Im}{Q}(z)
\method{Re}{Q}(z)
\method{Im}{Q}(z,...)<- value
\method{Re}{Q}(z,...)<- value
}
\arguments{
  \item{z}{The Quaternion}
  \item{value}{a numeric for Re and a quaternion vor Im}
  \item{\dots}{Further generic arguments}
}
\value{
  a numeric or complex in case of for \code{Re.Q} and a pure imaginary Quaternion
  for  \code{Im.Q}
}
\details{
  A quaternion has a real and an imaginary part. The real part
  correponds to its first entry and the imaginary part is a quaternion
  with the first entry 0 and all others as before.
}
\examples{
x <- Q(1,2,3,4)
x
Re(x)
Im(x)
Re(x) <- 4
Im(x) <- 1:3
x
}
\author{K.Gerald van den Boogaart}
\keyword{arith}