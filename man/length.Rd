\name{length}
\alias{length.Q}
\title{Length of a quaternion vector}
\description{
  Returns the length of the quaternion vector.
}
\usage{
 \method{length}{Q}(x)
}
\arguments{
  \item{x}{a quaternion}
}
\value{
  an integer giving the number of quaternions stored in the object.
}
\details{
  Quaternions try to simulate the behavior of vectors and not of
  lists. Thus the length is really the length of the vector simulate by
  this object.
}
\seealso{
  \code{\link{Q}}
}
\examples{
x <- Q(1:4,0,0,0)
length(x)
}
\author{K.Gerald van den Boogaart}
\keyword{arith}