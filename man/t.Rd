\name{t}
\alias{t.Q}
\title{Transposed Quaternion matrix}
\description{
  Computes the transposed of a quaternion matrix
}
\usage{
\method{t}{Q}(x)
}
\arguments{
  \item{x}{The quaternion matrix}
}
\value{
  a transposed matrix
}
\details{
  a transposed matrix of quaternions is a formal transposition without
  conjugating elements. 
}
\seealso{
  \code{\link{Conj.Q}}, \code{\link{adj.Q}} 
}
\examples{
A<- matrix(Q(1:4,5:8,9:12,13:16),nrow=2)
A
t(A)
}
\author{K.Gerald van den Boogaart}
\keyword{arith}