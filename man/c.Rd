\name{c}
\alias{c.Q}
\alias{[.Q}
\alias{[<-.Q}
\alias{rep.Q}
\alias{names.Q}
\alias{names<-.Q}
\title{Concatnation of quaternions}
\description{
  Concatnated vectors of quaternions.
}
\usage{
\method{c}{Q}(x,...)
\method{names}{Q}(x)
\method{names}{Q}(x)<-value
## x[...]<-value
## x[...]
\method{rep}{Q}(x,times,...)
}
\arguments{
  \item{x}{a quaternion vector}
  \item{\dots}{further quaternion vectors or indices}
  \item{value}{the names of the elements}
  \item{times}{number of repetition}
}
\value{
   a quaternion vector
  
}
\details{
  The vectors are concatenated.
}
\seealso{
  \code{\link{Q}}
}
\examples{
q1 <- Q(1:3,4:5)
q2 <- Q(1:3,4:5)
c(q1,q2)
}
\author{K.Gerald van den Boogaart}
\keyword{arith}