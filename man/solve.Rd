\name{solve}
\alias{\%*\%.Q}
\alias{\%*\%}
\alias{\%*\%.default}
\alias{solve.Q}
\alias{norm.Q}
\alias{norm}
\title{Linear Algebra with quaternions}
\description{
  The linear algebra can be applied to quaternions. 
}
\usage{
norm(x,...)
\S3method{norm}{Q}(x,...,each=FALSE)
\S3method{solve}{Q}(a,b,...)
## a \%*\% b
}
\arguments{
  \item{x}{a quaternion or numeric vector}
  \item{a}{a quaternion or numeric vector}
  \item{b}{a quaternion or numeric vector}
  \item{\dots}{additional arguments}
  \item{each}{a boolean, if true the norm is computed for each element
    of the vectors, otherwise it is computed for the whole vector}
}
\value{
  The solution of the given calculation as a quaterion or numeric.
}
\details{
  \itemize{
    \item{\%*\%}{
      The matrix product of quaternions is defined as classical matrix
      products, where every multiplication is defined according to ordinary
      quaternion multiplication. 
    }
    \item{norm}{
      The norm is defined as its vector norm.
    }
    \item{solve}{
      Returns the solution of \code{a\%*\% x== b}
    }
  }
}
\seealso{
  \code{\link{add.Q}}, \code{\link{Q}}, \code{\link{base}{solve}}
}
\examples{
A <- matrix(Q(rnorm(25),rnorm(25),rnorm(25),rnorm(25)),nrow=5)
x <- matrix(Q(rnorm(25),rnorm(25),rnorm(25),rnorm(25)),nrow=5)
A \%*\% x
solve(A, A \%*\% x) - x
norm(x)
}
\author{K.Gerald van den Boogaart}
\keyword{arith}