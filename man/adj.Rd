\name{adj}
\alias{adj}
\alias{adj.numeric}
\alias{adj.complex}
\alias{adj.Q}
\title{Compute adjugated matrix}
\description{
  Computes the adjungated or transposed matrix. 
}
\usage{
adj(x,...)
\method{adj}{numeric}(x,...)
\method{adj}{complex}(x,...)
\method{adj}{Q}(x,...)
}
\arguments{
\item{x}{The matrix to be adjungated}
\item{\dots}{The matrix to be adjungated}
}
\value{
  The adjungated matrix.
}
\details{
  The adjungated matrix is the conjugated transposed. 
}
\seealso{
  \code{\link{Q}}
}
\examples{
# Adjungation of real matrices
A <- matrix(1:4,nrow=2)
adj(A)
# Adjungation of complex matrices
A <- matrix(c(1+1i,1+2i,2+1i,2+2i),nrow=2)
A
adj(A)
# Adjungation of quaternions
q <- matrix(Q(1:4,1:4),nrow=2)
q
adj(q)
# Full Adjungation of complex quaternions
r <- rnorm(9)+1i*rnorm(9)
i <- rnorm(9)+1i*rnorm(9)
j <- rnorm(9)+1i*rnorm(9)
k <- rnorm(9)+1i*rnorm(9)
A <- matrix(Q(r,i,j,k),nrow=3)
r <- rnorm(3)+1i*rnorm(3)
i <- rnorm(3)+1i*rnorm(3)
j <- rnorm(3)+1i*rnorm(3)
k <- rnorm(3)+1i*rnorm(3)
q1 <- Q(r,i,j,k)
r <- rnorm(3)+1i*rnorm(3)
i <- rnorm(3)+1i*rnorm(3)
j <- rnorm(3)+1i*rnorm(3)
k <- rnorm(3)+1i*rnorm(3)
q2 <- Q(r,i,j,k)
adj(q1) \%*\% (A \%*\% q2) - adj(A \%*\% q1) \%*\%  q2  

}
\author{K.Gerald van den Boogaart}
\keyword{arith}