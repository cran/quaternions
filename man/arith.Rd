\name{arith}
\alias{+.Q}
\alias{-.Q}
\alias{*.Q}
\alias{/.Q}
\alias{add.Q}
\alias{inv.Q}
\title{Quaternion Arithmetik}
\description{
  Arithmetic operations of Cliffords Quaternions
}
\usage{
## Methods for Q
## x + y
## x - y
## - x
## x * y
## x / y
## x \%*\% y
inv.Q(x,...)
add.Q(x,y)
}
\arguments{
  \item{x}{a quaternion or numeric vector}
  \item{y}{a quaternion or numeric vector}
  \item{\dots}{additional arguments}
}
\value{
  a numeric giving the results
}
\details{
  For details on quaternion arithmetics see \code{\link{Q}}. The
  operations try to simulate the parallel operations on vectors as
  common in R for real and complex vectors in a way that quaternions can
  be handled just in full analogy to full vectors.
  \itemize{
  \item{+}{The sum of quaternions defined element by element:}
    \deqn{(a,b,c,d)+(r,x,y,z)=(a+r,b+x,c+y,d+z)}{(a,b,c,d)+(r,x,y,z)=(a+r,b+x,c+y,d+z)}
  \item{-}{The difference of quaternions} 
    \deqn{(a,b,c,d)-(r,x,y,z)=(a-r,b-x,c-y,d-z)}{(a,b,c,d)-(r,x,y,z)=(a-r,b-x,c-y,d-z)}
  \item{*}{The product of the quaternions}
    \deqn{(a,b,c,d)(w,x,y,z)=(aw+by+cz+dx,ax+bw+cz-dy,ay+wc+bz-dx,az+dw+by-wc)}{(a,b,c,d)(w,x,y,z)=(aw+by+cz+dx,ax+bw+cz-dy,ay+wc+bz-dx,az+dw+by-wc)
    }
 \item{/}{The product of the quaternions}
   \deqn{(a,b,c,d)(w,x,y,z)=(aw+by+cz+dx,ax+bw+cz-dy,ay+wc+bz-dx,az+dw+by-wc)}{
     (a,b,c,d)(w,x,y,z)=(aw+by+cz+dx,ax+bw+cz-dy,ay+wc+bz-dx,az+dw+by-wc)
   }
 \item{/}{The product of the quaternions}
   \deqn{(a,b,c,d)(w,x,y,z)=(a,b,c,d)(w,-x,-y,-z)/(|w|^2+|x|^2+|y|^2+|z|^2)}{
     (a,b,c,d)(w,x,y,z)=(a,b,c,d)(w,-x,-y,-z)/(|w|^2+|x|^2+|y|^2+|z|^2)
   }
 \item{inv}{Inversion is given by}
   \deqn{(a,b,c,d)^{-1}=(a,-b,-c,-d)/(|a|^2+|b|^2+|c|^2+|d|^2)}{
     (a,b,c,d)^{-1}=(a,-b,-c,-d)/(|a|^2+|b|^2+|c|^2+|d|^2)
   }
 }
}
\seealso{
  \code{\link{Q}}
}
\examples{
q1 <- Q(1,2,3,4)
q2 <- Q(1,2,0,0)
q1+q2
q1-q2
q1*q2
q1/q2
A <- matrix(Q(1:4,rnorm(4)),nrow=2)
B <- matrix(Q(rnorm(4),rnorm(4)),nrow=2)
A \%*\% B
norm(q1)
norm(A)
}
\author{K.Gerald van den Boogaart}
\keyword{arith}