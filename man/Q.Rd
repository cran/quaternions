\name{Q}
\alias{Q}
\alias{Q.Q}
\alias{Qr}
\alias{Qi}
\alias{Qj}
\alias{Qk}
\alias{Q.default}
\title{Constructs a quaternion}
\description{
  A quaternion is an element of a skew field, quit similar to the
  complex number, but with three different instead of only one imaginary
  unit. 
}
\usage{
Q(\dots)
\method{Q}{Q}(r,...)
\method{Q}{default}(r=0,i=0,j=0,k=0,...)
Qr
Qi
Qj
Qr
}
\arguments{
  \item{r}{a quaternion or  a numeric vector or tensor, giving the real part}
  \item{i}{typically a numeric or complex vector or tensor giving the first
    imaginary dimension}
  \item{j}{typically a numeric vector or tensor giving the second
    imaginary dimension}
  \item{k}{typically a numeric vector or tensor giving the third
    imaginary dimension}
  \item{\dots}{only used for generics}
  }
\value{
  an quaternion object of class \code{"Q"} containing the given values.
  }
  \details{
    
    The quaternion class tries to simulate the behavior of numeric vectors
    and matrices as good as possible for quaternions.
    There are two types of quaternions: real quaternios, which are
    isomorphic as real vector space to \eqn{R^4} and complex
    quaternions, which are isomorphic as complex vector
    spaces to \eqn{C^4}. Additionally two vectors can be multiplied
    according to the rule:
    \deqn{(a,b,c,d)(w,x,y,z)=(aw+by+cz+dx,ax+bw+cz-dy,ay+wc+bz-dx,az+dw+by-wc)}
    Real or complex numbers \eqn{a} are identifed with the vectors:
    \eqn{(a,0,0,0)}
    further we define the units  \eqn{1=(1,0,0,0)},
    \eqn{i=(0,1,0,0)}, \eqn{j=(0,0,1,0)}, and \eqn{k=(0,0,0,1)} such
    that we can write any quaternion \eqn{(a,b,c,d)} as
    \deqn{a+bi+cj+dk}
    And accordingly using the rules \eqn{i^2=j^2=k^2=-1} and
    \eqn{ij=k=-ji}, \eqn{jk=i=-kj}, \eqn{ki=j=-ik} we can do
    multiplication directly on the new notation. These vectors are
    accessible as the constants
    \code{Qr}, \code{Q.i}, \code{Q.j}, \code{Q.k}
    \cr
    Together with vector
    addtion this new
    multiplication the quaternions form a skew field. This means that
    all axioms of field execept the commutativity of the multiplication
    are given. Vectors of quaternions are stored as a list with entries
    r, x, y, z each containing a vector giving the one of the entries of
    each of
    the quaternions.
    \cr
    For quaternion matrices the dimension attribute is stored in each of
    the list elements.
}
\seealso{
  \link{quaternions}
}
\examples{
Q(1,2,3,4)
q<- Q(1:4,2:5,3:6,4:7)
q
q$r
q$i
q$i <- 1:3
}
\keyword{arith}