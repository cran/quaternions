\docType{package}
\name{quaternions}
\alias{quaternions}
\alias{quaternions-package}
\alias{quaternions}
\title{A package for quaternion arithmetics and linear algebra}
\description{
  Quaternions form a skew field and are a generalisation of real and
  complex numbers.
}
\details{
  \tabular{ll}{
    Package: \tab qunaternions\cr
    Type: \tab Package\cr
    Version: \tab 0.1\cr
    Date: \tab 2006-07-20\cr
    License: \tab GPL Version 2 or newer\cr
  }
  The quaternion package is made for arithmetic and linear algebra with
  the skew field of quaternions.
}
\seealso{
  \code{\link{Q}}, \code{\link{\%*\%.Q}}, \code{\link{+.Q}}
}
\examples{
 q1 <- Q(i=1)
 q2 <- Q(j=1)
 q1 * q2  # 1k
}
\author{K.Gerald van den Boogaart <boogaart@uni-greifswald.de}
\keyword{package}