\name{print.Q}
\alias{print.Q}
\title{Output of quaternions}
\description{
  This prints quaternions in the form of quaternion strings, readly
  arrange in the right way.
}
\usage{
\S3method{print}{Q}(x,...)
}
\arguments{
  \item{x}{a quaternion vector/tensor}
  \item{\dots}{further arguments passed to or from further methods}
  }
\value{
  the input invisible
}
\details{
  Quaternions are printed in their mathematical form
  r\code{+}i\code{i+}j\code{j+}k\code{k}
}
\seealso{
  \code{\link{Q}}
}
\examples{
print( Q(1:100,cumprod(0.1)) )
print( Q(1:100*1i,cumprod(0.1)) )

}
\keyword{arith}