\name{FGSN.EM}
\alias{FGSN.EM}
\title{FGSN.EM function}
\usage{
FGSN.EM(y, w, xi, s, la, iter.max=200, tol=10^-6)
}
\description{
Fit the FGSN distribution using EM-algorithm
}
\examples{
y <- r.FGSN(n=1000, xi=5, s=2, la=c(-5,3) )          # Simulating samples from FGSN distribution.
FGSN.EM(y, xi=5, s=2, la=c(-5,3) )                   # Fitted FGSN distribution
   }
