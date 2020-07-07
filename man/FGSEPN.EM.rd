\name{FGSEPN.EM}
\alias{FGSEPN.EM}
\title{FGSEPN.EM function}
\usage{
FGSEPN.EM(y, w, xi, s, la, nu, iter.max=200, tol=10^-6)
}
\description{
Fit the FGSEPN distribution using EM-algorithm
}
\examples{
y <- r.FGSEPN(n=1000, xi=5, s=2, la=c(-5,3), nu=0.7 )          # Simulating samples from FGSEPN distribution.
FGSEPN.EM(y, xi=5, s=2, la=c(-5,3), nu=0.7 )                   # Fitted FGSEPN distribution
   }
