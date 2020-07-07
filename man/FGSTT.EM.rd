\name{FGSTT.EM}
\alias{FGSTT.EM}
\title{FGSTT.EM function}
\usage{
FGSTT.EM(y, w, xi, s, la, nu, iter.max=200, tol=10^-6, equal=FA:SE)
}
\description{
Fit the FGSTT distribution using EM-algorithm
}
\examples{
y <- r.FGSTT(n=1000, xi=5, s=2, la=c(-5,3), nu1=3, nu2=3 )          # Simulating samples from FGSTTe distribution.
FGSTT.EM(y, xi=5, s=2, la=c(-5,3), nu1=3, nu2=3, equal=TRUE )                   # Fitted FGSTTe distribution
   }
