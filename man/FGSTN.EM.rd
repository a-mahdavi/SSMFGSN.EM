\name{FGSTN.EM}
\alias{FGSTN.EM}
\title{FGSTN.EM function}
\usage{
FGSTN.EM(y, w, xi, s, la, nu, iter.max=200, tol=10^-6)
}
\description{
Fit the FGSTN distribution using EM-algorithm
}
\examples{
y <- r.FGSTN(n=1000, xi=5, s=2, la=c(-5,3), nu=3 )          # Simulating samples from FGSTN distribution.
FGSTN.EM(y, xi=5, s=2, la=c(-5,3), nu=3 )                   # Fitted FGSTN distribution
   }
