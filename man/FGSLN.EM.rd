\name{FGSLN.EM}
\alias{FGSLN.EM}
\title{FGSLN.EM function}
\usage{
FGSLN.EM(y, w, xi, s, la, iter.max=200, tol=10^-6)
}
\description{
Fit the FGSLN distribution using EM-algorithm
}
\examples{
y <- r.FGSLN(n=1000, xi=5, s=2, la=c(-5,3) )          # Simulating samples from FGSLN distribution.
FGSLN.EM(y, xi=5, s=2, la=c(-5,3) )                   # Fitted FGSLN distribution
   }
