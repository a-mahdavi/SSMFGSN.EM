\name{FGSSN.EM}
\alias{FGSSN.EM}
\title{FGSSN.EM function}
\usage{
FGSSN.EM(y, w, xi, s, la, nu, iter.max=200, tol=10^-6)
}
\description{
Fit the FGSSN distribution using EM-algorithm
}
\examples{
y <- r.FGSSN(n=1000, xi=5, s=2, la=c(-5,3), nu=3 )          # Simulating samples from FGSSN distribution.
FGSSN.EM(y, xi=5, s=2, la=c(-5,3), nu=3 )                   # Fitted FGSSN distribution
   }
