\name{FGSCN.EM}
\alias{FGSCN.EM}
\title{FGSCN.EM function}
\usage{
FGSCN.EM(y, w, xi, s, la, nu1, nu2, iter.max=200, tol=10^-6, equal=FALSE)
}
\description{
Fit the FGSCN distribution using EM-algorithm
}
\examples{
y <- r.FGSCN(n=1000, xi=5, s=2, la=c(-5,3), nu1=0.7, nu2=0.7 )          # Simulating samples from FGSCNe distribution.
FGSCN.EM(y, xi=5, s=2, la=c(-5,3), nu1=0.7 , nu2=0.7, equal=TRUE )                   # Fitted FGSCNe distribution
   }
