\name{FGST.EM}
\alias{FGST.EM}
\title{FGST.EM function}
\usage{
FGST.EM(y, w, xi, s, la, nu, iter.max=200, tol=10^-6)
}
\description{
Fit the FGST distribution using EM-algorithm
}
\examples{
y <- r.FGST(n=1000, xi=5, s=2, la=c(-5,3), nu=3 )          # Simulating samples from FGST distribution.
FGST.EM(y, xi=5, s=2, la=c(-5,3), nu=3 )                   # Fitted FGST distribution
   }
