 r.FGSLN <- function(n , xi, s, la){
	 m <- length(la) ; a=seq(1,2*m-1,by=2)
	 y <- 0
	 for(i in 1:n){
	 x0 <- rnorm(1)*sqrt(rgamma(1,1,1/2))
	 a.x0 <- t(outer(x0,a,'^'))
	 x1<- rnorm(1)
	 if( x1 < as.numeric(t(la)%*%a.x0) )
	 y[i] <- x0
	 else
	 y[i] <- -x0
	 }
	 return(xi + s*y)
	 }