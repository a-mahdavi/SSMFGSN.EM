r.FGSEPN <- function(n , xi, s, la, nu){
	require(stabledist)
		m <- length(la) ; a=seq(1,2*m-1,by=2)
		y <- 0
		for(i in 1:n){
		x0 <- rnorm(1)
		tau <- rstable(1, alpha=nu, beta=1, gamma = .1, delta = 0, pm = 1)
		x0 <- tau^(-1/2)*x0
		a.x0 <- t(outer(x0,a,'^'))
		x1<- rnorm(1)
		if( x1 < as.numeric(t(la)%*%a.x0) )
		y[i] <- x0
		else
		y[i] <- -x0
		}
		return(xi + s*y)
		}