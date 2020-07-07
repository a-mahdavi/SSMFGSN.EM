r.FGST <- function(n , xi, s, la, nu){
		m <- length(la) ; a=seq(1,2*m-1,by=2)
		y <- 0
		for(i in 1:n){
		x0 <- rnorm(1)
		x1<- rnorm(1)
		tau <- rgamma(1,nu/2,nu/2)
		x0 <- tau^(-1/2)*x0 ; x1 <- tau^(-1/2)*x1
		a.x0 <- t(outer(x0,a,'^'))
		if( x1 < as.numeric(t(la)%*%a.x0) )
		y[i] <- x0
		else
		y[i] <- -x0
		}
		return(xi + s*y)
		}