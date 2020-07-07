FGSCN.EM <- function(y, xi, s, la, nu1, nu2, iter.max=200, tol=10^-6, equal=TRUE){
	if(equal==TRUE){
	nu <- nu1
	m <- length(la) ; a=seq(1,2*m-1,by=2)
	n <- length(y)
        dif <- 1
        count <- 0
	alp <- as.matrix(la/s^a)
	y.xi <- y-xi
	a.y.xi <- t(outer(y.xi,a,'^'))
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi);
	alp.a.y.xi[which(alp.a.y.xi<(-30))] <- (-30)
	while ((dif > tol) && (count <= iter.max)) {
# E step
	s1 <- (1-nu+nu*nu^(3/2)*exp((1-nu)*(y.xi/s)^2/2))/(1-nu+nu*nu^(1/2)*exp((1-nu)*(y.xi/s)^2/2))
	s3 <- alp.a.y.xi + dnorm(alp.a.y.xi)/pnorm(alp.a.y.xi)
	LL <- sum(log(2/s*(nu*nu^(1/2)*dnorm(nu^(1/2)*y.xi/s)+(1-nu)*dnorm(y.xi/s))*pnorm(alp.a.y.xi))) # log-likelihood function
# MCE steps
	xi <-  optim(xi,function(x){
	y.xi <- y-x
	a.y.xi <- t(outer(y.xi,a,'^'))
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	return(-sum(log(2/s*(nu*nu^(1/2)*dnorm(nu^(1/2)*y.xi/s)+(1-nu)*dnorm(y.xi/s))*pnorm(alp.a.y.xi))))
		},method= "BFGS")$par
	y.xi <- y-xi
	s <- sqrt(1/n*sum(s1*y.xi^2))
	a.y.xi <- t(outer(y.xi,a,'^'))
	s3.a.y.xi <- matrix(0,m,n)
	for(i in 1:n)
	s3.a.y.xi[,i] <- s3[i]*a.y.xi[,i]
	alp <- solve(a.y.xi%*%t(a.y.xi))%*%rowSums(s3.a.y.xi)
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	alp.a.y.xi[which(alp.a.y.xi<(-30))] <- (-30)
	nu <- optim(nu,function(x){
		nu=x
		-sum(log(2/s*(nu*nu^(1/2)*dnorm(nu^(1/2)*y.xi/s)+(1-nu)*dnorm(y.xi/s))*pnorm(alp.a.y.xi)))
		},method="L-BFGS-B",lower=.0001,upper=1)$par
	LL.new <- sum(log(2/s*(nu*nu^(1/2)*dnorm(nu^(1/2)*y.xi/s)+(1-nu)*dnorm(y.xi/s))*pnorm(alp.a.y.xi))) # log-likelihood function
	count <- count +1
	dif <- abs(LL.new/LL-1)
	}
	nu1 <- nu2 <- nu
	aic <- -2 * LL.new + 2 * (2+m+1)
	bic <- -2 * LL.new + log(n) * (2+m+1)
	}
	else{
	m <- length(la) ; a=seq(1,2*m-1,by=2)
	n <- length(y)
        dif <- 1
        count <- 0
	alp <- as.matrix(la/s^a)
	y.xi <- y-xi
	a.y.xi <- t(outer(y.xi,a,'^'))
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi);
	alp.a.y.xi[which(alp.a.y.xi<(-30))] <- (-30)
	while ((dif > tol) && (count <= iter.max)) {
# E step
	s1 <- (1-nu1+nu1*nu2^(3/2)*exp((1-nu2)*(y.xi/s)^2/2))/(1-nu1+nu1*nu2^(1/2)*exp((1-nu2)*(y.xi/s)^2/2))
	s3 <- alp.a.y.xi + dnorm(alp.a.y.xi)/pnorm(alp.a.y.xi)
	LL <- sum(log(2/s*(nu1*nu2^(1/2)*dnorm(nu2^(1/2)*y.xi/s)+(1-nu1)*dnorm(y.xi/s))*pnorm(alp.a.y.xi))) # log-likelihood function
# MCE steps
	xi <-  optim(xi,function(x){
	y.xi <- y-x
	a.y.xi <- t(outer(y.xi,a,'^'))
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	return(-sum(log(2/s*(nu1*nu2^(1/2)*dnorm(nu2^(1/2)*y.xi/s)+(1-nu1)*dnorm(y.xi/s))*pnorm(alp.a.y.xi))))
		},method= "BFGS")$par
	y.xi <- y-xi
	s <- sqrt(1/n*sum(s1*y.xi^2))
	a.y.xi <- t(outer(y.xi,a,'^'))
	s3.a.y.xi <- matrix(0,m,n)
	for(i in 1:n)
	s3.a.y.xi[,i] <- s3[i]*a.y.xi[,i]
	alp <- solve(a.y.xi%*%t(a.y.xi))%*%rowSums(s3.a.y.xi)
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	alp.a.y.xi[which(alp.a.y.xi<(-30))] <- (-30)
	nu <- optim(c(nu1,nu2),function(x){
		nu1=x[1] ; nu2=x[2]
		-sum(log(2/s*(nu1*nu2^(1/2)*dnorm(nu2^(1/2)*y.xi/s)+(1-nu1)*dnorm(y.xi/s))*pnorm(alp.a.y.xi)))
		},method="L-BFGS-B",lower=.0001,upper=1)$par
	nu1 <- nu[1] ; nu2 <- nu[2]
	LL.new <- sum(log(2/s*(nu1*nu2^(1/2)*dnorm(nu2^(1/2)*y.xi/s)+(1-nu1)*dnorm(y.xi/s))*pnorm(alp.a.y.xi))) # log-likelihood function
	count <- count +1
	dif <- abs(LL.new/LL-1)
	}
	aic <- -2 * LL.new + 2 * (2+m+2)
	bic <- -2 * LL.new + log(n) * (2+m+2)
	}
list(xi=xi, sigma=s, la=alp*s^a, nu1=nu1, nu2=nu2, loglik=LL.new, AIC=aic, BIC=bic, iter=count)
}
