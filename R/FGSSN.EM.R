FGSSN.EM <- function(y, xi, s, la, nu, iter.max=200, tol=10^-6){
	f.s <- function(eta,nu)
		nu*gamma(nu+1/2)*2^nu/(s*sqrt(pi))*pgamma(eta^2/2,nu+1/2)/abs(eta)^(2*nu+1)
	m <- length(la) ; a=seq(1,2*m-1,by=2)
	n <- length(y)
        dif <- 1
        count <- 0
	alp <- as.matrix(la/s^a)
	y.xi <- y-xi
	a.y.xi <- t(outer(y.xi,a,'^'))
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	while ((dif > tol) && (count <= iter.max)) {
# E step
	eta <- y.xi/s
	s1 <- (2*nu+1)/eta^2*pgamma(eta^2/2,nu+3/2)/pgamma(eta^2/2,nu+1/2)
	s3 <- alp.a.y.xi + dnorm(alp.a.y.xi)/pnorm(alp.a.y.xi)
	LL <- sum(log(2*f.s(eta,nu)*pnorm(alp.a.y.xi))) # log-likelihood function
# MCE steps
	xi <-  optim(xi,function(x){
	y.xi <- y-x
	a.y.xi <- t(outer(y.xi,a,'^'))
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	eta <- y.xi/s
	return(-sum(log(2*f.s(eta,nu)*pnorm(alp.a.y.xi))))
		},method= "BFGS")$par
	y.xi <- y-xi
	s <- sqrt(1/n*sum(s1*y.xi^2))
	a.y.xi <- t(outer(y.xi,a,'^'))
	s3.a.y.xi <- matrix(0,m,n)
	for(i in 1:n)
	s3.a.y.xi[,i] <- s3[i]*a.y.xi[,i]
	alp <- solve(a.y.xi%*%t(a.y.xi))%*%rowSums(s3.a.y.xi)
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	eta <- y.xi/s
	f.int2=function(nu){
	f <-0
	for(i in 1:n)
	f[i] <- integrate(function(x) log(x)*x^(nu-1/2)*exp(-x),lower=0,upper=eta[i]^2/2 )$value
	f
	}
	s4 <- log(2/eta^2)+1/(pgamma(eta^2/2,nu+1/2)*gamma(nu+1/2))*f.int2(nu)
	nu <- optim(nu,function(x){
		-sum(log(2*f.s(eta,x)*pnorm(alp.a.y.xi)))
		},method="L-BFGS-B",lower=.001,upper=Inf)$par
	LL.new <- sum(log(2*f.s(eta,nu)*pnorm(alp.a.y.xi)))   # log-likelihood function
	count <- count +1
	dif <- abs(LL.new/LL-1)
	}
	aic <- -2 * LL.new + 2 * (2+m+1)
	bic <- -2 * LL.new + log(n) * (2+m+1)
list(xi=xi, sigma=s, la=alp*s^a, nu=nu, loglik=LL.new, AIC=aic, BIC=bic, iter=count)
}
