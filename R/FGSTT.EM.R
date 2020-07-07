FGSTT.EM <- function(y, xi, s, la, nu1, nu2, iter.max=200, tol=10^-6, equal=TRUE){
	if(equal==TRUE){
	nu <- nu1
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
	s1 <- (nu+1)/(nu+(y.xi/s)^2)
	s2 <- pt(alp.a.y.xi*sqrt((nu+2)/nu),nu+2)/pt(alp.a.y.xi,nu)
	s3 <- 1/pt(alp.a.y.xi,nu)*(alp.a.y.xi*pt(alp.a.y.xi*sqrt((nu+2)/nu),nu+2) +dt(alp.a.y.xi,nu) )
	LL <- sum(log(2/s*dt(y.xi/s,nu)*pt(alp.a.y.xi,nu))) # log-likelihood function
# MCE steps
	xi <-  optim(xi,function(x){
	y.xi <- y-x
	a.y.xi <- t(outer(y.xi,a,'^'))
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	return(-sum(log(2/s*dt(y.xi/s,nu)*pt(alp.a.y.xi,nu))))
		},method= "BFGS")$par
	y.xi <- y-xi
	s <- sqrt(1/n*sum(s1*y.xi^2))
	a.y.xi <- t(outer(y.xi,a,'^'))
	s3.a.y.xi <- s2.a.y.xi <- matrix(0,m,n)
	for(i in 1:n){
	s2.a.y.xi[,i] <- s2[i]*a.y.xi[,i]
	s3.a.y.xi[,i] <- s3[i]*a.y.xi[,i]
		}
	alp <- solve(s2.a.y.xi%*%t(a.y.xi))%*%rowSums(s3.a.y.xi)
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	nu <- optim(nu,function(x){
		if( nu>0  )
		-sum(log(2/s*dt(y.xi/s,x)*pt(alp.a.y.xi,x)))
		else NULL
		},method="L-BFGS-B",lower=.001,upper=Inf)$par
	LL.new <- sum(log(2/s*dt(y.xi/s,nu)*pt(alp.a.y.xi,nu))) # log-likelihood function
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
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	while ((dif > tol) && (count <= iter.max)) {
# E step
	s1 <- (nu1+1)/(nu1+(y.xi/s)^2)
	s2 <- pt(alp.a.y.xi*sqrt((nu2+2)/nu2),nu2+2)/pt(alp.a.y.xi,nu2)
	s3 <- 1/pt(alp.a.y.xi,nu2)*(alp.a.y.xi*pt(alp.a.y.xi*sqrt((nu2+2)/nu2),nu2+2) +dt(alp.a.y.xi,nu2) )
	LL <- sum(log(2/s*dt(y.xi/s,nu1)*pt(alp.a.y.xi,nu2))) # log-likelihood function
# MCE steps
	xi <-  optim(xi,function(x){
	y.xi <- y-x
	a.y.xi <- t(outer(y.xi,a,'^'))
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	return(-sum(log(2/s*dt(y.xi/s,nu1)*pt(alp.a.y.xi,nu2))))
		},method= "BFGS")$par
	y.xi <- y-xi
	s <- sqrt(1/n*sum(s1*y.xi^2))
	a.y.xi <- t(outer(y.xi,a,'^'))
	s3.a.y.xi <- s2.a.y.xi <- matrix(0,m,n)
	for(i in 1:n){
	s2.a.y.xi[,i] <- s2[i]*a.y.xi[,i]
	s3.a.y.xi[,i] <- s3[i]*a.y.xi[,i]
		}
	alp <- solve(s2.a.y.xi%*%t(a.y.xi))%*%rowSums(s3.a.y.xi)
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	nu <- optim(c(nu1,nu2),function(x){
		nu1=x[1] ; nu2=x[2]
		-sum(log(2/s*dt(y.xi/s,nu1)*pt(alp.a.y.xi,nu2)))
		},method="L-BFGS-B",lower=c(0.001,0.001),upper=c(Inf,Inf))$par
	nu1 <- nu[1] ; nu2 <- nu[2]
	LL.new <- sum(log(2/s*dt(y.xi/s,nu1)*pt(alp.a.y.xi,nu2))) # log-likelihood function
	count <- count +1
	dif <- abs(LL.new/LL-1)
	}
	aic <- -2 * LL.new + 2 * (2+m+2)
	bic <- -2 * LL.new + log(n) * (2+m+2)
	}
list(xi=xi, sigma=s, la=alp*s^a, nu1=nu1, nu2=nu2, loglik=LL.new, AIC=aic, BIC=bic, iter=count)
}
