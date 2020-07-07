FGST.EM <- function(y, xi, s, la, nu, iter.max=200, tol=10^-6){
	m <- length(la) ; a=seq(1,2*m-1,by=2)
	n <- length(y)
        dif <- 1
        count <- 0
	alp <- as.matrix(la/s^a)
	y.xi <- y-xi
	a.y.xi <- t(outer(y.xi,a,'^'))
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	M <- alp.a.y.xi*sqrt((nu+1)/((y.xi/s)^2+nu))
	while ((dif > tol) && (count <= iter.max)) {
# E step
	s1 <- s2 <- (nu+1)/(nu+(y.xi/s)^2)*pt(M*sqrt((nu+3)/(nu+1)),nu+3)/pt(M,nu+1)
	s3 <- 1/pt(M,nu+1)*(M*sqrt((nu+1)/(nu+(y.xi/s)^2))*pt(M*sqrt((nu+3)/(nu+1)),nu+3)+((y.xi/s)^2+nu)^(-1/2)*gamma(nu/2+1)/(gamma(1/2)*gamma((nu+1)/2))*(1+alp.a.y.xi^2/((y.xi/s)^2+nu))^(-nu/2-1))
	LL <- sum(log(2/s*dt(y.xi/s,nu)*pt(M,nu+1))) # log-likelihood function
# MCE steps
	xi <-  optim(xi,function(x){
	y.xi <- y-x
	a.y.xi <- t(outer(y.xi,a,'^'))
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	M <- alp.a.y.xi*sqrt((nu+1)/((y.xi/s)^2+nu))
	return(-sum(log(2/s*dt(y.xi/s,nu)*pt(M,nu+1))))
		},method= "BFGS")$par
	y.xi <- y-xi
	s <- sqrt(1/n*sum(s1*y.xi^2))
	a.y.xi <- t(outer(y.xi,a,'^'))
	s2.a.y.xi <-0 ; s3.a.y.xi <- matrix(0,m,n)
	for(i in 1:n){
	s3.a.y.xi[,i] <- s3[i]*a.y.xi[,i]
	s2.a.y.xi <- s2.a.y.xi + s2[i]*a.y.xi[,i]%*%t(a.y.xi[,i])
			}
	alp <- solve(s2.a.y.xi)%*%rowSums(s3.a.y.xi)
	alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	nu <- optim(nu,function(x){
		M <- alp.a.y.xi*sqrt((x+1)/((y.xi/s)^2+x))
		-sum(log(2/s*dt(y.xi/s,x)*pt(M,x+1)))
		},method="L-BFGS-B",lower=.001,upper=Inf)$par
		M <- alp.a.y.xi*sqrt((nu+1)/((y.xi/s)^2+nu))
	LL.new <- sum(log(2/s*dt(y.xi/s,nu)*pt(M,nu+1)))  # log-likelihood function
	count <- count +1
	dif <- abs(LL.new/LL-1)
	}
	aic <- -2 * LL.new + 2 * (2+m+1)
	bic <- -2 * LL.new + log(n) * (2+m+1)
list(xi=xi, sigma=s, la=alp*s^a, nu=nu, loglik=LL.new, AIC=aic, BIC=bic, iter=count)
}
