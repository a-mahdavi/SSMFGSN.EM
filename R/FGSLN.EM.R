FGSLN.EM <- function(y, xi, s, la, iter.max=200, tol=10^-6){
	 f.L=function(x)
	 exp(-abs(x))/2
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
	 s1 <- 1/abs(eta)*besselK(abs(eta),1/2)/besselK(abs(eta),-1/2)
	 s3 <- alp.a.y.xi + dnorm(alp.a.y.xi)/pnorm(alp.a.y.xi)
	 LL <- sum(log(2/s*f.L(y.xi/s)*pnorm(alp.a.y.xi))) # log-likelihood function
	 # MCE steps
	 xi <-  optim(xi,function(x){
	 y.xi <- y-x
	 a.y.xi <- t(outer(y.xi,a,'^'))
	 alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	 return(-sum(log(2/s*f.L(y.xi/s)*pnorm(alp.a.y.xi))))
		},method= "BFGS")$par
	 y.xi <- y-xi
	 s <- sqrt(1/n*sum(s1*y.xi^2))
	 a.y.xi <- t(outer(y.xi,a,'^'))
	 s3.a.y.xi <- matrix(0,m,n)
	 for(i in 1:n)
	 s3.a.y.xi[,i] <- s3[i]*a.y.xi[,i]
	 alp <- solve(a.y.xi%*%t(a.y.xi))%*%rowSums(s3.a.y.xi)
	 alp.a.y.xi <- as.numeric(t(alp)%*%a.y.xi)
	 LL.new <- sum(log(2/s*f.L(y.xi/s)*pnorm(alp.a.y.xi)))  # log-likelihood function
	 count <- count +1
	 dif <- abs(LL.new/LL-1)
	 }
	 aic <- -2 * LL.new + 2 * (2+m+1)
	 bic <- -2 * LL.new + log(n) * (2+m+1)
	 list(xi=xi, sigma=s, la=alp*s^a, loglik=LL.new, AIC=aic, BIC=bic, iter=count)
	 }
