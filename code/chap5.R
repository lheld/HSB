#################################################################
### R code for Chapter 5
#################################################################

### Encoding: UTF-8

###################################################
### code chunk number 1: normLikelihood
###################################################
getOption("SweaveHooks")[["fig"]]()
xq <- round(mean(iten[["tf500a1"]]),2)
s <- round(sd(iten[["tf500a1"]]),2)
n <- nrow(iten)

                                        # Farben?
farbe.profil = ifelse(getOption ("myColour"), 2, 1)
farbe.est = ifelse(getOption ("myColour"), 3, 1)

normal.log.likelihood <- function(n, xq, s, mu, sigma2){

    term1 <- -n/2*log(sigma2)
    term2 <- -1/(2*sigma2)*((n-1)*s^2+n*(xq-mu)^2)

    return(term1+term2)
}

mugrid <- seq(round(xq-3*s/sqrt(n)),round(xq+3*s/sqrt(n)), 
              length=200)
sigma2grid <- seq(round(s^2/1.3),round(s^2*1.3), 
                  length=200)

nlikelihood <- matrix(nrow=length(mugrid), ncol=length(sigma2grid), NA)

for(i in 1:length(mugrid))
    for(j in 1:length(sigma2grid))
        nlikelihood[i,j] <- normal.log.likelihood(n, xq, s, mugrid[i],sigma2grid[j])

nlikelihood <- exp(nlikelihood - max(nlikelihood))
muhat <- xq
sigma2hat <- s^2/n*(n-1)

contour(mugrid, sigma2grid, nlikelihood, xlab=math (mu), ylab=math
        (sigma^2), xaxs = "i", yaxs = "i", labcex=1)
points(muhat, sigma2hat, pch=19)


###################################################
### code chunk number 2: trinomialLikelihood
###################################################
    x <- c(233, 385, 129)
    n <- sum(x)

    # Farben?
    farbe.profil = ifelse(getOption ("myColour"), 2, 1)
    farbe.est = ifelse(getOption ("myColour"), 3, 1)

    gridSize <- 500
    eps <- 1e-5
    p1grid <- seq(eps,1-eps,length=gridSize)
    p2grid <- seq(eps,1-eps,length=gridSize)

    likelihood2 <- matrix(NA, nrow=length(p1grid), ncol=length(p2grid))

    for(i in 1:length(p1grid))
    for(j in 1:length(p2grid))
  if((p1grid[i]+p2grid[j]) < 1)
  likelihood2[i,j] <- dmultinom(x, size=n, prob=c(p1grid[i],p2grid[j], 1-p1grid[i]-p2grid[j]), log=T)
likelihood2 <- likelihood2-max(likelihood2, na.rm=T)
p1hat <- x[1]/n
p2hat <- x[2]/n


###################################################
### code chunk number 3: trinomialLikelihoodPlot
###################################################
getOption("SweaveHooks")[["fig"]]()
    contour(p1grid, p2grid, likelihood2, nlevels=50, xlab=math (pi[1]),
    ylab=math (pi[2]), xaxs = "i", yaxs = "i", labcex=1)
    points(p1hat, p2hat, pch=19)


###################################################
### code chunk number 4: normProfilLikeA
###################################################
getOption("SweaveHooks")[["fig"]]()
     xq <- round(mean(iten[["tf500a1"]]),2)
    s <- round(sd(iten[["tf500a1"]]),2)
    n <- nrow(iten)

      # Farben?
      farbe.profil = ifelse(getOption ("myColour"), 2, 1)
      farbe.est = ifelse(getOption ("myColour"), 3, 1)

      profile <- function(mu, xq, s, n, log=FALSE){
        term <- ((xq-mu)^2 + s^2/n*(n-1))
        if(log)
            return(-n/2 * log(term))
        else
            return(exp(-n/2*log(term)))
      }

      estimated <- function(mu, xq, s, n, log=FALSE){
        term <- (n*(xq-mu)^2 + s^2*(n-1))
        sigma2hat <- s^2/n*(n-1)
        if(log)
            return(-term/(2*sigma2hat))
        else
            return(exp(-term/(2*sigma2hat)))
      }

      profile.target <- function(mu, xq, s, n, maxlogprofile, c){

        return(exp(profile(mu=mu, xq, s, n, log=TRUE) - maxlogprofile)-c)
      }

      muhat <- xq
      sigma2hat <- s^2/n*(n-1)

      profile.mu <- rep(NA, length(mugrid))
      estimated.mu <- rep(NA, length(mugrid))

      for(i in 1:length(mugrid)){
        profile.mu[i] <- profile(mugrid[i], xq, s, n, log=TRUE)
        estimated.mu[i]<- estimated(mugrid[i], xq, s, n, log=TRUE)
      }
      estimated.mu <- exp(estimated.mu - max(estimated.mu))
      maxlogprofile <- max(profile.mu)
      profile.mu <- exp(profile.mu - maxlogprofile)

      contour(mugrid, sigma2grid, nlikelihood, xlab=math (mu), ylab=math
      (sigma^2), xaxs = "i", yaxs = "i", labcex=1)
      points(muhat, sigma2hat, pch=19)
      lines(mugrid, rep(sigma2hat,length(mugrid)),col=farbe.est, lty = 2)

      sigma2hat.profile <- rep(NA, length(mugrid))

      for(i in 1:length(mugrid))
      sigma2hat.profile[i] <- mean((xq-mugrid[i])^2 + s^2/n*(n-1))
      lines(mugrid, sigma2hat.profile,col=farbe.profil)


###################################################
### code chunk number 5: normProfilLikeB
###################################################
getOption("SweaveHooks")[["fig"]]()
      matplot(mugrid, cbind(profile.mu, estimated.mu),type="l",lty=c(1,2),col=c(farbe.profil, farbe.est),
      xlab=math (mu),ylab=math(tilde(L)(mu)))
      myc <- exp(-.5*qchisq(0.95,df=1))
      lines(c(min(mugrid),max(mugrid)), rep(myc,2),lty=4)

      legend("topright", legend=c("Profile", "Estimated", "95% threshold"), lty=c(1,2,4),
      col=c(farbe.profil, farbe.est, 1), bty = "n")

      mylower <- muhat-3*sqrt(sigma2hat/length(data))
      myupper <- muhat

      result <- uniroot(profile.target, xq=xq, s=s, n=n, maxlogprofile = maxlogprofile, c=myc, lower=mylower, upper=myupper)

      mulower <- result$root

      mylower <- muhat
      myupper <- muhat+3*sqrt(sigma2hat/length(data))

      result <- uniroot(profile.target, xq=xq, s=s, n=n, maxlogprofile = maxlogprofile, c=myc, lower=mylower, upper=myupper)

      muupper <- result$root

       tlower <- xq-qnorm(0.975)*s/sqrt(n)
       tupper <- xq+qnorm(0.975)*s/sqrt(n)



###################################################
### code chunk number 6: logOddsLike
###################################################
getOption("SweaveHooks")[["fig"]]()
      # farbe?
      farbe.profil = ifelse(getOption ("myColour"), 2, 1)
      loglik <- function(eta, psi, x, n){
        psi*x[1]+eta*sum(x)-n[1]*log(1+exp(psi+eta))-n[2]*log(1+exp(eta))
      }
      loglik.full <- function(theta, x, n){
        loglik(theta[1], theta[2], x, n)
      }
      psi.profile.location <- function(psi, x, n){

        result <- optim(-2, loglik, psi=psi, x=x, n=n, method="BFGS", control=list(fnscale=-1))

        return(result$par)
      }
      x <- c(6, 2)
      n <- c(108, 103)

      psi.grid = seq(-1, 4, length = 100)
      eta.grid = seq(-7, -2, length = 100)
      whole.grid = expand.grid(eta.grid, psi.grid)
      werte = apply(whole.grid, 1, loglik.full, x = x, n = n)
      erg <- optim(rep(1, 2), loglik.full, method="BFGS", control=list(fnscale=-1), x=x, n=n, hessian=T)
      maxLogLike = erg$value
      contour(psi.grid, eta.grid,
              matrix(werte, ncol = length(eta.grid), byrow = TRUE) - maxLogLike,levels = -c(7:0)/2,
              xlab = math (psi), ylab = math (eta), xaxs = "i", yaxs = "i", labcex=1)
      lines(psi.grid, sapply(psi.grid, psi.profile.location, x = x, n = n), col = farbe.profil)
      points(erg$par[2],erg$par[1], pch = 19)


###################################################
### code chunk number 7: logOddsProfil
###################################################
getOption("SweaveHooks")[["fig"]]()
      # farbe for approx?
      farbe.approx = ifelse(getOption ("myColour"), 3, 1)
      psi.profile <- function(psi, x, n){

        result <- optim(-2, loglik, psi=psi, x=x, n=n, method="BFGS", control=list(fnscale=-1))

        return(result$value)
      }
      x <- c(6, 2)
      n <- c(108, 103)

      result <- optim(rep(1, 2), loglik.full, method="BFGS",
      control=list(fnscale=-1), x=x, n=n, hessian=T)
      ml <- result$par
      se <- sqrt(diag(solve(-result$hessian)))

      y <- n-x

      ml2 <- c(log(x[2]/y[2]), log(x[1]*y[2]/(x[2]*y[1])))
      se2 <- c(sqrt((1/x[2]+1/y[2])), sqrt(sum(1/x+1/y)))

      # Wald-KI for psi ist:
      wald = ml2[2] + c(-1,1)*se2[2]*1.96

      # Profil-Likelihood-KI for psi berechnen:
      # relative Profil-log-likelihood for psi:
      norm.psi.profile = function(psi, x, n){
        y = n-x
        psi.profile(psi, x, n) - psi.profile(log(x[1]*y[2]/(x[2]*y[1])), x, n)
      }
      # Nullstellen mit Threshold finden:
      lower = uniroot(function(psi){norm.psi.profile(psi,x,n)+1.92}, c(-100, ml2[2]))
      upper = uniroot(function(psi){norm.psi.profile(psi,x,n)+1.92}, c(ml2[2], 100))
      lik.int = c(lower$root, upper$root)

      profile <- rep(NA, length(psi.grid))

      for(i in 1:length(psi.grid))
      profile[i] <- psi.profile(psi.grid[i], x=x, n=n)

      profile.psi <- profile <- profile - max(profile)
      approx <- -.5*(psi.grid-ml2[2])^2/(se2[2])^2

      matplot(psi.grid, cbind(profile, approx),
              type="l", xlab = math(psi), ylab = math(tilde(l)(psi)),
              col = c(farbe.profil, farbe.approx), lty = c(1,2), ylim=c(-4,0))
      abline(h=-(qchisq(0.95, df=1))/2, lty = 4)
      abline(v = 0, lty=3)
      legend("topright", legend=c("Profile", "Quadr. approx.", "95% threshold"), lty=c(1,2,4),
      col=c(farbe.profil, farbe.approx,1), bty = "n")


###################################################
### code chunk number 8: weibullVektor
###################################################
weibullLik <- function(mualpha, log = TRUE)
{
    mu <- mualpha[1]
    alpha <- mualpha[2]
    loglik <- with(pbcTreat,
                   sum(d * dweibull (time, alpha, mu, log = TRUE) +
                       (1 - d) * pweibull(time, alpha, mu, lower.tail = FALSE,
                                          log.p = TRUE))) 
    if(log)
        return(loglik)
    else
        return(exp(loglik))
}
start <- c(1000, 1)
result <- optim(start, weibullLik, control=list(fnscale=-1), hessian=TRUE)
(ml <- result$par)
(observedFisher <- - result$hessian)
(observedFisherInv <- solve(observedFisher))
(se <- sqrt(diag(observedFisherInv)))


###################################################
### code chunk number 9: alphaWaldKiBerechnen
###################################################
  alphaKI <- ml[2] + c(-1, 1) * 1.96 * se[2]


###################################################
### code chunk number 10: weibullProfil
###################################################
getOption("SweaveHooks")[["fig"]]()
target.mu <- function(mu, alpha, ...){ # 3-Punkte-Arg. for log
      mualpha <- c(mu, alpha)
      weibullLik (mualpha, ...)
    }

    profile.alpha <- function(alpha, mu=ml[1], ...){
      result <- optim(mu, target.mu, alpha = alpha,
      control=list(fnscale=-1), ...)
      result$value
    }

    alpha <- seq(0.7, 1.3, length = 100)
    profile <- sapply (alpha, profile.alpha, log = TRUE)
    profile <- profile - max(profile)   # Normierung

    plot(alpha, profile, type="l", xlab = math (alpha), ylab = math (tilde(l)[p](alpha)))
    lines(c(min(alpha), max(alpha)), rep(-1.92, 2), col=1, lty=4)
    legend("topright", legend = "95% threshold", col = 1, lty = 4, bty = "n")

    ## compute profile ci
    cutpoint <- -1.92
    fn <- function (alpha) profile.alpha (alpha) - profile.alpha (ml[2]) - cutpoint
    alphaLower <- uniroot (fn, c (alpha[1], ml[2]))$root
    alphaUpper <- uniroot (fn, c (ml[2], alpha[length (alpha)]))$root


###################################################
### code chunk number 11: gammaLikelihood
###################################################
  gammaLik <- function(muphi, log = TRUE){
    mu <- muphi[1]
    phi <- muphi[2]
    alpha <- mu / phi
    beta <- 1 / phi
    loglik <- with (pbcTreat,
    sum (d * dgamma (time, alpha, beta, log = TRUE) +
    (1 - d) * pgamma (time, alpha, beta, lower.tail = FALSE, log.p = TRUE)
    )
    )
    if (log)
    return (loglik)
    else
    return (exp (loglik))
  }
  start <- c (3000, 3000)
  result2 <- optim(start, gammaLik, control=list(fnscale=-1), hessian=TRUE)
  ml2 <- result2$par
  observedFisher2 <- - result2$hessian
  se2 <- sqrt(diag(solve(observedFisher2)))

  ## for Beispiel mit Deltamethode rechnen:
  jacobi = function(mu, phi){
    return(matrix(c(1/phi, 0, -mu/phi^2, -1/phi^2), nrow = 2))
  }

  estimate = jacobi(ml2[1], ml2[2]) %*% solve(observedFisher2) %*% t(jacobi(ml2[1],ml2[2]))
  se3 = sqrt(diag(estimate))

  # as.95%-confidence interval for alpha also
  gammaAlphaKI <- ml2[1] / ml2[2] + c(-1,1) * 1.96 * se3[1]

  ## for späteres Beispiel Loglik, AIC, BIC speichern:
  weibullValues <- c (result$val, 2 * (result$val - 2), 2 * result$val - log (nrow (pbcTreat)) * 2)
  gammaValues <- c (result2$val, 2 * (result2$val - 2), 2 * result2$val - log (nrow (pbcTreat)) * 2)


###################################################
### code chunk number 12: darmkrebs2A
###################################################
      # Farben?
      farbe.profil = ifelse(getOption ("myColour"), 2, 1)
      farbe.est = ifelse(getOption ("myColour"), 3, 1)

      data <- c(37, 22, 25, 29, 34, 49)
      #data <- c(0, 0,1,1,7,34)[6:1]

      k <- length(data)
      n <- sum(data)

      ####################################################
      # probability function of beta-binomial distribution
      ####################################################
      dbetabinom <- function(t, a, b, k){

        res <- exp(lgamma(k+1)-lgamma(t+1)-lgamma(k-t+1)+lbeta(a+t,b+k-t)-lbeta(a,b))
        return(res)
      }

      # function to calculate log-likelihood of truncated beta-binomial model
      # subject to y>0
      # INPUT:
      # theta = c(a, b): Parameter of beta-binomial model
      # counts: vector with counts of observations with y=1, y=2, ...

      ####################################################
      # log-likelihood of truncated beta-binomial distribution
      # with parameters a and b
      ####################################################
      trunc.beta.binom.loglik.a.b <- function(theta, counts){

        a <- theta[1]
        b <- theta[2]

        n <- sum(counts)
        k <- length(counts)
        vec <- c(1:k)
        result <- sum(counts*(lbeta(a+vec,b+k-vec)-lbeta(a,b))) - n*log(1-exp(lbeta(a,b+k)-lbeta(a,b)))
        return(result)
      }

      ####################################################
      # log-likelihood of truncated beta-binomial distribution
      # with parameters xi and r
      ####################################################
      trunc.beta.binom.loglik.xi.r <- function(theta, counts){

        target.repar <- function(b, aplusb, k, xi){

          term1 <- exp(lgamma(b+k)-lgamma(b)+lgamma(aplusb)-lgamma(aplusb+k))
          return(term1 - xi)
        }

        xi <- theta[1]
        r <- theta[2]
        k <- length(counts)

        aplusb <- (1-r)/r
        eps <- 1E-9
        b <- uniroot(target.repar, aplusb=aplusb, k=k, xi=xi, lower=0+eps, upper=1000)
        a <- aplusb - b$root

        return(trunc.beta.binom.loglik.a.b(c(a,b$root), counts))
      }

      ####################################################
      # log-likelihood of truncated beta-binomial distribution
      # with parameter r for fixed xi
      # index 2 to distinguish it from function in mu r parametrisation
      ####################################################

      trunc.beta.binom.loglik.r2 <- function(r, xi, counts){

        return(trunc.beta.binom.loglik.xi.r(c(xi,r), counts))
      }

      ####################################################
      # profile log-likelihood of truncated beta-binomial distribution
      # with parameter xi
      ####################################################

      trunc.beta.binom.profile.loglik.xi <- function(xi, counts){

        r <-  optim(par=c(.5), fn=trunc.beta.binom.loglik.r2, xi=xi, counts=data, control=mycontrol, method = "L-BFGS-B", lower=0+eps, upper=1-eps)$par[1]

        return(trunc.beta.binom.loglik.xi.r(c(xi,r), counts))

      }

      trunc.beta.binom.repar.loglik.sep <- function(r, xi, counts){

        k <- length(counts)

        aplusb <- (1-r)/r
        eps <- 1E-9
        #print(c(r, xi, k))
        b <- uniroot(target.repar, aplusb=aplusb, k=k, xi=xi, lower=0+eps, upper=1000)
        a <- aplusb - b$root

        return(trunc.beta.binom.loglik.a.b(c(a,b$root), counts))
      }

      fnf <- function(a, b, k){
        return(exp(lbeta(a,b+k)-lbeta(a,b)))
      }

      target.repar <- function(b, aplusb, k, xi){
        term1 <- exp(lgamma(b + k)-lgamma(b)+lgamma(aplusb)-lgamma(aplusb + k))
        return(term1 - xi)
      }


      eps <- 1E-9
      mycontrol <- list(fnscale=-1, maxit=100)

      eps <- .01
      bb.ml2 <- optim(par=c(.5,.5), fn=trunc.beta.binom.loglik.xi.r, control=mycontrol, counts=data, method = "L-BFGS-B", lower=0+eps, upper=1-eps, hessian=T)
      xigrid <- bb.ml2$par[1] + c(-30:100)/200

      xi.profile.loglik <- rep(NA, length(xigrid))
      xi.profile.rho.values <- rep(NA, length(xigrid))
      xi.empirical.loglik <- rep(NA, length(xigrid))

      for(i in 1:length(xigrid)){

        temp <- optim(par=.5, fn=trunc.beta.binom.repar.loglik.sep, xi =xigrid[i], control=mycontrol, counts=data, method = "L-BFGS-B", lower=0+eps, upper=1-eps)

        xi.profile.loglik[i] = temp$val - bb.ml2$val
        xi.profile.rho.values[i] = temp$par
        xi.empirical.loglik[i] <-trunc.beta.binom.repar.loglik.sep(bb.ml2$par[2], xigrid[i],counts=data) - bb.ml2$val
      }

      rgrid <- c(33:67)*2/200

      loglik <- matrix(ncol=length(rgrid), nrow=length(xigrid), NA)

      # Plot der gemeinsamen relativen log-likelihood von xi and rho:
      for(i in c(1:length(xigrid)))
      for(j in c(1:length(rgrid)))
      loglik[i,j] <- trunc.beta.binom.loglik.xi.r(c(xigrid[i], rgrid[j]), counts=data)
      loglik <- loglik-max(loglik)


###################################################
### code chunk number 13: darmkrebs2APlot
###################################################
getOption("SweaveHooks")[["fig"]]()
      contour(xigrid, rgrid, loglik, levels = -c(0,0.25,0.5,1.0,1.5,2.0,2.5),
      xlab=math (xi), ylab=math (rho), xaxs = "i", yaxs = "i", labcex=1)
      points(bb.ml2$par[1],bb.ml2$par[2], pch = 19)
      # Kontur der geschätzten and Profil-log-likelihood einzeichnen:
      abline(h = bb.ml2$par[2], lty = 2, col = farbe.est)
      lines(xigrid, xi.profile.rho.values, lty = 1, col = farbe.profil)
      # Legende hinzufügen
      legend("topleft", lty = c(1,2), legend = c("Profile", "Estimated"), bty = "n", col = c(farbe.profil, farbe.est))


###################################################
### code chunk number 14: darmkrebs2B
###################################################
getOption("SweaveHooks")[["fig"]]()
      # Plot der relativen Profil and geschätzten Log"=Likelihood von xi:
matplot(xigrid, cbind(xi.profile.loglik,xi.empirical.loglik),
type="l", ylim=c(-3, 0), ylab="", xlab = math (xi), col = 
c(farbe.profil,farbe.est))
abline(h = -0.5*qchisq(0.95, 1), col = "gray")
legend("topright", lty=c(1,2,1), col = c(farbe.profil,farbe.est,"gray"),
legend=c(expression(paste("Profile", math(tilde(l)[p](xi)))),
expression(paste("Estimated", math(tilde(l)[e](xi)))), "95% threshold"), bty
= "n")
drawml(bb.ml2$par[1], 0, down = TRUE)

      # Berechnung der Anpassungsstatistiken for Beispiel 3 der Darmkrebsserie:
      ml <- optim(par=c(1,1), fn=trunc.beta.binom.loglik.a.b, control=mycontrol, counts=data, method = "L-BFGS-B", lower=0+eps, upper=1e10, hessian=T)

      ml.a <- ml$par[1]
      ml.b <- ml$par[2]

      probs = dbetabinom(1:6, ml.a, ml.b, 6)
      probs = probs/sum(probs)
      fit = 196*probs

      chi3 = sum((data-fit)^2 / fit)
      g3 = 2*sum(data * log(data/fit))


###################################################
### code chunk number 15: weibullKI
###################################################
getOption("SweaveHooks")[["fig"]]()
data <- unc
  n=length(data)
veca=seq(from=0.78, to=1.68, by=0.01)
vecmu=seq(from=780, to=1780, by=3)
  values=matrix(nrow=length(veca), ncol=length(vecmu))
  p=prod(data)
  grid = expand.grid(veca, vecmu)
  values = matrix(apply(grid, 1, function(one){
    a = one[1]
    mu = one[2]
    return((a^n)/(mu^(n*a))*p^(a-1)*exp(-sum((data/mu)^a)))
  }), nrow = length(veca))
  # normieren:
  values = values / max(values)
  # plotten:
  contour(x=veca, y=vecmu, z=values,
          levels = c(2,5, 10, 25, 50, 80, 95)/100,
          xlab=math (alpha), ylab=math (mu), xaxs="i", yaxs="i", labcex=1)
  # KI schraffieren:
  region = contourLines(x = veca, y = vecmu, z = values, levels = 0.05)
  region = region[[1]]
  polygon(region$x, region$y, density = 10, angle = 45, border = NA, col = "gray")


######################################################################
### code chunk number 16: alcohol-study-compute-two-sample-t-test
######################################################################

x1 <- iten$tf500a1[iten$Sex == "w"]
x2 <- iten$tf500a1[iten$Sex == "m"]
res <- t.test(x1, x2, var.equal = TRUE)
p1 <- res$p.value
t1 <- res$statistic
n <- length(x1)+length(x2)
w <- n*log(1+t1^2/(n-2))
p2 <- 1- pchisq(w, df=1)


#######################################################
### code chunk number 17: anpassungstestBerechnungen
#######################################################
x <- c(233,385,129)
q <- (x[1]+.5*x[2])/sum(x)
lq <- function(q, x){
  p <- c(q*q, 2*q*(1-q), (1-q)*(1-q))
  return(sum(x*log(p)))
}
l1 <- lq(q, x)
ls <- function(x){
  p <- x/sum(x)
  return(sum(x*log(p)))
}
l2 <- ls(x)

chi <- function(q, x){
  e <- sum(x)*c(q*q, 2*q*(1-q), (1-q)*(1-q))
  return(sum((x-e)^2/e))
}
chi2 <- chi(q, x)


###############################################################
### code chunk number 18: screeningtestBspMitAnpassungstest
###############################################################
ml <- 0.6240842
data <- c(37, 22, 25, 29, 34, 49)
# chi^2 Anpassungsstatistik
probs <- dbinom(c(1:6), prob=ml, size=6)
# calculate expected cases; normalise it!
fit <- probs*196/sum(probs)
chi1 <- sum((data-fit)^2/fit)
# G^2 Anpassungsstatistik
g1 = 2*sum(data*log(data/fit))


#########################################################################################
### code chunk number 19: compare-profile-with-conditional-log-lik-in-binomial-example
#########################################################################################
getOption("SweaveHooks")[["fig"]]()
## the data
x <- c(6, 2)
n <- c(108, 103)

## build on results from the previous example in the
## profile likelihood section.
farbe.conditional <- ifelse(getOption ("myColour"), 3, 1)

## compute the relative conditional log-likelihood values:
conditional <- numeric(length(psi.grid))
for(i in seq_along(psi.grid))
{
    conditional[i] <- log(dnoncenhypergeom(x=x[1],
                                           n1=n[1],
                                           n2=n[2],
                                           m1=sum(x),
                                           psi=exp(psi.grid[i])))
}
conditional <- conditional - max(conditional)

## plot the difference
matplot(psi.grid, conditional - profile.psi,
     type="l", xlab=math(psi), ylab = 
     math(tilde(l)[c](psi) - tilde(l)[p](psi)))
abline(v=log(x[1]*y[2]/(x[2]*y[1])), col="gray", lty=2)



