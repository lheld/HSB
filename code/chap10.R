#################################################################
### R code for Chapter 10
#################################################################

### Encoding: UTF-8

# Remark: Running chunks numbers 24-28 requires the package INLA,
# which depends on numerous other packages. It may take several minutes to install INLA
# and its dependencies.
# To install INLA, use the following code:
# install.packages("INLA", repos=c(getOption("repos"), 
#                 INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
  

###################################################
### code chunk number 1
###################################################
getOption("SweaveHooks")[["fig"]]()

rem.data <- scan("../data/rem.txt")
  
plot.mc <- function(x, probs=NA, ...){
  if(is.na(sum(probs))){
    plot(x, pch=19, axes=F, ...)
    axis(1)
    axis(2, at=c(1,2), as.character(c("1","2")))
    box()
  }
  else{
    lim <- (c(1,2)-1.5)*1.05+.5
    par(mar=rep(5, 4))
    plot(probs, axes=FALSE, type="b", ..., ylim=lim)
    points(c(1:length(x)), (x-1.5)*1.05+.5, pch=19)
    lines(c(1,length(x)), rep(0.5, 2), lty=2)
    axis(1)
    axis(2, at=lim, as.character(c("1","2")))
    axis(4, at=c(0:10)/10, as.character(c(0:10)/10))
    mtext("posterior probability", side=4, line=2, las=0)
    box()
  }

}

plot.mc(rem.data, xlab="Time (minutes)", ylab="Sleep status")



###################################################
### code chunk number 2 
###################################################
log.lik <- function(theta, x, stationary=TRUE){
  P <- matrix(c(theta[1],1-theta[1],1-theta[2],theta[2]), nrow=2, ncol=2, byrow=T)
  pi <- rep(1,2)%*%solve(diag(2) - P + 1)
  if(stationary==TRUE)
    log.lik <- log(pi[x[1]])
  else
    log.lik <- 0
  for(i in 2:length(x))
    log.lik <- log.lik + log(P[x[i-1], x[i] ])
  return(log.lik)
}

theta <- c(0.5, 0.5)
eps <- 0.0001

res1 <- optim(theta, log.lik, method = "L-BFGS-B", lower = rep(eps, 2), upper = rep(1-eps, 2),
           control = list(fnscale=-1), x=rem.data, stationary=FALSE, hessian=TRUE)
res2 <- optim(theta, log.lik, method = "L-BFGS-B", lower = rep(eps, 2), upper = rep(1-eps, 2),
           control = list(fnscale=-1), x=rem.data, hessian=TRUE)
par1 <- res1$par
par2 <- res2$par
par11 <- c(par1[1], par2[1])
par22 <- c(par1[2], par2[2])

se1 <- diag(sqrt(-solve(res1$hessian)))
se2 <- diag(sqrt(-solve(res2$hessian)))
se11 <- c(se1[1], se2[1])
se22 <- c(se1[2], se2[2])
stat <- c("conditional", "full")

res <- data.frame(stat, par11, se11, par22, se22)
colnames(res) <- c("Likelihood", "$\\hat{p}_{11}$", "$\\se(\\hat{p}_{11})$", "$\\hat{p}_{22}$", "$\\se(\\hat{p}_{11})$")

print(xtable(res, digits=3, label="tab:rem.ml", caption="Conditional and full ML estimates and standard errors of the diagonal entries of $\\mathbf{P}$.", fig=FALSE),include.rownames=FALSE, sanitize.text.function=function(x){x})


###################################################
### code chunk number 3
###################################################
P <- matrix(NA, 2, 2)
P[1,1] <- par2[1]
P[1,2] <- 1-par2[1]
P[2,1] <- 1-par2[2]
P[2,2] <- par2[2]
P2 <- P%*%P
P3 <- P2%*%P
P4 <- P3%*%P
P5 <- P4%*%P
P6 <- P5%*%P
P7 <- P6%*%P
P8 <- P7%*%P
P9 <- P8%*%P
P10 <- P9%*%P



stationary <- function(P){
  K <- nrow(P)
  pi <- rep(1,K)%*%solve(diag(K) - P + 1)
  return(pi)
}
pi <- stationary(P)


###################################################
### code chunk number 4
###################################################
getOption("SweaveHooks")[["fig"]]()

z2 <- c(P[2,1],P2[2,1],P3[2,1],P4[2,1],P5[2,1],P6[2,1],P7[2,1],P8[2,1],P9[2,1],P10[2,1])
z1 <- c(P[1,1],P2[1,1],P3[1,1],P4[1,1],P5[1,1],P6[1,1],P7[1,1],P8[1,1],P9[1,1],P10[1,1])
matplot(c(1:10), cbind(z1, z2), type="b", xlab=math(k), ylab=expression(paste(math(k), "-step forecast probability", sep="")), col=1, lty=1)
lines(c(1,10), rep(pi[1],2), col=1, lty=2)



###################################################
### code chunk number 5
###################################################

library(MASS)
attach(beav2)
hours <- 24*(day-307) + trunc(time/100) + (time%%100)/60

temp0 <- temp[activ==0]
temp1 <- temp[activ==1]
hours0 <- hours[activ==0]
hours1 <- hours[activ==1]

z <- function(alpha){
  return(atanh(alpha))
}

z.inv <- function(z){
  return(tanh(z))
}

z.inv.deriv <- function(z){
  return(1-tanh(z)^2)
}


target <- function(theta, x, cond=TRUE){
  mu <- theta[1]
  alpha <- z.inv(theta[2])
  sigma2 <- exp(theta[3])
  n <- length(x)
  term1 <- -log(sigma2)*(n-1)/2-sum((x[-1]-mu-alpha*(x[-n]-mu))^2)/(2*sigma2)
  if(cond==TRUE)
    term2 <- 0
  else 
    term2 <- (log(1-alpha^2)-log(sigma2))/2-(1-alpha^2)*(x[1]-mu)^2/(2*sigma2)
    return(term1+term2)

}

start <- function(temp){
  alpha <- acf(temp, plot=FALSE)$acf[2]
  res <- c(mean(temp), z.inv(alpha), log(var(temp)*(1-alpha^2)))
  return(res)
}

res0.cond <- optim(start(temp0), target, method="BFGS", control=list(fnscale=-1), x=temp0, cond=TRUE, hessian=TRUE)
res0.full <- optim(start(temp0), target, method="BFGS", control=list(fnscale=-1), x=temp0, cond=FALSE, hessian=TRUE)
res1.cond <- optim(start(temp1), target, method="BFGS", control=list(fnscale=-1), x=temp1, cond=TRUE, hessian=TRUE)
res1.full <- optim(start(temp1), target, method="BFGS", control=list(fnscale=-1), x=temp1, cond=FALSE, hessian=TRUE)

results <- function(res){
  printout <- matrix(NA, nrow=2, ncol=3)
  colnames(printout) <- c("mu", "alpha", "sigma.squared")
  printout[1,] <- (round(c(res$par[1],z.inv(res$par[2]),exp(res$par[3])),3))
  se <- sqrt(diag(solve(-res$hessian)))
  se[2] <- se[2]*z.inv.deriv(res$par[2])
  se[3] <- se[3]*exp(res$par[3])
  printout[2,] <- (round(se,3))
  rownames(printout) <- c("MLE", "SE")
  return(printout)
}


###################################################
### code chunk number 6
###################################################
getOption("SweaveHooks")[["fig"]]()

plot(hours[1:38], temp[1:38], type="l", xlab="Time (hours)", ylab="Temperature", xlim=c(min(hours),max(hours)), ylim=c(min(temp),max(temp)))
lines(rep(15+9/12, 2),c(25,40), lty=2, col=1)
lines(hours[39:length(temp)], temp[39:length(temp)], type="l")



###################################################
### code chunk number 7
###################################################
activity <- c("inside", "inside", "outside", "outside")
lik <- c("conditional", "full","conditional", "full") 
mu <- c(results(res0.cond)[1,1],results(res0.full)[1,1],results(res1.cond)[1,1],results(res1.full)[1,1])
mu.se <- c(results(res0.cond)[2,1],results(res0.full)[2,1],results(res1.cond)[2,1],results(res1.full)[2,1])
alpha <- c(results(res0.cond)[1,2],results(res0.full)[1,2],results(res1.cond)[1,2],results(res1.full)[1,2])
alpha.se <- c(results(res0.cond)[2,2],results(res0.full)[2,2],results(res1.cond)[2,2],results(res1.full)[2,2])
sigma2 <- c(results(res0.cond)[1,3],results(res0.full)[1,3],results(res1.cond)[1,3],results(res1.full)[1,3])
sigma2.se <- c(results(res0.cond)[2,3],results(res0.full)[2,3],results(res1.cond)[2,3],results(res1.full)[2,3])

res <- data.frame(activity, lik, mu, mu.se, alpha, alpha.se, sigma2, sigma2.se)
colnames(res) <- c("Activity", "Likelihood", "$\\hat{\\mu}$", "$\\se(\\hat{\\mu})$", "$\\hat{\\alpha}$", "$\\se(\\hat{\\alpha})$", "$\\hat{\\sigma}^2$", "$\\se(\\hat{\\sigma}^2)$")

print(xtable(res, digits=3, label="tab:beaver.ml", caption="Conditional and full ML estimates with standard errors from separate analyses of beaver body temperature inside and outside the retreat."),include.rownames=FALSE, sanitize.text.function=function(x){x})


###################################################
### code chunk number 8
###################################################

model <- arima(temp, order=c(1,0,0), xreg=activ)
model2 <- arima(temp, order=c(1,0,0))
model3 <- arima(temp, order=c(0,0,0), xreg=activ)



###################################################
### code chunk number 9
###################################################
covariate <- c("yes", "no", "yes")
autoregression <- c("yes", "yes", "no")

se1 <- sqrt(diag((model[[3]])))
se2 <- sqrt(diag((model2[[3]])))
se3 <- sqrt(diag((model3[[3]])))

mu <- c(model[[1]][2],model2[[1]][2],model3[[1]][1])
mu.se <- c(se1[2],se2[2],se3[1])
ar1 <- c(model[[1]][1],model2[[1]][1],NA)
ar1.se <- c(se1[1],se2[1],NA)
activityEffect <- c(model[[1]][3],NA,model3[[1]][2])
activityEffect.se <- c(se1[3],NA,se3[2])
aic <- c(model$aic,model2$aic,model3$aic) 

res <- data.frame(covariate, autoregression, mu, mu.se, ar1, ar1.se, activityEffect, activityEffect.se, aic)
colnames(res) <- c("Covariate", "Autoregression", "$\\hat{\\mu}$", "$\\se(\\hat{\\mu})$",  "$\\hat{\\alpha}$", "$\\se(\\hat{\\alpha})$", "$\\hat{\\beta}$", "$\\se(\\hat{\\beta})$", "AIC")

print(xtable(res, digits=2, label="tab:beaver.ml2", caption="ML estimates and standard errors of parameters describing beaver body temperature with activity as binary covariate."),include.rownames=FALSE, sanitize.text.function=function(x){x})



###################################################
### code chunk number 10: predictBeaver
###################################################
## predict the next four hours in 10 min intervals
  n.ahead <- 6*4
  p <- predict(model, newxreg=rep(1, n.ahead), n.ahead=n.ahead)
  pred <- p$pred
  pred.se <- p$se
  round(pred, 3)
  round(pred.se, 3)


###################################################
### code chunk number 11
###################################################
getOption("SweaveHooks")[["fig"]]()

plot(hours[1:38], temp[1:38], type="l", xlab="Time (hours)", ylab="Temperature", xlim=c(min(hours),max(hours)+n.ahead/6), ylim=c(min(temp),max(temp)))
lines(rep(15+9/12, 2),c(25,40), lty=2, col=1)
lines(hours[39:length(temp)], temp[39:length(temp)], type="l")
lines(rep(max(hours)+1/12, 2),c(25,40), lty=2, col=1)

lines(max(hours)+c(1:n.ahead)/6, pred, lty=6, type="l")
lines(max(hours)+c(1:n.ahead)/6, pred+1.96*pred.se, lty=2, type="l")
lines(max(hours)+c(1:n.ahead)/6, pred-1.96*pred.se, lty=2, type="l")

lines(c(min(hours),15+8/12),rep(model[[1]][2],2), col="gray60", lty=2)
lines(c(15+10/12,max(hours)+n.ahead/6),rep(model[[1]][2]+model[[1]][3],2), col="gray60", lty=2)



###################################################
### code chunk number 12
###################################################


auxiliary <- function(y, m, theta){
  n <- length(y)
  probs <- matrix(NA, ncol=m, nrow=n)
  for(i in 1:m)
    probs[,i] <- theta[i,y]
  return(probs)
}

## computes MAP estimate for a hidden Markov model with transition matrix P 
## and misclassification probabilities theta based on time series y with 
## number of states K and initial distribution delta
## if delta==NULL the stationary distribution implied by P is used

viterbi <- function(y, K, theta, P, delta=NULL){
    if(is.null(delta)) 
        delta <- solve(t(diag(K)-P+1),rep(1,K))
    n <- length(y)
    probs <- auxiliary(y, K, theta)
    xi <- matrix(0, n, K)
    foo <- delta*probs[1,]
    xi[1,] <- foo/sum(foo)
    for(i in 2:n){
        foo <- apply(xi[i-1,]*P,2,max)*probs[i,]
        xi[i,] <- foo/sum(foo)
    }
    map <- numeric(n)
    map[n] <- which.max(xi[n,])
    for(i in (n-1):1){
        map[i] <- which.max(P[,map[i+1]]*xi[i,])
    }
    return(map)
}

noisy.channel <- scan("../data/noisyChannel.txt")
P <- matrix(c(0.75,0.25,0.25,0.75), ncol=2, byrow=T)
eps <- 0.2
theta <- matrix(c(1-eps, eps, eps, 1-eps), ncol=2, byrow=T)
(map <- viterbi(noisy.channel, K=2, theta, P))


###################################################
### code chunk number 13
###################################################

logit <- function(x) log(x/(1-x))
expit <- function(x) 1/(1+exp(-x))

HMM.loglik <- function(parvect, y, restrict=FALSE, logit=FALSE, ...){
  n <- length(y) 
  P <- matrix(0, 2, 2)
  pr <- matrix(0, 2, 2)
    if(logit==TRUE)
      parvect <- expit(parvect)
  for(i in 1:2){
    P[i,1] <- parvect[i] 
    P[i,2] <- 1 - P[i,1]
  }
  if(restrict==FALSE){
    for(i in 1:2){
      pr[i,1] <- parvect[2+i]
      pr[i,2] <- 1-pr[i,1]
      }
  }
  else{
    element <-  parvect[3]
    pr[1,1] <- pr[2,2] <- element
    pr[1,2] <- pr[2,1] <- 1-element
  }

  delta <- solve(t(diag(2)-P+1),rep(1,2))
  lscale <- 0
  foo <- delta
  for(i in 1:n){
    foo <- foo%*%P*pr[,y[i]]
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
  }
  result <- -lscale
  return(result)
}


HMM.max <- function(y, theta0, P0, restrict=FALSE, ...){
  if(restrict==FALSE)
    parvect0 <- c(as.vector(t(P0[,1])), as.vector(t(theta0[,1])))
  if(restrict==TRUE)
    parvect0 <- c(as.vector(t(P0[,1])), theta0)
  res <- optim(parvect0, HMM.loglik, y=y, method="BFGS", logit=TRUE, restrict=restrict, ...)
  return(res)
}




###################################################
### code chunk number 14
###################################################


lalphabeta <- function(y, m, theta, P, delta=NULL, report.lbeta=TRUE){
  if(is.null(delta)) 
    delta <- solve(t(diag(m)-P+1),rep(1,m))
  n <- length(y)
  lalpha <- matrix(NA, m, n)
  probs <- auxiliary(y, m, theta)
  foo <- delta*probs[1,]
  sumfoo <- sum(foo) 
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  lalpha[,1] <- log(foo) + lscale
  for(i in 2:n){
    foo <- foo%*%P*probs[i,] 
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo/sumfoo
     lalpha[,i] <- log(foo) + lscale
  }
  if(report.lbeta==FALSE)
    return(list(la=lalpha))    
  if(report.lbeta==TRUE){
    lbeta <- matrix(NA, m, n)
    lbeta[,n] <- rep(0, m)
    foo <- rep(1/m, m)
    lscale <- log(m)
    for(i in (n-1):1){
      foo <- P%*%(probs[i+1,]*foo)
      lbeta[,i] <- log(foo)+lscale
      sumfoo <- sum(foo)
      foo <- foo/sumfoo
      lscale <- lscale + log(sumfoo)
    }
    return(list(la=lalpha, lb=lbeta))    
  }
}

posterior.probs <- function(y, m, theta, P, delta=NULL){
  if(is.null(delta)) 
    delta <- solve(t(diag(m)-P+1),rep(1,m))
  n <- length(y)
  fb <- lalphabeta(y, m, theta, P, delta=delta)
  la <- fb$la
  lb <- fb$lb
  c <- max(la[,n])
  llk <- c + log(sum(exp(la[,n]-c)))
  stateprobs <- matrix(NA, ncol=n, nrow=m)
  for(i in 1:n) 
    stateprobs[,i] <- exp(la[,i]+lb[,i]-llk)
  return(stateprobs)
}
 

pp <- posterior.probs(noisy.channel, m=2, theta, P)



###################################################
### code chunk number 15
###################################################
getOption("SweaveHooks")[["fig"]]()

plot.mc(noisy.channel, probs=pp[2,], ylab="State", xlab="Time")


###################################################
### code chunk number 16
###################################################
res0 <- HMM.max(y=noisy.channel, m=2, theta0=cbind(theta, 1-theta), P0=P, restrict=FALSE, hessian=TRUE)
res0b <- HMM.max(y=noisy.channel, m=2, theta0=theta[1], P0=P, restrict=TRUE, hessian=TRUE)

p11 <- expit(res0$par[1])
p21 <- expit(res0$par[2])
theta11 <- 1-expit(res0$par[3])
theta21 <- expit(res0$par[4])

p11.0 <- expit(res0b$par[1])
p21.0 <- expit(res0b$par[2])
theta11.0 <- 1-expit(res0b$par[3])


###################################################
### code chunk number 17
###################################################
getOption("SweaveHooks")[["fig"]]()

P2 <- matrix(c(p11.0, 1-p11.0, p21.0, 1-p21.0), 2, 2, byrow=T)
theta2 <- matrix(c(theta11.0, 1-theta11.0, 1-theta11.0, theta11.0), 2, 2, byrow=T)


pp2 <- posterior.probs(noisy.channel, m=2, theta2, P2)
plot.mc(noisy.channel, probs=pp2[1,], ylab="State", xlab="Time")



###################################################
### code chunk number 18
###################################################

P <- matrix(c(0.75,0.25,0.25,0.75), ncol=2)
theta <- c(0.9,0.1)
res <- HMM.max(y=rem.data, m=2, theta0=cbind(theta, 1-theta), P0=P, restrict=FALSE)

res2 <- HMM.max(y=rem.data, m=2, theta0=theta[1], P0=P, restrict=TRUE)

P.ML <- matrix(NA, 2, 2)
P.ML[1,1] <- expit(res$par[1])
P.ML[1,2] <- 1-P.ML[1,1]
P.ML[2,1] <- expit(res$par[2])
P.ML[2,2] <- 1-P.ML[2,1]

theta.ML <- matrix(NA, 2, 2)
theta.ML[1,1] <- expit(res$par[3])
theta.ML[1,2] <- 1-theta.ML[1,1]
theta.ML[2,1] <- expit(res$par[4])
theta.ML[2,2] <- 1-theta.ML[2,1]



P.ML.0 <- matrix(NA, 2, 2)
P.ML.0[1,1] <- expit(res2$par[1])
P.ML.0[1,2] <- 1-P.ML.0[1,1]
P.ML.0[2,1] <- expit(res2$par[2])
P.ML.0[2,2] <- 1-P.ML.0[2,1]

theta.ML.0 <- matrix(NA, 2, 2)
theta.ML.0[1,1] <- expit(res2$par[3])
theta.ML.0[1,2] <- 1-theta.ML.0[1,1]
theta.ML.0[2,1] <- theta.ML.0[1,2]
theta.ML.0[2,2] <- theta.ML.0[1,1]

ll.diff <- res2$value-res$value
pvalue <- 1-pchisq(ll.diff, df=1)



###################################################
### code chunk number 19
###################################################
getOption("SweaveHooks")[["fig"]]()

P.ML <- matrix(NA, 2, 2)
P.ML[1,1] <- expit(res2$par[1])
P.ML[1,2] <- 1-P.ML[1,1]
P.ML[2,1] <- expit(res2$par[2])
P.ML[2,2] <- 1-P.ML[2,1]

theta.ML <- matrix(NA, 2, 2)
theta.ML[1,1] <- expit(res2$par[3])
theta.ML[1,2] <- 1-theta.ML[1,1]
theta.ML[2,1] <- 1-theta.ML[1,1]
theta.ML[2,2] <- theta.ML[1,1]
pp2 <- posterior.probs(rem.data, m=2, theta.ML, P.ML)
map.P2 <- viterbi(rem.data, K=2, theta.ML, P.ML)
plot.mc(rem.data, probs=pp2[2,], ylab="State", xlab="Time")


###################################################
### code chunk number 20
###################################################
burnin <- 100
nsamples <- 1000


###################################################
### code chunk number 21
###################################################

update.P <- function(x, prior.P){
  counts <- matrix(0.0, nrow=2, ncol=2)
  n <- length(x)
  for(i in 2:n)
    counts[x[i-1],x[i]] <- counts[x[i-1],x[i]] + 1
  P <- matrix(NA, 2, 2)
  P[1,1] <- rbeta(1, prior.P[1,1]+counts[1,1], prior.P[1,2]+counts[1,2])
  P[1,2] <- 1 - P[1,1]
  P[2,1] <- rbeta(1, prior.P[2,1]+counts[2,1], prior.P[2,2]+counts[2,2])
  P[2,2] <- 1 - P[2,1]
  return(P)
}

update.theta <- function(x, y, prior.theta){
   counts <- matrix(0.0, nrow=2, ncol=2)
  n <- length(x)
  for(i in 1:n)
    counts[x[i],y[i]] <- counts[x[i],y[i]] + 1
  theta <- matrix(NA, 2, 2)
  theta[1,1] <- rbeta(1, prior.theta[1,1]+counts[1,1], prior.theta[1,2]+counts[1,2])
  theta[1,2] <- 1 - theta[1,1]
  theta[2,1] <- rbeta(1, prior.theta[2,1]+counts[2,1], prior.theta[2,2]+counts[2,2])
  theta[2,2] <- 1 - theta[2,1]
  return(theta)
 }


update.x <- function(y, P, theta){
  n <- length(y)
  fb <- lalphabeta(y, m=2, theta, P, report.lbeta=FALSE)
  la <- fb$la
  probs <- exp(la[,n])
  probs <- probs/sum(probs)
  x <- rep(NA, n)
   x[n] <- rbinom(1, size=1, prob=probs[2])+1
   for(i in (n-1):1){
     probs <- exp(la[,i])
     probs <- probs*P[,x[i+1]]
     probs <- probs/sum(probs)
     x[i] <- rbinom(1, size=1, prob=probs[2])+1
   }
   return(x)
 }


x <- rem.data
y <- rem.data
prior.P <- matrix(c(9,1,1,9),2,2)
prior.theta <- matrix(c(9,1,1,9),2,2)
theta <- matrix(rep(0.5,4), 2, 2)
P <- matrix(rep(0.5,4), 2, 2)

n <- length(y)
x.sum <- rep(0.0,n)
P.sum <- matrix(0, 2, 2)
theta.sum <- matrix(0, 2, 2)
P11 <- rep(NA, nsamples)
P22 <- rep(NA, nsamples)
theta11 <- rep(NA, nsamples)
theta22 <- rep(NA, nsamples)


for(i in -burnin:nsamples){
  theta <- update.theta(x, y, prior.theta)
  P <- update.P(x, prior.P)
  x <- update.x(y, P, theta)
  if(i > 0){
  x.sum <- x.sum+x
  P.sum <- P.sum+P
  theta.sum <- theta.sum+theta
  P11[i] <- P[1,1]
  P22[i] <- P[2,2]
  theta11[i] <- theta[1,1]
  theta22[i] <- theta[2,2]
}
}


###################################################
### code chunk number 22
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(2,2))
truehist(P11, xlim=c(0.5,1), xlab=math(p[1][1]), col="gray")
truehist(P22, xlim=c(0.5,1), xlab=math(p[2][2]), col="gray")
truehist(1-theta11, xlim=c(0.0,0.5), xlab=math(theta[1]), col="gray")
truehist(1-theta22, xlim=c(0.0,0.5), xlab=math(theta[2]), col="gray")


###################################################
### code chunk number 23
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow=c(1,1))

plot.mc(rem.data, probs=x.sum/nsamples-1, ylab="State", xlab="Time")


###################################################
### code chunk number 24
###################################################

library("INLA")

len <- length(rem.data)
time <- c(1:len)
mydata <- list(y=rem.data-1, t=time)

formula1 <- y ~ f(t, model="ar1", hyper=list(rho=list(initial=5, prior="normal", c(0, 0.15)), prec=list(prior="loggamma", param=c(1.0, .005))))
model1 <- inla(formula1, family="binomial", data=mydata, control.predictor=list(compute=TRUE))

npred <- 10
mydata.pred <- list(y=c(rem.data-1, rep(NA, npred)), t=c(time, c((len+1):(len+npred))))
model1.pred <- inla(formula1, family="binomial", data=mydata.pred, control.predictor=list(compute=TRUE), Ntrials=rep(1, len+npred))



###################################################
### code chunk number 25
###################################################
getOption("SweaveHooks")[["fig"]]()

# desired summary estimates
quantities <- c("mean", "0.025quant", "0.975quant")

matplot(time, model1$summary.linear.predictor[, quantities], 
  type="l", lty=rep(c(1,2,2),2), col=rep(c(1,2),each=3), xlab="Time", ylab=expression(x))



###################################################
### code chunk number 26
###################################################
getOption("SweaveHooks")[["fig"]]()
# transformed linear predictor with inverse link

matplot(time, model1$summary.fitted.values[, quantities], 
  type="l", lty=rep(c(1,2,2),2), col=rep(c(1,2),each=3), xlab="Time", ylab="Probability", ylim=c(0,1))
x <- rem.data
points(c(1:length(x)), (x-1.5)*1.05+.5, pch=19)



###################################################
### code chunk number 27
###################################################
getOption("SweaveHooks")[["fig"]]()

prior.density <- function(rho){
  eta <- log((1+rho)/(1-rho))
  return(2/((1+rho)*(1-rho))*dnorm(eta, 0, sd=sqrt(1/0.15)))
}

# fit a spline through the marg. posterior from INLA to get a smoother line
plot(inla.smarginal(model1$marginals.hyperpar[[2]]), type="l", xlim=c(0.75, 1.0), 
     xlab=math(alpha), ylab="Density")
grid <- c(75000:99999)/100000
lines(grid, prior.density(grid), lty=2, col=1)
# adapt formatting of legend to Fig. 6.1
legend("topleft", legend=rev(c("prior", "posterior")), lty=rev(c(2,1)), bty="n")


###################################################
### code chunk number 28
###################################################
getOption("SweaveHooks")[["fig"]]()
# number of timepoints to predict
npred <- 10
time <- c(1:(len+npred))
mydata <- list(y=c(rem.data-1, rep(NA,npred)), t=time)

formula1 <- y ~ f(t, model="ar1", hyper=list(rho=list(initial=5)))
model1 <- inla(formula1, family="binomial", data=mydata, Ntrials=rep(1, npred+len), control.predictor=list(compute=TRUE))

formula2 <- y ~ f(t, model="rw1")
model2 <- inla(formula2, family="binomial", data=mydata, Ntrials=rep(1, npred+len), control.predictor=list(compute=TRUE))

# desired summary estimates
quantities <- c("mean", "0.025quant", "0.975quant")

# fitted values need to be generated explicitly for predictions
fitted.values1 <- matrix(NA, ncol=3, nrow=len+npred)
fitted.values2 <- matrix(NA, ncol=3, nrow=len+npred)
fitted.values1[1:len,] <- as.matrix(model1$summary.fitted.values[1:len, quantities])
fitted.values2[1:len,] <- as.matrix(model2$summary.fitted.values[1:len, quantities])
for(i in 1:npred){
  transMarg1 <- inla.tmarginal(function(x){exp(x)/(1 + exp(x))}, model1$marginals.linear.predictor[[len+i]])
  fitted.values1[len+i,] <- c(inla.emarginal(function(x) x,  transMarg1), inla.qmarginal(c(0.025,0.975),transMarg1))
  transMarg2 <- inla.tmarginal(function(x){exp(x)/(1 + exp(x))}, model2$marginals.linear.predictor[[len+i]])
  fitted.values2[len+i,] <- c(inla.emarginal(function(x) x,  transMarg2), inla.qmarginal(c(0.025,0.975),transMarg2))
}

matplot(time, fitted.values1, 
  type="l", lty=rep(c(1,2,2),2), col=rep(c(1,2),each=3), xlab="Time", ylab="Probability", ylim=c(0,1.0))
x <- rem.data
points(c(1:length(x)), (x-1.5)*1.05+.5, pch=19)
lines(rep(120.5, 2), c(0,1), col=1, lty=3)




