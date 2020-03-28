#################################################################
### R code for Chapter 7
#################################################################

### Encoding: UTF-8

###############################################################
### code chunk number 1: kriterienVergleichUeberlebenszeiten
###############################################################
## Weibull: numerisch
weibullLik <- function(mualpha, log = TRUE, subset=seq_len(nrow(pbcTreat))){
    mu <- mualpha[1]
    alpha <- mualpha[2]
    loglik <- with (pbcTreat[subset,,drop=FALSE],
                    sum (d * dweibull (time, alpha, mu, log = TRUE) +
                         (1 - d) * pweibull (time, alpha, mu, lower.tail = FALSE, log.p = TRUE)
                         )
                    )
    if (log)
        return (loglik)
    else
        return (exp (loglik))
}
start <- c (1000, 1)
resultWeibull <- optim(start, weibullLik, control=list(fnscale=-1), hessian=TRUE)
stopifnot (resultWeibull$convergence == 0)

aicWeibull <- - 2*resultWeibull$value + 2 * 2
bicWeibull <- - 2*resultWeibull$value + log(pbcTreat.nObs) * 2

## Gamma: numerisch
gammaLik2 <- function(muphi, log = TRUE, subset=seq_len(nrow(pbcTreat))){
    mu <- muphi[1]
    phi <- muphi[2]
    loglik <- with (pbcTreat[subset,,drop=FALSE],
                    sum (d * dgamma (time, shape = mu / phi, scale = phi, log = TRUE) +
                         (1 - d) * pgamma (time, shape = mu / phi, scale = phi, lower.tail = FALSE, log.p = TRUE)
                         )
                    )
    if (log)
        return (loglik)
    else
        return (exp (loglik))
}
start <- c (1000, 1000)
resultGamma <- optim(start, gammaLik2, control=list(fnscale=-1), hessian=TRUE)
stopifnot (resultGamma$convergence == 0)

aicGamma <- - 2 * resultGamma$value + 2 * 2
bicGamma <- - 2 * resultGamma$value + log (nrow (pbcTreat)) * 2

## Exponential: analytisch
lambdaMl <- with (pbcTreat, sum (d) / sum (time))
expLik <- function(lambda, log = TRUE, subset=seq_len(nrow(pbcTreat))){
    loglik <- with (pbcTreat[subset,,drop=FALSE],
                    sum (d * dexp (time, lambda, log = TRUE) +
                         (1 - d) * pexp (time, lambda, lower.tail = FALSE, log.p = TRUE)
                         )
                    )
    if (log)
        return (loglik)
    else
        return (exp (loglik))
}
expMaxLogLik <- expLik (lambdaMl)

bicExp <- - 2 * expMaxLogLik + log (nrow (pbcTreat)) * 1
aicExp <- - 2 * expMaxLogLik + 2 * 1

result <- data.frame (
                      model = Cs (Exp, Weibull, Gamma),
                      AIC = c (aicExp, aicWeibull, aicGamma),
                      BIC = c (bicExp, bicWeibull, bicGamma)
                      )


##########################################################
### code chunk number 2: model-selection-weibull-log-lik
##########################################################
start <- c(1000, 1)
resultWeibull <- optim(start, weibullLik, log=TRUE,
                       control=list(fnscale=-1), hessian=TRUE)
logLikWeibull <- resultWeibull$value
logLikWeibull


#########################################################
### code chunk number 3: model-selection-lr-test-stats
#########################################################
testWeibullVsExp <- 2 * (logLikWeibull - expMaxLogLik)
pvalWeibullVsExp <- pchisq(q=testWeibullVsExp, df=1, lower.tail=FALSE)


########################################################
### code chunk number 4: model-selection-gamma-log-lik
########################################################
## implement the (log) likelihood in the gamma model:
gammaLik <- function(muphi, log = TRUE){
    mu <- muphi[1]
    phi <- muphi[2]
    loglik <- with (pbcTreat,
                    sum(d * dgamma (time, 
                                    shape = mu / phi, 
                                    scale = phi, 
                                    log = TRUE) +
                        (1 - d) * pgamma (time, 
                                          shape = mu / phi, 
                                          scale = phi, 
                                          lower.tail = FALSE, 
                                          log.p = TRUE)))
    if(log)
    {        
        return(loglik)
    } else {
        return(exp(loglik))
    }
}
start <- c (1000, 1000)
resultGamma <- optim(start, gammaLik, control=list(fnscale=-1), hessian=TRUE)
logLikGamma <- resultGamma$value
logLikGamma


######################################################################
### code chunk number 5: model-selection-lr-test-stats-gamma-vs-exp
######################################################################
testGammaVsExp <- 2 * (logLikGamma - expMaxLogLik)
pvalGammaVsExp <- pchisq(q=testGammaVsExp, df=1, lower.tail=FALSE)


###################################################
### code chunk number 6
###################################################
q <- c(1:8)
p <- pchisq(2*q, df=q, lower.tail=FALSE)
p <- formatPval(p)

mat <- rbind(q, p)
s1 <- c("$\\boldsymbol{q}$","$\\boldsymbol{P}$\\textbf{-value}")
myoutput <- cbind(s1, mat)


###################################################
### code chunk number 7: latexMyOutput2
###################################################
latexMatrix (myoutput)


############################################################################
### code chunk number 8: model-selection-survival-show-equiv-of-cv-and-aic
############################################################################
## Weibull
cvLogliksWeibull <- numeric(pbcTreat.nObs)
for(i in seq_len(pbcTreat.nObs))
{
    subsetMle <- optim(c(1000, 1),
                       weibullLik,
                       log=TRUE, 
                       subset=
                       setdiff(x=seq_len(pbcTreat.nObs),
                               y=i),
                       control=list(fnscale=-1), 
                       hessian=TRUE)
    stopifnot(subsetMle$convergence==0)
    subsetMle <- subsetMle$par
    
    cvLogliksWeibull[i] <- weibullLik(subsetMle,
                                      log=TRUE,
                                      subset=i)
}
meanCvLoglikWeibull <- mean(cvLogliksWeibull)
scaledAicWeibull <- aicWeibull / (- 2 * pbcTreat.nObs)

## gamma
cvLogliksGamma <- numeric(pbcTreat.nObs)
for(i in seq_len(pbcTreat.nObs))
{
    subsetMle <- optim(c(1000, 1),
                       gammaLik2,
                       log=TRUE, 
                       subset=
                       setdiff(x=seq_len(pbcTreat.nObs),
                               y=i),
                       control=list(fnscale=-1), 
                       hessian=TRUE)
    stopifnot(subsetMle$convergence==0)
    subsetMle <- subsetMle$par
    
    cvLogliksGamma[i] <- gammaLik2(subsetMle,
                                      log=TRUE,
                                      subset=i)
}
meanCvLoglikGamma <- mean(cvLogliksGamma)
scaledAicGamma <- aicGamma / (- 2 * pbcTreat.nObs)

## exponential
cvLogliksExp <- numeric(pbcTreat.nObs)
for(i in seq_len(pbcTreat.nObs))
{
    subsetMle <- with (pbcTreat[setdiff(x=seq_len(pbcTreat.nObs),
                                        y=i),,drop=FALSE], 
                       sum (d) / sum (time))
        
    cvLogliksExp[i] <- expLik(subsetMle,
                              log=TRUE,
                              subset=i)
}
meanCvLoglikExp <- mean(cvLogliksExp)
scaledAicExp <- aicExp / (- 2 * pbcTreat.nObs)


###########################################################
### code chunk number 9: aic-and-crossvalidation-example
###########################################################
## Copy function gammaLik to gammaLik2, add an argument "subset" and
## replace "pbcTreat" by "pbcTreat[subset, , drop=FALSE]" in the function. 
cvLogliksGamma <- numeric(nrow(pbcTreat))
for(i in seq_len(nrow(pbcTreat)))
{
    ## compute the MLE from the training sample:
    subsetMle <- optim(c(1000, 1),
                       gammaLik2,
                       log=TRUE, 
                       subset=
                       setdiff(x=seq_len(nrow(pbcTreat)),
                               y=i),
                       control=list(fnscale=-1), 
                       hessian=TRUE)
    stopifnot(subsetMle$convergence==0)
    subsetMle <- subsetMle$par
    ## compute the log-likelihood for the validation observation:
    cvLogliksGamma[i] <- gammaLik2(subsetMle,
                                   log=TRUE,
                                   subset=i)
}
meanCvLoglikGamma <- mean(cvLogliksGamma)


###################################################
### code chunk number 10: BF-interpretation
###################################################
getOption("SweaveHooks")[["fig"]]()
n <- 100
N <- 100
cor.range <- seq(-1, 1, length.out = n)
mat <- matrix(rep(cor.range, N), ncol = N, byrow = TRUE)
cols <- colorRampPalette(c("grey20", "white"))
                                      
names <- rev(c(1, 3, 20, 150))

y <- seq(0, 1, length = 5)
posy <- y[-5] + 0.125
x <- 0.1

par(mar = c(1, 5, 1, 1))
image(z = mat, col = cols(n), axes = FALSE, bty = "n", ylab = math(BF[12],font=1), main = "")
axis(2, at = y[-1], labels = names, las = 1, lwd = 0, lwd.ticks = 1)

f.lines <- function(pos, x = 0.1, const = 0.01)
{
    lines(rep(x, 2), pos, col = "white")
    lines(c(x-const, x+const), rep(pos[1], 2), col = "white")
    lines(c(x-const, x+const), rep(pos[2], 2), col = "white")
}

eps <- c(0.005, -0.005)
f.lines(y[1:2] + eps, x = x-0.05)
f.lines(y[2:3] + eps, x = x-0.05)
f.lines(y[3:4] + eps, x = x-0.05)
f.lines(y[4:5] + eps, x = x-0.05)

text(x, posy[4], "slight evidence", col = "grey50", cex = 1.5, font = 2, adj = 0)
text(x, posy[3], "positive evidence", col = "white", cex = 1.5, font = 2, adj = 0)
text(x, posy[2], "strong evidence", col = "white", cex = 1.5, font = 2, adj = 0)
text(x, posy[1], "very strong evidence", col = "white", cex = 1.5, font = 2, adj = 0)


###################################################
### code chunk number 11: ModellwahlHardyWeinberg
###################################################
  log.marg.lik.binom <- function(x, alpha){
    n <- sum(x)
    logzaehler <- lgamma(n+1)+lgamma(2*x[1]+x[2]+alpha[1])+lgamma(2*x[3]+x[2]+alpha[2])+x[2]*log(2.)+lgamma(sum(alpha))
    lognenner <- sum(lgamma(x+1))+sum(lgamma(alpha))+lgamma(2*n+sum(alpha))
    return((logzaehler-lognenner))
  }

  # allgemeine Trinomial
  log.marg.lik.multinom <- function(x, alpha){
    n <- sum(x)
    alphastar <- alpha+x
    result <- lgamma(sum(alpha))-lgamma(sum(alphastar))
    result <- result + sum(lgamma(alphastar)) - sum(lgamma(alpha))
    result <- result +lgamma(n+1) - sum(lgamma(x+1))
    return(result)
  }

  # Hardy"=Weinberg mit \upsilon = 1/2
  log.lik <- function(x){
    n <- sum(x)
    logzaehler <- lgamma(n+1)+x[2]*log(2)+2*n*log(1/2)
    lognenner <- sum(lgamma(x+1))
    return((logzaehler-lognenner))
  }

  # data from Mourant 1954
  x <- c(233,385,129)

  lm1 <- rep(NA, 3)
  lm1[1] <- log.marg.lik.binom(x, rep(1, 2))
  lm1[2] <-log.marg.lik.multinom(x, alpha=rep(1, length(x)))
  lm1[3] <- log.lik(x)

  lbf1 <- matrix(0, nrow=3, ncol=3)
  for(i in 1:3)
  for(j in 1:3)
  lbf1[i,j] <- lm1[i]-lm1[j]


  lm2 <- rep(NA, 3)
  lm2[1] <- log.marg.lik.binom(x, rep(.5, 2))
  lm2[2] <-log.marg.lik.multinom(x, alpha=rep(.5, length(x)))
  lm2[3] <- log.lik(x)

  lbf2 <- matrix(0, nrow=3, ncol=3)
  for(i in 1:3)
  for(j in 1:3)
  lbf2[i,j] <- lm2[i]-lm2[j]


###################################################
### code chunk number 12
###################################################
sigma <- as.numeric(sd(iten[["tf500a1"]]))
kappa <- 1/sigma^2


###################################################
### code chunk number 13: m2LogMarginal
###################################################
x <- iten$tf500a1[iten$Sex == "w"]
y <- iten$tf500a1[iten$Sex == "m"]

## Funktion logmarglik1 berechnet den Logarithmus der
## marginalen Likelihood im Modell M_1

logmarglik1 <- function(x, kappa, lambda, nu)
{
    SS <- sum((x-mean(x))^2)
    n <- length(x)
    
    -n/2*log(2*pi) + 
        .5*(n*log(kappa)+log(lambda)-log(n*kappa+lambda)) - 
            .5*(n*kappa*lambda)/(n*kappa+lambda)*(mean(x)-nu)^2-.5*kappa*SS
}

## Funktion postmean1 berechnet den posteriori
## Erwartungswert im Modell M_1
postmean1 <- function(x, kappa, lambda, nu)
{
    n <- length(x)
    (lambda*nu+n*kappa*mean(x))/(lambda+n*kappa)
}

## Funktion logmarglik2 berechnet den Logarithmus der
## marginalen Likelihood im Modell M_2
logmarglik2 <- function(x1, x2, kappa, lambda, nu)
{
    logmarglik1(x=x1, 
                kappa=kappa, lambda=lambda, nu=nu) +
        logmarglik1(x=x2, 
                    kappa=kappa, lambda=lambda, nu=nu)
}

## Funktion postmean berechnet den posteriori
## Erwartungswert im Modell M_2 
postmean2 <- function(x1, x2, kappa, lambda, nu)
{
    c(postmean1(x=x1, 
                kappa=kappa, lambda=lambda, nu=nu),
      postmean1(x=x2, 
                kappa=kappa, lambda=lambda, nu=nu))    
}


## Funktion logmarglik3 berechnet den Logarithmus der
## marginalen Likelihood im Modell M_3

logmarglik3 <- function(x, mu, kappa)
{
    sum(dnorm(x=x, mean=mu, sd=sqrt(1 / kappa), log=TRUE))
}

## Hyperparameter

lambda <- 200^(-2)
nu <- 2000

## Berechnung von posteriori mean und log marginaler
## Likelihood im Modell M_1
xAndy <- c(x, y)
lml1 <- logmarglik1(xAndy, kappa, lambda, nu)
muhat1 <- postmean1(xAndy, kappa, lambda, nu)

## Berechnung von posteriori mean und log marginaler
## Likelihood im Modell M_2
lml2 <- logmarglik2(x, y, kappa, lambda, nu)
muhat2 <- postmean2(x, y, kappa, lambda, nu)

## Berechnung von log marginaler Likelihood im Modell M_3
lml3 <- logmarglik3(xAndy, mu=nu, kappa=kappa)

bf21 <- exp(lml2-lml1)
bf23 <- exp(lml2-lml3)

minlm <- min(lml1,lml2)
rlml1 <- lml1-minlm
rlml2 <- lml2-minlm

p <- exp(c(rlml1,rlml2))
p <- p/sum(p)
save(p, file="../data/postprob.RData")

muhat <- matrix(NA, nrow=2, ncol=3)
muhat[,1] <- rep(muhat1, 2)
muhat[,2] <- c(muhat2[1], muhat2[2])
muhat[,3] <- muhat[,1] * p[1] + muhat[,2]*p[2] 
muhat <- formatRound(muhat,1)

s1 <- c("Female","Male")

myoutput <- cbind(s1, muhat)


###################################################
### code chunk number 14: alco-lindleys-paradox
###################################################
lambda.new <- 10^(-10)
lml1.new <- logmarglik1(xAndy, kappa, lambda.new, nu)
lml2.new <- logmarglik2(x, y, kappa, lambda.new, nu)

minlm.new <- min(lml1.new,lml2.new)
rlml1.new <- lml1.new-minlm.new
rlml2.new <- lml2.new-minlm.new

p.new <- exp(c(rlml1.new,rlml2.new))
p.new <- p.new/sum(p.new)


#####################################################
### code chunk number 15: alco-model-comparison-bic
#####################################################
xAndy <- c(x, y)
n <- length(xAndy)
BIC1 <- - 2 * sum(dnorm(xAndy, mean=mean(xAndy), sd=sqrt(1/kappa), log=TRUE)) + 1*log(n)
AIC1 <- BIC1 - 1 * log(n) + 2

BIC2 <- - 2 * sum(dnorm(x, mean=mean(x), sd=sqrt(1/kappa), log=TRUE)) - 
    2 * sum(dnorm(y, mean=mean(y), sd=sqrt(1/kappa), log=TRUE)) + 2 * log(n)
AIC2 <- BIC2 - 2 * log(n) + 4

bml1 <- -BIC1/2
bml2 <- -BIC2/2

minbml <- min(bml1,bml2)

rbml1 <- bml1-minbml
rbml2 <- bml2-minbml

p2 <- exp(c(rbml1,rbml2))
p2 <- p2/sum(p2)




###################################################
### code chunk number 16: latexMyOutput
###################################################
latexMatrix (myoutput)


