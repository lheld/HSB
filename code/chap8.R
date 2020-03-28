#################################################################
### R code for Chapter 8
#################################################################

### Encoding: UTF-8

# Remark 1: The first 6 chunks contain pseudo-code only!
# Remark 2: Code chunk number 53 takes about 3 minutes to run.


##########################################################
### code chunk number 1: pseudo-lik-prior (eval = FALSE)
##########################################################
getOption("SweaveHooks")[["fig"]]()
## log.likelihood <- function(theta, data){...}
## log.prior <- function(theta){...}


##########################################################
### code chunk number 2: pseudo-optim-post (eval = FALSE)
##########################################################
getOption("SweaveHooks")[["fig"]]()
## log.unnorm.posterior <- function(theta, data)
##     log.likelihood(theta, data) + log.prior(theta)
## log.unnorm.posterior <- Vectorize(log.unnorm.posterior, "theta")
## result.opt <-  optimize(log.unnorm.posterior, maximum=TRUE, data=...,
##                         lower=..., upper=...)
## post.mode <- result.opt$maximum
## ordinate <- result.opt$objective


###############################################################
### code chunk number 3: pseudo-norm-posterior (eval = FALSE)
###############################################################
getOption("SweaveHooks")[["fig"]]()
## unnorm.posterior <- function(theta, data)
##     exp(log.unnorm.posterior(theta, data) - ordinate)
## 
## norm.const <- integrate(unnorm.posterior, data=..., 
##                         lower=..., upper=...)$value
## norm.posterior <- function(theta, data)
##     unnorm.posterior(theta, data) / norm.const


###################################################
### code chunk number 4: pseudo-mean (eval = FALSE)
###################################################
getOption("SweaveHooks")[["fig"]]()
## post.mean <- integrate(function(theta) theta * norm.posterior(theta, data=...), 
##                        lower=..., upper=...)$value


####################################################################
### code chunk number 5: pseudo-cdf-and-quantile-fun (eval = FALSE)
####################################################################
getOption("SweaveHooks")[["fig"]]()
## post.cdf <- function(x, data)
##     integrate(norm.posterior, data=data, 
##               lower=..., upper=x)$value
## post.quantile <- function(q, data)
##     uniroot(function(x) post.cdf(x, data) - q,
##             lower=..., upper=...)$root


###################################################################
### code chunk number 6: pseudo-median-and-equi-ci (eval = FALSE)
###################################################################
getOption("SweaveHooks")[["fig"]]()
## post.median <- post.quantile(0.5, data=...)
## 
## post.ci <- function(gamma, data)
##     c(post.quantile((1 - gamma) / 2, data),
##       post.quantile((1 + gamma) / 2, data))
## 
## post.95ci <- post.ci(0.95, data=...)


###################################################
### code chunk number 7: darmkrebs4-get-data
###################################################
counts <- colonCancer[-1]


#########################################################
### code chunk number 8: darmkrebs4-trunc-binom-loglik
#########################################################
## Truncated binomial log-likelihood function
## pi: the parameter, the probability of a positive test result
## data: vector with counts Z_1, ..., Z_N
log.likelihood <- function(pi, data)
{
    n <- sum(data)
    k <- length(data)
    vec <- seq_len(k)
    result <- 
        sum(data * (vec * log(pi) + (k - vec) * log(1 - pi))) - 
            n * log(1 - (1 - pi)^k)
    return(result)
}       
log.likelihood <- Vectorize(log.likelihood, "pi")


###########################################################
### code chunk number 9: darmkrebs4-prior-and-unnorm-post
###########################################################
log.prior <- function(pi)
    dbeta(pi, 0.5, 0.5, log=TRUE)

log.unnorm.posterior <- function(pi, data)
    log.likelihood(pi, data) + log.prior(pi)


###################################################
### code chunk number 10: darmkrebs4-norm-post
###################################################
## the data:
counts

## get posterior mode and its density ordinate:
result.opt <-  optimize(log.unnorm.posterior, maximum=TRUE, data=counts,
                        lower=0, upper=1)
post.mode <- result.opt$maximum
post.mode
ordinate <- result.opt$objective

## use that to compute the normalised posterior density:
unnorm.posterior <- function(pi, data)
    exp(log.unnorm.posterior(pi, data) - ordinate)
norm.const <- integrate(unnorm.posterior, data=counts, 
                        lower=0, upper=1)$value
norm.posterior <- function(pi, data)
    unnorm.posterior(pi, data) / norm.const


###########################################################
### code chunk number 11: darmkrebs4-post-mean-and-median
###########################################################
## posterior mean calculation as in the pseudo code:
post.mean <- integrate(function(pi) pi * norm.posterior(pi, data=counts), 
                       lower=0, upper=1)$value
post.mean

## likewise for the cdf and quantile functions, and hence the median:
post.cdf <- function(x, data)
    integrate(norm.posterior, data=data, 
              lower=0, upper=x)$value

## numerical problems occur if we go exactly to the boundaries here,
## therefore go away some small epsilon:
eps <- 1e-10
post.quantile <- function(q, data)
    uniroot(function(x) post.cdf(x, data) - q,
            lower=0 + eps, upper=1 - eps)$root

post.median <- post.quantile(0.5, data=counts)
post.median


###################################################
### code chunk number 12: darmkrebs4-95-ci
###################################################
post.ci <- function(gamma, data)
    c(post.quantile((1 - gamma) / 2, data),
      post.quantile((1 + gamma) / 2, data))

post.95ci <- post.ci(0.95, data=counts)
post.95ci


###################################################
### code chunk number 13: dummy
###################################################
gamma <- (1-post.median)^6
Z0 <- 196*gamma/(1-gamma)


###################################################
### code chunk number 14: darmkrebs4Posteriori
###################################################
getOption("SweaveHooks")[["fig"]]()
mypi <- c(550:700)/1000
post <- norm.posterior(pi=mypi, data=counts)
plot(mypi, post, type="l", xlab=math (pi), ylab = math (f(pi~plain("|")~z)))
drawml(post.mode, 
       norm.posterior(pi=post.mode, data=counts), digi = 4)


###########################################################
### code chunk number 15: darmkrebs5GemeinsamePosteriori
###########################################################
    mycontrol <- list(fnscale=-1, maxit=100)
    eps <- sqrt (.Machine$double.eps)

    #########################################################
    # log-likelihood of truncated beta-binomial distribution
    # with parameters a and b
    #########################################################
    trunc.beta.binom.loglik.a.b <- function(theta, counts){

      a <- theta[1]
      b <- theta[2]

      n <- sum(counts)
      k <- length(counts)
      vec <- c(1:k)
      result <- sum(counts*(lbeta(a+vec,b+k-vec)-lbeta(a,b))) - n*log(1-exp(lbeta(a,b+k)-lbeta(a,b)))
      return(result)
    }


    ##########################################################
    # log-likelihood of truncated beta-binomial distribution
    # with parameters mu and r
    ##########################################################

    trunc.beta.binom.loglik.mu.r <- function(theta, counts){

      mu <- theta[1]
      r <- theta[2]
      k <- length(counts)

      a <- (1-r)/r*mu
      b<- (1-r)/r*(1-mu)

      return(trunc.beta.binom.loglik.a.b(c(a,b), counts))
    }


    data <- colonCancer[-1]


    # compute unnormalised log posterior
    ########################################################################
    log.unnorm.posterior <- function(theta, counts){

      mu <- theta[1]
      r <- theta[2]

      result <- trunc.beta.binom.loglik.mu.r(theta, counts) + dbeta(mu, 1, 1, log=T) + dbeta(r, 1, 1, log=T)

      return(result)
    }

    ########################################################################
    # compute unnormalised posterior
    ########################################################################
    unnorm.posterior <- function(theta, counts){

      post.mode.value <- optim(rep(0.5,2),log.unnorm.posterior, counts=counts, lower=rep(0,2), upper=rep(1,2), max=T)$value

      return(exp(log.unnorm.posterior(theta, counts=counts)-post.mode.value))
    }
    ########################################################################
    # compute normalised posterior
    ########################################################################
    posterior.norm <- function(theta, counts, norm=NA){

      if(is.na(norm))
      {
         integrand <- function(theta)
         {
            unnorm.posterior(theta, counts=data)
         }
         norm <- adaptIntegrate(integrand, 
                                tol=1e-09,
                                lowerLimit=c(0,0), 
                                upperLimit=c(1,1))$integral
      }
      
      return(unnorm.posterior(theta, counts=counts)/norm)
    }
    ########################################################################

    # compute marginal posterior for mu
    ########################################################################
    posterior.norm.mu <- function(mu, myr, counts, norm=NA){
      result <- rep(NA, length(mu))
      for(i in c(1:length(mu))){
        theta <- c(mu[i], myr)
        result[i] <- posterior.norm(theta, counts=counts, norm=norm)
      }
      return(result)
    }


    # compute marginal posterior for r
    ########################################################################
    posterior.norm.r <- function(r, mu, counts, norm=NA){
      result <- rep(NA, length(r))
      for(i in c(1:length(r))){
        theta <- c(mu, r[i])
        result[i] <- posterior.norm(theta, counts=counts, norm=norm)
      }
      return(result)
    }


    ########################################################################

     unnorm.posterior <- function(theta, counts, post.mode.value=NA){

      mu <- theta[1]
      r <- theta[2]

      result <- trunc.beta.binom.loglik.mu.r(theta, counts) + dbeta(mu, 1, 1, log=T) + dbeta(r, 1, 1, log=T)
      return(exp(result))
    }

    mugrid <- seq(0.01, 0.99, length = 200)
    rgrid <- seq(0.01, 0.99, length = 200)

    posterior <- matrix(nrow=length(mugrid), ncol=length(rgrid), NA)

    integrand <- function(theta)
    {
       unnorm.posterior(theta, counts=data)
    }
    norm <- adaptIntegrate(integrand, 
                           tol=1e-09,
                           lowerLimit=c(0,0), 
                           upperLimit=c(1,1))$integral 

    for(i in 1:length(mugrid)){
      for(j in 1:length(rgrid)){
        posterior[i,j] <- posterior.norm(c(mugrid[i],rgrid[j]), counts=data, norm=norm)
      }
    }

save (log.unnorm.posterior,
      trunc.beta.binom.loglik.mu.r,
      trunc.beta.binom.loglik.a.b,
      file = "../data/posterioriDarmkrebs.RData"
      )


#############################################################
### code chunk number 16: darmkrebs5GemeinsamePosterioriPlot
#############################################################
getOption("SweaveHooks")[["fig"]]()
    ### Plot
    contour(mugrid, rgrid, posterior, xlab=math (mu), ylab=math (rho),
      xlim=c(0.17,0.6), ylim=c(0.3,0.7), labcex=1)


##########################################################
### code chunk number 17: darmkrebs5MarginalePosterioriB
##########################################################
posterior.r <- numeric(length(rgrid))
for(j in seq_along(rgrid))
{
    posterior.r[j] <- integrate(posterior.norm.mu, myr=rgrid[j], 
                                norm=norm, counts=data, 
                                lower=0, upper=1, rel.tol=1e-6)[["value"]]
}


#########################################################
### code chunk number 18: darmkrebs5MarginalePosterioriA
#########################################################
posterior.mu <- rep(NA, length(mugrid))
      for(i in 1:length(mugrid)){
        posterior.mu[i] <- integrate(posterior.norm.r, mu=mugrid[i],
                                     norm=norm,counts=data, lower=0,
                                     upper=1, rel.tol=1e-6)$value
      }


##############################################################
### code chunk number 19: darmkrebs5MarginalePosterioriAPlot
##############################################################
getOption("SweaveHooks")[["fig"]]()
      plot(mugrid, posterior.mu, type="l", xlab=math (mu), ylab=math (f*(mu~plain("|")~x)))


##########################################################
### code chunk number 20: darmkrebs5MarginalePosterioriB
##########################################################
getOption("SweaveHooks")[["fig"]]()
      plot(rgrid, posterior.r, type="l", xlab=math (rho), ylab=math (f*(rho~plain("|")~x)))


###################################################
### code chunk number 21: vergleichLaplaceTable
###################################################
#Tabelle  erzeugen mit n, x, LA, wahrem Wert und relativem Fehler 
results = function(n, x){ 
  ewert = (x+1/2)/(n+1) 
  approx = exp((n+1/2)*log(n-1) + (x+1)*log(x+1/2) - (n+3/2)*log(n) - x*log(x-1/2)) 
  relfehler = (approx - ewert)/ewert 
  approx2 = exp((n+1/2)*log(n) + (x+1)*log(x+1) - (n+3/2)*log(n+1) - x*log(x))
  relfehler2 = (approx2 - ewert)/ewert

    
    data.frame(xOverN = as.character(x/n), n = as.character(n), Exact = ewert, Laplace = approx, Error = relfehler, Laplace2 = approx2, Error2 = relfehler2)
  }

xOverN <- rep(c(0.6, 0.8, 1), each = 3)
n <- rep(c(5,20,100), 3)
x <- n * xOverN

  ergebnisse <- results(n, x)

  ergebnisse[, 4] <- paste(formatRound(ergebnisse[, 4], 4),
                           " \\, \\scriptsize{ (",
                           ifelse(round(ergebnisse[, 5],4) == 0, "-", ""),
                           formatRound(ergebnisse[, 5], 4),
                           ")}",
                           sep="")
ergebnisse[, 6] <- paste(formatRound(ergebnisse[, 6], 4),
                           " \\, \\scriptsize{ (",
                         ifelse(round(ergebnisse[, 7],4) == 0, "-", ""),
                           formatRound(ergebnisse[, 7], 4),
                           ")}",
                           sep="")
ergebnisse <- ergebnisse[, -c(5, 7)]

  names(ergebnisse) <- 
    paste ("$\\boldsymbol{", 
           c("x/n", "n", 
             "\\E(\\pi \\given x)",
             "\\hat{\\E}_{1}(\\pi \\given x)",
             "\\hat{\\E}_{2}(\\pi \\given x)"), 
           "}$", sep = "")
  #in Latex konvertieren
w <-
    latex(ergebnisse,                          # was konvertieren?
          file = latexTempFile,
          label = "tab:vergleichLaplace",  # labeln
          cgroup = c("Observation", "Posterior Mean"),    # um Spalten zusammenzufassen
          n.cgroup = c(2,3),
          rowname = NULL,                  # um keine Zeilennamen zu haben
          booktabs = TRUE,                 # um booktabs zu benutzen
           caption = "Comparison of two Laplace approximations $\\hat{\\E}_{1}(\\pi \\given x)$ and $\\hat{\\E}_{2}(\\pi \\given x)$ with the true posterior mean $\\E(\\pi \\given x)$ in a binomial experiment with Jeffreys'
    prior. The corresponding relative error of the approximation is printed in brackets.",
          center = "none",            # statt \begin{center}
          dec = 4,                         # fuer 2 Dezimalstellen
          cgroupTexCmd = "bfseries",
          colnamesTexCmd = "footnotesize\\bfseries",
          collabel.just = Cs(r, r, L, L, L),       # Spaltenueberschriften ausrichten
          col.just = c (Cs (r, r), rep (Cs (D), 3)),
          numeric.dollar = FALSE,
          where = "tbp"                 # default Werte
    )
postLatex (w, widthFactor = 0.8)


###################################################
### code chunk number 22: darmkrebsLaplace1
###################################################
optimObj <- optim(c(0.5, 0.5), log.unnorm.posterior, counts = data,
                  control = list(fnscale = -1), hessian = TRUE)
(mode <- optimObj$par)
curvature <- optimObj$hessian
(logDetCurvature <- as.numeric(determinant(curvature)$modulus))


###################################################
### code chunk number 23: darmkrebsLaplace2
###################################################
log.mu.times.unnorm.posterior <- function (theta, counts)
    log (theta[1]) + log.unnorm.posterior (theta, counts)
muOptimObj <- optim(c(0.5, 0.5), log.mu.times.unnorm.posterior, 
                    counts = data, control = list (fnscale = -1), 
                    hessian = TRUE)
(muMode <- muOptimObj$par)
muCurvature <- muOptimObj$hessian
(muLogDetCurvature <- as.numeric(determinant(muCurvature)$modulus))


###################################################
### code chunk number 24: darmkrebsLaplace3
###################################################
logPosteriorExpectationMu <-
    1/2 * logDetCurvature - 1/2 * muLogDetCurvature +
    muOptimObj$value - optimObj$value
(posteriorExpectationMu <- exp(logPosteriorExpectationMu))


###################################################
### code chunk number 25: darmkrebsLaplace4
###################################################
log.rho.times.unnorm.posterior <- function (theta, counts)
    log (theta[2]) + log.unnorm.posterior (theta, counts)
rhoOptimObj <- optim (c (0.5, 0.5), log.rho.times.unnorm.posterior, counts = data,
                   control = list (fnscale = -1), hessian = TRUE)
rhoMode <- rhoOptimObj$par
rhoCurvature <- rhoOptimObj$hessian
rhoLogDetCurvature <- as.numeric(determinant(rhoCurvature)$modulus)
logPosteriorExpectationRho <-
    1/2 * logDetCurvature - 1/2 * rhoLogDetCurvature +
    rhoOptimObj$value - optimObj$value
posteriorExpectationRho <- exp(logPosteriorExpectationRho)


###################################################
### code chunk number 26: setseed
###################################################
set.seed(09052003)


###################################################
### code chunk number 27: betapostEwertMC1
###################################################

M <- 10000
  theta <- rbeta(M, 4.5, 1.5)
  (Etheta <- mean(theta))

  (se.Etheta <- sqrt(var(theta)/M))

(Ptheta <- mean(theta<0.5))
  (se.Ptheta <- sqrt(var(theta<0.5)/M))



###################################################
### code chunk number 28: diagnostic3
###################################################
M <- 10000
## prev: samples from Beta(1.5, 99.5) distribution
prev <- rbeta(n, 1.5, 104)
## first use fixed values for sensitivity and specificity
sens <- 0.9
spec <- 0.9
## and calculate positive predictive value (PPV)
ppv <- sens * prev / (sens * prev + (1-spec)*(1-prev))
## now assume distributions for sensitivity and specificity
sens <- rbeta(n, 36.5, 4.5)
spec <- rbeta(n, 36.5, 4.5)
## and calculate the resulting samples for PPV
ppv2 <- sens * prev / (sens * prev + (1-spec)*(1-prev))


###################################################
### code chunk number 29: diagnostic4
###################################################
getOption("SweaveHooks")[["fig"]]()
par(pty="s")
boxplot(ppv, ppv2, names=c("(a)", "(b)"), ylab="positive predictive value")


###################################################
### code chunk number 30: alco-gender-comparison
###################################################
x <- iten$tf500a1[iten$Sex == "w"]
y <- iten$tf500a1[iten$Sex == "m"]

location1 <- mean (x)
scale1 <- sd (x) / sqrt (length (x))
df1 <- length (x) - 1

location2 <- mean (y)
scale2 <- sd (y) / sqrt (length (y))
df2 <- length (y) - 1

M <- 1e+4                             # Anzahl Ziehungen
Mformatted <- formatBig(trunc(M))


##########################################################
### code chunk number 31: fig-alco-simulated-differences
##########################################################
getOption("SweaveHooks")[["fig"]]()
set.seed (151)
mu1 <- rst(M, xi = location1, omega = scale1, nu = df1)
mu2 <- rst(M, xi = location2, omega = scale2, nu = df2)
theta2 <- mu1 - mu2
truehist (theta2,
          xlab = math (theta), ylab = math (hat (f)*(theta~plain("|")~x)),
          col = "gray")


###################################################
### code chunk number 32: dic-hardy-weinberg
###################################################
## data:
x <- c(233, 385, 129)
n <- sum(x)

## log-likelihoods:
triLoglik <- function(pi)
{
    if(is.vector(pi))
    {
        pi <- t(pi)
    }
    
    apply(pi, 1, 
          FUN=function(onePi) sum(x * log(onePi)))      
}

hwLoglik <- function(upsilon)
{
    pi <- cbind(upsilon * upsilon,
                2 * upsilon * (1 - upsilon),
                (1 - upsilon) * (1 - upsilon))
    triLoglik(pi)
}

## sample from the Hardy-Weinberg posterior:
aPost <- 1 + 2 * x[1] + x[2]
bPost <- 1 + x[2] + 2 * x[3]
upsilonSamples <- rbeta(n=10000, aPost, bPost)

## calculate DIC:
upsilonBar <- aPost / (aPost + bPost)
upsilonBar - mean(upsilonSamples) ## check: OK
hwDf <- mean(2 * (hwLoglik(upsilonBar) - hwLoglik(upsilonSamples)))
hwDf
hwDic <- - 2 * hwLoglik(upsilonBar) + 2 * hwDf
hwDic

## sample from the Dirichlet posterior:
alphaPost <- 1 + x
piSamples <- rdirichlet(n=10000,
                        alphaPost)
head(piSamples)

## calculate DIC:
piBar <- alphaPost / sum(alphaPost)
piBar - colMeans(piSamples) ## check: OK
triDf <- mean(2 * (triLoglik(piBar) - triLoglik(piSamples)))
triDf
triDic <- - 2 * triLoglik(piBar) + 2 * triDf
triDic


###################################################
### code chunk number 33: dic-hardy-weinberg-save
###################################################
save(hwDf, hwDic, 
     triDf, triDic,
     file="../data/hardyWeinbergDic.RData")


###################################################
### code chunk number 34: betapostEwertMC2
###################################################
## sort parameter samples
thetaorder <- theta[order(theta)]
## determine number and sizes of all possible credible intervals
M <- 10000
level <- 0.95
n.cis <- round(M * (1-level)) + 1
size <- numeric(n.cis)
for(i in seq_len(n.cis)){
    lower <- thetaorder[i]
    upper <- thetaorder[M - n.cis + i]
    size[i] <- upper - lower
}
## get the one with smallest size: the HPD interval
size.min <- which.min(size)
HPD.lower <- thetaorder[size.min]
HPD.upper <- thetaorder[M - n.cis + size.min]
## also compute the equi-tailed interval
ET.lower <- thetaorder[(M * (1-level)) / 2]
ET.upper <- thetaorder[M - (M * (1-level)) / 2]
## compare the results:
c(HPD.lower, HPD.upper)
c(ET.lower, ET.upper)


###################################################
### code chunk number 35: betapostHpd
###################################################
outerdens <- function(h, alpha, beta)
{
    ## compute the mode of this beta distribution
    mode <- max(min((alpha-1)/(alpha + beta - 2),
                    1), 0)
    
    ## to compute the intersection points of the height h
    ## and the density function, only go up to epsilon
    ## to the mode:
    eps <- 1e-15        
    
    ## compute the lower intersection point
    lower <- 
        if(mode <= 0)
        {
            0                           # the density is montonically decreasing 
        } else {
            uniroot(function(x){dbeta(x, alpha, beta) - h}, 
                    interval = c(0, mode - eps))$root
        }
    
    ## compute the upper intersection point
    upper <-
        if(mode >= 1)
        {
            1                           # the density is monotonically increasing 
        } else {
            uniroot(function(x){dbeta(x, alpha, beta) - h}, 
                    interval = c(mode + eps, 1))$root
        }
    
    ## compute the probability inside the interval (lower, upper):
    prob <- pbeta(lower, alpha, beta) + 
        pbeta(upper, alpha, beta, lower.tail = FALSE)
    
    ## return everything
    return(c(prob=prob, 
             lower=lower,
             upper=upper))
}

betaHpd <- function(alpha, 
                    beta, 
                    level=0.95)
{
    ## compute the mode of this beta distribution,
    ## but go epsilon away from 0 or 1
    eps <- 1e-15
    mode <- max(min((alpha-1)/(alpha + beta - 2),
                    1-eps), 0+eps)
    
    ## determine h_opt:
    result <- uniroot(function(h){outerdens(h, alpha, beta)["prob"] - 
                                      (1 - level)}, 
                      ## search in the interval (eps, f(mode) - eps)
                      interval = 
                      c(eps, 
                        dbeta(mode, alpha, beta) - eps))  
    height <- result$root
    
    ## this gives the HPD interval
    hpd <- outerdens(height, alpha, beta)[c("lower", "upper")]
    hpd 
}

alpha <- 4.5
beta <- 1.5
hpd <- betaHpd(alpha, beta)
hpd


###################################################
### code chunk number 36: importanceSampling
###################################################
u <- runif(M)
w <- dbeta(u, 4.5, 1.5)
(sum(w))
(Etheta.u <- sum(u * w) / sum(w))
(se.Etheta.u <- sqrt(sum((u - Etheta.u)^2 * w^2))) / sum(w)
(Ptheta.u <- sum((u < 0.5) * w) / sum(w))
(se.Ptheta.u <- sqrt(sum(((u < 0.5) - Ptheta.u)^2 * w^2))) / sum(w)


###################################################
### code chunk number 37: rejectionSampling
###################################################
## posterior parameters and mode:
alpha <- 4.5
beta <- 1.5
mode <- (alpha - 1) / (alpha + beta - 2)
a <- dbeta(mode, alpha, beta)
## number of samples to be produced:
M <- 10000
## vector where the samples will be stored:
theta <- numeric(M)
## also save the number of trials for each sample:
trials <- numeric(M)
## for each sample:
for(m in seq_along(theta))
{
    k <- 0
    while(TRUE)
    {
        k <- k + 1
        ## sample random variables
        u <- runif(1)
        z <- runif(1)
        ## check for acceptance, then exit the loop
        if(u <= dbeta(z, alpha, beta) / a)
            break
    }
    ## save the z realization as a theta sample
    theta[m] <- z
    ## and the number of trials
    trials[m] <- k
}
## average number of trials required:
mean(trials)
## estimate posterior mean of theta:
(Etheta <- mean(theta))
(se.Etheta <- sqrt(var(theta) / M))
## estimate P(theta < 0.5 | y):
(Ptheta <- mean(theta < 0.5))
(se.Ptheta <- sqrt(var(theta < 0.5) / M))


###################################################
### code chunk number 38: gibbsSamplerDarmkrebs
###################################################
## data set:
fulldata <- colonCancer
k <- 0:6
n <- sum(fulldata[-1])
## impute start value for Z0 (first element):
fulldata[1] <- 10
## MCMC settings:
nburnin <- 100
niter <- 10000
## where the samples will be saved:
pisamples <- Z0samples <- numeric(niter)
## set the random seed:
set.seed(920)
## do the sampling:
for(i in (- nburnin + 1):niter)
{
    ## draw pi from full conditional
    pi <- rbeta(1, 0.5 + sum(k * fulldata), 0.5 + sum((6 - k) * fulldata))
    ## draw Z0 from full conditional
    fulldata[1] <- rnbinom(1, size=n, prob=1-(1-pi)^6)
    ## if the burn-in iterations have passed, save the samples
    if(i>0)
    {
        pisamples[i] <- pi
        Z0samples[i] <- fulldata[1]
    }
}


###################################################
### code chunk number 39: gibbsPathsA
###################################################
getOption("SweaveHooks")[["fig"]]()

plot(pisamples, type = "p", pch = ".", ylim = c(0.5, 0.7), cex = 0.5,
     xlab = math (m), ylab = math (pi^(m))
     )


###################################################
### code chunk number 40: gibbsPathsB
###################################################
getOption("SweaveHooks")[["fig"]]()
      bereich = min(Z0samples):max(Z0samples)
      n = length(Z0samples)
      relfreq = matrix(NA, n, length(bereich))
      for(i in seq(along = bereich)){
        relfreq[,i] = cumsum(Z0samples == bereich[i])/(1:n)
      }

      if(getOption ("myColour")){
        matplot(relfreq, type = "l", ylim = c(0, max(relfreq[-c(1:50)])),
        xlab = math (m), ylab = "Cumulative relative frequency", lty = 1)
      } else
      {
        matplot(relfreq, type = "l", ylim = c(0, max(relfreq[-c(1:50)])),
        xlab = math (m), ylab = "Cumulative relative frequency", col = 1)
      }

      text(par("usr")[2], relfreq[n,], labels = bereich, xpd = TRUE, pos = 2)


#######################################################
### code chunk number 41: summariesDarmkrebsSampling
#######################################################
  summary(pisamples)
  summary(Z0samples)


###################################################
### code chunk number 42: gibbsDistributionsA
###################################################
getOption("SweaveHooks")[["fig"]]()
      truehist(pisamples, xlab = math (pi), col = "gray", ylab=math (hat (f)*(pi~plain("|")~z)))


###################################################
### code chunk number 43: gibbsDistributionsB
###################################################
getOption("SweaveHooks")[["fig"]]()
      barplot(relfreq[n,],
              xlab = math (k), ylab = math (hat(plain(Pr))(Z[0] == k~plain("|")~z)),
              ylim = c(0, max(relfreq[n,])+0.1), names.arg = bereich)


###################################################
### code chunk number 44: dqPrioriA
###################################################
getOption("SweaveHooks")[["fig"]]()
      x <- c(233, 385, 129)

      lower <- function(q){
        lower <- pmax(-q^2, -(1-q)^2)
        return(lower)
      }

      upper <- function(q){
        upper <- q*(1-q)
        return(upper)
      }

      # zunaechst zur Illustration samples aus der Priori

      alpha <- 1
      beta <- 1
      nsamples <- 10000
      qprior <- rbeta(nsamples, alpha, beta)
      dprior <-  runif(nsamples, min=lower(qprior), max=upper(qprior))

      # gemeinsame Verteilung:
      plot(qprior, dprior, xlab = math (q), ylab = math (delta), cex = .25, col = "darkgrey")
      # Grenzen einzeichnen:
      grenzen.wd = 2
      grenzen.col = 1
      qvec = seq(min(qprior), max(qprior), length = 1e3)
      lines(qvec, lower(qvec), col = grenzen.col, lwd = grenzen.wd)
      lines(qvec, upper(qvec), col = grenzen.col, lwd = grenzen.wd)


###################################################
### code chunk number 45: dqPrioriB
###################################################
getOption("SweaveHooks")[["fig"]]()
      truehist(dprior, xlab = math (delta), col = "gray", ylab=math (hat (f)*(delta)))


###################################################
### code chunk number 46: metropolisWithinGibbs
###################################################
## data set
x <- c(233, 385, 129)
## limits for delta as functions of v
lower <- function(v)
{
    lower <- pmax(-v^2, -(1-v)^2)
    return(lower)
}
upper <- function(v)
{
    upper <- v * (1 - v)
    return(upper)
}
## function to compute probabilities from v and d
myprob <- function(v, d)
{
    p1 <- v^2 + d
    p2 <- 2 * v * (1-v) - 2 * d
    p3 <- (1-v)^2 + d
    p <- c(p1, p2, p3)
    ## if the result is valid, return it,
    ## otherwise NAs:
    if(all(p >= 0) & all(p <= 1))
    {
        return(p)
    } else {
        return(rep(NA, 3))
    }
}
## use a uniform prior on v:
alpha <- 1
beta <- 1
## MCMC sampling setup
scale <- 0.03
niter <- 10000
nburnin <- 100
## Initialise samples and counters
vsamples <- numeric(niter)
dsamples <- numeric(niter)
v <- 0.5
d <- 0
vyes <- 0
dyes <- 0
## start MCMC sampling
for(i in (-nburnin + 1):niter)
{
    ## proposal for v:
    first <- alpha + 2 * x[1] + x[2]
    second <- beta + x[2] + 2 * x[3]
    vstar <- rbeta(1, first, second)
    ## is this a valid combination of v and d?
    valid <- (d >= lower(vstar)) && (d <= upper(vstar))
    ## compute the log-posterior ratio
    logPostRatio <-
        if(valid)
        {
            ## from the likelihood
            dmultinom(x, prob=myprob(vstar, d), log=TRUE) - 
                dmultinom(x, prob=myprob(v, d), log=TRUE) +
                    ## from the marginal prior on v
                    dbeta(vstar, alpha, beta, log=TRUE) - 
                        dbeta(v, alpha, beta, log=TRUE) +
                            ## from the prior on d given v
                            dunif(x=d, lower(vstar), upper(vstar), log=TRUE) -
                                dunif(x=d, min=lower(v), max=upper(v), log=TRUE)
        } else {
            ## if the combination is not valid, then the likelihood of 
            ## the proposal is zero.
            - Inf
        }
    ## compute the log proposal ratio
    logPropRatio <- dbeta(v, first, second, log=TRUE) - 
        dbeta(vstar, first, second, log=TRUE)
    ## hence we obtain the log acceptance probability
    logAcc <- logPostRatio + logPropRatio
    ## decide acceptance
    if(log(runif(1)) <= logAcc)
    {
        v <- vstar
        ## count acceptances
        if(i > 0)
        {
            vyes <- vyes + 1
        }
    }

    ## proposal for d:
    first <- max(d - scale, lower(v))
    second <- min(d + scale, upper(v))
    dstar <- runif(1, min=first, max=second)
    ## compute the log posterior ratio
    logPostRatio <- 
        ## from the likelihood
        dmultinom(x, prob=myprob(v, dstar), log=TRUE) - 
            dmultinom(x, prob=myprob(v, d), log=TRUE) +
                ## from the prior on d given v
                dunif(x=dstar, lower(v), upper(v), log=TRUE) -
                    dunif(x=d, min=lower(v), max=upper(v), log=TRUE)
    ## compute the log proposal ratio
    logPropRatio <- dunif(d, first,second, log=TRUE) - 
        dunif(dstar, first, second, log=TRUE)
    ## hence we obtain the log acceptance probability
    logAcc <- logPostRatio + logPropRatio
    ## decide acceptance
    if(log(runif(1)) <= logAcc)
    {
        d <- dstar
        ## count acceptances
        if(i > 0)
        {
            dyes <- dyes + 1
        }
    }
    ## if burnin was passed, save the samples
    if(i > 0)
    {
        vsamples[i] <- v
        dsamples[i] <- d
    }
}


###################################################
### code chunk number 47: hwungleichLikelihood
###################################################
  # ML-Analyse:

  loglik = function(theta){
    v = theta[1]
    d = theta[2]
    p = myprob(v,d)
    if(any(is.na(p)))
    -Inf
    else
    dmultinom(x, prob=myprob(v,d), log=T)
  }
  erg = optim(c(0.6,0), loglik, method = "BFGS", control = list(fnscale = -1), hessian = TRUE)
  ml = erg$par
  se = sqrt(diag(solve(-erg$hessian)))


###########################################################
### code chunk number 48: metropolisWithinGibbsPosterioriA
###########################################################
getOption("SweaveHooks")[["fig"]]()
      # Plot der Posteriori
      plot(vsamples, dsamples, xlab = math (upsilon), ylab = math (delta), cex=0.25, col = "darkgrey")
      # zweidim. Kerndichteschaetzer:
      vd <- cbind(vsamples, dsamples)
      est <- bkde2D(vd, bandwidth=c(0.005,0.005), gridsize = rep(100,2))
      contour(est$x1, est$x2, est$fhat, add = TRUE, nlevels = 10, labcex=1)


###########################################################
### code chunk number 49: metropolisWithinGibbsPosterioriB
###########################################################
getOption("SweaveHooks")[["fig"]]()
      truehist(dsamples, xlab = math (delta), col = "gray", ylab=math (hat (f)*(delta~plain("|")~x)))


########################################################
### code chunk number 50: darmkrebs_PosterioriSampling
########################################################
## load data
 load ("../data/posterioriDarmkrebs.RData") 

counts <- colonCancer[-1]

  ## determine proposal pars
  postOptim <- optim (rep (0.5, 2), fn = log.unnorm.posterior, counts = counts,
  control = list (fnscale = -1), method = "BFGS",
  hessian = TRUE)

  postMode <- postOptim$par
  postVar <- solve (- postOptim$hessian)

  varFactor <- 3
  proposalVar <- varFactor * postVar

  ## create markov chain
  M <- 1e+4
  numAccepted <- 0

  thetaChain <- matrix (nrow = M, ncol = length (postMode))

  set.seed (100)
  current <- postMode
  for (i in seq_len (M)){
    proposal <- mvrnorm (n = 1, mu = current, Sigma = proposalVar)
    logAcceptanceProb <- log.unnorm.posterior (proposal, counts) - log.unnorm.posterior (current, counts)
    if (!is.na (logAcceptanceProb) && log (runif (1)) < logAcceptanceProb){
      current <- proposal
      numAccepted <- numAccepted + 1
    }
    thetaChain[i,] <- current
  }

  ## discard burn-in
  thetaChain <- thetaChain[- (1:1e3), ]

  ## compute false negative fraction samples
  alphaChain <- thetaChain[,1] * (1 / thetaChain[,2] - 1)
  betaChain <- (1 - thetaChain[,1]) * (1 / thetaChain[,2] - 1)

  gammaChain <- exp (
  lbeta (alphaChain, betaChain + length (counts)) -
  lbeta (alphaChain, betaChain)
  )

# compute Z0 samples

Z0Chain <- sum(counts)*gammaChain/(1-gammaChain)


empMode <- function(x){
  mylevel <- 0.01
  hpd <- empiricalHpd (x, mylevel[1])
  return(mean(hpd))
}

empHpd <- empiricalHpd (gammaChain, level = 0.95)

empMode <- empMode(gammaChain)



###################################################
### code chunk number 51: darmkrebs_gammaHistogramm
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (density (gammaChain, from = 0, to = 1, kernel = "epan", adjust = 2),
      xlab = math (xi), ylab = math (hat (f)*(xi~plain("|")~z)),
      main = ""
      )
## mark hpd
rect (empHpd[1], par ("usr")[3], empHpd[2], 0, density = NA, col = "gray")


#############################################################
### code chunk number 52: metropolisWithinGibbsPosterioriZ0
#############################################################
getOption("SweaveHooks")[["fig"]]()
      truehist(Z0Chain, xlab = math (Z[0]), col = "gray", xlim=c(0,300), h=5, ylab = math (plain(Pr)(Z[0] == k~plain("|")~z)))


################################################################
### code chunk number 53: MCMC fuer GMRF-Beispiel mit SIR-Raten
################################################################
      
  
      ## Speicherort fuer Kette?
      chainFile <- "../data/scotland/gmrfChainBinary"

      chainConnection <- file (chainFile, "wb")
      open (chainConnection)

      x <- scotlandData$x
      e <- scotlandData$e

      n <- length (x)

      ## Priori-Parameter
      alpha <- 1
      beta <- 0.01
      prioriPars <- c (alpha = alpha, beta = beta)

      ## Nachbarschaftsdaten laden
      adjacency <- as.matrix (read.table ("../data/scotland/adjacency.txt"))
      colnames (adjacency) <- c ("numberOfNeighbours", paste ("index", 1:(ncol (adjacency)-1), sep = ""))

      ## 0-F?lle rauswerfen
      # x[x == 0] <- 1

      ## Startwerte und Kettenl?nge w?hlen
      lambda <- x / e
      lambda[x == 0] <- 1

      sigma2 <- 1

      M <- 1e+5
      indexVector <- seq_len (n)
      acceptanceVector <- integer (n)

      ## Gibbs-Sampler starten
      set.seed (233)
      inverseGammaShape <- alpha + (n - 1) / 2

      for (t in seq_len (M)){
        for (i in indexVector){ ## Metropolis-Hastings-Schritt f?r ite SIR
          lambdaiOld <- lambda[i]

          ni <- adjacency[i, 1]
          etaiBar <- sum (log (lambda [adjacency[i, -1]])) / ni

          ## eigentlich:

          ## meanLogNormal <- exp (etaiBar + sigma2 / (2 * ni))
          ## varianceLogNormal <- (exp (sigma2 / ni) - 1) * meanLogNormal^2

          ## und daraus dann Berechnung von gammaShape und gammaRate
          ## numerisch besser:

          meanSigma2 <- sigma2 / ni
          if(0){
            logMeanLogNormal <- - log (exp (meanSigma2) - 1)
            logVarianceLogNormal <- logMeanLogNormal - etaiBar - meanSigma2 / 2

            gammaShape <- x[i] + exp (logMeanLogNormal)
            gammaRate <- e[i] + exp (logVarianceLogNormal)
          }
          if(0){
            logMeanLogNormal <-  etaiBar + 0.5*meanSigma2
            logVarianceLogNormal <- meanSigma2 + log(exp(meanSigma2)-1)+2*etaiBar

            gammaShape <- x[i] + exp (2*logMeanLogNormal-logVarianceLogNormal)
            gammaRate <- e[i] + exp (logMeanLogNormal-logVarianceLogNormal)
          }
          if(1){
            MeanLogNormal <-  exp(etaiBar + 0.5*meanSigma2)
            VarianceLogNormal <- exp(meanSigma2)*(exp(meanSigma2)-1)*exp(2*etaiBar)

            gammaShape <- x[i] + MeanLogNormal^2/VarianceLogNormal
            gammaRate <- e[i] + MeanLogNormal/VarianceLogNormal
          }
          lambdaiNew <- rgamma (1, gammaShape, gammaRate)

          logAcceptanceProb <- dlnorm(lambdaiNew, meanlog = etaiBar, sdlog = sqrt(meanSigma2), log=T) -
            dlnorm(lambdaiOld, meanlog = etaiBar, sdlog = sqrt(meanSigma2), log=T) +
            dgamma(lambdaiOld, gammaShape, gammaRate, log=T) -
            dgamma(lambdaiNew, gammaShape, gammaRate, log=T) +
            dpois(x[i], e[i]*lambdaiNew, log=T) -
            dpois(x[i], e[i]*lambdaiOld, log=T)

          if (log (runif (1)) < logAcceptanceProb){
            lambda[i] <- lambdaiNew
            acceptanceVector[i] <- acceptanceVector[i] + 1
          }
        }
        ## Varianz
        eta <- log (lambda)
        pairedEtasSquaredDifferenceSum <-   # hier ein 1/2 hinzugef?gt
          sum (sapply (indexVector, function (i) sum ((eta[i] - eta[adjacency[i, -1]])^2))) / 2

        inverseGammaRate <- beta + pairedEtasSquaredDifferenceSum / 2

        sigma2 <- 1 / rgamma (1, inverseGammaShape, inverseGammaRate)

        ## speichern
        writeBin (c (lambda, sigma2), chainConnection)
      }
      close (chainConnection)

      chainMatrix <- matrix (data = readBin (chainFile, double (0), n = M * (n + 1)),
                             nrow = M, ncol = n + 1,
                             byrow = TRUE)

      # plot (chainMatrix[, 53])                # ein Fall mit x == 0
      # plot (chainMatrix[, 54])                # ein normaler Fall

      ## discard burnin
      chainMatrix <- chainMatrix[-c (1:(M/10)), ]

      lambdaPosteriorMeans <- colMeans (chainMatrix[, -(n+1)])
      sigma2Mean <- mean (chainMatrix[, n+1])


####################################################################
### code chunk number 54: PoissonVollBayes_PosterioriEwerteVonSIRs
####################################################################
getOption("SweaveHooks")[["fig"]]()
############################################################################
##                                                                        ##
##       Schottlandlandkarte grau, 1000 Farbstufen, Log-Skala             ##
##                                                                        ##
############################################################################
## mit Posteriori-Erwartungswerten von lambda_i

scot <- scotlandCoordinates
scotland <- list(length=56)
for(i in 1:56) scotland[[i]] <- scot[scot[,1]==i,]

xm <- (min(scotland[[54]][,2],na.rm=T) + max(scotland[[54]][,2],na.rm=T))/2
xd <- xm - 407.5
ym <- (min(scotland[[54]][,3],na.rm=T) + max(scotland[[54]][,3],na.rm=T))/2
yd <- ym - 940
scotland[[54]][,2] <- scotland[[54]][,2]-xd
scotland[[54]][,3] <- scotland[[54]][,3]-yd

xm <- (min(scotland[[55]][,2],na.rm=T) + max(scotland[[55]][,2],na.rm=T))/2
xd <- xm - 500
ym <- (min(scotland[[55]][,3],na.rm=T) + max(scotland[[55]][,3],na.rm=T))/2
yd <- ym - 902.5
scotland[[55]][,2] <- scotland[[55]][,2]-xd
scotland[[55]][,3] <- scotland[[55]][,3]-yd

## hier einzige AEnderung
sir <- lambdaPosteriorMeans


unten <- 0.25
oben <- 4


# Definieren der Farbmatrix

buntcol <- matrix(nrow=1001, ncol=3)
vari <- 80000
mit <- 650
buntcol[,3] <- 0.0
buntcol[1,1] <- 0.0
buntcol[1,2] <- 0.0
for(i in 1:1000){
buntcol[i+1,1] <- exp(-(i-mit)*(i-mit)/(2*vari))
buntcol[1001-i+1,2] <- exp(-(i-mit)*(i-mit)/(2*vari))
}

mycol <- gray(seq(1,0,length=1001))

farb <- 1:1001

eps <- 1e-10
grenz <- 1:1002
grenz[1] <- -eps
grenz[2] <- unten
grenz[1001] <- oben
grenz[1002] <- 1000000
diff <- log(grenz[1001])-log(grenz[2])
schritt <- diff/999
gleich <- 1:1000
for(i in 1:1000){
gleich[i] <- log(grenz[2]) + (i-1)*schritt
}
for(i in 1:998){
grenz[i+2] <- exp(gleich[i+1])
}

farben <- cut(sir,grenz)
for(i in 1:56){
if(farben[i]==1) farben[i] <- 2}

breite <- c(50,550)
hoehe <- c(500,1000)

par(pty="s", mar = rep (0, 4))
plot(breite,hoehe,type="n",col=0,xlab="",ylab="",axes=F)
for(k in 1:length(scotland)){
polygon(scotland[[k]][,2],scotland[[k]][,3],col=mycol[farben[k]],border=F)
polygon(scotland[[k]][,2],scotland[[k]][,3],density=0,lwd=.3,col=1)}
polygon(c(365,365,450,450),c(890,990,990,890),density=0,lwd=.3,col=1)
polygon(c(450,450,550,550),c(815,990,990,815),density=0,lwd=.3,col=1)

for(i in 1:1001){
polygon(c(440,440,455,455),c(550+.2*(i-1),550+.2*(i),550+.2*(i),550+.2*(i-1)),col=mycol[i],border=F)
if(grenz[i] < .25 && grenz[i+2] >= .25){
text(475,550+.2*(i+1),"0.250",cex=.7,col=1)
lines(c(455,458),c(550+.2*(i+1),550+.2*(i+1)),lwd=.3,col=1001)}
if(grenz[i] < .4 && grenz[i+2] >= .4){
text(475,550+.2*(i+1),"0.400",cex=.7,col=1)
lines(c(455,458),c(550+.2*(i+1),550+.2*(i+1)),lwd=.3,col=1001)}
if(grenz[i] < .625 && grenz[i+2] >= .625){
text(475,550+.2*(i+1),"0.625",cex=.7,col=1)
lines(c(455,458),c(550+.2*(i+1),550+.2*(i+1)),lwd=.3,col=1001)}
if(grenz[i] < 1 && grenz[i+2] >= 1){
text(475,550+.2*(i+1),"1.000",cex=.7,col=1)
lines(c(455,458),c(550+.2*(i+1),550+.2*(i+1)),lwd=.3,col=1001)}
if(grenz[i] < 1.6 && grenz[i+2] >= 1.6){
text(475,550+.2*(i+1),"1.600",cex=.7,col=1)
lines(c(455,458),c(550+.2*(i+1),550+.2*(i+1)),lwd=.3,col=1001)}
if(grenz[i] < 2.5 && grenz[i+2] >= 2.5){
text(475,550+.2*(i+1),"2.500",cex=.7,col=1)
lines(c(455,458),c(550+.2*(i+1),550+.2*(i+1)),lwd=.3,col=1001)}
if(grenz[i] < 4 && grenz[i+2] >= 4){
text(475,550+.2*(i+1),"4.000",cex=.7,col=1)
lines(c(455,458),c(550+.2*(i+1),550+.2*(i+1)),lwd=.3,col=1001)}
}
polygon(c(440,440,455,455),c(550,750,750,550),density=0,col=1,lwd=.3)


###############################################################
### code chunk number 55: ModellwahlBayesfaktorenHWodernicht
###############################################################
  x <- c(233, 385, 129)
  alpha <- 1
  beta <- 1
  # Forderungen an q
  lower <- function(q){
    lower <- pmax(-q^2, -(1-q)^2)
    return(lower)
  }

  upper <- function(q){
    upper <- q*(1-q)
    return(upper)
  }
  # Wkeiten aus q und d berechnen
  myprob <- function(q, d){
    p1 <- q^2+d
    p2 <- 2*q*(1-q)-2*d
    p3 <- (1-q)^2+d
    p <- c(p1,p2,p3)
    if(all(p>=0) & all(p<=1))
    return(p)
    else
    return(rep(NA, 3))
  }

  n = sum(x)
  # Marginale Likelihood im HW-Modell:
  hw.marglik = exp(
  lfactorial(n)-sum(lfactorial(x))
  + x[2]*log(2) - lbeta(alpha, beta)
  + lbeta(alpha + 2*x[1] + x[2], beta + x[2] + 2*x[3])
  )
  # Marginale Likelihood im HWU-Modell:

  # über numerische Integration:
  # bei Priori-Gleichverteilung von q und d gegeben q gleich auf eingeschränktem Bereich ergibt sich

  dintegral = function(q){
    erg = integrate(function(d){
      sapply(d, function(one.d){
        dmultinom(x, prob = myprob(q, one.d))
      })
    }, lower = lower(q), upper = upper(q))
    erg$value
  }
  erg.gesamt = integrate(function(q){
    sapply(q, function(one.q){
      exp(-log(upper(one.q)-lower(one.q))+log(dintegral(one.q)))
    })
  }, lower = 0, upper = 1)
  hwu.marglik = erg.gesamt$value

  # Bayes-Faktor ist also:
  bayesfactor = hw.marglik / hwu.marglik


###################################################
### code chunk number 56: darmkrebsMargLik1
###################################################
  ## Berechnungen für das Beispiel darmkrebsMargLik1
  # Daten laden
  Zbeob = colonCancer[-1] # ohne Z0
  n = sum(Zbeob) # 196
  N = length(Zbeob) # 6
  k = 1:N
  # beobachtete Loglikelihood
  darm.beob.loglikelihood = function(p){
    one.p = function(p){
      sum(Zbeob*dbinom(k, N, p, log = TRUE)) - n*log(1-(1-p)^N)
    }
    sapply(p, one.p)
  }
  # Priori-Samples
  M = 3*1e4
  wh = 5
  margLogLik1.matrix = lik.matrix = matrix(exp(darm.beob.loglikelihood(rbeta(wh*M, .5, .5))), ncol = wh)
  # sequentielle Mittelwerte
  for(i in 1:wh){
    margLogLik1.matrix[,i] = log(cumsum(lik.matrix[,i]) / (1:M))
  }
  margLogLik1.final = margLogLik1.matrix[M,]
#  Mformatted <- format (trunc(M), big.mark = ",\\\\\\\\")
Mformatted <- formatBig(M)


###################################################
### code chunk number 57: MCmargLik1Plot
###################################################
# Remark: Due to randomness, this figure may look a bit different
# than Figure 8.12a) in the book
getOption("SweaveHooks")[["fig"]]()
      marg.lim = -c(441, 440)
      matplot(margLogLik1.matrix, type = "l",
      col = ifelse(rep(getOption ("myColour"), wh), rainbow(wh), grey(seq(0,0.8,length = wh))),
      lty = 1,
      xlab = math (M), ylab = math (plain(log)~"{"~hat(f)(x)~"}"),
      ylim = marg.lim
      )


#################################################################
### code chunk number 58: darmkrebsVergleichJeffreyLikelihood
#################################################################
getOption("SweaveHooks")[["fig"]]()
      eps = 1e-3
      pvec = seq(eps,1-eps, length = 2500)
      prior = dbeta(pvec, .5, .5)
      prior = prior/max(prior)
      likelihood = exp(darm.beob.loglikelihood(pvec))
      likelihood = likelihood/max(likelihood)
      plot(pvec, prior, type = "l", xlab = math (pi), ylab = "", yaxt = "n",
           ylim = c(min(min(likelihood),min(prior)),1))
      lines(pvec, likelihood, lty=5)


################################################################
### code chunk number 59: darmkrebsMarginaleLikelihoodSamples
################################################################
## wegen großer Variabilität der Schätzungen:
  set.seed(3)

  M = 1e5
  kfull = 0:6
  wh = 5
  ## Ketten neu erzeugen
  pisamples.matrix = matrix(nrow = M, ncol = wh)
  for(j in 1:wh) {
    Z0 <- 10
    Z <- c(Z0, Zbeob)
    pi <- 0.5

    nburnin <- 100
    niter <- M + nburnin
    pisamples <- Z0samples <- rep(NA, niter)

    for(i in 1:niter){
      pi <- rbeta(1, 0.5+sum(kfull*Z), 0.5+sum((6-kfull)*Z))
      Z[1] <- rnbinom(1, size=n, prob=1-(1-pi)^6)
      pisamples[i] <- pi
      Z0samples[i] <- Z[1]
    }

    pisamples <- pisamples[-c(1:nburnin)]
    Z0samples <- Z0samples[-c(1:nburnin)]
    pisamples.matrix[,j] = pisamples
  }

  ## auf Log-Ebene rechnen um genauer zu sein
  margLogLik2.matrix = apply(pisamples.matrix, 2, function(one.pisamples){
    invers.lik = 1/exp(darm.beob.loglikelihood(one.pisamples))
    log(1:M)-log(cumsum(invers.lik))}
  )
  margLogLik2.final = margLogLik2.matrix[M,]


############################################################
### code chunk number 60: darmkrebsMarginaleLikelihoodPlot
############################################################
getOption("SweaveHooks")[["fig"]]()
matplot(margLogLik2.matrix, 
        type = "l",
        col = 
        ifelse(rep(getOption ("myColour"), wh), 
               rainbow(wh), grey(seq(0,0.8,length = wh))),
        lty = 1,
        xlab = math (M), 
        ylab = math (plain (log)~hat(f)(x)), 
        ylim = c(-439, -437.5),
        xaxt="n")

## plot the x axis labels separately for beauty
xticks <- c(0, 2, 4, 6, 8, 10) * 10^4
axis(side=1, at=xticks,
     labels=
     expression(0, 2 %.% 10^4, 4 %.% 10^4, 6 %.% 10^4, 8 %.% 10^4, 10^5))


#####################################################################################
### code chunk number 61: darmkrebsMarginaleLikelihoodBerechnungenFuerDrittenAnsatz
#####################################################################################
# Berechnungen für zweite Version:
post.mean.pi <- round(mean(pisamples),3)
post.pi.atmean.est <- mean(sapply(Z0samples, 
                                  function(one.Z0){
                                      Z[1] = one.Z0
                                      dbeta(post.mean.pi, 0.5+sum(kfull*Z), 
                                            0.5+sum((6-kfull)*Z))}))
likelihood.at.postmean <- exp(darm.beob.loglikelihood(post.mean.pi))
priori.at.postmean <- dbeta(post.mean.pi, .5, .5)
log.marglik2 <- log(likelihood.at.postmean) + log(priori.at.postmean) - 
    log(post.pi.atmean.est)


