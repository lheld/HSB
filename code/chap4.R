#################################################################
### R code for Chapter 4
#################################################################

### Encoding: UTF-8

###################################################
### code chunk number 1
###################################################
getOption("SweaveHooks")[["fig"]]()
T1 <- function(xq, eq, n, lambda0){
    return(sqrt(n)*(xq - eq*lambda0)/sqrt(eq*lambda0))
}
T2 <- function(xq, eq, n, lambda0){
    return(sqrt(n)*(xq - eq*lambda0)/sqrt(xq))
}
T3 <- function(xq, eq, n, lambda0){

    result <- n*xq*(log(xq/eq)-log(lambda0)-1)+lambda0*n*eq
    result <- ifelse(is.na(result), lambda0*n*eq, result)
    result <- result*2
    result <- sqrt(result)*sign(xq/eq-lambda0)
    return(result)
}

scotland <- scotlandData

x <- scotland[,1]
e <- scotland[,2]

t1 <- rep(NA, length(x))
t2 <- rep(NA, length(x))
t3 <- rep(NA, length(x))

for(i in 1:length(x)){
    t1[i] <- T1(x[i], e[i], 1, 1)
    t2[i] <- T2(x[i], e[i], 1, 1)
    t3[i] <- T3(x[i], e[i], 1, 1)
}
mint <- min(t1, t2[abs(t2)!=Inf])
maxt <- max(t1, t2[abs(t2)!=Inf])
range <- c(mint, maxt)

plot(t1, t2, na.rm=T, xlim=range, ylim=range,
     xlab = math (T[1]), ylab = math (T[2])
     )
lines(range, range, lty=2)

normalQuantiles <- c (-1, 1) * qnorm (0.975)
abline (v = normalQuantiles, h = normalQuantiles, col = "gray")

select <- t1[abs(t2)==Inf]


###################################################
### code chunk number 2: score-test-poisson-hist-t1
###################################################
getOption("SweaveHooks")[["fig"]]()
p1 <- 2 * pnorm(abs(t1), lower.tail=FALSE)
truehist(p1, 
         breaks=c(0:10)/10, 
         ylim=c(0,5),
         xlab=math(p),
         prob=TRUE,
         col="gray")
abline(h=1,
       lty=2)


###################################################
### code chunk number 3: score-test-poisson-hist-t2
###################################################
getOption("SweaveHooks")[["fig"]]()
p2 <- 2 * pnorm(abs(t2), lower.tail=FALSE)
truehist(p2, 
         breaks=c(0:10)/10, 
         ylim=c(0,5),
         xlab=math(p),
         prob=TRUE,
         col="gray")
abline(h=1,
       lty=2)


###################################################
### code chunk number 4: scoreKIPoisson
###################################################
getOption("SweaveHooks")[["fig"]]()
ci1 <- function(xq, eq, n, level=0.95){
    q <- qnorm(1-((1-level)/2))
    center <- (xq+q^2/(2*n))/eq
    term <- q/(2*eq)*sqrt(4*xq/n+q^2/n^2)
    lower <- center - term
    upper <- center + term
    signif <- (lower>1)|(upper<1)
    return(c(center, lower, upper, signif))
}
ci <- matrix(NA, ncol=4, nrow=length(x))
for(i in 1:length(x))
    ci[i,] <- ci1(x[i], e[i], n=1)

signif <- ci[,4]
par (font.lab = 1)

scotOrder <- order(ci[,1])
save(scotOrder, 
     file="../data/scotland/scotOrder.RData")

gplots::plotCI(x=ci[scotOrder,1], li=ci[scotOrder,2], ui=ci[scotOrder,3],
               err="y", ylab="Relative risk", xlab="Region",
               col="black",
               pch=19,
               gap=0.25,
               ylim=c(-0.5, 13))
abline(h=1, lty=2)


###################################################
### code chunk number 5
###################################################
x1 = 233
x2 = 385
x3 = 129
n <- x1+x2+x3
qml <- (2*x1+x2)/(2*n)
se.qml <- sqrt(qml*(1-qml)/(2*n))


###################################################
### code chunk number 6
###################################################
x1 <- 11
e1 <- 3.04
ml <- x1/e1
se.ml <- sqrt(ml/e1)
lower <- ml - 1.96*se.ml
upper <- ml + 1.96*se.ml

se2.ml <- sqrt(1/(4*e1))
lower2 <- sqrt(ml) - 1.96*se2.ml
upper2 <- sqrt(ml) + 1.96*se2.ml

lower3 <- lower2^2
upper3 <- upper2^2



###################################################
### code chunk number 7: varstabWaldKIPoisson
###################################################
getOption("SweaveHooks")[["fig"]]()
ci2 <- function(xq, eq, n, level=0.95)  # approximate ci based on varstab transform
{
    q <- qnorm(1-((1-level)/2))
    center <- sqrt (xq / eq)
    term <- q / sqrt (4 * n * eq)
    lower <- center - term
    upper <- center + term
    signif <- (lower>1)|(upper<1)
    return(c(center^2, lower^2, upper^2, signif))
}
ci <- matrix(NA, ncol=4, nrow=length(x))
for(i in 1:length(x))
    ci[i,] <- ci2(x[i], e[i], n=1)

signif2 <- ci[,4]
par (font.lab = 1)
gplots::plotCI(x=ci[scotOrder,1], li=ci[scotOrder,2], ui=ci[scotOrder,3],
               err="y", ylab="Relative risk", xlab="Region",
               col="black",
               pch=19,
               gap=0.25,
               ylim=c(-0.5, 13))
abline(h=1, lty=2)


###################################################
### code chunk number 8
###################################################
  h <- function(pi)
  asin(sqrt(pi))
  hinv <- function(x)
    sin(x)^2
  n <- 100
  x <- 2
  pid <- x/n
  hpid <- h(pid)
  hlow <- hpid - 1.96*sqrt(1/(4*n))
  hupp <- hpid + 1.96*sqrt(1/(4*n))
  low <- hinv(hlow)
  upp <- hinv(hupp)


###################################################
### code chunk number 9: FisherszTrafo
###################################################
x <- iten$BAK * 2100
y <- iten$L_500a_1
n <- length(x)
xq <- mean(x)
yq <- mean(y)
xycor <- cor(x, y)
sexycor <- (1 - xycor^2) / sqrt(n)
corxyci <- xycor + c(-1, 1) * 1.96 * sexycor
zxy <- atanh(xycor)
zxyci <- zxy + c(-1, 1) * 1.96 / sqrt(n)
corxyci2 <- tanh(zxyci)


###################################################
### code chunk number 10
###################################################
getOption("SweaveHooks")[["fig"]]()
  mint <- min(t1, t3)
  maxt <- max(t1, t3)
  range <- c(mint, maxt)
  plot(t1, t3, na.rm=T, xlim=range, ylim=range,
  xlab = math (T[1]), ylab = math (T[4]))
  lines(range, range, lty=2)
abline (v = normalQuantiles, h = normalQuantiles, col = "gray")


#######################################################
### code chunk number 11: score-test-poisson-hist-t1
#######################################################
getOption("SweaveHooks")[["fig"]]()
p1 <- 2 * pnorm(abs(t3), lower.tail=FALSE)
truehist(p1, 
         breaks=c(0:10)/10, 
         ylim=c(0,5),
         xlab=math(p),
         prob=TRUE,
         col="gray")
abline(h=1,
       lty=2)


###################################################
### code chunk number 12: LQKIPoissonSchottland
###################################################
## define a general function which computes likelihood confidence intervals
likelihood.ci <- function(gamma, ## the confidence level
                          loglik, ## the log-likelihood function
                          theta.hat, ## the MLE
                          lower, ## lower bound of parameter space
                          upper, ## upper bound of parameter space
                          comp.lower = TRUE, ## compute lower bound of CI?
                          comp.upper = TRUE, ## compute upper bound of CI?
                          ...) ## additional arguments for the log-likelihood 
{
    ## target function, such that f(theta)=0 gives CI limits
    f <- function(theta, ...)
    {
        loglik(theta, ...) - loglik(theta.hat, ...) + 1/2 * qchisq(gamma, df=1)
    }
    
    ## compute lower and upper bounds of CI
    ret <- c()
    if(comp.lower)
    {
        hl.lower <- uniroot(f,
                            interval = c(lower, theta.hat),
                            ...)$root
        ret <- c(ret, hl.lower)
    }
    if(comp.upper) 
    {
        hl.upper <- uniroot(f,
                            interval = c(theta.hat, upper),
                            ...)$root
        ret <- c(ret, hl.upper)
    }
    return(ret)
}

## the log-likelihood of lambda in the Poisson model
loglik.poisson <- function(lambda, x, e, log = TRUE)
{
    dpois(x = x, lambda = lambda * e, log = log)
}

## get the data
x <- scotlandData[, 1]
e <- scotlandData[, 2]

## here we will save the bounds of the CIs
likCI <- matrix(nrow=length(x), ncol=2)

## confidence level
confLevel <- 0.95 

## small positive value
eps <- sqrt(.Machine$double.eps)

## now process all regions
for(i in seq_along(x)) 
{
    res <- likelihood.ci(gamma = confLevel,
                         loglik = loglik.poisson, 
                         theta.hat = x[i] / e[i],
                         lower = eps, upper = 1/eps, 
                         x = x[i], e = e[i], 
                         comp.lower = x[i] > 0)
    if (x[i] == 0)
    {
        res <- c (0, res)
    }
    
    likCI[i, ] <- res
}


########################################################
### code chunk number 13: LQKIPoissonSchottlandFigure
########################################################
getOption("SweaveHooks")[["fig"]]()
plotCI (x=(x/e)[scotOrder], 
        li=likCI[scotOrder, 1], ui=likCI[scotOrder, 2],
        err="y", ylab="Relative risk", xlab="Region",
        col="black",
        pch=19,
        gap=0.25, ylim=c(-0.5, 13))
        
abline(h=1, lty=2)


###################################################
### code chunk number 14: LQWaldKI
###################################################
getOption("SweaveHooks")[["fig"]]()
  beob <- 11
  erw <- 3.04
  theta.ml = beob/erw
  theta.vector = seq(from = 1.5, to = 7, length = 200)
  rel.loglike = function(theta){
    return( sum(dpois(beob, theta*erw, log = TRUE)) - sum(dpois(beob, theta.ml*erw, log = TRUE)))
  }
  like.werte = sapply(theta.vector, rel.loglike)
  matplot(theta.vector, like.werte, xlab = math (theta), ylab = math (tilde(l)(theta)), type = "l", yaxt = "n")
  drawml(theta.ml, 0)
  quadrApprox = function(theta){
    return(-0.5*beob/(theta.ml^2)*(theta - theta.ml)^2)
  }
  approx.werte = quadrApprox(theta.vector)
  lines(theta.vector, approx.werte, lty = 2)
  c = -0.5*qchisq(1-0.05,1)
  axis(side = 2, at = c(0, c), labels = c(0, "c"))
  abline(h = c, lty = 2, col = "gray")
  wald.unten = theta.vector[min(which(approx.werte>=c))]
  wald.oben = theta.vector[max(which(approx.werte>=c))]
  like.unten = theta.vector[min(which(like.werte>=c))]
  like.oben = theta.vector[max(which(like.werte>=c))]
  ki.klammer = function(unten, oben, height, label = "KI", abstand = 1, ...){
    segments(unten, height, unten, height + abstand, ...) # links
    segments(oben, height, oben, height + abstand, ...) # rechts
    segments(unten, height + abstand, oben, height + abstand, ...) # horizontal
    text(x = mean(c(unten, oben)), y = height, labels = label, pos = ifelse(abstand > 0, 3, 1), ...)
  }
  if (getOption ("myColour")){
    farben = c(2,3)
  }else{
    farben = c("black","darkgrey")
  }
  ki.klammer(wald.unten, wald.oben, c, "Wald confidence interval", abstand = -0.2, col = farben[1] )
  ki.klammer(like.unten, like.oben, c, "Likelihood ratio confidence interval", abstand = +0.2, col = farben[1])


###################################################
### code chunk number 15: anteilZahlenBerechnung
###################################################
## Confidence Intervals for Binomial Experiment ********************

## Setting
data <- data.frame (
                    n = c (rep (10, 3), rep (100, 3)),
                    x = c (0, 1, 5, 0, 10, 50),
                    x.wald = c (0.5, 1, 5, 0.5, 10, 50)
                    )
prob = 0.95 # confidence level?
q <- qnorm(1-(1-prob)/2)
q.onesided = qnorm(prob)
q.wald = q.onesided*(data$x==0)+q*(data$x!=0)

## Wald
wald.se <- with (data, sqrt (x.wald/n * (1 - x.wald/n) / n))
wald.lower <- with (data, x/n - q.wald * wald.se)
wald.lower <- ifelse (wald.lower <= 0, 0, wald.lower)
wald.upper <- with (data, x/n + q.wald * wald.se)
wald.lower <- ifelse(data$x==0, 0, wald.lower)

## Wald for logit
wald2.se = with (data, sqrt(1/x.wald + 1/(n-x.wald)))
wald2.lower = with (data, plogis (log(x.wald/(n-x.wald)) - q.wald * wald2.se))
wald2.upper = with (data, plogis (log(x.wald/(n-x.wald)) + q.wald * wald2.se))
wald2.lower <- ifelse(data$x==0, 0, wald2.lower)


## varianzstabilisiertes Wald:
varstabWald.se <- with (data, sqrt (1 / (4*n)))
varstabWald.lower <- with (data, asin (sqrt (x / n)) - q.wald * varstabWald.se)
varstabWald.lower <- ifelse (varstabWald.lower <= 0, 0, sin (varstabWald.lower)^2)
varstabWald.upper <- with (data, asin (sqrt (x / n)) + q.wald * varstabWald.se)
varstabWald.upper <- ifelse (varstabWald.upper >= pi / 2, 1, sin (varstabWald.upper)^2)

## Wilson
pseudo.est = with (data, ( x + q^2 / 2 ) / ( n + q^2 ))
pseudo.se = with (data, sqrt( x/n*(1 - x/n) / n + q^2 / n /(4*n) ) / (1 + q^2 / n))
wilson.lower = pseudo.est - q * pseudo.se
wilson.upper = pseudo.est + q * pseudo.se

## Likelihood
loglik.binom <- function(prob, x, size, log = TRUE) # passendes Format for CI-Funktion
  dbinom (x = x, size = size, prob = prob, log = log)

eps <- sqrt (.Machine$double.eps)

keinTreffer <- with (data, x == 0)
nurTreffer <- with (data, x == n)
lik.lower <- lik.upper <- double (nrow (data))

for (i in 1:nrow (data)) {
  if (nurTreffer[i]){
    res <- with (data,
                 likelihood.ci(gamma = prob, loglik.binom, theta.hat = 1 - eps,
                               lower = eps, upper = 1 - eps, x = x[i], size = n[i], comp.upper = FALSE)
                 )
    res <- c (res, 1)
  }
  else if (keinTreffer[i]){
    res <- with (data,
                 likelihood.ci(gamma = prob, loglik.binom, theta.hat = eps,
                               lower = eps, upper = 1 - eps, x = x[i], size = n[i], comp.lower = FALSE)
                 )
    res <- c (0, res)
  }
  else{
    res <- with (data,
                 likelihood.ci (gamma = prob, loglik.binom, theta.hat = x[i] / n[i],
                                lower = eps, upper = 1 - eps, x = x[i], size = n[i])
                 )
  }
  lik.lower[i] <- res[1]
  lik.upper[i] <- res[2]
}

## Clopper-Pearson: in settings.R definierte Funktion benutzen
clopperP <- with (data, clopperPearson (x, n, prob))

## prepare table
bounds <- cbind (wald.lower, wald.upper, wald2.lower, wald2.upper, varstabWald.lower, varstabWald.upper, wilson.lower, wilson.upper, lik.lower, lik.upper, clopperP)
bounds <- format (round (bounds, 3), digits = 2)

wald <- paste (bounds[,1], " to ", bounds[,2], sep="")
wald2 <- paste (bounds[,3], " to ", bounds[,4], sep="")
varstabWald <- paste (bounds[,5], " to ", bounds[,6], sep="")
wilson <- paste (bounds[,7], " to ", bounds[,8], sep="")
lik <- paste (bounds[,9], " to ", bounds[,10], sep="")
clopperP <- paste (bounds[,11], " to ", bounds[,12], sep="")

data <- data.frame (n = data$n, x = data$x, 
               wald = wald, wald2 = wald2, varstabWald = varstabWald,
               wilson = wilson, lik = lik, clopperP = clopperP)


names(data) <- c("$\\boldsymbol{n}$", "$\\boldsymbol{x}$",
                 "1) Wald for $\\boldsymbol{\\pi}$", "2) Wald for $\\boldsymbol{\\logit (\\pi)}$",
                 "3) Wald for $\\boldsymbol{\\arcsin (\\sqrt{\\pi})}$",
                 "4) Wilson", "5) Likelihood", "6) Clopper--Pearson"
                 )


###################################################
### code chunk number 16: anteilZahlenA
###################################################
#in Latex konvertieren
w <-
latex(data[, c (1, 2, 3, 4, 5)],                             # was konvertieren?
      file = latexTempFile,             # um direkt auf Ausgabe zu schreiben
      table.env = FALSE,
      n.cgroup = c (2, 3),              # wieviele Spalten jeweils zusammenfassen?
      rowname = NULL,                  # um keine Zeilennamen zu haben
      booktabs = TRUE,                  # um booktabs zu benutzen
      colnamesTexCmd = "bfseries",
      col.just = Cs(r, r, C, C, C,),  # Spalten ausrichten
      center = "none",
      multicol = FALSE,
      numeric.dollar = FALSE
     )
postLatex (w, widthFactor = 1.0, minipage = FALSE)


###################################################
### code chunk number 17: anteilZahlenB
###################################################
w <-
latex(data[, c (1, 2, 6, 7, 8)],                             # was konvertieren?
      file = latexTempFile,             # um direkt auf Ausgabe zu schreiben
      table.env = FALSE,
      n.cgroup = c (2, 3),              # wieviele Spalten jeweils zusammenfassen?
      rowname = NULL,                  # um keine Zeilennamen zu haben
      booktabs = TRUE,                  # um booktabs zu benutzen
      colnamesTexCmd = "bfseries",
      col.just = Cs(r, r, C, C, C,),  # Spalten ausrichten
      center = "none",
      multicol = FALSE,
      numeric.dollar = FALSE
     )
postLatex (w, widthFactor = 1.0, minipage = FALSE)


###################################################
### code chunk number 18: anteilKIvergleich_Wald
###################################################
getOption("SweaveHooks")[["fig"]]()
    # Szenario festlegen:
    n1 = 50
    x1 = 0:n1
    x1.wald <- x1
    x1.wald[1] <- 0.5
    x1.wald[n1+1] <- 49.5


    prob = 0.95 # confidence level?
    q = qnorm(1-(1-prob)/2)
    ygrenzen = c(0.8, 1) # y-Skala?
    q.onesided = qnorm(prob)
    q.wald =  c(q.onesided, rep(q, 49), q.onesided)


    # Wald-Intervalle berechnen
    se = sqrt(x1.wald/n1*(1-x1.wald/n1)/n1)
    wald.lower = x1/n1 - q.wald*se
    wald.upper = x1/n1 + q.wald*se
    wald.lower <- ifelse (wald.lower <= 0, 0, wald.lower)
    wald.upper <- ifelse (wald.upper >= 1, 1, wald.upper)


    plot.coverage(cbind(wald.lower, wald.upper), ylim = ygrenzen, smooth = 0.05)

    wald.widths = wald.upper - wald.lower

    

#######################################################
### code chunk number 19: anteilKIvergleich_logitWald
#######################################################
getOption("SweaveHooks")[["fig"]]()
    # Wald-Intervalle for logit(pi) berechnen und dann rücktransformieren
    se.logit = sqrt(1/x1.wald + 1/(n1-x1.wald))
    wald2.lower = plogis( log(x1.wald/(n1-x1.wald)) - q.wald*se.logit )
    wald2.upper = plogis( log(x1.wald/(n1-x1.wald)) + q.wald*se.logit )
    # Nach unserer Def. geht Intervall in beiden Extremfällen von 0 bis 1
    wald2.lower[c(1)] = 0
    wald2.upper[c(n1 + 1)] = 1
    # plotten

    wald2.widths = wald2.upper - wald2.lower
    wald2.widths[c(1, n1 + 1)] = wald2.widths[c(2, n1 + 1)]
    wald2.lower[c(1, n1 + 1)] = 0
    wald2.upper[c(1, n1 + 1)] = wald2.widths[c(1, n1 + 1)]

    plot.coverage(cbind(wald2.lower, wald2.upper), ylim = ygrenzen, smooth = 0.05)



#########################################################
### code chunk number 20: anteilKIvergleich_varstabWald
#########################################################
getOption("SweaveHooks")[["fig"]]()
    ## varianzstabilisiertes Wald:
    varstabWald.se <- sqrt (1 / (4*n1))
    varstabWald.lower <- asin (sqrt (x1 / n1)) - q.wald * varstabWald.se
    varstabWald.upper <- asin (sqrt (x1 / n1)) + q.wald * varstabWald.se
    varstabWald.lower <- ifelse (varstabWald.lower <= 0, 0, sin (varstabWald.lower)^2)
    varstabWald.upper <- ifelse (varstabWald.upper >= pi / 2, 1, sin (varstabWald.upper)^2)
    varstabWald.lower[c(1)] = 0
    varstabWald.upper[c(n1 + 1)] = 1

    ## plotten
    plot.coverage(cbind(varstabWald.lower, varstabWald.upper), ylim = ygrenzen, smooth = 0.05)

    varstabWald.widths = varstabWald.upper - varstabWald.lower


###################################################
### code chunk number 21: anteilKIvergleich_Wilson
###################################################
getOption("SweaveHooks")[["fig"]]()
    # Wilson-Intervalle for pi berechnen
    pseudo.est = ( x1 + q^2 / 2 ) / ( n1 + q^2 )
    pseudo.se = sqrt( x1/n1*(1 - x1/n1) / n1 + q^2 / n1 /(4*n1) ) / (1 + q^2 / n1)
    wilson.lower = pseudo.est - q*pseudo.se
    wilson.upper = pseudo.est + q*pseudo.se
    # plotten
    plot.coverage(cbind(wilson.lower, wilson.upper), ylim = ygrenzen, smooth = 0.05)
    wilson.widths = wilson.upper - wilson.lower


########################################################
### code chunk number 22: anteilKIvergleich_Likelihood
########################################################
getOption("SweaveHooks")[["fig"]]()
    # Log-Likelihood for Binomialfall (ohne Konstante)
    loglik = function(pi, x, n){
      return( x*log(pi) + (n-x)*log(1-pi) )
    }

    # außer an Extremstellen so berechnen:
    eps <- 1e-12
    ci.like.middle <- sapply(x1[-c(1, n1 + 1)], function(one){
      likelihood.ci(gamma=prob, loglik, theta.hat=one / n1, lower=eps, upper=1-eps, x = one, n = n1)}
    )
    ci.like.middle = t(ci.like.middle)

    # an Extremstellen jeweils "mittlere" Grenze berechnen:
    # Fall x1 = 0
    ci.like.0 = c(0, likelihood.ci(gamma=prob, loglik, theta.hat = eps, lower=eps, upper=1-eps, x = 0, n = n1, comp.lower = FALSE))
    # Fall x1 = 1
    ci.like.1 = c(likelihood.ci(gamma=prob, loglik, theta.hat = 1-eps, lower=eps, upper=1-eps, x = n1, n = n1, comp.upper = FALSE), 1)

    ci.like = rbind(ci.like.0, ci.like.middle, ci.like.1)
    # plotten
    plot.coverage(ci.like, ylim = ygrenzen, smooth = 0.05)
    like.widths = ci.like[,2] - ci.like[,1]


############################################################
### code chunk number 23: anteilKIVergleich_ClopperPearson
############################################################
getOption("SweaveHooks")[["fig"]]()
    ci.clopperP <- clopperPearson (x1, n1, prob)
    plot.coverage (ci.clopperP, ylim = ygrenzen, smooth = 0.05)
    clopperP.widths <- ci.clopperP[, "upper"] - ci.clopperP[, "lower"]


###################################################
### code chunk number 24: anteilKIwidths
###################################################
getOption("SweaveHooks")[["fig"]]()
  breiten = (cbind(wald.widths, wald2.widths, varstabWald.widths,
                   wilson.widths, like.widths, clopperP.widths
                   )[1:(n1/2+1),]
             )
  par (font.lab = 1)
  matplot(x1[1:(n1/2+1)], breiten, xlab = math (x, font = 3), ylab = "Width", lty = c(1,2,4,1,2,3),
  col = c(rep(c("black", "gray"),c(3,2)), "black"), type = "l", ylim = c(0, max(breiten[breiten!=1])+0.01),
  xlim = c(0,n1/2), lwd = 2)
  legend("bottomright", lty = c(4, 3, 2, 1, 2, 1), col = c("black", "black", "black", "gray", "gray", "black"),
         legend = c("variance stab. Wald", "Clopper-Pearson","Wald via logit",
         "Wilson", "Likelihood", "Wald"), lwd = 2, bty = "n")


