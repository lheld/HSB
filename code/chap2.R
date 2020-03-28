#################################################################
### R code for Chapter 2
#################################################################

### Encoding: UTF-8

###################################################
### code chunk number 1: binom1
###################################################
getOption("SweaveHooks")[["fig"]]()
vec=seq(from=0, to=1, by=0.01)
values=dbinom(x=2, size=10, prob=vec)
matplot(vec, values, xlab=math(pi), ylab=math(L(pi)), type="l")
drawml(0.2, dbinom(x=2, size=10, prob=0.2))


###################################################
### code chunk number 2: binom2
###################################################
getOption("SweaveHooks")[["fig"]]()
values=dbinom(x=0, size=10, prob=vec)
matplot(vec, values, xlab=math(pi), ylab=math(L(pi)), type="l")
drawml(0.0, dbinom(0, 10, 0.0), dezi=1)


###################################################
### code chunk number 3: dhyper1
###################################################
M <- 13
n <- 10
x <- 5
ml <- c(25, 26)
(dhyper(x=x, m=M, n=ml-M, k=n))


###################################################
### code chunk number 4: capturerecapture
###################################################
getOption("SweaveHooks")[["fig"]]()
    plot.caprecap=function(M=26, n=63, y=5, right=1000){
      vec=seq(from=max(n, M+n-y), to=right, by=1)
      values=dhyper(x=y, m=M, n=vec-M, k=n)
      matplot(vec, values, xlab=math (N), ylab=math(L(N)), type="l")
      mlindex=which.max(values)
      drawml(vec[mlindex], max(values), down = TRUE)
    }
    plot.caprecap(right=1500)


###################################################
### code chunk number 5: crctable
###################################################
M <- c(26, 25, 25, 13)
n <- c(63, 30, 30, 10)
x <- c(5, 11, 10, 5)
ml <- trunc(M*n/x)
lik <- dhyper(x,M,ml-M,n)


###################################################
### code chunk number 6: exponential-lambda
###################################################
getOption("SweaveHooks")[["fig"]]()
vec=seq(from=0.0001, to=0.002, by=0.000005)
values=(vec^uncN)*exp(-vec*uncSum)
faktor <- ceiling (abs (log (max (values), 10)))
matplot(vec, values*(10^faktor), xlab=math(lambda),
        ylab=substitute(italic(L)(lambda) %*% 10^faktor, list (faktor = faktor)), type="l")
drawml(1/uncMean, (1/uncMean)^uncN*exp(-1/uncMean*uncSum)*10^faktor, digi=2, down=T)


###################################################
### code chunk number 7: exponential-mu
###################################################
getOption("SweaveHooks")[["fig"]]()
matplot(vec, values*(10^faktor), xlab=math(mu),
        ylab=substitute(italic(L)(mu)%*% 10^faktor, list (faktor = faktor)), 
        type="l",
        xaxt="n")
tickPoints <- c(0.0005, 0.001, 0.0015, 0.002)
axis(side=1, at=tickPoints, labels=formatRound(1/tickPoints))

x <- 1/uncMean
y <- (1/uncMean)^uncN*exp(-1/uncMean*uncSum)*10^faktor
digi <- 3
dezi <- 0
segments(x, par("usr")[3], x, y, lty="dashed", col = "black")
axis(side=1, at=x, font=2, labels=format(1/x, nsmall=dezi, digits=digi), line=1, tick=F)


###################################################
### code chunk number 8: weibull
###################################################
weibullML <- optim (c (2, 1000),
                    function (amu) exp (sum (dweibull (unc, shape = amu[1], scale = amu[2], log = TRUE))),
                    control = list (fnscale = -1e-171), method = "L-BFGS-B", lower = rep (1e-10, 2), hessian = TRUE
                    )
## important: correct scaling of the function!!
veca=seq(from=0.8, to=1.6, by=0.01)
vecmu=seq(from=800, to=1700, by=3)
values=matrix(nrow=length(veca), ncol=length(vecmu))
p=prod(unc)
for(i in 1:length(veca)){
 for(j in 1:length(vecmu)){
   a=veca[i]
   mu=vecmu[j]
   values[i,j]=exp (sum (dweibull (unc, a, mu, log = TRUE)))
 }
}


###################################################
### code chunk number 9: weibullPlot
###################################################
getOption("SweaveHooks")[["fig"]]()
contour(x=veca, y=vecmu, z=values*10^faktor, nlevels=20, xlab=math(alpha), ylab=math(mu),
        xaxs="i", yaxs="i", labcex=1)
abline (v = 1, lty = 2, col = "black")


###################################################
### code chunk number 10: gamma
###################################################
veca=seq(from=0.7, to=1.9, by=0.005)
vecb=seq(from=0.0005, to=0.0018, by=0.00002)
values=matrix(data = 0, nrow=length(veca), ncol=length(vecb))
## für Berechnung der Berechnungsdiagonale Steigung der beiden Geraden un einen Punkt unten und oben angeben,
## die Geradenpunkte sein sollen
steigung <- 0.0002 / 0.2
pu <- c (1.6, 0.001)
po <- c (1, 0.0014)
abu <- pu[2] - pu[1] * steigung
abo <- po[2] - po[1] * steigung

for(i in 1:length(veca)){
  a = veca[i]
  # Iterationen sparen, nur "dickere Diagonale" berechnen
  for(j in which(vecb >= (abu + steigung*a) & vecb <= (abo + steigung*a))){
    a=veca[i]
    b=vecb[j]
    values[i,j] <- exp(sum(dgamma(unc, shape=a, rate=b, log = TRUE)))
  }
}


###################################################
### code chunk number 11: gammaPlot
###################################################
getOption("SweaveHooks")[["fig"]]()
contour(x=veca, y=vecb, z=values*10^faktor, nlevels=15, xlab=math(alpha), ylab=math(beta),
        xaxs="i", yaxs="i", labcex=1)
abline (v = 1, lty = 2, col = "black")


###################################################
### code chunk number 12: gamma2
###################################################

## MLE bestimmen, um Grafik zentrieren zu können
gammaML <- optim (c (uncMean, uncSd^2 / uncMean),
                  function (muphi){
                    - sum(dgamma(unc, shape=muphi[1]/muphi[2], scale=muphi[2], log = TRUE)) # -loglik
                  }
                  )
gammaML <- gammaML$par
vecmu=seq(from=800, to=1600, length = 200)
vecphi=seq(from=500, to=1700, length = 200)
values=matrix(0, nrow=length(vecmu), ncol=length(vecphi)) 

for(i in 1:length(vecmu)){
  for(j in 1:length(vecphi)){
    m=vecmu[i]
    p=vecphi[j]
    values[i,j]=exp(sum(dgamma(unc, shape=m/p, scale=p, log = TRUE)))
  }
}


###################################################
### code chunk number 13: gamma2Plot
###################################################
getOption("SweaveHooks")[["fig"]]()
contour(x=vecmu, y=vecphi, z=values*10^faktor, nlevels=15, xlab=math(mu), ylab=math(phi), xaxs="i", yaxs="i", labcex=1)


###################################################
### code chunk number 14: likelihoodfunktionen1
###################################################
getOption("SweaveHooks")[["fig"]]()
vec=seq(from=0.001, to=0.799, by=0.01)
values=dbinom(x=2, size=10, prob=vec)
matplot(vec, values, xlab=math(pi), ylab=math(L(pi)), type="l")
drawml(0.2, dbinom(2,10,prob=0.2))


###################################################
### code chunk number 15: likelihoodfunktionen2
###################################################
getOption("SweaveHooks")[["fig"]]()
values=dbinom(x=2, size=10, prob=vec)/dbinom(2,10,prob=0.2)
matplot(vec, values, xlab=math(pi), ylab=math(tilde(L)(pi)), type="l")
drawml(0.2, 1)


###################################################
### code chunk number 16: likelihoodfunktionen3
###################################################
getOption("SweaveHooks")[["fig"]]()
values=dbinom(x=2, size=10, prob=vec, log=T)
matplot(vec, values, xlab=math(pi), ylab=math(l(pi)), type="l")
drawml(0.2, dbinom(2,10,0.2,log=T))


###################################################
### code chunk number 17: likelihoodfunktionen4
###################################################
getOption("SweaveHooks")[["fig"]]()
values=dbinom(x=2, size=10, prob=vec, log=T)-dbinom(2,10,0.2,log=T)
matplot(vec, values, xlab=math(pi), ylab=math(tilde(l)(pi)), type="l")
drawml(0.2, 0)


###################################################
### code chunk number 18
###################################################

x1 = 233
x2 = 385
x3 = 129
n <- x1+x2+x3
qml <- (2*x1+x2)/(2*n)
se.qml <- sqrt(qml*(1-qml)/n)

pi1ml <- qml^2
pi2ml <- 2*qml*(1-qml)
pi3ml <- (1-qml)^2


###################################################
### code chunk number 19: pbcTreatLambdaMlBerechnung
###################################################
plaSum <- with(pbcTreat, sum (time))
lambdaml <- with (pbcTreat, sum (d) / sum (time))


###################################################
### code chunk number 20: darmkrebs1
###################################################
## Truncated binomial log-likelihood function
## pi: the parameter, the probability of a positive test result
## data: vector with counts  Z_1, ..., Z_N
log.likelihood <- function(pi, data)
{
    n <- sum(data)
    k <- length(data)
    vec <- seq_len(k)
    result <- sum(data * (vec * log(pi) + (k-vec) * log(1-pi))) - 
           n * log(1 - (1-pi)^k)

}
data <- c(37, 22, 25, 29, 34, 49)
eps <- 1e-10
result <- optim(0.5, log.likelihood, data = data, 
                method = "L-BFGS-B", lower = eps, upper = 1-eps, 
                control = list(fnscale = -1), hessian = TRUE) 
(ml <- result$par)


###################################################
### code chunk number 21: beobfisherPrinten
###################################################
(observed.fisher <- - result$hessian)


###################################################
### code chunk number 22: emalgorithmImplementieren
###################################################
## data set
fulldata <- colonCancer
k <- 0:6
n <- sum(fulldata[-1])

## impute start value for Z0 (first element)
## and initialise some different old value
fulldata[1] <- 10
Z0old <- 9

## the EM algorithm
while(abs(Z0old - fulldata[1]) >= 1e-7)
{
    Z0old <- fulldata[1]
    pi <- sum(fulldata * k) / sum(fulldata) / 6
    xi <- (1-pi)^6
    fulldata[1] <- n * xi / (1-xi)
}


###################################################
### code chunk number 23: emErgebnisseTable
###################################################
data <- colonCancer
n <- sum(data, na.rm = TRUE)
Z0 <- 10

fulldata <- c(Z0, 37, 22, 25, 29, 34, 49)
k <- c(0:6)
pi <- 0.5


# start EM algorithm

diff <- 1

iter.matrix = matrix(nrow = 1, ncol = 3)
iter = 0

while(diff >= 1e-7){
  iter = iter + 1
  old.pi <- pi
  pi <- sum(fulldata*k)/sum(fulldata)/6
  xi <- (1-pi)^6
  Z0 <- n*xi/(1-xi)
  iter.matrix = rbind(iter.matrix, c(iter, pi, Z0))
  diff <- abs(old.pi - pi)
  fulldata[1] <- Z0
}

iter.table <- as.data.frame(iter.matrix[-1,])
names(iter.table) <- c("\\textbf{Iteration}", "$\\boldsymbol{\\hat{\\pi}}$", "$\\boldsymbol{\\hat{Z}_{0}}$")
iter.table[,1] <- as.character (iter.table[,1]) # um unschöne 1.0000 zu vermeiden

#in Latex konvertieren
w <- latex(iter.table,
           booktabs = TRUE,
           file = latexTempFile,
           label = "tab:EMergebnisse",
           rowname = NULL,
           caption = "Values of the EM algorithm until convergence (Difference between old and new estimate $\\hat{\\pi}$ smaller than $\\mathsf{10^{-7}}$).",
           center = "none",
           dec = 7,
           numeric.dollar = FALSE,
           collabel.just = Cs(c,c,c),
           col.just = Cs(A,A,A),
           where = NULL
           )
postLatex (w, widthFactor = 0.7)


###################################################
### code chunk number 24
###################################################
beob <- 11
erwartet <- 3.04


###################################################
### code chunk number 25: loglikePoisson
###################################################
getOption("SweaveHooks")[["fig"]]()
    theta.ml = beob/erwartet
    theta.vector = seq(from = 1.5, to = 7, length = 200)
    rel.loglike = function(theta){
      return( sum(dpois(beob, theta*erwartet, log = TRUE)) - sum(dpois(beob, theta.ml*erwartet, log = TRUE)))
    }
    matplot(theta.vector, sapply(theta.vector, rel.loglike), xlab = math (lambda), ylab = math (tilde(l)(lambda)), type = "l")
    drawml(theta.ml, 0)
    quadrApprox = function(theta){
      return(-0.5*beob/(theta.ml^2)*(theta - theta.ml)^2)
    }
    lines(theta.vector, quadrApprox(theta.vector), lty = 2)


###################################################
### code chunk number 26: approxBinom
###################################################
getOption("SweaveHooks")[["fig"]]()
binomexample=function(x, n, xlimits){
vec=seq(from=max(x/n-0.2, 0.01), to=min(x/n+0.2,1), by=0.001)
values=dbinom(x, size=n, prob=vec, log=T)-dbinom(x,n,x/n,log=T)
matplot(vec, values, xlab=math(pi), ylab=math(tilde(l)(pi)), type="l", xlim = xlimits, ylim = c(-5, 0))
drawml(x/n, 0)
approx=-1/2*n/(x/n*(1-x/n))*(vec-x/n)^2
lines(vec, approx, lty=2)
}
binomexample(8,10, xlimits = c(0.6, 1))


###################################################
### code chunk number 27: approxBinomB
###################################################
getOption("SweaveHooks")[["fig"]]()
binomexample(40, 50, xlimits = c(0.6, 1))


###################################################
### code chunk number 28: approxBinomC
###################################################
getOption("SweaveHooks")[["fig"]]()
binomexample(160, 200, xlimits = c(0.6, 1))


###################################################
### code chunk number 29: approxBinomD
###################################################
getOption("SweaveHooks")[["fig"]]()
binomexample(800, 1000, xlimits = c(0.6, 1))


###################################################
### code chunk number 30: uniform-likelihood
###################################################
getOption("SweaveHooks")[["fig"]]()
max=7
vecn=c(5,6,7)
vectheta=seq(from=max, to=max+10, by=0.1)

likelihood=function(theta, num){
return(theta^-num)
}

#Linie von 0 bis max zeichnen
upper=likelihood(max,vecn[1])*10^5
matplot(c(0,max),c(0,0), xlim=c(0,max+10),ylim=c(0,upper), type="l", xlab=math(theta), ylab=math(L(theta)%*%10^5))

#likelihood vervollständigen
for(i in 1:length(vecn)){
  values=sapply(vectheta, likelihood, num=vecn[i])
  lines(vectheta, values*10^5, lty=i)
}

drawml(max, upper)


###################################################
### code chunk number 31: uniform-loglikelihood
###################################################
getOption("SweaveHooks")[["fig"]]()
max=7
vecn=c(5,6,7,10,30,60)
vectheta=seq(from=max, to=max+10, by=0.1)

loglikelihood=function(theta, num){
  return(log(theta)*(-num))
}
#Grenzen für y-Achse
upper=loglikelihood(max,num=vecn[1])
lower=loglikelihood(max+10,num=vecn[length(vecn)])
#Startplot
values=sapply(vectheta, loglikelihood, num=vecn[1])
matplot(vectheta, values, xlab=math(theta), ylab=math(l(theta)), ylim=c(lower,upper), type="l")

#für weiter n's:
for(i in 2:length(vecn)){
  values=sapply(vectheta, loglikelihood, num=vecn[i])
  lines(vectheta, values, lty=i)
}

drawml(max, upper)


