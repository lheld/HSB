#################################################################
### R code for Chapter 9
#################################################################

### Encoding: UTF-8

# Remark: code chunk number 8 depends on chap7.R, run that file before chunk 8


###################################################
### code chunk number 1
###################################################
mu.hat <- as.numeric(mean(iten[["tf500a1"]]))
sd.hat <- as.numeric(sd(iten[["tf500a1"]]))
fac <- as.numeric(sqrt(1+1/length(iten[["tf500a1"]])))
q0 <- pnorm(2000, mu.hat, sd.hat)


###################################################
### code chunk number 2: poisson1
###################################################

n <- 1
x <- 11
ex <- 3.04
ey <- 3.04
thetahat <- x/ex

  predlikun = function(yv){
    sapply(yv, function(y){
      gesamtsumme = sum(c(y, x))
      logerg = y*log(ey) + gesamtsumme * log (gesamtsumme/(ex+ey)) - lfactorial(y) - gesamtsumme
      exp(logerg)
    })
  }

support = 0:100
fhat <- dpois(support, thetahat*ex)
pt <- predlikun(support)
normierung <- sum(pt)
fsupport <- pt/normierung
pred.ewert <- sum(support * fsupport)
pred.var <- sum ((support-pred.ewert)^2 * fsupport)



###################################################
### code chunk number 3: poisson2
###################################################
## Data:
x <- 11
ex <- 3.04
ey <- 3.04
## MLE:
(lambdahat <- x/ex)
## bootstrap prediction:
set.seed(1)
m <- 10000
lambdasample <- rpois(m, lambda=lambdahat * ex) / ex
support <- 0:100
ghat <- rep(0, length(support))
for(i in 1:m)
{
    ghat <- ghat + dpois(support, lambda=lambdasample[i] * ey)
}
ghat <- ghat/m
## empirical moments of the predictive distribution:
(e.g <- sum(ghat*support))
e.g2 <- sum(ghat*support^2)
(var.g <- e.g2 - e.g^2)


###################################################
### code chunk number 4: PrognosePoisson
###################################################
getOption("SweaveHooks")[["fig"]]()
farben = c("black", "darkgray", "grey")
                                        # links mehr Platz für ylab schaffen
oldmar = par("mar")
par(mar = oldmar * c(1, 2, 1, 1))
matplot(cbind(support-0.18, support, support + 0.18), cbind(fhat, fsupport, ghat), type = "h", xlim = c(0, 30),
        col = farben, lwd = 4, lty = 1,
        xlab = math (x),
        ylab = expression ("Predictive probability mass function")
        )
legend("topright", legend = Cs("plug-in", "likelihood", "bootstrap"), col = farben, lwd = 3, lty = 1, bty = "n")
par(mar = oldmar)


###################################################
### code chunk number 5
###################################################

nu <- 2000
delta <- 1/(200^2)

sigma <- as.numeric(sd(iten[["tf500a1"]]))
kappa <- 1/sigma^2


pred.mean <- function(x, k0, mean0, sigma){
    mean = mean(x)
    n = length(x)
    k = 1 / sigma^2
    postvar <- (n * k + k0)^(-1)
    postmean <- (n * k * mean + k0 * mean0) * postvar
    return(postmean)
}

pred.sd <- function(x, k0, sigma){
    mean = mean(x)
    sd = sigma
    n = length(x)
    k = 1 / sd(x)^2
    num <- k*sd^2+n+1
    den <- k*sd^2+n
    return(sd*sqrt(num/den))
}

pm.total <- pred.mean(x=iten$tf500a1,  k0=1/(200^2), mean0=2000, sigma=sigma)
pm.male <- pred.mean(x=iten$tf500a1[iten$Sex=="m"],  k0=1/(200^2), mean0=2000, sigma=sigma)
pm.female <- pred.mean(x=iten$tf500a1[iten$Sex=="w"],  k0=1/(200^2), mean0=2000, sigma=sigma)
    
psd.total <- pred.sd(x=iten$tf500a1,  k0=1/(200^2), sigma=sigma)
psd.male <- pred.sd(x=iten$tf500a1[iten$Sex=="m"],  k0=1/(200^2), sigma=sigma)
psd.female <- pred.sd(x=iten$tf500a1[iten$Sex=="w"],  k0=1/(200^2), sigma=sigma)
    



###################################################
### code chunk number 6
###################################################
a <- 6
b <- 2
set.seed(25012013)
pi <- 0.55
m <- 11
data <- rbinom(m, prob=pi, size=1)

a.cum <- a+cumsum(data)
b.cum <- b+c(1:m)-cumsum(data)
ab.cum <- a.cum + b.cum
pred.prob <- a.cum/ab.cum


###################################################
### code chunk number 7: poisson3
###################################################
prognosis <- 
    data.frame(dist = 
               c ("Plug-in", "Likelihood", "Bootstrap", "Bayesian"),  
               exp = c(lambdahat*ex, pred.ewert, 11, 11.5),
               var = c(lambdahat*ex, pred.var, 22, 23))

names(prognosis) <- c("Predictive distribution",
                      "Mean",
                      "Variance")

 #in Latex konvertieren
w <-
 latex(prognosis, # was konvertieren?
       file = latexTempFile,
       label = "tab:Prognoseverteilungen", # labeln
       rowname = NULL, # um keine Zeilennamen zu haben
       booktabs = TRUE, # um booktabs zu benutzen
       caption = "Mean and variance of different predictive distributions for $x=\\mathsf{11}$ and $e_x=e_y=\\mathsf{3.04}$. Exact values are given for the bootstrap prediction based on \\eqref{eq:boot1}
and \\eqref{eq:boot2}.",
       center = "none", # statt \\begin{center}
       dec = 3, # für Dezimalstellen
       colnamesTexCmd = "bfseries",
       numeric.dollar = FALSE,
       collabel.just = Cs(L, C, C), # Spaltenüberschriften ausrichten
       col.just = Cs(L, A, A), # Spalten ausrichten
       where = "tbp" # default-Werte wiederherstellen
      )
postLatex (w, widthFactor = 0.9)


###################################################
### code chunk number 8
###################################################
# load results from chap7.R
load("../data/postprob.RData")


###################################################
### code chunk number 9: fussball1
###################################################
## Daten erzeugen
prognosis <- c ("Prevalence", "Expert", "Oracle", "Inverted oracle")
auc <- c(0.5, 1, 1, 0)
s <- c(0, 0.04, 0, 1)
m <- c(0, 0.25, 0.25, 0.25)
bs <- 0.25 + s - m

fb <- cbind (prognosis, auc, s, m, bs)

## ausgeben
latexMatrix (fb)


###################################################
### code chunk number 10: PIT1
###################################################
    set.seed(1)
    M <- 10000
    x <- rnorm(M)
    y <- rnorm(M)
    pit.plugin <- pnorm(y, mean=x, sd=1)
    pit.bayes <- pnorm(y, mean=x, sd=sqrt(2))


###################################################
### code chunk number 11: PIT2
###################################################
logscore <- function(yobs, mu, sigma){
      return(-dnorm(yobs, mu, sigma, log=T))
    }

    crps <- function(yobs, mu, sigma){
      yt <- (yobs-mu)/sigma
      term <- yt*((2*pnorm(yt)-1)+2*dnorm(yt)-1/sqrt(pi))
      return(term*sigma)
    }

    mse <- function(yobs, mu){
      return((yobs-mu)^2)
    }

    plugin.ls <- mean(logscore(y, x, 1))
    bayes.ls <- mean(logscore(y, x, sqrt(2)))
    oracle.ls <- mean(logscore(y, 0, 1))


    plugin.crps <- mean(crps(y, x, 1))
    bayes.crps <- mean(crps(y, x, sqrt(2)))
    oracle.crps <- mean(crps(y, 0, 1))

    plugin.mse <- mean(mse(y, x))
    bayes.mse <- mean(mse(y, x))
    oracle.mse <- mean(mse(y, 0))



###################################################
### code chunk number 12: PIThistogramme_plugin
###################################################
getOption("SweaveHooks")[["fig"]]()
histbreaks <- seq (0, 1, by = 0.05) # if omitted, truehist plots something greater than 1
truehist(pit.plugin, breaks = histbreaks, col = "gray", ymax=2.5, ylab = "", xlab = "PIT")


###################################################
### code chunk number 13: PIThistogramme_bayes
###################################################
getOption("SweaveHooks")[["fig"]]()
truehist(pit.bayes, breaks = histbreaks, col = "gray", ymax=2.5, ylab = "", xlab = "PIT")


###################################################
### code chunk number 14: fussball1
###################################################
## Daten erzeugen
plugin <- c(plugin.ls, plugin.crps, plugin.mse)
bayes <- c(bayes.ls, bayes.crps, bayes.mse)
oracle <- c(oracle.ls, oracle.crps, oracle.mse)

data <- cbind (c ("Plug-in", "Bayesian", "Oracle"),
               formatRound (rbind (plugin, bayes, oracle), 4)
               )

## ausgeben
latexMatrix (data)


