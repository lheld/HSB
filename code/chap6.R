#################################################################
### R code for Chapter 6
#################################################################

### Encoding: UTF-8

# Remark: code chunk number 31 depends on chap4.R, run that file before chunk 31

###################################################
### code chunk number 1
###################################################

sens <- 0.9
spec <- 0.9
prev <- 0.01
pt <- sens*prev+(1-spec)*(1-prev)
pdt <- sens*prev/pt

prior.odds <- 0.01/0.99
lr.pos <- sens/(1-spec)
posterior.odds.pos <- prior.odds*lr.pos
PPV <- posterior.odds.pos/(1+posterior.odds.pos)
lr.neg <- (1-sens)/spec
posterior.odds.neg <- prior.odds*lr.neg
NPV <- 1-(posterior.odds.neg/(1+posterior.odds.neg))


###################################################
### code chunk number 2: binomialPosteriori
###################################################
getOption("SweaveHooks")[["fig"]]()
    alpha = 3
    beta = 2
    n = 10
    x = 8
    alpha.p = alpha + x
    beta.p = beta + n - x
    grid = seq(0, 1, length = 200)
    modus = (alpha.p - 1)/(alpha.p + beta.p - 2)
    # Priori-Dichte startplot
    plot(grid, dbeta(grid, alpha, beta),
         xlab = math (pi), ylab = math (f * (pi~plain("|")~x)),
         type = "l", lty = 2, ylim = c(0, dbeta(modus, alpha.p, beta.p))
         )
    # Masse grau markieren (zuerst, damit Dichte drüberkommt)
    bounds = qbeta(c(0.025, 0.975), alpha.p, beta.p)
    # links
    left.grid = seq(0, bounds[1], length = 300)
    lines(left.grid, dbeta(left.grid, alpha.p, beta.p), type = "h", col = "lightgray")
    # rechts
    right.grid = seq(bounds[2], 1, length = 300)
    lines(right.grid, dbeta(right.grid, alpha.p, beta.p), type = "h", col = "lightgray")
    # Posteriori Dichte
    lines(grid, dbeta(grid, alpha.p, beta.p), lty = 1)
    mean = (alpha.p)/(alpha.p + beta.p)
    points(mean,par("usr")[3], pch = 19, xpd = TRUE)
    meanFormatted <- format(mean, digits = 3)
    text(mean, par("usr")[3],
         labels = substitute(E*(paste(italic (pi)," | ", italic (x))) == m, list(m = meanFormatted)),
         pos = 3)
    drawml(modus, dbeta(modus, alpha + x, beta + n - x), down = TRUE)
    modusFormatted <- format(modus, digits=3)
    legend("topleft", lty = c(1,2), legend = c("posterior", "prior"), bty = "n")
    # gleichendiges KI einzeichnen
    funline(bounds[1], dbeta, shape1 = alpha.p, shape2 = beta.p, down = TRUE)
    funline(bounds[2], dbeta, shape1 = alpha.p, shape2 = beta.p, down = TRUE)



###################################################
### code chunk number 3
###################################################
  #Tabelle erzeugen
  post = function(one){
    n = one[1]
    x = one[2]
    return(c((x+1)/(n+2), x/n, qbeta(0.5, 1+x, 1+n-x), qbeta(0.025, 1+x, 1+n-x), qbeta(0.975, 1+x, 1+n-x)))
  }
  nx = matrix(c(rep(10,3), rep(100,3), 0, 1, 5, 0, 10, 50), ncol = 2)
  result = apply(nx, 1, post)
  output = cbind(nx, t(result))
  output <- as.data.frame (output)
  names(output) <- c("$\\boldsymbol{n}$", "$\\boldsymbol{x}$",
                     "\\textbf{Mean}", "\\textbf{Mode}", "\\textbf{Median}",
                     "\\textbf{2.5\\% quantile}", "\\textbf{97.5\\% quantile}")
  for (i in 1:2)
    output[,i] <- as.character (output[,i])
  #in Latex konvertieren
w <-
    latex(output,                          # was konvertieren?
          file = latexTempFile,                 # um direkt auf Ausgabe zu schreiben
          where = "tbp",
          booktabs = TRUE,
          label = "tab:posterioriGleich",  # labeln
          cgroup = c("Observation", "Posterior characteristics"),    # um Spalten zusammenzufassen
          n.cgroup = c(2,5),
          rowname = NULL, # um keine Zeilennamen zu haben
          caption =
  "Summary characteristics of the posterior distribution of $\\pi$ under a uniform prior distribution in the binomial model. The 95\\% credible interval based on the 2.5\\% and 97.5\\% quantiles should be compared with the 95\\% confidence intervals from Table~\\ref{tab:anteilZahlen}.", 
          col.just = c (rep ("R", 2), rep ("C", 7)),
          collabel.just = c (rep ("R", 2), rep ("C", 7)),
          center = "none",
          dec = 2,                      # for 2 Dezimalstellen
          numeric.dollar = FALSE,
          cgroupTexCmd = "bfseries",
          colnamesTexCmd = "footnotesize"
          )
postLatex (w, widthFactor = 1.0)


###################################################
### code chunk number 4
###################################################

a <- 0.5
b <- 5
x <- 1
n <- 100
apost <- a+x
bpost <- b+n-x

pmean <- apost/(apost+bpost)
pmode <- (apost-1)/(apost+bpost-2)
pmedian <- qbeta(0.5, apost, bpost)
plower <- qbeta(0.025, apost, bpost)
pupper <- qbeta(0.975, apost, bpost)
  


###################################################
### code chunk number 5: diagnosticTesting2-ppv
###################################################
getOption("SweaveHooks")[["fig"]]()
post.ppv <- function(theta, a, b, LR){
  c <- a * LR / b
  x <- c * (1 / theta - 1)
  result <- c * theta^(-2) * df(x, df1=2*b, df2=2*a)
  return(result)
}

target.ppv <- function(x, a, b, LR)
  return(post.ppv(theta=x, a=a, b=b, LR=LR) * x)

ppv.mean <- integrate(target.ppv, 
                      lower=0, upper=1, 
                      a=apost, b=bpost, LR=9)$value

ppv.mode <- optimize(post.ppv, 
                     lower=0, upper=1, 
                     a=apost, b=bpost, LR=9, 
                     maximum=TRUE)$maximum

# nsample: sample size
nsample <- 10^5
# prev: samples from beta(1.5, 99.5) distribution
prev <- rbeta(nsample, apost, bpost)
# values for sensitivity and specificity
sens <- 0.9
spec <- 0.9
ppv <- sens*prev/(sens*prev+(1-spec)*(1-prev))
# make a nice histogram
truehist(ppv[ppv<0.5], xlab="Positive predictive value", col="gray")

grid <- c(1:1000)/2000
f.ppv <- post.ppv(grid, a=apost, b=bpost, LR=9)
lines(grid, f.ppv, type="l", col="red")


###################################################
### code chunk number 6: diagnosticTesting2-npv
###################################################
getOption("SweaveHooks")[["fig"]]()
post.npv <- function(tau, a, b, LR){
    d <- b * LR / a
    x <- d * (1 / tau - 1)
    result <- d * tau^(-2) * df(x, df1=2*a, df2=2*b)
    return(result)
}

target.npv <- function(x, a, b, LR)
  return(post.npv(tau=x, a=a, b=b, LR=LR) * x)

npv.mean <- integrate(target.npv, 
                      lower=0, upper=1, 
                      a=apost, b=bpost, LR=9)$value

npv.mode <- optimize(post.npv, 
                     lower=0, upper=1, 
                     a=apost, b=bpost, LR=9, 
                     maximum=TRUE)$maximum
  
## calculate npv samples
npv <- spec * (1 - prev) / ((1 - sens) * prev + spec * (1-prev))
# make a nice histogram
truehist(npv[npv>0.99], xlab="Negative predictive value", col="gray")#, xlim=c(0,0.5),ylim=c(0,50))
grid.npv <- seq(from=0.99, to=0.999999, length=200)
f.npv <- post.npv(grid.npv, a=apost, b=bpost, LR=9)
lines(grid.npv, f.npv, type="l", col="red")


###################################################
### code chunk number 7
###################################################
## prior parameters
a <- 0.5
b <- 5
## data
x <- 1
n <- 100
## posterior parameters
apost <- a+x
bpost <- b+n-x
## sample size
nsample <- 10^5
## prevalence values sampled from Be(1.5, 99.5) distribution
prev <- rbeta(nsample, apost, bpost)
## set values for sensitivity and specificity
sens <- 0.9
spec <- 0.9
## compute resulting positive and negative predictive value samples
ppv <- sens * prev / (sens * prev + (1 - spec) * (1 - prev))
npv <- spec * (1 - prev) / ((1 - sens) * prev + spec * (1 - prev))


###################################################
### code chunk number 8: vergleichHPD
###################################################
getOption("SweaveHooks")[["fig"]]()
  alpha = 3
  beta = 2
  n = 10
  x = 8
  p1 = alpha + x
  p2 = beta + n - x
  grid = seq(0, 1, length = 400)
  modus = (alpha + x - 1)/(alpha + beta + n - 2)
  plot(grid, dbeta(grid, p1, p2), xlab = math (pi),xlim=c(0.2,1.0), 
  ylab = math (f * (pi~plain("|")~x)), type = "l", xaxs = "i", yaxs = "i",
  ylim = c(0, dbeta(modus, p1, p2) + 0.2)
  )
  drawml(modus, dbeta(modus, alpha + x, beta + n - x), down = TRUE)

  # gleichendiges Kredibilitätsintervall

  kr.l = qbeta(0.025, p1, p2)
  kr.u = qbeta(0.975, p1, p2)
  kr.lim =c(kr.l, kr.u)

  axis(side = 1, at = kr.lim, labels = FALSE)
  segments(kr.lim, rep(par("usr")[3],2), kr.lim, dbeta(kr.lim, p1, p2), lty = 2)
  # Flächen füllen
  x.krl = c(grid[grid <= kr.l], kr.l, 0)
  y.krl = c(dbeta(grid[grid <= kr.l], p1, p2), 0, 0)
  polygon(x.krl, y.krl, col = "gray", border = 1, lty = 2)

  x.kru = c(grid[grid >= kr.u], 1, kr.u)
  y.kru = c(dbeta(grid[grid >= kr.u], p1, p2), 0, 0)
  polygon(x.kru, y.kru, col = "lightgray", border = 1, lty = 2)

  # HPD-Intervall:

  result = uniroot(function(h){outerdens(h, p1, p2)[1] - 0.05}, interval = c(1e-15, dbeta(modus, p1, p2)-1e-15))
  height = result$root
  abline(h = height, lty = 2)

  hpd.lim = outerdens(height, p1, p2)[c(2,3)]
  hpd.l = hpd.lim[1]
  hpd.u = hpd.lim[2]

  axis(side = 1, at = hpd.lim, labels = FALSE)
  # Flächen schraffieren
  x.hpdl = c(grid[grid <= hpd.l], hpd.l, 0)
  y.hpdl = c(dbeta(grid[grid <= hpd.l], p1, p2), 0, 0)
  polygon(x.hpdl, y.hpdl, density = 10, col = 1, border = 1, lty = 1)

  x.hpdu = c(grid[grid >= hpd.u], 1, hpd.u)
  y.hpdu = c(dbeta(grid[grid >= hpd.u], p1, p2), 0, 0)
  polygon(x.hpdu, y.hpdu, density = 10, col = 1, border = 1, lty = 1)




###################################################
### code chunk number 9: CapRecapBayes
###################################################

  # Plot von priori und posteriori for den Parameter N im
  # capture-recapture experiment

  # input parameter:
  # M: Anzahl der markierten Fische
  # n: Größe der Stichprobe
  # x: Anzahl der markierten Fische in Stichprobe
  # gamma: Parameter zur Festlegung der Priori distribution:
  # gamma = 0 --> Gleichverteilung
  # gamma > 0 --> (gestutzte) geometrische distribution
  # maxN: Obere Grenze for den Träger der priori-distribution


  plot.posterior <- function(M, n, x, gamma=0.0, maxN=5*M*n/x, level=0.95, noplot=F, maxNplot=5*M*n/x, ylim=NA){

    if (maxN == Inf) maxN <- 10000

    # prior, definiert auf 1:maxN
    if(gamma>0)
      prior <- gamma*(1-gamma)^{c(1:maxN)-1}
    if(gamma==0) prior <- 1^{c(1:maxN)-1}
    
    # prior ist null for N < M
    prior[1:(M-1)] <- prior[1:(M-1)]*0

    # Normalisieren der priori
    prior <- prior/sum(prior)

    # N: Vektor zum Berechnen und Plotten der Priori und Posteriori

    N <- c(1:maxN)

    # posterior, definiert auf 1:maxN
    posterior <- vector(mode="numeric",length=maxN)
    posterior <- posterior*0

    traeger <- c(max(n, M+n-x):maxN)

    # Berechnung der Posteriori
    posterior[traeger] <- prior[traeger]*dhyper(x, M, traeger-M, n)

    # Normalisieren der posterior
    posterior <- posterior/sum(posterior)

    # Posteriori Modus
    pmodus <- N[posterior==max(posterior)]

    # Posteriori Erwartungswert
    perwartung <- round(sum(posterior[traeger]*traeger), 1)

    # Posteriori Median
    cumposterior <- cumsum(posterior)
    pmedian <- which.max(cumposterior>=.5)

    # Berechnung der HPD Region

    postorder <- order(posterior, decreasing = T)
    postsum <- cumsum(posterior[postorder])

    hpd <- N[postorder[1:which.max(postsum>=level)]]
    
    hpdmin <- min(hpd)
    hpdmax <- max(hpd)

    if(noplot==F){
      par(mfrow=c(2,1))
      maxpriorposterior <- max(prior,posterior)
      if(is.na(ylim))
        ylim <- c(0,maxpriorposterior)
      
      # plot der priori
      plot(N, prior,type="p", cex=0.25, pch=19, xlab=math (N), ylab = math (f * (N)), xlim=c(0, maxNplot), ylim=ylim)

      # plot der posteriori
      plot(N, posterior,type="p", cex=0.25,pch=19, xlab=math (N), ylab = math (f * (N~plain("|")~x)), xlim=c(0, maxNplot), ylim=ylim)



      # Eintragen der Punktschätzer
      abline(v = c(pmodus), lty = 2)
      abline(v = c(perwartung), lty = 2)
      abline(v = c(pmedian), lty = 2)


      lines(x=c(hpdmin,hpdmax), y=c(0,0), lty=1)
      lines(x=c(hpdmin,hpdmin), y=c(0,maxpriorposterior), lty=1)
      lines(x=c(hpdmax,hpdmax), y=c(0,maxpriorposterior), lty=1)
    }

    return(c(pmodus, pmedian, perwartung, hpdmin, hpdmax))

  }
  res <- plot.posterior(M=26, n=63, x=5, gamma=0, maxN=1500, noplot=T)



###################################################
### code chunk number 10: CapRecapBayesA
###################################################
getOption("SweaveHooks")[["fig"]]()
mygamma <- 0.001
      res <- plot.posterior(M=26, n=63, x=5, gamma=0, maxN=1500, maxNplot=1500,ylim=c(0,0.003))
      res <- plot.posterior(M=26, n=63, x=5, gamma=mygamma, maxN=Inf, noplot=T,ylim=c(0,0.003))


###################################################
### code chunk number 11: CapRecapBayesB
###################################################
getOption("SweaveHooks")[["fig"]]()
      res2 <- plot.posterior(M=26, n=63, x=5, gamma=mygamma, maxN=Inf, maxNplot=1500,ylim=c(0,0.003))


###################################################
### code chunk number 12: normal-posterior-setup
###################################################
 var0 <- 200^2
 k0 <- 1/var0
 mean0 <- 2000

 x <- iten$tf500a1
 alcoData <- list(x = x, mean = mean(x), sd = sd(x), n = length(x),
                  k = 1 / sd(x)^2)

 postvar <- with(alcoData, (n * k + k0)^(-1))
 postmean <- with(alcoData, (n * k * mean + k0 * mean0) * postvar)


###################################################
### code chunk number 13: normal-posterior
###################################################
 getOption("SweaveHooks")[["fig"]]()
 par(mfrow=c(1,1))
 margin <- postvar
 gridRange <- range(x)
 grid <- seq(gridRange[1], gridRange[2], length = 10000)

 plot(grid, dnorm(grid, postmean, sqrt(postvar)),
      xlab = math (mu),
      ylab = math (hat (f)*(mu~plain("|")~x[1*plain(":")*n])),
      type = "l")

 lines(grid, dnorm(grid, mean0, sqrt(var0)), lty = 2)

 rug(x)
    drawml(postmean, dnorm(postmean, postmean, sqrt(postvar)), digi = 5, down = TRUE)
    legend("topleft", lty = c(1,2), legend = c("posterior", "prior"), bty = "n")


###################################################
### code chunk number 14: semiBayes
###################################################

 a <- 6
 b <- 2
 c <- 102
 d <- 101

 data.mean <- log((a*d)/(b*c))
 data.var <- 1/a+1/b+1/c+1/d
 data.ef <- exp(1.96*sqrt(data.var))

 prior.mean <- 0
 prior.var <- 1

 posterior.var <- 1/(1/data.var+1/prior.var)
 posterior.mean <- posterior.var*data.mean/data.var
 posterior.ef <- exp(1.96*sqrt(posterior.var))




###################################################
### code chunk number 15
###################################################
  ## setting
  sizes <- c (rep(10, 3), rep(100, 3))
  data <- data.frame (n = sizes,
  x = c(0,1,5,0,10,50)
  )
  prob = 0.95                             # Konfidenzniveau?
  q = qnorm(1-(1-prob)/2)                 # nötiges Normalverteilungsquantil (bei 0.95 : 1.96...)

  ## Wald

  wald.se <- with (data, sqrt (x/n * (1 - x/n) / n))
  wald.lower <- with (data, x/n - q * wald.se)
  wald.upper <- with (data, x/n + q * wald.se)

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
      likelihood.ci(alpha = 1 - prob, loglik.binom, theta.hat = 1 - eps,
      lower = eps, upper = 1 - eps, x = x[i], size = n[i], comp.upper = FALSE)
      )
      res <- c (res, 1)
    }
    else if (keinTreffer[i]){
      res <- with (data,
      likelihood.ci(alpha = 1 - prob, loglik.binom, theta.hat = eps,
      lower = eps, upper = 1 - eps, x = x[i], size = n[i], comp.lower = FALSE)
      )
      res <- c (0, res)
    }
    else{
      res <- with (data,
      likelihood.ci (alpha = 1 - prob, loglik.binom, theta.hat = x[i] / n[i],
      lower = eps, upper = 1 - eps, x = x[i], size = n[i])
      )
    }
    lik.lower[i] <- res[1]
    lik.upper[i] <- res[2]
  }

  ## Hpd mit Jeffreys' Priori

  hpd.jeffreys.interval = function(x, n){
    # Parameter der Posteriori beta:
    p1 = 0.5 + x
    p2 = 0.5 + n - x
    modus = (p1-1)/(p1 + p2 - 2)
    if(x==0){
      hpd.l <- 0
      hpd.u <- qbeta(prob, p1, p2)
    }
    else{
      if(x==n){
        hpd.l <- qbeta(1-prob, p1, p2)
        hpd.u <- 1
      }
      else{
        eps = 1e-15
                                        # Funktion outerdens liefert for den Wert h die Wahrscheinlichkeit aller
                                        # Werte pi, deren Dichte kleiner als h ist, sowie die Schnittpunkte.
        result = uniroot(function(h){outerdens(h, p1, p2)[1] - (1-prob)}, interval = c(eps, dbeta(modus, p1, p2)-10*eps))
        height = result$root
        hpd.l = outerdens(height, p1, p2)[2]
        hpd.u = outerdens(height, p1, p2)[3]
      }
    }
    hpd.lim = c(hpd.l, hpd.u)
    return(hpd.lim)
 }

  HPD.int <- t (apply (data, 1, function (one) hpd.jeffreys.interval (one["x"], one["n"])))

  ## gleichendiges KI
  gleich.int.lower <- with (data,
  qbeta ((1 - prob) / 2,
  0.5 + x,
  0.5 + n - x)
  )
  gleich.int.upper <- with (data,
  qbeta ((1 - prob) / 2 + prob,
  0.5 + x,
  0.5 + n - x)
  )

  ##  alle zusammenführen

  bounds <- origBounds <-
  cbind (gleich.int.lower, gleich.int.upper, HPD.int, wald.lower, wald.upper, lik.lower, lik.upper)
  bounds <- format (round (bounds, 3), 2) # format um immer 3 Dezimalziffern zu haben

form <- function(x)
{
    paste("$\\mathsf{", x, "}$", sep="")
}

  gleich <- paste (form(bounds[,1]), " to ", form(bounds[,2]), sep = "")
  hpd <- paste (form(bounds[,3]), " to ", form(bounds[,4]), sep = "")
  wald <- paste (form(bounds[,5]), " to ", form(bounds[,6]), sep = "")
  lik <- paste (form(bounds[,7]), " to ", form(bounds[,8]), sep = "")

  tableData <- cbind (data, gleich = gleich, hpd = hpd, wald = wald, lik = lik)

  names(tableData) <- c("$\\boldsymbol{n}$", "$\\boldsymbol{x}$", "Equi-tailed", "HPD", "Wald", "Likelihood")
names (tableData) <- paste ("\\textbf{", names (tableData), "}", sep = "")

  #in Latex konvertieren
w <-
    latex(tableData,                             # was konvertieren?
          file = latexTempFile,                        # um direkt auf Ausgabe zu schreiben
          label = "tab:vergleichKRmitHPD",       # labeln
          booktabs = TRUE,
          cgroup = c ("Observation",
          paste ("Interval at ", formatRound (prob * 100), "\\% level of type", sep = "")),
          n.cgroup = c (2, 4),              # wieviele Spalten jeweils zusammenfassen?
          rowname = NULL,                  # um keine Zeilennamen zu haben       booktabs = TRUE,                  # um booktabs zu benutzen
          caption = "A comparison of different credible and confidence intervals for $X\\sim \\Bin(n, \\pi)$ and $\\pi \\sim
  \\Be({1}/{2}, {1}/{2})$. If $n=\\mathsf{100}$ and $x = \\mathsf{50}$, the four approaches yield
nearly identical intervals.",
          cgroupTexCmd = "bfseries",            # um "Posteriori etc" nicht fett zu haben
          colnamesTexCmd = "footnotesize",
          numeric.dollar = FALSE,
          col.just = Cs(R, R, r, r, r, r),  # Spalten ausrichten
          collabel.just = Cs(R, R, c, c, c, c),
          where = "tbp",                   # Floatparameter
          center = "none"              # for centering statt center-Umgebung
          )
postLatex (w, widthFactor = 1.0)

  ## Tabelle mit Breiten:

  breiten <- sapply (seq (1, ncol (origBounds) - 1, by = 2),
  function (lowerInd) diff (t (origBounds[, lowerInd+(0:1)]))
  )
  breiten <- as.data.frame (format (round (breiten, 3), 2))
  breiten <- cbind (data, breiten)
  names (breiten) <- c ("$\\boldsymbol{n}$", "$\\boldsymbol{x}$", "Equi-tailed", "HPD", "Wald", "Likelihood")
names (breiten) <- paste ("\\textbf{", names (breiten), "}", sep = "")
  #in Latex konvertieren
w <-
    latex(breiten,                             # was konvertieren?
          file = latexTempFile,                        # um direkt auf Ausgabe zu schreiben
          booktabs = TRUE,
          label = "tab:vergleichKRmitHPDBreiten",       # labeln
          cgroup = c ("Observation", paste ("Width of ", formatRound (prob * 100),"\\% interval of type", sep = "")),
          n.cgroup = c (2, 4),              # wieviele Spalten jeweils zusammenfassen?
          rowname = NULL,                  # um keine Zeilennamen zu haben       booktabs = TRUE,                  # um booktabs zu benutzen
          caption = "Comparison of the widths of the different confidence and credible intervals from 
Table~\\ref{tab:vergleichKRmitHPD}",
          cgroupTexCmd = "bfseries",
          colnamesTexCmd = "footnotesize",
          col.just = Cs(r, r, C, C, C, C),  # Spalten ausrichten
          collabel.just = Cs(r, r, C, C, C, C),
          where = "tbp",                   # Floatparameter
          center = "none"              # for centering statt center-Umgebung
          )
postLatex (w, widthFactor = 1.0)


###################################################
### code chunk number 16: anteilJeffreysHPD
###################################################
getOption("SweaveHooks")[["fig"]]()
      # Szenario festlegen:
      n1 = 50
      x1 = 0:n1

      prob = 0.95 # Konfidenzniveau?
      q = qnorm(1-(1-prob)/2)
      ygrenzen = c(0.8, 1) # y-Skala?

      n2 = 100
      x2 = 0:n2

      # HPD-Intervalle berechnen

      # Funktion hpd.jeffreys.interval berechnet for Beobachtung x bei festem n das HPD-Intervall bei Jeffreys' Priori.
      # Dazu darf x nicht 0 oder n sein! (sonst kein Modus)
      hpd.jeffreys.interval = function(x, n){
        # Parameter der Posteriori beta:
        p1 = 0.5 + x
        p2 = 0.5 + n - x
        modus = (p1-1)/(p1 + p2 - 2)
        eps = 1e-15
        # Funktion outerdens liefert for den Wert h die Wahrscheinlichkeit aller
        # Werte pi, deren Dichte kleiner als h ist, sowie die Schnittpunkte.
        result = uniroot(function(h){outerdens(h, p1, p2)[1] - (1-prob)}, interval = c(eps, dbeta(modus, p1, p2)-10*eps))
        height = result[["root"]]
        hpd.l = outerdens(height, p1, p2)[2]
        hpd.u = outerdens(height, p1, p2)[3]
        hpd.lim = c(hpd.l, hpd.u)
        return(hpd.lim)
      }

      # außer an den Endpunkten unimodale Dichte, deshalb dort:
      HPD.int.middle = t(sapply(x1[-c(1, n1+1)], hpd.jeffreys.interval, n = n1))
      # an den Endpunkten HPD berechnen:
      HPD.int.x10 = c(0, qbeta(prob, 0.5 + 0, 0.5 + n1 - 0))
      HPD.int.x1n1 = c(qbeta(1-prob, 0.5 + n1, 0.5 + n1 - n1),1)

      HPD.int = rbind(HPD.int.x10, HPD.int.middle, HPD.int.x1n1)

      plot.coverage(round(HPD.int,10), ylim = ygrenzen, smooth = 0.05)


###################################################
### code chunk number 17: anteilJeffreysKI
###################################################
getOption("SweaveHooks")[["fig"]]()
      KI.int = t(sapply(x1, function(one.x){
        alpha = 1-prob
        return(qbeta(c(alpha/2, 1-alpha/2), 0.5 + one.x, 0.5 + n1 - one.x))
      })
      )
      plot.coverage(KI.int, ylim = ygrenzen, smooth = 0.05)


###################################################
### code chunk number 18: normalGamma1-par-settings
###################################################
alpha <- 2.0
beta <- 1.2
nu <- 0
lambda <- 0.5


###################################################
### code chunk number 19: normalGamma1
###################################################
getOption("SweaveHooks")[["fig"]]()
    mu.grid = seq(-7, 7, length = 100)
    kappa.grid = seq(1e-11, 6.5, length = 100)
    grid = expand.grid(mu = mu.grid, kappa = kappa.grid)

    dens = function(param){
      mu = param[1]
      kappa = param[2]
      return(dgamma(kappa, shape = alpha, rate = beta)*dnorm(mu, nu, sqrt((lambda*kappa)^-1)))
    }

    werte = matrix(apply(grid, 1, dens), nrow = length(mu.grid))

    contour(mu.grid, kappa.grid, werte, xaxs = "i", yaxs = "i", xlab = math
(mu), ylab = math (kappa), nlevels = 12, ylim=c(0,5.5), xlim=c(-5,5), labcex=1)


###################################################
### code chunk number 20: normalGamma2
###################################################
getOption("SweaveHooks")[["fig"]]()
mydgamma <- dgamma(kappa.grid[kappa.grid>0.01], alpha, beta)
      plot(kappa.grid[kappa.grid>0.01], mydgamma, xlab = math (kappa), ylab = math (f(kappa)), type = "l",  ylim=c(0, max(mydgamma)))


###################################################
### code chunk number 21: normalGamma3
###################################################
getOption("SweaveHooks")[["fig"]]()
      # for die Dichte der t-distribution, arbeitet vektorwertig.
      my.t = function(x, mu, sigma2, alpha){
        const = gamma(0.5*(alpha + 1))/(gamma(0.5*alpha)*gamma(0.5))*(sigma2*alpha)^(-1/2)
        return(const*(1+(x-mu)^2/(alpha*sigma2))^(-(alpha+1)/2))
      }
myt <- my.t(mu.grid, nu, beta / (alpha * lambda), 2 * alpha)
plot(mu.grid, myt, xlab = math (mu), ylab = math (f(mu)), type = "l", ylim=c(0, max(myt)))


###################################################################
### code chunk number 22: alco-posterior-when-both-are-unknown-mu
###################################################################
getOption("SweaveHooks")[["fig"]]()
x <- iten$tf500a1
                                                     
location <- with (alcoData, mean)
scale <- with (alcoData, sd / sqrt (n))
df <- with (alcoData, n-1)

grid <- seq (min (x), max (x), length = 10000)
densVec <- dst (grid, xi = location, omega = scale, nu = df)

plot (grid, densVec, type = "l",
      xlab = math (mu), ylab = math (f*(mu~plain("|")~x[1*plain(":")*n])))
rug(x)


#########################################################################
### code chunk number 23: alco-posterior-when-both-are-unknown-sigma2
#########################################################################
getOption("SweaveHooks")[["fig"]]()
n <- length(x)
a <- (n-1)/2
b <- sum((x - mean(x))^2)/2

grid <- seq(40000, 80000, length = 200)
plot(grid, dinvgamma(grid, a, b),
     xlab = math (sigma^2), ylab = math (f*(sigma^2~plain("|")~x[1*plain(":")*n])),
     type = "l", yaxt="n")

yticks <- (0:7) * 10^(-5)
axis(side=2, at=yticks,
     labels=
     expression(0, 
                1 %.% 10^-5,
                2 %.% 10^-5,
                3 %.% 10^-5,
                4 %.% 10^-5,
                5 %.% 10^-5,
                6 %.% 10^-5,
                7 %.% 10^-5))

funline(b/(a+1), dinvgamma, shape = a, scale = b, col="black", label="")
funline(var(x)*(n-1)/n, dinvgamma, shape = a, scale = b, col="black", label="")
funline(var(x), dinvgamma, shape = a, scale = b, col="black", label="")
funline(b/(a-1), dinvgamma, shape = a, scale = b, col="black", label="")


###################################################
### code chunk number 24: discreteAsymA
###################################################
# Remark: Due to randomness, the top panel of this figure may look a bit different
# than in Figure 6.11a) in the book
getOption("SweaveHooks")[["fig"]]()
    thetat <- 0.25
    theta <- seq(0.05:0.95, by=0.1)
    prior <- rep(.1, 10)

    n <- c(10,100,1000)
    x <- sapply(n, function(one){rbinom(1, size=one, prob=thetat)})

    posterior = t(apply(cbind(n,x), 1, function(one){dbinom(one[2], prob=theta, size=one[1])*prior}))
    posterior <- posterior/rowSums(posterior)


    par(mfrow=c(2,1))
    barplot(posterior, names.arg=theta, xlab=math (theta),
            ylab=math (f*(theta~plain("|")~x)), ylim=c(0,1), beside = TRUE)
    legend("topright", bty = "n", col = gray.colors(length(n)),
           fill = gray.colors(length(n)), title = "n = ", legend = n)

    KLdistance <- function(fx, true){

      return(sum(true*(log(true)-log(fx))))


    }
    distance <- rep(NA, 10)

    for(i in 1:10)
    distance[i] <- KLdistance(dbinom(c(0:1), prob=thetat, size=1), dbinom(c(0:1), prob=theta[i], size=1))

    plot(theta, distance, type="b", ylim=c(0, max(distance)), xaxt = "n",
         ylab="KL discrepancy", xlab=math (theta))
    axis (1, theta)


###################################################
### code chunk number 25: discreteAsymB
###################################################
# Remark: Due to randomness, the top panel of this figure may look a bit different
# than in Figure 6.11b) in the book
getOption("SweaveHooks")[["fig"]]()
    thetat <- 0.33
    theta <- seq(0.05:0.95, by=0.1)
    prior <- rep(.1, 10)

    n <- c(10,100,1000)
    x <- sapply(n, function(one){rbinom(1, size=one, prob=thetat)})

    posterior = t(apply(cbind(n,x), 1, function(one){dbinom(one[2], prob=theta, size=one[1])*prior}))
    posterior <- posterior/rowSums(posterior)


    par(mfrow=c(2,1))
    barplot(posterior, names.arg=theta, xlab=math (theta),
            ylab=math (f*(theta~plain("|")~x)), ylim=c(0,1), beside = TRUE)
    legend("topright", bty = "n", col = gray.colors(length(n)),
           fill = gray.colors(length(n)), title = "n = ", legend = n)

    KLdistance <- function(fx, true){

      return(sum(true*(log(true)-log(fx))))


    }
    distance <- rep(NA, 10)

    for(i in 1:10)
    distance[i] <- KLdistance(dbinom(c(0:1), prob=thetat, size=1), dbinom(c(0:1), prob=theta[i], size=1))

    plot(theta, distance, type="b", ylim=c(0, max(distance)), xaxt = "n",
         ylab="KL discrepancy", xlab=math (theta))
    axis (1, theta)


###################################################
### code chunk number 26: continuousAsymA
###################################################
getOption("SweaveHooks")[["fig"]]()
      # Funktion contplot for stetige Asymptotik
      contplot = function(n, x){
        par(mfrow=c(1,1))

        thetat <- 0.1

        alpha <- .5
        beta <- .5

        mle <- x/n
        fisher <- n / (mle * (1 - mle))
        
        alpha2 <-alpha+x
        beta2 <- beta+n-x



        pmean <- alpha2/(alpha2+beta2)
        pmode <- (alpha2-1)/(alpha2+beta2-2)
        pvar <- alpha2*beta2/((alpha2+beta2)^2*(alpha2+beta2+1))

        lower <- max(0.01, pmean-3*sqrt(pvar))
        upper <- min(0.99, pmean+3*sqrt(pvar))

        grid <- seq(lower, upper, by=sqrt(pvar)/600)
        posterior <- dbeta(grid, alpha2, beta2)

        target <- function(pi, ...){
          return(-dbeta(pi, alpha2, beta2, log=T))

        }


        # exakte Bestimmung der Krümmung am Maximum
        pcurv = (alpha + beta + n - 2)^3 / (alpha + x - 1) / (beta + n - x - 1)

        ## compute all three approximations
        approx2 <- dnorm(grid, mean=mle, sd=sqrt(1/fisher))
        approx3 <- dnorm(grid, mean=pmode, sd=sqrt(1/pcurv))
        approx4 <- dnorm(grid, mean=pmean, sd=sqrt(pvar))

        linecolours = c(1, 1, 1, "darkgrey")
        linetypes = c(1,2,3,2)
        linewidths = 2

        matplot(grid, cbind(posterior, approx2, approx3, approx4), type="l",
                ylab=math (f*(pi~plain("|")~x)), xlab=math (pi), 
                col = linecolours, lty=linetypes, lwd = linewidths)

        legend("topright", 
               lty = linetypes, 
               legend=
               c("posterior", "approximation 2",
                 "approximation 3", "approximation 4"), 
               col = linecolours, bty = "n", lwd = linewidths)
      }
      set.seed(1133)
      contplot(10, 1)


###################################################
### code chunk number 27: continuousAsymB
###################################################
getOption("SweaveHooks")[["fig"]]()
      set.seed(123)
      contplot(100, 10)


###################################################
### code chunk number 28: continuousAsymC
###################################################
getOption("SweaveHooks")[["fig"]]()
      set.seed(123)
      contplot(1000, 100)


###################################################
### code chunk number 29: continuousAsymD
###################################################
getOption("SweaveHooks")[["fig"]]()
      set.seed(123)
      contplot(10000, 1000)


###################################################
### code chunk number 30: PoissonEmpirischerBayes
###################################################
  ## loglikelihood without normalizing constants
  poissonGammaLogLik <- function (theta, x, e)
  {
    n <- length (x)
    alpha <- theta[1]
    beta <- theta[2]

    n * (alpha * log (beta) - lgamma (alpha)) +
    sum (lgamma (alpha + x) - (alpha + x) * log (beta + e))
  }
  scotlandMargMle <- with (scotlandData,
  optim (c (1, 1), poissonGammaLogLik, x = x, e = e,
  method = "L-BFGS-B", lower = rep (eps, 2),
  control = list (fnscale = -1))
  )
alphaMl <- scotlandMargMle[[1]][1]
betaMl <- scotlandMargMle[[1]][2]


##############################################################
### code chunk number 31: PoissonEmpirischerBayes_Abbildung
##############################################################
# Hint: run chap4.R before this chunk!

getOption("SweaveHooks")[["fig"]]()
posteriorPars <- with (scotlandData,
                       cbind (alphaMl + x, betaMl + e)
                       )
ci <- cbind(posteriorPars[,1] / posteriorPars[,2],
            qgamma (0.025, posteriorPars[,1], posteriorPars[,2]),
            qgamma (0.975, posteriorPars[,1], posteriorPars[,2])
            )

## order in the same way as in chap4.R
load("../data/scotland/scotOrder.RData")
ci2 <- ci[scotOrder,]

signif <- ci2[,2] > 1 | ci2[,3] < 1

plotCI(x=ci2[,1], li=ci2[,2], ui=ci2[,3],
       err="y", ylab="Relative risk", xlab="Region",
       col="black", # ifelse(signif==0, "grey", "black"),
       gap=0.25)

with (scotlandData,
      points (1:56, (x / e)[scotOrder], pch = 19, col = "black"))
      # ifelse(signif==0, "grey", "black"))
      

abline (h = 1, lty=2)
abline (h = alphaMl / betaMl, lty=3)


###################################################
### code chunk number 32: preeclampsiaMetaLogOR
###################################################

diuretics <- with (preeclampsia,
                   c (sum (Diuretic[Preeclampsia == "yes"]), sum (Diuretic[Preeclampsia == "no"]))
                   )
controls <- with (preeclampsia,
                   c (sum (Control[Preeclampsia == "yes"]), sum (Control[Preeclampsia == "no"]))
                   )
logOddsRatio <- log (diuretics[1] * controls[2] / diuretics[2] / controls[1])
standardError <- sqrt (sum (1/diuretics) + sum (1/controls))
waldKiLogOddsRatio <- logOddsRatio + c (-1, 1) * 1.96 * standardError
waldData <- list (logOddsRatio, waldKiLogOddsRatio)
save (waldData, file = "../data/waldKiEklampsie.RData")

## compute data
oddsRatio <- function (square)
    (square[1,1] * square[2,2]) / (square[1,2] * square[2,1])

variance <- function (square)
    sum (1 / square)

groups <- split (subset (preeclampsia, select = c (Diuretic, Control)),
                 preeclampsia[["Trial"]])

logOddsRatios <- log (sapply (groups, oddsRatio))
variances <- sapply (groups, variance)

logOddsRatio2 <- sum(logOddsRatios/variances)/sum(1/variances)
standardError2 <- sqrt (1/sum(1/variances))
waldKiLogOddsRatio2 <- logOddsRatio2 + c (-1, 1) * 1.96 * standardError2
waldData2 <- list (logOddsRatio2, waldKiLogOddsRatio2)


###################################################
### code chunk number 33: preeclampsiaEmpricialBayes
###################################################
getOption("SweaveHooks")[["fig"]]()

## study effects
postSd <- sqrt (variances)
postLower <- logOddsRatios - postSd * 1.96
postUpper <- logOddsRatios + postSd * 1.96

## plot intervals with point estimates

panel.ci <- function(x, y, lx, ux, subscripts, horLine=NULL, ...)
{
    x <- as.numeric(x)
    y <- as.numeric(y)

    lx <- as.numeric(lx[subscripts])
    ux <- as.numeric(ux[subscripts])

    panel.dotplot(x, y, lty = 2, ...)            # normal dotplot for point estimates
    panel.arrows(lx, y, ux, y,          # draw intervals
                 length = 0.1, unit = "native",
                 angle = 90,            # deviation from line
                 code = 3,              # left and right whisker
                 ...)
    panel.abline (v = 0, lty = 2)       # reference line
    
    if(! is.null(horLine))
    {
        panel.abline(h=horLine, lty=1)
    }
}

studyNames <- c (names (postLower))
studyNames <- ordered (studyNames, levels = rev (studyNames)) # levels important for order!

ciData <- data.frame (low = postLower, 
                      up = postUpper,
                      mid = logOddsRatios,
                      names = studyNames
                      )
ciData[["signif"]] <- with (ciData,
                            up < 0 | low > 0)
ciData[["farbe"]] <- "black"

randomEffectsCiPlot <- with (ciData,
                             dotplot (names ~ mid,
                                      panel = panel.ci,
                                      lx = low, ux = up,
                                      pch = 19, col = farbe,
                                      xlim = c (-3, 3), xlab = "log odds ratio",
                                      scales = list (cex = 1)
                                      )
                             )
print (randomEffectsCiPlot)


###################################################
### code chunk number 34
###################################################
getOption("SweaveHooks")[["fig"]]()
## compute MLE
nuMle <- function (tau2)
{
    precisions <- 1 / (variances + tau2)
    weighted.mean (logOddsRatios, precisions)
}

profilLogLikTau2 <- function (tau2)
{
    margVars <- variances + tau2
    arg <- log (margVars) + (logOddsRatios - nuMle (tau2))^2 / margVars
    - 1/2 * sum(arg)
}


tau2Ml <- optimize (profilLogLikTau2, c (eps, 1/eps), maximum = TRUE)$maximum
nuMl <- nuMle (tau2Ml)

## empirical Bayes
weights <- variances / (variances + tau2Ml)

postExpectation <- weights * nuMl + (1 - weights) * logOddsRatios
postSd <- sqrt ((1 - weights) * variances)
postLower <- postExpectation - postSd * 1.96
postUpper <- postExpectation + postSd * 1.96

## profile log lik ci for nu
tau2Mle <- function (nu)
{
    scoreTau2 <- function (tau2)
    {
        margVars <- variances + tau2
        arg <- (logOddsRatios - nu)^2 / margVars - 1
        arg <- arg / margVars
        1/2 * sum (arg)
    }
    res <- uniroot (scoreTau2, c (eps, 1/eps))
    return(res$root)
}

profilLogLikNu <- function (nu)
{
    margVars <- variances + tau2Mle (nu)
    arg <- log (margVars) + (logOddsRatios - nu)^2 / margVars
    - 1/2 * sum (arg)
}

normProfilLogLikNu <- function (nu)
    profilLogLikNu (nu) - profilLogLikNu (nuMl)

likCiNu <- likelihood.ci (alpha = 0.05, normProfilLogLikNu, nuMl, lower = -10, upper = +10)

## data from previous fixed effect example
load ("../data/waldKiEklampsie.RData")    # loads list "waldData"

## plot intervals with point estimates

studyNames <- c (names (postLower), "Random effects", "Fixed effect")
studyNames <- ordered (studyNames, levels = rev (studyNames)) # levels important for order!

ciData <- data.frame (low = c (postLower, likCiNu[1], waldData2[[2]][1]),
                      up = c (postUpper, likCiNu[2], waldData2[[2]][2]),
                      mid = c (postExpectation, nuMl, waldData2[[1]]),
                      names = studyNames
                      )
ciData[["signif"]] <- with (ciData,
                            up < 0 | low > 0)
ciData[["farbe"]] <- "black"

randomEffectsCiPlot <- with (ciData,
                             dotplot (names ~ mid,
                                      panel = panel.ci,
                                      lx = low, ux = up,
                                      pch = 19, col = farbe,
                                      xlim = c (-1.75, 1.75), xlab = "log odds ratio",
                                      scales = list (cex = 1),
                                      horLine=2.5
                                      )
                             )
print (randomEffectsCiPlot)


