#################################################################
### Settings for the other R files
#################################################################

# Remark: Please run this file before running any other R file.


#################################################
### load required packages
#################################################

# for tables etc
library(Hmisc)
library(xtable)

# for Huber estimator, truehist
library(MASS)

# for 2-dim. kernel density estimator
library(KernSmooth)

## for dst
library(sn)

## for nice confidence interval plots
library(gplots)
library(cubature)
library(lattice)

## for the noncentral hypergeometric distribution,
## and rdirichlet
library(MCMCpack)

## for summaryBy
library(doBy)

## for survival analysis
library(survival)

## for formatPval
tmp <- library("biostatUZH", logical.return = TRUE)
if(! tmp)
{
  install.packages("biostatUZH", repos = "http://R-Forge.R-project.org")   
  # Mac OS X users should add:  type="source"
  library("biostatUZH")
}


#################################################
### data preparation
#################################################


## colon cancer ---------------------------------------------------------------------------------------

colonCancer <- c(NA, 37, 22, 25, 29, 34, 49)


## PBC survival times ----------------------------------------------------------------------------

pbcTreat <- read.table ("../data/pbcTreat.txt", header = TRUE)

pbcTreatSurv <- Surv(pbcTreat$time,pbcTreat$d) # d ist death-Variable, d.h. 1 heißt
# beobachtet, 0 zensiert.
pbcTreat.nObs <- nrow(pbcTreat)

uncensored <- subset (pbcTreat, d == 1)   # tatsächliche Üzeiten in Treatmentgruppe
unc <- uncensored$time
uncN <- length (unc)
uncMean <- mean (unc)
uncSd <- sd (unc)
uncSum <- sum (unc)


## lip cancer in Scotland -----------------------------------------------------------------------

scotlandData <- read.table ("../data/scotland/scotland.txt", header = TRUE)

scotlandCoordinates <-  matrix (scan ("../data/scotland/scotkoord.txt"),
                                ncol = 3, byrow = T)


## preeclampsia ------------------------------------------------------------------------------------

preeclampsia <- read.table ("../data/preeclampsia.txt", header = TRUE)


## blood alcohol content --------------------------------------------------

iten <- read.table("../data/tbsTotal.csv", header=TRUE)
iten$tf500a1 <- iten$BAK / iten$L_500a_1 * 2100



#################################################
### general settings and functions
#################################################


## general options -----------------------------------------------------------------------------
options(encoding = "utf8",            # determine encoding
        SweaveHooks =
        list(                           # every graphics is formatted as follows:
             fig = function(){
                 par(las = 1, # always horizontal axis labels
                     cex.lab = 1.35, # increase font height for labels
                     mar = c(4, 4.5, 0.2, 0.5)+0.1, # reduce top and right margins
                     xaxs="r", yaxs="r", # regular axis styles
                     font.lab = 1       # label font unchanged
                     )}
             ),
        myColour = FALSE,               # no colour in graphics
        width=60,                       # only 60 characters in R output
        continue="  ",                  # use two blanks instead of "+" for
                                        # continuation
        SweaveSyntax="SweaveSyntaxNoweb"
        )

## temporary file for LaTeX output:
latexTempFile <- tempfile ()

## fonts -------------------------------------------------------------------------------------------
CM <- Type1Font("CM",
                file.path("Schriften", "metrics",
                          c("fcsr8a.afm", "fcsb8a.afm", "fcmri8a.afm", "fcmbi8a.afm", "cmsyase.afm"))
                )
pdfFonts(CM = CM)
postscriptFonts(CM = CM)

ps.options (family = "CM")

## functions  -------------------------------------------------------------------------------------


## drop the last num elements 
allButLast <- function (vector, num = 1)
    vector[-(length (vector) - (0:(num-1)))]

## half-normal distribution, i.e distribution of |Normal(mean, sd^2)| : _hnorm

rhnorm <- function (n, mean = 0, sd = 1)
{
    abs (rnorm (n = n, mean = mean, sd = sd))
}

dhnorm <- function (x, mean = 0, sd = 1, log = FALSE)
{
    val <- dnorm (x, mean = mean, sd = sd, log = FALSE) +
        dnorm (-x, mean = mean, sd = sd, log = FALSE)

    if (log)
        return (log (val))
    else
        return (val)
}

phnorm <- function (q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
{
    val <- ifelse (q > 0,
                   pnorm (q, mean = mean, sd = sd) - pnorm (-q, mean = mean, sd = sd),
                   0)
    if (!lower.tail)
        val <- 1 - val
    if (log.p)
        return (log (val))
    else
        return (val)
}

qhnorm <- function (p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
{
    if (log.p)
        p <- exp (p)

    target <- function (q)
        phnorm (q, mean = mean, sd = sd, lower.tail = lower.tail, log.p = FALSE) - p

    unirootObj <- uniroot (target, lower = 0, upper = 1/sqrt (.Machine$double.eps))
    unirootObj$root
}

meanHnorm <- function (mean = 0, sd = 1)
{
    tmp <- mean / sd

    val <- sqrt (2) * sd * exp (-1/2 * tmp^2) / sqrt (pi) +
        mean * (2 * pnorm (tmp) - 1)

    val
}

varHnorm <- function (mean = 0, sd = 1)
{
    val <- sd^2 + mean^2 - meanHnorm (mean = mean, sd = sd)^2
    val
}

modeHnorm <- function (mean = 0, sd = 1)
{
    if (abs (mean) <= sd){
        0
    } else {
        gradient <- function (x, mean, sd){
            val <- dnorm ((x - mean) / sd) * (mean - x) -
                dnorm ((x + mean) / sd) * (mean + x)
            val / sd^3
        }
        eps <- sqrt (.Machine$double.eps)
        optimObj <- optim (mean/2, fn = dhnorm, gr = gradient,
                           mean = mean, sd = sd,
                           method = "L-BFGS-B", lower = eps, upper = mean + eps,
                           control = list (fnscale = -1))
        stopifnot (optimObj$convergence == 0)
        optimObj$par
    }
}

## this function returns all values pi, for which the beta density
## with parameters p1 and p2 is smaller than h
outerdens= function(h, p1, p2){
  modus = (p1-1)/(p1 + p2 - 2)
  if(x==0)
    schnitt.l <- 0
  else
    schnitt.l = uniroot(function(x){dbeta(x, p1, p2) - h}, interval = c(0, modus))$root
  if(x==n)
    schnitt.u <- 1
  else
    schnitt.u = uniroot(function(x){dbeta(x, p1, p2) - h}, interval = c(modus, 1))$root
  return(c(pbeta(schnitt.l, p1, p2) + pbeta(schnitt.u, p1, p2, lower.tail = FALSE), schnitt.l, schnitt.u))
}

## abbreviation for mathematical expression in font.lab font
## multiple arguments are possible
math <- function (..., font = 3)
{
    ## unpack dot-dot-dot argument
    dots <- as.list (substitute (...[]))
    dots <- dots[2:(length (dots) - 1)]

    ## font extension
    font <- switch (font, "plain", "bold", "italic", "bolditalic")
    dots <- paste (font, "(", dots, ")")

    ## parse and build expression
    do.call ("expression", as.list (parse (text = dots)))
}


## add MLE in figure
drawml <- function(x,                   # MLE
                   y,                   # value of L(MLE)
                   dezi=0,              # show at least dezi digits
                   digi=3,              # how many digits
                   down=F)              # should the label be shifted downwards?
{
    segments(x, par("usr")[3], x, y, lty="dashed", col = "black")
    if (down) { axis(side=1, at=x, font=2, labels=format(x, nsmall=dezi, digits=digi), line=1, tick=F) }
    else { axis(side=1, at=x, font=2, labels=format(x, nsmall=dezi, digits=digi)) }
}

## draws a vertical (gray if col==NULL) dashed line from x to fun(x)
funline <- function(x,                  # abscissa
                    fun,                # function
                    col = "gray",       # color of line
                    label = NULL,       # character label at x
                    dezi = 0,           # if label = NULL, x is added (see drawml)
                    digi = 3,           # dito (see drawml)
                    down = FALSE,       # dito (see drawml)
                    ...                 # additional arguments for fun
                    )
{
    segments(x, par("usr")[3], x, fun(x, ...), col = col, lty = 2)
    label <- ifelse (is.null (label), format(x, nsmall=dezi, digits=digi), label)
    axis(side = 1, at = x, labels = rep (label, length.out = length (x)),
         font = 2, padj = ifelse(down, 1, 0))
}

## plot the coverage probabilities for the binomial distribution
plot.coverage <- function(confidence.intervals, # matrix with 2 columns, which contains
                                                # the lower and upper limits of the CI, and n+1 rows
                          add = FALSE,  # should the plot be added?
                          smooth = NULL, # if not NULL, a local mean of the coverage probs
                                         # is added with the epsilon parameter smooth
                          faktor = 20,  # determines the number of p-values for which the coverage
                                        # is computed as n*faktor -1
                          prob = 0.95,  # nominal confidence level
                          ...
                          )
{
  # avoid too high precision (Wilson)
  confidence.intervals = round(confidence.intervals, 12)
  n = nrow(confidence.intervals)-1
  x = 0:n
  # which column contains the lower limits?
  lower.ind = which.min(confidence.intervals[1,])
  upper.ind = setdiff(c(1,2), lower.ind)
  # compute pvector: without 0 and 1
  pvector = seq(0, 1, length = faktor*n + 1)[-c(1, faktor*n+1)]
  gerade = !(length(pvector) %% 2)
  pvector.halb = pvector[1:ceiling(length(pvector)/2)]
  # compute coverage probabilities:
  # just for one half due to symmetry around p = 0.5
  coverage.probs = sapply(pvector.halb,
    function(p){
      sum(dbinom(x, n, p)[which(p >= confidence.intervals[,lower.ind] & p <= confidence.intervals[,upper.ind])])
    }
    )
  if(gerade)
    coverage.probs = c(coverage.probs, rev(coverage.probs))
  else
    coverage.probs = c(coverage.probs, rev(coverage.probs)[-1])
  # plot
  if(!add){
    plot(pvector, coverage.probs, type = "l", xlab = expression(paste("True ", pi)),
         ylab = "Coverage probability", col = "darkgrey", ...)
    abline(h = prob, lty = 2)
  }
  else
    {
      lines(pvector, coverage.probs, type = "l", col = "darkgrey", ...)
  }
  # draw smooth
  if(!is.null(smooth)){
    # general a function
    a = function(p){
      if(p==0) NA
      else if (p==1) NA
      else
      if(p <= smooth)
        NA#(p*(1-p)*p^(-2) - 1)*p#1 - 2*smooth
      else if(p >= 1 - smooth)
        NA#(p*(1-p)*(1-p)^(-2) - 1)*p #1/smooth - 3 + 2*smooth
      else
        (p*(1-p)*smooth^(-2) - 1)*p
    }
    # function to compute the local coverage
    local.coverage = function(p){
      ap = a(p)
      a1mp = a(1-p)
      alpha = ap + x
      beta = a1mp + n - x
      values.gamma = (lchoose(n, x)
                    + lgamma(ap + a1mp) - lgamma(ap) - lgamma(a1mp)
                    + lgamma(ap + x) + lgamma(a1mp + n - x) - lgamma(ap + a1mp + n)
                      )
      values.integral = log(
        pbeta(confidence.intervals[,upper.ind], alpha, beta) - pbeta(confidence.intervals[,lower.ind], alpha, beta)
        )
      sum(exp(values.gamma + values.integral))
    }

    # compute and draw:
    coverage.average = sapply(pvector.halb, local.coverage)

    if(gerade)
      coverage.average = c(coverage.average, rev(coverage.average))
    else
      coverage.average = c(coverage.average, rev(coverage.average)[-1])

    lines(pvector, coverage.average, lwd = 2)
  }
}


## Function to compute likelihood based confidence interval, basically
## the two solutions to
##            f(\theta) = l(\theta)-l(\hat{theta)) + 1/2 dchisq(1-alpha,df=1)=0
## are found. (by Michael Höhle)
likelihood.ci <- function(alpha=0.05,   # confidence level (see Equation 2.6 in Pawitan (2003))
                          loglik,       # Loglikelihoodfunktion
                          theta.hat,    # the MLE
                          lower,        # search interval [lower,theta.hat] for f=0
                          upper,        # search interval [theta.hat,upper] for f=0
                          comp.lower = TRUE, # should the lower boundary be computed?
                          comp.upper = TRUE, # should the upper boundary be computed?
                          ...
                          )
{
  # Highest Likelihood intervall -- target function
  f <- function(theta,...) {
    loglik(theta,...) - loglik(theta.hat,...) + 1/2*qchisq(1-alpha, df=1)
  }

  # Compute upper and lower boundary numerically, if desired only one boundary
  if(comp.lower) hl.lower <- uniroot(f,interval=c(lower,theta.hat),...)$root
  if(comp.upper) hl.upper <- uniroot(f,interval=c(theta.hat,upper),...)$root

  ret <- NULL
  if(comp.lower) ret <- c (ret, hl.lower)
  if(comp.upper) ret <- c (ret, hl.upper)

  return (ret)
}

## Compute Clopper-Pearson interval for binomial success probability.
clopperPearson <- function (x,          # number of successes
                            n,          # total number of runs
                            level = 0.95# wished (minimum) confidence level
                            )
{
    alpha <- 1 - level

    lower <- qbeta (alpha / 2, x, n - x + 1)
    upper <- qbeta (1 - alpha / 2, x + 1, n - x)

    return (cbind (lower, upper))
}


## post-processing of LaTeX output 
postLatex <- function (latexOutput,     # latex function output
                       widthFactor = 0.7,# tabular will have width (widthFactor * \textwidth)
                       minipage = TRUE,  # use minipage environment?
                       dropColTitles = FALSE # drop column titles?
                       )
{
    ## get string, skip comment lines
    x <- scan (file = latexOutput$file, what = character (0), sep = '\n', skip = 2, quiet = TRUE)

    ## optionally drop the column titles and the \midrule command.
    ## Both follow the \toprule command.
    if(dropColTitles)
    {
        colTitleIndex <- grep("toprule", x, fixed=TRUE) + 1
        x <- x[- c(colTitleIndex, colTitleIndex + 1)]
    }
    
    ## wrap around minipage
    if (minipage){
        x[1] <- paste (x[1], " \n ", "\\begin{minipage}{", widthFactor, "\\textwidth}", sep = "")
        x[length (x)] <- paste ("\\end{minipage} \n", x[length (x)])
    }

    ## replacements
    repl <- matrix (ncol = 2, byrow = TRUE, dimnames = list (NULL, c ("pattern", "replacement")),
                    data = c (
                                        # begin tabularx environment
                    "\\\\begin{tabular}",
                    paste ("\\\\begin{tabularx}{", ifelse (minipage, "", widthFactor), "\\\\textwidth}"),
                                        # end tabularx environment
                    "\\\\end{tabular}",
                    "\\\\end{tabularx}",
                                        # rule after column headings
                    "\\\\midrule",
                    "\\\\headendrule",
                                        # rule between column heading lines
                    "\\\\cline{.*}",
                    "\\\\midrule",
                                        # superfluous linebreak
                    "\\\\\\\\[ \t]*\\\\\\\\",
                    "\\\\\\\\"
                    ))

    ## process string
    for (i in seq_len (nrow (repl))){
        x <- gsub (repl[i, "pattern"],
                  repl[i, "replacement"],
                  x, perl = TRUE)
    }

    ## output
    cat(x, sep = '\n' )
}

## number formatting: add space after multiples of 10^3
formatBig <- function (x)
{
    mark <- ifelse (log10 (x) < 4, "", "\\\\\\\\,")
    format (x, big.mark = mark, scientific=FALSE)
}

## prints the number in scientific format
formatSc <- function(x, digits=3, backslash="\\\\")
{
    ## check which are exactly zero
    isZero <- x == 0.0

    ## start processing
    x <- signif(x, digits)
    exponent <- floor(log10(abs(x)))
    mantissa <- x / (10^exponent)

    ## format the mantissa
    mantissa <- format(mantissa, digits = digits)
    mantissa <-
        ifelse(mantissa == "1",
               "",
               paste(mantissa, backslash, "cdot ",
                     sep=""))

    result <- paste(mantissa,
                    "10^{", exponent,"}",
                    sep="")

    ## set to zero which were exactly zero
    result[isZero] <- "0"

    ## and return the result
    return(result)
}

## for rounding, such that formatRound(123.00, 1) is printed as 123.0 and not
## just 123.
formatRound <- function(x, digits=0,...)
{
   format(round(x, digits=digits),nsmall=digits, ...)
}


## convert matrix/data frame to LaTeX output
latexMatrix <- function (matrix)
    write.table (matrix, sep = " & ", eol="\\\\\n", quote=FALSE, col.names=FALSE, row.names = FALSE)

## empirical HPD interval
empiricalHpd <- function (theta,        # sample vector of parameter
                          level         # credibility level
                          )
{
    M <- length (theta)
    thetaorder <- sort.int (theta, method = "quick")

    alpha <- 1 - level
    maxSize <- ceiling (M * alpha)
    ind <- seq_len (maxSize)

    lower <- thetaorder[ind]
    upper <- thetaorder[M - maxSize + ind]
    size <- upper - lower
    sizeMin <- which.min(size)

    HPDLower <- thetaorder[sizeMin]
    HPDUpper <- thetaorder[M - maxSize + sizeMin]
    return (c (HPDLower, HPDUpper))
}


