#################################################################
### R code for Chapter 3
#################################################################

### Encoding: UTF-8

###################################################
### code chunk number 1
###################################################
## weitere Berechnungen in settings
uncQ <- qgamma (c (0.025, 0.975), uncN, 1)
lambdaKI <- uncQ / uncSum
ewertKI <- rev (1 / lambdaKI)


###################################################
### code chunk number 2: alcohol-exact
###################################################
data <- iten$tf500a1
lower <- mean(data) - qt(0.975, df=length(data)-1)*sd(data)/sqrt(length(data))
upper <- mean(data) + qt(0.975, df=length(data)-1)*sd(data)/sqrt(length(data))

lower.sd <- sqrt((length(data) - 1) * var(data) / qchisq(0.975, df=length(data)-1))
upper.sd <- sqrt((length(data) - 1) * var(data) / qchisq(0.025, df=length(data)-1))


###################################################
### code chunk number 3
###################################################
pmterm <- 1.96*uncSd / sqrt (uncN)
apprKI <- uncMean + c(-1,1)*pmterm


###################################################
### code chunk number 4
###################################################
n <- 100
x <- 2
pihat <- x/n
se.pihat <- sqrt(pihat*(1-pihat)/n)
pi.lower <- pihat - 1.96*se.pihat
pi.upper <- pihat + 1.96*se.pihat

logit <- function(x) log(x/(1-x))
expit <- function(x) exp(x)/(1+exp(x))
vphat <- logit(pihat)
se.vphat <- sqrt(1/x+1/(n-x))
vp.lower <- vphat - 1.96*se.vphat
vp.upper <- vphat + 1.96*se.vphat


###################################################
### code chunk number 5
###################################################
data <- iten$tf500a1

set.seed(090503)
nboot <- 10000
boot.data <- matrix(NA, nrow=nboot, ncol=length(data))
for(i in 1:nboot) boot.data[i,] <- sample(data, size=length(data),replace=TRUE)
boot.mean <- apply(boot.data, 1, mean)
boot.median <- apply(boot.data, 1, median)
cv <- function(data) sd(data)/mean(data)
boot.cv <- apply(boot.data, 1, cv)

boot.ci.mean <- quantile(boot.mean, c(0.025, 0.975))
boot.ci.median <- quantile(boot.median, c(0.025, 0.975))
boot.ci.cv <- quantile(boot.cv, c(0.025, 0.975))

lower <- mean(data) - qt(0.975, df=length(data)-1)*sd(data)/sqrt(length(data))
upper <- mean(data) + qt(0.975, df=length(data)-1)*sd(data)/sqrt(length(data))


###################################################
### code chunk number 6: boot-hist-means
###################################################
getOption("SweaveHooks")[["fig"]]()
truehist(boot.mean, 
         xlab=math(hat(mu)),
         prob=TRUE,
         col="gray",
         nbins=20)
abline(v=
       c(mean(boot.mean),
         mean(data)),
       lty=c(1, 2))


###################################################
### code chunk number 7: boot-hist-coef-variation
###################################################
getOption("SweaveHooks")[["fig"]]()
truehist(boot.cv, 
         xlab=math(hat(phi)),
         prob=TRUE,
         col="gray",
         nbins=20)
abline(v=
       c(mean(boot.cv),
         cv(data)),
       lty=c(1, 2))


###################################################
### code chunk number 8
###################################################
tmp <- library("MBESS", logical.return = TRUE)
if(! tmp)
{
  install.packages("MBESS", dependencies = "Depends")   
  # Mac OS X users should add:  type="source"
  library("MBESS")
}
res <- ci.cv(data=data, conf.level=.95)
res$Lower.Limit.CofV
res$Upper.Limit.CofV


###################################################
### code chunk number 9: boot-bca-example
###################################################
bca <- function(bootsamples,
                gamma,
                estimate,
                a,
                c)
{
    qfun <- function(alpha)
    {
        c + (c + qnorm(alpha)) / (1 - a * (c + qnorm(alpha)))
    }
    
    return(c(quantile(bootsamples, 
                      probs=
                      pnorm(c(qfun((1 - gamma) / 2),
                              qfun((1 + gamma) / 2))))))
}

jn.means <- sapply(seq_along(data),
                   function(i){mean(data[-i])})
a.means <- 1/6 * sum((mean(jn.means) - jn.means)^3) / 
    (sum((mean(jn.means) - jn.means)^2))^(3/2)
c.means <- qnorm(mean(boot.mean <= mean(data)))
bca.means <- bca(boot.mean,
                 gamma=0.95,
                 estimate=mean(data),
                 a=a.means,
                 c=c.means)

jn.cv <- sapply(seq_along(data),
                function(i){cv(data[-i])})
a.cv <- 1/6 * sum((mean(jn.cv) - jn.cv)^3) / 
    (sum((mean(jn.cv) - jn.cv)^2))^(3/2)
c.cv <- qnorm(mean(boot.cv <= cv(data)))
bca.cv <- bca(boot.cv,
              gamma=0.95,
              estimate=cv(data),
              a=a.cv,
              c=c.cv)

bca.cv.quantiles <- 
    as.numeric(gsub(pattern="%", replacement="", x=names(bca.cv), fixed=TRUE)) 


###################################################
### code chunk number 10: pvalue-evidence
###################################################
getOption("SweaveHooks")[["fig"]]()
n <- 100
N <- 100
cor.range <- seq(-1, 1, length.out = n)
mat <- matrix(rep(cor.range, N), ncol = N, byrow = TRUE)
cols <- colorRampPalette(c("grey20", "white"))
                                      
names <- rev(c(1, 0.1, 0.01, 0.001, format(0.0001, scientific = FALSE)))

y <- seq(0, 1, length = 5)
posy <- y[-5] + 0.125
x <- 0.1

par(mar = c(1, 4, 1, 1))
image(z = mat, col = cols(n), axes = FALSE, bty = "n", ylab = math(p), main = "")
axis(2, at = y, labels = names, las = 1, lwd = 0, lwd.ticks = 1)

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

text(x, posy[4], "no evidence", col = "grey50", cex = 1.5, font = 2, adj = 0)
text(x, posy[3], "weak evidence", col = "white", cex = 1.5, font = 2, adj = 0)
text(x, posy[2], "substantial evidence", col = "white", cex = 1.5, font = 2, adj = 0)
text(x, posy[1], "strong evidence", col = "white", cex = 1.5, font = 2, adj = 0)


###################################################
### code chunk number 11
###################################################
tervila <- matrix(c(6, 2, 102, 101), nrow = 2)
f.test <- fisher.test(tervila, alternative = "greater")


###################################################
### code chunk number 12: help.chunk0
###################################################
t0 <- exp(lchoose(108,0)+lchoose(103,8)-lchoose(211,8))
t1 <- exp(lchoose(108,1)+lchoose(103,7)-lchoose(211,8))
t2 <- exp(lchoose(108,2)+lchoose(103,6)-lchoose(211,8))
t3 <- exp(lchoose(108,3)+lchoose(103,5)-lchoose(211,8))

t6 <- exp(lchoose(108,6)+lchoose(103,2)-lchoose(211,8))
t7 <- exp(lchoose(108,7)+lchoose(103,1)-lchoose(211,8))
t8 <- exp(lchoose(108,8)+lchoose(103,0)-lchoose(211,8))


###################################################
### code chunk number 13: help.chunk
###################################################
t <- uncSum/1000
n <- uncN


###################################################
### code chunk number 14: gamma-pivot-p-value
###################################################
t
n
pgamma(t, shape=n, rate=1, lower.tail=FALSE)


###################################################
### code chunk number 15
###################################################
S0 <- sum(unc-1000)^2/(uncN)
T0 <- (1130.8-1000)/(sqrt(S0)/sqrt(47))
P0 <- pnorm(T0,lower.tail=FALSE)


