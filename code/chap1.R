#################################################################
### R code for Chapter 1
#################################################################

### Encoding: UTF-8

###################################################
### code chunk number 1: deFinetti
###################################################
getOption("SweaveHooks")[["fig"]]()
################################## Zeichnen von de Finetti-Diagramm ################################

## Setup vom Diagramm

s = 2/sqrt(3)
aa = c(0, 0)
AA = c(s, 0)
Aa = c(s/2, 1)

## grays = gray(c(0.7, 0, 0.4))
grays = gray(rep(0, 3))
names(grays) = Cs(aa, AA, Aa)
mylty = 2

## returns: nichts
finetti.setup = function(){
  triangle = cbind(aa, AA, Aa)

  ## Dreieck (aspect ratio = 1, um Gleichseitigkeit zu sehen)
  outerspace = 0.1
  plot(triangle[1,c(1,2,3,1)], triangle[2,c(1,2,3,1)], type = "l",
       main = "", xlab = "", ylab = "", xlim = c(-outerspace, s + outerspace), ylim=c(-0.1,1.25),
       lwd = 1, asp = 1)
                                        
  ## Beschriftung
  text(triangle[1,], triangle[2,], labels = Cs(aa, AA, Aa), col = grays, cex = 1.3, pos = c(2, 4, 3))
  ## HW-Parabel
  pvalues = seq(0, s, length = 100)
  heterozygote.fraction = 2*pvalues/s*(1-pvalues/s)
  lines(pvalues, heterozygote.fraction, lty=5, lwd=2)
  legend("topright", lty=5, legend="Hardy-Weinberg equilibrium", lwd=2, bty="n")
}

## Punkte einzeichnen
## a, b, sind die Anteile der Genotypen aa \bzw AA in der Population aus (0,1),
## labels kann optional ein Character-Vektor der Länge 4 sein mit Beschriftungen der p, a, b, c-Werte, sonst
## werden die Anteilswerte hingeschrieben.
## name: Name des Punktes ## returns: nichts
finetti.point = function(a, b, labels = NULL, name = NULL){                                         # Punkt im Koordinatensystem bestimmen
  c = 1-a-b
  pvalue = (2*a + c)/2*s
  punkt = c(pvalue, c)
                                        # Allelhäufigkeit von A ist p
  p = pvalue/s
                                        # Fußpunkte von a, b und c bestimmen
  ## point: ein Vektor der Länge n
  ## line: Matrix mit zwei Spalten, deren Vektoren der Länge n die Gerade bestimmen, auf
  ## die der point projeziert wird
  ## returns: gibt die Koordinaten des Lotfußpunktes zurück.
  projection = function(point, line){
    direction = line[,2] - line[,1]
    foot = (direction %*% (point-line[,1]))/sum(direction^2) * direction + line[,1]
    return(foot)
  }
  a.foot = projection(punkt, cbind(Aa, aa))
  b.foot = projection(punkt, cbind(Aa, AA))
  c.foot = c(pvalue, 0)
  foots = cbind(a.foot, b.foot, c.foot)
  ## Einzeichnen
                                        # Lote
  segments(rep(punkt[1], 3), rep(punkt[2], 3),
           foots[1,], foots[2,],
           col = grays[Cs(AA, aa, Aa)],
           lty = mylty
           )
                                        # Punkt selber
  points(pvalue, c, pch = 19)
  ## Beschriftung
  if(is.null(labels)) labels = formatRound(c(p, a, b, c), 3)
  middles = cbind(c.foot, (foots + punkt)/2)
  text(middles[1,],
       middles[2,],
       labels = labels,
       cex = 1, col = c(1, grays[Cs(AA, aa, Aa)]), pos = c(1, 2, 4, 4)
       )
  if(!is.null(name)) text(punkt[1], punkt[2], labels = name, pos = 3) }
## zeichnet Punkt im HWE ein mit A-Häuf. p
finetti.hwpoint = function(p, ...){
  finetti.point(p^2, (1-p)^2, ...)
}

finetti.setup()
a.hwu = 0.3
b.hwu = 0.6
finetti.point(a.hwu, b.hwu, labels = Cs(q, a, b, c), name = Cs(F)) # extra Fußpunkt von c benennen und markieren:
s = 2/sqrt(3)
c.hwu = 1-a.hwu-b.hwu
pvalue = (2*a.hwu + c.hwu)/2*s
points(pvalue, 0, pch = 19)
text(pvalue-0.01, 0.01, labels = Cs(Q), adj = c(1,0))
par (mar = c (1, 1.5, 0, 0))
finetti.hwpoint(0.5, name = Cs(G))


###################################################
### code chunk number 2: Lippenkrebs
###################################################

############################################################################
##                                                                        ##
##       Schottlandlandkarte grau, 1000 Farbstufen, Log-Skala             ##
##                                                                        ##
############################################################################

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


data <- scotlandData
sir <- data[,1]/data[,2]


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


###################################################
### code chunk number 3: LippenkrebsPlot
###################################################
getOption("SweaveHooks")[["fig"]]()
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


###################################################
### code chunk number 4: Lippenkrebs2
###################################################
getOption("SweaveHooks")[["fig"]]()


plot(sqrt(data[,2]), sqrt(sir),xlab="Expected number",ylab="SIR",axes=F, xlim=c(0, 10), ylim=c(0,sqrt(7)), pch=19)
lines(c(0,sqrt(100)), rep(1,2), lty=2)

x <- c(0, 0.25, 1, 2.5, 4.5, 7.0, 10)*10
y <- x[-length(x)]/10
axis(1, at=sqrt (x), labels=x)
axis(2, at=sqrt (y), labels=y)
box()




###################################################
### code chunk number 5: iten-table
###################################################
itenTable <- summaryBy(tf500a1 ~ Sex, data=iten, FUN=list(length, mean, sd))


###################################################
### code chunk number 6: pbcTreatDaten
###################################################
pbcTreatDaten <-
    with (pbcTreat,
          paste ("{\\footnotesize", 
                 formatBig (time), 
                 ifelse (d, "\\hphantom{+}", "+"), 
                 "}",
                 sep = ""))

## in Matrixform bringen
pbcTreatN <- length (pbcTreatDaten)
tableDim <- n2mfrow (pbcTreatN)
pbcTreatDaten <- matrix (c (pbcTreatDaten, rep ("", prod (tableDim) - pbcTreatN)), tableDim[1], byrow = TRUE)
## und als Tabelle ausgeben
w <-
latex(pbcTreatDaten,                     # was konvertieren?
      file = latexTempFile,             
      label = "tab:pbcTreatDaten",       # labeln
      booktabs = TRUE,                  # um booktabs zu benutzen
      colheads = NULL,                  # keine Spaltenüberschriften
      caption = paste("Survival times of", nrow (pbcTreat), "patients under Azathioprine treatment in days. Censored observations are marked with a plus sign."),
      col.just = rep ("R", tableDim[2]),
      center = "none",            # statt \\begin{center}
      where = NULL
      )
postLatex (w, 
           1.0, 
           dropColTitles=TRUE)
## für anschließende Grafik die erste Spalte auswählen:
selection <- seq (from = 1, by = ncol (pbcTreatDaten), length = nrow (pbcTreatDaten))


###################################################
### code chunk number 7: pbcTreatGrafik
###################################################
getOption("SweaveHooks")[["fig"]]()
    ## Funktion zum Plotten von Überlebenszeiten. observed ist Indikator für nicht-zensierte Beob.
    ## obsLabels kann Vektor mit y-Labels für die Zeiten sein.
    timesPlot <- function (times, observed, obsLabels = NULL, ...){
        n <- length (times)
        ## order times
        timesOrder <- rev (order (times))     # rev um lange Dauern unten zu haben
        times <- times[timesOrder]
        observed <- observed[timesOrder]
        obsLabels <- obsLabels[timesOrder]
        ## plot points
        plot (times, 1:n, pch = ifelse (observed, 16, 3), cex=1.75, yaxt = "n", ...)
        legend("topright", legend=c("censored","death"), pch=c(3,19), bty="n") 
       leftlim <- par ("usr")[1]
        ## plot corresponding lines
        segments (leftlim, 1:n, times, 1:n)
        ## optional labeling
        if (!is.null (obsLabels)){
            axis (side = 2, at = 1:n, labels = obsLabels, tick = FALSE)
        }
    }
exampleCases <- pbcTreat[selection,]

par (font.lab = 1)
with(exampleCases,
     timesPlot (time, d, xlab = "Days", ylab="Patient ID",
                obsLabels = selection)
     )


