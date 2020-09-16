###############
# DEPENDENCIES
###############
library(extrafont)
library(fontcm)
loadfonts(quiet = TRUE)
library(plyr)
library(lattice)
library(latticeExtra)
library(wesanderson)

source('./R/functions-Simulations.R')

#######################
# AUXILLIARY FUNCTIONS
#######################

toPdf <- function(expr, filename, ...) {
  toDev(expr, pdf, filename, ...)
}

figPath  <-  function(name) {
  file.path('output/figs', name)
}

toDev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf('Creating %s\n', filename))
  dev(filename, family='CM Roman', ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}



####################
# PLOTTING FUNCTIONS
####################

#' Plot text or points according to relative axis position.
#'
#' @title Plot text or points according to relative axis position
#' @param px Relative x-axis position (in proportion) where character is to be plotted.
#' @param py Relative y-axis position (in proportion) where character is to be plotted.
#' @param lab Plotted text. Works if argument \code{\link[graphics]{text}} is \code{TRUE}.
#' @param adj See argument of same name in R base function \code{\link[graphics]{par}}.
#' @param text Logical. Should text or points be plotted?
#' @param log Used if the original plot uses the argument log, e.g. \code{log='x'}, \code{log='y'} or \code{log='xy'}.
#' @param ... Additional arguments to R base function \code{\link[graphics]{text}}.
#' @export
proportionalLabel <- function(px, py, lab, adj=c(0, 1), text=TRUE, log=FALSE, ...) {
    usr  <-  par('usr')
    x.p  <-  usr[1] + px*(usr[2] - usr[1])
    y.p  <-  usr[3] + py*(usr[4] - usr[3])
    if(log=='x') {
        x.p<-10^(x.p)
    }
    if(log=='y') {
        y.p<-10^(y.p)
    }
    if(log=='xy') {
        x.p<-10^(x.p)
        y.p<-10^(y.p)
    }
    if(text){
        text(x.p, y.p, lab, adj=adj, ...)
    } else {
        points(x.p, y.p, ...)
    }
}

#' Draw equally-spaced white lines on plot window.
#'
#' @title Equally-spaced white lines on plot window
#' @param ... Additional arguments to internal function \code{\link{proportionalLabel}}.
#' @author Diego Barneche
#' @export
plotGrid  <-  function(lineCol='white',...) {
    proportionalLabel(rep(0.2, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.4, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.6, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.8, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.2, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.4, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.6, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.8, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
}


#' Internal. Create nice rounded numbers for plotting.
#'
#' @title Rounded numbers for plotting
#' @param value A numeric vector.
#' @param precision Number of rounding digits.
#' @return A character vector.
#' @author Diego Barneche.
rounded  <-  function(value, precision=1) {
  sprintf(paste0('%.', precision, 'f'), round(value, precision))
}


#' Creates transparent colours
#'
#' @title Creates transparent colours
#' @param col Colour.
#' @param opacity equivalent to alpha transparency parameter
#' @export
transparentColor <- function(col, opacity=0.5) {
    if (length(opacity) > 1 && any(is.na(opacity))) {
        n        <-  max(length(col), length(opacity))
        opacity  <-  rep(opacity, length.out=n)
        col      <-  rep(col, length.out=n)
        ok       <-  !is.na(opacity)
        ret      <-  rep(NA, length(col))
        ret[ok]  <-  Recall(col[ok], opacity[ok])
        ret
    } else {
        tmp  <-  col2rgb(col)/255
        rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
    }
}

fibonacci.scale  <-  function(n) {
    fibs  <-  c(0,1)
    for(i in 2:n) {
        fibs  <-  c(fibs, (fibs[i] + fibs[i-1]))
    }
    (fibs/max(fibs))[-2]
}



##############################################################
##############################################################
##  Preliminary Figs


#' Figs showing Kidwell Funnel plots with simulation results
#' and 1-locus pop. gen. invasion conditions
FunnelPlots  <-  function(df1, df2, df3, df4) {

# Make filenames for import from df names
    fName1  <-  paste('./output/simData/', df1, '.csv', sep="")
    fName2  <-  paste('./output/simData/', df2, '.csv', sep="")
    fName3  <-  paste('./output/simData/', df3, '.csv', sep="")
    fName4  <-  paste('./output/simData/', df4, '.csv', sep="")

# Extract plotting parameter values from df names
    d1   <-  strsplit(df1, '_')[[1]][5:8]
    parsA  <-  list(
                    "hf"     =  as.numeric(strsplit(d1[1],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d1[2],'m')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d1[3],'C')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d1[4],'a')[[1]][2])
                    )
    d2   <-  strsplit(df2, '_')[[1]][5:8]
    parsB  <-  list(
                    "hf"     =  as.numeric(strsplit(d1[1],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d1[2],'m')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d1[3],'C')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d1[4],'a')[[1]][2])
                    )
    d3   <-  strsplit(df3, '_')[[1]][5:8]
    parsC  <-  list(
                    "hf"     =  as.numeric(strsplit(d1[1],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d1[2],'m')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d1[3],'C')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d1[4],'a')[[1]][2])
                    )
    d4   <-  strsplit(df4, '_')[[1]][5:8]
    parsD  <-  list(
                    "hf"     =  as.numeric(strsplit(d1[1],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d1[2],'m')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d1[3],'C')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d1[4],'a')[[1]][2])
                    )

# Import data sets
    panelA  <-  read.csv(fName1, head=TRUE)
    panelB  <-  read.csv(fName2, head=TRUE)
    panelC  <-  read.csv(fName3, head=TRUE)
    panelD  <-  read.csv(fName4, head=TRUE)

# Color scheme
    COLS  <-  list(
                    "invLine"  =  transparentColor('#252525', opacity=1),
                    "poly"     =  transparentColor('green2', opacity=0.6),
                    "aFix"     =  transparentColor('purple', opacity=0.6),
                    "AFix"     =  transparentColor('dodgerblue2', opacity=0.6),
                    "extinct"  =  transparentColor('red', opacity=0.7)
                )
#  Create vector of male selection coefficiets for invasion functions
smLine  <-  seq(0,0.15,length=100)

# Set plot layout
layout.mat <- matrix(c(1,2,3,4), nrow=2, ncol=2, byrow=TRUE)
layout <- layout(layout.mat,respect=TRUE)

##  Row 1: Additive allelic effects

    ##  Panel A: C = 0; hf = hm = 1/2
        dat  <-  panelA
        # Calculate additional variables for plotting
        inv_A  <-  popGen_A_invade(hf=dat$hf[1], hm=dat$hm[1], sm=smLine, C=dat$C[1])
        inv_a  <-  popGen_a_invade(hf=dat$hf[1], hm=dat$hm[1], sm=smLine, C=dat$C[1])
        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        points(sf ~ sm, pch=21, bg=COLS$aFix, data=fix_A)
        points(sf ~ sm, pch=21, bg=COLS$AFix, data=fix_a)
        points(sf[extinct == 1] ~ sm[extinct == 1], pch=21, bg=COLS$extinct, data=dat)
        points(sf[poly == 1 & extinct == 0] ~ sm[poly == 1 & extinct == 0], pch=21, bg=COLS$poly, data=dat)
        #  pop. gen. invasion conditions
        lines(inv_A[inv_A <= 0.15] ~ smLine[inv_A <= 0.15], lwd=2, col=COLS$invLine)
        lines(inv_a[inv_a <= 0.15] ~ smLine[inv_a <= 0.15], lwd=2, col=COLS$invLine)
#        polygon(c(smLine,rev(smLine)), c(inv_a_polygon, rev(inv_A)), col=transparentColor('grey80', 0.6), border='grey70')
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.45, 0.5, substitute(paste(italic(h), " = ", hh), list(hh = parsA$hf)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.3, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, 1.15, substitute(paste(italic(C), ' = ', CC, ", ", delta, ' = ', dd), list(CC=parsA$C, dd=parsA$delta)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)


    ##  Panel B: C = 0.5; hf = hm = 1/2
        dat  <-  panelB
        # Calculate additional variables for plotting
        inv_A  <-  popGen_A_invade(hf=dat$hf[1], hm=dat$hm[1], sm=smLine, C=dat$C[1])
        inv_a  <-  popGen_a_invade(hf=dat$hf[1], hm=dat$hm[1], sm=smLine, C=dat$C[1])
        inv_a_polygon  <-  inv_a
        inv_a_polygon[inv_a > 0.15]  <-  0.15
        x      <-  dat[round(dat$Eq_paa, digits=5)==1,]
        fix_a  <-  x[x$extinct==0,]
        x      <-  dat[round(dat$Eq_pAA, digits=5)==1,]
        fix_A  <-  x[x$extinct==0,]
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2, bty='L')
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        points(sf ~ sm, pch=21, bg=COLS$aFix, data=fix_A)
        points(sf ~ sm, pch=21, bg=COLS$AFix, data=fix_a)
        points(sf[extinct == 1] ~ sm[extinct == 1], pch=21, bg=COLS$extinct, data=dat)
        points(sf[poly == 1 & extinct == 0] ~ sm[poly == 1 & extinct == 0], pch=21, bg=COLS$poly, data=dat)
        #  pop. gen. invasion conditions
        lines(inv_A[inv_A <= 0.15] ~ smLine[inv_A <= 0.15], lwd=2, col=COLS$invLine)
        lines(inv_a[inv_a <= 0.15] ~ smLine[inv_a <= 0.15], lwd=2, col=COLS$invLine)
#        polygon(c(smLine,rev(smLine)), c(inv_a_polygon, rev(inv_A)), col=transparentColor('grey80', 0.6), border='grey70')
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.15, substitute(paste(italic(C), ' = ', CC, ", ", delta, ' = ', dd),list(CC=parsB$C, dd=parsB$delta)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

par(xpd=NA)
legend(#"topright",inset=c(-0.2,0),
            x       =  usr[2]*1.65,
            y       =  usr[4],
        #    title   =  expression(paste(Outcome~of~invasion~analysis)),
            legend  =  c(
                        expression(paste("Extinct")),
                        expression(paste("Female-benefit allele Fixes")),
                        expression(paste("Male-benefit allele Fixes")),
                        expression(paste("Polymorphism"))),
            pch     =  21,
            pt.bg   =  c(COLS$extinct,
                         COLS$aFix,
                         COLS$AFix,
                         COLS$poly),
            cex     =  0.75,
            xjust   =  1,
            yjust   =  1,
            bty     =  'n',
            border  =  NA)

    ##  Panel C: C = 0; hf = hm = 1/4
        dat  <-  panelC
        # Calculate additional variables for plotting
        inv_A  <-  popGen_A_invade(hf=dat$hf[1], hm=dat$hm[1], sm=smLine, C=dat$C[1])
        inv_a  <-  popGen_a_invade(hf=dat$hf[1], hm=dat$hm[1], sm=smLine, C=dat$C[1])
        inv_a_polygon  <-  inv_a
        inv_a_polygon[inv_a > 0.15]  <-  0.15
        x      <-  dat[round(dat$Eq_paa, digits=5)==1,]
        fix_a  <-  x[x$extinct==0,]
        x      <-  dat[round(dat$Eq_pAA, digits=5)==1,]
        fix_A  <-  x[x$extinct==0,]
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        points(sf ~ sm, pch=21, bg=COLS$aFix, data=fix_A)
        points(sf ~ sm, pch=21, bg=COLS$AFix, data=fix_a)
        points(sf[extinct == 1] ~ sm[extinct == 1], pch=21, bg=COLS$extinct, data=dat)
        points(sf[poly == 1 & extinct == 0] ~ sm[poly == 1 & extinct == 0], pch=21, bg=COLS$poly, data=dat)
        #  pop. gen. invasion conditions
        lines(inv_A[inv_A <= 0.15] ~ smLine[inv_A <= 0.15], lwd=2, col=COLS$invLine)
        lines(inv_a[inv_a <= 0.15] ~ smLine[inv_a <= 0.15], lwd=2, col=COLS$invLine)
#        polygon(c(smLine,rev(smLine)), c(inv_a_polygon, rev(inv_A)), col=transparentColor('grey80', 0.6), border='grey70')
        # axes        
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.45, 0.5,substitute(paste(italic(h), " = ", hh), list(hh = parsC$hf)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.3, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)



    ##  Panel D: C = 0.5; hf = hm = 1/2
        dat  <-  panelD
        # Calculate additional variables for plotting
        inv_A  <-  popGen_A_invade(hf=dat$hf[1], hm=dat$hm[1], sm=smLine, C=dat$C[1])
        inv_a  <-  popGen_a_invade(hf=dat$hf[1], hm=dat$hm[1], sm=smLine, C=dat$C[1])
        inv_a_polygon  <-  inv_a
        inv_a_polygon[inv_a > 0.15]  <-  0.15
        x      <-  dat[round(dat$Eq_paa, digits=5)==1,]
        fix_a  <-  x[x$extinct==0,]
        x      <-  dat[round(dat$Eq_pAA, digits=5)==1,]
        fix_A  <-  x[x$extinct==0,]
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        points(sf ~ sm, pch=21, bg=COLS$aFix, data=fix_A)
        points(sf ~ sm, pch=21, bg=COLS$AFix, data=fix_a)
        points(sf[extinct == 1] ~ sm[extinct == 1], pch=21, bg=COLS$extinct, data=dat)
        points(sf[poly == 1 & extinct == 0] ~ sm[poly == 1 & extinct == 0], pch=21, bg=COLS$poly, data=dat)
        #  pop. gen. invasion conditions
        lines(inv_A[inv_A <= 0.15] ~ smLine[inv_A <= 0.15], lwd=2, col=COLS$invLine)
        lines(inv_a[inv_a <= 0.15] ~ smLine[inv_a <= 0.15], lwd=2, col=COLS$invLine)
#        polygon(c(smLine,rev(smLine)), c(inv_a_polygon, rev(inv_A)), col=transparentColor('grey80', 0.6), border='grey70')
        # axes        
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.03, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)


}













#' Figs showing Kidwell Funnel plots with simulation results
#' and 1-locus pop. gen. invasion conditions
FunnelEigSimCompare  <-  function(df1, df2, df3, df4) {

# Make filenames for import from df names
    fName1  <-  paste('./output/simData/', df1, '.csv', sep="")
    fName2  <-  paste('./output/simData/', df2, '.csv', sep="")
    fName3  <-  paste('./output/simData/', df3, '.csv', sep="")
    fName4  <-  paste('./output/simData/', df4, '.csv', sep="")

# Extract plotting parameter values from df names
    d1   <-  strsplit(df1, '_')[[1]][5:8]
    parsA  <-  list(
                    "hf"     =  as.numeric(strsplit(d1[1],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d1[2],'m')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d1[3],'C')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d1[4],'a')[[1]][2])
                    )
    d2   <-  strsplit(df2, '_')[[1]][5:8]
    parsB  <-  list(
                    "hf"     =  as.numeric(strsplit(d2[1],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d2[2],'m')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d2[3],'C')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d2[4],'a')[[1]][2])
                    )
    d3   <-  strsplit(df3, '_')[[1]][5:8]
    parsC  <-  list(
                    "hf"     =  as.numeric(strsplit(d3[1],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d3[2],'m')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d3[3],'C')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d3[4],'a')[[1]][2])
                    )
    d4   <-  strsplit(df4, '_')[[1]][5:8]
    parsD  <-  list(
                    "hf"     =  as.numeric(strsplit(d4[1],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d4[2],'m')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d4[3],'C')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d4[4],'a')[[1]][2])
                    )

# Import data sets
    panelA  <-  read.csv(fName1, head=TRUE)
    panelB  <-  read.csv(fName2, head=TRUE)
    panelC  <-  read.csv(fName3, head=TRUE)
    panelD  <-  read.csv(fName4, head=TRUE)

# Color scheme
    COLS  <-  list(
                    "invLine"  =  transparentColor('#252525', opacity=1),
                    "agree"    =  transparentColor('green', opacity=0.6),
                    "sim"      =  transparentColor('purple', opacity=0.6),
                    "eig"      =  transparentColor('dodgerblue2', opacity=0.6),
                    "extinct"  =  transparentColor('red', opacity=0.7)
                )
#  Create vector of male selection coefficiets for invasion functions
    smLine  <-  seq(0,0.15,length=100)

# Set plot layout
    layout.mat  <- matrix(c(1,2,3,4), nrow=2, ncol=2, byrow=TRUE)
    layout      <- layout(layout.mat,respect=TRUE)

##  Row 1: Additive allelic effects

    ##  Panel A: C = 0; hf = hm = 1/2
        dat  <-  panelA
        # Calculate additional variables for plotting
        inv_A  <-  popGen_A_invade(hf=dat$hf[1], hm=dat$hm[1], sm=smLine, C=dat$C[1])
        inv_a  <-  popGen_a_invade(hf=dat$hf[1], hm=dat$hm[1], sm=smLine, C=dat$C[1])
        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        points(sf[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 1] ~ sm[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 1], pch=21, bg=COLS$agree, data=dat)
        points(sf[dat$zeta_AA < 1 & dat$poly == 1| dat$zeta_aa < 1 & dat$poly == 1] ~ sm[dat$zeta_AA < 1 & dat$poly == 1| dat$zeta_aa < 1 & dat$poly == 1], pch=21, bg=COLS$sim, data=dat)
        points(sf[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 0] ~ sm[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 0], pch=21, bg=COLS$eig, data=dat)
        points(sf[extinct == 1] ~ sm[extinct == 1], pch=21, bg=COLS$extinct, data=dat)
        #  pop. gen. invasion conditions
        lines(inv_A[inv_A <= 0.15] ~ smLine[inv_A <= 0.15], lwd=2, col=COLS$invLine)
        lines(inv_a[inv_a <= 0.15] ~ smLine[inv_a <= 0.15], lwd=2, col=COLS$invLine)
#        polygon(c(smLine,rev(smLine)), c(inv_a_polygon, rev(inv_A)), col=transparentColor('grey80', 0.6), border='grey70')
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.45, 0.5, substitute(paste(italic(h), " = ", hh), list(hh = parsA$hf)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.3, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, 1.15, substitute(paste(italic(C), ' = ', CC, ", ", delta, ' = ', dd), list(CC=parsA$C, dd=parsA$delta)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)


    ##  Panel B: C = 0.5; hf = hm = 1/2
        dat  <-  panelB
        # Calculate additional variables for plotting
        inv_A  <-  popGen_A_invade(hf=dat$hf[1], hm=dat$hm[1], sm=smLine, C=dat$C[1])
        inv_a  <-  popGen_a_invade(hf=dat$hf[1], hm=dat$hm[1], sm=smLine, C=dat$C[1])
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2, bty='L')
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        points(sf[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 1] ~ sm[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 1], pch=21, bg=COLS$agree, data=dat)
        points(sf[dat$zeta_AA < 1 & dat$poly == 1| dat$zeta_aa < 1 & dat$poly == 1] ~ sm[dat$zeta_AA < 1 & dat$poly == 1| dat$zeta_aa < 1 & dat$poly == 1], pch=21, bg=COLS$sim, data=dat)
        points(sf[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 0] ~ sm[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 0], pch=21, bg=COLS$eig, data=dat)
        points(sf[extinct == 1] ~ sm[extinct == 1], pch=21, bg=COLS$extinct, data=dat)
        #  pop. gen. invasion conditions
        lines(inv_A[inv_A <= 0.15] ~ smLine[inv_A <= 0.15], lwd=2, col=COLS$invLine)
        lines(inv_a[inv_a <= 0.15] ~ smLine[inv_a <= 0.15], lwd=2, col=COLS$invLine)
#        polygon(c(smLine,rev(smLine)), c(inv_a_polygon, rev(inv_A)), col=transparentColor('grey80', 0.6), border='grey70')
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.15, substitute(paste(italic(C), ' = ', CC, ", ", delta, ' = ', dd),list(CC=parsB$C, dd=parsB$delta)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

par(xpd=NA)
legend(#"topright",inset=c(-0.2,0),
            x       =  usr[2]*1.65,
            y       =  usr[4],
        #    title   =  expression(paste(Outcome~of~invasion~analysis)),
            legend  =  c(
                        expression(paste("Agree")),
                        expression(paste("Sim. Polymorphic Only")),
                        expression(paste(zeta[i]," Polymorphic Only")),
                        expression(paste("Extinct"))),
            pch     =  21,
            pt.bg   =  c(COLS$agree,
                         COLS$sim,
                         COLS$eig,
                         COLS$extinct),
            cex     =  0.75,
            xjust   =  1,
            yjust   =  1,
            bty     =  'n',
            border  =  NA)

    ##  Panel C: C = 0; hf = hm = 1/4
        dat  <-  panelC
        # Calculate additional variables for plotting
        inv_A  <-  popGen_A_invade(hf=dat$hf[1], hm=dat$hm[1], sm=smLine, C=dat$C[1])
        inv_a  <-  popGen_a_invade(hf=dat$hf[1], hm=dat$hm[1], sm=smLine, C=dat$C[1])
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        points(sf[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 1] ~ sm[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 1], pch=21, bg=COLS$agree, data=dat)
        points(sf[dat$zeta_AA < 1 & dat$poly == 1| dat$zeta_aa < 1 & dat$poly == 1] ~ sm[dat$zeta_AA < 1 & dat$poly == 1| dat$zeta_aa < 1 & dat$poly == 1], pch=21, bg=COLS$sim, data=dat)
        points(sf[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 0] ~ sm[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 0], pch=21, bg=COLS$eig, data=dat)
        points(sf[extinct == 1] ~ sm[extinct == 1], pch=21, bg=COLS$extinct, data=dat)
        #  pop. gen. invasion conditions
        lines(inv_A[inv_A <= 0.15] ~ smLine[inv_A <= 0.15], lwd=2, col=COLS$invLine)
        lines(inv_a[inv_a <= 0.15] ~ smLine[inv_a <= 0.15], lwd=2, col=COLS$invLine)
#        polygon(c(smLine,rev(smLine)), c(inv_a_polygon, rev(inv_A)), col=transparentColor('grey80', 0.6), border='grey70')
        # axes        
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.45, 0.5,substitute(paste(italic(h), " = ", hh), list(hh = parsC$hf)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.3, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)



    ##  Panel D: C = 0.5; hf = hm = 1/2
        dat  <-  panelD
        # Calculate additional variables for plotting
        inv_A  <-  popGen_A_invade(hf=dat$hf[1], hm=dat$hm[1], sm=smLine, C=dat$C[1])
        inv_a  <-  popGen_a_invade(hf=dat$hf[1], hm=dat$hm[1], sm=smLine, C=dat$C[1])
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        points(sf[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 1] ~ sm[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 1], pch=21, bg=COLS$agree, data=dat)
        points(sf[dat$zeta_AA < 1 & dat$poly == 1| dat$zeta_aa < 1 & dat$poly == 1] ~ sm[dat$zeta_AA < 1 & dat$poly == 1| dat$zeta_aa < 1 & dat$poly == 1], pch=21, bg=COLS$sim, data=dat)
        points(sf[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 0] ~ sm[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 0], pch=21, bg=COLS$eig, data=dat)
        points(sf[extinct == 1] ~ sm[extinct == 1], pch=21, bg=COLS$extinct, data=dat)
        #  pop. gen. invasion conditions
        lines(inv_A[inv_A <= 0.15] ~ smLine[inv_A <= 0.15], lwd=2, col=COLS$invLine)
        lines(inv_a[inv_a <= 0.15] ~ smLine[inv_a <= 0.15], lwd=2, col=COLS$invLine)
#        polygon(c(smLine,rev(smLine)), c(inv_a_polygon, rev(inv_A)), col=transparentColor('grey80', 0.6), border='grey70')
        # axes        
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.03, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)


}