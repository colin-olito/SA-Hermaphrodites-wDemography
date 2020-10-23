###############
# DEPENDENCIES
###############
library(extrafont)
loadfonts(quiet = TRUE)
library(fontcm)
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
    d1   <-  strsplit(df1, '_')[[1]][c(2,4:7)]
    parsA  <-  list(
                    "sMax"   =  as.numeric(strsplit(d1[1],'x')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d1[1],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d1[2],'m')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d1[3],'C')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d1[4],'a')[[1]][2])
                    )
    d2   <-  strsplit(df2, '_')[[1]][c(2,4:7)]
    parsB  <-  list(
                    "sMax"   =  as.numeric(strsplit(d2[1],'x')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d2[1],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d2[2],'m')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d2[3],'C')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d2[4],'a')[[1]][2])
                    )
    d3   <-  strsplit(df3, '_')[[1]][c(2,4:7)]
    parsC  <-  list(
                    "sMax"   =  as.numeric(strsplit(d3[1],'x')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d3[1],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d3[2],'m')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d3[3],'C')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d3[4],'a')[[1]][2])
                    )
    d4   <-  strsplit(df4, '_')[[1]][c(2,4:7)]
    parsD  <-  list(
                    "sMax"   =  as.numeric(strsplit(d4[1],'x')[[1]][2]),
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
                    "poly"     =  transparentColor('green2', opacity=0.6),
                    "aFix"     =  transparentColor('purple', opacity=0.6),
                    "AFix"     =  transparentColor('dodgerblue2', opacity=0.6),
                    "extinct"  =  transparentColor('red', opacity=0.7)
                )
#  Create vector of male selection coefficients for invasion functions
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
        fix_a  <-  dat[dat$zeta_AA < 1 & dat$zeta_aa > 1 & dat$extinct == 0,]
        fix_A  <-  dat[dat$zeta_AA > 1 & dat$zeta_aa < 1 & dat$extinct == 0,]
        poly   <-  dat[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$extinct == 0,]
        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,parsA$sMax), ylim = c(0,parsA$sMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        points(sf ~ sm, pch=21, bg=COLS$aFix, data=fix_A)
        points(sf ~ sm, pch=21, bg=COLS$AFix, data=fix_a)
        points(sf[extinct == 1] ~ sm[extinct == 1], pch=21, bg=COLS$extinct, data=dat)
        points(sf[extinct == 0] ~ sm[extinct == 0], pch=21, bg=COLS$poly, data=poly)
        #  pop. gen. invasion conditions
        lines(inv_A[inv_A <= parsA$sMax] ~ smLine[inv_A <= parsA$sMax], lwd=2, col=COLS$invLine)
        lines(inv_a[inv_a <= parsA$sMax] ~ smLine[inv_a <= parsA$sMax], lwd=2, col=COLS$invLine)
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
        fix_a  <-  dat[dat$zeta_AA < 1 & dat$zeta_aa > 1 & dat$extinct == 0,]
        fix_A  <-  dat[dat$zeta_AA > 1 & dat$zeta_aa < 1 & dat$extinct == 0,]
        poly   <-  dat[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$extinct == 0,]
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,parsB$sMax), ylim = c(0,parsB$sMax), ylab='', xlab='', cex.lab=1.2, bty='L')
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        points(sf ~ sm, pch=21, bg=COLS$aFix, data=fix_A)
        points(sf ~ sm, pch=21, bg=COLS$AFix, data=fix_a)
        points(sf[extinct == 1] ~ sm[extinct == 1], pch=21, bg=COLS$extinct, data=dat)
        points(sf[extinct == 0] ~ sm[extinct == 0], pch=21, bg=COLS$poly, data=poly)
        #  pop. gen. invasion conditions
        lines(inv_A[inv_A <= parsB$sMax] ~ smLine[inv_A <= parsB$sMax], lwd=2, col=COLS$invLine)
        lines(inv_a[inv_a <= parsB$sMax] ~ smLine[inv_a <= parsB$sMax], lwd=2, col=COLS$invLine)
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
        fix_a  <-  dat[dat$zeta_AA < 1 & dat$zeta_aa > 1 & dat$extinct == 0,]
        fix_A  <-  dat[dat$zeta_AA > 1 & dat$zeta_aa < 1 & dat$extinct == 0,]
        poly   <-  dat[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$extinct == 0,]
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,parsC$sMax), ylim = c(0,parsC$sMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        points(sf ~ sm, pch=21, bg=COLS$aFix, data=fix_A)
        points(sf ~ sm, pch=21, bg=COLS$AFix, data=fix_a)
        points(sf[extinct == 1] ~ sm[extinct == 1], pch=21, bg=COLS$extinct, data=dat)
        points(sf[extinct == 0] ~ sm[extinct == 0], pch=21, bg=COLS$poly, data=poly)
        #  pop. gen. invasion conditions
        lines(inv_A[inv_A <= parsC$sMax] ~ smLine[inv_A <= parsC$sMax], lwd=2, col=COLS$invLine)
        lines(inv_a[inv_a <= parsC$sMax] ~ smLine[inv_a <= parsC$sMax], lwd=2, col=COLS$invLine)
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
        fix_a  <-  dat[dat$zeta_AA < 1 & dat$zeta_aa > 1 & dat$extinct == 0,]
        fix_A  <-  dat[dat$zeta_AA > 1 & dat$zeta_aa < 1 & dat$extinct == 0,]
        poly   <-  dat[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$extinct == 0,]
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,parsD$sMax), ylim = c(0,parsD$sMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        points(sf ~ sm, pch=21, bg=COLS$aFix, data=fix_A)
        points(sf ~ sm, pch=21, bg=COLS$AFix, data=fix_a)
        points(sf[extinct == 1] ~ sm[extinct == 1], pch=21, bg=COLS$extinct, data=dat)
        points(sf[extinct == 0] ~ sm[extinct == 0], pch=21, bg=COLS$poly, data=poly)
        #  pop. gen. invasion conditions
        lines(inv_A[inv_A <= parsD$sMax] ~ smLine[inv_A <= parsD$sMax], lwd=2, col=COLS$invLine)
        lines(inv_a[inv_a <= parsD$sMax] ~ smLine[inv_a <= parsD$sMax], lwd=2, col=COLS$invLine)
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
    d1   <-  strsplit(df1, '_')[[1]][c(2,4:7)]
    parsA  <-  list(
                    "sMax"   =  as.numeric(strsplit(d1[1],'x')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d1[1],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d1[2],'m')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d1[3],'C')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d1[4],'a')[[1]][2])
                    )
    d2   <-  strsplit(df2, '_')[[1]][c(2,4:7)]
    parsB  <-  list(
                    "sMax"   =  as.numeric(strsplit(d2[1],'x')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d2[1],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d2[2],'m')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d2[3],'C')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d2[4],'a')[[1]][2])
                    )
    d3   <-  strsplit(df3, '_')[[1]][c(2,4:7)]
    parsC  <-  list(
                    "sMax"   =  as.numeric(strsplit(d3[1],'x')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d3[1],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d3[2],'m')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d3[3],'C')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d3[4],'a')[[1]][2])
                    )
    d4   <-  strsplit(df4, '_')[[1]][c(2,4:7)]
    parsD  <-  list(
                    "sMax"   =  as.numeric(strsplit(d4[1],'x')[[1]][2]),
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
#  Create vector of male selection coefficients for invasion functions
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
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,parsA$sMax), ylim = c(0,parsA$sMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        points(sf[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 1] ~ sm[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 1], pch=21, bg=COLS$agree, data=dat)
        points(sf[dat$zeta_AA < 1 & dat$poly == 1 | dat$zeta_aa < 1 & dat$poly == 1] ~ sm[dat$zeta_AA < 1 & dat$poly == 1| dat$zeta_aa < 1 & dat$poly == 1], pch=21, bg=COLS$sim, data=dat)
        points(sf[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 0] ~ sm[dat$zeta_AA > 1 & dat$zeta_aa > 1 & dat$poly == 0], pch=21, bg=COLS$eig, data=dat)
        points(sf[extinct == 1] ~ sm[extinct == 1], pch=21, bg=COLS$extinct, data=dat)
        #  pop. gen. invasion conditions
        lines(inv_A[inv_A <= parsA$sMax] ~ smLine[inv_A <= parsA$sMax], lwd=2, col=COLS$invLine)
        lines(inv_a[inv_a <= parsA$sMax] ~ smLine[inv_a <= parsA$sMax], lwd=2, col=COLS$invLine)
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
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,parsB$sMax), ylim = c(0,parsB$sMax), ylab='', xlab='', cex.lab=1.2, bty='L')
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
        lines(inv_A[inv_A <= parsB$sMax] ~ smLine[inv_A <= parsB$sMax], lwd=2, col=COLS$invLine)
        lines(inv_a[inv_a <= parsB$sMax] ~ smLine[inv_a <= parsB$sMax], lwd=2, col=COLS$invLine)
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
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,parsC$sMax), ylim = c(0,parsC$sMax), ylab='', xlab='', cex.lab=1.2)
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
        lines(inv_A[inv_A <= parsC$sMax] ~ smLine[inv_A <= parsC$sMax], lwd=2, col=COLS$invLine)
        lines(inv_a[inv_a <= parsC$sMax] ~ smLine[inv_a <= parsC$sMax], lwd=2, col=COLS$invLine)
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
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,parsD$sMax), ylim = c(0,parsD$sMax), ylab='', xlab='', cex.lab=1.2)
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
        lines(inv_A[inv_A <= parsD$sMax] ~ smLine[inv_A <= parsD$sMax], lwd=2, col=COLS$invLine)
        lines(inv_a[inv_a <= parsD$sMax] ~ smLine[inv_a <= parsD$sMax], lwd=2, col=COLS$invLine)
#        polygon(c(smLine,rev(smLine)), c(inv_a_polygon, rev(inv_A)), col=transparentColor('grey80', 0.6), border='grey70')
        # axes        
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.03, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)


}



#' Fig showing Proportion of polymorphic parameter space
#' 
polySpaceFig  <-  function(df1, df2, df3, df4, df5, df6) {

# Make filenames for import from df names
    fName1  <-  paste('./output/simData/', df1, '.csv', sep="")
    fName2  <-  paste('./output/simData/', df2, '.csv', sep="")
    fName3  <-  paste('./output/simData/', df3, '.csv', sep="")
    fName4  <-  paste('./output/simData/', df4, '.csv', sep="")
    fName5  <-  paste('./output/simData/', df5, '.csv', sep="")
    fName6  <-  paste('./output/simData/', df6, '.csv', sep="")

# Extract plotting parameter values from df names
    d1   <-  strsplit(df1, '_')[[1]][c(2,4:10)]
    pars1  <-  list(
                    "sMax"   =  as.numeric(strsplit(d1[1],'x')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d1[2],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d1[3],'m')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d1[4],'a')[[1]][2]),
                    "dj"     =  as.numeric(strsplit(d1[5],'j')[[1]][2]),
                    "da"     =  as.numeric(strsplit(d1[6],'a')[[1]][2]),
                    "dg"     =  as.numeric(strsplit(d1[7],'g')[[1]][2]),
                    "f"      =  as.numeric(strsplit(d1[8],'f')[[1]][2])
                    )
    d2   <-  strsplit(df2, '_')[[1]][c(2,4:10)]
    pars2  <-  list(
                    "sMax"   =  as.numeric(strsplit(d2[1],'x')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d2[2],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d2[3],'m')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d2[4],'a')[[1]][2]),
                    "dj"     =  as.numeric(strsplit(d2[5],'j')[[1]][2]),
                    "da"     =  as.numeric(strsplit(d2[6],'a')[[1]][2]),
                    "dg"     =  as.numeric(strsplit(d2[7],'g')[[1]][2]),
                    "f"      =  as.numeric(strsplit(d2[8],'f')[[1]][2])
                    )
    d3   <-  strsplit(df3, '_')[[1]][c(2,4:10)]
    pars3  <-  list(
                    "sMax"   =  as.numeric(strsplit(d3[1],'x')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d3[2],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d3[3],'m')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d3[4],'a')[[1]][2]),
                    "dj"     =  as.numeric(strsplit(d3[5],'j')[[1]][2]),
                    "da"     =  as.numeric(strsplit(d3[6],'a')[[1]][2]),
                    "dg"     =  as.numeric(strsplit(d3[7],'g')[[1]][2]),
                    "f"      =  as.numeric(strsplit(d3[8],'f')[[1]][2])
                    )
    d4   <-  strsplit(df4, '_')[[1]][c(2,4:10)]
    pars4  <-  list(
                    "sMax"   =  as.numeric(strsplit(d4[1],'x')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d4[2],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d4[3],'m')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d4[4],'a')[[1]][2]),
                    "dj"     =  as.numeric(strsplit(d4[5],'j')[[1]][2]),
                    "da"     =  as.numeric(strsplit(d4[6],'a')[[1]][2]),
                    "dg"     =  as.numeric(strsplit(d4[7],'g')[[1]][2]),
                    "f"      =  as.numeric(strsplit(d4[8],'f')[[1]][2])
                    )
    d5   <-  strsplit(df5, '_')[[1]][c(2,4:10)]
    pars5  <-  list(
                    "sMax"   =  as.numeric(strsplit(d5[1],'x')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d5[2],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d5[3],'m')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d5[4],'a')[[1]][2]),
                    "dj"     =  as.numeric(strsplit(d5[5],'j')[[1]][2]),
                    "da"     =  as.numeric(strsplit(d5[6],'a')[[1]][2]),
                    "dg"     =  as.numeric(strsplit(d5[7],'g')[[1]][2]),
                    "f"      =  as.numeric(strsplit(d5[8],'f')[[1]][2])
                    )
    d6   <-  strsplit(df6, '_')[[1]][c(2,4:10)]
    pars6  <-  list(
                    "sMax"   =  as.numeric(strsplit(d6[1],'x')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d6[2],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d6[3],'m')[[1]][2]),
                    "delta"  =  as.numeric(strsplit(d6[4],'a')[[1]][2]),
                    "dj"     =  as.numeric(strsplit(d6[5],'j')[[1]][2]),
                    "da"     =  as.numeric(strsplit(d6[6],'a')[[1]][2]),
                    "dg"     =  as.numeric(strsplit(d6[7],'g')[[1]][2]),
                    "f"      =  as.numeric(strsplit(d6[8],'f')[[1]][2])
                    )

# Import data sets
    plt1  <-  read.csv(fName1, head=TRUE)
    plt2  <-  read.csv(fName2, head=TRUE)
    plt3  <-  read.csv(fName3, head=TRUE)
    plt4  <-  read.csv(fName4, head=TRUE)
    plt5  <-  read.csv(fName5, head=TRUE)
    plt6  <-  read.csv(fName6, head=TRUE)

# Color scheme
    COLS  <-  list(
                    "PG"   =  transparentColor('#252525', opacity=1),
                    "low"  =  transparentColor('dodgerblue4', opacity=0.6),
                    "med"  =  transparentColor('darkolivegreen4', opacity=0.6),
                    "hi"   =  transparentColor('tomato', opacity=0.6),
                    "low2"  =  transparentColor('dodgerblue4', opacity=1),
                    "med2"  =  transparentColor('darkolivegreen4', opacity=1),
                    "hi2"   =  transparentColor('tomato', opacity=1)                )
#  Create vector of selfing rates for pop gen function.
    CLine          <-  seq(0,0.9,length=100)
    addPGSpace     <-  c()
    domRevPGSpace  <-  c()
    for(i in 1:length(CLine)) {
        addPGSpace[i]     <-  popGen_PolySpace(hf=pars1$hf, hm=pars1$hm, C=CLine[i], sMax=pars1$sMax)
        domRevPGSpace[i]  <-  popGen_PolySpace(hf=pars4$hf, hm=pars4$hm, C=CLine[i], sMax=pars4$sMax)
    }

# Set plot layout
    layout.mat  <- matrix(c(1,2), nrow=2, ncol=1, byrow=TRUE)
    layout      <- layout(layout.mat,respect=TRUE)

    ##  Panel A: hf = hm = 1/2
        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(plt2$Cs)), ylim = c(0,0.105), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(addPGSpace ~ CLine, lwd=2, col=COLS$PG)
        points(eigPolyViable ~ Cs, pch=21, bg=COLS$low, col=COLS$low2, data=plt1)
        points(eigPolyViable ~ Cs, pch=21, bg=COLS$med, col=COLS$med2, data=plt2)
        points(eigPolyViable ~ Cs, pch=21, bg=COLS$hi, col=COLS$hi2, data=plt3)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.1, substitute(paste(italic(h), " = ", hh), list(hh = pars2$hf)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.5, 1.1, expression(paste(italic(h), " = ", 1/2)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(-0.3, 0.5, expression(paste("Proportion of parameter space")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=90)

        #Legend
        legend( x       =  usr[2],
                y       =  usr[4],
                legend  =  c(
                             expression(paste("Pop. Gen.")),
                             substitute(paste("High fertility (", italic(f), " = ", ff, ")"), list(ff = pars3$f)),
                             substitute(paste("Med. fertility (", italic(f), " = ", ff, ".0)"), list(ff = pars2$f)),
                             substitute(paste("Low  fertility (", italic(f), " = ", ff, ")"), list(ff = pars1$f))),
                 lty     =  c(1,NA,NA,NA),
                 lwd     =  c(2,NA,NA,NA),
                 col     =  c(COLS$PG,
                              COLS$hi2,
                              COLS$med2,
                              COLS$low2),
                 pch     =  c(NA,21,21,21),
                 pt.bg   =  c(NA,
                              COLS$hi,
                              COLS$med,
                              COLS$low),
                 cex     =  0.65,
                 pt.cex  =  0.75,
                 xjust   =  1,
                 yjust   =  1,
                 bty     =  'n',
                 border  =  NA)

    ##  Panel B: hf = hm = 1/2
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(plt5$Cs)), ylim = c(0,(max(plt5$popGenPoly)*1.05)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(domRevPGSpace ~ CLine, lwd=2, col=COLS$PG)
        points(eigPolyViable ~ Cs, pch=21, bg=COLS$low, col=COLS$low2, data=plt4)
        points(eigPolyViable ~ Cs, pch=21, bg=COLS$med, col=COLS$med2, data=plt5)
        points(eigPolyViable ~ Cs, pch=21, bg=COLS$hi, col=COLS$hi2, data=plt6)
        # axes        
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(0.03, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.1, substitute(paste(italic(h), " = ", hh), list(hh = pars5$hf)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.5, 1.1, expression(paste(italic(h), " = ", 1/4)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.5, -0.3, expression(paste("Selfing rate (",italic(C), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(-0.3, 0.5, expression(paste("Proportion of parameter space")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=90)

        #Legend
        legend( x       =  usr[2],
                y       =  usr[4],
                legend  =  c(
                             expression(paste("Pop. Gen.")),
                             substitute(paste("High fertility (", italic(f), " = ", ff, ")"), list(ff = pars6$f)),
                             substitute(paste("Med. fertility (", italic(f), " = ", ff, ")"), list(ff = pars5$f)),
                             substitute(paste("Low  fertility (", italic(f), " = ", ff, ")"), list(ff = pars4$f))),
                 lty     =  c(1,NA,NA,NA),
                 lwd     =  c(2,NA,NA,NA),
                 col     =  c(COLS$PG,
                              COLS$hi2,
                              COLS$med2,
                              COLS$low2),
                 pch     =  c(NA,21,21,21),
                 pt.bg   =  c(NA,
                              COLS$hi,
                              COLS$med,
                              COLS$low),
                 cex     =  0.65,
                 pt.cex  =  0.75,
                 xjust   =  1,
                 yjust   =  1,
                 bty     =  'n',
                 border  =  NA)        
}













deltaEffectsPolySpaceFig  <-  function(df1, df2, df3, df4, df5, df6) {

# Make filenames for import from df names
    fName1  <-  paste('./output/simData/', df1, '.csv', sep="")
    fName2  <-  paste('./output/simData/', df2, '.csv', sep="")
    fName3  <-  paste('./output/simData/', df3, '.csv', sep="")
    fName4  <-  paste('./output/simData/', df4, '.csv', sep="")
    fName5  <-  paste('./output/simData/', df5, '.csv', sep="")
    fName6  <-  paste('./output/simData/', df6, '.csv', sep="")

# Extract plotting parameter values from df names
    d1   <-  strsplit(df1, '_')[[1]][c(2,4:7)]
    pars1  <-  list(
                    "sMax"   =  as.numeric(strsplit(d1[1],'x')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d1[2],'C')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d1[3],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d1[4],'m')[[1]][2]),
                    "f"      =  as.numeric(strsplit(d1[5],'f')[[1]][2])
                    )
    d2   <-  strsplit(df2, '_')[[1]][c(2,4:7)]
    pars2  <-  list(
                    "sMax"   =  as.numeric(strsplit(d2[1],'x')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d2[2],'C')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d2[3],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d2[4],'m')[[1]][2]),
                    "f"      =  as.numeric(strsplit(d2[5],'f')[[1]][2])
                    )
    d3   <-  strsplit(df3, '_')[[1]][c(2,4:7)]
    pars3  <-  list(
                    "sMax"   =  as.numeric(strsplit(d3[1],'x')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d3[2],'C')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d3[3],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d3[4],'m')[[1]][2]),
                    "f"      =  as.numeric(strsplit(d3[5],'f')[[1]][2])
                    )
    d4   <-  strsplit(df4, '_')[[1]][c(2,4:7)]
    pars4  <-  list(
                    "sMax"   =  as.numeric(strsplit(d4[1],'x')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d4[2],'C')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d4[3],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d4[4],'m')[[1]][2]),
                    "f"      =  as.numeric(strsplit(d4[5],'f')[[1]][2])
                    )
    d5   <-  strsplit(df5, '_')[[1]][c(2,4:7)]
    pars5  <-  list(
                    "sMax"   =  as.numeric(strsplit(d5[1],'x')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d5[2],'C')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d5[3],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d5[4],'m')[[1]][2]),
                    "f"      =  as.numeric(strsplit(d5[5],'f')[[1]][2])
                    )
    d6   <-  strsplit(df6, '_')[[1]][c(2,4:7)]
    pars6  <-  list(
                    "sMax"   =  as.numeric(strsplit(d6[1],'x')[[1]][2]),
                    "C"      =  as.numeric(strsplit(d6[2],'C')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d6[3],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d6[4],'m')[[1]][2]),
                    "f"      =  as.numeric(strsplit(d6[5],'f')[[1]][2])
                    )

# Import data sets
    plt1  <-  read.csv(fName1, head=TRUE)
    plt2  <-  read.csv(fName2, head=TRUE)
    plt3  <-  read.csv(fName3, head=TRUE)
    plt4  <-  read.csv(fName4, head=TRUE)
    plt5  <-  read.csv(fName5, head=TRUE)
    plt6  <-  read.csv(fName6, head=TRUE)

# Color scheme
    COLS  <-  list(
                    "PG"     =  transparentColor('#252525', opacity=1),
                    "dSim"   =  transparentColor('#252525', opacity=0.6),
                    "d_j"    =  transparentColor('dodgerblue4', opacity=0.6),
                    "d_a"    =  transparentColor('dodgerblue', opacity=0.6),
                    "d_g"    =  transparentColor('tomato', opacity=0.6),
                    "dSim2"  =  transparentColor('#252525', opacity=1),
                    "d_j2"   =  transparentColor('dodgerblue4', opacity=1),
                    "d_a2"   =  transparentColor('dodgerblue', opacity=1),
                    "d_g2"   =  transparentColor('tomato', opacity=1)
                    )
#  Create vector of delta values for pop gen predictions.
    dLine          <-  seq(0,0.9,length=100)

# Set plot layout
    layout.mat  <- matrix(c(1:6), nrow=2, ncol=3, byrow=TRUE)
    layout      <- layout(layout.mat,respect=TRUE)

## Row 1: Additive SA (hf = hm = 1/2)
    ##  Panel A: C = 1/4
        PGSpace     <-  c()
        for(i in 1:length(dLine)) {
            PGSpace[i]     <-  popGen_PolySpace_Delta_Add(C=pars1$C, delta=dLine[i], sMax=pars1$sMax)
        }
        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(plt1$d)), ylim = c(0,0.0925), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(PGSpace ~ dLine, lwd=2, col=COLS$PG)
        points(d_eigPolyViable ~ d, pch=21, bg=COLS$dSim, col=COLS$dSim2, data=plt1)
        points(d_j_eigPolyViable ~ d_j, pch=21, bg=COLS$d_j, col=COLS$d_j2, data=plt1)
        points(d_a_eigPolyViable ~ d_a, pch=21, bg=COLS$d_a, col=COLS$d_a2, data=plt1)
        points(d_g_eigPolyViable ~ d_g, pch=21, bg=COLS$d_g, col=COLS$d_g2, data=plt1)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, substitute(paste(italic(C), " = ", cc), list(cc = pars1$C)), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(-0.5, 0.5, expression(paste(italic(h), " = ", 1/2)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.35, -0.15, expression(paste("Proportion of parameter space")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)



    ##  Panel B: C = 1/2
        PGSpace     <-  c()
        for(i in 1:length(dLine)) {
            PGSpace[i]     <-  popGen_PolySpace_Delta_Add(C=pars2$C, delta=dLine[i], sMax=pars2$sMax)
        }
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(plt2$d)), ylim = c(0,0.0925), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(PGSpace ~ dLine, lwd=2, col=COLS$PG)
        points(d_eigPolyViable ~ d, pch=21, bg=COLS$dSim, col=COLS$dSim2, data=plt2)
        points(d_j_eigPolyViable ~ d_j, pch=21, bg=COLS$d_j, col=COLS$d_j2, data=plt2)
        points(d_a_eigPolyViable ~ d_a, pch=21, bg=COLS$d_a, col=COLS$d_a2, data=plt2)
        points(d_g_eigPolyViable ~ d_g, pch=21, bg=COLS$d_g, col=COLS$d_g2, data=plt2)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.03, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, substitute(paste(italic(C), " = ", cc), list(cc = pars2$C)), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)


    ##  Panel C: C = 3/4
        PGSpace     <-  c()
        for(i in 1:length(dLine)) {
            PGSpace[i]     <-  popGen_PolySpace_Delta_Add(C=0.75, delta=dLine[i], sMax=pars3$sMax)
        }
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(plt3$d)), ylim = c(0,0.0925), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(PGSpace ~ dLine, lwd=2, col=COLS$PG)
        points(d_eigPolyViable ~ d, pch=21, bg=COLS$dSim, col=COLS$dSim2, data=plt3)
        points(d_j_eigPolyViable ~ d_j, pch=21, bg=COLS$d_j, col=COLS$d_j2, data=plt3)
        points(d_a_eigPolyViable ~ d_a, pch=21, bg=COLS$d_a, col=COLS$d_a2, data=plt3)
        points(d_g_eigPolyViable ~ d_g, pch=21, bg=COLS$d_g, col=COLS$d_g2, data=plt3)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.03, 1.075, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, substitute(paste(italic(C), " = ", cc), list(cc = pars3$C)), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)

        #Legend
        legend( x       =  usr[2]*0.45,
                y       =  usr[4],
                legend  =  c(
                             expression(paste("Pop. Gen. ", italic(delta))),
                             expression(paste(delta)),
                             expression(paste(delta[italic(j)])),
                             expression(paste(delta[italic(a)])),
                             expression(paste(delta[gamma]))),
                 lty     =  c(1,NA,NA,NA,NA),
                 lwd     =  c(2,NA,NA,NA,NA),
                 col     =  c(COLS$PG,
                              COLS$dSim,
                              COLS$d_j,
                              COLS$d_a,
                              COLS$d_g),
                 pch     =  c(NA,21,21,21,21),
                 pt.bg   =  c(NA,
                              COLS$dSim,
                              COLS$d_j,
                              COLS$d_a,
                              COLS$d_g),
                 cex     =  0.75,
                 pt.cex  =  0.75,
                 xjust   =  1,
                 yjust   =  1,
                 bty     =  'n',
                 border  =  NA)
    ##  DOMINANCE REVERSAL: hf = hm = 1/4
    ##  Panel D: C = 1/4
        PGSpace     <-  c()
        for(i in 1:length(dLine)) {
            PGSpace[i]     <-  popGen_PolySpace_Delta_DomRev(C=pars4$C, delta=dLine[i], sMax=pars4$sMax)
        }
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(plt4$d)), ylim = c(0,max(PGSpace)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(PGSpace ~ dLine, lwd=2, col=COLS$PG)
        points(d_eigPolyViable ~ d, pch=21, bg=COLS$dSim, col=COLS$dSim2, data=plt4)
        points(d_j_eigPolyViable ~ d_j, pch=21, bg=COLS$d_j, col=COLS$d_j2, data=plt4)
        points(d_a_eigPolyViable ~ d_a, pch=21, bg=COLS$d_a, col=COLS$d_a2, data=plt4)
        points(d_g_eigPolyViable ~ d_g, pch=21, bg=COLS$d_g, col=COLS$d_g2, data=plt4)
        # axes        
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(0.03, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.1, substitute(paste(italic(C), " = ", cc), list(cc = pars4$C)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(-0.5, 0.5, expression(paste(italic(h), " = ", 1/4)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)

    ##  Panel E: C = 1/2
        PGSpace     <-  c()
        for(i in 1:length(dLine)) {
            PGSpace[i]     <-  popGen_PolySpace_Delta_DomRev(C=pars5$C, delta=dLine[i], sMax=pars5$sMax)
        }
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(plt5$d)), ylim = c(0,0.65), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(PGSpace ~ dLine, lwd=2, col=COLS$PG)
        points(d_eigPolyViable ~ d, pch=21, bg=COLS$dSim, col=COLS$dSim2, data=plt5)
        points(d_j_eigPolyViable ~ d_j, pch=21, bg=COLS$d_j, col=COLS$d_j2, data=plt5)
        points(d_a_eigPolyViable ~ d_a, pch=21, bg=COLS$d_a, col=COLS$d_a2, data=plt5)
        points(d_g_eigPolyViable ~ d_g, pch=21, bg=COLS$d_g, col=COLS$d_g2, data=plt5)
        # axes        
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.03, 1.075, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.1, substitute(paste(italic(C), " = ", cc), list(cc = pars5$C)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.5, -0.35, expression(paste("Inbreeding depression (", delta[italic(i)], ")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)

    ##  Panel F: C = 3/4
        PGSpace     <-  c()
        for(i in 1:length(dLine)) {
            PGSpace[i]     <-  popGen_PolySpace_Delta_DomRev(C=0.75, delta=dLine[i], sMax=pars6$sMax)
        }
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(plt6$d)), ylim = c(0,max(plt6$d_popGenPoly)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(PGSpace ~ dLine, lwd=2, col=COLS$PG)
        points(d_eigPolyViable ~ d, pch=21, bg=COLS$dSim, col=COLS$dSim2, data=plt6)
        points(d_j_eigPolyViable ~ d_j, pch=21, bg=COLS$d_j, col=COLS$d_j2, data=plt6)
        points(d_a_eigPolyViable ~ d_a, pch=21, bg=COLS$d_a, col=COLS$d_a2, data=plt6)
        points(d_g_eigPolyViable ~ d_g, pch=21, bg=COLS$d_g, col=COLS$d_g2, data=plt6)
        # axes        
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.03, 1.075, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.1, substitute(paste(italic(C), " = ", cc), list(cc = pars6$C)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)

}





df1  <-  "deltaSelfingSimPolySpace_sMax0.15_nSamples1000_dStar0.8_hf0.5_hm0.5_f7"
df2  <-  "deltaSelfingSimPolySpace_sMax0.15_nSamples1000_dStar0.8_hf0.25_hm0.25_f7"
deltaSelfingLoadPolySpaceFig  <-  function(df1, df2) {

# Make filenames for import from df names
    fName1  <-  paste('./output/simData/', df1, '.csv', sep="")
    fName2  <-  paste('./output/simData/', df2, '.csv', sep="")

# Extract plotting parameter values from df names
    d1   <-  strsplit(df1, '_')[[1]][c(2,4:7)]
    pars1  <-  list(
                    "sMax"   =  as.numeric(strsplit(d1[1],'x')[[1]][2]),
                    "dStar"  =  as.numeric(strsplit(d1[2],'r')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d1[3],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d1[4],'m')[[1]][2]),
                    "f"      =  as.numeric(strsplit(d1[5],'f')[[1]][2])
                    )
    d2   <-  strsplit(df2, '_')[[1]][c(2,4:7)]
    pars2  <-  list(
                    "sMax"   =  as.numeric(strsplit(d2[1],'x')[[1]][2]),
                    "dStar"  =  as.numeric(strsplit(d2[2],'r')[[1]][2]),
                    "hf"     =  as.numeric(strsplit(d2[3],'f')[[1]][2]),
                    "hm"     =  as.numeric(strsplit(d2[4],'m')[[1]][2]),
                    "f"      =  as.numeric(strsplit(d2[5],'f')[[1]][2])
                    )

# Import data sets
    plt1  <-  read.csv(fName1, head=TRUE)
    plt2  <-  read.csv(fName2, head=TRUE)

# Color scheme
    COLS  <-  list(
                    "PG"     =  transparentColor('#252525', opacity=1),
                    "dSim"   =  transparentColor('#252525', opacity=0.6),
                    "d_j"    =  transparentColor('dodgerblue4', opacity=0.6),
                    "d_a"    =  transparentColor('dodgerblue', opacity=0.6),
                    "d_g"    =  transparentColor('tomato', opacity=0.6),
                    "dSim2"  =  transparentColor('#252525', opacity=1),
                    "d_j2"   =  transparentColor('dodgerblue4', opacity=1),
                    "d_a2"   =  transparentColor('dodgerblue', opacity=1),
                    "d_g2"   =  transparentColor('tomato', opacity=1)
                    )
#  Create vector of delta values for pop gen predictions.
    dStar  <-  pars1$dStar
    CLine  <-  seq(0,0.9,length=100)
    dLine  <-  predDelta(dStar=dStar, b=1/2, a=0.2, C=CLine) 


# Set plot layout
    layout.mat  <- matrix(c(1:2), nrow=2, ncol=1, byrow=TRUE)
    layout      <- layout(layout.mat,respect=TRUE)

    ## Panel A: Additive SA (hf = hm = 1/2)
        PGSpace     <-  c()
        for(i in 1:length(CLine)) {
            PGSpace[i]     <-  popGen_PolySpace_Delta_Add(C=CLine[i], delta=dLine[i], sMax=pars1$sMax)
        }
        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.925), ylim = c(0,0.12), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(PGSpace ~ CLine, lwd=2, col=COLS$PG)
        points(d_eigPolyViable   ~ C, pch=21, bg=COLS$dSim, col=COLS$dSim2, data=plt1)
        points(d_j_eigPolyViable ~ C, pch=21, bg=COLS$d_j, col=COLS$d_j2, data=plt1)
        points(d_a_eigPolyViable ~ C, pch=21, bg=COLS$d_a, col=COLS$d_a2, data=plt1)
        points(d_g_eigPolyViable ~ C, pch=21, bg=COLS$d_g, col=COLS$d_g2, data=plt1)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.1, expression(paste(italic(h), " = ", 1/2)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(-0.35, -0.15, expression(paste("Proportion of parameter space")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)

      #Legend
        legend( x       =  usr[2]*0.95,
                y       =  usr[4],
                legend  =  c(
                             expression(paste("Pop. Gen. ", italic(delta))),
                             expression(paste(delta)),
                             expression(paste(delta[italic(j)])),
                             expression(paste(delta[italic(a)])),
                             expression(paste(delta[gamma]))),
                 lty     =  c(1,NA,NA,NA,NA),
                 lwd     =  c(2,NA,NA,NA,NA),
                 col     =  c(COLS$PG,
                              COLS$dSim,
                              COLS$d_j,
                              COLS$d_a,
                              COLS$d_g),
                 pch     =  c(NA,21,21,21,21),
                 pt.bg   =  c(NA,
                              COLS$dSim,
                              COLS$d_j,
                              COLS$d_a,
                              COLS$d_g),
                 cex     =  0.75,
                 pt.cex  =  0.75,
                 xjust   =  1,
                 yjust   =  1,
                 bty     =  'n',
                 border  =  NA)
    ## Panel B: Dominance Reversal SA (hf = hm = 1/4)
        CLine2  <-  CLine[-1]
        dLine2  <-  dLine[-1]
        PGSpace     <-  c()
        for(i in 1:length(CLine2)) {
            PGSpace[i]     <-  popGen_PolySpace_Delta_DomRev(C=CLine2[i], delta=dLine2[i], sMax=pars2$sMax)
        }
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='', xlim = c(0,0.925), ylim = c(0,0.8), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(PGSpace ~ CLine2, lwd=2, col=COLS$PG)
        points(d_eigPolyViable ~ C, pch=21, bg=COLS$dSim, col=COLS$dSim2, data=plt2)
        points(d_j_eigPolyViable ~ C, pch=21, bg=COLS$d_j, col=COLS$d_j2, data=plt2)
        points(d_a_eigPolyViable ~ C, pch=21, bg=COLS$d_a, col=COLS$d_a2, data=plt2)
        points(d_g_eigPolyViable ~ C, pch=21, bg=COLS$d_g, col=COLS$d_g2, data=plt2)
        # axes        
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(0.03, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.1, expression(paste(italic(h), " = ", 1/4)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.5, -0.35, expression(paste("Selfing Rate (", italic(C), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)

}









# Simple inv6 Fig
MimulusInv6Fig  <-  function() {

    sfs  <-  seq(0,1, by=0.001)
    sms  <-  seq(0,1, by=0.001)

    hf  <-  1/2
    hm  <-  0.35

    sfInv6  <-  0.308
    smInv6  <-  0.88

    inv_A  <-  popGen_A_invade(hf=hf, hm=hm, sm=sms, C=0)
    inv_a  <-  popGen_a_invade(hf=hf, hm=hm, sm=sms, C=0)
    inv_a[inv_a > 1]  <-  1

# Color scheme
    COLS  <-  list(
                    "PG"     =  transparentColor('#252525', opacity=1),
                    "fill"   =  transparentColor('#252525', opacity=0.25),
                    "inv6"   =  transparentColor('tomato', opacity=1),
                    "inv6bg" =  transparentColor('tomato', opacity=0.6)
                    )
# Set plot layout
    layout.mat  <- matrix(c(1:2), nrow=1, ncol=2, byrow=TRUE)
    layout      <- layout(layout.mat,respect=TRUE)

    ## Panel A: Outcrossing

        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1.05), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        polygon(c(rev(sms),sms), c(rev(inv_A), inv_a), col=transparentColor('grey80', 0.6), border=NA)
        lines(inv_A ~ sms, lwd=2, col=COLS$PG)
        lines(inv_a[1:436] ~ sms[1:436], lwd=2, col=COLS$PG)
        points(sfInv6 ~ smInv6, pch=21, col=COLS$inv6, bg=COLS$inv6bg, cex=1.25)
        # axes        
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel((smInv6-0.05), (sfInv6-0.07), expression(italic('inv6')), cex=1, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.1, expression(paste(italic(C), " = ", 0)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(-0.35, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.35, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)


    ## Panel B: partial selfing

    inv_A  <-  popGen_A_invade(hf=hf, hm=hm, sm=sms, C=0.24) #, C=0.38)
    inv_a  <-  popGen_a_invade(hf=hf, hm=hm, sm=sms, C=0.24) #, C=0.38)
    inv_a[inv_a > 1]  <-  1

        # Make the plot
        plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1.05), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        polygon(c(rev(sms),sms), c(rev(inv_A), inv_a), col=transparentColor('grey80', 0.6), border=NA)
        lines(inv_A ~ sms, lwd=2, col=COLS$PG)
        lines(inv_a[1:572] ~ sms[1:572], lwd=2, col=COLS$PG)
        points(sfInv6 ~ smInv6, pch=21, col=COLS$inv6, bg=COLS$inv6bg, cex=1.25)
        # axes        
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.03, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.1, expression(paste(italic(C), " = ", 0.24)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.5, -0.35, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel((smInv6-0.05), (sfInv6-0.07), expression(italic('inv6')), cex=1, adj=c(0.5, 0.5), xpd=NA)



}


