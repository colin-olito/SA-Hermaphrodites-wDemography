###############
# DEPENDENCIES
###############
library(extrafont)
library(fontcm)
loadfonts(quiet=TRUE)
library(plyr)
library(lattice)
library(latticeExtra)
library(wesanderson)
library(MASS)
library(raster)
library(akima)

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



proportionalArrows <- function(px1, py1, px2, py2, adj=c(0, 1), log=FALSE, length=length, ...) {
    usr  <-  par('usr')
    x.p1  <-  usr[1] + px1*(usr[2] - usr[1])
    y.p1  <-  usr[3] + py1*(usr[4] - usr[3])
    x.p2  <-  usr[1] + px2*(usr[2] - usr[1])
    y.p2  <-  usr[3] + py2*(usr[4] - usr[3])
    if(log=='x') {
        x.p1  <-  10^(x.p1)
        x.p2  <-  10^(x.p2)
    }
    if(log=='y') {
        y.p1  <-  10^(y.p1)
        y.p2  <-  10^(y.p2)
    }
    if(log=='xy') {
        x.p1  <-  10^(x.p1)
        y.p1  <-  10^(y.p1)
        x.p2  <-  10^(x.p2)
        y.p2  <-  10^(y.p2)
    }
    arrows(x0=x.p1, y0=y.p1, x1=x.p2, y1=y.p2, length=length,...)
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
##  Final Figs for Paper
##############################################################
##############################################################

demViablePolySpaceFig  <-  function() {

    # Import data sets
    invA       <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.5_hm0.5_C0_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extA_f5.8  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.5_hm0.5_C0_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extA_f6.0  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.5_hm0.5_C0_delta0_dj0_da0_dg0_f6.csv", head=TRUE)
    extA_f6.2  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.5_hm0.5_C0_delta0_dj0_da0_dg0_f6.2.csv", head=TRUE)

    invB       <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.5_hm0.5_C0.25_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extB_f5.8  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.5_hm0.5_C0.25_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extB_f6.0  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.5_hm0.5_C0.25_delta0_dj0_da0_dg0_f6.csv", head=TRUE)
    extB_f6.2  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.5_hm0.5_C0.25_delta0_dj0_da0_dg0_f6.2.csv", head=TRUE)

    invC       <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.5_hm0.5_C0.5_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extC_f5.8  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.5_hm0.5_C0.5_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extC_f6.0  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.5_hm0.5_C0.5_delta0_dj0_da0_dg0_f6.csv", head=TRUE)
    extC_f6.2  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.5_hm0.5_C0.5_delta0_dj0_da0_dg0_f6.2.csv", head=TRUE)

    invD       <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.25_hm0.25_C0_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extD_f5.8  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.25_C0_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extD_f5.9  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.25_C0_delta0_dj0_da0_dg0_f5.9.csv", head=TRUE)
    extD_f6.0  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.25_C0_delta0_dj0_da0_dg0_f6.csv", head=TRUE)

    invE       <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.25_hm0.25_C0.25_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extE_f5.8  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.25_C0.25_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extE_f5.9  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.25_C0.25_delta0_dj0_da0_dg0_f5.9.csv", head=TRUE)
    extE_f6.0  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.25_C0.25_delta0_dj0_da0_dg0_f6.csv", head=TRUE)

    invF       <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.25_hm0.25_C0.5_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extF_f5.8  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.25_C0.5_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extF_f5.9  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.25_C0.5_delta0_dj0_da0_dg0_f5.9.csv", head=TRUE)
    extF_f6.0  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.25_C0.5_delta0_dj0_da0_dg0_f6.csv", head=TRUE)

    # Clean up a couple missing values from data
    extA_f5.8$sfThreshold[16]  <-  mean(extA_f5.8$sfThreshold[15], extA_f5.8$sfThreshold[17])
    extA_f6.0$sfThreshold[87]  <-  mean(extA_f6.0$sfThreshold[86], extA_f6.0$sfThreshold[88])
    extA_f6.2$sfThreshold[88]  <-  mean(extA_f6.2$sfThreshold[87], extA_f6.2$sfThreshold[89])   
    extA_f6.2$sfThreshold[95]  <-  mean(extA_f6.2$sfThreshold[94], extA_f6.2$sfThreshold[96])   
    extA_f6.2$sfThreshold[56]  <-  extA_f6.2$sfThreshold[57]
    extD_f5.8$sfThreshold[33]  <-  mean(extD_f5.8$sfThreshold[32], extD_f5.8$sfThreshold[34])
    extD_f5.8$sfThreshold[37]  <-  mean(extD_f5.8$sfThreshold[36], extD_f5.8$sfThreshold[38])
    extD_f5.8$sf[11]           <-  mean(extD_f5.8$sf[11], extD_f5.8$sfThreshold[21])
    extD_f5.8$smThreshold[11]  <-  mean(extD_f5.8$smThreshold[12], extD_f5.8$sms[21])
    extE_f5.8$sfThreshold[67]  <-  mean(extE_f5.8$sfThreshold[66], extE_f5.8$sfThreshold[68])
    extE_f5.8$sf[9]            <-  mean(extE_f5.8$sf[9], extE_f5.8$sfThreshold[29])
    extE_f5.8$smThreshold[11]  <-  mean(extE_f5.8$smThreshold[10], extE_f5.8$sms[29])
    extE_f5.9$sfThreshold[87]  <-  mean(extE_f5.9$sfThreshold[86], extE_f5.9$sfThreshold[88])
    extD_f6.0$sfThreshold[83]  <-  (extD_f6.0$sfThreshold[82] + extD_f6.0$sfThreshold[84])/2

# Color scheme
    COLS  <-  list("line"     =  transparentColor('#252525', opacity=1),
                   "extinct"  =  transparentColor('red', opacity=0.15))

# Set plot layout
layout.mat <- matrix(c(1:6), nrow=2, ncol=3, byrow=TRUE)
layout     <- layout(layout.mat,respect=TRUE)

##  Row 1: Additive allelic effects

    ##  Panel A: C = 0; hf = hm = 1/2
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s',xpd=TRUE)
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Extinction Thresholds
        # f = 5.8 polygon
            t    <-  cbind(extA_f5.8$sf[!is.na(extA_f5.8$smThreshold)], extA_f5.8$smThreshold[!is.na(extA_f5.8$smThreshold)])
            t    <-  rbind(t, c(0.15, 0.15))
            b    <-  cbind(extA_f5.8$sms[!is.na(extA_f5.8$sfThreshold)], extA_f5.8$sfThreshold[!is.na(extA_f5.8$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extA_f5.8)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extA_f5.8)   
        proportionalLabel(0.85, 0.14, expression(paste(italic(f), " = 5.8")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # f = 6.0 polygon
            t    <-  cbind(extA_f6.0$sf[!is.na(extA_f6.0$smThreshold)], extA_f6.0$smThreshold[!is.na(extA_f6.0$smThreshold)])
            t    <-  rbind(t, c(0.15, 0.15))
            b    <-  cbind(extA_f6.0$sms[!is.na(extA_f6.0$sfThreshold)], extA_f6.0$sfThreshold[!is.na(extA_f6.0$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extA_f6.0)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extA_f6.0)   
        proportionalLabel(0.91, 0.35, expression(paste("6.0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # f = 6.2 polygon
            t    <-  cbind(extA_f6.2$sf[!is.na(extA_f6.2$smThreshold)], extA_f6.2$smThreshold[!is.na(extA_f6.2$smThreshold)])
            t    <-  rbind(t, c(0.15, 0.15))
            b    <-  cbind(extA_f6.2$sms[!is.na(extA_f6.2$sfThreshold)], extA_f6.2$sfThreshold[!is.na(extA_f6.2$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extA_f6.2)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extA_f6.2)   
        proportionalLabel(0.91, 0.54, expression(paste("6.2")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # Invasion Boundaries
        lines(aInvBound[1:88] ~ sms[1:88], lty=1, lwd=1.5, col=COLS$line, data=invA)
        lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invA)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.5, 0.5, expression(paste(italic(h), " = 1/2")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.3, -0.2, expression(paste("Selection on female reproductive function (", italic(s[f]), ")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
#        proportionalLabel(-0.3, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = 0')), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    ##  Panel B: C = 1/4; hf = hm = 1/2
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2, bty='L')
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Extinction Thresholds
        # f = 5.8 polygon
            t    <-  cbind(extB_f5.8$sf[!is.na(extB_f5.8$smThreshold)], extB_f5.8$smThreshold[!is.na(extB_f5.8$smThreshold)])
            b    <-  cbind(extB_f5.8$sms[!is.na(extB_f5.8$sfThreshold)], extB_f5.8$sfThreshold[!is.na(extB_f5.8$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extB_f5.8)
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extB_f5.8)
        proportionalLabel(0.85, 0.14, expression(paste(italic(f), " = 5.8")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # f = 6.0 polygon
            t    <-  cbind(extB_f6.0$sf[!is.na(extB_f6.0$smThreshold)], extB_f6.0$smThreshold[!is.na(extB_f6.0$smThreshold)])
            b    <-  cbind(extB_f6.0$sms[!is.na(extB_f6.0$sfThreshold)], extB_f6.0$sfThreshold[!is.na(extB_f6.0$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extB_f6.0)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extB_f6.0)   
        proportionalLabel(0.91, 0.35, expression(paste("6.0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # f = 6.2 polygon
            t    <-  cbind(extB_f6.2$sf[!is.na(extB_f6.2$smThreshold)], extB_f6.2$smThreshold[!is.na(extB_f6.2$smThreshold)])
            b    <-  cbind(extB_f6.2$sms[!is.na(extB_f6.2$sfThreshold)], extB_f6.2$sfThreshold[!is.na(extB_f6.2$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extB_f6.2)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extB_f6.2)   
#        proportionalLabel(0.9, 0.44, expression(paste("6.2")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(1.1, 0.525, expression(paste("6.2")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalArrows(px1=1.1, py1=0.525, px2=0.97, py2=0.525, length=0.03)
        # Invasion Boundaries
        lines(aInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invB)
        lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invB)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = 1/4')), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    ##  Panel C: C = 1/2; hf = hm = 1/2
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2, bty='L')
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Extinction Thresholds
        # f = 5.8 polygon
            t    <-  cbind(extC_f5.8$sf[!is.na(extC_f5.8$smThreshold)], extC_f5.8$smThreshold[!is.na(extC_f5.8$smThreshold)])
            b    <-  cbind(extC_f5.8$sms[!is.na(extC_f5.8$sfThreshold)], extC_f5.8$sfThreshold[!is.na(extC_f5.8$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extC_f5.8)
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extC_f5.8)
        proportionalLabel(0.85, 0.14, expression(paste(italic(f), " = 5.8")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # f = 6.0 polygon
            t    <-  cbind(extC_f6.0$sf[!is.na(extC_f6.0$smThreshold)], extC_f6.0$smThreshold[!is.na(extC_f6.0$smThreshold)])
            b    <-  cbind(extC_f6.0$sms[!is.na(extC_f6.0$sfThreshold)], extC_f6.0$sfThreshold[!is.na(extC_f6.0$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extC_f6.0)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extC_f6.0)   
#        proportionalLabel(0.85, 0.22, expression(paste("6.0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
#        proportionalArrows(0.875, 0.26, 0.95, 0.3, length=0.051, lwd=1)
        proportionalLabel(1.1, 0.32, expression(paste("6.0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalArrows(px1=1.1, py1=0.32, px2=0.97, py2=0.32, length=0.03)
        # Invasion Boundaries
        lines(aInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invC)
        lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invC)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = 1/2')), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.075, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)


##  Row 2: Dominance Reversal
  ##  Panel D: C = 0; hf = hm = 1/4
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2, bty='L')
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Extinction Thresholds
        # f = 5.8 polygon
            t    <-  cbind(extD_f5.8$sf[!is.na(extD_f5.8$smThreshold)], extD_f5.8$smThreshold[!is.na(extD_f5.8$smThreshold)])
            t    <-  rbind(t, c(0.15, 0.15))
            b    <-  cbind(extD_f5.8$sms[!is.na(extD_f5.8$sfThreshold)], extD_f5.8$sfThreshold[!is.na(extD_f5.8$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extD_f5.8)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extD_f5.8)   
        proportionalLabel(0.45, 0.8, expression(paste(italic(f), " = 5.8")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=75)
        # f = 5.9 polygon
            t    <-  cbind(extD_f5.9$sf[!is.na(extD_f5.9$smThreshold)], extD_f5.9$smThreshold[!is.na(extD_f5.9$smThreshold)])
            t    <-  rbind(t, c(0.15, 0.15))
            b    <-  cbind(extD_f5.9$sms[!is.na(extD_f5.9$sfThreshold)], extD_f5.9$sfThreshold[!is.na(extD_f5.9$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extD_f5.9)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extD_f5.9)   
#        proportionalLabel(0.85, 0.25, expression(paste("5.9")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.64, 0.86, expression(paste("5.9")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=75)
        # f = 6.0 polygon
            t    <-  cbind(extD_f6.0$sf[!is.na(extD_f6.0$smThreshold)], extD_f6.0$smThreshold[!is.na(extD_f6.0$smThreshold)])
            t    <-  rbind(t, c(0.15, 0.15))
            b    <-  cbind(extD_f6.0$sms, extD_f6.0$sfThreshold)
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extD_f6.0)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extD_f6.0)   
#        proportionalLabel(0.85, 0.7, expression(paste("6.0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.83, 0.86, expression(paste("6.0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=75)
        # Invasion Boundaries
        lines(aInvBound[1:33] ~ sms[1:33], lty=1, lwd=1.5, col=COLS$line, data=invD)
        lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invD)
        # axes        
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.5, 0.5, expression(paste(italic(h), " = 1/4")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
#        proportionalLabel(-0.3, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.04, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)


  ##  Panel E: C = 1/4; hf = hm = 1/4
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2, bty='L')
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Extinction Thresholds
        # f = 5.8 polygon
            t    <-  cbind(extE_f5.8$sf[!is.na(extE_f5.8$smThreshold)], extE_f5.8$smThreshold[!is.na(extE_f5.8$smThreshold)])
            t    <-  rbind(t, c(0.15, 0.15))
            b    <-  cbind(extE_f5.8$sms[!is.na(extE_f5.8$sfThreshold)], extE_f5.8$sfThreshold[!is.na(extE_f5.8$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extE_f5.8)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extE_f5.8)   
#        proportionalLabel(0.85, 0.14, expression(paste(italic(f), " = 5.8")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.8, 0.85, expression(paste(italic(f), " = 5.8")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=60)
        # f = 5.9 polygon
            t    <-  cbind(extE_f5.9$sf[!is.na(extE_f5.9$smThreshold)], extE_f5.9$smThreshold[!is.na(extE_f5.9$smThreshold)])
            b    <-  cbind(extE_f5.9$sms[!is.na(extE_f5.9$sfThreshold)], extE_f5.9$sfThreshold[!is.na(extE_f5.9$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extE_f5.9)
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extE_f5.9)
#        proportionalLabel(0.85, 0.33, expression(paste("5.9")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.9, 0.63, expression(paste("5.9")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=60)
        # Invasion Boundaries
        lines(aInvBound[1:68] ~ sms[1:68], lty=1, lwd=1.5, col=COLS$line, data=invE)
        lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invE)
        # axes        
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.04, 1.075, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.3, expression(paste("Selection on male reproductive function (", italic(s[m]), ")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
#        proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)

  ##  Panel F: C = 1/2; hf = hm = 1/4
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2, bty='L')
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Extinction Thresholds
        # f = 5.8 polygon
            t    <-  cbind(extF_f5.8$sf[!is.na(extF_f5.8$smThreshold)], extF_f5.8$smThreshold[!is.na(extF_f5.8$smThreshold)])
            t    <-  rbind(t, c(0.15, 0.15))
            b    <-  cbind(extF_f5.8$sms[!is.na(extF_f5.8$sfThreshold)], extF_f5.8$sfThreshold[!is.na(extF_f5.8$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extF_f5.8)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extF_f5.8)   
#        proportionalLabel(0.85, 0.14, expression(paste(italic(f), " = 5.8")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.85, 0.36, expression(paste(italic(f), " = 5.8")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=30)
        # f = 5.9 polygon
            t    <-  cbind(extF_f5.9$sf[!is.na(extF_f5.9$smThreshold)], extF_f5.9$smThreshold[!is.na(extF_f5.9$smThreshold)])
            b    <-  cbind(extF_f5.9$sms[!is.na(extF_f5.9$sfThreshold)], extF_f5.9$sfThreshold[!is.na(extF_f5.9$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extF_f5.9)
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extF_f5.9)
#        proportionalLabel(0.85, 0.33, expression(paste("5.9")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(1.1, 0.24, expression(paste("5.9")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalArrows(px1=1.1, py1=0.24, px2=0.97, py2=0.24, length=0.03)
        # Invasion Boundaries
        lines(aInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invF)
        lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invF)
        # axes        
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.04, 1.075, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
}

# Function to create color scalebar for image() heatmap
# from https://menugget.blogspot.com/2011/08/adding-scale-to-image-plot.html#more
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
 if(!missing(breaks)){
  if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
 }
 if(missing(breaks) & !missing(zlim)){
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
 }
 if(missing(breaks) & missing(zlim)){
  zlim <- range(z, na.rm=TRUE)
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 poly <- vector(mode="list", length(col))
 for(i in seq(poly)){
  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
 }
 xaxt <- ifelse(horiz, "s", "n")
 yaxt <- ifelse(horiz, "n", "s")
 if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
 if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
 if(missing(xlim)) xlim=XLIM
 if(missing(ylim)) ylim=YLIM
 plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
 for(i in seq(poly)){
  if(horiz){
   polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
  }
  if(!horiz){
   polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
  }
 }
}


#' Revised Fig. 1 showing heatmap of lambda 
#' instead of extinction thresholds

lambdaPolySpaceFig  <-  function() {
 
    # Import data sets
    invA        <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.5_hm0.5_C0_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    heatA_f5.8  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_alpha0_hf0.5_hm0.5_C0_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    heatA_f6.0  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_alpha0_hf0.5_hm0.5_C0_delta0_dj0_da0_dg0_f6.csv", head=TRUE)
    heatA_f6.2  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_alpha0_hf0.5_hm0.5_C0_delta0_dj0_da0_dg0_f6.2.csv", head=TRUE)

    invB        <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.5_hm0.5_C0.25_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    heatB_f5.8  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.5_hm0.5_C0.25_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    heatB_f6.0  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.5_hm0.5_C0.25_delta0_dj0_da0_dg0_f6.csv", head=TRUE)
    heatB_f6.2  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.5_hm0.5_C0.25_delta0_dj0_da0_dg0_f6.2.csv", head=TRUE)

    invC        <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.5_hm0.5_C0.5_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    heatC_f5.8  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.5_hm0.5_C0.5_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    heatC_f6.0  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.5_hm0.5_C0.5_delta0_dj0_da0_dg0_f6.csv", head=TRUE)
    heatC_f6.2  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.5_hm0.5_C0.5_delta0_dj0_da0_dg0_f6.2.csv", head=TRUE)

    invD        <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.25_hm0.25_C0_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    heatD_f5.8  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.25_hm0.25_C0_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    heatD_f5.9  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.25_hm0.25_C0_delta0_dj0_da0_dg0_f5.9.csv", head=TRUE)
    heatD_f6.0  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.25_hm0.25_C0_delta0_dj0_da0_dg0_f6.csv", head=TRUE)

    invE        <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.25_hm0.25_C0.25_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    heatE_f5.8  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.25_hm0.25_C0.25_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    heatE_f5.9  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.25_hm0.25_C0.25_delta0_dj0_da0_dg0_f5.9.csv", head=TRUE)
    heatE_f6.0  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.25_hm0.25_C0.25_delta0_dj0_da0_dg0_f6.csv", head=TRUE)

    invF        <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.25_hm0.25_C0.5_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    heatF_f5.8  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.25_hm0.25_C0.5_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    heatF_f5.9  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.25_hm0.25_C0.5_delta0_dj0_da0_dg0_f5.9.csv", head=TRUE)
    heatF_f6.0  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.25_hm0.25_C0.5_delta0_dj0_da0_dg0_f6.csv", head=TRUE)


    # Color scheme
    COLS  <-  list("line"     =  transparentColor('#252525', opacity=1),
                   "extinct"  =  transparentColor('red', opacity=0.15))

    
    # heatmap colors modified from wes_anderson palette "Zissou1"
    Z1 = rev(c("dodgerblue4","dodgerblue4","#3B9AB2","#3B9AB2", "#EBCC2A", "#EBCC2A", "#E1AF00", "#E1AF00", "#F21A00", "#F21A00"))
    blu     <-  colorRampPalette(Z1[length(Z1)])(2)
    notBlu  <-  colorRampPalette(Z1, interpolate = c("linear"))(60)
    HEAT    <-  c(notBlu, blu)

    # resolution for heatmaps
    resolution = (heatA_f6.2$sf[2] - heatA_f6.2$sf[1])
    # breaks for heatmap colors
    breaks <- seq(min(c(heatA_f5.8$lambda_sim,
                      heatB_f5.8$lambda_sim,
                      heatC_f5.8$lambda_sim,
                      heatD_f5.8$lambda_sim,
                      heatE_f5.8$lambda_sim,
                      heatF_f5.8$lambda_sim)), 
                max(c(heatA_f6.2$lambda_sim,
                      heatB_f6.2$lambda_sim,
                      heatC_f6.2$lambda_sim,
                      heatD_f6.0$lambda_sim,
                      heatE_f6.0$lambda_sim,
                      heatF_f6.0$lambda_sim)),length.out=(length(HEAT)+1))
#    image.scale(A_f6.2, col=HEAT, breaks=breaks, horiz=FALSE, yaxt="n")
#    axis(4, las=2)

    # Set plot layout
    layout.mat <- matrix(c(
                           1,1, 2,2, 3,3, 10,
                           1,1, 2,2, 3,3, 10,
                           4,4, 5,5, 6,6, 10,
                           4,4, 5,5, 6,6, 10,
                           7,7, 8,8, 9,9, 10,
                           7,7, 8,8, 9,9, 10,
                           11,11, 11,11, 11,11, 11,
                           12,12, 13,13, 14,14, 21,
                           12,12, 13,13, 14,14, 21,
                           15,15, 16,16, 17,17, 21,
                           15,15, 16,16, 17,17, 21,
                           18,18, 19,19, 20,20, 21,
                           18,18, 19,19, 20,20, 21), 
                        nrow=13, ncol=7, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)



##  Row 1: Fertility = 6.2

    # Panel 1: h = 1/2,  C = 0, , f = 6.2
    A_f6.2 <- interp(x=heatA_f6.2$sm, y=heatA_f6.2$sf, z=heatA_f6.2$lambda_sim, 
                     xo=seq(min(heatA_f6.2$sm),max(heatA_f6.2$sm),by=resolution), 
                     yo=seq(min(heatA_f6.2$sf),max(heatA_f6.2$sf),by=resolution), duplicate="mean")
    # heatmap of lambda
    par(omi=c(1, 1.5, 0.5, 1.25), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s', xpd=TRUE)
    image(A_f6.2, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("0","0.05","0.10", "0.15"), cex=1.2)
    # Invasion Boundaries
    lines(aInvBound[1:88] ~ sms[1:88], lty=1, lwd=1.5, col=COLS$line, data=invA)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invA)
    # Labels/Annotations
    proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = 0')), cex=1.75, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.04, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(-0.5, 0.5, expression(paste(italic(f), " = 6.2")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)

    # Panel 2: h = 1/2,  C = 1/4, f = 6.2
    B_f6.2 <- interp(x=heatB_f6.2$sm, y=heatB_f6.2$sf, z=heatB_f6.2$lambda_sim, 
                xo=seq(min(heatB_f6.2$sm),max(heatB_f6.2$sm),by=resolution), 
                yo=seq(min(heatB_f6.2$sf),max(heatB_f6.2$sf),by=resolution), duplicate="mean")
    image(B_f6.2, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    # Invasion Boundaries
    lines(aInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invB)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invB)
    # Labels/Annotations
#    proportionalLabel(0.5, 1.35, expression(paste(italic(h), " = 1/2")), cex=2.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
    proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = 1/4')), cex=1.75, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.04, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # Panel 3: h = 1/2,  C = 1/2, f = 6.2
    C_f6.2 <- interp(x=heatC_f6.2$sm, y=heatC_f6.2$sf, z=heatC_f6.2$lambda_sim, 
                xo=seq(min(heatC_f6.2$sm),max(heatC_f6.2$sm),by=resolution), 
                yo=seq(min(heatC_f6.2$sf),max(heatC_f6.2$sf),by=resolution), duplicate="mean")
    image(C_f6.2, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    # Invasion Boundaries
    lines(aInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invC)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invC)
    # Labels/Annotations
    proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = 1/2')), cex=1.75, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.04, 1.075, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

##  Row 2: Fertility = 6.0

    # Panel 4: h = 1/2,  C = 0, , f = 6.0
    A_f6.0 <- interp(x=heatA_f6.0$sm, y=heatA_f6.0$sf, z=heatA_f6.0$lambda_sim, 
                     xo=seq(min(heatA_f6.0$sm),max(heatA_f6.0$sm),by=resolution), 
                     yo=seq(min(heatA_f6.0$sf),max(heatA_f6.0$sf),by=resolution), duplicate="mean")
    # heatmap of lambda
    image(A_f6.0, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("0","0.05","0.10", "0.15"), cex.axis=1.2)
    # Invasion Boundaries
    lines(aInvBound[1:88] ~ sms[1:88], lty=1, lwd=1.5, col=COLS$line, data=invA)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invA)
    # Labels/Annotations
    proportionalLabel(0.04, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(-0.3, 0.5, expression(paste("Selection on female reproductive function (", italic(s[f]), ")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
    proportionalLabel(-0.75, 0.5, expression(paste(italic(h), " = 1/2")), cex=2.0, adj=c(0.5, 0.5), xpd=NA, srt=90)
    proportionalLabel(-0.5, 0.5, expression(paste(italic(f), " = 6.0")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)

    # Panel 5: h = 1/2,  C = 1/4, f = 6.0
    B_f6.0 <- interp(x=heatB_f6.0$sm, y=heatB_f6.0$sf, z=heatB_f6.0$lambda_sim, 
                xo=seq(min(heatB_f6.0$sm),max(heatB_f6.0$sm),by=resolution), 
                yo=seq(min(heatB_f6.0$sf),max(heatB_f6.0$sf),by=resolution), duplicate="mean")
    image(B_f6.0, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    # Invasion Boundaries
    lines(aInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invB)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invB)
    # Labels/Annotations
    proportionalLabel(0.04, 1.075, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # Panel 6: h = 1/2,  C = 1/2, f = 6.0
    C_f6.0 <- interp(x=heatC_f6.0$sm, y=heatC_f6.0$sf, z=heatC_f6.0$lambda_sim, 
                xo=seq(min(heatC_f6.0$sm),max(heatC_f6.0$sm),by=resolution), 
                yo=seq(min(heatC_f6.0$sf),max(heatC_f6.0$sf),by=resolution), duplicate="mean")
    image(C_f6.0, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    # Invasion Boundaries
    lines(aInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invC)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invC)
    # Labels/Annotations
    proportionalLabel(0.04, 1.075, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

##  Row 3: Fertility = 5.8

    # Panel 7: h = 1/2,  C = 0, , f = 5.8
    A_f5.8 <- interp(x=heatA_f5.8$sm, y=heatA_f5.8$sf, z=heatA_f5.8$lambda_sim, 
                     xo=seq(min(heatA_f5.8$sm),max(heatA_f5.8$sm),by=resolution), 
                     yo=seq(min(heatA_f5.8$sf),max(heatA_f5.8$sf),by=resolution), duplicate="mean")
    # heatmap of lambda
    image(A_f5.8, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("0","0.05","0.10", "0.15"), cex.axis=1.2)
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("0","0.05","0.10", "0.15"), cex.axis=1.2)
    # Invasion Boundaries
    lines(aInvBound[1:88] ~ sms[1:88], lty=1, lwd=1.5, col=COLS$line, data=invA)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invA)
    # Labels/Annotations
    proportionalLabel(0.04, 1.075, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(-0.5, 0.5, expression(paste(italic(f), " = 5.8")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)

    # Panel 8: h = 1/2,  C = 1/4, f = 5.8
    B_f5.8 <- interp(x=heatB_f5.8$sm, y=heatB_f5.8$sf, z=heatB_f5.8$lambda_sim, 
                xo=seq(min(heatB_f5.8$sm),max(heatB_f5.8$sm),by=resolution), 
                yo=seq(min(heatB_f5.8$sf),max(heatB_f5.8$sf),by=resolution), duplicate="mean")
    image(B_f5.8, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("0","0.05","0.10", "0.15"), cex.axis=1.2)
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    # Invasion Boundaries
    lines(aInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invB)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invB)
    # Labels/Annotations
    proportionalLabel(0.04, 1.075, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.5, -0.3, expression(paste("Selection on male reproductive function (", italic(s[m]), ")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)

    # Panel 9: h = 1/2,  C = 1/2, f = 5.8
    C_f5.8 <- interp(x=heatC_f5.8$sm, y=heatC_f5.8$sf, z=heatC_f5.8$lambda_sim, 
                xo=seq(min(heatC_f5.8$sm),max(heatC_f5.8$sm),by=resolution), 
                yo=seq(min(heatC_f5.8$sf),max(heatC_f5.8$sf),by=resolution), duplicate="mean")
    image(C_f5.8, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("0","0.05","0.10", "0.15"), cex.axis=1.2)
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    # Invasion Boundaries
    lines(aInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invC)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invC)
    # Labels/Annotations
    proportionalLabel(0.04, 1.075, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # Color ScaleBar
    image.scale(A_f6.2, col=HEAT, breaks=breaks, horiz=FALSE, yaxt="n")
    axis(4, las=2, cex.axis=1.5)
    proportionalLabel(2.5, 0.5, expression(paste("Population intrinsic growth rate (", italic(lambda), ")")), cex=2.5, adj=c(0.5, 0.5), xpd=NA, srt=270)


##  Blank Space between Plots
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')

##  Row 4: Fertility = 6.0

    # Panel 10: h = 1/2,  C = 0, , f = 6.0
    D_f6.0 <- interp(x=heatD_f6.0$sm, y=heatD_f6.0$sf, z=heatD_f6.0$lambda_sim, 
                     xo=seq(min(heatD_f6.0$sm),max(heatD_f6.0$sm),by=resolution), 
                     yo=seq(min(heatD_f6.0$sf),max(heatD_f6.0$sf),by=resolution), duplicate="mean")
    # heatmap of lambda
    image(D_f6.0, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("0","0.05","0.10", "0.15"), cex.axis=1.2)
    # Invasion Boundaries
    lines(aInvBound[1:33] ~ sms[1:33], lty=1, lwd=1.5, col=COLS$line, data=invD)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invD)
    # Labels/Annotations
    proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = 0')), cex=1.75, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.04, 1.075, 'J', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(-0.5, 0.5, expression(paste(italic(f), " = 6.0")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)

    # Panel 11: h = 1/2,  C = 1/4, f = 6.0
    E_f6.0 <- interp(x=heatE_f6.0$sm, y=heatE_f6.0$sf, z=heatE_f6.0$lambda_sim, 
                xo=seq(min(heatE_f6.0$sm),max(heatE_f6.0$sm),by=resolution), 
                yo=seq(min(heatE_f6.0$sf),max(heatE_f6.0$sf),by=resolution), duplicate="mean")
    image(E_f6.0, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    # Invasion Boundaries
    lines(aInvBound[1:68] ~ sms[1:68], lty=1, lwd=1.5, col=COLS$line, data=invE)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invE)
    # Labels/Annotations
#    proportionalLabel(0.5, 1.35, expression(paste(italic(h), " = 1/4")), cex=2.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
    proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = 1/4')), cex=1.75, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.04, 1.075, 'K', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # Panel 12: h = 1/2,  C = 1/2, f = 6.0
    F_f6.0 <- interp(x=heatF_f6.0$sm, y=heatF_f6.0$sf, z=heatF_f6.0$lambda_sim, 
                xo=seq(min(heatF_f6.0$sm),max(heatF_f6.0$sm),by=resolution), 
                yo=seq(min(heatF_f6.0$sf),max(heatF_f6.0$sf),by=resolution), duplicate="mean")
    image(F_f6.0, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    # Invasion Boundaries
    lines(aInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invF)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invF)
    # Labels/Annotations
    proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = 1/2')), cex=1.75, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.04, 1.075, 'L', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

##  Row 5: Fertility = 5.9
    # Panel 13: h = 1/4,  C = 0, , f = 6.0
    D_f5.9 <- interp(x=heatD_f5.9$sm, y=heatD_f5.9$sf, z=heatD_f5.9$lambda_sim, 
                     xo=seq(min(heatD_f5.9$sm),max(heatD_f5.9$sm),by=resolution), 
                     yo=seq(min(heatD_f5.9$sf),max(heatD_f5.9$sf),by=resolution), duplicate="mean")
    # heatmap of lambda
    image(D_f5.9, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("0","0.05","0.10", "0.15"), cex.axis=1.2)
    # Invasion Boundaries
    lines(aInvBound[1:33] ~ sms[1:33], lty=1, lwd=1.5, col=COLS$line, data=invD)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invD)
    # Labels/Annotations
    proportionalLabel(0.04, 1.075, 'M', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(-0.75, 0.5, expression(paste(italic(h), " = 1/4")), cex=2.0, adj=c(0.5, 0.5), xpd=NA, srt=90)
    proportionalLabel(-0.3, 0.5, expression(paste("Selection on female reproductive function (", italic(s[f]), ")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
    proportionalLabel(-0.5, 0.5, expression(paste(italic(f), " = 5.9")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)

    # Panel 14: h = 1/4,  C = 1/4, f = 5.9
    E_f5.9 <- interp(x=heatE_f5.9$sm, y=heatE_f5.9$sf, z=heatE_f5.9$lambda_sim, 
                xo=seq(min(heatE_f5.9$sm),max(heatE_f5.9$sm),by=resolution), 
                yo=seq(min(heatE_f5.9$sf),max(heatE_f5.9$sf),by=resolution), duplicate="mean")
    image(E_f5.9, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    # Invasion Boundaries
    lines(aInvBound[1:68] ~ sms[1:68], lty=1, lwd=1.5, col=COLS$line, data=invE)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invE)
    # Labels/Annotations
    proportionalLabel(0.04, 1.075, 'N', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # Panel 15: h = 1/2,  C = 1/2, f = 5.9
    F_f5.9 <- interp(x=heatF_f5.9$sm, y=heatF_f5.9$sf, z=heatF_f5.9$lambda_sim, 
                xo=seq(min(heatF_f5.9$sm),max(heatF_f5.9$sm),by=resolution), 
                yo=seq(min(heatF_f5.9$sf),max(heatF_f5.9$sf),by=resolution), duplicate="mean")
    image(F_f5.9, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    # Invasion Boundaries
    lines(aInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invF)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invF)
    # Labels/Annotations
    proportionalLabel(0.04, 1.075, 'O', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

##  Row 6: h = 1/4, Fertility = 6.0
    # Panel 16: h = 1/4,  C = 0, f = 6.0
    D_f5.8 <- interp(x=heatD_f5.8$sm, y=heatD_f5.8$sf, z=heatD_f5.8$lambda_sim, 
                     xo=seq(min(heatD_f5.8$sm),max(heatD_f5.8$sm),by=resolution), 
                     yo=seq(min(heatD_f5.8$sf),max(heatD_f5.8$sf),by=resolution), duplicate="mean")
    # heatmap of lambda
    image(D_f5.8, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("0","0.05","0.10", "0.15"), cex.axis=1.2)
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("0","0.05","0.10", "0.15"), cex.axis=1.2)
    # Invasion Boundaries
    lines(aInvBound[1:33] ~ sms[1:33], lty=1, lwd=1.5, col=COLS$line, data=invD)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invD)
    # Labels/Annotations
    proportionalLabel(0.04, 1.075, 'P', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(-0.5, 0.5, expression(paste(italic(f), " = 5.8")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)



    # Panel 17: h = 1/4,  C = 1/4, f = 5.8
    E_f5.8 <- interp(x=heatE_f5.8$sm, y=heatE_f5.8$sf, z=heatE_f5.8$lambda_sim, 
                xo=seq(min(heatE_f5.8$sm),max(heatE_f5.8$sm),by=resolution), 
                yo=seq(min(heatE_f5.8$sf),max(heatE_f5.8$sf),by=resolution), duplicate="mean")
    image(E_f5.8, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("0","0.05","0.10", "0.15"), cex.axis=1.2)
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    # Invasion Boundaries
    lines(aInvBound[1:68] ~ sms[1:68], lty=1, lwd=1.5, col=COLS$line, data=invE)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invE)
    # Labels/Annotations
    proportionalLabel(0.04, 1.075, 'Q', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(0.5, -0.3, expression(paste("Selection on male reproductive function (", italic(s[m]), ")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)

    # Panel 18: h = 1/2,  C = 1/2, f = 6.2
    F_f5.8 <- interp(x=heatF_f5.8$sm, y=heatF_f5.8$sf, z=heatF_f5.8$lambda_sim, 
                xo=seq(min(heatF_f5.8$sm),max(heatF_f5.8$sm),by=resolution), 
                yo=seq(min(heatF_f5.8$sf),max(heatF_f5.8$sf),by=resolution), duplicate="mean")
    image(F_f5.8, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
    axis(1, at=c(0,0.05,0.10, 0.1485), labels=c("0","0.05","0.10", "0.15"), cex.axis=1.2)
    axis(2, at=c(0,0.05,0.10, 0.1485), labels=c("","","", ""))
    # Invasion Boundaries
    lines(aInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invF)
    lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invF)
    # Labels/Annotations
    proportionalLabel(0.04, 1.075, 'R', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#    proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)

    # Color ScaleBar
    image.scale(A_f6.2, col=HEAT, breaks=breaks, horiz=FALSE, yaxt="n")
    axis(4, las=2, cex.axis=1.5)
    proportionalLabel(2.5, 0.5, expression(paste("Population intrinsic growth rate (", italic(lambda), ")")), cex=2.5, adj=c(0.5, 0.5), xpd=NA, srt=270)

}




resultsIllustrationFig  <-  function() {

    # Import Data
    invD        <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.25_hm0.25_C0_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    heatD_f5.8  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.25_hm0.25_C0_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    heatD_f5.9  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.25_hm0.25_C0_delta0_dj0_da0_dg0_f5.9.csv", head=TRUE)
    heatD_f6.0  <-  read.csv("./output/simData/lambdaHeatMapData_sMax0.15_len100_hf0.25_hm0.25_C0_delta0_dj0_da0_dg0_f6.csv", head=TRUE)

    extD_f5.8  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.25_C0_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extD_f5.9  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.25_C0_delta0_dj0_da0_dg0_f5.9.csv", head=TRUE)
    extD_f6.0  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.25_C0_delta0_dj0_da0_dg0_f6.csv", head=TRUE)

    # Clean up a couple missing values from data
    extD_f5.8$sfThreshold[33]  <-  mean(extD_f5.8$sfThreshold[32], extD_f5.8$sfThreshold[34])
    extD_f5.8$sfThreshold[37]  <-  mean(extD_f5.8$sfThreshold[36], extD_f5.8$sfThreshold[38])
    extD_f5.8$sf[11]           <-  mean(extD_f5.8$sf[11], extD_f5.8$sfThreshold[21])
    extD_f5.8$smThreshold[11]  <-  mean(extD_f5.8$smThreshold[12], extD_f5.8$sms[21])
    extD_f6.0$sfThreshold[83]  <-  (extD_f6.0$sfThreshold[82] + extD_f6.0$sfThreshold[84])/2

    # Color scheme
    COLS  <-  list("line"     =  transparentColor('#252525', opacity=1),
                   "extinct"  =  transparentColor('red', opacity=0.15))

    # heatmap colors modified from wes_anderson palette "Zissou1"
    Z1 = rev(c("dodgerblue4","dodgerblue4","#3B9AB2","#3B9AB2", "#EBCC2A", "#EBCC2A", "#E1AF00", "#E1AF00", "#F21A00", "#F21A00", "#F21A00"))
    blu     <-  colorRampPalette(Z1[length(Z1)])(10)
    notBlu  <-  colorRampPalette(Z1, interpolate = c("linear"))(60)
    HEAT    <-  c(notBlu, blu)
#    HEAT    <-  colorRampPalette(Z1, interpolate = c("linear"))(60)

    # resolution for heatmaps
    resolution = (heatD_f6.0$sf[2] - heatD_f6.0$sf[1])
    # breaks for heatmap colors
    breaks <- seq(min(heatD_f6.0$lambda_sim), 
                  max(heatD_f6.0$lambda_sim),length.out=(length(HEAT)+1))
#    image.scale(A_f6.2, col=HEAT, breaks=breaks, horiz=FALSE, yaxt="n")
#    axis(4, las=2)

    # Set plot layout
    layout.mat <- matrix(c(
                           1,1,1,1, 2,2,2,2, 3,
                           1,1,1,1, 2,2,2,2, 3,
                           1,1,1,1, 2,2,2,2, 3,
                           1,1,1,1, 2,2,2,2, 3), 
                            nrow=4, ncol=9, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

# Make the plot
    # Panel A: Population Genetic Predictions
    par(omi=c(0.5, 0.5, 0.5, 1.2), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s', xpd=TRUE)
    plot(NA, axes=FALSE, type='n', main='', xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
    # Invasion Boundaries
        lines(aInvBound[1:33] ~ sms[1:33], lty=1, lwd=2, col=COLS$line, data=invD)
        lines(AInvBound ~ sms, lty=1, lwd=2, col=COLS$line, data=invD)
    # axes        
        axis(1, at=c(0,0.05,0.10, 0.15), labels=c("0","0.05","0.10", "0.15"))
        axis(2, at=c(0,0.05,0.10, 0.15), labels=c("0","0.05","0.10", "0.15"))
        proportionalLabel(-0.28, 0.5, expression(paste("Selection on female")), cex=1.4, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.2, 0.5, expression(paste("reproductive function (", italic(s[f]), ")")), cex=1.4, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(1.125, -0.2, expression(paste("Selection on male reproductive function (", italic(s[m]), ")")), cex=1.4, adj=c(0.5, 0.5), xpd=NA, srt=0)
    # Labels/Annotations
#        proportionalLabel(0.03, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.1, 'Population Genetic Outcome', cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.125, 0.65, expression(paste(italic(A)," allele fixes")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=70)
        proportionalLabel(0.6, 0.6, expression(paste("Polymorphism")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.65, 0.1, expression(paste(italic(a)," allele fixes")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)

    # Panel B: Demographic Predictions
        D_f6.0 <- interp(x=heatD_f6.0$sm, y=heatD_f6.0$sf, z=heatD_f6.0$lambda_sim, 
                     xo=seq(min(heatD_f6.0$sm),max(heatD_f6.0$sm),by=resolution), 
                     yo=seq(min(heatD_f6.0$sf),max(heatD_f6.0$sf),by=resolution), duplicate="mean")
    # heatmap of lambda
        image(D_f6.0, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
        axis(1, at=c(0.0015,0.05,0.10, 0.15), labels=c("0","0.05","0.10", "0.15"))
        axis(2, at=c(0.0015,0.05,0.10, 0.15), labels=c("","","", ""))
    # Invasion Boundaries
#        lines(aInvBound[1:33] ~ sms[1:33], lty=4, lwd=1, col=COLS$line, data=invD)
#        lines(AInvBound ~ sms, lty=4, lwd=1, col=COLS$line, data=invD)
    # Extinction Boundary
        lines(sf ~ smThreshold, lty=1, lwd=1, col=COLS$line, data=extD_f6.0)   
        lines(sfThreshold ~ sms, lty=1, lwd=1, col=COLS$line, data=extD_f6.0)   
    # Labels/Annotations
        proportionalLabel(0.5, 1.1, 'Demographic Outcome', cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.03, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.9, 0.7, expression(paste("Extinction (", lambda, " < 1)")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=85)
#        proportionalLabel(0.6, 0.9, expression(paste("Extinction")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
#        proportionalArrows(px1=0.6, py1=0.85, px2=0.9, py2=0.6, length=0.1)
#        proportionalLabel(0.4, 0.35, expression(paste("Exp. Growth (", lambda," > 1)")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.55, 0.425, expression(paste("Exp. Growth")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.55, 0.35, expression(paste("(", lambda," > 1)")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)

    # Color ScaleBar
    image.scale(A_f6.2, col=HEAT, breaks=breaks, horiz=FALSE, yaxt="n")
    axis(4, las=2, cex.axis=1.1)
    proportionalLabel(4, 0.5, expression(paste("Population intrinsic growth rate (", italic(lambda), ")")), cex=1.25, adj=c(0.5, 0.5), xpd=NA, srt=270)

}











#' Fig showing Proportion of polymorphic parameter space
#' (w/out inbreeding depression)
polySpaceFigTitrate  <-  function(df = "dataPolySpaceFig_sMax0.15_res0.003_delta0_dj0_da0_dg0") {

    # Make filenames for import from df names
    fName  <-  paste('./output/simData/', df, '.csv', sep="")

    # Extract plotting parameter values from df names
    df1  <-  strsplit(df, '_')[[1]][c(2:7)]
    pars  <-  list(
                    "sMax"  =  as.numeric(strsplit(df1[1],'x')[[1]][2]),
                    "res"   =  as.numeric(strsplit(df1[2],'s')[[1]][2]),
                    "d"     =  as.numeric(strsplit(df1[3],'a')[[1]][2]),
                    "dj"    =  as.numeric(strsplit(df1[4],'j')[[1]][2]),
                    "da"    =  as.numeric(strsplit(df1[5],'a')[[1]][2]),
                    "dg"    =  as.numeric(strsplit(df1[6],'g')[[1]][2])
                    )

    # import data
    data  <-  read.csv(file=fName, header=TRUE)

    # clean data set & quantify parameter space
    dat   <-  quantPolySpace(data = data, pars = pars)

    # Color scheme
    COLS  <-  list(
                    "PG"    =  transparentColor('#252525', opacity=1),
                    "low"   =  transparentColor('dodgerblue4', opacity=0.6),
                    "med"   =  transparentColor('darkolivegreen4', opacity=0.6),
                    "hi"    =  transparentColor('tomato', opacity=0.6),
                    "low2"  =  transparentColor('dodgerblue4', opacity=1),
                    "med2"  =  transparentColor('darkolivegreen4', opacity=1),
                    "hi2"   =  transparentColor('tomato', opacity=1)
                   )

#  Create vector of selfing rates for pop gen function.
    CLine          <-  seq(0,0.9,length=100)
    addPGSpace     <-  c()
    domRevPGSpace  <-  c()
    for(i in 1:length(CLine)) {
        addPGSpace[i]     <-  popGen_PolySpace(hf=0.5, hm=0.5, C=CLine[i], sMax=pars$sMax)
        domRevPGSpace[i]  <-  popGen_PolySpace(hf=0.25, hm=0.25, C=CLine[i], sMax=pars$sMax)
    }

# Set plot layout
    layout.mat  <- matrix(c(1,2), nrow=2, ncol=1, byrow=TRUE)
    layout      <- layout(layout.mat,respect=TRUE)

    ##  Panel A: hf = hm = 1/2
        # subset data
        d  <-  dat[dat$h == 0.5,]
        dlow  <-  d[d$f == 5.8,]
        dmed  <-  d[d$f == 6.0,]
        dhi   <-  d[d$f == 6.5,]
        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(d$C)), ylim = c(0,0.105), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(addPGSpace ~ CLine, lwd=2, col=COLS$PG)
        points(PrViaPoly ~ C, pch=21, bg=COLS$low, col=COLS$low2, data=dlow)
        points(PrViaPoly ~ C, pch=21, bg=COLS$med, col=COLS$med2, data=dmed)
        points(PrViaPoly ~ C, pch=21, bg=COLS$hi, col=COLS$hi2, data=dhi)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.1, substitute(paste(italic(h), " = ", hh), list(hh = pars2$hf)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.5, 1.1, expression(paste(italic(h), " = ", 1/2)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(-0.3, -0.16, expression(paste("Proportion viable polymorphic parameter space")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=90)

        #Legend
        legend( x       =  usr[2],
                y       =  usr[4],
                legend  =  c(
                             expression(paste("Pop. Gen.")),
                             substitute(paste("High fertility (", italic(f), " = ", ff, ")"), list(ff = dhi$f[1])),
                             substitute(paste("Med. fertility (", italic(f), " = ", ff, ".0)"), list(ff = dmed$f[1])),
                             substitute(paste("Low  fertility (", italic(f), " = ", ff, ")"), list(ff = dlow$f[1]))),
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
        # subset data
        d  <-  dat[dat$h == 0.25,]
        dlow  <-  d[d$f == 5.8,]
        dmed  <-  d[d$f == 5.9,]
        dhi   <-  d[d$f == 6.5,]
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(d$C)), ylim = c(0,(max(d$PrViaPoly)*1.05)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(domRevPGSpace ~ CLine, lwd=2, col=COLS$PG)
        points(PrViaPoly ~ C, pch=21, bg=COLS$low, col=COLS$low2, data=dlow)
        points(PrViaPoly ~ C, pch=21, bg=COLS$med, col=COLS$med2, data=dmed)
        points(PrViaPoly ~ C, pch=21, bg=COLS$hi, col=COLS$hi2, data=dhi)
        # axes        
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(0.03, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.1, substitute(paste(italic(h), " = ", hh), list(hh = pars5$hf)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.5, 1.1, expression(paste(italic(h), " = ", 1/4)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.5, -0.3, expression(paste("Selfing rate (",italic(C), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
#        proportionalLabel(-0.3, 0.5, expression(paste("Prop. viable polymorphic space")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=90)

        #Legend
        legend( x       =  usr[2],
                y       =  usr[4],
                legend  =  c(
                             expression(paste("Pop. Gen.")),
                             substitute(paste("High fertility (", italic(f), " = ", ff, ")"), list(ff = dhi$f[1])),
                             substitute(paste("Med. fertility (", italic(f), " = ", ff, ")"), list(ff = dmed$f[1])),
                             substitute(paste("Low  fertility (", italic(f), " = ", ff, ")"), list(ff = dlow$f[1]))),
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





#' Fig showing Effect of Inbreeding Depression on 
#' Proportion of polymorphic parameter space

deltaSelfingLoadPolySpaceFigTitrate  <-  function(df1 = "dataDeltaPolySpaceFig_sMax0.15_res0.003_dStar0.8_f6.5",
                                                  df2 = "dataDeltaPolySpaceFig_sMax0.15_res0.003_dStar0.8_f7.5",
                                                  df3 = "dataDeltaPolySpaceFig_sMax0.15_res0.003_dStar0.8_f8.5") {

    # Make filenames for import from df names
    fName1  <-  paste('./output/simData/', df1, '.csv', sep="")
    fName2  <-  paste('./output/simData/', df2, '.csv', sep="")
    fName3  <-  paste('./output/simData/', df3, '.csv', sep="")

    # import data
    data1  <-  read.csv(file=fName1, header=TRUE)
    data2  <-  read.csv(file=fName2, header=TRUE)
    data3  <-  read.csv(file=fName3, header=TRUE)

    # Extract plotting parameter values from df names
    d1   <-  strsplit(df1, '_')[[1]][c(2:5)]
    pars1  <-  list(
                    "sMax"  =  as.numeric(strsplit(d1[1],'x')[[1]][2]),
                    "res"   =  as.numeric(strsplit(d1[2],'s')[[1]][2]),
                    "dStar" =  as.numeric(strsplit(d1[3],'r')[[1]][2]),
                    "f"     =  as.numeric(strsplit(d1[4],'f')[[1]][2])
                    )
    d2   <-  strsplit(df2, '_')[[1]][c(2:5)]
    pars2  <-  list(
                    "sMax"  =  as.numeric(strsplit(d2[1],'x')[[1]][2]),
                    "res"   =  as.numeric(strsplit(d2[2],'s')[[1]][2]),
                    "dStar" =  as.numeric(strsplit(d2[3],'r')[[1]][2]),
                    "f"     =  as.numeric(strsplit(d2[4],'f')[[1]][2])
                    )
    d3   <-  strsplit(df3, '_')[[1]][c(2:5)]
    pars3  <-  list(
                    "sMax"  =  as.numeric(strsplit(d2[1],'x')[[1]][2]),
                    "res"   =  as.numeric(strsplit(d2[2],'s')[[1]][2]),
                    "dStar" =  as.numeric(strsplit(d2[3],'r')[[1]][2]),
                    "f"     =  as.numeric(strsplit(d2[4],'f')[[1]][2])
                    )
    hLev  <-  unique(data1$h)
    dLev  <-  unique(data1$Delta)
    CLev  <-  unique(data1$C)
    nHs   <-  length(hLev)
    nDs   <-  length(dLev)
    nCs   <-  length(CLev)

    # clean data set & quantify parameter space
    dat1   <-  quantDeltaPolySpace(data = data1, pars = pars1)
    dat1[dat1 < 0]  <-  0
    dat2   <-  quantDeltaPolySpace(data = data2, pars = pars2)
    dat2[dat2 < 0]  <-  0
    dat3   <-  quantDeltaPolySpace(data = data3, pars = pars3)
    dat3[dat3 < 0]  <-  0

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
    ##          early-acting delta
        PGSpace     <-  c()
        for(i in 1:length(CLine)) {
            PGSpace[i]     <-  popGen_PolySpace_Delta_Add(C=CLine[i], delta=dLine[i], sMax=pars1$sMax)
        }
        d1  <-  dat1[dat1$h == 0.5,]
        d2  <-  dat2[dat2$h == 0.5,]
        d2$PrViaPoly[d2$Delta == "d_g"][c(18)]  <-  as.character((as.numeric(d2$PrViaPoly[d2$Delta == "d_g"][c(17)])  + as.numeric(d2$PrViaPoly[d2$Delta == "d_g"][c(19)]) )/2)
        d3  <-  dat3[dat3$h == 0.5,]
        d3$PrViaPoly[d3$Delta == "d_g"][c(18)]  <-  as.character((as.numeric(d3$PrViaPoly[d3$Delta == "d_g"][c(17)])  + as.numeric(d3$PrViaPoly[d3$Delta == "d_g"][c(19)]) )/2)

        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,4,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.925), ylim = c(0,0.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        # f = 6.5
        points(PrViaPoly[Delta == "d"][c(1:9,23)] ~ C[Delta == "d"][c(1:9,23)], pch=21, bg=COLS$dSim, col=COLS$dSim2, data=d1)
        points(PrViaPoly[Delta == "d_j"][c(1:7,23)] ~ C[Delta == "d_j"][c(1:7,23)], pch=21, bg=COLS$d_j, col=COLS$d_j2, data=d1)
        points(PrViaPoly[Delta == "d_a"][c(1:12,23)] ~ C[Delta == "d_a"][c(1:12,23)], pch=21, bg=COLS$d_a, col=COLS$d_a2, data=d1)
        points(PrViaPoly[Delta == "d_g"][c(1:10,23)] ~ C[Delta == "d_g"][c(1:10,23)], pch=21, bg=COLS$d_g, col=COLS$d_g2, data=d1)
        # f = 7.5
        points(PrViaPoly[Delta == "d"][c(1:21,37)] ~ C[Delta == "d"][c(1:21,37)], pch=3, col=COLS$dSim2, data=d2)
        points(PrViaPoly[Delta == "d_j"][c(1:15,37)] ~ C[Delta == "d_j"][c(1:15,37)], pch=3, col=COLS$d_j2, data=d2)
        points(PrViaPoly[Delta == "d_a"][c(1:24,37)] ~ C[Delta == "d_a"][c(1:24,37)], pch=3, col=COLS$d_a2, data=d2)
        points(PrViaPoly[Delta == "d_g"][c(1:22,37)] ~ C[Delta == "d_g"][c(1:22,37)], pch=3, col=COLS$d_g2, data=d2)
        # f = 8.5
        points(PrViaPoly[Delta == "d"][c(1:32,37)] ~ C[Delta == "d"][c(1:32,37)], pch=4, col=COLS$dSim2, data=d3)
        points(PrViaPoly[Delta == "d_j"][c(1:21,37)] ~ C[Delta == "d_j"][c(1:21,37)], pch=4, col=COLS$d_j2, data=d3)
        points(PrViaPoly[Delta == "d_a"][c(1:35,37)] ~ C[Delta == "d_a"][c(1:35,37)], pch=4, col=COLS$d_a2, data=d3)
        points(PrViaPoly[Delta == "d_g"][c(1:33,37)] ~ C[Delta == "d_g"][c(1:33,37)], pch=4, col=COLS$d_g2, data=d3)
        # PG expectation
        lines(PGSpace ~ CLine, lwd=2, col=COLS$PG)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Labels/annotations
        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(italic(h), " = ", 1/2)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(-0.35, -0.16, expression(paste("Proportion viable polymorphic parameter space")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        legend( x       =  usr[2],
                y       =  usr[4],
                legend  =  c(
                             expression(paste(italic(delta), " (pop. gen.)")),
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
        legend( x       =  usr[2]*0.93,
                y       =  usr[4]*0.62,
                legend  =  c(
                             expression(paste(italic(f), " = ", 8.5)),
                             expression(paste(italic(f), " = ", 7.5)),
                             expression(paste(italic(f), " = ", 6.5))),
                 col     =  c(COLS$PG),
                 pch     =  c(4,3,21),
                 pt.bg   =  c(NA),
                 cex     =  0.75,
                 pt.cex  =  0.75,
                 xjust   =  1,
                 yjust   =  1,
                 bty     =  'n',
                 border  =  NA)


    ## Panel B: Dominance Reversal SA (hf = hm = 1/4)
    ##          early-acting delta
        CLine2  <-  CLine[-1]
        dLine2  <-  dLine[-1]
        PGSpace     <-  c()
        for(i in 1:length(CLine2)) {
            PGSpace[i]     <-  popGen_PolySpace_Delta_DomRev(C=CLine2[i], delta=dLine2[i], sMax=pars2$sMax)
        }
        d1  <-  dat1[dat1$h == 0.25,]
        d2  <-  dat2[dat2$h == 0.25,]
        d3  <-  dat3[dat3$h == 0.25,]
        d1$PrViaPoly[30]  <-  "0"
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='', xlim = c(0,0.925), ylim = c(0,0.8), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        # f = 6.5
        points(PrViaPoly[Delta == "d"][c(1:10,23)] ~ C[Delta == "d"][c(1:10,23)], pch=21, bg=COLS$dSim, col=COLS$dSim2, data=d1)
        points(PrViaPoly[Delta == "d_j"][c(1:7,23)] ~ C[Delta == "d_j"][c(1:7,23)], pch=21, bg=COLS$d_j, col=COLS$d_j2, data=d1)
        points(PrViaPoly[Delta == "d_a"][c(1:12,23)] ~ C[Delta == "d_a"][c(1:12,23)], pch=21, bg=COLS$d_a, col=COLS$d_a2, data=d1)
        points(PrViaPoly[Delta == "d_g"][c(1:10,23)] ~ C[Delta == "d_g"][c(1:10,23)], pch=21, bg=COLS$d_g, col=COLS$d_g2, data=d1)
        # f = 7.5
        points(PrViaPoly[Delta == "d"][c(1:22,37)] ~ C[Delta == "d"][c(1:22,37)], pch=3, bg=COLS$dSim, col=COLS$dSim2, data=d2)
        points(PrViaPoly[Delta == "d_j"][c(1:15,37)] ~ C[Delta == "d_j"][c(1:15,37)], pch=3, bg=COLS$d_j, col=COLS$d_j2, data=d2)
        points(PrViaPoly[Delta == "d_a"][c(1:25,37)] ~ C[Delta == "d_a"][c(1:25,37)], pch=3, bg=COLS$d_a, col=COLS$d_a2, data=d2)
        points(PrViaPoly[Delta == "d_g"][c(1:23,37)] ~ C[Delta == "d_g"][c(1:23,37)], pch=3, bg=COLS$d_g, col=COLS$d_g2, data=d2)        # axes        
        # f = 8.5
        points(PrViaPoly[Delta == "d"][c(1:32,37)] ~ C[Delta == "d"][c(1:32,37)], pch=4, col=COLS$dSim2, data=d3)
        points(PrViaPoly[Delta == "d_j"][c(1:21,37)] ~ C[Delta == "d_j"][c(1:21,37)], pch=4, col=COLS$d_j2, data=d3)
        points(PrViaPoly[Delta == "d_a"][c(1:35,37)] ~ C[Delta == "d_a"][c(1:35,37)], pch=4, col=COLS$d_a2, data=d3)
        points(PrViaPoly[Delta == "d_g"][c(1:34,37)] ~ C[Delta == "d_g"][c(1:34,37)], pch=4, col=COLS$d_g2, data=d3)
        # Pop. Gen Prediction
        lines(PGSpace ~ CLine2, lwd=2, col=COLS$PG)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Labels/annotations
        proportionalLabel(0.03, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(italic(h), " = ", 1/4)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
#        proportionalLabel(-0.35, 0.5, expression(paste("Prop. viable polymorphic space")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.35, expression(paste("Selfing Rate (", italic(C), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)

}





# Simple inv6 Fig
MimulusInv6Fig  <-  function() {

    # import data
    path  <-  './output/simData/'
    inv1  <-  read.csv(paste(path, "invBoundMimulus_sMax0.99_res0.01_hf0.5_hm0.5_C0_deltaF_useCompadreFALSE_EM", '.csv', sep=""), header=TRUE)
    inv2  <-  read.csv(paste(path, "invBoundMimulus_sMax0.99_res0.01_hf0.5_hm0.5_C0.29_deltaT_useCompadreFALSE_EM", '.csv', sep=""), header=TRUE)
    inv3  <-  read.csv(paste(path, "invBoundMimulus_sMax0.99_res0.01_hf0.5_hm0.5_C0_deltaF_useCompadreFALSE_LEP", '.csv', sep=""), header=TRUE)
    inv4  <-  read.csv(paste(path, "invBoundMimulus_sMax0.99_res0.01_hf0.5_hm0.5_C0.29_deltaT_useCompadreFALSE_LEP", '.csv', sep=""), header=TRUE)
    ext3  <-  read.csv(paste(path, "extThresholdMimulus_SfxSm_sMax0.99_res0.01_hf0.5_hm0.5_C0_useCompadreFALSE_IDFALSE_LEP", '.csv', sep=""), header=TRUE)
    ext4  <-  read.csv(paste(path, "extThresholdMimulus_SfxSm_sMax0.99_res0.01_hf0.5_hm0.5_C0.29_useCompadreFALSE_IDTRUE_LEP", '.csv', sep=""), header=TRUE)

#    ext4$sf[46]  <-  mean(ext4$sf[45], ext4$sf[47])
    ext4$smThreshold[46]  <-  mean(c(ext4$smThreshold[45], ext4$smThreshold[47]))
    # Selection coefficients for Inv6 (including average flower # effect)
    sfInv6  <-  0.31
    smInv6  <-  0.30

# Color scheme
    COLS  <-  list(
                    "line"   =  transparentColor('#252525', opacity=1),
                    "fill"   =  transparentColor('#252525', opacity=0.5),
                    "inv6"   =  transparentColor('tomato', opacity=1),
                    "inv6bg" =  transparentColor('tomato', opacity=0.6),
                    "extinct"  =  transparentColor('red', opacity=0.15)
                    )
# Set plot layout
    layout.mat  <- matrix(c(1:2), nrow=1, ncol=2, byrow=TRUE)
    layout      <- layout(layout.mat,respect=TRUE)

    ## Panel A: Outcrossing

        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,4,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
#        polygon(c(rev(sms),sms), c(rev(inv_A), inv_a), col=transparentColor('grey80', 0.6), border=NA)
        lines(AInvBound ~ sms, lwd=2, col=COLS$line, data=inv1)
        lines(aInvBound[aInvBound < max(sms)] ~ sms[aInvBound < max(sms)], lwd=2, col=COLS$line, data=inv1)
        lines(AInvBound ~ sms, lty=2, lwd=2, col=COLS$line, data=inv2)
        lines(c(inv2$aInvBound[inv2$aInvBound < max(inv2$sms)],max(inv2$sms)) ~ c(inv2$sms[inv2$aInvBound < max(inv2$sms)],0.50), lty=2, lwd=2, col=COLS$line, data=inv2)
        # no lines to plot for extinction threshold (all polymorphic space is viable)
        points(sfInv6 ~ smInv6, pch=21, col=COLS$line, bg=COLS$fill, cex=1.25)
        # axes        
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel((smInv6-0.1), (sfInv6+0.05), expression('inv6'), cex=1, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, 'Local Demographic Rates', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.25, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.25, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
      #Legend
        legend( x       =  usr[2],
                y       =  usr[4]*0.98,
#        legend( x       =  usr[2]*0.45,
#                y       =  usr[4],
                legend  =  c(
                             expression(paste(italic(C)," = 0")),
                             expression(paste(italic(C)," = 0.29, I.D."))),
                 lty     =  c(1,2),
                 lwd     =  c(2,2),
                 col     =  c(COLS$line),
                 seg.len =  c(2),
                 cex     =  0.75,
                 xjust   =  1,
                 yjust   =  1,
                 bty     =  'n',
                 border  =  NA)

        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
#        polygon(c(rev(sms),sms), c(rev(inv_A), inv_a), col=transparentColor('grey80', 0.6), border=NA)
        lines(AInvBound ~ sms, lwd=2, col=COLS$line, data=inv3)
        lines(aInvBound[aInvBound < max(sms)] ~ sms[aInvBound < max(sms)], lwd=2, col=COLS$line, data=inv3)
        lines(AInvBound ~ sms, lty=2, lwd=2, col=COLS$line, data=inv4)
        lines(c(aInvBound[aInvBound < max(sms)],0.99) ~ c(sms[aInvBound < max(sms)],0.562), lty=2, lwd=2, col=COLS$line, data=inv4)
        # PLOT LINES FOR EXTINCTION THRESHOLDS
        t    <-  cbind(ext3$sf[!is.na(ext3$smThreshold)], ext3$smThreshold[!is.na(ext3$smThreshold)])
        t    <-  rbind(t, c(0.99, 0.99))
        b    <-  cbind(ext3$sms[!is.na(ext3$sfThreshold)], ext3$sfThreshold[!is.na(ext3$sfThreshold)])
        polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lwd=1.5, col=COLS$line, data=ext3)
        lines(sfThreshold ~ sms, lwd=1.5, col=COLS$line, data=ext3)        

        t    <-  cbind(ext4$sf[!is.na(ext4$smThreshold)], ext4$smThreshold[!is.na(ext4$smThreshold)])
        t    <-  rbind(t, c(0.99, 0.99))
        b    <-  cbind(ext4$sms[!is.na(ext4$sfThreshold)], ext4$sfThreshold[!is.na(ext4$sfThreshold)])
        polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lwd=1.5, lty=2, col=COLS$line, data=ext4)
        lines(sfThreshold ~ sms, lwd=1.5, lty=2, col=COLS$line, data=ext4)
        # Plot inv6
        points(sfInv6 ~ smInv6, pch=21, col=COLS$line, bg=COLS$fill, cex=1.25)
        # axes        
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel((smInv6-0.1), (sfInv6+0.05), expression('inv6'), cex=1, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, 'Non-local Demographic Rates', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.15, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.25, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)

}


# Mimulus fig w/ heatmap for lambda + inv6 
MimulusInv6LambdaFig  <-  function() {

    # import data
    path  <-  './output/simData/'
    inv1  <-  read.csv(paste(path, "invBoundMimulus_sMax0.99_res0.01_hf0.5_hm0.5_C0_deltaF_useCompadreFALSE_EM", '.csv', sep=""), header=TRUE)
    inv2  <-  read.csv(paste(path, "invBoundMimulus_sMax0.99_res0.01_hf0.5_hm0.5_C0_deltaF_useCompadreFALSE_LEP", '.csv', sep=""), header=TRUE)
    inv3  <-  read.csv(paste(path, "invBoundMimulus_sMax0.99_res0.01_hf0.5_hm0.5_C0.29_deltaT_useCompadreFALSE_EM", '.csv', sep=""), header=TRUE)
    inv4  <-  read.csv(paste(path, "invBoundMimulus_sMax0.99_res0.01_hf0.5_hm0.5_C0.29_deltaT_useCompadreFALSE_LEP", '.csv', sep=""), header=TRUE)
    heat_1  <-  read.csv("./output/simData/lambdaHeatMapMimulusData_sMax0.99_len100_alpha0_hf0.5_hm0.5_C0_useCompadreFALSE_IDFALSE_EM.csv", head=TRUE)
    heat_2  <-  read.csv("./output/simData/lambdaHeatMapMimulusData_sMax0.99_len100_alpha0_hf0.5_hm0.5_C0_useCompadreFALSE_IDFALSE_LEP.csv", head=TRUE)
    heat_3  <-  read.csv("./output/simData/lambdaHeatMapMimulusData_sMax0.99_len100_alpha0_hf0.5_hm0.5_C0.29_useCompadreFALSE_IDTRUE_EM.csv", head=TRUE)
    heat_4  <-  read.csv("./output/simData/lambdaHeatMapMimulusData_sMax0.99_len100_alpha0_hf0.5_hm0.5_C0.29_useCompadreFALSE_IDTRUE_LEP.csv", head=TRUE)

    # Clean up abberrant lambda_sim values
        oddOnes  <-  which(heat_1$lambda_sim < 1e-14)
        heat_1$lambda_sim[oddOnes]  <-  (heat_1$lambda_sim[oddOnes+1] + heat_1$lambda_sim[oddOnes-1])/2    
        # two more odd results: 
        oddOnes  <-  which(heat_1$lambda_sim < 1)
        heat_1$lambda_sim[oddOnes]  <-  (heat_1$lambda_sim[oddOnes+2] + heat_1$lambda_sim[oddOnes-1])/2    

        oddOnes  <-  which(heat_2$lambda_sim < 1e-14)
        heat_2$lambda_sim[oddOnes]  <-  (heat_2$lambda_sim[oddOnes+1] + heat_2$lambda_sim[oddOnes-1])/2    

        oddOnes  <-  which(heat_3$lambda_sim < 1e-14)
        heat_3$lambda_sim[oddOnes]  <-  (heat_3$lambda_sim[oddOnes+1] + heat_3$lambda_sim[oddOnes-1])/2    

        oddOnes  <-  which(heat_4$lambda_sim < 1e-14)
        heat_4$lambda_sim[oddOnes]  <-  (heat_4$lambda_sim[oddOnes+1] + heat_4$lambda_sim[oddOnes-1])/2    


    # heatmap colors modified from wes_anderson palette "Zissou1"
    Z1 = rev(c("dodgerblue4","dodgerblue4","#3B9AB2","#3B9AB2", "#EBCC2A", "#EBCC2A", "#E1AF00", "#E1AF00", "#F21A00", "#F21A00", "#F21A00"))
    blu     <-  colorRampPalette(Z1[length(Z1)])(10)
    notBlu  <-  colorRampPalette(Z1, interpolate = c("linear"))(60)
    HEAT    <-  c(notBlu, blu)
#    HEAT    <-  colorRampPalette(Z1, interpolate = c("linear"))(60)

    # resolution for heatmaps
    resolution = (heat_1$sf[2] - heat_1$sf[1])
    # breaks for heatmap colors
#    breaks <- seq(min(heat_1$lambda_sim,
#                      heat_2$lambda_sim,
#                      heat_3$lambda_sim,
#                      heat_4$lambda_sim), 
#                  max(heat_1$lambda_sim,
#                      heat_2$lambda_sim,
#                      heat_3$lambda_sim,
#                      heat_4$lambda_sim),length.out=(length(HEAT)+1))
    breaks <- c(seq(min(heat_1$lambda_sim,
                      heat_2$lambda_sim,
                      heat_3$lambda_sim,
                      heat_4$lambda_sim), 
                  1.1,length.out=(length(HEAT))),1.1)
    #  relabel large lambdas so color bar is pretty
        heat_1$lambda_sim[heat_1$lambda_sim > 1.1]  <- 1.1
        heat_3$lambda_sim[heat_3$lambda_sim > 1.1]  <- 1.1

    # Set plot layout
    layout.mat <- matrix(c(
                           1,1,1,1, 2,2,2,2, 3,
                           1,1,1,1, 2,2,2,2, 3,
                           1,1,1,1, 2,2,2,2, 3,
                           1,1,1,1, 2,2,2,2, 3,
                           4,4,4,4, 5,5,5,5, 6,
                           4,4,4,4, 5,5,5,5, 6,
                           4,4,4,4, 5,5,5,5, 6,
                           4,4,4,4, 5,5,5,5, 6), 
                            nrow=8, ncol=9, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Selection coefficients for Inv6 (including average flower # effect)
    sfInv6  <-  0.31
    smInv6  <-  0.30

# Color scheme
    COLS  <-  list(
                    "line"   =  transparentColor('#252525', opacity=1),
                    "fill"   =  transparentColor('#252525', opacity=0.5),
                    "inv6"   =  transparentColor('tomato', opacity=1),
                    "inv6bg" =  transparentColor('tomato', opacity=0.6),
                    "extinct"  =  transparentColor('red', opacity=0.15)
                    )

## Panel A: Eagle Meadows Population, Obligate Outcrossing
    # Make the plot
    par(omi=c(0.5, 1, 0.5, 1.2), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
    EM_C0 <- interp(x=heat_1$sm, y=heat_1$sf, z=heat_1$lambda_sim, 
                     xo=seq(min(heat_1$sm),max(heat_1$sm),by=resolution), 
                     yo=seq(min(heat_1$sf),max(heat_1$sf),by=resolution), duplicate="mean")
    # heatmap of lambda
        image(EM_C0, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
        axis(1, at=c(0.01, 0.25, 0.5, 0.75, 0.99), labels=c("","","","",""))
        axis(2, at=c(0.01, 0.25, 0.5, 0.75, 0.99), labels=c("0","0.25","0.5", "0.75", "1.0"))
    # Invasion Boundaries
        lines(AInvBound ~ sms, lwd=2, col=COLS$line, data=inv1)
        lines(aInvBound[aInvBound < max(sms)] ~ sms[aInvBound < max(sms)], lwd=2, col=COLS$line, data=inv1)
    # Plot inv6
        points(sfInv6 ~ smInv6, pch=21, col=COLS$line, bg='white', cex=1.25)
        proportionalLabel((smInv6-0.1), (sfInv6+0.05), expression('inv6'), cex=1, adj=c(0.5, 0.5), xpd=NA, col='white')
    # Labels/Annotations
        proportionalLabel(0.03, 1.04, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, 'Local Demographic Rates', cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.4, 0.5, expression(italic(C) == 0), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)

## Panel B: Low-Elevation Perennials, Obligate Outcrossing
        # Make the plot
        LEP_C0 <- interp(x=heat_2$sm, y=heat_2$sf, z=heat_2$lambda_sim, 
                     xo=seq(min(heat_2$sm),max(heat_2$sm),by=resolution), 
                     yo=seq(min(heat_2$sf),max(heat_2$sf),by=resolution), duplicate="mean")
    # heatmap of lambda
        image(LEP_C0, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
        axis(1, at=c(0.01, 0.25, 0.5, 0.75, 0.99), labels=c("","","","",""))
        axis(2, at=c(0.01, 0.25, 0.5, 0.75, 0.99), labels=c("","","","",""))
    # Invasion Boundaries
        lines(AInvBound ~ sms, lwd=2, col=COLS$line, data=inv2)
        lines(aInvBound[aInvBound < max(sms)] ~ sms[aInvBound < max(sms)], lwd=2, col=COLS$line, data=inv2)
    # Plot inv6
        points(sfInv6 ~ smInv6, pch=21, col=COLS$line, bg='white', cex=1.25)
        proportionalLabel((smInv6-0.1), (sfInv6+0.05), expression('inv6'), cex=1, adj=c(0.5, 0.5), xpd=NA, col='white')
    # Labels/Annotations
        proportionalLabel(0.03, 1.04, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, 'Non-local Demographic Rates', cex=1.5, adj=c(0.5, 0.5), xpd=NA)

    # Color ScaleBar
        image.scale(LEP_C0, col=HEAT, breaks=breaks, horiz=FALSE, yaxt="n", xpd=TRUE)
        axis(4, at=c(0.98,1.00,1.02,1.04,1.06,1.08,1.10), labels=c("0.98","1.00","1.02","1.04","1.06","1.08","> 1.10"), las=2, cex.axis=1.1)
        proportionalLabel(5, 0.5, expression(paste("Population intrinsic growth rate (", italic(lambda), ")")), cex=1.25, adj=c(0.5, 0.5), xpd=NA, srt=270)



## Panel C: Eagle Meadows Population, Selfing w/ I.D>
    EM_C0.29 <- interp(x=heat_3$sm, y=heat_3$sf, z=heat_3$lambda_sim, 
                     xo=seq(min(heat_3$sm),max(heat_3$sm),by=resolution), 
                     yo=seq(min(heat_3$sf),max(heat_3$sf),by=resolution), duplicate="mean")
    # heatmap of lambda
        image(EM_C0.29, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
        axis(1, at=c(0.01, 0.25, 0.5, 0.75, 0.99), labels=c("0","0.25","0.5", "0.75", "1.0"))
        axis(2, at=c(0.01, 0.25, 0.5, 0.75, 0.99), labels=c("0","0.25","0.5", "0.75", "1.0"))
    # Invasion Boundaries
        lines(AInvBound ~ sms, lwd=2, col=COLS$line, data=inv3)
        lines(aInvBound[aInvBound < max(sms)] ~ sms[aInvBound < max(sms)], lwd=2, col=COLS$line, data=inv3)
    # Plot inv6
        points(sfInv6 ~ smInv6, pch=21, col=COLS$line, bg='white', cex=1.25)
        proportionalLabel((smInv6-0.1), (sfInv6+0.05), expression('inv6'), cex=1, adj=c(0.5, 0.5), xpd=NA, col='white')
    # Labels/Annotations
        proportionalLabel(0.03, 1.04, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(C), " = 0.29, w/ I.D.")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)

        proportionalLabel(-0.25, 1.125, expression(paste("Selection on female reproductive function (", italic(s[f]), ")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(1.125, -0.25, expression(paste("Selection on male reproductive function (", italic(s[m]), ")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)

## Panel D: Low-Elevation Perennials, Obligate Outcrossing
        # Make the plot
        LEP_C0.29 <- interp(x=heat_4$sm, y=heat_4$sf, z=heat_4$lambda_sim, 
                     xo=seq(min(heat_4$sm),max(heat_4$sm),by=resolution), 
                     yo=seq(min(heat_4$sf),max(heat_4$sf),by=resolution), duplicate="mean")
    # heatmap of lambda
        image(LEP_C0.29, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
        axis(1, at=c(0.01, 0.25, 0.5, 0.75, 0.99), labels=c("0","0.25","0.5", "0.75", "1.0"))
        axis(2, at=c(0.01, 0.25, 0.5, 0.75, 0.99), labels=c("","","","",""))
    # Invasion Boundaries
        lines(AInvBound ~ sms, lwd=2, col=COLS$line, data=inv4)
        lines(aInvBound[aInvBound < max(sms)] ~ sms[aInvBound < max(sms)], lwd=2, col=COLS$line, data=inv4)
    # Plot inv6
        points(sfInv6 ~ smInv6, pch=21, col=COLS$line, bg='white', cex=1.25)
        proportionalLabel((smInv6-0.1), (sfInv6+0.05), expression('inv6'), cex=1, adj=c(0.5, 0.5), xpd=NA, col='black')
    # Labels/Annotations
        proportionalLabel(0.03, 1.04, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # Color ScaleBar
        image.scale(LEP_C0, col=HEAT, breaks=breaks, horiz=FALSE, yaxt="n", xpd=TRUE)
        axis(4, at=c(0.98,1.00,1.02,1.04,1.06,1.08,1.10), labels=c("0.98","1.00","1.02","1.04","1.06","1.08","> 1.10"), las=2, cex.axis=1.1)
        proportionalLabel(5, 0.5, expression(paste("Population intrinsic growth rate (", italic(lambda), ")")), cex=1.25, adj=c(0.5, 0.5), xpd=NA, srt=270)


}







# Mimulus fig w/ heatmap for lambda + inv6 
MimulusInv6NstarFig  <-  function() {

    # import data
    path    <-  './output/simData/'
    inv1    <-  read.csv(paste(path, "invBoundMimulus_sMax0.99_res0.01_hf0.5_hm0.5_C0_deltaF_useCompadreFALSE_EM", '.csv', sep=""), header=TRUE)
    inv2    <-  read.csv(paste(path, "invBoundMimulus_sMax0.99_res0.01_hf0.5_hm0.5_C0_deltaF_useCompadreFALSE_LEP", '.csv', sep=""), header=TRUE)
    inv3    <-  read.csv(paste(path, "invBoundMimulus_sMax0.99_res0.01_hf0.5_hm0.5_C0.29_deltaT_useCompadreFALSE_EM", '.csv', sep=""), header=TRUE)
    inv4    <-  read.csv(paste(path, "invBoundMimulus_sMax0.99_res0.01_hf0.5_hm0.5_C0.29_deltaT_useCompadreFALSE_LEP", '.csv', sep=""), header=TRUE)
    heat_1  <-  read.csv("./output/simData/lambdaHeatMapMimulusData_sMax0.99_len100_alpha1e-04_hf0.5_hm0.5_C0_useCompadreFALSE_IDFALSE_EM.csv", head=TRUE)
    heat_2  <-  read.csv("./output/simData/lambdaHeatMapMimulusData_sMax0.99_len100_alpha1e-04_hf0.5_hm0.5_C0_useCompadreFALSE_IDFALSE_LEP.csv", head=TRUE)
    heat_3  <-  read.csv("./output/simData/lambdaHeatMapMimulusData_sMax0.99_len100_alpha1e-04_hf0.5_hm0.5_C0.29_useCompadreFALSE_IDTRUE_EM.csv", head=TRUE)
    heat_4  <-  read.csv("./output/simData/lambdaHeatMapMimulusData_sMax0.99_len100_alpha1e-04_hf0.5_hm0.5_C0.29_useCompadreFALSE_IDTRUE_LEP.csv", head=TRUE)

    # Clean up extremely small values so plots scale reasonably well.
        smallOnes  <-  which(heat_2$Nstar < 1e-5)
        heat_2$Nstar[smallOnes]  <-  1e-5    

        smallOnes  <-  which(heat_4$Nstar < 1e-5)
        heat_4$Nstar[smallOnes]  <-  1e-5

    # Calculate log equilibrium frequency
    heat_1$lNstar  <-  log(heat_1$Nstar)
    heat_2$lNstar  <-  log(heat_2$Nstar)
    heat_3$lNstar  <-  log(heat_3$Nstar)
    heat_4$lNstar  <-  log(heat_4$Nstar)

    # heatmap colors modified from wes_anderson palette "Zissou1"
    Z1 = rev(c("dodgerblue4","#3B9AB2", "#EBCC2A", "#E1AF00", "#F21A00"))
#    blu     <-  colorRampPalette(Z1[length(Z1)])(0)
#    notBlu  <-  colorRampPalette(Z1, interpolate = c("linear"))(60)
#    HEAT    <-  c(notBlu, blu)
    HEAT  <-  colorRampPalette(Z1, interpolate = c("linear"))(60)

    # resolution for heatmaps
    resolution = (heat_1$sf[2] - heat_1$sf[1])

    # breaks for heatmap colors
    breaks <- seq(min(heat_2$lNstar,
                      heat_4$lNstar), 
                  max(heat_1$lNstar,
                      heat_3$lNstar),length.out=(length(HEAT)+1))
    # fix aberrant point corresponding to min(heat_4$lNstar)
    heat_4$lNstar[heat_4$sf == 0.2970 & heat_4$sm == 0.4455]  <-  -11.5129 # original value: -11.51293
    
    # Set plot layout
    layout.mat <- matrix(c(
                           1,1,1,1, 2,2,2,2, 3,
                           1,1,1,1, 2,2,2,2, 3,
                           1,1,1,1, 2,2,2,2, 3,
                           1,1,1,1, 2,2,2,2, 3,
                           4,4,4,4, 5,5,5,5, 6,
                           4,4,4,4, 5,5,5,5, 6,
                           4,4,4,4, 5,5,5,5, 6,
                           4,4,4,4, 5,5,5,5, 6), 
                            nrow=8, ncol=9, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Selection coefficients for Inv6 (including average flower # effect)
    sfInv6  <-  0.31
    smInv6  <-  0.30

# Color scheme
    COLS  <-  list(
                    "line"   =  transparentColor('#252525', opacity=1),
                    "fill"   =  transparentColor('#252525', opacity=0.5),
                    "inv6"   =  transparentColor('tomato', opacity=1),
                    "inv6bg" =  transparentColor('tomato', opacity=0.6),
                    "extinct"  =  transparentColor('red', opacity=0.15)
                    )

## Panel A: Eagle Meadows Population, Obligate Outcrossing
    # Make the plot
    par(omi=c(0.5, 1, 0.5, 1.2), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
    EM_C0 <- interp(x=heat_1$sm, y=heat_1$sf, z=heat_1$lNstar, 
                     xo=seq(min(heat_1$sm),max(heat_1$sm),by=resolution), 
                     yo=seq(min(heat_1$sf),max(heat_1$sf),by=resolution), duplicate="mean")
    # heatmap of lambda
        image(EM_C0, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
        axis(1, at=c(0.01, 0.25, 0.5, 0.75, 0.99), labels=c("","","","",""))
        axis(2, at=c(0.01, 0.25, 0.5, 0.75, 0.99), labels=c("0","0.25","0.5", "0.75", "1.0"))
    # Invasion Boundaries
        lines(AInvBound ~ sms, lwd=2, col=COLS$line, data=inv1)
        lines(aInvBound[aInvBound < max(sms)] ~ sms[aInvBound < max(sms)], lwd=2, col=COLS$line, data=inv1)
    # Plot inv6
        points(sfInv6 ~ smInv6, pch=21, col=COLS$line, bg='white', cex=1.25)
        proportionalLabel((smInv6-0.1), (sfInv6+0.05), expression('inv6'), cex=1, adj=c(0.5, 0.5), xpd=NA, col='white')
    # Labels/Annotations
        proportionalLabel(0.03, 1.04, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, 'Local Demographic Rates', cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.4, 0.5, expression(italic(C) == 0), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)

## Panel B: Low-Elevation Perennials, Obligate Outcrossing
        # Make the plot
        LEP_C0 <- interp(x=heat_2$sm, y=heat_2$sf, z=heat_2$lNstar, 
                     xo=seq(min(heat_2$sm),max(heat_2$sm),by=resolution), 
                     yo=seq(min(heat_2$sf),max(heat_2$sf),by=resolution), duplicate="mean")
    # heatmap of lambda
        image(LEP_C0, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
        axis(1, at=c(0.01, 0.25, 0.5, 0.75, 0.99), labels=c("","","","",""))
        axis(2, at=c(0.01, 0.25, 0.5, 0.75, 0.99), labels=c("","","","",""))
    # Invasion Boundaries
        lines(AInvBound ~ sms, lwd=2, col=COLS$line, data=inv2)
        lines(aInvBound[aInvBound < max(sms)] ~ sms[aInvBound < max(sms)], lwd=2, col=COLS$line, data=inv2)
    # Plot inv6
        points(sfInv6 ~ smInv6, pch=21, col=COLS$line, bg='white', cex=1.25)
        proportionalLabel((smInv6-0.1), (sfInv6+0.05), expression('inv6'), cex=1, adj=c(0.5, 0.5), xpd=NA, col='white')
    # Labels/Annotations
        proportionalLabel(0.03, 1.04, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, 'Non-local Demographic Rates', cex=1.5, adj=c(0.5, 0.5), xpd=NA)

    # Color ScaleBar
        image.scale(LEP_C0, col=HEAT, breaks=breaks, horiz=FALSE, yaxt="n", xpd=TRUE)
        axis(4, at=c(-10, -5, 0, 5), labels=c("< -10","-5","0","5"), las=2, cex.axis=1.1)
        proportionalLabel(5, 0.5, expression(paste("Log Equilibrium Density (ln[", italic(N)^{"*"}, "])")), cex=1.25, adj=c(0.5, 0.5), xpd=NA, srt=270)



## Panel C: Eagle Meadows Population, Selfing w/ I.D>
    EM_C0.29 <- interp(x=heat_3$sm, y=heat_3$sf, z=heat_3$lNstar, 
                     xo=seq(min(heat_3$sm),max(heat_3$sm),by=resolution), 
                     yo=seq(min(heat_3$sf),max(heat_3$sf),by=resolution), duplicate="mean")
    # heatmap of lambda
        image(EM_C0.29, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
        axis(1, at=c(0.01, 0.25, 0.5, 0.75, 0.99), labels=c("0","0.25","0.5", "0.75", "1.0"))
        axis(2, at=c(0.01, 0.25, 0.5, 0.75, 0.99), labels=c("0","0.25","0.5", "0.75", "1.0"))
    # Invasion Boundaries
        lines(AInvBound ~ sms, lwd=2, col=COLS$line, data=inv3)
        lines(aInvBound[aInvBound < max(sms)] ~ sms[aInvBound < max(sms)], lwd=2, col=COLS$line, data=inv3)
    # Plot inv6
        points(sfInv6 ~ smInv6, pch=21, col=COLS$line, bg='white', cex=1.25)
        proportionalLabel((smInv6-0.1), (sfInv6+0.05), expression('inv6'), cex=1, adj=c(0.5, 0.5), xpd=NA, col='white')
    # Labels/Annotations
        proportionalLabel(0.03, 1.04, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(C), " = 0.29, w/ I.D.")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)

        proportionalLabel(-0.25, 1.125, expression(paste("Selection on female reproductive function (", italic(s[f]), ")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(1.125, -0.25, expression(paste("Selection on male reproductive function (", italic(s[m]), ")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)

## Panel D: Low-Elevation Perennials, Obligate Outcrossing
        # Make the plot
        LEP_C0.29 <- interp(x=heat_4$sm, y=heat_4$sf, z=heat_4$lNstar, 
                     xo=seq(min(heat_4$sm),max(heat_4$sm),by=resolution), 
                     yo=seq(min(heat_4$sf),max(heat_4$sf),by=resolution), duplicate="mean")
    # heatmap of lambda
        image(LEP_C0.29, col=HEAT, breaks=breaks, axes=FALSE)
    # axes
        axis(1, at=c(0.01, 0.25, 0.5, 0.75, 0.99), labels=c("0","0.25","0.5", "0.75", "1.0"))
        axis(2, at=c(0.01, 0.25, 0.5, 0.75, 0.99), labels=c("","","","",""))
    # Invasion Boundaries
        lines(AInvBound ~ sms, lwd=2, col=COLS$line, data=inv4)
        lines(aInvBound[aInvBound < max(sms)] ~ sms[aInvBound < max(sms)], lwd=2, col=COLS$line, data=inv4)
    # Plot inv6
        points(sfInv6 ~ smInv6, pch=21, col=COLS$line, bg='white', cex=1.25)
        proportionalLabel((smInv6-0.1), (sfInv6+0.05), expression('inv6'), cex=1, adj=c(0.5, 0.5), xpd=NA, col='white')
    # Labels/Annotations
        proportionalLabel(0.03, 1.04, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # Color ScaleBar
        image.scale(LEP_C0, col=HEAT, breaks=breaks, horiz=FALSE, yaxt="n", xpd=TRUE)
#        axis(4, las=2, cex.axis=1.1)
        axis(4, at=c(-10, -5, 0, 5), labels=c("< -10","-5","0","5"), las=2, cex.axis=1.1)
        proportionalLabel(5, 0.5, expression(paste("Log Equilibrium Density (ln[", italic(N)^{"*"}, "])")), cex=1.25, adj=c(0.5, 0.5), xpd=NA, srt=270)


}

##############################################################
##############################################################
##  Figures for Supplementary Material
##############################################################
##############################################################


demViablePolySpace_SexSpecDominance_Fig  <-  function() {

    # Import data sets
    invA       <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.25_hm0.75_C0_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extA_f5.8  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.75_C0_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extA_f6.0  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.75_C0_delta0_dj0_da0_dg0_f6.csv", head=TRUE)
    extA_f6.2  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.75_C0_delta0_dj0_da0_dg0_f6.2.csv", head=TRUE)

    invB       <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.25_hm0.75_C0.25_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extB_f5.8  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.75_C0.25_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extB_f6.0  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.75_C0.25_delta0_dj0_da0_dg0_f6.csv", head=TRUE)
    extB_f6.2  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.75_C0.25_delta0_dj0_da0_dg0_f6.2.csv", head=TRUE)

    invC       <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.25_hm0.75_C0.5_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extC_f5.8  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.75_C0.5_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extC_f6.0  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.75_C0.5_delta0_dj0_da0_dg0_f6.csv", head=TRUE)
    extC_f6.2  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.25_hm0.75_C0.5_delta0_dj0_da0_dg0_f6.2.csv", head=TRUE)

    invD       <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.75_hm0.25_C0_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extD_f5.8  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.75_hm0.25_C0_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extD_f6.0  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.75_hm0.25_C0_delta0_dj0_da0_dg0_f6.csv", head=TRUE)
    extD_f6.2  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.75_hm0.25_C0_delta0_dj0_da0_dg0_f6.2.csv", head=TRUE)

    invE       <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.75_hm0.25_C0.25_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extE_f5.8  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.75_hm0.25_C0.25_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extE_f6.0  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.75_hm0.25_C0.25_delta0_dj0_da0_dg0_f6.csv", head=TRUE)
    extE_f6.2  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.75_hm0.25_C0.25_delta0_dj0_da0_dg0_f6.2.csv", head=TRUE)

    invF       <-  read.csv("./output/simData/invasionBoundaries_sMax0.15_res0.0015_hf0.75_hm0.25_C0.5_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extF_f5.8  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.75_hm0.25_C0.5_delta0_dj0_da0_dg0_f5.8.csv", head=TRUE)
    extF_f6.0  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.75_hm0.25_C0.5_delta0_dj0_da0_dg0_f6.csv", head=TRUE)
    extF_f6.2  <-  read.csv("./output/simData/extThreshold_SfxSm_sMax0.15_res0.0015_hf0.75_hm0.25_C0.5_delta0_dj0_da0_dg0_f6.2.csv", head=TRUE)

    # Clean up some anomalous values from titrated thresholds
    extA_f5.8$sfThreshold[64]  <-  mean(c(extA_f5.8$sfThreshold[63], extA_f5.8$sfThreshold[65]))
    diff  <-  (extA_f5.8$smThreshold[92]-extA_f5.8$smThreshold[88])/4
    for (i in 89:91) {
        extA_f5.8$smThreshold[i]   <-  extA_f5.8$smThreshold[i-1] + diff
    }
    extA_f6.0$sfThreshold[46]  <-  mean(c(extA_f6.0$sfThreshold[45], extA_f6.0$sfThreshold[47]))
    extA_f6.0$sfThreshold[69]  <-  mean(c(extA_f6.0$sfThreshold[68], extA_f6.0$sfThreshold[70]))
    extA_f6.0$sfThreshold[76]  <-  mean(c(extA_f6.0$sfThreshold[75], extA_f6.0$sfThreshold[77]))
    extA_f6.0$smThreshold[88:91]
    diff  <-  (extA_f6.0$smThreshold[91]-extA_f6.0$smThreshold[88])/3
    for (i in 89:90) {
        extA_f6.0$smThreshold[i]   <-  extA_f6.0$smThreshold[i-1] + diff
    }
    extA_f6.2$sfThreshold[56]  <-  extA_f6.2$sfThreshold[57]
    extA_f6.2$sfThreshold[60]  <-  mean(c(extA_f6.2$sfThreshold[59], extA_f6.2$sfThreshold[61]))
    extA_f6.2$sfThreshold[95]  <-  mean(c(extA_f6.2$sfThreshold[94], extA_f6.2$sfThreshold[96]))
    diff  <-  (extA_f6.2$smThreshold[90]-extA_f6.2$smThreshold[88])/2
    extA_f6.2$smThreshold[89]  <-  extA_f6.2$smThreshold[88] + diff
    extB_f5.8$sfThreshold[29]  <-  mean(c(extB_f5.8$sfThreshold[28], extB_f5.8$sfThreshold[30]))
    extB_f5.8$smThreshold[57]  <-  mean(c(extB_f5.8$smThreshold[56], extB_f5.8$smThreshold[58]))
    extC_f5.8$sfThreshold[45]  <-  mean(c(extC_f5.8$sfThreshold[44], extC_f5.8$sfThreshold[46]))
    extC_f5.8$smThreshold[33]  <-  mean(c(extC_f5.8$smThreshold[32], extC_f5.8$smThreshold[34]))
    extD_f5.8$sfThreshold[69]  <-  mean(c(extD_f5.8$sfThreshold[68], extD_f5.8$sfThreshold[70]))
    extD_f5.8$sfThreshold[56]  <-  mean(c(extD_f5.8$sfThreshold[55], extD_f5.8$sfThreshold[57]))
    extD_f5.8$smThreshold[11]  <-  extD_f5.8$sfThreshold[10] + 0.001
    diff  <-  (extD_f5.8$smThreshold[31] - extD_f5.8$smThreshold[11])/20
    for (i in 12:30) {
        extD_f5.8$smThreshold[i]   <-  extD_f5.8$smThreshold[i-1] + diff
    }
    extD_f5.8$smThreshold[10]  <-  extD_f5.8$sfThreshold[11]
    diff  <-  (extD_f5.8$smThreshold[52] - extD_f5.8$smThreshold[47])/5
    for (i in 48:51) {
        extD_f5.8$smThreshold[i]   <-  extD_f5.8$smThreshold[i-1] + diff
    }
    diff  <-  (extD_f5.8$smThreshold[97] - extD_f5.8$smThreshold[88])/9
    for (i in 89:96) {
        extD_f5.8$smThreshold[i]   <-  extD_f5.8$smThreshold[i-1] + diff
    }
    extD_f6.0$smThreshold[78]  <-  mean(c(extD_f6.0$smThreshold[77],extD_f6.0$smThreshold[79]))
    diff  <-  (extD_f6.0$smThreshold[94] - extD_f6.0$smThreshold[84])/10
    for (i in 85:93) {
        extD_f6.0$smThreshold[i]   <-  extD_f6.0$smThreshold[i-1] + diff
    }
    extD_f6.2$sfThreshold[60]  <-  mean(c(extD_f6.2$sfThreshold[59],extD_f6.2$sfThreshold[61]))
    extD_f6.2$sfThreshold[95]  <-  mean(c(extD_f6.2$sfThreshold[94],extD_f6.2$sfThreshold[96]))
    diff  <-  (extD_f6.2$smThreshold[92] - extD_f6.2$smThreshold[87])/5
    for (i in 88:91) {
        extD_f6.2$smThreshold[i]   <-  extD_f6.2$smThreshold[i-1] + diff
    }
    diff  <-  (extE_f5.8$smThreshold[60] - extE_f5.8$smThreshold[55])/5
    for (i in 56:59) {
        extE_f5.8$smThreshold[i]   <-  extE_f5.8$smThreshold[i-1] + diff
    }
    extE_f6.0$smThreshold[57]  <-  mean(c(extE_f6.0$smThreshold[58],extE_f6.0$smThreshold[60]))
    diff  <-  (extE_f6.0$smThreshold[60] - extE_f6.0$smThreshold[55])/5
    for (i in 56:59) {
        extE_f6.0$smThreshold[i]   <-  extE_f6.0$smThreshold[i-1] + diff
    }
    extF_f5.8$sfThreshold[45]  <-  mean(c(extF_f5.8$sfThreshold[44], extF_f5.8$sfThreshold[46]))
    extF_f5.8$smThreshold[33]  <-  mean(c(extF_f5.8$smThreshold[32], extF_f5.8$smThreshold[34]))


# Color scheme
    COLS  <-  list("line"     =  transparentColor('#252525', opacity=1),
                   "extinct"  =  transparentColor('red', opacity=0.15))

# Set plot layout
layout.mat <- matrix(c(1:6), nrow=2, ncol=3, byrow=TRUE)
layout     <- layout(layout.mat,respect=TRUE)

##  Row 1: Additive allelic effects

    ##  Panel A: C = 0; hf = 3/4, hm = 1/4
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s',xpd=TRUE)
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Extinction Thresholds
        # f = 5.8 polygon
            t    <-  cbind(extA_f5.8$sf[!is.na(extA_f5.8$smThreshold)], extA_f5.8$smThreshold[!is.na(extA_f5.8$smThreshold)])
            t    <-  rbind(t, c(0.15, 0.15))
            b    <-  cbind(extA_f5.8$sms[!is.na(extA_f5.8$sfThreshold)], extA_f5.8$sfThreshold[!is.na(extA_f5.8$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extA_f5.8)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extA_f5.8)   
        proportionalLabel(0.85, 0.14, expression(paste(italic(f), " = 5.8")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # f = 6.0 polygon
            t    <-  cbind(extA_f6.0$sf[!is.na(extA_f6.0$smThreshold)], extA_f6.0$smThreshold[!is.na(extA_f6.0$smThreshold)])
            t    <-  rbind(t, c(0.15, 0.15))
            b    <-  cbind(extA_f6.0$sms[!is.na(extA_f6.0$sfThreshold)], extA_f6.0$sfThreshold[!is.na(extA_f6.0$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extA_f6.0)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extA_f6.0)   
        proportionalLabel(0.91, 0.35, expression(paste("6.0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # f = 6.2 polygon
            t    <-  cbind(extA_f6.2$sf[!is.na(extA_f6.2$smThreshold)], extA_f6.2$smThreshold[!is.na(extA_f6.2$smThreshold)])
            t    <-  rbind(t, c(0.15, 0.15))
            b    <-  cbind(extA_f6.2$sms[!is.na(extA_f6.2$sfThreshold)], extA_f6.2$sfThreshold[!is.na(extA_f6.2$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extA_f6.2)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extA_f6.2)   
        proportionalLabel(0.91, 0.54, expression(paste("6.2")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # Invasion Boundaries
        lines(aInvBound[1:88] ~ sms[1:88], lty=1, lwd=1.5, col=COLS$line, data=invA)
        lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invA)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.5, 0.5, expression(paste(italic(h[f]), " = 1/4, ", italic(h[m]), " = 3/4")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.3, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = 0')), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    ##  Panel B: C = 1/4; hf = 3/4, hm = 1/4
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2, bty='L')
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Extinction Thresholds
        # f = 5.8 polygon
            t    <-  cbind(extB_f5.8$sf[!is.na(extB_f5.8$smThreshold)], extB_f5.8$smThreshold[!is.na(extB_f5.8$smThreshold)])
            b    <-  cbind(extB_f5.8$sms[!is.na(extB_f5.8$sfThreshold)], extB_f5.8$sfThreshold[!is.na(extB_f5.8$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extB_f5.8)
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extB_f5.8)
        proportionalLabel(0.85, 0.14, expression(paste(italic(f), " = 5.8")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # f = 6.0 polygon
            t    <-  cbind(extB_f6.0$sf[!is.na(extB_f6.0$smThreshold)], extB_f6.0$smThreshold[!is.na(extB_f6.0$smThreshold)])
            b    <-  cbind(extB_f6.0$sms[!is.na(extB_f6.0$sfThreshold)], extB_f6.0$sfThreshold[!is.na(extB_f6.0$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extB_f6.0)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extB_f6.0)   
        proportionalLabel(0.91, 0.35, expression(paste("6.0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # f = 6.2 polygon
            t    <-  cbind(extB_f6.2$sf[!is.na(extB_f6.2$smThreshold)], extB_f6.2$smThreshold[!is.na(extB_f6.2$smThreshold)])
            b    <-  cbind(extB_f6.2$sms[!is.na(extB_f6.2$sfThreshold)], extB_f6.2$sfThreshold[!is.na(extB_f6.2$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extB_f6.2)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extB_f6.2)   
        proportionalLabel(1.1, 0.525, expression(paste("6.2")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalArrows(px1=1.1, py1=0.525, px2=0.97, py2=0.525, length=0.03)
        # Invasion Boundaries
        lines(aInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invB)
        lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invB)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = 1/4')), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    ##  Panel C: C = 1/2; hf = 3/4, hm = 1/4
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2, bty='L')
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Extinction Thresholds
        # f = 5.8 polygon
            t    <-  cbind(extC_f5.8$sf[!is.na(extC_f5.8$smThreshold)], extC_f5.8$smThreshold[!is.na(extC_f5.8$smThreshold)])
            b    <-  cbind(extC_f5.8$sms[!is.na(extC_f5.8$sfThreshold)], extC_f5.8$sfThreshold[!is.na(extC_f5.8$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extC_f5.8)
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extC_f5.8)
        proportionalLabel(0.85, 0.14, expression(paste(italic(f), " = 5.8")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # f = 6.0 polygon
            t    <-  cbind(extC_f6.0$sf[!is.na(extC_f6.0$smThreshold)], extC_f6.0$smThreshold[!is.na(extC_f6.0$smThreshold)])
            b    <-  cbind(extC_f6.0$sms[!is.na(extC_f6.0$sfThreshold)], extC_f6.0$sfThreshold[!is.na(extC_f6.0$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extC_f6.0)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extC_f6.0)   
        proportionalLabel(1.1, 0.32, expression(paste("6.0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalArrows(px1=1.1, py1=0.32, px2=0.97, py2=0.32, length=0.03)
        # Invasion Boundaries
        lines(aInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invC)
        lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invC)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = 1/2')), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.075, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)


##  Row 2:
  ##  Panel D: C = 0; hf = 1/4, hm = 3/4
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2, bty='L')
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Extinction Thresholds
        # f = 5.8 polygon
            t    <-  cbind(extD_f5.8$sf[!is.na(extD_f5.8$smThreshold)], extD_f5.8$smThreshold[!is.na(extD_f5.8$smThreshold)])
            t    <-  rbind(t, c(0.15, 0.15))
            b    <-  cbind(extD_f5.8$sms[!is.na(extD_f5.8$sfThreshold)], extD_f5.8$sfThreshold[!is.na(extD_f5.8$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extD_f5.8)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extD_f5.8)   
        proportionalLabel(0.85, 0.14, expression(paste(italic(f), " = 5.8")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # f = 6.0 polygon
            t    <-  cbind(extD_f6.0$sf[!is.na(extD_f6.0$smThreshold)], extD_f6.0$smThreshold[!is.na(extD_f6.0$smThreshold)])
            t    <-  rbind(t, c(0.15, 0.15))
            b    <-  cbind(extD_f6.0$sms, extD_f6.0$sfThreshold)
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extD_f6.0)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extD_f6.0)   
        proportionalLabel(0.91, 0.35, expression(paste("6.0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # f = 6.2 polygon
            t    <-  cbind(extD_f6.2$sf[!is.na(extD_f6.2$smThreshold)], extD_f6.2$smThreshold[!is.na(extD_f6.2$smThreshold)])
            t    <-  rbind(t, c(0.15, 0.15))
            b    <-  cbind(extD_f6.2$sms[!is.na(extD_f6.2$sfThreshold)], extD_f6.2$sfThreshold[!is.na(extD_f6.2$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extD_f6.2)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extD_f6.2)   
        proportionalLabel(0.91, 0.54, expression(paste("6.2")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # Invasion Boundaries
        lines(aInvBound[1:88] ~ sms[1:88], lty=1, lwd=1.5, col=COLS$line, data=invD)
        lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invD)
        # axes        
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.5, 0.5, expression(paste(italic(h[f]), " = 3/4, ", italic(h[m]), " = 1/4")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.3, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.04, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)


  ##  Panel E: C = 1/4; hf = 1/4, hm = 3/4
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2, bty='L')
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Extinction Thresholds
        # f = 5.8 polygon
            t    <-  cbind(extE_f5.8$sf[!is.na(extE_f5.8$smThreshold)], extE_f5.8$smThreshold[!is.na(extE_f5.8$smThreshold)])
            t    <-  rbind(t, c(0.102, 0.15))
            b    <-  cbind(extE_f5.8$sms[!is.na(extE_f5.8$sfThreshold)], extE_f5.8$sfThreshold[!is.na(extE_f5.8$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extE_f5.8)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extE_f5.8)   
        proportionalLabel(0.85, 0.14, expression(paste(italic(f), " = 5.8")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # f = 6.0 polygon
            t    <-  cbind(extE_f6.0$sf[!is.na(extE_f6.0$smThreshold)], extE_f6.0$smThreshold[!is.na(extE_f6.0$smThreshold)])
            t    <-  rbind(t, c(0.15, 0.15))
            b    <-  cbind(extE_f6.0$sms, extE_f6.0$sfThreshold)
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extE_f6.0)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extE_f6.0)   
        proportionalLabel(0.91, 0.35, expression(paste("6.0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # f = 6.2 polygon
            t    <-  cbind(extE_f6.2$sf[!is.na(extE_f6.2$smThreshold)], extE_f6.2$smThreshold[!is.na(extE_f6.2$smThreshold)])
            b    <-  cbind(extE_f6.2$sms[!is.na(extE_f6.2$sfThreshold)], extE_f6.2$sfThreshold[!is.na(extE_f6.2$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extE_f6.2)
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extE_f6.2)
        proportionalLabel(1.1, 0.525, expression(paste("6.2")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalArrows(px1=1.1, py1=0.525, px2=0.97, py2=0.525, length=0.03)
        # Invasion Boundaries
        lines(aInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invE)
        lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invE)
        # axes        
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.04, 1.075, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)

  ##  Panel F: C = 1/2; hf = 1/4, hm = 3/4
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2, bty='L')
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Extinction Thresholds
        # f = 5.8 polygon
            t    <-  cbind(extF_f5.8$sf[!is.na(extF_f5.8$smThreshold)], extF_f5.8$smThreshold[!is.na(extF_f5.8$smThreshold)])
            t    <-  rbind(t, c(0.056, 0.15))
            b    <-  cbind(extF_f5.8$sms[!is.na(extF_f5.8$sfThreshold)], extF_f5.8$sfThreshold[!is.na(extF_f5.8$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extF_f5.8)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extF_f5.8)   
        proportionalLabel(0.85, 0.14, expression(paste(italic(f), " = 5.8")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        # f = 6.0 polygon
            t    <-  cbind(extF_f6.0$sf[!is.na(extF_f6.0$smThreshold)], extF_f6.0$smThreshold[!is.na(extF_f6.0$smThreshold)])
            t    <-  rbind(t, c(0.048, 0.15))
            b    <-  cbind(extF_f6.0$sms[!is.na(extF_f6.0$sfThreshold)], extF_f6.0$sfThreshold[!is.na(extF_f6.0$sfThreshold)])
            polygon(c(t[,2],rev(b[,1])), c(t[,1],rev(b[,2])), col=COLS$extinct, border=NA)
        lines(sf ~ smThreshold, lty=1, lwd=0.5, col=COLS$line, data=extF_f6.0)   
        lines(sfThreshold ~ sms, lty=1, lwd=0.5, col=COLS$line, data=extF_f6.0)   
        proportionalLabel(1.1, 0.32, expression(paste("6.0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalArrows(px1=1.1, py1=0.32, px2=0.97, py2=0.32, length=0.03)
        # Invasion Boundaries
        lines(aInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invF)
        lines(AInvBound ~ sms, lty=1, lwd=1.5, col=COLS$line, data=invF)
        # axes        
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.04, 1.075, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
}






#' Fig showing Proportion of polymorphic parameter space
#' (w/out inbreeding depression)
polySpaceFigTitrateSexSpec  <-  function(df = "dataPolySpaceFig_hf0.25_0.75_hm0.75_0.25_sMax0.15_res0.003_delta0_dj0_da0_dg0") {

    # Make filenames for import from df names
    fName  <-  paste('./output/simData/', df, '.csv', sep="")

    # Extract plotting parameter values from df names
    df1  <-  strsplit(df, '_')[[1]][c(2:11)]
    pars1  <-  list(
                    "hf"    =  as.numeric(strsplit(df1[1],'f')[[1]][2]),
                    "hm"    =  as.numeric(strsplit(df1[3],'m')[[1]][2]),
                    "sMax"  =  as.numeric(strsplit(df1[5],'x')[[1]][2]),
                    "res"   =  as.numeric(strsplit(df1[6],'s')[[1]][2]),
                    "d"     =  as.numeric(strsplit(df1[7],'a')[[1]][2]),
                    "dj"    =  as.numeric(strsplit(df1[8],'j')[[1]][2]),
                    "da"    =  as.numeric(strsplit(df1[9],'a')[[1]][2]),
                    "dg"    =  as.numeric(strsplit(df1[10],'g')[[1]][2])
                    )
    df2  <-  strsplit(df, '_')[[1]][c(2:11)]
    pars2  <-  list(
                    "hf"    =  as.numeric(df2[2]),
                    "hm"    =  as.numeric(df2[4]),
                    "sMax"  =  as.numeric(strsplit(df2[5],'x')[[1]][2]),
                    "res"   =  as.numeric(strsplit(df2[6],'s')[[1]][2]),
                    "d"     =  as.numeric(strsplit(df2[7],'a')[[1]][2]),
                    "dj"    =  as.numeric(strsplit(df2[8],'j')[[1]][2]),
                    "da"    =  as.numeric(strsplit(df2[9],'a')[[1]][2]),
                    "dg"    =  as.numeric(strsplit(df2[10],'g')[[1]][2])
                    )

    # import data
    data  <-  read.csv(file=fName, header=TRUE)

    # clean data set & quantify parameter space
    dat1   <-  quantPolySpaceSexSpec(data = data, pars = pars1)
    dat2   <-  quantPolySpaceSexSpec(data = data, pars = pars2)

    # Color scheme
    COLS  <-  list(
                    "PG"    =  transparentColor('#252525', opacity=1),
                    "low"   =  transparentColor('dodgerblue4', opacity=0.6),
                    "med"   =  transparentColor('darkolivegreen4', opacity=0.6),
                    "hi"    =  transparentColor('tomato', opacity=0.6),
                    "low2"  =  transparentColor('dodgerblue4', opacity=1),
                    "med2"  =  transparentColor('darkolivegreen4', opacity=1),
                    "hi2"   =  transparentColor('tomato', opacity=1)
                   )

#  Create vector of selfing rates for pop gen function.
    CLine     <-  seq(0,0.9,length=100)
    FrecMdom  <-  c()
    FdomMrec  <-  c()
    for(i in 1:length(CLine)) {
        FrecMdom[i]  <-  popGen_PolySpace(hf=0.25, hm=0.75, C=CLine[i], sMax=pars1$sMax)
        FdomMrec[i]  <-  popGen_PolySpace(hf=0.75, hm=0.25, C=CLine[i], sMax=pars1$sMax)
    }

# Set plot layout
    layout.mat  <- matrix(c(1,2), nrow=2, ncol=1, byrow=TRUE)
    layout      <- layout(layout.mat,respect=TRUE)

    ##  Panel A: hf = 1/4; hm = 3/4
        # subset data
        d  <-  dat1[dat1$hf == 0.25,]
        dlow  <-  d[d$f == 5.8,]
        dmed  <-  d[d$f == 6.0,]
        dhi   <-  d[d$f == 6.5,]
        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(d$C)), ylim = c(0,0.105), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(FrecMdom ~ CLine, lwd=2, col=COLS$PG)
        points(PrViaPoly ~ C, pch=21, bg=COLS$low, col=COLS$low2, data=dlow)
        points(PrViaPoly ~ C, pch=21, bg=COLS$med, col=COLS$med2, data=dmed)
        points(PrViaPoly ~ C, pch=21, bg=COLS$hi, col=COLS$hi2, data=dhi)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.1, substitute(paste(italic(h), " = ", hh), list(hh = pars2$hf)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.5, 1.1, expression(paste(italic(h[f]), " = ", 1/4, "; ", italic(h[m]), " = ", 3/4)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(-0.3, -0.16, expression(paste("Proportion viable polymorphic parameter space")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=90)
        #Legend
        legend( x       =  usr[2],
                y       =  usr[4],
                legend  =  c(
                             expression(paste("Pop. Gen.")),
                             substitute(paste("High fertility (", italic(f), " = ", ff, ")"), list(ff = dhi$f[1])),
                             substitute(paste("Med. fertility (", italic(f), " = ", ff, ".0)"), list(ff = dmed$f[1])),
                             substitute(paste("Low  fertility (", italic(f), " = ", ff, ")"), list(ff = dlow$f[1]))),
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
        # subset data
        d  <-  dat1[dat1$hf == 0.75,]
        dlow  <-  d[d$f == 5.8,]
        dmed  <-  d[d$f == 5.9,]
        dhi   <-  d[d$f == 6.5,]
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(d$C)), ylim = c(0,0.105), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(FdomMrec ~ CLine, lwd=2, col=COLS$PG)
        points(PrViaPoly ~ C, pch=21, bg=COLS$low, col=COLS$low2, data=dlow)
        points(PrViaPoly ~ C, pch=21, bg=COLS$med, col=COLS$med2, data=dmed)
        points(PrViaPoly ~ C, pch=21, bg=COLS$hi, col=COLS$hi2, data=dhi)
        # axes        
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(0.03, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.1, expression(paste(italic(h[f]), " = ", 3/4, "; ", italic(h[m]), " = ", 1/4)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.5, -0.3, expression(paste("Selfing rate (",italic(C), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
#        proportionalLabel(-0.3, 0.5, expression(paste("Prop. viable polymorphic space")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=90)
}





#' Fig showing Effect of Inbreeding Depression on 
#' Proportion of polymorphic parameter space

deltaSelfingLoadPolySpaceFigTitrateSexSpec  <-  function(df1 = "dataDeltaPolySpaceFigSexSpecDom_sMax0.15_res0.003_dStar0.8_f6.5",
                                                         df2 = "dataDeltaPolySpaceFigSexSpecDom_sMax0.15_res0.003_dStar0.8_f7.5",
                                                         df3 = "dataDeltaPolySpaceFigSexSpecDom_sMax0.15_res0.003_dStar0.8_f8.5") {

    # Make filenames for import from df names
    fName1  <-  paste('./output/simData/', df1, '.csv', sep="")
    fName2  <-  paste('./output/simData/', df2, '.csv', sep="")
    fName3  <-  paste('./output/simData/', df3, '.csv', sep="")

    # import data
    data1  <-  read.csv(file=fName1, header=TRUE)
    data2  <-  read.csv(file=fName2, header=TRUE)
    data3  <-  read.csv(file=fName3, header=TRUE)

    # Extract plotting parameter values from df names
    d1   <-  strsplit(df1, '_')[[1]][c(2:5)]
    pars1  <-  list(
                    "sMax"  =  as.numeric(strsplit(d1[1],'x')[[1]][2]),
                    "res"   =  as.numeric(strsplit(d1[2],'s')[[1]][2]),
                    "dStar" =  as.numeric(strsplit(d1[3],'r')[[1]][2]),
                    "f"     =  as.numeric(strsplit(d1[4],'f')[[1]][2])
                    )
    d2   <-  strsplit(df2, '_')[[1]][c(2:5)]
    pars2  <-  list(
                    "sMax"  =  as.numeric(strsplit(d2[1],'x')[[1]][2]),
                    "res"   =  as.numeric(strsplit(d2[2],'s')[[1]][2]),
                    "dStar" =  as.numeric(strsplit(d2[3],'r')[[1]][2]),
                    "f"     =  as.numeric(strsplit(d2[4],'f')[[1]][2])
                    )
    d3   <-  strsplit(df3, '_')[[1]][c(2:5)]
    pars3  <-  list(
                    "sMax"  =  as.numeric(strsplit(d2[1],'x')[[1]][2]),
                    "res"   =  as.numeric(strsplit(d2[2],'s')[[1]][2]),
                    "dStar" =  as.numeric(strsplit(d2[3],'r')[[1]][2]),
                    "f"     =  as.numeric(strsplit(d2[4],'f')[[1]][2])
                    )
    hLev  <-  unique(data1$hf)
    dLev  <-  unique(data1$Delta)
    CLev  <-  unique(data1$C)
    nHs   <-  length(hLev)
    nDs   <-  length(dLev)
    nCs   <-  length(CLev)

    # clean data set & quantify parameter space
    dat1   <-  quantDeltaPolySpace(data = data1, pars = pars1, sexSpec = TRUE)
    dat1[dat1 < 0]  <-  0
    dat2   <-  quantDeltaPolySpace(data = data2, pars = pars2, sexSpec = TRUE)
    dat2[dat2 < 0]  <-  0
    dat3   <-  quantDeltaPolySpace(data = data3, pars = pars3, sexSpec = TRUE)
    dat3[dat3 < 0]  <-  0

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
    ##          early-acting delta
        PGSpace     <-  c()
        for(i in 1:length(CLine)) {
            PGSpace[i]     <-  popGen_PolySpace_Delta_Add(C=CLine[i], delta=dLine[i], sMax=pars1$sMax)
        }
        d1  <-  dat1[dat1$h == hLev[1],]
        d2  <-  dat2[dat2$h == hLev[1],]
        d2$PrViaPoly[d2$Delta == "d_g"][c(18)]  <-  as.character((as.numeric(d2$PrViaPoly[d2$Delta == "d_g"][c(17)])  + as.numeric(d2$PrViaPoly[d2$Delta == "d_g"][c(19)]) )/2)
        d3  <-  dat3[dat3$h == hLev[1],]
        d3$PrViaPoly[d3$Delta == "d_g"][c(18)]  <-  as.character((as.numeric(d3$PrViaPoly[d3$Delta == "d_g"][c(17)])  + as.numeric(d3$PrViaPoly[d3$Delta == "d_g"][c(19)]) )/2)

        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,4,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.925), ylim = c(0,0.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        # f = 6.5
        points(d1$PrViaPoly[d1$Delta == "d"][c(1:9,23)] ~ C[Delta == "d"][c(1:9,23)], pch=21, bg=COLS$dSim, col=COLS$dSim2, data=d1)
        points(PrViaPoly[Delta == "d_j"][c(1:7,23)] ~ C[Delta == "d_j"][c(1:7,23)], pch=21, bg=COLS$d_j, col=COLS$d_j2, data=d1)
        points(PrViaPoly[Delta == "d_a"][c(1:12,23)] ~ C[Delta == "d_a"][c(1:12,23)], pch=21, bg=COLS$d_a, col=COLS$d_a2, data=d1)
        points(PrViaPoly[Delta == "d_g"][c(1:10,23)] ~ C[Delta == "d_g"][c(1:10,23)], pch=21, bg=COLS$d_g, col=COLS$d_g2, data=d1)
        # f = 7.5
        points(PrViaPoly[Delta == "d"][c(1:21,37)] ~ C[Delta == "d"][c(1:21,37)], pch=3, col=COLS$dSim2, data=d2)
        points(PrViaPoly[Delta == "d_j"][c(1:15,37)] ~ C[Delta == "d_j"][c(1:15,37)], pch=3, col=COLS$d_j2, data=d2)
        points(PrViaPoly[Delta == "d_a"][c(1:24,37)] ~ C[Delta == "d_a"][c(1:24,37)], pch=3, col=COLS$d_a2, data=d2)
        points(PrViaPoly[Delta == "d_g"][c(1:22,37)] ~ C[Delta == "d_g"][c(1:22,37)], pch=3, col=COLS$d_g2, data=d2)
        # f = 8.5
        points(PrViaPoly[Delta == "d"][c(1:32,37)] ~ C[Delta == "d"][c(1:32,37)], pch=4, col=COLS$dSim2, data=d3)
        points(PrViaPoly[Delta == "d_j"][c(1:21,37)] ~ C[Delta == "d_j"][c(1:21,37)], pch=4, col=COLS$d_j2, data=d3)
        points(PrViaPoly[Delta == "d_a"][c(1:35,37)] ~ C[Delta == "d_a"][c(1:35,37)], pch=4, col=COLS$d_a2, data=d3)
        points(PrViaPoly[Delta == "d_g"][c(1:33,37)] ~ C[Delta == "d_g"][c(1:33,37)], pch=4, col=COLS$d_g2, data=d3)
        # PG expectation
        lines(PGSpace ~ CLine, lwd=2, col=COLS$PG)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Labels/annotations
        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(italic(h[f]), " = ", 1/4, ", ", italic(h[m]), " = ", 3/4)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(-0.35, -0.16, expression(paste("Proportion viable polymorphic parameter space")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        legend( x       =  usr[2],
                y       =  usr[4],
                legend  =  c(
                             expression(paste(italic(delta), " (pop. gen.)")),
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
        legend( x       =  usr[2]*0.93,
                y       =  usr[4]*0.62,
                legend  =  c(
                             expression(paste(italic(f), " = ", 8.5)),
                             expression(paste(italic(f), " = ", 7.5)),
                             expression(paste(italic(f), " = ", 6.5))),
                 col     =  c(COLS$PG),
                 pch     =  c(4,3,21),
                 pt.bg   =  c(NA),
                 cex     =  0.75,
                 pt.cex  =  0.75,
                 xjust   =  1,
                 yjust   =  1,
                 bty     =  'n',
                 border  =  NA)


    ## Panel B: Dominance Reversal SA (hf = hm = 1/4)
    ##          early-acting delta
        CLine2  <-  CLine[-1]
        dLine2  <-  dLine[-1]
        PGSpace     <-  c()
        for(i in 1:length(CLine2)) {
            PGSpace[i]     <-  popGen_PolySpace_Delta_Add(C=CLine[i], delta=dLine[i], sMax=pars1$sMax)
        }
        d1  <-  dat1[dat1$h == 0.25,]
        d2  <-  dat2[dat2$h == 0.25,]
        d3  <-  dat3[dat3$h == 0.25,]
        d1$PrViaPoly[30]  <-  "0"
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='', xlim = c(0,0.925), ylim = c(0,0.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        # f = 6.5
        points(PrViaPoly[Delta == "d"][c(1:10,23)] ~ C[Delta == "d"][c(1:10,23)], pch=21, bg=COLS$dSim, col=COLS$dSim2, data=d1)
        points(PrViaPoly[Delta == "d_j"][c(1:7,23)] ~ C[Delta == "d_j"][c(1:7,23)], pch=21, bg=COLS$d_j, col=COLS$d_j2, data=d1)
        points(PrViaPoly[Delta == "d_a"][c(1:12,23)] ~ C[Delta == "d_a"][c(1:12,23)], pch=21, bg=COLS$d_a, col=COLS$d_a2, data=d1)
        points(PrViaPoly[Delta == "d_g"][c(1:10,23)] ~ C[Delta == "d_g"][c(1:10,23)], pch=21, bg=COLS$d_g, col=COLS$d_g2, data=d1)
        # f = 7.5
        points(PrViaPoly[Delta == "d"][c(1:22,37)] ~ C[Delta == "d"][c(1:22,37)], pch=3, bg=COLS$dSim, col=COLS$dSim2, data=d2)
        points(PrViaPoly[Delta == "d_j"][c(1:15,37)] ~ C[Delta == "d_j"][c(1:15,37)], pch=3, bg=COLS$d_j, col=COLS$d_j2, data=d2)
        points(PrViaPoly[Delta == "d_a"][c(1:25,37)] ~ C[Delta == "d_a"][c(1:25,37)], pch=3, bg=COLS$d_a, col=COLS$d_a2, data=d2)
        points(PrViaPoly[Delta == "d_g"][c(1:23,37)] ~ C[Delta == "d_g"][c(1:23,37)], pch=3, bg=COLS$d_g, col=COLS$d_g2, data=d2)        # axes        
        # f = 8.5
        points(PrViaPoly[Delta == "d"][c(1:32,37)] ~ C[Delta == "d"][c(1:32,37)], pch=4, col=COLS$dSim2, data=d3)
        points(PrViaPoly[Delta == "d_j"][c(1:21,37)] ~ C[Delta == "d_j"][c(1:21,37)], pch=4, col=COLS$d_j2, data=d3)
        points(PrViaPoly[Delta == "d_a"][c(1:35,37)] ~ C[Delta == "d_a"][c(1:35,37)], pch=4, col=COLS$d_a2, data=d3)
        points(PrViaPoly[Delta == "d_g"][c(1:34,37)] ~ C[Delta == "d_g"][c(1:34,37)], pch=4, col=COLS$d_g2, data=d3)
        # Pop. Gen Prediction
        lines(PGSpace ~ CLine2, lwd=2, col=COLS$PG)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Labels/annotations
        proportionalLabel(0.03, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(italic(h[f]), " = ", 3/4, ", ", italic(h[m]), " = ", 1/4)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
#        proportionalLabel(-0.35, 0.5, expression(paste("Prop. viable polymorphic space")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.35, expression(paste("Selfing Rate (", italic(C), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)

}












suppPolySpaceThresholdFigs  <-  function(df = "dataPolySpaceFig_sMax0.15_res0.003_delta0_dj0_da0_dg0") {

    # Make filenames for import from df names
    fName  <-  paste('./output/simData/', df, '.csv', sep="")

    # import data
    data  <-  read.csv(file=fName, header=TRUE)

    # Clean up aberrant smExt value
    data  <-  cleanPolySpaceData(data = data)
    data$smExt[90]  <-  mean(c(data$smExt[89], data$smExt[91]))

    # Extract plotting parameter values from df names
    d1   <-  strsplit(df, '_')[[1]][c(2:7)]
    pars  <-  list(
                    "sMax"  =  as.numeric(strsplit(d1[1],'x')[[1]][2]),
                    "res"   =  as.numeric(strsplit(d1[2],'s')[[1]][2]),
                    "d"     =  as.numeric(strsplit(d1[3],'a')[[1]][2]),
                    "dj"    =  as.numeric(strsplit(d1[4],'j')[[1]][2]),
                    "da"    =  as.numeric(strsplit(d1[5],'a')[[1]][2]),
                    "dg"    =  as.numeric(strsplit(d1[6],'g')[[1]][2])
                    )

    hLev  <-  unique(data$h)
    fLev  <-  unique(data$f)
    CLev  <-  unique(data$C)
    nHs   <-  length(hLev)
    nCs   <-  length(CLev)

    # Color scheme
    COLS  <-  list("line"     =  transparentColor('#252525', opacity=1),
                   "extinct"  =  transparentColor('red', opacity=0.15))

    # Set plot layout
    nMat        <-  nHs*nCs
    layout.mat  <- matrix(c(1:nMat), nrow=nHs, ncol=nCs, byrow=TRUE)
    layout      <- layout(layout.mat,respect=TRUE)

    # loop through and generate plots
    for(i in 1:nHs) {
        # Subset data
        dd  <-  data[data$h == hLev[i],]            
        
        for(j in 1:nCs) {
            
            d  <-  dd[dd$C == CLev[j],]            
            if(i == 1 && j == 1){
                par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s',xpd=TRUE)
            }
            plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2)
            usr  <-  par('usr')
            rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
            plotGrid(lineCol='grey80')
            box()
            # Plot Thresholds and Invasion boundaries
            for(k in 1:length(unique(d$f))) {
                # Extinction Thresholds
                lines(sf ~ smExt, lty=1, lwd=0.75, col=COLS$line, data=d[d$f == fLev[k],])
                lines(sfExt ~ sm, lty=1, lwd=0.75, col=COLS$line, data=d[d$f == fLev[k],])
                fTextLoc  <-  ((min(d$sfExt[d$f == fLev[k]][pars$sMax/pars$res])*0.97)/(pars$sMax+abs(2*usr[1])))
                proportionalLabel(0.85, fTextLoc, substitute(paste(italic(f), " = ", ff), list(ff = fLev[k])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
                if(k == 1) {
                    # Invasion Boundaries
                    lines(aInv[aInv < pars$sMax] ~ smInv[aInv < pars$sMax], lty=1, lwd=2, col=COLS$line, data=d[d$f == fLev[k],])
                    lines(AInv ~ smInv, lty=1, lwd=2, col=COLS$line, data=d[d$f == fLev[k],])
                }
            }
            axis(1, las=1, labels=NA)
            axis(2, las=1, labels=NA)
            # labels & annotations
            if(i == 1 && j == 1){
                axis(2, las=1)
                proportionalLabel(-0.3, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
            }
            if(i == 2 && j == 1){
                axis(2, las=1)
                proportionalLabel(-0.3, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
            }
            if(i == 1) {
                proportionalLabel(0.5, 1.15, substitute(paste(italic(C), " = ", CC), list(CC = CLev[j])), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
            }
            if(j==1) {
                proportionalLabel(-0.5, 0.5, substitute(paste(italic(h), " = ", hh), list(hh = hLev[i])), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
            }
            if(i == 2) {
                axis(1, las=1)
                proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
            }
        }
    }
 
}




#' Fig showing Effect of Inbreeding Depression on 
#' Proportion of polymorphic parameter space

suppDeltaSelfingLoadPolySpaceFigTitrate  <-  function(df1 = "dataDeltaPolySpaceFigSexSpecDom_sMax0.15_res0.003_dStar0.8_f6.5",
                                                      df2 = "dataDeltaPolySpaceFigSexSpecDom_sMax0.15_res0.003_dStar0.8_f7.5",
                                                      df3 = "dataDeltaPolySpaceFigSexSpecDom_sMax0.15_res0.003_dStar0.8_f8.5") {

    # Make filenames for import from df names
    fName1  <-  paste('./output/simData/', df1, '.csv', sep="")
    fName2  <-  paste('./output/simData/', df2, '.csv', sep="")
    fName3  <-  paste('./output/simData/', df3, '.csv', sep="")

    # import data
    data1  <-  read.csv(file=fName1, header=TRUE)
    data2  <-  read.csv(file=fName2, header=TRUE)
    data3  <-  read.csv(file=fName3, header=TRUE)

    # Extract plotting parameter values from df names
    d1   <-  strsplit(df1, '_')[[1]][c(2:5)]
    pars1  <-  list(
                    "sMax"  =  as.numeric(strsplit(d1[1],'x')[[1]][2]),
                    "res"   =  as.numeric(strsplit(d1[2],'s')[[1]][2]),
                    "dStar" =  as.numeric(strsplit(d1[3],'r')[[1]][2]),
                    "f"     =  as.numeric(strsplit(d1[4],'f')[[1]][2])
                    )
    d2   <-  strsplit(df2, '_')[[1]][c(2:5)]
    pars2  <-  list(
                    "sMax"  =  as.numeric(strsplit(d2[1],'x')[[1]][2]),
                    "res"   =  as.numeric(strsplit(d2[2],'s')[[1]][2]),
                    "dStar" =  as.numeric(strsplit(d2[3],'r')[[1]][2]),
                    "f"     =  as.numeric(strsplit(d2[4],'f')[[1]][2])
                    )
    d3   <-  strsplit(df3, '_')[[1]][c(2:5)]
    pars3  <-  list(
                    "sMax"  =  as.numeric(strsplit(d2[1],'x')[[1]][2]),
                    "res"   =  as.numeric(strsplit(d2[2],'s')[[1]][2]),
                    "dStar" =  as.numeric(strsplit(d2[3],'r')[[1]][2]),
                    "f"     =  as.numeric(strsplit(d2[4],'f')[[1]][2])
                    )
    hLev  <-  unique(data1$h)
    dLev  <-  unique(data1$Delta)
    CLev  <-  unique(data1$C)
    nHs   <-  length(hLev)
    nDs   <-  length(dLev)
    nCs   <-  length(CLev)

    # clean data set & quantify parameter space
    dat1   <-  quantPolySpaceSexSpec(data = data1, pars = pars1)
    dat1[dat1 < 0]  <-  0
    dat2   <-  quantPolySpaceSexSpec(data = data2, pars = pars2)
    dat2[dat2 < 0]  <-  0
    dat3   <-  quantPolySpaceSexSpec(data = data3, pars = pars3)
    dat3[dat3 < 0]  <-  0

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
    ##          early-acting delta
        PGSpace     <-  c()
        for(i in 1:length(CLine)) {
            PGSpace[i]     <-  popGen_PolySpace_Delta_Add(C=CLine[i], delta=dLine[i], sMax=pars1$sMax)
        }
        d1  <-  dat1[dat1$h == 0.5,]
        d2  <-  dat2[dat2$h == 0.5,]
        d2$PrViaPoly[d2$Delta == "d_g"][c(18)]  <-  as.character((as.numeric(d2$PrViaPoly[d2$Delta == "d_g"][c(17)])  + as.numeric(d2$PrViaPoly[d2$Delta == "d_g"][c(19)]) )/2)
        d3  <-  dat3[dat3$h == 0.5,]
        d3$PrViaPoly[d3$Delta == "d_g"][c(18)]  <-  as.character((as.numeric(d3$PrViaPoly[d3$Delta == "d_g"][c(17)])  + as.numeric(d3$PrViaPoly[d3$Delta == "d_g"][c(19)]) )/2)

        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,4,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.925), ylim = c(0,0.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        # f = 6.5
        points(PrViaPoly[Delta == "d"][c(1:9,23)] ~ C[Delta == "d"][c(1:9,23)], pch=21, bg=COLS$dSim, col=COLS$dSim2, data=d1)
        points(PrViaPoly[Delta == "d_j"][c(1:7,23)] ~ C[Delta == "d_j"][c(1:7,23)], pch=21, bg=COLS$d_j, col=COLS$d_j2, data=d1)
        points(PrViaPoly[Delta == "d_a"][c(1:12,23)] ~ C[Delta == "d_a"][c(1:12,23)], pch=21, bg=COLS$d_a, col=COLS$d_a2, data=d1)
        points(PrViaPoly[Delta == "d_g"][c(1:10,23)] ~ C[Delta == "d_g"][c(1:10,23)], pch=21, bg=COLS$d_g, col=COLS$d_g2, data=d1)
        # f = 7.5
        points(PrViaPoly[Delta == "d"][c(1:21,37)] ~ C[Delta == "d"][c(1:21,37)], pch=3, col=COLS$dSim2, data=d2)
        points(PrViaPoly[Delta == "d_j"][c(1:15,37)] ~ C[Delta == "d_j"][c(1:15,37)], pch=3, col=COLS$d_j2, data=d2)
        points(PrViaPoly[Delta == "d_a"][c(1:24,37)] ~ C[Delta == "d_a"][c(1:24,37)], pch=3, col=COLS$d_a2, data=d2)
        points(PrViaPoly[Delta == "d_g"][c(1:22,37)] ~ C[Delta == "d_g"][c(1:22,37)], pch=3, col=COLS$d_g2, data=d2)
        # f = 8.5
        points(PrViaPoly[Delta == "d"][c(1:32,37)] ~ C[Delta == "d"][c(1:32,37)], pch=4, col=COLS$dSim2, data=d3)
        points(PrViaPoly[Delta == "d_j"][c(1:21,37)] ~ C[Delta == "d_j"][c(1:21,37)], pch=4, col=COLS$d_j2, data=d3)
        points(PrViaPoly[Delta == "d_a"][c(1:35,37)] ~ C[Delta == "d_a"][c(1:35,37)], pch=4, col=COLS$d_a2, data=d3)
        points(PrViaPoly[Delta == "d_g"][c(1:33,37)] ~ C[Delta == "d_g"][c(1:33,37)], pch=4, col=COLS$d_g2, data=d3)
        # PG expectation
        lines(PGSpace ~ CLine, lwd=2, col=COLS$PG)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Labels/annotations
        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(italic(h), " = ", 1/2)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(-0.35, -0.16, expression(paste("Proportion viable polymorphic parameter space")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        legend( x       =  usr[2],
                y       =  usr[4],
                legend  =  c(
                             expression(paste(italic(delta), " (pop. gen.)")),
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
        legend( x       =  usr[2]*0.93,
                y       =  usr[4]*0.62,
                legend  =  c(
                             expression(paste(italic(f), " = ", 8.5)),
                             expression(paste(italic(f), " = ", 7.5)),
                             expression(paste(italic(f), " = ", 6.5))),
                 col     =  c(COLS$PG),
                 pch     =  c(4,3,21),
                 pt.bg   =  c(NA),
                 cex     =  0.75,
                 pt.cex  =  0.75,
                 xjust   =  1,
                 yjust   =  1,
                 bty     =  'n',
                 border  =  NA)


    ## Panel B: Dominance Reversal SA (hf = hm = 1/4)
    ##          early-acting delta
        CLine2  <-  CLine[-1]
        dLine2  <-  dLine[-1]
        PGSpace     <-  c()
        for(i in 1:length(CLine2)) {
            PGSpace[i]     <-  popGen_PolySpace_Delta_DomRev(C=CLine2[i], delta=dLine2[i], sMax=pars2$sMax)
        }
        d1  <-  dat1[dat1$h == 0.25,]
        d2  <-  dat2[dat2$h == 0.25,]
        d3  <-  dat3[dat3$h == 0.25,]
        d1$PrViaPoly[30]  <-  "0"
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='', xlim = c(0,0.925), ylim = c(0,0.8), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        # f = 6.5
        points(PrViaPoly[Delta == "d"][c(1:10,23)] ~ C[Delta == "d"][c(1:10,23)], pch=21, bg=COLS$dSim, col=COLS$dSim2, data=d1)
        points(PrViaPoly[Delta == "d_j"][c(1:7,23)] ~ C[Delta == "d_j"][c(1:7,23)], pch=21, bg=COLS$d_j, col=COLS$d_j2, data=d1)
        points(PrViaPoly[Delta == "d_a"][c(1:12,23)] ~ C[Delta == "d_a"][c(1:12,23)], pch=21, bg=COLS$d_a, col=COLS$d_a2, data=d1)
        points(PrViaPoly[Delta == "d_g"][c(1:10,23)] ~ C[Delta == "d_g"][c(1:10,23)], pch=21, bg=COLS$d_g, col=COLS$d_g2, data=d1)
        # f = 7.5
        points(PrViaPoly[Delta == "d"][c(1:22,37)] ~ C[Delta == "d"][c(1:22,37)], pch=3, bg=COLS$dSim, col=COLS$dSim2, data=d2)
        points(PrViaPoly[Delta == "d_j"][c(1:15,37)] ~ C[Delta == "d_j"][c(1:15,37)], pch=3, bg=COLS$d_j, col=COLS$d_j2, data=d2)
        points(PrViaPoly[Delta == "d_a"][c(1:25,37)] ~ C[Delta == "d_a"][c(1:25,37)], pch=3, bg=COLS$d_a, col=COLS$d_a2, data=d2)
        points(PrViaPoly[Delta == "d_g"][c(1:23,37)] ~ C[Delta == "d_g"][c(1:23,37)], pch=3, bg=COLS$d_g, col=COLS$d_g2, data=d2)        # axes        
        # f = 8.5
        points(PrViaPoly[Delta == "d"][c(1:32,37)] ~ C[Delta == "d"][c(1:32,37)], pch=4, col=COLS$dSim2, data=d3)
        points(PrViaPoly[Delta == "d_j"][c(1:21,37)] ~ C[Delta == "d_j"][c(1:21,37)], pch=4, col=COLS$d_j2, data=d3)
        points(PrViaPoly[Delta == "d_a"][c(1:35,37)] ~ C[Delta == "d_a"][c(1:35,37)], pch=4, col=COLS$d_a2, data=d3)
        points(PrViaPoly[Delta == "d_g"][c(1:34,37)] ~ C[Delta == "d_g"][c(1:34,37)], pch=4, col=COLS$d_g2, data=d3)
        # Pop. Gen Prediction
        lines(PGSpace ~ CLine2, lwd=2, col=COLS$PG)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Labels/annotations
        proportionalLabel(0.03, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(italic(h), " = ", 1/4)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
#        proportionalLabel(-0.35, 0.5, expression(paste("Prop. viable polymorphic space")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.35, expression(paste("Selfing Rate (", italic(C), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)

}



















suppDeltaPolySpaceThresholdFigs  <-  function(df = "dataDeltaPolySpaceFig_sMax0.15_res0.003_dStar0.8_f6.5") {

    # Make filenames for import from df names
    fName  <-  paste('./output/simData/', df, '.csv', sep="")

    # import data
    data  <-  read.csv(file=fName, header=TRUE)

    # Clean up aberrant smExt value
    data  <-  cleanDeltaPolySpaceData(data = data)
    data$smExt[6530:6544]  <-  0
    data$sfExt[6578:6579]  <-  0

    # Extract plotting parameter values from df names
    d1   <-  strsplit(df, '_')[[1]][c(2:5)]
    pars  <-  list(
                    "sMax"  =  as.numeric(strsplit(d1[1],'x')[[1]][2]),
                    "res"   =  as.numeric(strsplit(d1[2],'s')[[1]][2]),
                    "dStar" =  as.numeric(strsplit(d1[3],'r')[[1]][2]),
                    "f"     =  as.numeric(strsplit(d1[4],'f')[[1]][2])
                    )

    hLev  <-  unique(data$h)
    dLev  <-  unique(data$Delta)
    CLev  <-  unique(data$C)
    nHs   <-  length(hLev)
    nDs   <-  length(dLev)
    nCs   <-  length(CLev)

    # Color scheme
    COLS  <-  list("line"     =  transparentColor('#252525', opacity=1),
                   "extinct"  =  transparentColor('red', opacity=0.15))

    # Set plot layout
    layout.mat  <- rbind(
                         matrix(c(1:(4*nCs)), nrow=4, ncol=nCs, byrow=TRUE),
                         rep((4*nCs+1), times=nCs),
                         matrix(c(((4*nCs+1)+(1:(4*nCs)))), nrow=4, ncol=nCs, byrow=TRUE)
                        )

    layout      <- layout(layout.mat,respect=TRUE)


        # Subset data by dominance coefficient
        i = 1
        d  <-  data[data$h == hLev[i],]            

        for(j in 1:nDs) {
            # Subset data by d_i
            dd  <-  d[d$Delta == dLev[j],]            

            # Subset data by Selfing Rate/d_i value
            for(k in 1:length(unique(dd$C))) {
                ddd  <-  dd[dd$C == CLev[k],]            

                # Initate plotting environment on first plot
                if(i == 1 && j == 1 && k == 1){
                    par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s',xpd=TRUE)
                }

                # Make dominance label 
                if(i == 1  & j == 1 & k == 13) {
                    proportionalLabel(0.5, 2, substitute(paste(italic(h), " = ", hh), list(hh = hLev[1])), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=0)
                }

                # Make plots
                plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2)
                usr  <-  par('usr')
                rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
                plotGrid(lineCol='grey80')
                box()
                # Extinction Thresholds
                lines(sf ~ smExt, lty=1, lwd=0.75, col=COLS$line, data=ddd)
                lines(sfExt ~ sm, lty=1, lwd=0.75, col=COLS$line, data=ddd)
                # Invasion Boundaries
                if(j == 1) {
                    lines(aInv[aInv < pars$sMax] ~ smInv[aInv < pars$sMax], lty=1, lwd=2, col=COLS$line, data=ddd)
                    lines(AInv ~ smInv, lty=1, lwd=2, col=COLS$line, data=ddd)
                }
                if(j > 1) {
                    lines(aInv[aInv <= pars$sMax] ~ smInv[aInv <= pars$sMax], lty=1, lwd=2, col=COLS$line, data=ddd)
                    lines(AInv ~ smInv, lty=1, lwd=2, col=COLS$line, data=ddd)
                }
                # axes
                axis(1, las=1, labels=NA)
                axis(2, las=1, labels=NA)
                # labels & annotations
                if(k == 1) {
                    axis(2, las=1)
                }
                if(j == 4) {
                    axis(1, las=1)
                }
                if(j == 1 && k == 1){
                    proportionalLabel(-0.3, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
                    proportionalLabel(-0.7, 0.5, expression(paste(italic(delta))), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=90)
                }
                if(j==2 & k==1) {
                    proportionalLabel(-0.7, 0.5, expression(paste(italic(delta[J]))), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=90)
                }
                if(j==3 & k==1) {
                    proportionalLabel(-0.7, 0.5, expression(paste(italic(delta[A]))), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=90)
                }
                if(j==4 & k==1) {
                    proportionalLabel(-0.7, 0.5, expression(paste(italic(delta[gamma]))), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=90)
                }

                if(i == 1 && j == 1) {
                    proportionalLabel(0.5, 1.15, substitute(paste(italic(C), " = ", CC), list(CC = CLev[k])), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
                }
                if(i == 2 && j == 4) {
                    axis(1, las=1)
                    proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
                }
            }
        }
    
        # Make empty plot for Dominance Label
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2)
        proportionalLabel(0.5, 0.5, substitute(paste(italic(h), " = ", hh), list(hh = hLev[2])), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=0)

        # Subset data by dominance coefficient
        i = 2
        d  <-  data[data$h == hLev[i],]            

        for(j in 1:nDs) {
            # Subset data by d_i
            dd  <-  d[d$Delta == dLev[j],]            

            # Subset data by Selfing Rate/d_i value
            for(k in 1:length(unique(dd$C))) {
                ddd  <-  dd[dd$C == CLev[k],]            

                # Initate plotting environment on first plot
                if(i == 1 && j == 1 && k == 1){
                    par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s',xpd=TRUE)
                }

                # Make dominance label 
                if(i == 1  & j == 1 & k == 12) {
                    proportionalLabel(0.5, 0.5, substitute(paste(italic(h), " = ", hh), list(hh = hLev[1])), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=0)
                }

                # Make plots
                plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.15), ylim = c(0,0.15), ylab='', xlab='', cex.lab=1.2)
                usr  <-  par('usr')
                rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
                plotGrid(lineCol='grey80')
                box()
                # Extinction Thresholds
                lines(sf ~ smExt, lty=1, lwd=0.75, col=COLS$line, data=ddd)
                lines(sfExt ~ sm, lty=1, lwd=0.75, col=COLS$line, data=ddd)
                # Invasion Boundaries
                if(j == 1) {
                    lines(aInv[aInv < pars$sMax] ~ smInv[aInv < pars$sMax], lty=1, lwd=2, col=COLS$line, data=ddd)
                    lines(AInv ~ smInv, lty=1, lwd=2, col=COLS$line, data=ddd)
                }
                if(j > 1) {
                    lines(aInv[aInv <= pars$sMax] ~ smInv[aInv <= pars$sMax], lty=1, lwd=2, col=COLS$line, data=ddd)
                    lines(AInv ~ smInv, lty=1, lwd=2, col=COLS$line, data=ddd)
                }
                # axes
                axis(1, las=1, labels=NA)
                axis(2, las=1, labels=NA)
                # labels & annotations
                if(k == 1) {
                    axis(2, las=1)
                }
                if(j == 4) {
                    axis(1, las=1)
                }
                if(j == 1 && k == 1){
                    proportionalLabel(-0.3, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
                    proportionalLabel(-0.7, 0.5, expression(paste(italic(delta))), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=90)
                }
                if(j==2 & k==1) {
                    proportionalLabel(-0.7, 0.5, expression(paste(italic(delta[J]))), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=90)
                }
                if(j==3 & k==1) {
                    proportionalLabel(-0.7, 0.5, expression(paste(italic(delta[A]))), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=90)
                }
                if(j==4 & k==1) {
                    proportionalLabel(-0.7, 0.5, expression(paste(italic(delta[gamma]))), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=90)
                }

                if(i == 1 && j == 1) {
                    proportionalLabel(0.5, 1.15, substitute(paste(italic(C), " = ", CC), list(CC = CLev[k])), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
                }
                if(i == 2 && j == 4) {
                    axis(1, las=1)
                    proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
                }
            }
        }
}




##############################################################
##############################################################
##  Preliminary Figs
##############################################################
##############################################################


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







