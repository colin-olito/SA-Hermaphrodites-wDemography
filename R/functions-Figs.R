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
        proportionalLabel(-0.3, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
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
        proportionalLabel(-0.3, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.04, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)


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
        proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)

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
        proportionalLabel(0.5, -0.3, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
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
        proportionalLabel(-0.3, 0.5, expression(paste("Prop. viable polymorphic space")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=90)

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
        proportionalLabel(-0.3, 0.5, expression(paste("Prop. viable polymorphic space")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=90)

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

deltaSelfingLoadPolySpaceFigTitrate  <-  function(df = "dataDeltaPolySpaceFig_sMax0.15_res0.003_dStar0.8_f6.5") {

    # Make filenames for import from df names
    fName  <-  paste('./output/simData/', df, '.csv', sep="")

    # import data
    data  <-  read.csv(file=fName, header=TRUE)

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

    # clean data set & quantify parameter space
    dat   <-  quantDeltaPolySpace(data = data, pars = pars)


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
    dStar  <-  pars$dStar
    CLine  <-  seq(0,0.9,length=100)
    dLine  <-  predDelta(dStar=dStar, b=1/2, a=0.2, C=CLine) 


# Set plot layout
    layout.mat  <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
    layout      <- layout(layout.mat,respect=TRUE)

    ## Panel A: Additive SA (hf = hm = 1/2)
    ##          early-acting delta
        PGSpace     <-  c()
        for(i in 1:length(CLine)) {
            PGSpace[i]     <-  popGen_PolySpace_Delta_Add(C=CLine[i], delta=dLine[i], sMax=pars$sMax)
        }
        d  <-  dat[dat$h == 0.5,]
        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,4,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.925), ylim = c(0,0.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(PGSpace ~ CLine, lwd=2, col=COLS$PG)
        points(PrViaPoly[Delta == "d"][c(1:9,21:23)] ~ C[Delta == "d"][c(1:9,21:23)], pch=21, bg=COLS$dSim, col=COLS$dSim2, data=d)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.2, expression(paste("Early-acting Inbreeding Depression")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(-0.5, 0.5, expression(paste(italic(h), " = ", 1/2)), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.3, 0.5, expression(paste("Prop. viable polymorphic space")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)

      #Legend
        legend( x       =  usr[2]*0.95,
                y       =  usr[4],
                legend  =  c(
                             expression(paste(italic(delta), " (Pop. Gen.)")),
                             expression(paste(delta))),
                 lty     =  c(1,NA),
                 lwd     =  c(2,NA),
                 col     =  c(COLS$PG,
                              COLS$dSim),
                 pch     =  c(NA,21),
                 pt.bg   =  c(NA,
                              COLS$dSim),
                 cex     =  0.75,
                 pt.cex  =  0.75,
                 xjust   =  1,
                 yjust   =  1,
                 bty     =  'n',
                 border  =  NA)
    
    ## Panel B: Additive SA (hf = hm = 1/2)
    ##          Late-acting delta
        plot(NA, axes=FALSE, type='n', main='', xlim = c(0,0.925), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(PGSpace ~ CLine, lwd=2, col=COLS$PG)
        points(PrViaPoly[Delta == "d_j"][c(1:8,21:23)] ~ C[Delta == "d_j"][c(1:8,21:23)], pch=21, bg=COLS$d_j, col=COLS$d_j2, data=d)
        points(PrViaPoly[Delta == "d_a"][c(1:12,21:23)] ~ C[Delta == "d_a"][c(1:12,21:23)], pch=21, bg=COLS$d_a, col=COLS$d_a2, data=d)
        points(PrViaPoly[Delta == "d_g"][c(1:10,21:23)] ~ C[Delta == "d_g"][c(1:10,21:23)], pch=21, bg=COLS$d_g, col=COLS$d_g2, data=d)
        # axes        
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.03, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.2, expression(paste("Late-acting Inbreeding Depression")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        #Legend
        legend( x       =  usr[2]*0.95,
                y       =  usr[4],
                legend  =  c(
                             expression(paste(italic(delta), " (Pop. Gen.)")),
                             expression(paste(delta[italic(j)])),
                             expression(paste(delta[italic(a)])),
                             expression(paste(delta[gamma]))),
                 lty     =  c(1,NA,NA,NA),
                 lwd     =  c(2,NA,NA,NA),
                 col     =  c(COLS$PG,
                              COLS$d_j,
                              COLS$d_a,
                              COLS$d_g),
                 pch     =  c(NA,21,21,21),
                 pt.bg   =  c(NA,
                              COLS$d_j,
                              COLS$d_a,
                              COLS$d_g),
                 cex     =  0.75,
                 pt.cex  =  0.75,
                 xjust   =  1,
                 yjust   =  1,
                 bty     =  'n',
                 border  =  NA)


    ## Panel C: Dominance Reversal SA (hf = hm = 1/4)
    ##          early-acting delta
        CLine2  <-  CLine[-1]
        dLine2  <-  dLine[-1]
        PGSpace     <-  c()
        for(i in 1:length(CLine2)) {
            PGSpace[i]     <-  popGen_PolySpace_Delta_DomRev(C=CLine2[i], delta=dLine2[i], sMax=pars$sMax)
        }
        d  <-  dat[dat$h == 0.25,]
        # Make the plot
        plot(NA, axes=FALSE, type='n', main='', xlim = c(0,0.925), ylim = c(0,0.8), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(PGSpace ~ CLine2, lwd=2, col=COLS$PG)
        points(PrViaPoly[Delta == "d"][c(1:9,21:23)] ~ C[Delta == "d"][c(1:9,21:23)], pch=21, bg=COLS$dSim, col=COLS$dSim2, data=d)
        # axes        
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(0.03, 1.075, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.5, 0.5, expression(paste(italic(h), " = ", 1/4)), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.3, 0.5, expression(paste("Prop. viable polymorphic space")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.3, expression(paste("Selfing Rate (", italic(C), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)

    ## Panel D: Dominance Reversal SA (hf = hm = 1/4)
    ##          Late-acting delta
        plot(NA, axes=FALSE, type='n', main='', xlim = c(0,0.925), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation Results
        lines(PGSpace ~ CLine2, lwd=2, col=COLS$PG)
        points(PrViaPoly[Delta == "d_j"][c(1:8,21:23)] ~ C[Delta == "d_j"][c(1:8,21:23)], pch=21, bg=COLS$d_j, col=COLS$d_j2, data=d)
        points(PrViaPoly[Delta == "d_a"][c(1:12,21:23)] ~ C[Delta == "d_a"][c(1:12,21:23)], pch=21, bg=COLS$d_a, col=COLS$d_a2, data=d)
        points(PrViaPoly[Delta == "d_g"][c(1:10,21:23)] ~ C[Delta == "d_g"][c(1:10,21:23)], pch=21, bg=COLS$d_g, col=COLS$d_g2, data=d)
        # axes        
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(0.03, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.3, expression(paste("Selfing Rate (", italic(C), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)

}


##############################################################
##############################################################
##  Figures for Supplementary Material
##############################################################
##############################################################

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
        #alternative selfing rate (C=0.38)
    inv_A  <-  popGen_A_invade(hf=hf, hm=hm, sm=sms, C=0.24)
    inv_a  <-  popGen_a_invade(hf=hf, hm=hm, sm=sms, C=0.24)
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


