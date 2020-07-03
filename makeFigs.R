#####################################################
#'  SA in hermaphrodites with demography 
#' and partial selfing
#'
#'  Functions to generate Figures for: 
#'    Article title goes here...
#'  
#'  Author: Colin Olito
#'
#'
#'  NOTES: Run this file, either from terminal using Rscript,
#'		  or interactively in R. This should create all the 
#'		  figures needed to correctly compile the mansucript
#'		  LaTeX file.  
#'          

rm(list=ls())
################
## Dependencies
source('R/functions-Figs.R')

###############
# PRELIM FIGS
###############

toPdf(FunnelPlots(), figPath(name='FunnelPlots.pdf'), width=10, height=7)
embed_fonts(figPath(name='FunnelPlots.pdf'))