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

toPdf(FunnelPlots(df1 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.5_hm0.5_C0_delta0',
				  df2 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.5_hm0.5_C0.5_delta0',
				  df3 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.25_hm0.25_C0_delta0',
				  df4 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.25_hm0.25_C0.5_delta0'), 
figPath(name='FunnelPlots_MEqAInv.pdf'), width=10, height=7)
embed_fonts(figPath(name='FunnelPlots_MEqAInv.pdf'))

toPdf(FunnelPlots(df1 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.5_hm0.5_C0_delta0',
				  df2 = 'demSimsSfxSm_MinorEqAlleleInv2_sMax0.15_nSamples1000_hf0.5_hm0.5_C0.5_delta0.125',
				  df3 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.25_hm0.25_C0_delta0',
				  df4 = 'demSimsSfxSm_MinorEqAlleleInv2_sMax0.15_nSamples1000_hf0.25_hm0.25_C0.5_delta0.125'), 
figPath(name='FunnelPlots_MEqAInv_ID.pdf'), width=10, height=7)
embed_fonts(figPath(name='FunnelPlots_MEqAInv.pdf'))






