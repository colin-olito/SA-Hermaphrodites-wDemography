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

toPdf(FunnelEigSimCompare(df1 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.5_hm0.5_C0_delta0',
						  df2 = 'demSimsSfxSm_MinorEqAlleleInv2_sMax0.15_nSamples1000_hf0.5_hm0.5_C0.5_delta0',
						  df3 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.25_hm0.25_C0_delta0',
						  df4 = 'demSimsSfxSm_MinorEqAlleleInv2_sMax0.15_nSamples1000_hf0.25_hm0.25_C0.5_delta0'), 
			figPath(name='Funnel_Compare_Eig_Sim.pdf'), width=10, height=7)
embed_fonts(figPath(name='Funnel_Compare_Eig_Sim.pdf'))









toPdf(polySpaceFig(df1 = "simPolySpace_sMax0.15_nSamples5000_hf0.5_hm0.5_delta0_dj0_da0_dg0_f5.8",
				   df2 = "simPolySpace_sMax0.15_nSamples5000_hf0.5_hm0.5_delta0_dj0_da0_dg0_f6",
				   df3 = "simPolySpace_sMax0.15_nSamples5000_hf0.5_hm0.5_delta0_dj0_da0_dg0_f6.5",
				   df4 = "simPolySpace_sMax0.15_nSamples1000_hf0.25_hm0.25_delta0_dj0_da0_dg0_f5.8",
				   df5 = "simPolySpace_sMax0.15_nSamples1000_hf0.25_hm0.25_delta0_dj0_da0_dg0_f6",
				   df6 = "simPolySpace_sMax0.15_nSamples1000_hf0.25_hm0.25_delta0_dj0_da0_dg0_f6.5"), 
			figPath(name='polymorphicSpace.pdf'), width=5, height=7)
embed_fonts(figPath(name='polymorphicSpace.pdf'))
