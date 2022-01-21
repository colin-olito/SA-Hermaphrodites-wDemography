#####################################################
#'  R code to generate .pdf figures for:
#' 
#'  Olito, C. and C. DeVries. 2022. Demographic consequences of 
#'  sexually antagonistic selection in partially selfing
#' populations. American Naturalist XX: XX--XX. doi: XXX
#'
#'  Author: Colin Olito
#' 
#'
#'  NOTES: Run this file, either from terminal using Rscript,
#'		   or interactively in R AFTER running the necessary
#'		   simulations provided in run-Simulations.R. This 
#' 		   should create .pdf files for each of the figures
#' 		   needed to correctly compile the manuscript LaTeX
#'		   file.  
#'         
#' 			To reproduce the figures for the paper and 
#' 			supplementary figures, run the following R code
#' 			up to L.108. This will generate output data as
#' 			.csv files saved to ./output/figs/
#' 			Be sure to create this directory locally prior
#' 			to running the simulations.

rm(list=ls())
################
## Dependencies
source('R/functions-Figs.R')
#source('./R/loadData-Compadre.R')



#########################################################
#' Functions to generate figures for the published paper
#########################################################

#' Fig. 1 - Didactic figure illustrating genetic and demographic predictions
#'			made by the model
toPdf(resultsIllustrationFig(), 
			figPath(name='didacticFig.pdf'), width=8, height=4)
embed_fonts(figPath(name='didacticFig.pdf'))


#' Fig. 2 - Illustration of invasion conditions & extinction thresholds
toPdf(demViablePolySpaceFig(), 
			figPath(name='extinctionThresholdsFig.pdf'), width=10, height=5)
embed_fonts(figPath(name='extinctionThresholdsFig.pdf'))

#' Fig. S1 - Illustration of population intrinsic growth rates (lambda)
#'				   for the same parameter sets used in Fig.1 
toPdf(lambdaPolySpaceFig(), 
			figPath(name='lambdaFig1_2.pdf'), width=10, height=16)
embed_fonts(figPath(name='lambdaFig1_2.pdf'))


#' Fig. 3 - Quantification of demographically viable parameter space, 
#' 			no inbreeding depression
toPdf(polySpaceFigTitrate(), 
			figPath(name='polymorphicSpaceTitrate.pdf'), width=5, height=7)
embed_fonts(figPath(name='polymorphicSpaceTitrate.pdf'))

#' Fig. 4 - Quantification of demographically viable parameter space, 
#'			showing effects of inbreeding depression
toPdf(deltaSelfingLoadPolySpaceFigTitrate(), 
			figPath(name='deltaPolymorphicSpaceTitrate.pdf'), width=5, height=7)
embed_fonts(figPath(name='deltaPolymorphicSpaceTitrate.pdf'))

#' Deprecated Fig. 5 - Mimulus Case study: polymorphic parameter space w/ 
#'			Inv6 plotted (using extinction thresholds)
#' toPdf(MimulusInv6Fig(), 
#'			figPath(name='MimulusInv6FigNewSel.pdf'), width=10, height=5)
#' embed_fonts(figPath(name='MimulusInv6FigNewSel.pdf'))


#' Fig. 5 - Mimulus Case study: lambda heatmap w/ Inv6
toPdf(MimulusInv6LambdaFig(), 
			figPath(name='MimulusInv6LambdaFig.pdf'), width=8, height=8)
embed_fonts(figPath(name='MimulusInv6LambdaFig.pdf'))


#' Fig. B1 - Density Dependent Dynamics
#' Mimulus Case study: Nstar heatmap w/ Inv6
toPdf(MimulusInv6NstarFig(), 
			figPath(name='MimulusInv6NstarFig.pdf'), width=8, height=8)
embed_fonts(figPath(name='MimulusInv6NstarFig.pdf'))



########################
# Supplementary Figures
########################

# Fig. S1 -- Illustration of Invaions & Extinction thresholds for sex-specific dominance
toPdf(demViablePolySpace_SexSpecDominance_Fig(), 
			figPath(name='extinctionThresholdsSexSpecDomFig.pdf'), width=10, height=5)
embed_fonts(figPath(name='extinctionThresholdsSexSpecDomFig.pdf'))

# Fig. S2 -- Quantification of demographically viable parameter space, no inbreeding depression

toPdf(polySpaceFigTitrateSexSpec(), 
			figPath(name='polymorphicSpaceTitrateSexSpec.pdf'), width=5, height=7)
embed_fonts(figPath(name='polymorphicSpaceTitrateSexSpec.pdf'))

# Fig. S3 -- Quantification of demographically viable parameter space, 
#			showing effects of inbreeding depression
toPdf(deltaSelfingLoadPolySpaceFigTitrateSexSpec(), 
			figPath(name='deltaPolymorphicSpaceTitrateCompSexSpec.pdf'), width=5, height=7)
embed_fonts(figPath(name='deltaPolymorphicSpaceTitrateCompSexSpec.pdf'))



####################################################################
# Preliminary Figures

#' The remainder of this file contains code to generate exploratory
#' figures & analyses written during the development of the study. 
#' Some  functions are buggy, and some are deprecated, having been 
#' replaced by better/more functional code for the final analyses.
####################################################################

#'  Funnel plots comparing population genetic model predictions
#'  with simulation results, distinguishing between different 
#'  outcomes (i.e., A fixes, a fixes, polymorphism, extinction)
toPdf(FunnelPlots(df1 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.5_hm0.5_C0_delta0',
				  df2 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.5_hm0.5_C0.5_delta0',
				  df3 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.25_hm0.25_C0_delta0',
				  df4 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.25_hm0.25_C0.5_delta0'), 
			figPath(name='FunnelPlots_MEqAInv.pdf'), width=10, height=7)
embed_fonts(figPath(name='FunnelPlots_MEqAInv.pdf'))


#'  Same as above, but using slightly different method for invading allele
toPdf(FunnelPlots(df1 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.5_hm0.5_C0_delta0',
				  df2 = 'demSimsSfxSm_MinorEqAlleleInv2_sMax0.15_nSamples1000_hf0.5_hm0.5_C0.5_delta0.125',
				  df3 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.25_hm0.25_C0_delta0',
				  df4 = 'demSimsSfxSm_MinorEqAlleleInv2_sMax0.15_nSamples1000_hf0.25_hm0.25_C0.5_delta0.125'), 
			figPath(name='FunnelPlots_MEqAInv_ID.pdf'), width=10, height=7)
embed_fonts(figPath(name='FunnelPlots_MEqAInv.pdf'))


#'  Funnel plots comparing simulation outcomes with the 
#'  analytic invasion criteria for the full matrix model
toPdf(FunnelEigSimCompare(df1 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.5_hm0.5_C0_delta0',
						  df2 = 'demSimsSfxSm_MinorEqAlleleInv2_sMax0.15_nSamples1000_hf0.5_hm0.5_C0.5_delta0',
						  df3 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.25_hm0.25_C0_delta0',
						  df4 = 'demSimsSfxSm_MinorEqAlleleInv2_sMax0.15_nSamples1000_hf0.25_hm0.25_C0.5_delta0'), 
			figPath(name='Funnel_Compare_Eig_Sim.pdf'), width=10, height=7)
embed_fonts(figPath(name='Funnel_Compare_Eig_Sim.pdf'))



#'  Polymorphic parameter space as a function of the selfing rate
#'  assuming no inbreeding depression
toPdf(polySpaceFig(df1 = "simPolySpace_sMax0.15_nSamples5000_hf0.5_hm0.5_delta0_dj0_da0_dg0_f5.8",
				   df2 = "simPolySpace_sMax0.15_nSamples5000_hf0.5_hm0.5_delta0_dj0_da0_dg0_f6",
				   df3 = "simPolySpace_sMax0.15_nSamples5000_hf0.5_hm0.5_delta0_dj0_da0_dg0_f6.5",
				   df4 = "simPolySpace_sMax0.15_nSamples1000_hf0.25_hm0.25_delta0_dj0_da0_dg0_f5.8",
				   df5 = "simPolySpace_sMax0.15_nSamples1000_hf0.25_hm0.25_delta0_dj0_da0_dg0_f6",
				   df6 = "simPolySpace_sMax0.15_nSamples1000_hf0.25_hm0.25_delta0_dj0_da0_dg0_f6.5"), 
			figPath(name='polymorphicSpace.pdf'), width=5, height=7)
embed_fonts(figPath(name='polymorphicSpace.pdf'))

toPdf(polySpaceFig(df1 = "simPolySpaceNew_sMax0.15_nSamples2000_hf0.5_hm0.5_delta0_dj0_da0_dg0_f5.8",
				   df2 = "simPolySpaceNew_sMax0.15_nSamples2000_hf0.5_hm0.5_delta0_dj0_da0_dg0_f6",
				   df3 = "simPolySpaceNew_sMax0.15_nSamples2000_hf0.5_hm0.5_delta0_dj0_da0_dg0_f6.5",
				   df4 = "simPolySpaceNew_sMax0.15_nSamples1000_hf0.25_hm0.25_delta0_dj0_da0_dg0_f5.8",
				   df5 = "simPolySpaceNew_sMax0.15_nSamples1000_hf0.25_hm0.25_delta0_dj0_da0_dg0_f5.9",
				   df6 = "simPolySpaceNew_sMax0.15_nSamples1000_hf0.25_hm0.25_delta0_dj0_da0_dg0_f6.5"), 
			figPath(name='polymorphicSpaceNew.pdf'), width=5, height=7)
embed_fonts(figPath(name='polymorphicSpaceNew.pdf'))


#'  Polymorphic parameter space as a function of inbreeding depression
#'  for different fixed selfing rates
toPdf(deltaEffectsPolySpaceFig(
			df1 = "deltaSimPolySpace_sMax0.15_nSamples1000_C0.25_hf0.5_hm0.5_f6.5",
			df2 = "deltaSimPolySpace_sMax0.15_nSamples1000_C0.5_hf0.5_hm0.5_f6.5",
			df3 = "deltaSimPolySpace_sMax0.15_nSamples1000_C0.75_hf0.5_hm0.5_f6.5",
			df4 = "deltaSimPolySpace_sMax0.15_nSamples1000_C0.25_hf0.25_hm0.25_f6.5",
			df5 = "deltaSimPolySpace_sMax0.15_nSamples1000_C0.5_hf0.25_hm0.25_f6.5",
			df6 = "deltaSimPolySpace_sMax0.15_nSamples1000_C0.75_hf0.25_hm0.25_f6.5"),
			figPath(name='deltaPolySpace.pdf'), width=7, height=5)
embed_fonts(figPath(name='deltaPolySpace.pdf'))


#'  Polymorphic parameter space as a function of the selfing rate
#'  using heuristic model of inbreeding depression ~ selfing
toPdf(deltaSelfingLoadPolySpaceFig(
			df1 = "deltaSelfingSimPolySpace_sMax0.15_nSamples1000_dStar0.8_hf0.5_hm0.5_f7",
			df2 = "deltaSelfingSimPolySpace_sMax0.15_nSamples1000_dStar0.8_hf0.25_hm0.25_f7"),
			figPath(name='deltaSelfingPolySpace.pdf'), width=5, height=7)
embed_fonts(figPath(name='deltaSelfingPolySpace.pdf'))





#'  Polymorphic parameter space as a function of the selfing rate
#'  using heuristic model of inbreeding depression ~ selfing
toPdf(MimulusInv6Fig(),
			figPath(name='mimulusInv6Fig.pdf'), width=7, height=5)
embed_fonts(figPath(name='mimulusInv6Fig.pdf'))



#' Ugly multi-panel figure showing invasion conditions and extinction
#' thresholds for the same parameter conditions that are summarized
#' in figures 3 & 4.
toPdf(suppPolySpaceThresholdFigs(), 
			figPath(name='SuppFig-polySpace.pdf'), width=20, height=8)
embed_fonts(figPath(name='SuppFig-polySpace.pdf'))


toPdf(suppDeltaPolySpaceThresholdFigs(), 
			figPath(name='SuppFig-deltaPolySpace.pdf'), width=40, height=20)
embed_fonts(figPath(name='SuppFig-deltaPolySpace.pdf'))

