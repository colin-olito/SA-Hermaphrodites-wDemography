################################################################
#' Load measurements for estimating selection coefficients for
#' inv6 from figures in Lee et al. (2016), back-calculate
#' selection coefficients, and save as objects to be called 
#' by ./makeFigs.R
#'	
#'  Author: Colin Olito
#' 
#' 	Notes:	Estimates for male and female fitness components
#' were back-calculated from direct measures of figures 2, 5,
#' & 6 in:
#' 
#' Lee, Y. W., L. Fishman, J. K. Kelly, and J. H. Willis. 2016.
#' A segregating inversion generates fitness variation in 
#' yellow monkeyflower (Mimulus guttatus). Genetics 202:1473–1484.
#' 
#' Raw measurements of the figures are provided as supplementary
#' files  inv6-measurements.pptx and inv6-measurements.odp, 
#' which are available on DRYAD at:
#' 
#' DRYAD Reviewer link (stable link to be provided when assigned)
#' https://datadryad.org/stash/share/81sAuXGEg8cSh-S9VVL0PfBCsl6YLkG1OIFBCvOefac 





###############
# Import data
###############

inv6Data  <-  read.csv('./data/inv6_Figure_Measurements.csv', header=TRUE)

#############################################
# Rescale back to original measurement units
#############################################

# Female fitness components
no_inv6_flowers   <-  (inv6Data$no_inv6_flowers_inches / inv6Data$pptInches[inv6Data$scaleUnits == 'flowers']) * inv6Data$scaleAxis[inv6Data$scaleUnits == 'flowers'] 
inv6_het_flowers  <-  (inv6Data$inv6_het_flowers_inches / inv6Data$pptInches[inv6Data$scaleUnits == 'flowers']) * inv6Data$scaleAxis[inv6Data$scaleUnits == 'flowers'] 

no_inv6_fruits   <-  (inv6Data$no_inv6_fruits_inches / inv6Data$pptInches[inv6Data$scaleUnits == 'fruits']) * inv6Data$scaleAxis[inv6Data$scaleUnits == 'fruits'] 
inv6_het_fruits  <-  (inv6Data$inv6_het_fruits_inches / inv6Data$pptInches[inv6Data$scaleUnits == 'fruits']) * inv6Data$scaleAxis[inv6Data$scaleUnits == 'fruits'] 

no_inv6_seeds   <-  (inv6Data$no_inv6_seeds_inches / inv6Data$pptInches[inv6Data$scaleUnits == 'seeds']) * inv6Data$scaleAxis[inv6Data$scaleUnits == 'seeds'] 
inv6_het_seeds  <-  (inv6Data$inv6_het_seeds_inches / inv6Data$pptInches[inv6Data$scaleUnits == 'seeds']) * inv6Data$scaleAxis[inv6Data$scaleUnits == 'seeds'] 

# Pollen Viability
no_inv6_PrViablePollen   <-  (inv6Data$no_inv6_PrViablePollen_inches / inv6Data$pptInches[inv6Data$scaleUnits == 'PrViablePollen']) * inv6Data$scaleAxis[inv6Data$scaleUnits == 'PrViablePollen'] 
inv6_het_PrViablePollen  <-  (inv6Data$inv6_het_PrViablePollen_inches / inv6Data$pptInches[inv6Data$scaleUnits == 'PrViablePollen']) * inv6Data$scaleAxis[inv6Data$scaleUnits == 'PrViablePollen'] 

#############################################
# Calculate inv6_homozygote values, assuming
# additive fitness effects

inv6_hom_flowers         <-  inv6_het_flowers + (inv6_het_flowers - no_inv6_flowers)
inv6_hom_fruits          <-  inv6_het_fruits + (inv6_het_fruits - no_inv6_fruits)
inv6_hom_seeds           <-  inv6_het_seeds + (inv6_het_seeds - no_inv6_seeds)
inv6_hom_PrViablePollen  <-  inv6_het_PrViablePollen - (no_inv6_PrViablePollen - inv6_het_PrViablePollen)

#############################################
# Calculate geometric mean flower number 
# for (2012 & 2013)

gM_no_inv6_flowers   <-  exp(mean(log(c(no_inv6_flowers[inv6Data$year == 2012], no_inv6_flowers[inv6Data$year == 2013]))))
gM_inv6_het_flowers  <-  exp(mean(log(c(inv6_het_flowers[inv6Data$year == 2012], inv6_het_flowers[inv6Data$year == 2013]))))
gM_inv6_hom_flowers  <-  gM_inv6_het_flowers + (gM_inv6_het_flowers - gM_no_inv6_flowers)

#############################################
# Calculate metric of pollen production 
# (pollen viability * geom. mean flower number)

no_inv6_pollenProd   <-  no_inv6_PrViablePollen * gM_no_inv6_flowers
inv6_het_pollenProd  <-  inv6_het_PrViablePollen * gM_inv6_het_flowers
inv6_hom_pollenProd  <-  inv6_het_pollenProd - (no_inv6_pollenProd - inv6_het_pollenProd)

#############################################
# Calculate relative fitnesses

no_inv6_wf   <-  gM_no_inv6_flowers / gM_inv6_hom_flowers
inv6_het_wf  <-  gM_inv6_het_flowers / gM_inv6_hom_flowers
inv6_hom_wf  <-  gM_inv6_hom_flowers / gM_inv6_hom_flowers

no_inv6_wm   <-  (no_inv6_pollenProd / no_inv6_pollenProd)[1]
inv6_het_wm  <-  (inv6_het_pollenProd / no_inv6_pollenProd)[1]
inv6_hom_wm  <-  (inv6_hom_pollenProd / no_inv6_pollenProd)[1]


#############################################
# Calculate s_f & s_m

inv6_sf  <-  round(inv6_hom_wf - no_inv6_wf, digits=2)
inv6_sm  <-  round(no_inv6_wm - inv6_hom_wm, digits=2)
