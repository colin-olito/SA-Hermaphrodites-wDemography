################################################################
#'  R code to reproduce main results for:
#' 
#'  Olito, C. and C. DeVries. 2022. Demographic consequences of 
#'  sexually antagonistic selection in partially selfing
#' populations. American Naturalist XX: XX--XX. doi: XXX
#'
#'  Author: Colin Olito
#'
#'  NOTES:  To reproduce the analyses for the paper, run the 
#' 			following R code up to L.728. This will generate
#' 			output data as .csv files saved to ./output/simData/
#' 			Be sure to create this directory locally prior
#' 			to running the simulations.
#'		


#####################
##  Dependencies
rm(list=ls())
source('R/functions-Simulations.R')
source('R/loadData-Compadre.R')

###########################################
##  Import parameter values from tables in
##  Peterson et al. (2016) & Willis (1993)
petersonData  <-  read.csv('./data/Peterson_2016_Data.csv', header=TRUE)
WillisData    <-  read.csv('./data/Willis_1993_Data.csv', header=TRUE)
Peterson2016.EM.theta.list  <-  list(D = petersonData$D[petersonData$population == "EM"], 
						G = petersonData$G[petersonData$population == "EM"],
						F = petersonData$F[petersonData$population == "EM"],
						O = petersonData$O[petersonData$population == "EM"],
						A = petersonData$A[petersonData$population == "EM"],
						S = petersonData$S[petersonData$population == "EM"],
						R = petersonData$R[petersonData$population == "EM"],
						pop = petersonData$population[petersonData$population == "EM"]
						)
Peterson2016.LEP.theta.list  <-  list(D = petersonData$D[petersonData$population == "LEP"], 
						 G = petersonData$G[petersonData$population == "LEP"],
						 F = petersonData$F[petersonData$population == "LEP"],
						 O = petersonData$O[petersonData$population == "LEP"],
						 A = petersonData$A[petersonData$population == "LEP"],
						 S = petersonData$S[petersonData$population == "LEP"],
						 R = petersonData$R[petersonData$population == "LEP"],
						 pop = petersonData$population[petersonData$population == "LEP"]
						)
Willis1993.SelfingRate  <-  WillisData$C
Willis1993.delta.list   <-  list(delta_D = WillisData$delta_D,
								 delta_G = WillisData$delta_G,
								 delta_F = WillisData$delta_F,
								 delta_O = WillisData$delta_O,
								 delta_S = WillisData$delta_S
							 )

######################################
#' Simulations for the published paper
######################################

######################################
#' Create output directories if they
#' do not already exist
outputDirectoriesExist  <-  dir.exists("./output")

if(!outputDirectoriesExist) {
	dir.create("./output")
	dir.create("./output/simData")
	dir.create("./output/figs")
}

#############################################################
#' Figs. 1 & 2 - Invasion conditions & extinction thresholds

## Additive SA Fitness
## hf = hm = 1/2

	# Obligate Outcrossing
	# C = 0
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  0,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

titrateInvBoundaries(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^3, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6.2), theta_prime = 6.2, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)



	# Intermediate Selfing
	# C = 1/4
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  1/4,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

titrateInvBoundaries(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^3, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6.2), theta_prime = 6.2, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

	# Intermediate Selfing
	# C = 1/2
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  1/2,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

titrateInvBoundaries(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^3, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6.2), theta_prime = 6.2, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

## Partially Recessive SA Fitness
## hf = hm = 1/4

	# Obligate Outcrossing
	# C = 0
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  0,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

titrateInvBoundaries(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^3, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.9), theta_prime = 5.9, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6.0), theta_prime = 6.0, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

	# Intermediate Selfing
	# C = 1/4
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  1/4,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

titrateInvBoundaries(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^3, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
						 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
						 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
						 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
						 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
						 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.9), theta_prime = 5.9, 
						 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
						 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
						 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
						 om = 2, g = 3, theta = c(0.6,0.6,0.05,6.0), theta_prime = 6.0, 
						 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
						 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
						 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)



	# Intermediate Selfing
	# C = 1/2
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  1/2,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

titrateInvBoundaries(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^3, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.9), theta_prime = 5.9, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6.0), theta_prime = 6.0, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)



##########################################################
#' Figs. 1 & S1 - Make data to generate heatmaps of lambda 
#'				  for same parameter conditions presented
#' 				  in Figs. 1 & S1. 


## Additive SA Fitness
## hf = hm = 1/2

	# Obligate Outcrossing
	# C = 0
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  0,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)
# Export list for dopar function
allObj  <-  ls()
defObj  <-  c("i")
funs    <-  allObj[-which(allObj %in% defObj)]

# Note: you can set nCluster to what ever value you want (provided it makes sense
#       given how many cores your computer has available. If left as nCluster = 'NA')
#       the default will be to use 2*(detectCores() - 1)
makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4, alpha=0,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 4*10^5, eqThreshold = 1e-9,
					 nCluster = 10, funs = funs, writeFile = TRUE)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4, alpha=0,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 4*10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4, alpha=0,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6.2), theta_prime = 6.2, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 4*10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)




	# Intermediate Selfing
	# C = 1/4
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  1/4,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6.2), theta_prime = 6.2, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)


	# Intermediate Selfing
	# C = 1/2
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  1/2,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6.2), theta_prime = 6.2, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)


## Partially Recessive SA Fitness
## hf = hm = 1/4

	# Obligate Outcrossing
	# C = 0
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  0,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.9), theta_prime = 5.9, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)



	# Intermediate Selfing
	# C = 1/4
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  1/4,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.9), theta_prime = 5.9, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)





	# Intermediate Selfing
	# C = 1/2
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  1/2,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.9), theta_prime = 5.9, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)

makeLambdaHeatMapData(sMax=0.15, len=100, precision = 1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8,
					 nCluster = 12, funs=funs, writeFile = TRUE)



##########################################################
#' Fig. 3 - Run simulations to calculate demographically 
#'			viable polymorphic parameter space (no I.D.)

makeDataPolyParamSpace(sMax=0.15, res=0.003, precision = 1e-4,
						om = 2, g = 3, theta = c(0.6,0.6,0.05,NA), theta_prime = NA, 
						hfVals= c(1/2, 1/4), hmVals= c(1/2, 1/4), fVals = c(5.8, 5.9, 6.0, 6.5), 
						delta = 0, delta_j = 0, delta_a = 0, delta_gamma = 0,
						tlimit = 10^5, eqThreshold = 1e-8)



##########################################################
#' Fig. 4 - Run simulations to calculate demographically 
#'			viable polymorphic parameter space (w/ I.D.)

makeDataDeltaPolyParamSpace(sMax=0.15, res=0.003, precision = 1e-4,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = 6.5, 
							hfVals= c(1/2, 1/4), hmVals= c(1/2, 1/4), dStar = 0.8, 
							tlimit = 10^5, eqThreshold = 1e-8)

makeDataDeltaPolyParamSpace(sMax=0.15, res=0.003, precision = 1e-4,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,7), theta_prime = 7, 
							hfVals= c(1/2, 1/4), hmVals= c(1/2, 1/4), dStar = 0.8, 
							tlimit = 10^5, eqThreshold = 1e-8)

makeDataDeltaPolyParamSpace(sMax=0.15, res=0.003, precision = 1e-4,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,7.5), theta_prime = 7.5, 
							hfVals= c(1/2, 1/4), hmVals= c(1/2, 1/4), dStar = 0.8, 
							tlimit = 10^5, eqThreshold = 1e-8)

makeDataDeltaPolyParamSpace(sMax=0.15, res=0.003, precision = 1e-4,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,8.5), theta_prime = 8.5, 
							hfVals= c(1/2, 1/4), hmVals= c(1/2, 1/4), dStar = 0.8, 
							tlimit = 10^5, eqThreshold = 1e-8)
 



#############################################################
#' Fig. 5 - Comparison of Demographically viable polymorphic 
#'			parameter space using locally and non-locally 
#' 			adapted population demographic parameters 
#' 			for C =0 and C = empirical est. (w/ I.D.).
#' 			
#' 			Makes data to generate heatmaps of lambda
#' 			taking the average (2012 & 2013) flower #
#' 			effect into account when calculating male
#' 			selection coeficients and dominance for
#' 			inv6.
#' 
#' 			The empirical estimates for the demographic
#' 			rates used in this analysis are provided in
#' 			objects theta.list and delta.list below.
#' 			They are taken directly from Table S1 from
#' 			Peterson et al. (2016) (see references)
#' 			
#' 			Empirical estimates for the selfing rates
#' 			for Iron Mountain are taken directly from
#' 			Table 1 of Willis (1993). 
#' 			Inbreeding depression parameters are 
#' 			calculated directly from Table 2 of 
#' 			Willis (1993). 

# Export list for dopar function
allObj  <-  ls()
defObj  <-  c("i")
funs    <-  allObj[-which(allObj %in% defObj)]

# Eagle Meadows population (locally adapted)
datMat      <-  matList_Mg_EM # created in ./R/loadData-Compadre.R
theta.list  <-  Peterson2016.EM.theta.list

# inbreeding depression parameters set to 0 due to obligate outcrossing (C = 0)
delta.list  <-  list(delta_D = 0,
					 delta_G = 0,
					 delta_F = 0,
					 delta_O = 0,
					 delta_S = 0
					 )

## Obligate oucrossing
makeLambdaHeatMapMimulusData(sMax = 0.99, len = 100,
							 datMat = datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
							 hf = 1/2, hm = 1/2, C = 0,
							 tlimit = 4*10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
							 nCluster = 12, funs = funs, writeFile = TRUE)

## Empirical selfing rate estimates
## Empirical Inbreeding Depression (delta values calculated from Willis 1993)
delta.list  <-  Willis1993.delta.list

# Empirical selfing rate
makeLambdaHeatMapMimulusData(sMax = 0.99, len = 100,
							 datMat = datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
							 hf=1/2, hm = 1/2, C = Willis1993.SelfingRate, 
							 tlimit = 4*10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
							 nCluster = 12, funs = funs, writeFile = TRUE)


##########################
# Low-elevation perenials (Peterson et al. 2016; Correction Table 1)
theta.list  <-  Peterson2016.LEP.theta.list

# inbreeding depression parameters set to 0 due to obligate outcrossing (C = 0)
delta.list  <-  list(delta_D = 0,
					 delta_G = 0,
					 delta_F = 0,
					 delta_O = 0,
					 delta_S = 0
					 )

## Obligate oucrossing
makeLambdaHeatMapMimulusData(sMax = 0.99, len = 100,
							 datMat = datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
							 hf=1/2, hm = 1/2, C = 0, 
							 tlimit = 4*10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
							 nCluster = 12, funs = funs, writeFile = TRUE)

## Empirical selfing rate estimates
## Empirical Inbreeding Depression (delta values calculated from Willis 1993)
delta.list  <-  Willis1993.delta.list

	# Empirical selfing rate
makeLambdaHeatMapMimulusData(sMax = 0.99, len = 100,
							 datMat = datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
							 hf = 1/2, hm = 1/2, C = Willis1993.SelfingRate, 
							 tlimit = 4*10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
							 nCluster = 12, funs = funs, writeFile = TRUE)




###########################################################
#' Fig. B1 - Supplementary density dependence figure.
#' 			 Make data to generate heatmap of lambda 
#' 			 for same parameter conditions as for Fig. 5.

# Export list for dopar function
allObj  <-  ls()
defObj  <-  c("i")
funs    <-  allObj[-which(allObj %in% defObj)]

#################################################
#' Find extinction threshold for Mimulus example
#' TAKING AVERAGE (2012 & 2013) FLOWER # EFFECT 
#' INTO ACCOUNT WHEN CALCULATING MALE SELECTION
#' COEFFICIENTS AND DOMINANCE FOR INV6
# Eagle Meadows population  (Peterson et al. 2016)
datMat      <-  matList_Mg_EM
theta.list  <-  Peterson2016.EM.theta.list

# inbreeding depression parameters set to 0 due to obligate outcrossing (C = 0)
delta.list  <-  list(delta_D = 0,
					 delta_G = 0,
					 delta_F = 0,
					 delta_O = 0,
					 delta_S = 0
					 )
alpha  <-  1e-4

## Obligate oucrossing
makeLambdaHeatMapMimulusData(sMax = 0.99, len = 100,
							 datMat = datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
							 hf = 1/2, hm = 1/2, C = 0, alpha = alpha,
							 tlimit = 4*10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
							 nCluster = 10, funs = funs, writeFile = TRUE)

## Empirical selfing rate estimates
## Empirical Inbreeding Depression (delta values calculated from Willis 1993)
delta.list  <-  Willis1993.delta.list

# Empirical selfing rates
makeLambdaHeatMapMimulusData(sMax = 0.99, len = 100,
							 datMat = datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
							 hf=1/2, hm = 1/2, C = Willis1993.SelfingRate,  alpha = alpha,
							 tlimit = 4*10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
							 nCluster = 12, funs = funs, writeFile = TRUE)


##########################
# Low-elevation perenials (Peterson et al. 2016; Correction Table 1)
theta.list  <-  Peterson2016.LEP.theta.list

# inbreeding depression parameters set to 0 due to obligate outcrossing (C = 0)
delta.list  <-  list(delta_D = 0,
					 delta_G = 0,
					 delta_F = 0,
					 delta_O = 0,
					 delta_S = 0
					 ) 

## Obligate oucrossing
makeLambdaHeatMapMimulusData(sMax = 0.99, len = 100,
							 datMat = datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
							 hf=1/2, hm = 1/2, C = 0,  alpha = alpha,
							 tlimit = 4*10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
							 nCluster = 10, funs = funs, writeFile = TRUE)

## Empirical selfing rate estimates
## Empirical Inbreeding Depression (delta values calculated from Willis 1993)
delta.list  <-  Willis1993.delta.list

	# Empirical selfing rates
makeLambdaHeatMapMimulusData(sMax = 0.99, len = 100,
							 datMat = datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
							 hf = 1/2, hm = 1/2, C = Willis1993.SelfingRate,  alpha = alpha,
							 tlimit = 4*10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
							 nCluster = 12, funs = funs, writeFile = TRUE)






####################################################################
#' Preliminary/Exploratory figures
#' 
#' The remainder of this file contains code to generate exploratory
#' figures & analyses written during the development of the study. 
#' Some  functions are buggy, and some are deprecated, having been 
#' replaced by better/more functional code for the final analyses.
####################################################################


#' Scratch code to iron out effects of inbreeding depression on 
#' invasion boundaries.
# Cs[11]
# deltaSeq  <-  predDelta(dStar=dStar, b=1/2, a=0.2, C=Cs)
 #
# j=11
# invData_f8.5  <-  titrateInvBoundaries(sMax=0.15, res=0.003, precision = 1e-4, alpha = 0,
# 					 om = 2, g = 3, theta = c(0.6,0.6,0.05,8.5), theta_prime = 8.5, 
# 					 hf = 1/2, hm = 1/2, C = Cs[j], 
# 					 delta = 0, delta_j = deltaSeq[j], delta_a = 0, delta_gamma = 0,
# 					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=FALSE)
# invData_f6.5  <-  titrateInvBoundaries(sMax=0.15, res=0.003, precision = 1e-4, alpha = 0,
# 					 om = 2, g = 3, theta = c(0.6,0.6,0.05,8.5), theta_prime = 6.5, 
# 					 hf = 1/2, hm = 1/2, C = Cs[j], 
# 					 delta = 0, delta_j = deltaSeq[j], delta_a = 0, delta_gamma = 0,
# 					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=FALSE)
 #
# pgAinv  <-  popGen_A_invade_Delta_Add(C=Cs[j], delta=deltaSeq[j], sm=invData$sms)
# pgainv  <-  popGen_a_invade_Delta_Add(C=Cs[j], delta=deltaSeq[j], sm=invData$sms)
 #
# plot(pgAinv ~ invData_f8.5$sms, type='l', lwd=2, col=1, ylim=c(0,0.15))
# lines(pgainv ~ invData_f8.5$sms, lwd=2, col=1)
 #
# lines(invData_f8.5$AInvBound ~ invData_f8.5$sms, lwd=2, col=2)
# lines(invData_f8.5$aInvBound ~ invData_f8.5$sms, lwd=2, col=2)
 #
# lines(invData_f6.5$AInvBound ~ invData_f6.5$sms, lwd=2, col=4)
# lines(invData_f6.5$aInvBound ~ invData_f6.5$sms, lwd=2, col=4)


## Additive SA Fitness
## hf = hm = 1/2

	# Obligate Outcrossing
	# C = 0
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  0,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = 6, 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	# Intermediate Selfing
	# C = 1/4
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  1/4,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = 6, 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)


	# Intermediate Selfing
	# C = 1/2
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  1/2,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = 6, 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	# Intermediate Selfing w/ I.D.
	# C = 1/2, delta = 0.2
#	pars  <-  list(
#					"hf"     =  1/2,
#					"hm"     =  1/2,
#					"C"      =  1/2,
#					"delta"  =  1/5,
#					"dj"     =  0,
#					"da"     =  0,
#					"dg"     =  0
#					)

#	selLoop(sMax = 0.15, nSamples=100,
#			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = 6.5, 
#			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
#			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
#			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)




## Partially Recessive SA Fitness
## hf = hm = 1/4

	# Obligate Outcrossing
	# C = 0
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  0,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = 6, 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	# Intermediate Selfing
	# C = 1/4
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  1/4,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = 6, 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	# Intermediate Selfing
	# C = 1/2
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  1/2,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = 6, 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)








############################################
## Plots showing effects of early- vs. late-
## acting inbreeding depression

	# Intermediate Selfing
	# C = 1/2
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  1/2,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	pars$delta  =  0.2

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	pars$delta  =  0.21

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	pars$delta  =  0.22

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	pars$delta  =  0.25

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	pars$delta  =  0.275

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)



	pars$delta  =  0
	pars$dj     =  0.1

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	pars$dj     =  0.11

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	pars$dj     =  0.12

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	pars$dj     =  0.13

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)


	pars$dj  =  0
	pars$da  =  0.1

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	pars$da     =  0.11

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	pars$da     =  0.12

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	pars$da     =  0.13

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	pars$da  =  0
	pars$dg  =  0.2

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	pars$dg     =  0.22

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	pars$dg     =  0.24

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	pars$dg     =  0.26

	selLoop(sMax = 0.15, nSamples=100,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, eqThreshold = 1e-8, writeFile = TRUE)

	# Low Selfing
	# C = 1/4
	# Additive
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  1/4
					)

	deltaParamSpaceMakeData(sMax = 0.15, nSamples=1000,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
							hf = pars$hf, hm = pars$hm, C = pars$C,
							tlimit = 10^5, eqThreshold = 1e-8)

	# Intermediate Selfing
	# C = 1/2
	# Dominance Reversal
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  1/4
					)

	deltaParamSpaceMakeData(sMax = 0.15, nSamples=1000,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
							hf = pars$hf, hm = pars$hm, C = pars$C,
							tlimit = 10^5, eqThreshold = 1e-8)
								# Intermediate Selfing
	# C = 1/2
	# Additive
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  1/2
					)

	deltaParamSpaceMakeData(sMax = 0.15, nSamples=1000,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
							hf = pars$hf, hm = pars$hm, C = pars$C,
							tlimit = 10^5, eqThreshold = 1e-8)

	# Intermediate Selfing
	# C = 1/2
	# Dominance Reversal
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  1/2
					)

	deltaParamSpaceMakeData(sMax = 0.15, nSamples=1000,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
							hf = pars$hf, hm = pars$hm, C = pars$C,
							tlimit = 10^5, eqThreshold = 1e-8)


	# C = 3/4
	# Additive
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  3/4
					)

	deltaParamSpaceMakeData(sMax = 0.15, nSamples=1000,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
							hf = pars$hf, hm = pars$hm, C = pars$C,
							tlimit = 10^5, eqThreshold = 1e-8)

	# Intermediate Selfing
	# C = 1/2
	# Dominance Reversal
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  3/4
					)

	deltaParamSpaceMakeData(sMax = 0.15, nSamples=1000,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
							hf = pars$hf, hm = pars$hm, C = pars$C,
							tlimit = 10^5, eqThreshold = 1e-8)



############################################
#' Plots showing effects of inbreeding
#' depression, using heuristic model of 
#' delta ~ C

#' Additive SA fitness
deltaSelfingPolySpaceMakeData(sMax = 0.15, nSamples=1e+3,
							  om = 2, g = 3, theta = c(0.6,0.6,0.05,7), theta_prime = 7, 
							  hf = 1/2, hm = 1/2, dStar = 0.8, 
							  tlimit = 10^5, eqThreshold = 1e-8)

#' Dominance Reversal SA fitness
deltaSelfingPolySpaceMakeData(sMax = 0.15, nSamples=1e+3,
							  om = 2, g = 3, theta = c(0.6,0.6,0.05,7), theta_prime = 7, 
							  hf = 1/4, hm = 1/4, dStar = 0.8, 
							  tlimit = 10^5, eqThreshold = 1e-8)

############################################
## Plots showing change in polymorphic
## pararmeter space as a function of the
## selfing rate

## "on the edge of demographic viability" (f = 6)
## Additive SA Fitness
## hf = hm = 1/2

	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)
	polyParamSpaceMakeData( sMax = 0.15, nSamples=2000,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
							hf = pars$hf, hm = pars$hm, delta = pars$delta, 
							delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
							tlimit = 10^5, eqThreshold = 1e-8)


## Partially Recessive SA Fitness
## hf = hm = 1/4

	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)
	polyParamSpaceMakeData( sMax = 0.15, nSamples=1000,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
							hf = pars$hf, hm = pars$hm, delta = pars$delta, 
							delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
							tlimit = 10^5, eqThreshold = 1e-8)


## "Poor viability" (f = 5.8)
## Additive SA Fitness
## hf = hm = 1/2

	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)
	polyParamSpaceMakeData( sMax = 0.15, nSamples=2000,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 6, 
							hf = pars$hf, hm = pars$hm, delta = pars$delta, 
							delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
							tlimit = 10^5, eqThreshold = 1e-8)


## Partially Recessive SA Fitness
## hf = hm = 1/4

	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)
	polyParamSpaceMakeData( sMax = 0.15, nSamples=1000,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,5.9), theta_prime = 6, 
							hf = pars$hf, hm = pars$hm, delta = pars$delta, 
							delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
							tlimit = 10^5, eqThreshold = 1e-8)


## "High viability" (f = 6.5)
## Additive SA Fitness
## hf = hm = 1/2

	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)
	polyParamSpaceMakeData( sMax = 0.15, nSamples=2000,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = 6, 
							hf = pars$hf, hm = pars$hm, delta = pars$delta, 
							delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
							tlimit = 10^5, eqThreshold = 1e-8)


## Partially Recessive SA Fitness
## hf = hm = 1/4

	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)
	polyParamSpaceMakeData( sMax = 0.15, nSamples=1000,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = 6, 
							hf = pars$hf, hm = pars$hm, delta = pars$delta, 
							delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
							tlimit = 10^5, eqThreshold = 1e-8)





#################################################
#' Old version of what is now fig. 5, using 
#' extinciton thresholds instead of heatmaps for
#' lambda

#################################################
#' Find extinction threshold for Mimulus example
#' TAKING AVERAGE (2012 & 2013) FLOWER # EFFECT 
#' INTO ACCOUNT WHEN CALCULATING MALE SELECTION
#' COEFFICIENTS AND DOMINANCE FOR INV6
# Eagle Meadows population  (Peterson et al. 2016)
datMat  <-  matList_Mg_EM
theta.list  <-  list(D = 0.534, 
					 G = 0.469,
					 F = 0.64,
					 O = 614,
					 A = 6.7e-4,
					 S = 0.179,
					 R = 8.71,
					 pop = "EM")
delta.list  <-  list(delta_D = 0,
					 delta_G = 0,
					 delta_F = 0,
					 delta_O = 0,
					 delta_S = 0
					 )


## Obligate oucrossing
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm = 1/2, C=0, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 1/2, C = 0, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 makePlots=TRUE, writeFile=TRUE, verbose=TRUE) 

## Empirical selfing rate estimates
## Empirical Inbreeding Depression (delta values calculated from Willis 1993)
delta.list  <-  list(delta_D = 0,
					 delta_G = 0.085,
					 delta_F = 0.2,
					 delta_O = 0,
					 delta_S = 0.38
					 )
	# Low-end of Empirical selfing rates
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm = 1/2, C=0.29, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 1/2, C = 0.29, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 writeFile=TRUE, verbose=TRUE) 


##########################
# Low-elevation perenials (Peterson et al. 2016; Correction Table 1)
theta.list  <-  list(D = 0.534, 
					 G = 0.652,
					 F = 4.09,
					 O = 494,
					 A = 6.7e-4,
					 S = 0,
					 R = 0,
					 pop = "LEP")
delta.list  <-  list(delta_D = 0,
					 delta_G = 0,
					 delta_F = 0,
					 delta_O = 0,
					 delta_S = 0
					 )


## Obligate oucrossing
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm = 1/2, C=0, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 1/2, C = 0, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 makePlots=TRUE, writeFile=TRUE, verbose=TRUE) 

## Empirical selfing rate estimates
## Empirical Inbreeding Depression (delta values calculated from Willis 1993)
delta.list  <-  list(delta_D = 0,
					 delta_G = 0.085,
					 delta_F = 0.2,
					 delta_O = 0,
					 delta_S = 0.38
					 )
	# Empirical selfing rates
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm = 1/2, C=0.29, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 1/2, C = 0.29, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 writeFile=TRUE, verbose=TRUE) 


#################################################
#' Find extinction threshold for Mimulus example

# Eagle Meadows population  (Peterson et al. 2016)
datMat  <-  matList_Mg_EM
theta.list  <-  list(D = 0.534, 
					 G = 0.469,
					 F = 0.64,
					 O = 614,
					 A = 6.7e-4,
					 S = 0.179,
					 R = 8.71,
					 pop = "EM")
delta.list  <-  list(delta_D = 0,
					 delta_G = 0,
					 delta_F = 0,
					 delta_O = 0,
					 delta_S = 0
					 )


## Obligate oucrossing
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm=0.35, C=0, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.35, C = 0, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 makePlots=TRUE, writeFile=TRUE, verbose=TRUE) 

## Empirical selfing rate estimates
## No Inbreeding Depression
	# Low-end of Empirical selfing rates
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm=0.35, C=0.29, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.35, C = 0.29, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 makePlots=TRUE, writeFile=TRUE, verbose=TRUE) 

	# High-end of Empirical selfing rates
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm=0.35, C=0.75, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.35, C = 0.75, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 makePlots=TRUE, writeFile=TRUE, verbose=TRUE) 
## Empirical selfing rate estimates
## Empirical Inbreeding Depression (delta values calculated from Willis 1993)
delta.list  <-  list(delta_D = 0,
					 delta_G = 0.085,
					 delta_F = 0.2,
					 delta_O = 0,
					 delta_S = 0.38
					 )
	# Low-end of Empirical selfing rates
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm=0.35, C=0.29, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.35, C = 0.29, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 writeFile=TRUE, verbose=TRUE) 
	# High-end of Empirical selfing rates
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm=0.35, C=0.75, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.35, C = 0.75, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 writeFile=TRUE, verbose=TRUE) 


##########################
# Low-elevation perenials (Peterson et al. 2016; Correction Table 1)
theta.list  <-  list(D = 0.534, 
					 G = 0.652,
					 F = 4.09,
					 O = 494,
					 A = 6.7e-4,
					 S = 0,
					 R = 0,
					 pop = "LEP")
delta.list  <-  list(delta_D = 0,
					 delta_G = 0,
					 delta_F = 0,
					 delta_O = 0,
					 delta_S = 0
					 )


## Obligate oucrossing
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm=0.35, C=0, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.35, C = 0, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 makePlots=TRUE, writeFile=TRUE, verbose=TRUE) 

## Empirical selfing rate estimates
## No Inbreeding Depression
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm=0.35, C=0.29, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.35, C = 0.29, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 makePlots=TRUE, writeFile=TRUE, verbose=TRUE) 


## Empirical selfing rate estimates
## Empirical Inbreeding Depression (delta values calculated from Willis 1993)
delta.list  <-  list(delta_D = 0,
					 delta_G = 0.085,
					 delta_F = 0.2,
					 delta_O = 0,
					 delta_S = 0.38
					 )
	# Empirical selfing rates
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm=0.35, C=0.29, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.35, C = 0.29, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 writeFile=TRUE, verbose=TRUE) 
	





	





#################################################
#' Find extinction threshold for Mimulus example
#' TAKING 2013 FLOWER # EFFECT INTO ACCOUNT WHEN 
#' CALCULATING MALE SELECTION COEFFICIENTS.
# Eagle Meadows population  (Peterson et al. 2016)
datMat  <-  matList_Mg_EM
theta.list  <-  list(D = 0.534, 
					 G = 0.469,
					 F = 0.64,
					 O = 614,
					 A = 6.7e-4,
					 S = 0.179,
					 R = 8.71,
					 pop = "EM")
delta.list  <-  list(delta_D = 0,
					 delta_G = 0,
					 delta_F = 0,
					 delta_O = 0,
					 delta_S = 0
					 )


## Obligate oucrossing
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm = 0.05, C=0, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.05, C = 0, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 makePlots=TRUE, writeFile=TRUE, verbose=TRUE) 

## Empirical selfing rate estimates
## No Inbreeding Depression
	# Low-end of Empirical selfing rates
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm = 0.05, C=0.29, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.05, C = 0.29, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 makePlots=TRUE, writeFile=TRUE, verbose=TRUE) 

## Empirical selfing rate estimates
## Empirical Inbreeding Depression (delta values calculated from Willis 1993)
delta.list  <-  list(delta_D = 0,
					 delta_G = 0.085,
					 delta_F = 0.2,
					 delta_O = 0,
					 delta_S = 0.38
					 )
	# Low-end of Empirical selfing rates
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm = 0.05, C=0.29, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.05, C = 0.29, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 writeFile=TRUE, verbose=TRUE) 


##########################
# Low-elevation perenials (Peterson et al. 2016; Correction Table 1)
theta.list  <-  list(D = 0.534, 
					 G = 0.652,
					 F = 4.09,
					 O = 494,
					 A = 6.7e-4,
					 S = 0,
					 R = 0,
					 pop = "LEP")
delta.list  <-  list(delta_D = 0,
					 delta_G = 0,
					 delta_F = 0,
					 delta_O = 0,
					 delta_S = 0
					 )


## Obligate oucrossing
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm = 0.05, C=0, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.05, C = 0, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 makePlots=TRUE, writeFile=TRUE, verbose=TRUE) 

## Empirical selfing rate estimates
## No Inbreeding Depression
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm = 0.05, C=0.29, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.05, C = 0.29, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 makePlots=TRUE, writeFile=TRUE, verbose=TRUE) 


## Empirical selfing rate estimates
## Empirical Inbreeding Depression (delta values calculated from Willis 1993)
delta.list  <-  list(delta_D = 0,
					 delta_G = 0.085,
					 delta_F = 0.2,
					 delta_O = 0,
					 delta_S = 0.38
					 )
	# Empirical selfing rates
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm = 0.05, C=0.29, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.05, C = 0.29, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 writeFile=TRUE, verbose=TRUE) 
	










##########################
# Iron Mountain (Willis 1993; Table 2)
theta.list  <-  list(D = 0.534, 
					 G = 0.906,
					 F = 0.89,
					 O = 28.628,
					 A = 6.7e-4,
					 S = 0.435,
					 R = 0,
					 pop = "IM")
delta.list  <-  list(delta_D = 0,
					 delta_G = 0,
					 delta_F = 0,
					 delta_O = 0,
					 delta_S = 0
					 )


## Obligate oucrossing
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm=0.35, C=0, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.35, C = 0, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 makePlots=TRUE, writeFile=TRUE, verbose=TRUE) 

## Empirical selfing rate estimates
## No Inbreeding Depression
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm=0.35, C=0.29, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.35, C = 0.29, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 makePlots=TRUE, writeFile=TRUE, verbose=TRUE) 


## Empirical selfing rate estimates
## Empirical Inbreeding Depression (delta values calculated from Willis 1993)
delta.list  <-  list(delta_D = 0,
					 delta_G = 0.085,
					 delta_F = 0.2,
					 delta_O = 0,
					 delta_S = 0.38
					 )
	# Empirical selfing rates
titrateInvBoundMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
						datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = FALSE,
						hf=1/2, hm=0.35, C=0.29, 
						tlimit = 10^5, eqThreshold = 1e-8, writeFile=TRUE)
extinctThreshMimulus(sMax = 0.99, res=0.01, precision=1e-4, 
					 datMat=datMat, theta.list = theta.list, delta.list = delta.list, useCompadre = FALSE,
					 hf = 1/2, hm = 0.35, C = 0.29, 
					 tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
					 writeFile=TRUE, verbose=TRUE) 



##########################################################
##########################################################
# Supplementary Figures
##########################################################
##########################################################

## Illustrating Effects of sex-specific dominance
# Assume complementary dominance (hf = 1 - hm)
# Fig. SX - 

## hf 1/4
## hm = 3/4 

	# Obligate Outcrossing
	# C = 0
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  3/4,
					"C"      =  0,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

titrateInvBoundaries(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^3, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6.2), theta_prime = 6.2, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)


	# Intermediate Selfing
	# C = 1/4
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  3/4,
					"C"      =  1/4,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

titrateInvBoundaries(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^3, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6.2), theta_prime = 6.2, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

	# Intermediate Selfing
	# C = 1/2
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  3/4,
					"C"      =  1/2,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

titrateInvBoundaries(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^3, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6.2), theta_prime = 6.2, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)


## hf 3/4
## hm = 1/4 

	# Obligate Outcrossing
	# C = 0
	pars  <-  list(
					"hf"     =  3/4,
					"hm"     =  1/4,
					"C"      =  0,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

titrateInvBoundaries(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^3, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6.2), theta_prime = 6.2, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)


	# Intermediate Selfing
	# C = 1/4
	pars  <-  list(
					"hf"     =  3/4,
					"hm"     =  1/4,
					"C"      =  1/4,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

titrateInvBoundaries(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^3, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6.2), theta_prime = 6.2, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

	# Intermediate Selfing
	# C = 1/2
	pars  <-  list(
					"hf"     =  3/4,
					"hm"     =  1/4,
					"C"      =  1/2,
					"delta"  =  0,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

titrateInvBoundaries(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^3, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)

extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 5.8, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)
extinctThreshTitrate(sMax=0.15, res=0.0015, precision=1e-4,
					 om = 2, g = 3, theta = c(0.6,0.6,0.05,6.2), theta_prime = 6.2, 
					 hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
					 delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg,
					 tlimit = 10^5, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE)




##########################################################
# Fig. SX - Dem. viable polymorphic param. space (no I.D.)
#			Sex-specific dominance

makeDataPolyParamSpace(sMax=0.15, res=0.003, precision = 1e-4,
						om = 2, g = 3, theta = c(0.6,0.6,0.05,NA), theta_prime = NA, 
						hfVals= c(1/4, 3/4), hmVals= c(3/4, 1/4), fVals = c(5.8, 5.9, 6.0, 6.5), 
						delta = 0, delta_j = 0, delta_a = 0, delta_gamma = 0,
						tlimit = 10^5, eqThreshold = 1e-8)



##########################################################
# Fig. SX- Dem. viable polymorphic param. space (w/ I.D.)

makeDataDeltaPolyParamSpace(sMax=0.15, res=0.003, precision = 1e-4,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = 6.5, 
							hfVals= c(1/4, 3/4), hmVals= c(3/4, 1/4), dStar = 0.8, 
							tlimit = 10^5, eqThreshold = 1e-8)

makeDataDeltaPolyParamSpace(sMax=0.15, res=0.003, precision = 1e-4,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,7.5), theta_prime = 7.5, 
							hfVals= c(1/4, 3/4), hmVals= c(3/4, 1/4), dStar = 0.8, 
							tlimit = 10^5, eqThreshold = 1e-8)

makeDataDeltaPolyParamSpace(sMax=0.15, res=0.003, precision = 1e-4,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,8.5), theta_prime = 8.5, 
							hfVals= c(1/4, 3/4), hmVals= c(3/4, 1/4), dStar = 0.8, 
							tlimit = 10^5, eqThreshold = 1e-8)
