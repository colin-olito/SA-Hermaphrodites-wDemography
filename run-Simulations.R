################################################################
#  RUN SIMULATIONS FOR 
#
#  R code for forward simulations. Generates output data
#  as .csv files saved to ./output/simData/
#
#
#  Author: Colin Olito
#
#  NOTES:  
#		


#####################
##  Dependencies
rm(list=ls())
source('R/functions-Simulations.R')
#source('R/loadData-Compadre.R')



######################
# Preliminary figures
######################

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

	# Intermediate Selfing w/ I.D.
	# C = 1/2, delta = 0.5
#	pars  <-  list(
#					"hf"     =  1/4,
#					"hm"     =  1/4,
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

