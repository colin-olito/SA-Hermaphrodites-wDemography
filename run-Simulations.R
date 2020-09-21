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

	selLoop(sMax = 0.15, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = c(0.6,0.6,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, writeFile = TRUE)

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

	selLoop(sMax = 0.15, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = c(0.6,0.6,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, writeFile = TRUE)

	# Intermediate Selfing w/ I.D.
	# C = 1/2, delta = 0.2
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  1/2,
					"delta"  =  1/5,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

	selLoop(sMax = 0.15, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = c(0.6,0.6,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, writeFile = TRUE)




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

	selLoop(sMax = 0.15, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = c(0.6,0.6,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, writeFile = TRUE)

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

	selLoop(sMax = 0.15, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = c(0.6,0.6,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, writeFile = TRUE)

	# Intermediate Selfing w/ I.D.
	# C = 1/2, delta = 0.5
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  1/2,
					"delta"  =  1/5,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

	selLoop(sMax = 0.15, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = c(0.6,0.6,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, writeFile = TRUE)









## Same plots showing larger s_f x s_m parameter space

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

	selLoop(sMax = 1, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = c(0.6,0.6,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, writeFile = TRUE)

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

	selLoop(sMax = 1, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = c(0.6,0.6,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, writeFile = TRUE)

	# Intermediate Selfing w/ I.D.
	# C = 1/2, delta = 0.25
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  1/2,
					"delta"  =  1/4,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

	selLoop(sMax = 1, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = c(0.6,0.6,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, writeFile = TRUE)





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

	selLoop(sMax = 1, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = c(0.6,0.6,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, writeFile = TRUE)

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

	selLoop(sMax = 1, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = c(0.6,0.6,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, writeFile = TRUE)

	# Intermediate Selfing w/ I.D.
	# C = 1/2, delta = 0.25
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  1/2,
					"delta"  =  1/4,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

	selLoop(sMax = 1, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = c(0.6,0.6,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, writeFile = TRUE)



############################################
## Plots showing effects of early- vs. late-
## acting inbreeding depression

	# Intermediate Selfing
	# C = 1/2
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  1/2,
					"delta"  =  0.1,
					"dj"     =  0,
					"da"     =  0,
					"dg"     =  0
					)

	selLoop(sMax = 0.15, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, 
			delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
			tlimit = 10^5, writeFile = TRUE)


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
	polyParamSpaceMakeData( sMax = 0.15, nSamples=5000,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = c(0.6,0.6,0.05,6), 
							hf = pars$hf, hm = pars$hm, delta = pars$delta, 
							delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
							tlimit = 10^5)


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
							om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = c(0.6,0.6,0.05,6), 
							hf = pars$hf, hm = pars$hm, delta = pars$delta, 
							delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
							tlimit = 10^5)


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
	polyParamSpaceMakeData( sMax = 0.15, nSamples=5000,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = c(0.6,0.6,0.05,5.8), 
							hf = pars$hf, hm = pars$hm, delta = pars$delta, 
							delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
							tlimit = 10^5)


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
							om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = c(0.6,0.6,0.05,5.75), 
							hf = pars$hf, hm = pars$hm, delta = pars$delta, 
							delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
							tlimit = 10^5)


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
	polyParamSpaceMakeData( sMax = 0.15, nSamples=5000,
							om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
							hf = pars$hf, hm = pars$hm, delta = pars$delta, 
							delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
							tlimit = 10^5)


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
							om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = c(0.6,0.6,0.05,6.5), 
							hf = pars$hf, hm = pars$hm, delta = pars$delta, 
							delta_j = pars$dj, delta_a = pars$da, delta_gamma = pars$dg, 
							tlimit = 10^5)

