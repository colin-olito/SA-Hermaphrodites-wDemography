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
					"delta"  =  0
					)

	selLoop(sMax = 0.15, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.6,0.6,0.05,0.05,6), theta_prime = c(0.6,0.6,0.6,0.6,0.05,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, tlimit = 10^5, returnRes = FALSE)

	# Intermediate Selfing
	# C = 1/2
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  1/2,
					"delta"  =  0
					)

	selLoop(sMax = 0.15, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.6,0.6,0.05,0.05,6), theta_prime = c(0.6,0.6,0.6,0.6,0.05,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, tlimit = 10^5, returnRes = FALSE)

	# Intermediate Selfing w/ I.D.
	# C = 1/2, delta = 0.2
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  1/2,
					"delta"  =  1/5
					)

	selLoop(sMax = 0.15, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.6,0.6,0.05,0.05,6), theta_prime = c(0.6,0.6,0.6,0.6,0.05,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, tlimit = 10^5, returnRes = FALSE)




## Partially Recessive SA Fitness
## hf = hm = 1/4

	# Obligate Outcrossing
	# C = 0
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  0,
					"delta"  =  0
					)

	selLoop(sMax = 0.15, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.6,0.6,0.05,0.05,6), theta_prime = c(0.6,0.6,0.6,0.6,0.05,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, tlimit = 10^5, returnRes = FALSE)

	# Intermediate Selfing
	# C = 1/2
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  1/2,
					"delta"  =  0
					)

	selLoop(sMax = 0.15, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.6,0.6,0.05,0.05,6), theta_prime = c(0.6,0.6,0.6,0.6,0.05,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, tlimit = 10^5, returnRes = FALSE)

	# Intermediate Selfing w/ I.D.
	# C = 1/2, delta = 0.5
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  1/2,
					"delta"  =  1/5
					)

	selLoop(sMax = 0.15, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.6,0.6,0.05,0.05,6), theta_prime = c(0.6,0.6,0.6,0.6,0.05,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, tlimit = 10^5, returnRes = FALSE)









## Same plots showing larger s_f x s_m parameter space

## Additive SA Fitness
## hf = hm = 1/2

	# Obligate Outcrossing
	# C = 0
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  0,
					"delta"  =  0
					)

	selLoop(sMax = 1, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.6,0.6,0.05,0.05,6), theta_prime = c(0.6,0.6,0.6,0.6,0.05,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, tlimit = 10^5, returnRes = FALSE)

	# Intermediate Selfing
	# C = 1/2
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  1/2,
					"delta"  =  0
					)

	selLoop(sMax = 1, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.6,0.6,0.05,0.05,6), theta_prime = c(0.6,0.6,0.6,0.6,0.05,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, tlimit = 10^5, returnRes = FALSE)

	# Intermediate Selfing w/ I.D.
	# C = 1/2, delta = 0.25
	pars  <-  list(
					"hf"     =  1/2,
					"hm"     =  1/2,
					"C"      =  1/2,
					"delta"  =  1/4
					)

	selLoop(sMax = 1, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.6,0.6,0.05,0.05,6), theta_prime = c(0.6,0.6,0.6,0.6,0.05,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, tlimit = 10^5, returnRes = FALSE)





## Partially Recessive SA Fitness
## hf = hm = 1/4

	# Obligate Outcrossing
	# C = 0
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  0,
					"delta"  =  0
					)

	selLoop(sMax = 1, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.6,0.6,0.05,0.05,6), theta_prime = c(0.6,0.6,0.6,0.6,0.05,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, tlimit = 10^5, returnRes = FALSE)

	# Intermediate Selfing
	# C = 1/2
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  1/2,
					"delta"  =  0
					)

	selLoop(sMax = 1, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.6,0.6,0.05,0.05,6), theta_prime = c(0.6,0.6,0.6,0.6,0.05,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, tlimit = 10^5, returnRes = FALSE)

	# Intermediate Selfing w/ I.D.
	# C = 1/2, delta = 0.25
	pars  <-  list(
					"hf"     =  1/4,
					"hm"     =  1/4,
					"C"      =  1/2,
					"delta"  =  1/4
					)

	selLoop(sMax = 1, nSamples=1000,
			om = 2, g = 3, theta = c(0.6,0.6,0.6,0.6,0.05,0.05,6), theta_prime = c(0.6,0.6,0.6,0.6,0.05,0.05,6), 
			hf = pars$hf, hm = pars$hm, C = pars$C, delta = pars$delta, tlimit = 10^5, returnRes = FALSE)

