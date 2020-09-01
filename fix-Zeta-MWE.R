##  Minimum working example of disagreement b/w
##   Simulations and Zeta_i


# run simulations necessary to show disagreement
# (can take 10 min  or so)
rm(list=ls())
source('R/functions-Simulations.R')

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

# plotting functions
source('./R/functions-Figs.R')
#  Plot to show that when C > 0, zeta_i never predicts polymorphism
FunnelEigSimCompare(df1 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.5_hm0.5_C0_delta0',
					df2 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.5_hm0.5_C0.5_delta0',
					df3 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.25_hm0.25_C0_delta0',
					df4 = 'demSimsSfxSm_MinorEqAlleleInv_sMax0.15_nSamples1000_hf0.25_hm0.25_C0.5_delta0')



###################
# enter parameters so we can step through a simulation 
# the h and sf, sm values should land us in polymorphic
# parameter space
om = 2
g = 3
theta = c(0.58, 0.6, 0.6, 0.6, 0.05, 0.05, 5.9)
theta_prime = c(0.6, 0.6, 0.6, 0.6, 0.05, 0.05, 5.9)
hf = 1/4
hm = 1/4
sf = 0.03
sm = 0.1
C = 1/2
delta = 0
tlimit = 10^4
Ainvade = FALSE
intInit = FALSE


# Uncomment and run to confirm that zeta_i and simulations
# agree when C = 0
# sf = 0.05
# sm = 0.051
# C = 0

# Create identity and ones matrices
	Ig   <-  diag(g)
	Iom  <-  diag(om)
	eg   <-  ones(c(g,1))
	eom  <-  ones(c(om,1))
	K    <-  vecperm(om,g)

	####################################################
	# PARAMETERS AND INITIAL CONDITIONS 
	####################################################

	# BASELINE PARAMETERS for survival and 
	# fertility through Female and Male function
	# theta=[sigmaS_J, sigmaS_A, sigmaX_J, sigmaX_A, gammaS, gammaX, f_ii]
	theta       <-  rep.col(theta,3)
	theta_prime <- rep.col(theta_prime,3)

	# Fitness Expressions
	fii        <-  c(1, 1 - hf*sf, 1 - sf)
	fii_prime  <-  c(1 - sm, 1 - hm*sm, 1)

	##SELECTION DIFFERENTIAL FEMALE-FUNCTION
	theta[7,]  <-  theta[7,]*fii

	#SELECTION DIFFERENTIAL MALE-FUNCTION
	theta_prime[7,]  <-  theta_prime[7,]*fii_prime

	#Initial conditions for invasion: nzero  <-  [juv_AA, ad_AA, juv_Aa, ad_Aa, juv_aa, ad_aa]
	# All aa
		if(Ainvade) {
			nzero  <- round(c(C*(1 - delta)*100*c(0,0,0,0,141.3,18.2),(1 - C)*100*c(0,0,0,0,141.3,18.2)))
#			nzero  <- round(c((C*(1 - delta)/(1 - C*delta))*100*c(0,0,0,0,141.3,18.2),(1 - (C*(1 - delta)/(1 - C*delta)))*100*c(0,0,0,0,141.3,18.2)))
		} 
		#All AA
		if(!Ainvade) {
			nzero  <-  round(c(C*(1 - delta)*100*c(141.3,18.2,0,0,0,0),(1 - C)*100*c(141.3,18.2,0,0,0,0)))
#			nzero  <-  round(c((C*(1 - delta)/(1 - C*delta))*100*c(141.3,18.2,0,0,0,0),(1 - (C*(1 - delta)/(1 - C*delta)))*100*c(141.3,18.2,0,0,0,0)))
		}
		#All Aa
		if(intInit) {
			nzero  <-  round(c(C*(1 - delta)*100*c(0,0,141.3,18.2,0,0),(1 - C)*100*c(0,0,141.3,18.2,0,0)))
#			nzero  <-  round(c((C*(1 - delta)/(1 - C*delta))*100*c(0,0,141.3,18.2,0,0),(1 - (C*(1 - delta)/(1 - C*delta)))*100*c(0,0,141.3,18.2,0,0)))
		}


	####################################################
	# Population and female fertility perameters
	####################################################

	sigmaS_J   <-  theta[1,]
	sigmaS_A   <-  theta[2,]
	sigmaX_J   <-  theta[3,]
	sigmaX_A   <-  theta[4,]
	gammaS     <-  theta[5,]
	gammaX     <-  theta[6,]
	f          <-  theta[7,]
	f_prime    <-  theta_prime[7,]
	USi        <-  zeros(c(om,om,g))
	UXi        <-  zeros(c(om,om,g))
	Fii        <-  zeros(c(om,om,g))
	Fii_pr     <-  zeros(c(om,om,g))
	Nfun       <-  zeros(c(om,om,g))
	R          <-  zeros(c(om,om,g))
	lambda     <-  zeros(c(1,g))
	R_0        <-  zeros(c(1,g))
	eigs_temp  <-  zeros(c(g,2))

	# create genotype-specific survival and fertility submatrices
	for (i in 1:3){
	    USi[,,i]       <- rbind(c(sigmaS_J[i]*(1 - gammaS[i]), 0         ),
	                            c(sigmaS_J[i]*gammaS[i],       sigmaS_A[i]))
	    UXi[,,i]       <- rbind(c(sigmaX_J[i]*(1 - gammaX[i]), 0         ),
	                            c(sigmaX_J[i]*gammaX[i],       sigmaX_A[i]))
	    Fii[,,i]       <- rbind(c(0,f[i]),
						        c(0,0))
	    Fii_pr[,,i]    <- rbind(c(0,f_prime[i]),
								c(0,0))
    
		# genotype-specific population growth rate
		eigs_temp[i,]  <-  eigen((C*(1 - delta)*Fii[,,i] + USi[,,i]) + ((1 - C)*Fii[,,i] + UXi[,,i]), symmetric=FALSE, only.values = TRUE)$values
		lambda[i]      <-  max(eigs_temp)
		Nfun[,,i]      <-  solve((diag(2) - (C*(1 - delta)*USi[,,i] + (1 - C)*UXi[,,i])), diag(2))
#		Nfun[,,i]      <-  solve((diag(2) - ((C*(1 - delta)/(1 - C*delta))*USi[,,i] + (1 - (C*(1 - delta)/(1 - C*delta)))*UXi[,,i])), diag(2))
		R[,,i]         <-  (C*(1 - delta)*Fii[,,i] + (1 - C)*Fii[,,i]) %*% Nfun[,,i]
		#genotype-specific R0
		R_0[i] <- max(eigen(R[,,i],symmetric=FALSE, only.values = TRUE)$values)
	}

	Atilde_genotype  <-  zeros(c(2*om,2*om,g))
	lambda_full      <-  as.vector(zeros(c(1,g)))
	for (i in 1:3){
		Atilde_genotype[,,i]  <- rbind( cbind(USi[,,i] + C*(1 - delta)*Fii[,,i], C*(1 - delta)*Fii[,,i]), 
										cbind((1 - C)*Fii[,,i], UXi[,,i] + (1 - C)*Fii[,,i])
								   	   )
#		Atilde_genotype[,,i]  <- rbind( cbind(USi[,,i] + (C*(1 - delta)/(1 - C*delta))*Fii[,,i], (C*(1 - delta)/(1 - C*delta))*Fii[,,i]), 
#										cbind((1 - (C*(1 - delta)/(1 - C*delta)))*Fii[,,i], UXi[,,i] + (1 - (C*(1 - delta)/(1 - C*delta)))*Fii[,,i])
#								   	   )
		lambda_full[i]    <- max(eigen(Atilde_genotype[,,i],symmetric=FALSE, only.values = TRUE)$values)
	}

	# CREATE BLOCK DIAGONAL MATRICES
	d 			 <-  diag(g)
	blkD         <-  kronecker(Iom,d)
	blkUS        <-  zeros(c(g*om,g*om))
	blkUX        <-  zeros(c(g*om,g*om))
	blkFS        <-  zeros(c(g*om,g*om))
	blkFX        <-  zeros(c(g*om,g*om))
	blkFX_prime  <-  zeros(c(g*om,g*om))
	for (i in 1:g){
	    blkUS        <-  blkUS + kronecker(emat(g,g,i,i),USi[,,i])
	    blkUX        <-  blkUX + kronecker(emat(g,g,i,i),UXi[,,i])
	    blkFS        <-  blkFS + kronecker(emat(g,g,i,i),C*(1 - delta)*Fii[,,i])
	    blkFX        <-  blkFX + kronecker(emat(g,g,i,i),(1 - C)*Fii[,,i])
	    blkFX_prime  <-  blkFX_prime + kronecker(emat(g,g,i,i),(1 - C)*Fii_pr[,,i])
	}
	W_prime  <-  rbind(cbind(ones(c(1,om)),.5*ones(c(1,om)),0*ones(c(1,om))),
					   cbind(0*ones(c(1,om)),.5*ones(c(1,om)),ones(c(1,om))))
	W <- rbind( c(1,0.5,0),
				c(0,0.5,1)
				)
	Z <- rbind( c(1,0,0,0),
				c(0,1,1,0),
				c(0,0,0,1)
				)

	# INITIAL CONDITIONS
	nzero  <- c(nzero)
	n      <- t(t(nzero))
	nout   <- zeros(c((2*om*g),tlimit))
	pout   <- zeros(c(g,tlimit))

	###############################################################
	# Simulating the Stage X genotype dynamics - generation loop
	###############################################################
	i  <-  1
	tmp       <-  kronecker(diag(g),ones(c(1,om))) %*% (nzero[1:6] + nzero[7:12])
	p         <-  tmp/colSums(tmp)
	pDelta    <-  c(1,1,1)

		nout[,i]  <-  n

        #Creating the male gamete pool
        nX    <-  n[1:6] + n[7:12]
        ngam  <-  ones(c(1,2)) %*% W_prime %*% blkFX_prime %*% nX
        
        #Equation 5 in the manuscript
        q_prime <- (W_prime%*%blkFX_prime%*%nX)/ngam[1] 

        UtildeS  <-  blkUS
        UtildeX  <-  blkUX

        # Parent-Offspring genotype map - Selfing
        HS <- rbind(c(1, 1/4, 0),
					c(0, 1/2, 0),
					c(0, 1/4, 1))
        # Parent-Offspring genotype map - Outcrossing
        HX <- zeros(c(g,g))

        for (ii in 1:g){
            Pi  <-  Ig[,ii]
            qi  <-  W %*% Pi #allele frequencies in oocytes of genotype i

            # genotype frequencies in the offspring of mothers of genotype i
            # produced by outcrossing (equation 16 in de Vries and Caswell, 2018a (American Naturalist))
            Piprime  <-  Z %*% kronecker(qi,q_prime)

            HX[,ii]  <-  Piprime # the outcrossing parent-offspring matrix
        }

        blkHS    <-  kronecker(Iom,HS)
        blkHX    <-  kronecker(Iom,HX)

        FtildeS  <-  t(K) %*% blkHS %*% K %*% blkFS
        FtildeX  <-  t(K) %*% blkHX %*% K %*% blkFX
        Atilde   <-  rbind(cbind((UtildeS + FtildeS), FtildeS),
						   cbind(FtildeX            , (UtildeX + FtildeX)))

		# Project population into next generation
        nnext  <-  Atilde %*% n


########################################
# Below here is where we see the problem
########################################

# initial frequencies for the two boundaries 
# where AA genotype is fixed, and where aa genotype is fixed
nBoundAA  <-  round(c(C*(1 - delta)*100*c(141.3,18.2,0,0,0,0),(1 - C)*100*c(141.3,18.2,0,0,0,0)))
nBoundaa  <-  round(c(C*(1 - delta)*100*c(0,0,0,0,141.3,18.2),(1 - C)*100*c(0,0,0,0,141.3,18.2)))
pBoundAA  <-  nBoundAA/norm(as.matrix(nBoundAA), type="1")
pBoundaa  <-  nBoundaa/norm(as.matrix(nBoundaa), type="1")

# we rearrange the order of the phat values so they match
# the structure of our jacobian, which was ordered by genotype, then self/outcross
pHat_AA    <-  pBoundAA[c(1,2,7,8,3,4,9,10,5,6,11,12)]
pHat_aa    <-  pBoundaa[c(1,2,7,8,3,4,9,10,5,6,11,12)]

# calculate Atilde[ptilde]
testAtilde_AA  <-  AtildeBound(nBound=pBoundAA, W_prime=W_prime, blkFX_prime = blkFX_prime, blkUS=blkUS, blkUX=blkUX, W=W, Iom=Iom, HS=HS, HX=HX, K=K, blkFS=blkFS,blkFX=blkFX)
testAtilde_aa  <-  AtildeBound(nBound=pBoundaa, W_prime=W_prime, blkFX_prime=blkFX_prime, blkUS=blkUS, blkUX=blkUX, W=W, Iom=Iom, HS=HS, HX=HX, K=K, blkFS=blkFS,blkFX=blkFX)

# Calculate a naive genotype-specific eigenvalue for comparison
# note, that this gives the same value for lambda_i as we get 
# elsewhere in the simulation (see lambda_full)
testLambda_AA  <-  c(max(eigen(testAtilde_AA[c(1:4),c(1:4)],symmetric=FALSE, only.values = TRUE)$values),
				     max(eigen(testAtilde_AA[c(5:8),c(5:8)],symmetric=FALSE, only.values = TRUE)$values),
				     max(eigen(testAtilde_AA[c(9:12),c(9:12)],symmetric=FALSE, only.values = TRUE)$values)
					)
testLambda_aa  <-  c(max(eigen(testAtilde_aa[c(1:4),c(1:4)],symmetric=FALSE, only.values = TRUE)$values),
				     max(eigen(testAtilde_aa[c(5:8),c(5:8)],symmetric=FALSE, only.values = TRUE)$values),
				     max(eigen(testAtilde_aa[c(9:12),c(9:12)],symmetric=FALSE, only.values = TRUE)$values)
					)

# Calcuate lambda_AA and lambda_aa both using the as described in the notes file
# (Eq(17) in the linearization at the boundary equilibrium subsection)
lambda_i   <-  c(testLambda_AA[1], testLambda_aa[3])
lambda2_i  <-  c(t(ones(2*om*g)) %*% testAtilde_AA %*% pHat_AA,
				 t(ones(2*om*g)) %*% testAtilde_aa %*% pHat_aa)

# Calculate zeta_i, dividing be both of these versions of lambda_i, 
# and see that they do not agree... interestingly, the result using lambda2_i also
# doesn't work when C = 0
calcZeta(om=om, Fii=Fii, Fii_pr=Fii_pr, USi=USi, UXi=UXi,
		 pHat_AA=pHat_AA, pHat_aa=pHat_aa, C=C, delta=delta, lambda_AA_aa=lambda_i)
calcZeta(om=om, Fii=Fii, Fii_pr=Fii_pr, USi=USi, UXi=UXi,
		 pHat_AA=pHat_AA, pHat_aa=pHat_aa, C=C, delta=delta, lambda_AA_aa=lambda2_i)


# now run full simulation
sim  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
						hf = hf, hm = hm, sf = sf, sm = sm, C = C, delta = delta, tlimit = 10^4)
sim$poly
sim$pEq



## So, I'm not sure what's going on here. I probably just have a minor
## mistake or typo somewhere, but I've struggled to find it, and think
## I need a fresh pair of eyest to have a look... sorry Lotte.


