################################################################
#  Functions to simulate demographic model across selection
#  parameter space
#
#
#  Author: Colin Olito, adapted from L. DeVries
#
#  NOTES:  
#		


rm(list=ls())
#####################
##  Dependencies
source('R/functions-MatModels.R')


###########################
#' 1-locus Pop Gen Invasion
#' Conditions
#' 
#' Parameters:
#' 
popGen_A_invade  <-  function(hf, hm, sm, C) {
	(sm*(C - 1)*(2*hm*(C - 1) - C)) / (sm*(C - 1)*(2*hm*(C - 1) - C) + (C + 1)*(2 - C + 2*hf*(C - 1)))
}
popGen_a_invade  <-  function(hf, hm, sm, C) {
	(sm*(1 - C)*(2 - C+2*hm*(C - 1))) / ((C + 1)*(2*hf*(C - 1) - C)*(sm - 1))
}


###########################
#' Calculate proportion of sf x sm parameter space 
#' that is polymorphic by numerically integrating
#' invasion conditions
#' 
popGen_PolySpace  <-  function(hf, hm, C, sMax) {
	sms        <-  seq(0,sMax, length=10000)
	test_ainv  <-  popGen_a_invade(hf=hf, hm=hm, sm=sms, C=C)
	smCrit     <-  max(sms[test_ainv <= sMax])
	ainv       <-  function(x){(x*(1 - C)*(2 - C+2*hm*(C - 1))) / ((C + 1)*(2*hf*(C - 1) - C)*(x - 1))}
	Ainv       <-  function(x){(x*(C - 1)*(2*hm*(C - 1) - C)) / (x*(C - 1)*(2*hm*(C - 1) - C) + (C + 1)*(2 - C + 2*hf*(C - 1)))}
	part1      <-  integrate(ainv, lower=0, upper=smCrit)$value - integrate(Ainv, lower=0, upper=smCrit)$value
	part2      <-  ((sMax - smCrit)*sMax) - integrate(Ainv, lower=smCrit, upper=sMax)$value
	polySpace  <-  (part1 + part2)/sMax^2
	return(polySpace)
}

###########################################
#' 1-locus Pop Gen Equilibrium Frequencies
#' Conditions
#' 
#' Parameters:
#' 
popGen_qHat_Add  <-  function(sf, sm, C) {
	qHat  <-  (sf*(1 + C) - sm*(1 - sf)*(1 - C)) / (2*sm*sf)
	if(qHat < 0) {
		qHat  <-  0
	}
	if(qHat > 1) {
		qHat  <-  1
	}
	qHat
}
popGen_qHat_DomRev  <-  function(h, sf, sm, C) {
	qHat  <-  (2*sf*(1 + C ) - (sf*(1 + C ) + sm*(1 - C ))*(C + 2*h*(1 - C ))) / (2*(1 - C)*(1 - 2*h)*(sf*(1 + C)+sm*(1 - C)))
	if(qHat < 0) {
		qHat  <-  0
	}
	if(qHat > 1) {
		qHat  <-  1
	}
	qHat
}


###########################################
#' Eigenvalue Calculator for full demographic model
#' 
#(C*(1 - delta)/(1 - C*delta))
#(1 - (C*(1 - delta)/(1 - C*delta)))
calcZeta  <-  function(om, Fii, Fii_pr, USi, UXi, pHat_AA, pHat_aa, C, delta, lambda_AA_aa){
	n_AABound  <-  round(c(C*(1 - delta)*100*c(141.3,18.2,0,0,0,0),(1 - C)*100*c(141.3,18.2,0,0,0,0)))[c(1,2,7,8,3,4,9,10,5,6,11,12)]
	n_aaBound  <-  round(c(C*(1 - delta)*100*c(0,0,0,0,141.3,18.2),(1 - C)*100*c(0,0,0,0,141.3,18.2)))[c(1,2,7,8,3,4,9,10,5,6,11,12)]
	pHat_AA    <-  n_AABound/norm(as.matrix(n_AABound), type="1")
	pHat_aa    <-  n_aaBound/norm(as.matrix(n_aaBound), type="1")
	pn_AA      <-  ones(c(om)) %*% Fii_pr[,,1] %*% (pHat_AA[1:2] + pHat_AA[3:4])
	pn_aa      <-  ones(c(om)) %*% Fii_pr[,,3] %*% (pHat_aa[9:10] + pHat_aa[11:12])
	M22_AA  <-  rbind(
					cbind(USi[,,2] + C*(1 - delta)*Fii[,,2]/2, C*(1 - delta)*Fii[,,2]/2, zeros(c(2,2))),
					cbind(           ((1 - C)/2)*Fii[,,2] + kronecker(c(((1 - C)/(2*pn_AA)))*Fii[,,1]%*%(pHat_AA[1:2] + pHat_AA[3:4]), (ones(om)%*%Fii_pr[,,2])),
						  UXi[,,2] + ((1 - C)/2)*Fii[,,2] + kronecker(c(((1 - C)/(2*pn_AA)))*Fii[,,1]%*%(pHat_AA[1:2] + pHat_AA[3:4]), (ones(om)%*%Fii_pr[,,2])), 
						                 (1 - C)*Fii[,,3] + kronecker(c(((1 - C)/(  pn_AA)))*Fii[,,1]%*%(pHat_AA[1:2] + pHat_AA[3:4]), (ones(om)%*%Fii_pr[,,3]))),
					cbind(((C*(1 - delta))/4)*Fii[,,2], ((C*(1 - delta))/4)*Fii[,,2], USi[,,3] + C*(1 - delta)*Fii[,,3])
					)
	M22_aa  <-  rbind(
					cbind(USi[,,2] + C*(1 - delta)*Fii[,,2]/2, C*(1 - delta)*Fii[,,2]/2, zeros(c(2,2))),
					cbind(           ((1 - C)/2)*Fii[,,2] + kronecker(c(((1 - C)/(2*pn_aa)))*Fii[,,3]%*%(pHat_aa[9:10] + pHat_aa[11:12]), (ones(om)%*%Fii_pr[,,2])),
						  UXi[,,2] + ((1 - C)/2)*Fii[,,2] + kronecker(c(((1 - C)/(2*pn_aa)))*Fii[,,3]%*%(pHat_aa[9:10] + pHat_aa[11:12]), (ones(om)%*%Fii_pr[,,2])),
						                 (1 - C)*Fii[,,1] + kronecker(c(((1 - C)/(  pn_aa)))*Fii[,,3]%*%(pHat_aa[9:10] + pHat_aa[11:12]), (ones(om)%*%Fii_pr[,,1]))),
					cbind(((C*(1 - delta))/4)*Fii[,,2], ((C*(1 - delta))/4)*Fii[,,2], USi[,,1] + C*(1 - delta)*Fii[,,1])
					)
	zeta_i  <-  c(max(eigen(M22_AA, symmetric=FALSE, only.values = TRUE)$values),
				  max(eigen(M22_aa, symmetric=FALSE, only.values = TRUE)$values))
	zeta_i  <-  zeta_i/lambda_AA_aa
	return(zeta_i)
}



#########################
##  Calculate Atilde[ntilde] for 
##  boundary conditions
AtildeBound  <-  function(nBound, g, W_prime, blkFX_prime, blkUS, blkUX, W, Iom, Ig, HS, HX, K, Z, blkFS, blkFX) {

        #Creating the male gamete pool
        nX    <-  nBound[1:6] + nBound[7:12]
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
        AtildeCoexist  <-  Atilde[c(1,2,7,8,3,4,9,10,5,6,11,12),c(1,2,7,8,3,4,9,10,5,6,11,12)]
        return(AtildeCoexist)
}









###########################
#' Forward Simulation 
#'
#' Parameters:
#' dims: 	vector c(om, g), with om = # stages, g = # genotypes
#' theta: 	vector of length 4 (but actually becomes vector of length 7), c(sigmaS_J, sigmaS_A, sigmaX_J, sigmaX_A, gammaS, gammaX, f_ii)
#' theta: 	vector of length 7, c(sigmaS_J, sigmaS_A, gammaS, sigmaX_J, sigmaX_A, gammaX, f_ii)
fwdDemModelSim  <-  function(	om = 2, g = 3, theta = c(0.6, 0.6, 0.05, 6), theta_prime = c(0.6, 0.6, 0.6, 0.6, 0.05, 0.05, 5.9), 
								hf = 1/2, hm = 1/2, sf = 0.1, sm = 0.105, C = 0, delta = 0, 
								delta_j = 0, delta_a = 0, delta_gamma = 0,
								tlimit = 10^4, Ainvade = FALSE, intInit = FALSE, ...) {

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
	# fertility through Female and Male function, depending on whether individuals
	# are produce by selfing or outcrossing
	# theta=[sigmaS_J, sigmaS_A, gammaS, sigmaX_J, sigmaX_A, gammaX, f_ii]
	# Late-acting inbreeding depression
	theta        <-  c(0,0,0, theta)
	theta[1]     <-  theta[4]*(1 - delta_j)
	theta[2]     <-  theta[5]*(1 - delta_a)
	theta[3]     <-  theta[6]*(1 - delta_gamma)
	theta        <-  rep.col(theta,3)
	theta_prime  <-  rep.col(theta_prime,3)
	theta_prime  <-  rep.col(theta_prime,3)

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
	gammaS     <-  theta[3,]
	sigmaX_J   <-  theta[4,]
	sigmaX_A   <-  theta[5,]
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
	# begin generation loop
	while( sum(n) > 1  &&  any(pDelta > 1e-9) && i <= tlimit) {

		nout[,i]  <-  n

	    # Introduce new allele by introducing one heterozygote juvenile
	    if (i==10){
	      n[9]  <-  1
	    }

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

        # Regulate population size to avoid crashing simulation
		if (sum(nnext) > 1e+200) {
			nnext  <-  nnext * 1e-10
			n      <-  nnext
			tmp    <-  kronecker(diag(g),ones(c(1,om))) %*% (nnext[1:6] + nnext[7:12])
			pnext  <-  tmp/colSums(tmp)
			pcheck <-  p
			p      <-  pnext
		} else{
			n      <-  nnext
			tmp    <-  kronecker(diag(g),ones(c(1,om))) %*% (nnext[1:6] + nnext[7:12])
			pnext  <-  tmp/colSums(tmp)
			pcheck <-  p
			p      <-  pnext
		}

		if(i > 20) {
			pDelta  <-  abs(pnext - pcheck)
		}

		#  Calculate eigenvalues based on Analytic results
		if(i == 1) {
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
			testAtilde_AA  <-  AtildeBound(nBound=pBoundAA, g=g, W_prime=W_prime, blkFX_prime = blkFX_prime, blkUS=blkUS, blkUX=blkUX, W=W, Iom=Iom, Ig=Ig, HS=HS, HX=HX, K=K, Z=Z, blkFS=blkFS,blkFX=blkFX)
			testAtilde_aa  <-  AtildeBound(nBound=pBoundaa, g=g, W_prime=W_prime, blkFX_prime=blkFX_prime, blkUS=blkUS, blkUX=blkUX, W=W, Iom=Iom, Ig=Ig, HS=HS, HX=HX, K=K, Z=Z, blkFS=blkFS,blkFX=blkFX)
			
			# Calculate a naive genotype-specific eigenvalue for comparison
			# note, that this gives the same value for lambda_i as we get 
			# elsewhere in the simulation (see lambda_full)
			testLambda_AA  <-  c(max(eigen(testAtilde_AA[c(1:4), c(1:4)],  symmetric=FALSE, only.values = TRUE)$values),
							     max(eigen(testAtilde_AA[c(5:8), c(5:8)],  symmetric=FALSE, only.values = TRUE)$values),
							     max(eigen(testAtilde_AA[c(9:12),c(9:12)], symmetric=FALSE, only.values = TRUE)$values)
							  )
			testLambda_aa  <-  c(max(eigen(testAtilde_aa[c(1:4), c(1:4)],  symmetric=FALSE, only.values = TRUE)$values),
							     max(eigen(testAtilde_aa[c(5:8), c(5:8)],  symmetric=FALSE, only.values = TRUE)$values),
							     max(eigen(testAtilde_aa[c(9:12),c(9:12)], symmetric=FALSE, only.values = TRUE)$values)
							  )
			
			# Calcuate lambda_AA and lambda_aa both using the as described in the notes file
			# (Eq(17) in the linearization at the boundary equilibrium subsection)
			lambda_i   <-  c(testLambda_AA[1], testLambda_aa[3])
#			lambda2_i  <-  c(t(ones(2*om*g)) %*% testAtilde_AA %*% pHat_AA,
#							 t(ones(2*om*g)) %*% testAtilde_aa %*% pHat_aa)

			# Calculate coexistence conditions based on leading eigenvalue of the Jacobian
			zeta_i  <-  calcZeta(om=om, Fii=Fii, Fii_pr=Fii_pr, USi=USi, UXi=UXi,
								 pHat_AA=pHat_AA, pHat_aa=pHat_aa, C=C, delta=delta, 
								 lambda_AA_aa=lambda_i)
		}

		i  <-  i + 1
	}

	##################
	# results
	temp         <-  kronecker(diag(g),ones(c(1,om))) %*% (nout[1:6,] + nout[7:12,])
	p_genotypes  <-  sweep(temp,2,colSums(temp),'/')
	res  <-  list(
					"nout"         =  nout,
					"pout"         =  pout,
					"p_genotypes"  =  p_genotypes,
					"pEq"          =  p_genotypes[,(i-1)],
					"pDelta"       =  pDelta,
					"nzero"        =  nzero,
					"lambda"       =  lambda,
					"lambda_full"  =  lambda_full,
					"zeta_i"       =  zeta_i,
					"om"           =  om,
					"g"            =  g,
					"theta"        =  theta,
					"theta_prime"  =  theta_prime,
					"hf"           =  hf,
					"hm"           =  hm,
					"sf"           =  sf,
					"sm"           =  sm,
					"C"            =  C,
					"delta"        =  delta,
					"delta_j"      =  delta_j,
					"delta_a"      =  delta_a,
					"delta_gamma"  =  delta_gamma,
					"tlimit"       =  tlimit,
					"Ainvade"      =  Ainvade,
					"intInit"      =  intInit,
					"extinct"      =  sum(n) < 1,
					"polymorphism" =  !any(round(p_genotypes[,(i-1)], digits=5) == 1),
					"runtime"      =  (i-1)
					)
	return(res)

}

##############################
#' Loop over selection space 
#'
#' Parameters:
#' dims: 	vector c(om, g), with om = # stages, g = # genotypes
#' theta: 	vector of length 4, c(sigma_J, sigma_A, gamma, f_ii)
#' selPars:	vector of length 4, c(hf, hm, sf, sm,)
selLoop  <-  function(sMax = 0.15, nSamples=1e+2,
					  om = 2, g = 3, theta = c(0.6,0.6,0.05,6.1), theta_prime = c(0.6,0.6,0.05,6.1), 
					  hf = 1/2, hm = 1/2, C = 0, delta = 0, 
					  delta_j = 0, delta_a = 0, delta_gamma = 0,
					  tlimit = 10^5, intInit = FALSE, writeFile=TRUE, progBar = TRUE, ...) {
	
	sfs           <-  runif(min = 0, max=sMax, n=nSamples)
	sms           <-  runif(min = 0, max=sMax, n=nSamples)
	extinct       <-  rep(NA,times=nSamples)
	polymorphism  <-  rep(NA,times=nSamples)
	pEq           <-  matrix(0, ncol=3, nrow=nSamples)
	zeta_i        <-  matrix(0, ncol=2, nrow=nSamples)
	lambda_i      <-  matrix(0, ncol=3, nrow=nSamples)

	for(i in 1:nSamples) {
		if(hf == hm && hf == 1/2) {
			qHat  <-  popGen_qHat_Add(sf = sfs[i], sm = sms[i], C = C)
		}
		if(hf == hm && hf < 1/2) {
			qHat  <-  popGen_qHat_DomRev(h = hf, sf = sfs[i], sm = sms[i], C = C)
		}
			if(qHat < 1/2){
				Ainvade  <-  TRUE
			}
			if(qHat > 1/2){
				Ainvade  <-  FALSE
			}

			results  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
										hf = hf, hm = hm, sf = sfs[i], sm = sms[i], C = C, delta = delta, 
										delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
										tlimit = 10^5, Ainvade = Ainvade, intInit = intInit)
			intInit  <-  FALSE
			extinct[i]       <-  results$extinct
			polymorphism[i]  <-  results$polymorphism
			pEq[i,]          <-  results$pEq
			zeta_i[i,]       <-  results$zeta_i
			lambda_i[i,]     <-  results$lambda_full

		if(progBar){
			cat('\r', paste(100*(i/nSamples),'% Complete'))
		}

	}

	# compile results as data frame
	hfVec       <-  rep(hf, times=nSamples)
	hmVec       <-  rep(hm, times=nSamples)
	CVec        <-  rep(C, times=nSamples)
	deltaVec    <-  rep(delta, times=nSamples)
	results.df  <-  as.data.frame(cbind(sfs, 
										sms, 
										extinct, 
										polymorphism, 
										pEq, 
										zeta_i,
										lambda_i,
										hfVec, 
										hmVec, 
										CVec, 
										delta
										)
								 )
	colnames(results.df)  <-  c("sf",
								"sm",
								"extinct",
								"poly",
								"Eq_pAA",
								"Eqp_Aa",
								"Eq_paa",
								"zeta_AA",
								"zeta_aa",
								"lambda_AA",
								"lambda_Aa",
								"lambda_aa",
								"hf",
								"hm",
								"C",
								"delta"
								)

	# export data as .csv to ./output/data
	if(writeFile) {
			filename <-  paste("./output/simData/demSimsSfxSm", "_sMax", sMax, "_nSamples", nSamples, "_hf", hf, "_hm", hm, 
							"_C", C, "_delta", delta, "_dj", delta_j, "_da", delta_a, "_dg", delta_gamma, ".csv", sep="")
			write.csv(results.df, file=filename, row.names = FALSE)
	} else{
			return(results.df)
	}

}






##############################
#' Quantify proportion of polymorphic AND/OR demographically viable
#' Parameter space for different Selfing rates & dominance
#'
polyParamSpaceMakeData  <-  function(sMax = 0.15, nSamples=1e+2,
									 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = c(0.6,0.6,0.05,6), 
									 hf = 1/2, hm = 1/2, C = 0, delta = 0, 
									 delta_j = 0, delta_a = 0, delta_gamma = 0,
									 tlimit = 10^5) {
	
	Cs             <-  seq(0,0.9, by=0.1)
	extinct        <-  c()
	popGenPoly     <-  c()
	simPoly        <-  c()
	simPolyViable  <-  c()
	simAFix        <-  c()
	simAFixViable  <-  c()
	simaFix        <-  c()
	simaFixViable  <-  c()
	eigPoly        <-  c()
	eigPolyViable  <-  c()
	eigAFix        <-  c()
	eigAFixViable  <-  c()
	eigaFix        <-  c()
	eigaFixViable  <-  c()
	for(i in 1:length(Cs)) {
		res  <-  selLoop(sMax=sMax, nSamples=nSamples, 
						 om = om, g = g, theta = theta, theta_prime = theta_prime, 
						 hf = hf, hm = hm, C = Cs[i], delta = delta,
						 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
						 tlimit = tlimit, intInit = FALSE, writeFile = FALSE, progBar = FALSE)
	
		# Calculate & store results for plotting
		extinct[i]        <-  sum(res$extinct)/nSamples
		popGenPoly[i]     <-  popGen_PolySpace(hf=hf, hm=hm, C=Cs[i], sMax=sMax)
		simPoly[i]        <-  sum(res$poly)/nSamples
		simPolyViable[i]  <-  simPoly[i] - sum(res$extinct == 1 & res$poly == 1)/nSamples
		simAFix[i]        <-  sum(res$poly == 0 & round(res$Eq_pAA) == 1)/nSamples
		simAFixViable[i]  <-  simAFix[i] - sum(res$extinct[res$poly == 0 & round(res$Eq_pAA) == 1])/nSamples
		simaFix[i]        <-  sum(res$poly == 0 & round(res$Eq_paa) == 1)/nSamples
		simaFixViable[i]  <-  simaFix[i] - sum(res$extinct[res$poly == 0 & round(res$Eq_paa) == 1])/nSamples
		eigPoly [i]       <-  sum(res$zeta_AA > 1 & res$zeta_aa > 1)/nSamples
		eigPolyViable[i]  <-  eigPoly[i] - sum(res$extinct == 1 & res$zeta_AA > 1 & res$zeta_aa > 1)/nSamples
		eigAFix[i]        <-  sum(res$zeta_AA > 1 & res$zeta_aa < 1)/nSamples
		eigAFixViable[i]  <-  eigAFix[i] - sum(res$extinct[res$poly == 0 & res$zeta_AA > 1 & res$zeta_aa < 1])/nSamples
		eigaFix[i]        <-  sum(res$zeta_AA < 1 & res$zeta_aa > 1)/nSamples
		eigaFixViable[i]  <-  eigaFix[i] - sum(res$extinct[res$poly == 0 & res$zeta_AA < 1 & res$zeta_aa > 1])/nSamples
	
		cat('\r', paste(100*(i/length(Cs)),'% Complete'))
	}

	# compile results as data frame
	results.df  <-  as.data.frame(cbind(Cs, 
										extinct,
										popGenPoly,
										simPoly,
										simPolyViable,
										simAFix,
										simAFixViable,
										simaFix,
										simaFixViable,
										eigPoly,
										eigPolyViable,
										eigAFix,
										eigAFixViable,
										eigaFix,
										eigaFixViable
										)
								 )

	# Export data as .csv to ./output/data
	filename <-  paste("./output/simData/simPolySpace", "_sMax", sMax, "_nSamples", nSamples, "_hf", hf, "_hm", hm, 
					   "_delta", delta, "_dj", delta_j, "_da", delta_a, "_dg", delta_gamma, "_f", theta[4], ".csv", sep="")
	write.csv(results.df, file=filename, row.names = FALSE)

}


