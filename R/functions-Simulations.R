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

###########################
#' Forward Simulation 
#'
#' Parameters:
#' dims: 	vector c(om, g), with om = # stages, g = # genotypes
#' theta: 	vector of length 7, c(sigmaS_J, sigmaS_A, sigmaX_J, sigmaX_A, gammaS, gammaX, f_ii)
fwdDemModelSim  <-  function(	om = 2, g = 3, theta = c(0.6, 0.6, 0.6, 0.6, 0.05, 0.05, 5.9), theta_prime = c(0.6, 0.6, 0.6, 0.6, 0.05, 0.05, 5.9), 
								hf = 1/2, hm = 1/2, sf = 0.1, sm = 0.105, C = 0, delta = 0, tlimit = 10^4, Ainvade = FALSE, intInit = FALSE, ...) {

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
	FSi        <-  zeros(c(om,om,g))
	FXi        <-  zeros(c(om,om,g))
	FXi_prime  <-  zeros(c(om,om,g))
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
	    FSi[,,i]       <- rbind(c(0,C*f[i]*(1 - delta)),
						        c(0,0))
#	    FSi[,,i]       <- rbind(c(0,(C*(1 - delta)/(1 - C*delta))*f[i]),
#						        c(0,0))
	    FXi[,,i]       <- rbind(c(0,(1 - (C*(1 - delta)/(1 - C*delta)))*f[i]),
						        c(0,0))
	    FXi_prime[,,i] <- rbind(c(0,(1 - (C*(1 - delta)/(1 - C*delta)))*f_prime[i]),
								c(0,0))
    
		# genotype-specific population growth rate
		eigs_temp[i,]  <-  eigen((FSi[,,i] + USi[,,i]) + (FXi[,,i] + UXi[,,i]), symmetric=FALSE, only.values = TRUE)$values
		lambda[i]      <-  max(eigs_temp)
		Nfun[,,i]      <-  solve((diag(2) - (C*(1 - delta)*USi[,,i] + (1 - C)*UXi[,,i])), diag(2))
#		Nfun[,,i]      <-  solve((diag(2) - ((C*(1 - delta)/(1 - C*delta))*USi[,,i] + (1 - (C*(1 - delta)/(1 - C*delta)))*UXi[,,i])), diag(2))
		R[,,i]         <-  (FSi[,,i] + FXi[,,i]) %*% Nfun[,,i]
		#genotype-specific R0
		R_0[i] <- max(eigen(R[,,i],symmetric=FALSE, only.values = TRUE)$values)
	}
browser()
	Atilde_genotype  <-  zeros(c(2*om,2*om,g))
	lambda_full      <-  zeros(c(1,g))
	for (i in 1:3){
		Atilde_genotype[,,i]  <- rbind( cbind(USi[,,i] + FSi[,,i], zeros(c(2,2))), 
										cbind(FXi[,,i], UXi[,,i])
								   	   )
		lambda_full[i]    <- max(eigen(Atilde_genotype[,,i],symmetric=FALSE, only.values = TRUE)$values)
	}

# Coexistence Conditions Based on Leading Eigenvalue of the Jacobian
	n_AABound  <-  round(c(C*(1 - delta)*100*c(141.3,18.2,0,0,0,0),(1 - C)*100*c(141.3,18.2,0,0,0,0)))[c(1,2,7,8,3,4,9,10,5,6,11,12)]
	n_aaBound  <-  round(c(C*(1 - delta)*100*c(0,0,0,0,141.3,18.2),(1 - C)*100*c(0,0,0,0,141.3,18.2)))[c(1,2,7,8,3,4,9,10,5,6,11,12)]

	pHat_AA         <-  n_AABound/c(sum(n_AABound[1:2]), sum(n_AABound[1:2]), sum(n_AABound[3:4]), sum(n_AABound[3:4]), 1,1,1,1,1,1,1,1)
	pHat_aa         <-  n_aaBound/c(1,1,1,1,1,1,1,1, sum(n_aaBound[9:10]), sum(n_aaBound[9:10]), sum(n_aaBound[11:12]), sum(n_aaBound[11:12]))
	pn_AA           <-  ones(c(om))%*%FXi_prime[,,1]%*%(pHat_AA[1:2] + pHat_AA[3:4])
	pn_aa           <-  ones(c(om))%*%FXi_prime[,,3]%*%(pHat_aa[9:10] + pHat_aa[11:12])
	M22_AA  <-  rbind(
					cbind(USi[,,2] + C*(1 - delta)*FSi[,,2]/2, C*(1 - delta)*FSi[,,2]/2, zeros(c(2,2))),
					cbind(((1 - C)/2)*FSi[,,2] + c(((1 - C)/(2*pn_AA)))*FSi[,,1]%*%kronecker((pHat_AA[1:2] + pHat_AA[3:4]), (ones(om)%*%FXi_prime[,,2])),
						  UXi[,,2] + ((1 - C)/2)*FSi[,,2] + c(((1 - C)/(2*pn_AA)))*FSi[,,1]%*%kronecker((pHat_AA[1:2] + pHat_AA[3:4]), (ones(om)%*%FXi_prime[,,2])), c((1 - C))*FSi[,,3]),
					cbind(((C*(1 - delta))/4)*FSi[,,2], ((C*(1 - delta))/4)*FSi[,,2], USi[,,3] + C*(1 - delta)*FSi[,,3])
					)
	M22_aa  <-  rbind(
					cbind(USi[,,2] + C*(1 - delta)*FSi[,,2]/2, C*(1 - delta)*FSi[,,2]/2, zeros(c(2,2))),
					cbind(((1 - C)/2)*FSi[,,2] + c(((1 - C)/(2*pn_aa)))*FSi[,,3]%*%kronecker((pHat_aa[1:2] + pHat_aa[3:4]), (ones(om)%*%FXi_prime[,,2])),
						  UXi[,,2] + ((1 - C)/2)*FSi[,,2] + c(((1 - C)/(2*pn_aa)))*FSi[,,3]%*%kronecker((pHat_aa[1:2] + pHat_AA[3:4]), (ones(om)%*%FXi_prime[,,2])), c((1 - C))*FSi[,,1]),
					cbind(((C*(1 - delta))/4)*FSi[,,2], ((C*(1 - delta))/4)*FSi[,,2], USi[,,1] + C*(1 - delta)*FSi[,,1])
					)
	max(eigen(M22_AA, symmetric=FALSE, only.values = TRUE)$values)
	max(eigen(M22_aa, symmetric=FALSE, only.values = TRUE)$values)


	# CREATE BLOCK DIAGONAL MATRICES
	d 			 <- diag(g)
	blkD         <-  kronecker(Iom,d)
	blkUS        <-  zeros(c(g*om,g*om))
	blkUX        <-  zeros(c(g*om,g*om))
	blkFS        <-  zeros(c(g*om,g*om))
	blkFX        <-  zeros(c(g*om,g*om))
	blkFX_prime  <-  zeros(c(g*om,g*om))
	for (i in 1:g){
	    blkUS        <-  blkUS + kronecker(emat(g,g,i,i),USi[,,i])
	    blkUX        <-  blkUX + kronecker(emat(g,g,i,i),UXi[,,i])
	    blkFS        <-  blkFS + kronecker(emat(g,g,i,i),FSi[,,i])
	    blkFX        <-  blkFX + kronecker(emat(g,g,i,i),FXi[,,i])
	    blkFX_prime  <-  blkFX_prime + kronecker(emat(g,g,i,i),FXi_prime[,,i])
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
#	for (i in 1:tlimit){

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
            pi  <-  Ig[,ii]
            qi  <-  W %*% pi #allele frequencies in oocytes of genotype i

            # genotype frequencies in the offspring of mothers of genotype i
            # produced by outcrossing (equation 16 in de Vries and Caswell, 2018a (American Naturalist))
            piprime  <-  Z %*% kronecker(qi,q_prime)

            HX[,ii]  <-  piprime # the outcrossing parent-offspring matrix
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
#		cat('\r', i/tlimit,'%')
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
					  hf = 1/2, hm = 1/2, C = 0, delta = 0, tlimit = 10^5, intInit = FALSE, returnRes=FALSE) {
	
	sfs           <-  runif(min = 0, max=sMax, n=nSamples)
	sms           <-  runif(min = 0, max=sMax, n=nSamples)
	extinct       <-  rep(NA,times=nSamples)
	polymorphism  <-  rep(NA,times=nSamples)
	pEq           <-  matrix(0, ncol=3, nrow=nSamples)

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
							hf = hf, hm = hm, sf = sfs[i], sm = sms[i], C = C, delta = delta, tlimit = 10^5, 
							Ainvade = Ainvade, intInit = intInit)
			intInit  <-  FALSE
			extinct[i]       <-  results$extinct
			polymorphism[i]  <-  results$polymorphism
			pEq[i,]          <-  results$pEq

		cat('\r', paste(100*(i/nSamples),'% Complete'))

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
								"hf",
								"hm",
								"C",
								"delta"
								)

	# export data as .csv to ./output/data
	filename <-  paste("./output/simData/demSimsSfxSm_MinorEqAlleleInv", "_sMax", sMax, "_nSamples", nSamples, "_hf", hf, "_hm", hm, "_C", C, "_delta", delta, ".csv", sep="")
	write.csv(results.df, file=filename, row.names = FALSE)

	if(returnRes) {
		return(results.df)
	}
}


