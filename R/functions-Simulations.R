################################################################
#  Functions to simulate demographic model across selection
#  parameter space
#
#
#  Author: Colin Olito, adapted from L. DeVries
#
#  NOTES:  
#		


#rm(list=ls())
#####################
##  Dependencies
source('R/functions-MatModels.R')


###########################
#' 1-locus Pop Gen Invasion
#' Conditions
#' 
popGen_A_invade  <-  function(hf, hm, sm, C) {
	(sm*(C - 1)*(2*hm*(C - 1) - C)) / (sm*(C - 1)*(2*hm*(C - 1) - C) + (C + 1)*(2 - C + 2*hf*(C - 1)))
}
popGen_a_invade  <-  function(hf, hm, sm, C) {
	(sm*(1 - C)*(2 - C+2*hm*(C - 1))) / ((C + 1)*(2*hf*(C - 1) - C)*(sm - 1))
}

popGen_A_invade_Delta_Add  <-  function(C, delta, sm) {
	(sm*(1 - C)) / (1 + sm*(1 - C) + C*(1 - 2*delta))
}
popGen_a_invade_Delta_Add  <-  function(C, delta, sm) {
	(sm*(1 - C)) / ((1 - sm)*(1 + C*(1 - 2*delta)))
}

popGen_a_invade_Delta_DomRev  <-  function(C, delta, sm) {
	(2 - 2*sm + C*(4 - 8*delta + sm*(-9 + 13*delta)) + (C^2)*(2*(1 - 2*delta)^2 + 
	sm*(3 + 3*delta - 8*delta^2)) - sqrt(4*(-1 + sm)^2 - 4*C*(-1 + sm)*(4 - 8*delta + sm*(-3 + 7*delta)) + 
	(C^2)*(24*(1 - 2*delta)^2 - 4*sm*(19 - 59*delta + 52*delta^2) + (sm^2)*(77 - 190*delta + 137*delta^2)) - 
	2*(C^3)*(8*(-1 + 2*delta)^3 - 2*sm*(-9 + 65*delta - 124*delta^2 + 76*delta^3) + (sm^2)*(15 + 32*delta - 127*delta^2 + 88*delta^3)) + 
	(C^4)*(4*(1 - 2*delta)^4 - 4*sm*(-5 + 11*delta + 16*delta^2 - 52*delta^3 + 32*delta^4) + 
	(sm^2)*(1 + 26*delta - 7*delta^2 - 80*delta^3 + 64*delta^4)))) / (2*C*(-1 + sm)*(-1 + delta)*(-1 + C*(-1 + 2*delta)))
}
popGen_A_invade_Delta_DomRev  <-  function(C, delta, sm) {
	-((-6 - 4*C + 2*(C^2) - 4*sm + C*sm + 3*(C^2)*sm + 16*C*delta + 7*C*sm*delta - 7*(C^2)*sm*delta - 
	8*(C^2)*delta^2 + sqrt((-6 - 4*sm + C*(-4 + sm + 16*delta + 7*sm*delta) + (C^2)*(2 + sm*(3 - 7*delta) - 8*delta^2))^2 - 
	8*(-1 + C)*sm*(-1 + C*(-1 + 2*delta))*(2*(3 + sm) - C*(-7 + sm + 19*delta + 3*sm*delta) + 
	(C^2)*(1 - 9*delta + 14*delta^2 + sm*(-1 + 3*delta))))) / (4*(3 + sm) - 2*C*(-7 + sm + 19*delta + 3*sm*delta) + 
	2*(C^2)*(1 - 9*delta + 14*delta^2 + sm*(-1 + 3*delta))))
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

popGen_PolySpace_Delta_Add  <-  function(C, delta, sMax) {
	sms        <-  seq(0,sMax, length=10000)
	test_ainv  <-  popGen_a_invade_Delta_Add(C=C, delta = delta, sm=sms)
	smCrit     <-  max(sms[test_ainv <= sMax])
	ainv       <-  function(x){(x*(1 - C)) / ((1 - x)*(1 + C*(1 - 2*delta)))}
	Ainv       <-  function(x){(x*(1 - C)) / (1 + x*(1 - C) + C*(1 - 2*delta))}
	part1      <-  integrate(ainv, lower=0, upper=smCrit)$value - integrate(Ainv, lower=0, upper=smCrit)$value
	part2      <-  ((sMax - smCrit)*sMax) - integrate(Ainv, lower=smCrit, upper=sMax)$value
	polySpace  <-  (part1 + part2)/sMax^2
	return(polySpace)
}

popGen_PolySpace_Delta_DomRev  <-  function(C, delta, sMax) {
	sms        <-  seq(0,sMax, length=10000)
	test_ainv  <-  popGen_a_invade_Delta_DomRev(C=C, delta=delta, sm=sms)
	smCrit     <-  max(sms[test_ainv <= sMax])
	Ainv       <-  function(x){-((-6 - 4*C + 2*(C^2) - 4*x + C*x + 3*(C^2)*x + 16*C*delta + 7*C*x*delta - 7*(C^2)*x*delta - 
	8*(C^2)*delta^2 + sqrt((-6 - 4*x + C*(-4 + x + 16*delta + 7*x*delta) + (C^2)*(2 + x*(3 - 7*delta) - 8*delta^2))^2 - 
	8*(-1 + C)*x*(-1 + C*(-1 + 2*delta))*(2*(3 + x) - C*(-7 + x + 19*delta + 3*x*delta) + 
	(C^2)*(1 - 9*delta + 14*delta^2 + x*(-1 + 3*delta))))) / (4*(3 + x) - 2*C*(-7 + x + 19*delta + 3*x*delta) + 
	2*(C^2)*(1 - 9*delta + 14*delta^2 + x*(-1 + 3*delta))))}
	ainv       <-  function(x){(2 - 2*x + C*(4 - 8*delta + x*(-9 + 13*delta)) + (C^2)*(2*(1 - 2*delta)^2 + 
	x*(3 + 3*delta - 8*delta^2)) - sqrt(4*(-1 + x)^2 - 4*C*(-1 + x)*(4 - 8*delta + x*(-3 + 7*delta)) + 
	(C^2)*(24*(1 - 2*delta)^2 - 4*x*(19 - 59*delta + 52*delta^2) + (x^2)*(77 - 190*delta + 137*delta^2)) - 
	2*(C^3)*(8*(-1 + 2*delta)^3 - 2*x*(-9 + 65*delta - 124*delta^2 + 76*delta^3) + (x^2)*(15 + 32*delta - 127*delta^2 + 88*delta^3)) + 
	(C^4)*(4*(1 - 2*delta)^4 - 4*x*(-5 + 11*delta + 16*delta^2 - 52*delta^3 + 32*delta^4) + 
	(x^2)*(1 + 26*delta - 7*delta^2 - 80*delta^3 + 64*delta^4)))) / (2*C*(-1 + x)*(-1 + delta)*(-1 + C*(-1 + 2*delta)))}
	part1      <-  integrate(ainv, lower=0, upper=smCrit)$value - integrate(Ainv, lower=0, upper=smCrit)$value
	part2      <-  ((sMax - smCrit)*sMax) - integrate(Ainv, lower=smCrit, upper=sMax)$value
	polySpace  <-  (part1 + part2)/sMax^2
	return(polySpace)
}


###########################################
#' 1-locus Pop Gen Equilibrium Frequencies
#' Conditions
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

#############################################
#' Heuristic inbreeding depression functions
#' (see Olito & Connallon 2019, Appendix E)
#' 
#' Parameters
#' a: 	shape parameter determiing curvature of the line. 
#' 		Linear when a = 1, concave up when a < 1, concave
#' 		downward when a > 1
#' C:	Selfing rate
Load  <-  function(a, C) {
	a*(1 - C)/(C + a*(1 - C))
}
#' Parameters
#' dStar:	Hypothetical severity of inbreeding depression in 
#' 			obligately outcrossing population 
#' b:		shape parameter determining how far delta will 
#' 			decline under complete selfing (C = 1)
#' a:		As above.
#' C:		As above.
predDelta  <-  function(dStar, b=1/2, a=0.2, C) {
	dStar - dStar*b*(1 - Load(a = a, C = C))
}

###########################################
#' Eigenvalue Calculator for full demographic model
#' 
#calcZeta  <-  function(om, FSi, FXi, FXi_pr, USi, UXi, pHat_AA, pHat_aa, C, delta, lambda_i){
#
#	pn_AA   <-  ones(c(om)) %*% FXi_pr[,,1] %*% (pHat_AA[1:om] + pHat_AA[om+1:om])
#	pn_aa   <-  ones(c(om)) %*% FXi_pr[,,3] %*% (pHat_aa[(4*om+1):(5*om)] + pHat_aa[(5*om+1):(6*om)])
#	M22_AA  <-  rbind(
#					cbind(USi[,,2] + C*(1 - delta)*FSi[,,2]/2, C*(1 - delta)*FSi[,,2]/2, zeros(c(om,om))),
#					cbind(           ((1 - C)/2)*FXi[,,2] + kronecker(c(((1 - C)/(2*pn_AA)))*FXi[,,1]%*%(pHat_AA[1:om] + pHat_AA[om+1:om]), (ones(om)%*%FXi_pr[,,2])),
#						  UXi[,,2] + ((1 - C)/2)*FXi[,,2] + kronecker(c(((1 - C)/(2*pn_AA)))*FXi[,,1]%*%(pHat_AA[1:om] + pHat_AA[om+1:om]), (ones(om)%*%FXi_pr[,,2])), 
#						                 (1 - C)*FXi[,,3] + kronecker(c(((1 - C)/(  pn_AA)))*FXi[,,1]%*%(pHat_AA[1:om] + pHat_AA[om+1:om]), (ones(om)%*%FXi_pr[,,3]))),
#					cbind(((C*(1 - delta))/4)*FSi[,,2], ((C*(1 - delta))/4)*FSi[,,2], USi[,,3] + C*(1 - delta)*FSi[,,3])
#					)
#	M22_aa  <-  rbind(
#					cbind(USi[,,2] + C*(1 - delta)*FSi[,,2]/2, C*(1 - delta)*FSi[,,2]/2, zeros(c(om,om))),
#					cbind(           ((1 - C)/2)*FXi[,,2] + kronecker(c(((1 - C)/(2*pn_aa)))*FXi[,,3]%*%(pHat_aa[(4*om+1):(5*om)] + pHat_aa[(5*om+1):(6*om)]), (ones(om)%*%FXi_pr[,,2])),
#						  UXi[,,2] + ((1 - C)/2)*FXi[,,2] + kronecker(c(((1 - C)/(2*pn_aa)))*FXi[,,3]%*%(pHat_aa[(4*om+1):(5*om)] + pHat_aa[(5*om+1):(6*om)]), (ones(om)%*%FXi_pr[,,2])),
#						                 (1 - C)*FXi[,,1] + kronecker(c(((1 - C)/(  pn_aa)))*FXi[,,3]%*%(pHat_aa[(4*om+1):(5*om)] + pHat_aa[(5*om+1):(6*om)]), (ones(om)%*%FXi_pr[,,1]))),
#					cbind(((C*(1 - delta))/4)*FSi[,,2], ((C*(1 - delta))/4)*FSi[,,2], USi[,,1] + C*(1 - delta)*FSi[,,1])
#					)
#	zeta_i  <-  c(max(Re(eigen(M22_AA, symmetric=FALSE, only.values = TRUE)$values)),
#				  max(Re(eigen(M22_aa, symmetric=FALSE, only.values = TRUE)$values)))
#	zeta_i  <-  zeta_i/lambda_i[c(1,3)]
#	return(zeta_i)
#}



calcZeta  <-  function(om, FSi, FXi, FXi_pr, USi, UXi, pHat_AA, pHat_aa, C, delta, lambda_i){

	pn_AA   <-  ones(c(om)) %*% FXi_pr[,,1] %*% (pHat_AA[1:om] + pHat_AA[om+1:om])
	pn_aa   <-  ones(c(om)) %*% FXi_pr[,,3] %*% (pHat_aa[(4*om+1):(5*om)] + pHat_aa[(5*om+1):(6*om)])
	M22_AA  <-  rbind(
					cbind(USi[,,2] + FSi[,,2]/2, FSi[,,2]/2, zeros(c(om,om))),
					cbind(           (1/2)*FXi[,,2] + kronecker(c((1/(2*pn_AA)))*FXi[,,1]%*%(pHat_AA[1:om] + pHat_AA[om+1:om]), (ones(om)%*%FXi_pr[,,2])),
						  UXi[,,2] + (1/2)*FXi[,,2] + kronecker(c((1/(2*pn_AA)))*FXi[,,1]%*%(pHat_AA[1:om] + pHat_AA[om+1:om]), (ones(om)%*%FXi_pr[,,2])), 
						                   FXi[,,3] + kronecker(c((1/(  pn_AA)))*FXi[,,1]%*%(pHat_AA[1:om] + pHat_AA[om+1:om]), (ones(om)%*%FXi_pr[,,3]))),
					cbind((1/4)*FSi[,,2], (1/4)*FSi[,,2], USi[,,3] + FSi[,,3])
					)
	M22_aa  <-  rbind(
					cbind(USi[,,2] + FSi[,,2]/2, FSi[,,2]/2, zeros(c(om,om))),
					cbind(           (1/2)*FXi[,,2] + kronecker(c((1/(2*pn_aa)))*FXi[,,3]%*%(pHat_aa[(4*om+1):(5*om)] + pHat_aa[(5*om+1):(6*om)]), (ones(om)%*%FXi_pr[,,2])),
						  UXi[,,2] + (1/2)*FXi[,,2] + kronecker(c((1/(2*pn_aa)))*FXi[,,3]%*%(pHat_aa[(4*om+1):(5*om)] + pHat_aa[(5*om+1):(6*om)]), (ones(om)%*%FXi_pr[,,2])),
						                   FXi[,,1] + kronecker(c((1/(  pn_aa)))*FXi[,,3]%*%(pHat_aa[(4*om+1):(5*om)] + pHat_aa[(5*om+1):(6*om)]), (ones(om)%*%FXi_pr[,,1]))),
					cbind((1/4)*FSi[,,2], (1/4)*FSi[,,2], USi[,,1] + FSi[,,1])
					)
	zeta_i  <-  c(max(Re(eigen(M22_AA, symmetric=FALSE, only.values = TRUE)$values)),
				  max(Re(eigen(M22_aa, symmetric=FALSE, only.values = TRUE)$values)))
	zeta_i  <-  zeta_i/lambda_i[c(1,3)]
	return(zeta_i)
}



#########################
##  Evaluate Atilde[ntilde] for 
##  a given initial equilibrium
evalAtilde  <-  function(nHat, om, g, W_prime, blkFX_prime, blkUS, blkUX, W, Iom, Ig, HS, HX, K, Z, blkFS, blkFX, reorder=TRUE) {

        #Creating the male gamete pool
        nX    <-  nHat[1:g*om] + nHat[(g*om+1):(2*g*om)]
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

        if(reorder == FALSE) {
        	AtildeEval  <-  Atilde	
        } else{
        	AtildeEval  <-  Atilde[c(1:om, g*om+1:om, om+1:om, g*om+om+1:om, 2*om+1:om, g*om+2*om+1:om),
								   c(1:om, g*om+1:om, om+1:om, g*om+om+1:om, 2*om+1:om, g*om+2*om+1:om)]
		}
        return(AtildeEval)
}




#############################
#' Forward dynamics to Equil 
#'
#' Parameters:
#' om:		Number of stages in life-cycle
#' g:		number of genotypes (default 3 for 1-locus, 2-allele)
#' theta: 	vector of length 4 (but actually becomes vector of length 7), c(sigmaS_J, sigmaS_A, sigmaX_J, sigmaX_A, gammaS, gammaX, f_ii)
#' theta_prime: Pollen production (specific value is not super important, default is set to equal f)
#' hf,hm:	Dominance coefficient for SA fitness effects
#' sf,sm:	Selection coefficient for SA fitness effects
#' C:		Population selfing rate
#' delta:	Early-acting inbreeding depression (proportion ovules aborted due to I.D.)
#' delta_j:	Proportional decrease in juvenile survival rate for inbred offspring
#' delta_j:	Proportional decrease in adult survival rate for inbred offspring
#' delta_gamma: Proportional decrease in juv. --> adult transition rate for inbred offspring
#' tlimit:	Max # of generations for fwd simulation
#' Ainvade:	Should A allele start off rare
#' intInit:	Optional initial frequency of A allele
fwdDyn2Eq  <-  function(nzero, om, g, W_prime, blkFX_prime, blkUS, blkUX, W, Iom, Ig, K, Z, blkFS, blkFX, tlimit, eqThreshold, ...){

	# initial state vector, frequency vector, and empty storage matrix for time-series
	n     <-  t(t(nzero))
	p     <-  n/norm(as.matrix(n), type="1")
	tmp   <-  kronecker(diag(g),ones(c(1,om))) %*% (nzero[1:(om*g)] + nzero[(om*g)+1:(om*g)])
	p_g   <-  tmp/colSums(tmp)
	nout  <-  zeros(c((2*om*g),tlimit))
	
	# set arbitrary value for Euclidian distance 
	# used to check if eq. has been reached 
	pDist  <-  1
	
	# begin generation loop
	i  <-  1
	while( pDist > eqThreshold && i <= tlimit) {

		# Populate n timeseries for current generation
		nout[,i]  <-  n

        #Creating the male gamete pool
        nX    <-  n[1:(g*om)] + n[(g*om+1):(2*g*om)]
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
		if (sum(nnext) > 1e+100) {
			nnext  <-  nnext * 1e-20
			n      <-  nnext
			pnext  <-  n/norm(as.matrix(n), type="1")
			pcheck <-  p
			p      <-  pnext
		} else{
			n      <-  nnext
			pnext  <-  n/norm(as.matrix(n), type="1")
			pcheck <-  p
			p      <-  pnext
		}

		# check to see if equilibrium has been reached
		pDist  <-  eucDist(pnext, pcheck)

		# next generation
		i  <-  i + 1
	}

	# clean up time-series
	nout         <-  nout[,1:(i-1)]
	temp         <-  kronecker(diag(g),ones(c(1,om))) %*% (nout[1:(3*om),] + nout[(3*om+1):(6*om),])
	p_genotypes  <-  sweep(temp,2,colSums(temp),'/')

	# Return results
	results  <-  list(
					  "nout"         =  nout,
					  "n"            =  nout[,ncol(nout)],
					  "lambda_sim"   =  sum(nout[,ncol(nout)]) / sum(nout[,ncol(nout)-1]),
					  "p"            =  p,
					  "p_genotypes"  =  p_genotypes,
					  "eqReached"    =  pDist <= eqThreshold
					)
	return(results)
}




###########################
#' Forward Simulation 
#'
#' Parameters:
#' om:		Number of stages in life-cycle
#' g:		number of genotypes (default 3 for 1-locus, 2-allele)
#' theta: 	vector of length 4 (but actually becomes vector of length 7), c(sigmaS_J, sigmaS_A, sigmaX_J, sigmaX_A, gammaS, gammaX, f_ii)
#' theta_prime: Pollen production (specific value is not super important, default is set to equal f)
#' hf,hm:	Dominance coefficient for SA fitness effects
#' sf,sm:	Selection coefficient for SA fitness effects
#' C:		Population selfing rate
#' delta:	Early-acting inbreeding depression (proportion ovules aborted due to I.D.)
#' delta_j:	Proportional decrease in juvenile survival rate for inbred offspring
#' delta_j:	Proportional decrease in adult survival rate for inbred offspring
#' delta_gamma: Proportional decrease in juv. --> adult transition rate for inbred offspring
#' tlimit:	Max # of generations for fwd simulation
#' Ainvade:	Should A allele start off rare
#' intInit:	Optional initial frequency of A allele
fwdDemModelSim  <-  function(om = 2, g = 3, theta = c(0.6, 0.6, 0.05, 6), theta_prime = 5.9, 
							 hf = 1/2, hm = 1/2, sf = 0.1, sm = 0.105, C = 0, delta = 0, 
							 delta_j = 0, delta_a = 0, delta_gamma = 0, datMat = NA,
							 tlimit = 10^4, eqThreshold = 1e-8, Ainvade = FALSE, intInit = FALSE, ...) {

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
	theta     <-  c(0,0,0, theta)
	theta[1]  <-  theta[4]*(1 - delta_j)
	theta[2]  <-  theta[5]*(1 - delta_a)
	theta[3]  <-  theta[6]*(1 - delta_gamma)
	theta     <-  rep.col(theta,3)
	f_prime   <-  rep(theta_prime,3)	
	
	# Fitness Expressions
	fii        <-  c(1, 1 - hf*sf, 1 - sf)
	fii_prime  <-  c(1 - sm, 1 - hm*sm, 1)	
	
	##SELECTION DIFFERENTIAL FEMALE-FUNCTION
	theta[7,]  <-  theta[7,]*fii	
	
	#SELECTION DIFFERENTIAL MALE-FUNCTION
	f_prime  <-  f_prime*fii_prime	
	

	####################################################
	# Population and female fertility perameters
	####################################################	
	sigmaS_J  <-  theta[1,]
	sigmaS_A  <-  theta[2,]
	gammaS    <-  theta[3,]
	sigmaX_J  <-  theta[4,]
	sigmaX_A  <-  theta[5,]
	gammaX    <-  theta[6,]
	f         <-  theta[7,]
	USi       <-  zeros(c(om,om,g))
	UXi       <-  zeros(c(om,om,g))
	FXi       <-  zeros(c(om,om,g))
	FSi       <-  zeros(c(om,om,g))
	FXi_pr    <-  zeros(c(om,om,g))
	lambda_i  <-  zeros(g)
	
	# create genotype-specific survival and fertility submatrices
	for (i in 1:3){
		USi[,,i]       <- rbind(c(sigmaS_J[i]*(1 - gammaS[i]), 0         ),
		                        c(sigmaS_J[i]*gammaS[i],       sigmaS_A[i]))
		UXi[,,i]       <- rbind(c(sigmaX_J[i]*(1 - gammaX[i]), 0         ),
		                        c(sigmaX_J[i]*gammaX[i],       sigmaX_A[i]))
		FSi[,,i]       <- rbind(c(0,C*(1 - delta)*f[i]),
						        c(0,0))
		FXi[,,i]       <- rbind(c(0,(1 - C)*f[i]),
							    c(0,0))
		FXi_pr[,,i]    <- rbind(c(0,f_prime[i]),
								c(0,0))	
		
		# Calculate gentoype-specific pop. growth rates
		lambda_i[i]  <-  max(Re(eigen(FSi[,,i] + USi[,,i] + FXi[,,i], symmetric=FALSE, only.values = TRUE)$values))
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
	    blkFS        <-  blkFS + kronecker(emat(g,g,i,i),FSi[,,i])
	    blkFX        <-  blkFX + kronecker(emat(g,g,i,i),FXi[,,i])
	    blkFX_prime  <-  blkFX_prime + kronecker(emat(g,g,i,i),FXi_pr[,,i])
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

	###############################################################
	# Simulating the Stage X genotype dynamics - generation loop
	###############################################################

	# Set initial state vectors for each boundary (fixed for AA & aa)
	n0AA  <-  c(      C*c(c(100-(om-1), rep(1,times=(om-1))), rep(0,times = 2*om)),
				(1 - C)*c(c(100-(om-1), rep(1,times=(om-1))), rep(0,times = 2*om)))
	n0aa  <-  c(      C*c(rep(0,times = 2*om), c(100-(om-1), rep(1,times=(om-1)))), 
				(1 - C)*c(rep(0,times = 2*om), c(100-(om-1), rep(1,times=(om-1)))))

	# Simulate to demographic equilibrium for each boundary
	AAEq  <-  fwdDyn2Eq(nzero=n0AA, om=om, g=g, W_prime=W_prime, blkFX_prime=blkFX_prime, blkUS=blkUS, blkUX=blkUX, W=W, Iom=Iom, Ig=Ig, K=K, Z=Z, blkFS=blkFS, blkFX=blkFX, tlimit=10^2, eqThreshold=eqThreshold)
	aaEq  <-  fwdDyn2Eq(nzero=n0aa, om=om, g=g, W_prime=W_prime, blkFX_prime=blkFX_prime, blkUS=blkUS, blkUX=blkUX, W=W, Iom=Iom, Ig=Ig, K=K, Z=Z, blkFS=blkFS, blkFX=blkFX, tlimit=10^2, eqThreshold=eqThreshold)


	# use intermediate frequency to find Eq.?
	if(intInit) {
		initEq  <-  c(      C*c(rep(c(32,rep(4/3, times=om-1)), times=g)),
					  (1 - C)*c(rep(c(32,rep(4/3, times=om-1)), times=g)))
		} else{
			# pick boundary for invasion based on major allele frequency		
			if(Ainvade) {
				initEq  <-  aaEq$p*100
			}
			if(!Ainvade) {
				initEq  <-  AAEq$p*100
			}
			# Introduce rare allele in heterozygotes
			initEq[(om+1)]       <-  0.1
			initEq[(om*g+om+1)]  <-  0.1
 		}

	# Simulate to demographic equilibrium with rare allele
	invadeEq  <-  fwdDyn2Eq(nzero=initEq, om=om, g=g, W_prime=W_prime, blkFX_prime=blkFX_prime, blkUS=blkUS, blkUX=blkUX, W=W, Iom=Iom, Ig=Ig, K=K, Z=Z, blkFS=blkFS, blkFX=blkFX, tlimit=tlimit, eqThreshold=eqThreshold)

	# Calculate coexistence conditions based on leading eigenvalue of the Jacobian
	# NOTE: we rearrange the order of the pHat_i values so they match the structure
	# of the jacobian, which is ordered by genotype then self/outcross
	nHat_AA  <-  AAEq$n[c(1:om, g*om+1:om, om+1:om, g*om+om+1:om, 2*om+1:om, g*om+2*om+1:om)]
	nHat_aa  <-  aaEq$n[c(1:om, g*om+1:om, om+1:om, g*om+om+1:om, 2*om+1:om, g*om+2*om+1:om)]
	pHat_AA  <-  nHat_AA/norm(as.matrix(nHat_AA), type="1")
	pHat_aa  <-  nHat_aa/norm(as.matrix(nHat_aa), type="1")
	zeta_i   <-  calcZeta(om=om, FSi=FSi, FXi=FXi, FXi_pr=FXi_pr, USi=USi, UXi=UXi,
						 pHat_AA=pHat_AA, pHat_aa=pHat_aa, C=C, delta=delta, 
						 lambda_i=lambda_i)

	##################
	# results
	res  <-  list(
					"nout"         =  invadeEq$nout,
					"p_genotypes"  =  invadeEq$p_genotypes,
					"pEq"          =  invadeEq$p_genotypes[,ncol(invadeEq$p_genotypes)],
					"eqReached"    =  invadeEq$eqReached,
					"lambda_i"     =  lambda_i,
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
					"lambda_sim"   =  invadeEq$lambda_sim,
					"extinct"      =  invadeEq$lambda_sim < 1,
					"polymorphism" =  !any(round(invadeEq$p_genotypes[,ncol(invadeEq$p_genotypes)], digits=3) == 1)
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
					  om = 2, g = 3, theta = c(0.6,0.6,0.05,6.1), theta_prime = 6, 
					  hf = 1/2, hm = 1/2, C = 0, delta = 0, 
					  delta_j = 0, delta_a = 0, delta_gamma = 0,
					  tlimit = 10^5, intInit = FALSE, eqThreshold = 1e-6, writeFile=TRUE, progBar = TRUE, ...) {
	
	sfs           <-  runif(min = 0, max=sMax, n=nSamples)
	sms           <-  runif(min = 0, max=sMax, n=nSamples)
	extinct       <-  rep(NA,times=nSamples)
	polymorphism  <-  rep(NA,times=nSamples)
	pEq           <-  matrix(0, ncol=3, nrow=nSamples)
	zeta_i        <-  matrix(0, ncol=2, nrow=nSamples)
	lambda_i      <-  matrix(0, ncol=3, nrow=nSamples)
	lambda_sim    <-  rep(0, nrow=nSamples)

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
										tlimit = 10^5, eqThreshold = eqThreshold, Ainvade = Ainvade, intInit = intInit)
			intInit          <-  FALSE
			extinct[i]       <-  results$extinct
			polymorphism[i]  <-  results$polymorphism
			pEq[i,]          <-  results$pEq
			zeta_i[i,]       <-  results$zeta_i
			lambda_i[i,]     <-  results$lambda_i
			lambda_sim[i,]   <-  results$lambda_sim

		if(progBar){
			cat('\r', paste(100*(i/nSamples),'% Complete'))
		}

	}

	# compile results as data frame
	hfVec       <-  rep(hf,    times=nSamples)
	hmVec       <-  rep(hm,    times=nSamples)
	CVec        <-  rep(C,     times=nSamples)
	deltaVec    <-  rep(delta, times=nSamples)
	results.df  <-  as.data.frame(cbind(sfs, 
										sms, 
										extinct, 
										polymorphism, 
										pEq, 
										zeta_i,
										lambda_i,
										lambda_sim,
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
								"lambda_sim",
								"hf",
								"hm",
								"C",
								"delta"
								)

	# export data as .csv to ./output/data
	if(writeFile) {
			filename <-  paste("./output/simData/demSimsSfxSmNew", "_sMax", sMax, "_nSamples", nSamples, "_hf", hf, "_hm", hm, 
							"_C", C, "_delta", delta, "_dj", delta_j, "_da", delta_a, "_dg", delta_gamma, "_f", theta[4], ".csv", sep="")
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
									 om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
									 hf = 1/2, hm = 1/2, C = 0, delta = 0, 
									 delta_j = 0, delta_a = 0, delta_gamma = 0,
									 tlimit = 10^5, eqThreshold = 1e-8) {
	
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
						 tlimit = tlimit, eqThreshold = eqThreshold, writeFile = FALSE, progBar = FALSE)

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
	filename <-  paste("./output/simData/simPolySpaceNew", "_sMax", sMax, "_nSamples", nSamples, "_hf", hf, "_hm", hm, 
					   "_delta", delta, "_dj", delta_j, "_da", delta_a, "_dg", delta_gamma, "_f", theta[4], ".csv", sep="")
	write.csv(results.df, file=filename, row.names = FALSE)

}






##############################
#' Quantify proportion of polymorphic AND/OR demographically viable
#' Parameter space for different Selfing rates & dominance
#'
deltaParamSpaceMakeData  <-  function(sMax = 0.15, nSamples=1e+2,
									om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = 6.5, 
									hf = 1/2, hm = 1/2, C = 1/2,
									tlimit = 10^5, eqThreshold = 1e-8) {
	
#	d    =  c(seq(0,0.3, by=0.05), seq(0.3, 0.9, by=0.1))
#	d_j  =  c(0, seq(0.1,0.2, by=0.02), seq(0.3, 0.9, by=0.1))
#	d_a  =  c(0, seq(0.1,0.2, by=0.02), seq(0.3, 0.9, by=0.1))
#	d_g  =  c(0,0.1,seq(0.2,0.3, by=0.02), seq(0.4, 0.9, by=0.1))
	d    =  seq(0, 0.9, by=0.1)
	d_j  =  seq(0, 0.9, by=0.1)
	d_a  =  seq(0, 0.9, by=0.1)
	d_g  =  seq(0, 0.9, by=0.1)

	equalLen  <-  length(d) == length(d_j) & length(d_j) == length(d_a) & length(d_a) == length(d_g)
	if(!equalLen) {
		stop('Need delta_i vectors of equal length')
	}

	d_extinct        <-  c()
	d_popGenPoly     <-  c()
	d_simPoly        <-  c()
	d_simPolyViable  <-  c()
	d_simAFix        <-  c()
	d_simAFixViable  <-  c()
	d_simaFix        <-  c()
	d_simaFixViable  <-  c()
	d_eigPoly        <-  c()
	d_eigPolyViable  <-  c()
	d_eigAFix        <-  c()
	d_eigAFixViable  <-  c()
	d_eigaFix        <-  c()
	d_eigaFixViable  <-  c()
	for(i in 1:length(d)) {
		d_res  <-  selLoop(sMax=sMax, nSamples=nSamples, 
							om = om, g = g, theta = theta, theta_prime = theta_prime, 
							hf = hf, hm = hm, C = C, delta = d[i],
							delta_j = 0, delta_a = 0, delta_gamma = 0,
							tlimit = tlimit, eqThreshold = eqThreshold, intInit = FALSE, writeFile = FALSE, progBar = FALSE)
	
		# Calculate & store results for plotting
		d_extinct[i]        <-  sum(d_res$extinct)/nSamples
		if(hf == hm & hm == 1/2) {
			d_popGenPoly[i]     <-  popGen_PolySpace_Delta_Add(C=C, delta=d[i], sMax=sMax)
		}
		if(hf == hm & hm != 1/2) {
			d_popGenPoly[i]     <-  popGen_PolySpace_Delta_DomRev(C=C, delta=d[i], sMax=sMax)
		}
		d_simPoly[i]        <-  sum(d_res$poly)/nSamples
		d_simPolyViable[i]  <-  d_simPoly[i] - sum(d_res$extinct == 1 & d_res$poly == 1)/nSamples
		d_simAFix[i]        <-  sum(d_res$poly == 0 & round(d_res$Eq_pAA) == 1)/nSamples
		d_simAFixViable[i]  <-  d_simAFix[i] - sum(d_res$extinct[d_res$poly == 0 & round(d_res$Eq_pAA) == 1])/nSamples
		d_simaFix[i]        <-  sum(d_res$poly == 0 & round(d_res$Eq_paa) == 1)/nSamples
		d_simaFixViable[i]  <-  d_simaFix[i] - sum(d_res$extinct[d_res$poly == 0 & round(d_res$Eq_paa) == 1])/nSamples
		d_eigPoly [i]       <-  sum(d_res$zeta_AA > 1 & d_res$zeta_aa > 1)/nSamples
		d_eigPolyViable[i]  <-  d_eigPoly[i] - sum(d_res$extinct == 1 & d_res$zeta_AA > 1 & d_res$zeta_aa > 1)/nSamples
		d_eigAFix[i]        <-  sum(d_res$zeta_AA > 1 & d_res$zeta_aa < 1)/nSamples
		d_eigAFixViable[i]  <-  d_eigAFix[i] - sum(d_res$extinct[d_res$poly == 0 & d_res$zeta_AA > 1 & d_res$zeta_aa < 1])/nSamples
		d_eigaFix[i]        <-  sum(d_res$zeta_AA < 1 & d_res$zeta_aa > 1)/nSamples
		d_eigaFixViable[i]  <-  d_eigaFix[i] - sum(d_res$extinct[d_res$poly == 0 & d_res$zeta_AA < 1 & d_res$zeta_aa > 1])/nSamples
	
		cat('\r', paste("d gradient ", round(100*(i/length(d))),'% Complete'))
		flush.console()
	}

	d_j_extinct        <-  c()
	d_j_popGenPoly     <-  c()
	d_j_simPoly        <-  c()
	d_j_simPolyViable  <-  c()
	d_j_simAFix        <-  c()
	d_j_simAFixViable  <-  c()
	d_j_simaFix        <-  c()
	d_j_simaFixViable  <-  c()
	d_j_eigPoly        <-  c()
	d_j_eigPolyViable  <-  c()
	d_j_eigAFix        <-  c()
	d_j_eigAFixViable  <-  c()
	d_j_eigaFix        <-  c()
	d_j_eigaFixViable  <-  c()
	for(i in 1:length(d_j)) {
		d_j_res  <-  selLoop(sMax=sMax, nSamples=nSamples, 
							om = om, g = g, theta = theta, theta_prime = theta_prime, 
							hf = hf, hm = hm, C = C, delta = 0,
							delta_j = d_j[i], delta_a = 0, delta_gamma = 0,
							tlimit = tlimit, eqThreshold = eqThreshold, intInit = FALSE, writeFile = FALSE, progBar = FALSE)
	
		# Calculate & store results for plotting
		d_j_extinct[i]        <-  sum(d_j_res$extinct)/nSamples
		if(hf == hm & hm == 1/2) {
			d_j_popGenPoly[i]     <-  popGen_PolySpace_Delta_Add(C=C, delta=d[i], sMax=sMax)
		}
		if(hf == hm & hm != 1/2) {
			d_j_popGenPoly[i]     <-  popGen_PolySpace_Delta_DomRev(C=C, delta=d[i], sMax=sMax)
		}
		d_j_simPoly[i]        <-  sum(d_j_res$poly)/nSamples
		d_j_simPolyViable[i]  <-  d_j_simPoly[i] - sum(d_j_res$extinct == 1 & d_j_res$poly == 1)/nSamples
		d_j_simAFix[i]        <-  sum(d_j_res$poly == 0 & round(d_j_res$Eq_pAA) == 1)/nSamples
		d_j_simAFixViable[i]  <-  d_j_simAFix[i] - sum(d_j_res$extinct[d_j_res$poly == 0 & round(d_j_res$Eq_pAA) == 1])/nSamples
		d_j_simaFix[i]        <-  sum(d_j_res$poly == 0 & round(d_j_res$Eq_paa) == 1)/nSamples
		d_j_simaFixViable[i]  <-  d_j_simaFix[i] - sum(d_j_res$extinct[d_j_res$poly == 0 & round(d_j_res$Eq_paa) == 1])/nSamples
		d_j_eigPoly [i]       <-  sum(d_j_res$zeta_AA > 1 & d_j_res$zeta_aa > 1)/nSamples
		d_j_eigPolyViable[i]  <-  d_j_eigPoly[i] - sum(d_j_res$extinct == 1 & d_j_res$zeta_AA > 1 & d_j_res$zeta_aa > 1)/nSamples
		d_j_eigAFix[i]        <-  sum(d_j_res$zeta_AA > 1 & d_j_res$zeta_aa < 1)/nSamples
		d_j_eigAFixViable[i]  <-  d_j_eigAFix[i] - sum(d_j_res$extinct[d_j_res$poly == 0 & d_j_res$zeta_AA > 1 & d_j_res$zeta_aa < 1])/nSamples
		d_j_eigaFix[i]        <-  sum(d_j_res$zeta_AA < 1 & d_j_res$zeta_aa > 1)/nSamples
		d_j_eigaFixViable[i]  <-  d_j_eigaFix[i] - sum(d_j_res$extinct[d_j_res$poly == 0 & d_j_res$zeta_AA < 1 & d_j_res$zeta_aa > 1])/nSamples
	
		cat('\r', paste("d_j gradient ", round(100*(i/length(d))),'% Complete'))
		flush.console()
	}

	d_a_extinct        <-  c()
	d_a_popGenPoly     <-  c()
	d_a_simPoly        <-  c()
	d_a_simPolyViable  <-  c()
	d_a_simAFix        <-  c()
	d_a_simAFixViable  <-  c()
	d_a_simaFix        <-  c()
	d_a_simaFixViable  <-  c()
	d_a_eigPoly        <-  c()
	d_a_eigPolyViable  <-  c()
	d_a_eigAFix        <-  c()
	d_a_eigAFixViable  <-  c()
	d_a_eigaFix        <-  c()
	d_a_eigaFixViable  <-  c()
	for(i in 1:length(d_a)) {

		d_a_res  <-  selLoop(sMax=sMax, nSamples=nSamples, 
							om = om, g = g, theta = theta, theta_prime = theta_prime, 
							hf = hf, hm = hm, C = C, delta = 0,
							delta_j = 0, delta_a = d_a[i], delta_gamma = 0,
							tlimit = tlimit, eqThreshold = eqThreshold, intInit = FALSE, writeFile = FALSE, progBar = FALSE)
	
		# Calculate & store results for plotting
		d_a_extinct[i]        <-  sum(d_a_res$extinct)/nSamples
		if(hf == hm & hm == 1/2) {
			d_a_popGenPoly[i]     <-  popGen_PolySpace_Delta_Add(C=C, delta=d[i], sMax=sMax)
		}
		if(hf == hm & hm != 1/2) {
			d_a_popGenPoly[i]     <-  popGen_PolySpace_Delta_DomRev(C=C, delta=d[i], sMax=sMax)
		}
		d_a_simPoly[i]        <-  sum(d_a_res$poly)/nSamples
		d_a_simPolyViable[i]  <-  d_a_simPoly[i] - sum(d_a_res$extinct == 1 & d_a_res$poly == 1)/nSamples
		d_a_simAFix[i]        <-  sum(d_a_res$poly == 0 & round(d_a_res$Eq_pAA) == 1)/nSamples
		d_a_simAFixViable[i]  <-  d_a_simAFix[i] - sum(d_a_res$extinct[d_a_res$poly == 0 & round(d_a_res$Eq_pAA) == 1])/nSamples
		d_a_simaFix[i]        <-  sum(d_a_res$poly == 0 & round(d_a_res$Eq_paa) == 1)/nSamples
		d_a_simaFixViable[i]  <-  d_a_simaFix[i] - sum(d_a_res$extinct[d_a_res$poly == 0 & round(d_a_res$Eq_paa) == 1])/nSamples
		d_a_eigPoly [i]       <-  sum(d_a_res$zeta_AA > 1 & d_a_res$zeta_aa > 1)/nSamples
		d_a_eigPolyViable[i]  <-  d_a_eigPoly[i] - sum(d_a_res$extinct == 1 & d_a_res$zeta_AA > 1 & d_a_res$zeta_aa > 1)/nSamples
		d_a_eigAFix[i]        <-  sum(d_a_res$zeta_AA > 1 & d_a_res$zeta_aa < 1)/nSamples
		d_a_eigAFixViable[i]  <-  d_a_eigAFix[i] - sum(d_a_res$extinct[d_a_res$poly == 0 & d_a_res$zeta_AA > 1 & d_a_res$zeta_aa < 1])/nSamples
		d_a_eigaFix[i]        <-  sum(d_a_res$zeta_AA < 1 & d_a_res$zeta_aa > 1)/nSamples
		d_a_eigaFixViable[i]  <-  d_a_eigaFix[i] - sum(d_a_res$extinct[d_a_res$poly == 0 & d_a_res$zeta_AA < 1 & d_a_res$zeta_aa > 1])/nSamples
	
		cat('\r', paste("d_a gradient ", round(100*(i/length(d_a))),'% Complete'))
		flush.console()
	}

	d_g_extinct        <-  c()
	d_g_popGenPoly     <-  c()
	d_g_simPoly        <-  c()
	d_g_simPolyViable  <-  c()
	d_g_simAFix        <-  c()
	d_g_simAFixViable  <-  c()
	d_g_simaFix        <-  c()
	d_g_simaFixViable  <-  c()
	d_g_eigPoly        <-  c()
	d_g_eigPolyViable  <-  c()
	d_g_eigAFix        <-  c()
	d_g_eigAFixViable  <-  c()
	d_g_eigaFix        <-  c()
	d_g_eigaFixViable  <-  c()
	for(i in 1:length(d_g)) {
		d_g_res  <-  selLoop(sMax=sMax, nSamples=nSamples, 
							om = om, g = g, theta = theta, theta_prime = theta_prime, 
							hf = hf, hm = hm, C = C, delta = 0,
							delta_j = 0, delta_a = 0, delta_gamma = d_g[i],
							tlimit = tlimit, eqThreshold = eqThreshold, intInit = FALSE, writeFile = FALSE, progBar = FALSE)
	
		# Calculate & store results for plotting
		d_g_extinct[i]        <-  sum(d_g_res$extinct)/nSamples
		if(hf == hm & hm == 1/2) {
			d_g_popGenPoly[i]     <-  popGen_PolySpace_Delta_Add(C=C, delta=d[i], sMax=sMax)
		}
		if(hf == hm & hm != 1/2) {
			d_g_popGenPoly[i]     <-  popGen_PolySpace_Delta_DomRev(C=C, delta=d[i], sMax=sMax)
		}
		d_g_simPoly[i]        <-  sum(d_g_res$poly)/nSamples
		d_g_simPolyViable[i]  <-  d_g_simPoly[i] - sum(d_g_res$extinct == 1 & d_g_res$poly == 1)/nSamples
		d_g_simAFix[i]        <-  sum(d_g_res$poly == 0 & round(d_g_res$Eq_pAA) == 1)/nSamples
		d_g_simAFixViable[i]  <-  d_g_simAFix[i] - sum(d_g_res$extinct[d_g_res$poly == 0 & round(d_g_res$Eq_pAA) == 1])/nSamples
		d_g_simaFix[i]        <-  sum(d_g_res$poly == 0 & round(d_g_res$Eq_paa) == 1)/nSamples
		d_g_simaFixViable[i]  <-  d_g_simaFix[i] - sum(d_g_res$extinct[d_g_res$poly == 0 & round(d_g_res$Eq_paa) == 1])/nSamples
		d_g_eigPoly [i]       <-  sum(d_g_res$zeta_AA > 1 & d_g_res$zeta_aa > 1)/nSamples
		d_g_eigPolyViable[i]  <-  d_g_eigPoly[i] - sum(d_g_res$extinct == 1 & d_g_res$zeta_AA > 1 & d_g_res$zeta_aa > 1)/nSamples
		d_g_eigAFix[i]        <-  sum(d_g_res$zeta_AA > 1 & d_g_res$zeta_aa < 1)/nSamples
		d_g_eigAFixViable[i]  <-  d_g_eigAFix[i] - sum(d_g_res$extinct[d_g_res$poly == 0 & d_g_res$zeta_AA > 1 & d_g_res$zeta_aa < 1])/nSamples
		d_g_eigaFix[i]        <-  sum(d_g_res$zeta_AA < 1 & d_g_res$zeta_aa > 1)/nSamples
		d_g_eigaFixViable[i]  <-  d_g_eigaFix[i] - sum(d_g_res$extinct[d_g_res$poly == 0 & d_g_res$zeta_AA < 1 & d_g_res$zeta_aa > 1])/nSamples
	
		cat('\r', paste("d_gamma gradient ", round(100*(i/length(d))),'% Complete'))
		flush.console()
	}


	# compile results as data frame
	results.df  <-  as.data.frame(cbind(d,
										d_j,
										d_a,
										d_g,
										d_extinct,
										d_popGenPoly,
										d_simPoly,
										d_simPolyViable,
										d_simAFix,
										d_simAFixViable,
										d_simaFix,
										d_simaFixViable,
										d_eigPoly,
										d_eigPolyViable,
										d_eigAFix,
										d_eigAFixViable,
										d_eigaFix,
										d_eigaFixViable,
										d_j_extinct,
										d_j_popGenPoly,
										d_j_simPoly,
										d_j_simPolyViable,
										d_j_simAFix,
										d_j_simAFixViable,
										d_j_simaFix,
										d_j_simaFixViable,
										d_j_eigPoly,
										d_j_eigPolyViable,
										d_j_eigAFix,
										d_j_eigAFixViable,
										d_j_eigaFix,
										d_j_eigaFixViable,
										d_a_extinct,
										d_a_popGenPoly,
										d_a_simPoly,
										d_a_simPolyViable,
										d_a_simAFix,
										d_a_simAFixViable,
										d_a_simaFix,
										d_a_simaFixViable,
										d_a_eigPoly,
										d_a_eigPolyViable,
										d_a_eigAFix,
										d_a_eigAFixViable,
										d_a_eigaFix,
										d_a_eigaFixViable,
										d_g_extinct,
										d_g_popGenPoly,
										d_g_simPoly,
										d_g_simPolyViable,
										d_g_simAFix,
										d_g_simAFixViable,
										d_g_simaFix,
										d_g_simaFixViable,
										d_g_eigPoly,
										d_g_eigPolyViable,
										d_g_eigAFix,
										d_g_eigAFixViable,
										d_g_eigaFix,
										d_g_eigaFixViable
										)
								 )

	# Export data as .csv to ./output/data
	filename <-  paste("./output/simData/deltaSimPolySpace", "_sMax", sMax, "_nSamples", nSamples, "_C", C, "_hf", hf, "_hm", hm,
					   "_f", theta[4], ".csv", sep="")
	write.csv(results.df, file=filename, row.names = FALSE)

}






##############################
#' Quantify proportion of polymorphic AND/OR demographically viable
#' Parameter space for different Selfing rates & dominance
#'
deltaSelfingPolySpaceMakeData  <-  function(sMax = 0.15, nSamples=1e+2,
											om = 2, g = 3, theta = c(0.6,0.6,0.05,7), theta_prime = 7, 
											hf = 1/2, hm = 1/2, dStar = 0.8, 
											tlimit = 10^5, eqThreshold = 1e-8) {
	
	if(hf == 1/2 & hm == 1/2) {
		Cs        <-  seq(0, 0.9, by=0.05)
	}
	if(hf ==hm & hm < 1/2) {
		Cs        <-  c(0.01, seq(0.05, 0.9, by=0.05))
	}	
	deltaSeq  <-  predDelta(dStar=dStar, b=1/2, a=0.2, C=Cs) 

	d_extinct        <-  c()
	d_popGenPoly     <-  c()
	d_simPoly        <-  c()
	d_simPolyViable  <-  c()
	d_simAFix        <-  c()
	d_simAFixViable  <-  c()
	d_simaFix        <-  c()
	d_simaFixViable  <-  c()
	d_eigPoly        <-  c()
	d_eigPolyViable  <-  c()
	d_eigAFix        <-  c()
	d_eigAFixViable  <-  c()
	d_eigaFix        <-  c()
	d_eigaFixViable  <-  c()
	for(i in 1:length(deltaSeq)) {
		d_res  <-  selLoop(sMax=sMax, nSamples=nSamples, 
							om = om, g = g, theta = theta, theta_prime = theta_prime, 
							hf = hf, hm = hm, C = Cs[i], delta = deltaSeq[i],
							delta_j = 0, delta_a = 0, delta_gamma = 0,
							tlimit = tlimit, eqThreshold = eqThreshold, 
							intInit = FALSE, writeFile = FALSE, progBar = FALSE)
	
		# Calculate & store results for plotting
		d_extinct[i]        <-  sum(d_res$extinct)/nSamples
		if(hf == hm & hm == 1/2) {
			d_popGenPoly[i]     <-  popGen_PolySpace_Delta_Add(C=Cs[i], delta=deltaSeq[i], sMax=sMax)
		}
		if(hf == hm & hm != 1/2) {
			d_popGenPoly[i]     <-  popGen_PolySpace_Delta_DomRev(C=Cs[i], delta=deltaSeq[i], sMax=sMax)
		}
		d_simPoly[i]        <-  sum(d_res$poly)/nSamples
		d_simPolyViable[i]  <-  d_simPoly[i] - sum(d_res$extinct == 1 & d_res$poly == 1)/nSamples
		d_simAFix[i]        <-  sum(d_res$poly == 0 & round(d_res$Eq_pAA) == 1)/nSamples
		d_simAFixViable[i]  <-  d_simAFix[i] - sum(d_res$extinct[d_res$poly == 0 & round(d_res$Eq_pAA) == 1])/nSamples
		d_simaFix[i]        <-  sum(d_res$poly == 0 & round(d_res$Eq_paa) == 1)/nSamples
		d_simaFixViable[i]  <-  d_simaFix[i] - sum(d_res$extinct[d_res$poly == 0 & round(d_res$Eq_paa) == 1])/nSamples
		d_eigPoly [i]       <-  sum(d_res$zeta_AA > 1 & d_res$zeta_aa > 1)/nSamples
		d_eigPolyViable[i]  <-  d_eigPoly[i] - sum(d_res$extinct == 1 & d_res$zeta_AA > 1 & d_res$zeta_aa > 1)/nSamples
		d_eigAFix[i]        <-  sum(d_res$zeta_AA > 1 & d_res$zeta_aa < 1)/nSamples
		d_eigAFixViable[i]  <-  d_eigAFix[i] - sum(d_res$extinct[d_res$poly == 0 & d_res$zeta_AA > 1 & d_res$zeta_aa < 1])/nSamples
		d_eigaFix[i]        <-  sum(d_res$zeta_AA < 1 & d_res$zeta_aa > 1)/nSamples
		d_eigaFixViable[i]  <-  d_eigaFix[i] - sum(d_res$extinct[d_res$poly == 0 & d_res$zeta_AA < 1 & d_res$zeta_aa > 1])/nSamples
	
		cat('\r', paste("d gradient ", round(100*(i/length(deltaSeq))),'% Complete'))
		flush.console()
	}

	d_j_extinct        <-  c()
	d_j_popGenPoly     <-  c()
	d_j_simPoly        <-  c()
	d_j_simPolyViable  <-  c()
	d_j_simAFix        <-  c()
	d_j_simAFixViable  <-  c()
	d_j_simaFix        <-  c()
	d_j_simaFixViable  <-  c()
	d_j_eigPoly        <-  c()
	d_j_eigPolyViable  <-  c()
	d_j_eigAFix        <-  c()
	d_j_eigAFixViable  <-  c()
	d_j_eigaFix        <-  c()
	d_j_eigaFixViable  <-  c()
	for(i in 1:length(deltaSeq)) {
		d_j_res  <-  selLoop(sMax=sMax, nSamples=nSamples, 
							om = om, g = g, theta = theta, theta_prime = theta_prime, 
							hf = hf, hm = hm, C = Cs[i], delta = 0,
							delta_j = deltaSeq[i], delta_a = 0, delta_gamma = 0,
							tlimit = tlimit, eqThreshold = eqThreshold, 
							intInit = FALSE, writeFile = FALSE, progBar = FALSE)
	
		# Calculate & store results for plotting
		d_j_extinct[i]        <-  sum(d_j_res$extinct)/nSamples
		if(hf == hm & hm == 1/2) {
			d_j_popGenPoly[i]     <-  popGen_PolySpace_Delta_Add(C=Cs[i], delta=deltaSeq[i], sMax=sMax)
		}
		if(hf == hm & hm != 1/2) {
			d_j_popGenPoly[i]     <-  popGen_PolySpace_Delta_DomRev(C=Cs[i], delta=deltaSeq[i], sMax=sMax)
		}
		d_j_simPoly[i]        <-  sum(d_j_res$poly)/nSamples
		d_j_simPolyViable[i]  <-  d_j_simPoly[i] - sum(d_j_res$extinct == 1 & d_j_res$poly == 1)/nSamples
		d_j_simAFix[i]        <-  sum(d_j_res$poly == 0 & round(d_j_res$Eq_pAA) == 1)/nSamples
		d_j_simAFixViable[i]  <-  d_j_simAFix[i] - sum(d_j_res$extinct[d_j_res$poly == 0 & round(d_j_res$Eq_pAA) == 1])/nSamples
		d_j_simaFix[i]        <-  sum(d_j_res$poly == 0 & round(d_j_res$Eq_paa) == 1)/nSamples
		d_j_simaFixViable[i]  <-  d_j_simaFix[i] - sum(d_j_res$extinct[d_j_res$poly == 0 & round(d_j_res$Eq_paa) == 1])/nSamples
		d_j_eigPoly [i]       <-  sum(d_j_res$zeta_AA > 1 & d_j_res$zeta_aa > 1)/nSamples
		d_j_eigPolyViable[i]  <-  d_j_eigPoly[i] - sum(d_j_res$extinct == 1 & d_j_res$zeta_AA > 1 & d_j_res$zeta_aa > 1)/nSamples
		d_j_eigAFix[i]        <-  sum(d_j_res$zeta_AA > 1 & d_j_res$zeta_aa < 1)/nSamples
		d_j_eigAFixViable[i]  <-  d_j_eigAFix[i] - sum(d_j_res$extinct[d_j_res$poly == 0 & d_j_res$zeta_AA > 1 & d_j_res$zeta_aa < 1])/nSamples
		d_j_eigaFix[i]        <-  sum(d_j_res$zeta_AA < 1 & d_j_res$zeta_aa > 1)/nSamples
		d_j_eigaFixViable[i]  <-  d_j_eigaFix[i] - sum(d_j_res$extinct[d_j_res$poly == 0 & d_j_res$zeta_AA < 1 & d_j_res$zeta_aa > 1])/nSamples
	
		cat('\r', paste("d_j gradient ", round(100*(i/length(deltaSeq))),'% Complete'))
		flush.console()
	}

	d_a_extinct        <-  c()
	d_a_popGenPoly     <-  c()
	d_a_simPoly        <-  c()
	d_a_simPolyViable  <-  c()
	d_a_simAFix        <-  c()
	d_a_simAFixViable  <-  c()
	d_a_simaFix        <-  c()
	d_a_simaFixViable  <-  c()
	d_a_eigPoly        <-  c()
	d_a_eigPolyViable  <-  c()
	d_a_eigAFix        <-  c()
	d_a_eigAFixViable  <-  c()
	d_a_eigaFix        <-  c()
	d_a_eigaFixViable  <-  c()
	for(i in 1:length(deltaSeq)) {

		d_a_res  <-  selLoop(sMax=sMax, nSamples=nSamples, 
							om = om, g = g, theta = theta, theta_prime = theta_prime, 
							hf = hf, hm = hm, C = Cs[i], delta = 0,
							delta_j = 0, delta_a = deltaSeq[i], delta_gamma = 0,
							tlimit = tlimit, eqThreshold = eqThreshold, 
							intInit = FALSE, writeFile = FALSE, progBar = FALSE)
	
		# Calculate & store results for plotting
		d_a_extinct[i]        <-  sum(d_a_res$extinct)/nSamples
		if(hf == hm & hm == 1/2) {
			d_a_popGenPoly[i]     <-  popGen_PolySpace_Delta_Add(C=Cs[i], delta=deltaSeq[i], sMax=sMax)
		}
		if(hf == hm & hm != 1/2) {
			d_a_popGenPoly[i]     <-  popGen_PolySpace_Delta_DomRev(C=Cs[i], delta=deltaSeq[i], sMax=sMax)
		}
		d_a_simPoly[i]        <-  sum(d_a_res$poly)/nSamples
		d_a_simPolyViable[i]  <-  d_a_simPoly[i] - sum(d_a_res$extinct == 1 & d_a_res$poly == 1)/nSamples
		d_a_simAFix[i]        <-  sum(d_a_res$poly == 0 & round(d_a_res$Eq_pAA) == 1)/nSamples
		d_a_simAFixViable[i]  <-  d_a_simAFix[i] - sum(d_a_res$extinct[d_a_res$poly == 0 & round(d_a_res$Eq_pAA) == 1])/nSamples
		d_a_simaFix[i]        <-  sum(d_a_res$poly == 0 & round(d_a_res$Eq_paa) == 1)/nSamples
		d_a_simaFixViable[i]  <-  d_a_simaFix[i] - sum(d_a_res$extinct[d_a_res$poly == 0 & round(d_a_res$Eq_paa) == 1])/nSamples
		d_a_eigPoly [i]       <-  sum(d_a_res$zeta_AA > 1 & d_a_res$zeta_aa > 1)/nSamples
		d_a_eigPolyViable[i]  <-  d_a_eigPoly[i] - sum(d_a_res$extinct == 1 & d_a_res$zeta_AA > 1 & d_a_res$zeta_aa > 1)/nSamples
		d_a_eigAFix[i]        <-  sum(d_a_res$zeta_AA > 1 & d_a_res$zeta_aa < 1)/nSamples
		d_a_eigAFixViable[i]  <-  d_a_eigAFix[i] - sum(d_a_res$extinct[d_a_res$poly == 0 & d_a_res$zeta_AA > 1 & d_a_res$zeta_aa < 1])/nSamples
		d_a_eigaFix[i]        <-  sum(d_a_res$zeta_AA < 1 & d_a_res$zeta_aa > 1)/nSamples
		d_a_eigaFixViable[i]  <-  d_a_eigaFix[i] - sum(d_a_res$extinct[d_a_res$poly == 0 & d_a_res$zeta_AA < 1 & d_a_res$zeta_aa > 1])/nSamples
	
		cat('\r', paste("d_a gradient ", round(100*(i/length(deltaSeq))),'% Complete'))
		flush.console()
	}

	d_g_extinct        <-  c()
	d_g_popGenPoly     <-  c()
	d_g_simPoly        <-  c()
	d_g_simPolyViable  <-  c()
	d_g_simAFix        <-  c()
	d_g_simAFixViable  <-  c()
	d_g_simaFix        <-  c()
	d_g_simaFixViable  <-  c()
	d_g_eigPoly        <-  c()
	d_g_eigPolyViable  <-  c()
	d_g_eigAFix        <-  c()
	d_g_eigAFixViable  <-  c()
	d_g_eigaFix        <-  c()
	d_g_eigaFixViable  <-  c()
	for(i in 1:length(deltaSeq)) {
		d_g_res  <-  selLoop(sMax=sMax, nSamples=nSamples, 
							om = om, g = g, theta = theta, theta_prime = theta_prime, 
							hf = hf, hm = hm, C = Cs[i], delta = 0,
							delta_j = 0, delta_a = 0, delta_gamma = deltaSeq[i],
							tlimit = tlimit, eqThreshold = eqThreshold, 
							intInit = FALSE, writeFile = FALSE, progBar = FALSE)
	
		# Calculate & store results for plotting
		d_g_extinct[i]        <-  sum(d_g_res$extinct)/nSamples
		if(hf == hm & hm == 1/2) {
			d_g_popGenPoly[i]     <-  popGen_PolySpace_Delta_Add(C=Cs[i], delta=deltaSeq[i], sMax=sMax)
		}
		if(hf == hm & hm != 1/2) {
			d_g_popGenPoly[i]     <-  popGen_PolySpace_Delta_DomRev(C=Cs[i], delta=deltaSeq[i], sMax=sMax)
		}
		d_g_simPoly[i]        <-  sum(d_g_res$poly)/nSamples
		d_g_simPolyViable[i]  <-  d_g_simPoly[i] - sum(d_g_res$extinct == 1 & d_g_res$poly == 1)/nSamples
		d_g_simAFix[i]        <-  sum(d_g_res$poly == 0 & round(d_g_res$Eq_pAA) == 1)/nSamples
		d_g_simAFixViable[i]  <-  d_g_simAFix[i] - sum(d_g_res$extinct[d_g_res$poly == 0 & round(d_g_res$Eq_pAA) == 1])/nSamples
		d_g_simaFix[i]        <-  sum(d_g_res$poly == 0 & round(d_g_res$Eq_paa) == 1)/nSamples
		d_g_simaFixViable[i]  <-  d_g_simaFix[i] - sum(d_g_res$extinct[d_g_res$poly == 0 & round(d_g_res$Eq_paa) == 1])/nSamples
		d_g_eigPoly [i]       <-  sum(d_g_res$zeta_AA > 1 & d_g_res$zeta_aa > 1)/nSamples
		d_g_eigPolyViable[i]  <-  d_g_eigPoly[i] - sum(d_g_res$extinct == 1 & d_g_res$zeta_AA > 1 & d_g_res$zeta_aa > 1)/nSamples
		d_g_eigAFix[i]        <-  sum(d_g_res$zeta_AA > 1 & d_g_res$zeta_aa < 1)/nSamples
		d_g_eigAFixViable[i]  <-  d_g_eigAFix[i] - sum(d_g_res$extinct[d_g_res$poly == 0 & d_g_res$zeta_AA > 1 & d_g_res$zeta_aa < 1])/nSamples
		d_g_eigaFix[i]        <-  sum(d_g_res$zeta_AA < 1 & d_g_res$zeta_aa > 1)/nSamples
		d_g_eigaFixViable[i]  <-  d_g_eigaFix[i] - sum(d_g_res$extinct[d_g_res$poly == 0 & d_g_res$zeta_AA < 1 & d_g_res$zeta_aa > 1])/nSamples
	
		cat('\r', paste("d_gamma gradient ", round(100*(i/length(deltaSeq))),'% Complete'))
		flush.console()
	}


	# compile results as data frame
	results.df  <-  as.data.frame(cbind(Cs,
										deltaSeq,
										d_extinct,
										d_popGenPoly,
										d_simPoly,
										d_simPolyViable,
										d_simAFix,
										d_simAFixViable,
										d_simaFix,
										d_simaFixViable,
										d_eigPoly,
										d_eigPolyViable,
										d_eigAFix,
										d_eigAFixViable,
										d_eigaFix,
										d_eigaFixViable,
										d_j_extinct,
										d_j_popGenPoly,
										d_j_simPoly,
										d_j_simPolyViable,
										d_j_simAFix,
										d_j_simAFixViable,
										d_j_simaFix,
										d_j_simaFixViable,
										d_j_eigPoly,
										d_j_eigPolyViable,
										d_j_eigAFix,
										d_j_eigAFixViable,
										d_j_eigaFix,
										d_j_eigaFixViable,
										d_a_extinct,
										d_a_popGenPoly,
										d_a_simPoly,
										d_a_simPolyViable,
										d_a_simAFix,
										d_a_simAFixViable,
										d_a_simaFix,
										d_a_simaFixViable,
										d_a_eigPoly,
										d_a_eigPolyViable,
										d_a_eigAFix,
										d_a_eigAFixViable,
										d_a_eigaFix,
										d_a_eigaFixViable,
										d_g_extinct,
										d_g_popGenPoly,
										d_g_simPoly,
										d_g_simPolyViable,
										d_g_simAFix,
										d_g_simAFixViable,
										d_g_simaFix,
										d_g_simaFixViable,
										d_g_eigPoly,
										d_g_eigPolyViable,
										d_g_eigAFix,
										d_g_eigAFixViable,
										d_g_eigaFix,
										d_g_eigaFixViable
										)
								 )

	# Export data as .csv to ./output/data
	filename <-  paste("./output/simData/deltaSelfingSimPolySpaceNew", "_sMax", sMax, "_nSamples", nSamples, "_dStar", dStar, "_hf", hf, "_hm", hm,
					   "_f", theta[4], ".csv", sep="")
	write.csv(results.df, file=filename, row.names = FALSE)

}



midPoint  <-  function(x1,x2) {
	x1 + (x2-x1)/2
}


##############################
#' Identify extinction threshold across sf x sm parameter space
#' faster method for calculating proportions of sf x sm where
#' different dynamical outcomes happen(?)
#'
extinctThreshTitrate  <-  function(sMax=0.15, res=0.0015, 
								   om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
								   hf = 1/2, hm = 1/2, C = 0, delta = 0, 
								   delta_j = 0, delta_a = 0, delta_gamma = 0,
								   tlimit = 10^5, Ainvade=FALSE, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE) {
	
	# vector of sf values to titrate along
	sfs          <-  seq(res,(sMax-res),by=res)
	sms          <-  sfs
	smThreshold  <-  rep(NA, times=length(sfs))
	sfThreshold  <-  smThreshold

	# Three values we need to keep track of
	titrateStartVals  <-  matrix(rbind(c(sfs[1], sfs[length(sfs)]),
								  c(NA,NA)),nrow=2,ncol=2)
	rownames(titrateStartVals)  <-  c("sVals", "extinct")
	colnames(titrateStartVals)  <-  c("left", "right")

	# loop over sf values
	for(i in 1:length(sfs)) {

		titrateVals  <-  titrateStartVals

		# start at boundaries of [0, sMax]
		L  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
								 hf = hf, hm = hm, sf = sfs[i], sm = titrateVals[1,1], C = C, delta = delta, 
								 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
								 tlimit = 10^5, eqThreshold = eqThreshold, Ainvade = FALSE, intInit = TRUE)
		R  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
								 hf = hf, hm = hm, sf = sfs[i], sm = titrateVals[1,2], C = C, delta = delta, 
								 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
								 tlimit = 10^5, eqThreshold = eqThreshold, Ainvade = FALSE, intInit = TRUE)
		
		# check to see if extinction outcomes differ. If not, go to next sf value
		if(all(L$lambda_sim > 1 && R$lambda_sim > 1)) {
			if(verbose) {
					cat('\r', paste("sf grad.:", round(100*(i/length(sfs))), "% complete"))
					flush.console()
				}
			next
		} else {
			# if they do, titrate to estimate threshold
			titrateVals[2,]  <-  c(L$lambda_sim, R$lambda_sim)
			titrateDelta     <-  1

			while(!all(titrateVals[2,1] > 1 && titrateVals[2,2] > 1) & titrateDelta > res) {
				sm_prop       <-  midPoint(titrateVals[1,1],titrateVals[1,2])
				titrateDelta  <-  sm_prop - min(titrateVals[1,1], titrateVals[1,2])
			
				proposal  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
											 hf = hf, hm = hm, sf = sfs[i], sm = sm_prop, C = C, delta = delta, 
											 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
											 tlimit = 10^5, eqThreshold = eqThreshold, Ainvade = Ainvade, intInit = TRUE)
				# replace appropraite boundary value
				if(proposal$lambda_sim  < 1) {
					titrateVals[1,][titrateVals[2,] < 1]  <-  sm_prop
					titrateVals[2,][titrateVals[2,] < 1]  <-  proposal$lambda_sim
				}
				if(proposal$lambda_sim  > 1) {
					titrateVals[1,][titrateVals[2,] > 1]  <-  sm_prop
					titrateVals[2,][titrateVals[2,] > 1]  <-  proposal$lambda_sim
				}
				if(verbose) {
					cat('\r', paste("sf grad.:", round(100*(i/length(sfs))), "% complete; Titration ratio = ", round(titrateDelta/res, digits=2), "(done when < 1)"))
					flush.console()
				}

			}

		}
		# Record sm threshold value
		smThreshold[i]  <-  sm_prop
	}

	# loop over sm values
	for(i in 1:length(sms)) {

		titrateVals  <-  titrateStartVals

		# start at boundaries of [0, sMax]
		L  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
								 hf = hf, hm = hm, sf = titrateVals[1,1], sm = sms[i], C = C, delta = delta, 
								 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
								 tlimit = tlimit, eqThreshold = eqThreshold, Ainvade = Ainvade, intInit = TRUE)
		R  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
								 hf = hf, hm = hm, sf = titrateVals[1,2], sm = sms[i], C = C, delta = delta, 
								 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
								 tlimit = tlimit, eqThreshold = eqThreshold, Ainvade = Ainvade, intInit = TRUE)
		
		# check to see if extinction outcomes differ. If not, go to next sf value
		if(all(L$lambda_sim > 1 && R$lambda_sim > 1)) {
			if(verbose) {
					cat('\r', paste("sm grad.:", round(100*(i/length(sms))), "% complete"))
					flush.console()
				}
			next
		} else {
			# if they do, titrate to estimate threshold
			titrateVals[2,]  <-  c(L$lambda_sim, R$lambda_sim)
			titrateDelta     <-  1

			while(!all(titrateVals[2,1] > 1 && titrateVals[2,2] > 1) & titrateDelta > res) {
				sf_prop       <-  midPoint(titrateVals[1,1],titrateVals[1,2])
				titrateDelta  <-  sf_prop - min(titrateVals[1,1], titrateVals[1,2])
			
				proposal  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
											 hf = hf, hm = hm, sf = sf_prop, sm = sms[i], C = C, delta = delta, 
											 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
											 tlimit = tlimit, eqThreshold = eqThreshold, Ainvade = Ainvade, intInit = TRUE)
				# replace appropraite boundary value
				if(proposal$lambda_sim  < 1) {
					titrateVals[1,][titrateVals[2,] < 1]  <-  sf_prop
					titrateVals[2,][titrateVals[2,] < 1]  <-  proposal$lambda_sim
				}
				if(proposal$lambda_sim  > 1) {
					titrateVals[1,][titrateVals[2,] > 1]  <-  sf_prop
					titrateVals[2,][titrateVals[2,] > 1]  <-  proposal$lambda_sim
				}
				if(verbose) {
					cat('\r', paste("sm grad.:", round(100*(i/length(sms))), "% complete; Titrate ratio = ", round(titrateDelta/res, digits=2), "(done when < 1)"))
					flush.console()
				}

			}

		}
		# Record sf threshold value
		sfThreshold[i]  <-  sf_prop
	}

	# compile results as data frame
	results.df  <-  as.data.frame(cbind(sfs, 
										smThreshold, 
										sms, 
										sfThreshold
										)
								 )
	colnames(results.df)  <-  c("sf",
								"smThreshold",
								"sms",
								"sfThreshold"
								)

	# export data as .csv to ./output/data
	if(writeFile) {
			filename <-  paste("./output/simData/extThreshold_SfxSm", "_sMax", sMax, "_res", res, "_hf", hf, "_hm", hm, 
							"_C", C, "_delta", delta, "_dj", delta_j, "_da", delta_a, "_dg", delta_gamma, "_f", theta[4], ".csv", sep="")
			write.csv(results.df, file=filename, row.names = FALSE)
	} else{
			return(results.df)
	}
}




###################################################
#  ANALYSES USING REAL DEMOGRAPHIC DATA
###################################################

#' Parameters:
#' datMat:	List of Demography matrices from COMPADRE
#' theta.list: 	list of demographic parameters from Peterson et al. (2016)
#' 			theta = list(D, G, F, O, A, S, R), where
#' 			D = seed bank survival (0.534)
#' 			G = seed germination rate (0.469)
#' 			F = flower production (0.64)
#' 			O = ovules per flower (614)
#' 			A = seedling recruits proportional to clonal rosette recruits (6.7e-4)
#' 			S = overwinter survival (0.179)
#' 			R = rosette production (8.71)
#' delta.list: 	list of demographic parameters from Peterson et al. (2016)
#' 			delta_D
#' 			delta_G
#' 			delta_F
#' 			delta_O
#' 			delta_S
#' hf,hm:	Dominance coefficient for SA fitness effects
#' sf,sm:	Selection coefficient for SA fitness effects
#' C:		Population selfing rate
#' tlimit:	Max # of generations for fwd simulation
#' eqThreshold:	threshold per gen. change in allele frequencies for determining demographic equilibrium
#' Ainvade:	Should A allele start off rare
#' intInit:	Optional initial frequency of A allele
#' STILL NOT WORKING
fwdSimMimulusDat  <-  function(datMat, theta.list, delta.list, useCompadre = TRUE,
								hf = 1/2, hm = 1/2, sf = 0.025, sm = 0.03, C = 0,
								tlimit = 10^4, eqThreshold=1e-9, Ainvade = FALSE, intInit = FALSE, ...) {

	# calculate # life stages from data matrices
	g   <-  3
	om  <-  ncol(datMat$A)

	# Create identity and ones matrices
	Ig   <-  diag(g)
	Iom  <-  diag(om)
	eg   <-  ones(c(g,1))
	eom  <-  ones(c(om,1))
	K    <-  vecperm(om,g)

	# Fitness Expressions
	fii        <-  c(1, 1 - hf*sf, 1 - sf)
	fii_prime  <-  c(1 - sm, 1 - hm*sm, 1)

	# unpack demographic parameters
	D  <-  theta.list$D
	G  <-  theta.list$G
	F  <-  theta.list$F
	O  <-  theta.list$O
	A  <-  theta.list$A
	S  <-  theta.list$S
	R  <-  theta.list$R

	# unpack inbreeding depression parameters
	delta_D  <-  delta.list$delta_D
	delta_G  <-  delta.list$delta_G
	delta_F  <-  delta.list$delta_F
	delta_O  <-  delta.list$delta_O
	delta_S  <-  delta.list$delta_S

	if(useCompadre) {
	# unpack inbreeding depression parameters
	delta_O  <-  delta.list$delta_O
	delta_S  <-  delta.list$delta_S

	# create delta matrix 
	deltaMat  <-  matrix(rbind( c(1, (1 - delta_O), (1 - delta_O)),
								c(1, (1 - delta_O), (1 - delta_O)),
								c(1, (1 - delta_S), (1 - delta_S))), nrow=om, ncol=om)
	datMatF_S  <-  datMat$F*deltaMat
	datMatF_X  <-  datMat$F
	datMatU_S  <-  datMat$U*deltaMat
	datMatU_X  <-  datMat$U
	}

	if(!useCompadre) {
	parDatMatS  <-  list(
					 U = matrix(rbind(c(D*(1 - delta_D)*(1 - (G*(1 - delta_G))), 0, 0),
									  c(D*(1 - delta_D)*G*(1 - delta_G), F*(1 - delta_F)*O*(1 - delta_O)*A*G*(1 - delta_G), F*(1 - delta_F)*O*(1 - delta_O)*A*G*(1 - delta_G)),
									  c(0, S*(1 - delta_S)*R*(1 - delta_S), S*(1 - delta_S)*R*(1 - delta_S))), nrow=om, ncol=om),
					 F = matrix(rbind(c(0, F*(1 - delta_F)*O*(1 - delta_O)*A*(1 - G*(1 - delta_G)), F*(1 - delta_F)*O*(1 - delta_O)*A*(1 - G*(1 - delta_G))),
									  c(0, 0, 0),
									  c(0, 0, 0)), nrow=om, ncol=om)
					 )
	parDatMatX  <-  list(
					 U = matrix(rbind(c(D*(1 - G), 0, 0),
									  c(D*G, F*O*A*G, F*O*A*G),
									  c(0, S*R, S*R)), nrow=om, ncol=om),
					 F = matrix(rbind(c(0, F*O*A*(1 - G), F*O*A*(1 - G)),
									  c(0, 0, 0),
									  c(0, 0, 0)), nrow=om, ncol=om)
					 )
	}

	####################################################
	# Population and female fertility perameters
	####################################################
	USi       <-  zeros(c(om,om,g))
	UXi       <-  zeros(c(om,om,g))
	FXi       <-  zeros(c(om,om,g))
	FSi       <-  zeros(c(om,om,g))
	FXi_pr    <-  zeros(c(om,om,g))
	lambda_i  <-  zeros(g)

	# create genotype-specific survival and fertility submatrices
	for (i in 1:g){
	    
	    if(useCompadre) {
			USi[,,i]     <-  datMatU_S
			UXi[,,i]     <-  datMatU_X
			FSi[,,i]     <-  datMatF_S
			FXi[,,i]     <-  datMatF_X
			FXi_pr[,,i]  <-  datMatF_X
	    }
		if(!useCompadre) {
	    	USi[,,i]     <-  parDatMatS$U
	    	UXi[,,i]     <-  parDatMatX$U
	    	FSi[,,i]     <-  parDatMatS$F
	    	FXi[,,i]     <-  parDatMatS$F
	    	FXi_pr[,,i]  <-  parDatMatX$F
		}

		FSi[,,i][FSi[,,i] != 0]        <-  C*(FSi[,,i][FSi[,,i] != 0])*fii[i]
		FXi[,,i][FXi[,,i] != 0]        <-  (1 - C)*(FXi[,,i][FXi[,,i] != 0])*fii[i]
		FXi_pr[,,i][FXi_pr[,,i] != 0]  <-  (FXi_pr[,,i][FXi_pr[,,i] != 0])*fii_prime[i]

		# Calculate gentoype-specific pop. growth rates
		lambda_i[i]  <-  max(Re(eigen(FSi[,,i] + USi[,,i] + FXi[,,i], symmetric=FALSE, only.values = TRUE)$values))
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
		blkFS        <-  blkFS + kronecker(emat(g,g,i,i),FSi[,,i])
		blkFX        <-  blkFX + kronecker(emat(g,g,i,i),FXi[,,i])
		blkFX_prime  <-  blkFX_prime + kronecker(emat(g,g,i,i),FXi_pr[,,i])
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

	###############################################################
	# Simulating the Stage X genotype dynamics - generation loop
	###############################################################

	# Set initial state vectors for each boundary (fixed for AA & aa)
	n0AA  <-  c(      C*c(c(100-(om-1), rep(1,times=(om-1))), rep(0,times = 2*om)),
				(1 - C)*c(c(100-(om-1), rep(1,times=(om-1))), rep(0,times = 2*om)))
	n0aa  <-  c(      C*c(rep(0,times = 2*om), c(100-(om-1), rep(1,times=(om-1)))), 
				(1 - C)*c(rep(0,times = 2*om), c(100-(om-1), rep(1,times=(om-1)))))

	# Simulate to demographic equilibrium for each boundary
	AAEq  <-  fwdDyn2Eq(nzero=n0AA, om=om, g=g, W_prime=W_prime, blkFX_prime=blkFX_prime, blkUS=blkUS, blkUX=blkUX, W=W, Iom=Iom, Ig=Ig, K=K, Z=Z, blkFS=blkFS, blkFX=blkFX, tlimit=10^3, eqThreshold=eqThreshold)
	aaEq  <-  fwdDyn2Eq(nzero=n0aa, om=om, g=g, W_prime=W_prime, blkFX_prime=blkFX_prime, blkUS=blkUS, blkUX=blkUX, W=W, Iom=Iom, Ig=Ig, K=K, Z=Z, blkFS=blkFS, blkFX=blkFX, tlimit=10^3, eqThreshold=eqThreshold)


	# use intermediate frequency to find Eq.?
	if(intInit) {
		initEq  <-  c(      C*c(rep(c(32,rep(4/3, times=om-1)), times=g)),
					  (1 - C)*c(rep(c(32,rep(4/3, times=om-1)), times=g)))
		} else{
			# pick boundary for invasion based on major allele frequency		
			if(Ainvade) {
				initEq  <-  aaEq$p*10
			}
			if(!Ainvade) {
				initEq  <-  AAEq$p*10
			}
			# Introduce rare allele in heterozygotes
			initEq[(om+1)]       <-  0.1
			initEq[(om*g+om+1)]  <-  0.1
 		}

	# Simulate to demographic equilibrium with rare allele
	invadeEq  <-  fwdDyn2Eq(nzero=initEq, om=om, g=g, W_prime=W_prime, blkFX_prime=blkFX_prime, blkUS=blkUS, blkUX=blkUX, W=W, Iom=Iom, Ig=Ig, K=K, Z=Z, blkFS=blkFS, blkFX=blkFX, tlimit=tlimit, eqThreshold=eqThreshold)

	# Calculate coexistence conditions based on leading eigenvalue of the Jacobian
	# NOTE: we rearrange the order of the pHat_i values so they match the structure
	# of the jacobian, which is ordered by genotype then self/outcross
	nHat_AA  <-  AAEq$n[c(1:om, g*om+1:om, om+1:om, g*om+om+1:om, 2*om+1:om, g*om+2*om+1:om)]
	nHat_aa  <-  aaEq$n[c(1:om, g*om+1:om, om+1:om, g*om+om+1:om, 2*om+1:om, g*om+2*om+1:om)]
	pHat_AA  <-  nHat_AA/norm(as.matrix(nHat_AA), type="1")
	pHat_aa  <-  nHat_aa/norm(as.matrix(nHat_aa), type="1")
	zeta_i   <-  calcZeta(om=om, FSi=FSi, FXi=FXi, FXi_pr=FXi_pr, USi=USi, UXi=UXi,
						 pHat_AA=pHat_AA, pHat_aa=pHat_aa, C=C, delta=delta_O, 
						 lambda_i=lambda_i)

	##################
	# results
	res  <-  list(
					"nout"         =  invadeEq$nout,
					"p_genotypes"  =  invadeEq$p_genotypes,
					"pEq"          =  invadeEq$p_genotypes[,ncol(invadeEq$p_genotypes)],
					"eqReached"    =  invadeEq$eqReached,
					"lambda_i"     =  lambda_i,
					"zeta_i"       =  zeta_i,
					"om"           =  om,
					"g"            =  g,
					"hf"           =  hf,
					"hm"           =  hm,
					"sf"           =  sf,
					"sm"           =  sm,
					"C"            =  C,
					"delta.list"   =  delta.list,
					"tlimit"       =  tlimit,
					"Ainvade"      =  Ainvade,
					"intInit"      =  intInit,
					"extinct"      =  invadeEq$lambda_sim < 1,
					"polymorphism" =  !any(round(invadeEq$p_genotypes[,ncol(invadeEq$p_genotypes)], digits=3) == 1)
					)
	return(res)

}





##############################
#' Loop over selection space 
#' Using Compadre Data
#'
#' Parameters:
#' dims: 	vector c(om, g), with om = # stages, g = # genotypes
#' theta: 	vector of length 4, c(sigma_J, sigma_A, gamma, f_ii)
#' selPars:	vector of length 4, c(hf, hm, sf, sm,)
selLoopDatMat  <-  function(sMax = 0.15, nSamples = 1e+2,
							datMat, theta.list, delta.list, useCompadre = FALSE,
							hf = 1/2, hm = 1/2, C = 0, 
							tlimit = 10^5, intInit = FALSE, eqThreshold=1e-9, 
							writeFile=TRUE, progBar = TRUE, ...) {
	
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


			results  <-  fwdSimMimulusDat(datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = useCompadre,
										hf = hf, hm = hm, sf = sfs[i], sm = sms[i], C = C,
										tlimit = tlimit, eqThreshold = eqThreshold, Ainvade = Ainvade, intInit = intInit)

			intInit  <-  FALSE
			extinct[i]       <-  results$extinct
			polymorphism[i]  <-  results$polymorphism
			pEq[i,]          <-  results$pEq
			zeta_i[i,]       <-  results$zeta_i
			lambda_i[i,]     <-  results$lambda_i

		if(progBar){
			cat('\r', paste(100*(i/nSamples),'% Complete'))
		}

	}

	# compile results as data frame
	hfVec       <-  rep(hf, times=nSamples)
	hmVec       <-  rep(hm, times=nSamples)
	CVec        <-  rep(C, times=nSamples)
	results.df  <-  as.data.frame(cbind(sfs, 
										sms, 
										extinct, 
										polymorphism, 
										pEq, 
										zeta_i,
										lambda_i,
										hfVec, 
										hmVec, 
										CVec
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
								"C"
								)

	# export data as .csv to ./output/data
	if(writeFile) {
			filename <-  paste("./output/simData/Mimulus_SfxSm", "_sMax", sMax, "_nSamples", nSamples, 
							   "_hf", hf, "_hm", hm, "_C", C, "_Compadre", useCompadre, ".csv", sep="")
			write.csv(results.df, file=filename, row.names = FALSE)
	} else{
			return(results.df)
	}

}








##############################
#' Identify extinction threshold across sf x sm parameter space
#' faster method for calculating proportions of sf x sm where
#' different dynamical outcomes happen(?)
#'
extinctThreshTitrateMimulus  <-  function(sMax = 0.15, res=0.001, 
										  datMat, theta.list, delta.list, useCompadre = FALSE,
										  hf = 1/2, hm = 1/2, C = 0, 
										  tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
										  writeFile=TRUE, verbose=TRUE) {
	
	# vector of sf values to titrate along
	sfs          <-  seq(res,(sMax-res),by=res)
	sms          <-  sfs
	smThreshold  <-  rep(NA, times=length(sfs))
	sfThreshold  <-  smThreshold

	# Three values we need to keep track of
	titrateStartVals  <-  matrix(rbind(c(sfs[1], sfs[length(sfs)]),
								  c(NA,NA)),nrow=2,ncol=2)
	rownames(titrateStartVals)  <-  c("sVals", "extinct")
	colnames(titrateStartVals)  <-  c("left", "right")

	# loop over sf values
	for(i in 1:length(sfs)) {

		titrateVals  <-  titrateStartVals

		# start at boundaries of [0, sMax]
		L  <-  fwdSimMimulusDat(datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = TRUE,
								hf = hf, hm = hm, sf = sfs[i], sm = titrateVals[1,1], C = C,
								tlimit = tlimit, eqThreshold=eqThreshold, Ainvade = FALSE, intInit = intInit)
		R  <-  fwdSimMimulusDat(datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = TRUE,
								hf = hf, hm = hm, sf = sfs[i], sm = titrateVals[1,2], C = C,
								tlimit = tlimit, eqThreshold=eqThreshold, Ainvade = FALSE, intInit = intInit)
		
		# check to see if extinction outcomes differ. If not, go to next sf value
		if(all(L$lambda_sim > 1 && R$lambda_sim > 1)) {
			if(verbose) {
					cat('\r', paste("sf grad.:", round(100*(i/length(sfs))), "% complete"))
					flush.console()
				}
			next
		} else {
			# if they do, titrate to estimate threshold
			titrateVals[2,]  <-  c(L$lambda_sim, R$lambda_sim)
			titrateDelta     <-  1

			while(!all(titrateVals[2,1] > 1 && titrateVals[2,2] > 1) & titrateDelta > res) {
				sm_prop       <-  midPoint(titrateVals[1,1],titrateVals[1,2])
				titrateDelta  <-  sm_prop - min(titrateVals[1,1], titrateVals[1,2])
			
				proposal  <-  fwdSimMimulusDat(datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = TRUE,
											   hf = hf, hm = hm, sf = sfs[i], sm = sm_prop, C = C,
											   tlimit = tlimit, eqThreshold=eqThreshold, Ainvade = FALSE, intInit = intInit)

				# replace appropraite boundary value
				if(proposal$lambda_sim  < 1) {
					titrateVals[1,][titrateVals[2,] < 1]  <-  sm_prop
					titrateVals[2,][titrateVals[2,] < 1]  <-  proposal$lambda_sim
				}
				if(proposal$lambda_sim  > 1) {
					titrateVals[1,][titrateVals[2,] > 1]  <-  sm_prop
					titrateVals[2,][titrateVals[2,] > 1]  <-  proposal$lambda_sim
				}
				if(verbose) {
					cat('\r', paste("sf grad.:", round(100*(i/length(sfs))), "% complete; Titrate ratio = ", round(titrateDelta/res, digits=2), "(done when < 1)"))
					flush.console()
				}
			}
		}
		# Record sm threshold value
		smThreshold[i]  <-  sm_prop
	}

	# loop over sm values
	for(i in 1:length(sms)) {

		titrateVals  <-  titrateStartVals

		# start at boundaries of [0, sMax]
		L  <-  fwdSimMimulusDat(datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = TRUE,
								hf = hf, hm = hm, sf = titrateVals[1,1], sm = sms[i], C = C,
								tlimit = tlimit, eqThreshold=eqThreshold, Ainvade = FALSE, intInit = intInit)
		R  <-  fwdSimMimulusDat(datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = TRUE,
								hf = hf, hm = hm, sf = titrateVals[1,2], sm = sms[i], C = C,
								tlimit = tlimit, eqThreshold=eqThreshold, Ainvade = FALSE, intInit = intInit)
		
		# check to see if extinction outcomes differ. If not, go to next sf value
		if(all(L$lambda_sim > 1 && R$lambda_sim > 1)) {
			if(verbose) {
					cat('\r', paste("sm grad.:", round(100*(i/length(sms))), "% complete"))
					flush.console()
				}
			next
		} else {
			# if they do, titrate to estimate threshold
			titrateVals[2,]  <-  c(L$lambda_sim, R$lambda_sim)
			titrateDelta     <-  1

			while(!all(titrateVals[2,1] > 1 && titrateVals[2,2] > 1) & titrateDelta > res) {
				sf_prop       <-  midPoint(titrateVals[1,1],titrateVals[1,2])
				titrateDelta  <-  sf_prop - min(titrateVals[1,1], titrateVals[1,2])
			
				proposal  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
											 hf = hf, hm = hm, sf = sf_prop, sm = sms[i], C = C, delta = delta, 
											 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
											 tlimit = tlimit, eqThreshold = eqThreshold, Ainvade = Ainvade, intInit = TRUE)
				# replace appropraite boundary value
				if(proposal$lambda_sim  < 1) {
					titrateVals[1,][titrateVals[2,] < 1]  <-  sf_prop
					titrateVals[2,][titrateVals[2,] < 1]  <-  proposal$lambda_sim
				}
				if(proposal$lambda_sim  > 1) {
					titrateVals[1,][titrateVals[2,] > 1]  <-  sf_prop
					titrateVals[2,][titrateVals[2,] > 1]  <-  proposal$lambda_sim
				}
				if(verbose) {
					cat('\r', paste("sm grad.:", round(100*(i/length(sms))), "% complete; Titrate ratio = ", round(titrateDelta/res, digits=2), "(done when < 1)"))
					flush.console()
				}

			}

		}
		# Record sf threshold value
		sfThreshold[i]  <-  sf_prop
	}

	# compile results as data frame
	results.df  <-  as.data.frame(cbind(sfs, 
										smThreshold, 
										sms, 
										sfThreshold 
										)
								 )
	colnames(results.df)  <-  c("sf",
								"smThreshold",
								"sms",
								"sfThreshold"
								)

	# export data as .csv to ./output/data
	if(writeFile) {
			filename <-  paste("./output/simData/extThreshold_SfxSm", "_sMax", sMax, "_res", res, "_hf", hf, "_hm", hm, 
							"_C", C, "_delta", delta, "_dj", delta_j, "_da", delta_a, "_dg", delta_gamma, "_f", theta[4], ".csv", sep="")
			write.csv(results.df, file=filename, row.names = FALSE)
	} else{
			return(results.df)
	}
}

