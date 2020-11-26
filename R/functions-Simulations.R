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
	if(C == 0){
		M22_AA  <-  UXi[,,2] + (1/2)*FXi[,,2] + kronecker(c((1/(2*pn_AA)))*FXi[,,1]%*%(pHat_AA[1:om] + pHat_AA[om+1:om]), (ones(om)%*%FXi_pr[,,2]))
		M22_aa  <-  UXi[,,2] + (1/2)*FXi[,,2] + kronecker(c((1/(2*pn_aa)))*FXi[,,3]%*%(pHat_aa[(4*om+1):(5*om)] + pHat_aa[(5*om+1):(6*om)]), (ones(om)%*%FXi_pr[,,2]))
		zeta_i  <-  c(max(Re(eigen(M22_AA, symmetric=FALSE, only.values = TRUE)$values)),
						 max(Re(eigen(M22_aa, symmetric=FALSE, only.values = TRUE)$values)))
		zeta_i  <-  zeta_i/lambda_i[c(1,3)]
	}
	if(C > 0){
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
	}
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
		lambda_i[i]  <-  max(Re(eigen(C*(1 - delta)*(USi[,,i]) + FSi[,,i] + (1 - C)*UXi[,,i] + FXi[,,i], symmetric=FALSE, only.values = TRUE)$values))
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

pickInvadingAllele  <-  function(hf, hm, sf, sm, C) {
		if(hf == hm && hf == 1/2) {
			qHat  <-  popGen_qHat_Add(sf = sf, sm = sm, C = C)
		}
		if(hf == hm && hf < 1/2) {
			qHat  <-  popGen_qHat_DomRev(h = hf, sf = sf, sm = sm, C = C)
		}
			if(qHat < 1/2){
				Ainvade  <-  TRUE
			}
			if(qHat > 1/2){
				Ainvade  <-  FALSE
			}
	Ainvade
}

midPoint  <-  function(x1,x2) {
	x1 + (x2-x1)/2 + runif(1, min=-1, max=1)*(x2-x1)/10
}


##############################
#' Identify extinction threshold across sf x sm parameter space
#' faster method for calculating proportions of sf x sm where
#' different dynamical outcomes happen(?)
#'
extinctThreshTitrate  <-  function(sMax=0.15, res=0.015, precision = 1e-4,
								   om = 2, g = 3, theta = c(0.6,0.6,0.05,6), theta_prime = 6, 
								   hf = 1/2, hm = 1/2, C = 0, delta = 0, 
								   delta_j = 0, delta_a = 0, delta_gamma = 0,
								   tlimit = 10^5, Ainvade=FALSE, eqThreshold = 1e-8, 
								   makePlots=FALSE, verbose=TRUE, writeFile=TRUE) {
	
	# vector of sf values to titrate along
	sfs          <-  seq(0, sMax, by=res)
	sms          <-  sfs
	smThreshold  <-  rep(NA, times=length(sfs))
	sfThreshold  <-  smThreshold

	# Three values we need to keep track of
	titrateStartVals  <-  matrix(rbind(c(sfs[1], sfs[length(sfs)]),
								  c(NA,NA)),nrow=2,ncol=2)
	rownames(titrateStartVals)  <-  c("sVals", "extinct")
	colnames(titrateStartVals)  <-  c("left", "right")

if(makePlots) {
	dev.new()
	plot(NA, xlim=c(0,sMax), ylim=c(0,sMax))
}
	# loop over sf values
	for(i in 1:length(sfs)) {

		titrateVals   <-  titrateStartVals
		titrateDelta  <-  1

		# start tit. at boundaries of [0, sMax] OR at narrower window based on previous threshold
		L  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
							  hf = hf, hm = hm, sf = sfs[i], sm = titrateVals[1,1], C = C, delta = delta, 
							  delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
							  tlimit = tlimit, eqThreshold = eqThreshold, Ainvade = FALSE, intInit = FALSE)
		R  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
							  hf = hf, hm = hm, sf = sfs[i], sm = titrateVals[1,2], C = C, delta = delta, 
							  delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
							  tlimit = tlimit, eqThreshold = eqThreshold, Ainvade = TRUE, intInit = FALSE)
		
		# If populations are viable at both boundaries of sf x sm space, 
		# go to next sf value
		if(L$lambda_sim > 1 && R$lambda_sim > 1) {
			if(verbose) {
					cat('\r', paste("Titrating sf value ", i, "/", length(sfs)))
					flush.console()
				}
			next
		} 
		# If extinction occurs at both boundaries of sf x sm space,
		# set sm extinction threshold to 0
		if(L$lambda_sim < 1 && R$lambda_sim < 1) {
			smThreshold[i]  <-  0
if(makePlots) {
	points(sfs[i] ~ 0)
}
			next
		}
		
		# if extinction outcomes differ, titrate to find threshold
		titrateVals[2,]  <-  c(L$lambda_sim, R$lambda_sim)
		errCount  <-  0
		while(titrateDelta > precision && errCount < 21) {

			# Propose new tit. boundary midway between previous boundary values
			sm_prop       <-  midPoint(titrateVals[1,1],titrateVals[1,2])
			titrateDelta  <-  sm_prop - min(titrateVals[1,1], titrateVals[1,2])
			errCount  <-  errCount + 1
if(makePlots) {
	points(sfs[i] ~ sm_prop, cex=0.25)
}
			proposal  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
										 hf = hf, hm = hm, sf = sfs[i], sm = sm_prop, C = C, delta = delta, 
										 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
										 tlimit = tlimit, eqThreshold = eqThreshold, Ainvade = FALSE, intInit = TRUE)
			# replace appropriate boundary value
			if(proposal$lambda_sim  < 1) {
				titrateVals[1,][titrateVals[2,] < 1]  <-  sm_prop
				titrateVals[2,][titrateVals[2,] < 1]  <-  proposal$lambda_sim
			}
			if(proposal$lambda_sim  > 1) {
				titrateVals[1,][titrateVals[2,] > 1]  <-  sm_prop
				titrateVals[2,][titrateVals[2,] > 1]  <-  proposal$lambda_sim
			}
			if(verbose) {
				cat('\r', paste("Titrating sf value ", i, "/", length(sfs), "; Tit. ratio = ", round(titrateDelta/precision, digits=2), "(done when < 1)"))
				flush.console()
			}

		}

		# Record sm threshold value
		smThreshold[i]  <-  sm_prop
if(makePlots) {
	points(sfs[i] ~ smThreshold[i])
}
	}

	# Rename titrate Columns for top/bottom titration
	colnames(titrateStartVals)  <-  c("bott.", "top")

	# loop over sm values (top/bottom titration)
	for(i in 1:length(sms)) {

		titrateVals   <-  titrateStartVals
		titrateDelta  <-  1

		# if sms value is smaller than the smallest smThreshold value found
		# by L/R titration, go to next sm tit. value 
		if(!all(is.na(smThreshold)) && sms[i] < min(na.omit(smThreshold))) {
			if(verbose) {
					cat('\r', paste("Titrating sm value ", i, "/", length(sfs)))
					flush.console()
				}
			next
		}

		# if the tit. sm value is larger than the smallest smThreshold
		# found by L/R titration, narrow the titration window using 
		if(!all(is.na(smThreshold)) && sms[i] > min(na.omit(smThreshold))) {
			smDev  <-  smThreshold[!is.na(smThreshold)] - sms[i]
			titrateVals[1,]  <-  c( 0, 
									min(sfs[!is.na(smThreshold)][ smDev == max(smDev[ smDev < 0])], sMax))
		} 

		# start at boundaries of [0, sMax] OR at narrower window based on previous threshold
		B  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
								 hf = hf, hm = hm, sf = titrateVals[1,1], sm = sms[i], C = C, delta = delta, 
								 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
								 tlimit = tlimit, eqThreshold = eqThreshold, Ainvade = FALSE, intInit = TRUE)
		T  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
								 hf = hf, hm = hm, sf = titrateVals[1,2], sm = sms[i], C = C, delta = delta, 
								 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
								 tlimit = tlimit, eqThreshold = eqThreshold, Ainvade = FALSE, intInit = TRUE)
		
		# If extinction does not occur at either boundary of sf x sm space,
		# go to next sm tit. value
		if(B$lambda_sim > 1 && T$lambda_sim > 1) {
			if(verbose) {
					cat('\r', paste("Titrating sm value ", i, "/", length(sfs)))
					flush.console()
				}
			next
		}
		# If extinction occurs at both T & B,
		# set sf extinction threshold to 0
		if(B$lambda_sim < 1 && T$lambda_sim < 1) {
			sfThreshold[i]  <-  0
if(makePlots) {
	points(0 ~ sms[i])
}
			next
		}

		# Assign lambda vals. to tit. values
		titrateVals[2,]  <-  c(B$lambda_sim, T$lambda_sim)
		errCount  <-  0
		# Titrate!
		while(titrateDelta > precision && errCount < 21) {

			# If extinction outcome differs, propose new boundary
			# midway between previous boundary values
			sf_prop       <-  midPoint(titrateVals[1,1],titrateVals[1,2])
			errCount  <-  errCount + 1
			titrateDelta  <-  sf_prop - min(titrateVals[1,1], titrateVals[1,2])
if(makePlots) {
	points(sf_prop ~ sms[i], cex=0.5)
}
			proposal  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
										 hf = hf, hm = hm, sf = sf_prop, sm = sms[i], C = C, delta = delta, 
										 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
										 tlimit = tlimit, eqThreshold = eqThreshold, Ainvade = Ainvade, intInit = FALSE)
			
			# replace appropriate boundary value
			if(proposal$lambda_sim  < 1) {
				titrateVals[1,][titrateVals[2,] < 1]  <-  sf_prop
				titrateVals[2,][titrateVals[2,] < 1]  <-  proposal$lambda_sim
			}
			if(proposal$lambda_sim  > 1) {
				titrateVals[1,][titrateVals[2,] > 1]  <-  sf_prop
				titrateVals[2,][titrateVals[2,] > 1]  <-  proposal$lambda_sim
			}
			if(verbose) {
				cat('\r', paste("Titrating sf value ", i, "/", length(sfs), "; Tit. ratio = ", round(titrateDelta/precision, digits=2), "(done when < 1)"))
				flush.console()
			}

		}

		# Record sf threshold value
		sfThreshold[i]  <-  sf_prop
if(makePlots) {
	points(sfThreshold[i] ~ sms[i], pch=2)
}
	}

	# clean up any smThreshold values that fall below the largest sfThreshold value
	if(!all(is.na(sfThreshold))) {
		smThreshold[sfs < max(na.omit(sfThreshold))]  <-  NA
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






##############################
#' Identify extinction threshold across sf x sm parameter space
#' faster method for calculating proportions of sf x sm where
#' different dynamical outcomes happen(?)
#'
extinctThreshTitrateRotate  <-  function(sMax=0.15, res=0.0015, precision=1e-4,
								   om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 6, 
								   hf = 1/2, hm = 1/2, C = 0, delta = 0, 
								   delta_j = 0, delta_a = 0, delta_gamma = 0,
								   tlimit = 10^5, Ainvade=FALSE, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE) {
	
	# vector of sf values to titrate along
	sfs          <-  seq(0, sMax, by=res)
	sms          <-  sfs
	extinctionThreshold  <-  cbind(rep(NA, times=2*length(sfs)),
								   rep(NA, times=2*length(sfs)))

	# Three values we need to keep track of
	titrateVals  <-  matrix(rbind(c(0, 0, NA),
								  c(0, 0, NA)), nrow=2,ncol=3)
	rownames(titrateVals)  <-  c("left", "right")
	colnames(titrateVals)  <-  c("sm", "sf", "lambda")

dev.new()
plot(NA, xlim=c(0,sMax), ylim=c(0,sMax))
#plot(NA, xlim=c(0,1), ylim=c(0,1))

	# loop over titration start values for upper-left triangle
	# of sf x sm space
	for(i in 1:(length(sfs))) {

		# Step through starting point on upper-left edges of sf x sm space
		titrateVals[1,c(1:2)]  <-  c(0, rev(sfs)[i])
		titrateVals[2,c(1:2)]  <-  c(sms[i], sMax)
		titrateVals[,3]        <-  c(NA,NA)

points(	titrateVals[,2] ~ 	titrateVals[,1])
		# calculate population growth rate at these boundaries
		L  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
								 hf = hf, hm = hm, sf = titrateVals[1,2], sm = titrateVals[1,1], C = C, delta = delta, 
								 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
								 tlimit = 10^5, eqThreshold = eqThreshold, Ainvade = FALSE, intInit = TRUE)
		R  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
								 hf = hf, hm = hm, sf = titrateVals[2,2], sm = titrateVals[2,1], C = C, delta = delta, 
								 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
								 tlimit = 10^5, eqThreshold = eqThreshold, Ainvade = FALSE, intInit = TRUE)
		
		# check to see if extinction outcomes differ. If not, go to next sf value
		if(all(L$lambda_sim > 1 && R$lambda_sim > 1)) {
			if(verbose) {
					cat('\r', paste("Titration", round(100*(i/(2*length(sfs)))), "% complete"))
					flush.console()
				}
			next
		} 
		if(!all(L$lambda_sim > 1 && R$lambda_sim > 1)) {

			# if they do, titrate to find threshold
			titrateVals[,3]  <-  c(L$lambda_sim, R$lambda_sim)
			titrateDelta     <-  eucDist(titrateVals[1,c(1:2)], titrateVals[2,c(1:2)])

			# If extinction occurs at both edges of sf x sm space,
			# set extinction threshold to x-axis
			if(all(titrateVals[,3] < 1)) {
				extinctionThreshold[i,]  <-  c(titrateVals[2,1], 0)
				points(0 ~ titrateVals[2,1], pch=3, cex=2)
				next
			}

			while(!all(titrateVals[,3] > 1) && titrateDelta > precision) {
				s_prop       <-  c(midPoint(titrateVals[1,1],titrateVals[2,1]),
								   midPoint(titrateVals[1,2],titrateVals[2,2]))
# points(s_prop[2]~s_prop[2])
				titrateDelta  <-  eucDist(s_prop, titrateVals[1,c(1:2)])
				proposal  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
											 hf = hf, hm = hm, sf = s_prop[2], sm = s_prop[1], C = C, delta = delta, 
											 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
											 tlimit = 10^5, eqThreshold = eqThreshold, Ainvade = Ainvade, intInit = TRUE)
				# replace appropriate boundary value
				if(proposal$lambda_sim  < 1) {
					titrateVals[titrateVals[,3] < 1, ]  <-  c(s_prop, proposal$lambda_sim)
				}
				if(proposal$lambda_sim  > 1) {
					titrateVals[titrateVals[,3] > 1, ]  <-  c(s_prop, proposal$lambda_sim)
				}
				if(verbose) {
					cat('\r', paste("Titration", round(100*(i/(2*length(sfs)))), "% complete; Tit. ratio = ", round(titrateDelta/precision, digits=2), "(done when < 1)"))
					flush.console()
				}

			}
		# Record sm threshold value
		extinctionThreshold[i,]  <-  s_prop
		points(extinctionThreshold[i,2] ~ extinctionThreshold[i,1], pch=3, cex=2)

		}

	}

	# loop over titration start values for lower-right triangle
	# of sf x sm space
	for(i in 1:(length(sfs)-1)) {

		# Step through starting point on upper-left edges of sf x sm space
		titrateVals[1,c(1:2)]  <-  c(sfs[i+1], 0)
		titrateVals[2,c(1:2)]  <-  c(sMax, sms[length(sms)-i])
		titrateVals[,3]        <-  c(NA,NA)

points(	titrateVals[,2] ~ 	titrateVals[,1])
		# calculate population growth rate at these boundaries
		L  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
								 hf = hf, hm = hm, sf = titrateVals[1,2], sm = titrateVals[1,1], C = C, delta = delta, 
								 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
								 tlimit = 10^5, eqThreshold = eqThreshold, Ainvade = FALSE, intInit = TRUE)
		R  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
								 hf = hf, hm = hm, sf = titrateVals[2,2], sm = titrateVals[2,1], C = C, delta = delta, 
								 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
								 tlimit = 10^5, eqThreshold = eqThreshold, Ainvade = FALSE, intInit = TRUE)

		# check to see if extinction outcomes differ. If not, go to next sf value
		if(all(L$lambda_sim > 1 && R$lambda_sim > 1)) {
			if(verbose) {
					cat('\r', paste("Titration", round(100*((length(sfs)+i)/(2*length(sfs)))), "% complete"))
					flush.console()
				}
			next
		} 
		if(!all(L$lambda_sim > 1 && R$lambda_sim > 1)) {
			
			# if they do, titrate to find threshold
			titrateVals[,3]  <-  c(L$lambda_sim, R$lambda_sim)
			titrateDelta     <-  eucDist(titrateVals[1,c(1:2)], titrateVals[2,c(1:2)])

			# If extinction occurs at both edges of sf x sm space,
			# set extinction threshold to x-axis
			if(all(titrateVals[,3] < 1)) {
				extinctionThreshold[i,]  <-  c(titrateVals[2,1], 0)
				points(0 ~ titrateVals[2,1], pch=3, cex=2)
				next
			}

			while(!all(titrateVals[,3] > 1) && titrateDelta > precision) {
				s_prop       <-  c(midPoint(titrateVals[1,1],titrateVals[2,1]),
								   midPoint(titrateVals[1,2],titrateVals[2,2]))
				titrateDelta  <-  eucDist(s_prop, titrateVals[1,c(1:2)])
				proposal  <-  fwdDemModelSim(om = om, g = g, theta = theta, theta_prime = theta_prime, 
											 hf = hf, hm = hm, sf = s_prop[2], sm = s_prop[1], C = C, delta = delta, 
											 delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
											 tlimit = 10^5, eqThreshold = eqThreshold, Ainvade = Ainvade, intInit = TRUE)
				# replace appropriate boundary value
				if(proposal$lambda_sim  < 1) {
					titrateVals[titrateVals[,3] < 1, ]  <-  c(s_prop, proposal$lambda_sim)
				}
				if(proposal$lambda_sim  > 1) {
					titrateVals[titrateVals[,3] > 1, ]  <-  c(s_prop, proposal$lambda_sim)
				}
				if(verbose) {
					cat('\r', paste("Titration", round(100*((length(sfs)+i)/(2*length(sfs)))), "% complete; Tit. ratio = ", round(titrateDelta/precision, digits=2), "(done when < 1)"))
					flush.console()
				}

			}

		# Record sm threshold value
		extinctionThreshold[(length(sfs) + i),]  <-  s_prop
		points(extinctionThreshold[(length(sfs) + i),2] ~ extinctionThreshold[(length(sfs) + i),1], pch=3, cex=2)

		}

	}

#plot(NA, xlim=c(0,sMax), ylim=c(0,sMax))
#points(extinctionThreshold[,2] ~ extinctionThreshold[,1], pch=3, cex=2)

	# compile results as data frame
	results.df  <-  as.data.frame(extinctionThreshold)
	colnames(results.df)  <-  c("sm",
								"sf"
								)

	# export data as .csv to ./output/data
	if(writeFile) {
			filename <-  paste("./output/simData/extThresholdRotate_SfxSm", "_sMax", sMax, "_res", res, "_hf", hf, "_hm", hm, 
							"_C", C, "_delta", delta, "_dj", delta_j, "_da", delta_a, "_dg", delta_gamma, "_f", theta[4], ".csv", sep="")
			write.csv(results.df, file=filename, row.names = FALSE)
	} else{
			return(results.df)
	}
}


#####################################################
###########################
#' Find Boundary Eq. & Evaluate Zeta_i
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
findInvBoundZeta  <-  function(om = 2, g = 3, theta = c(0.6, 0.6, 0.05, 6), theta_prime = 5.9, 
							hf = 1/2, hm = 1/2, sf = 0.1, sm = 0.105, C = 0, delta = 0, 
							delta_j = 0, delta_a = 0, delta_gamma = 0, datMat = NA,
							tlimit = 10^2, eqThreshold = 1e-8, ...) {

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
		lambda_i[i]  <-  max(Re(eigen(C*(1 - delta)*(USi[,,i]) + FSi[,,i] + (1 - C)*UXi[,,i] + FXi[,,i], symmetric=FALSE, only.values = TRUE)$values))
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

	# Set initial state vectors for each boundary (fixed for AA & aa)
	n0AA  <-  c(      C*c(c(100-(om-1), rep(1,times=(om-1))), rep(0,times = 2*om)),
				(1 - C)*c(c(100-(om-1), rep(1,times=(om-1))), rep(0,times = 2*om)))
	n0aa  <-  c(      C*c(rep(0,times = 2*om), c(100-(om-1), rep(1,times=(om-1)))), 
				(1 - C)*c(rep(0,times = 2*om), c(100-(om-1), rep(1,times=(om-1)))))
#browser()
	# Simulate to demographic equilibrium for each boundary
	AAEq  <-  fwdDyn2Eq(nzero=n0AA, om=om, g=g, W_prime=W_prime, blkFX_prime=blkFX_prime, blkUS=blkUS, blkUX=blkUX, W=W, Iom=Iom, Ig=Ig, K=K, Z=Z, blkFS=blkFS, blkFX=blkFX, tlimit=10^2, eqThreshold=eqThreshold)
	aaEq  <-  fwdDyn2Eq(nzero=n0aa, om=om, g=g, W_prime=W_prime, blkFX_prime=blkFX_prime, blkUS=blkUS, blkUX=blkUX, W=W, Iom=Iom, Ig=Ig, K=K, Z=Z, blkFS=blkFS, blkFX=blkFX, tlimit=10^2, eqThreshold=eqThreshold)

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
					"lambda_i"     =  lambda_i,
					"zeta_i"       =  zeta_i
					)
	return(res)

}






##############################
#' Identify extinction threshold across sf x sm parameter space
#' faster method for calculating proportions of sf x sm where
#' different dynamical outcomes happen(?)
#'
titrateInvBoundaries  <-  function(sMax=0.15, res=0.0015, precision=1e-4,
								   om = 2, g = 3, theta = c(0.6,0.6,0.05,5.8), theta_prime = 6, 
								   hf = 1/2, hm = 1/2, C = 0, delta = 0, 
								   delta_j = 0, delta_a = 0, delta_gamma = 0,
								   tlimit = 10^2, eqThreshold = 1e-8, verbose=TRUE, writeFile=TRUE) {
	
	# vector of sf values to titrate along
	sms          <-  seq(res, sMax, by=res)
	aInvBound    <-  rep(NA, times=length(sms))
	AInvBound    <-  rep(NA, times=length(sms))

	## Titrate to find invasion boundary for a allele 
	# Three values we need to keep track of
	titrateVals  <-  matrix(rbind(c(0, 0, NA),
								  c(0, 0, NA)), nrow=2,ncol=3)
	rownames(titrateVals)  <-  c("left", "right")
	colnames(titrateVals)  <-  c("sm", "sf", "Zeta_AA")

# plot(NA, xlim=c(0,sMax), ylim=c(0,sMax))

	# loop over sf values, and titrate for sm value corresponding to
	# a allele invasion condition 
	for(i in 1:(length(sms))) {

		# Step through starting point on upper-left edges of sf x sm space
		titrateVals[1,c(1:2)]  <-  c(sms[i], sMax)
		titrateVals[2,c(1:2)]  <-  c(sms[i], 0)
		titrateVals[,3]        <-  c(NA,NA)

# points(	titrateVals[,2] ~ 	titrateVals[,1])
		# calculate population growth rate at these boundaries
		T  <-  findInvBoundZeta(om = om, g = g, theta = theta, theta_prime = theta_prime,
							hf = hf, hm = hm, sf = titrateVals[1,2], sm = titrateVals[1,1], C = C, delta = delta, 
							delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
							tlimit = tlimit, eqThreshold = eqThreshold)
		B  <-  findInvBoundZeta(om = om, g = g, theta = theta, theta_prime = theta_prime,
							hf = hf, hm = hm, sf = titrateVals[2,2], sm = titrateVals[2,1], C = C, delta = delta, 
							delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
							tlimit = tlimit, eqThreshold = eqThreshold)

		# Titrate to find threshold value where zeta_AA crosses 1
		titrateVals[,3]  <-  c(T$zeta_i[1], B$zeta_i[1])
		titrateDelta     <-  eucDist(titrateVals[1,c(1:2)], titrateVals[2,c(1:2)])

		# If aBound is > sMax, record sMax as top of inv. space, skip to next
		# sm value in loop
		if(all(titrateVals[,3] > 1)) {
			aInvBound[i]  <-  sMax
#			points(aInvBound[i] ~ sms[i], cex=1)
			next
		}

		while(titrateDelta > precision) {
			s_prop       <-  c(midPoint(titrateVals[1,1],titrateVals[2,1]),
							   midPoint(titrateVals[1,2],titrateVals[2,2]))

			proposal  <-  findInvBoundZeta(om = om, g = g, theta = theta, theta_prime = theta_prime,
							hf = hf, hm = hm, sf = s_prop[2], sm = s_prop[1], C = C, delta = delta, 
							delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
							tlimit = tlimit, eqThreshold = eqThreshold)

			titrateDelta  <-  eucDist(s_prop, titrateVals[1,c(1:2)])

			# replace appropriate boundary value
			if(proposal$zeta_i[1]  < 1) {
				titrateVals[titrateVals[,3] < 1, ]  <-  c(s_prop, proposal$zeta_i[1])
			}
			if(proposal$zeta_i[1]  > 1) {
				titrateVals[titrateVals[,3] > 1, ]  <-  c(s_prop, proposal$zeta_i[1])
			}
			if(proposal$zeta_i[1]  == 1) {
				titrateDelta  <-  0
			}
			if(verbose) {
				cat('\r', paste("Titration", round(100*(i/(2*length(sms)))), "% complete"))
				flush.console()
			}
		}

		# Record sm threshold value
		aInvBound[i]  <-  s_prop[2]
#		points(aInvBound[i] ~ sms[i], cex=1)
	}

	## Titrate to find invasion boundary for A allele 
	# Three values we need to keep track of
	titrateVals  <-  matrix(rbind(c(0, 0, NA),
								  c(0, 0, NA)), nrow=2,ncol=3)
	rownames(titrateVals)  <-  c("top", "bottom")
	colnames(titrateVals)  <-  c("sm", "sf", "Zeta_aa")


	# loop over sm values, and titrate for sf value corresponding to
	# a allele invasion condition 
	for(i in 1:(length(sms))) {

		# Step through starting point on upper-left edges of sf x sm space
		titrateVals[1,c(1:2)]  <-  c(sms[i], sMax)
		titrateVals[2,c(1:2)]  <-  c(sms[i], 0)
		titrateVals[,3]        <-  c(NA,NA)

# points(	titrateVals[,2] ~ 	titrateVals[,1])
		# calculate population growth rate at these boundaries
		T  <-  findInvBoundZeta(om = om, g = g, theta = theta, theta_prime = theta_prime,
							hf = hf, hm = hm, sf = titrateVals[1,2], sm = titrateVals[1,1], C = C, delta = delta, 
							delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
							tlimit = tlimit, eqThreshold = eqThreshold)
		B  <-  findInvBoundZeta(om = om, g = g, theta = theta, theta_prime = theta_prime,
							hf = hf, hm = hm, sf = titrateVals[2,2], sm = titrateVals[2,1], C = C, delta = delta, 
							delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
							tlimit = tlimit, eqThreshold = eqThreshold)
		
		# Titrate to find threshold value where zeta_AA crosses 1
		titrateVals[,3]  <-  c(T$zeta_i[2], B$zeta_i[2])
		titrateDelta     <-  eucDist(titrateVals[1,c(1:2)], titrateVals[2,c(1:2)])

		# If ABound is > sMax, record 0 as bottom of inv. space, skip to next
		# sm value in loop
		if(all(titrateVals[,3] > 1)) {
			AInvBound[i]  <-  0
#			points(AInvBound[i] ~ sms[i], cex=1)
			next
		}

		while(titrateDelta > precision) {
			s_prop       <-  c(midPoint(titrateVals[1,1],titrateVals[2,1]),
							   midPoint(titrateVals[1,2],titrateVals[2,2]))

			proposal  <-  findInvBoundZeta(om = om, g = g, theta = theta, theta_prime = theta_prime,
							hf = hf, hm = hm, sf = s_prop[2], sm = s_prop[1], C = C, delta = delta, 
							delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
							tlimit = tlimit, eqThreshold = eqThreshold)

			titrateDelta  <-  eucDist(s_prop, titrateVals[1,c(1:2)])

			# replace appropriate boundary value
			if(proposal$zeta_i[2]  < 1) {
				titrateVals[titrateVals[,3] < 1, ]  <-  c(s_prop, proposal$zeta_i[2])
			}
			if(proposal$zeta_i[2]  > 1) {
				titrateVals[titrateVals[,3] > 1, ]  <-  c(s_prop, proposal$zeta_i[2])
			}
			if(proposal$zeta_i[2]  == 1) {
				titrateDelta  <-  0
			}
			if(verbose) {
				cat('\r', paste("Titration", round(100*((length(sms)+i)/(2*length(sms)))), "% complete"))
				flush.console()
			}
		}

		# Record sm threshold value
		AInvBound[i]  <-  s_prop[2]
#		points(AInvBound[i] ~ sms[i], cex=1)
	}


	# compile results as data frame
	results.df  <-  as.data.frame(rbind(c(0,0,0),cbind(sms, aInvBound, AInvBound)))

	# export data as .csv to ./output/data
	if(writeFile) {
			filename <-  paste("./output/simData/invasionBoundaries", "_sMax", sMax, "_res", res, "_hf", hf, "_hm", hm, 
							"_C", C, "_delta", delta, "_dj", delta_j, "_da", delta_a, "_dg", delta_gamma, "_f", theta[4], ".csv", sep="")
			write.csv(results.df, file=filename, row.names = FALSE)
	} else{
			return(results.df)
	}
}



###################################################
# makeData functions using titration to quantify
# demographically viable polymorphic parameter 
# space 
###################################################

##############################
#' function to make data to quantify demographically viable 
#' polymorphic parameter space for different dominance, scenarios
#' selfing rates, and fertility values.
makeDataPolyParamSpace  <-  function(sMax=0.1, res=0.01, precision = 1e-4,
									 om = 2, g = 3, theta = c(0.6,0.6,0.05,NA), theta_prime = NA, 
									 hVals= c(1/2, 1/4), fVals = c(5.8, 5.9, 6.0, 6.5), 
									 delta = 0, delta_j = 0, delta_a = 0, delta_gamma = 0,
									 tlimit = 10^5, eqThreshold = 1e-8) {

	# Vector of selfing rates, storage dataframe
	Cs       <-  seq(0,0.9, by=0.1)
	dataSet  <-  rep(NA, times=10)

	# Setup progress bar
	print('Making Data For Fig.2')
	pb   <-  txtProgressBar(min=0, max=(length(hVals)*length(fVals)*length(Cs)), style=3)
	setTxtProgressBar(pb, 0)

	# Loop over dominance, fertility, selfing
	for(i in 1:length(hVals)) {
		for(j in 1:length(fVals)) {

			# assign fertility value to theta & theta_prime
			theta[4]     <-  fVals[j]
			theta_prime  <-  fVals[j]

			for(k in 1:length(Cs)) {

				# Find Invasion Boundaries
				invData  <-  titrateInvBoundaries(sMax=sMax, res=res, precision=precision,
												  om = om, g = g, theta = theta, theta_prime = theta_prime, 
												  hf = hVals[i], hm = hVals[i], C = Cs[k], 
												  delta = delta, delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
												  tlimit = tlimit, eqThreshold = eqThreshold, verbose=FALSE, writeFile=FALSE)

				# Find Extinction Threshold
				extData  <-  extinctThreshTitrate(sMax=sMax, res=res, precision=precision,
												  om = om, g = g, theta = theta, theta_prime = theta_prime, 
												  hf = hVals[i], hm = hVals[i], C = Cs[k], 
												  delta = delta, delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
												  tlimit = tlimit, eqThreshold = eqThreshold, 
												  makePlots=FALSE, verbose=FALSE, writeFile=FALSE)

				# Append new data to dataSet
				dataSet  <-  rbind(dataSet, cbind(rep(hVals[i], times=nrow(invData)),
												  rep(fVals[j], times=nrow(invData)),
												  rep(Cs[k], times=nrow(invData)),
												  invData, 
												  extData))

			# Update Progress Bar
			setTxtProgressBar(pb, (((i-1)*length(fVals)*length(Cs)) + (j-1)*length(Cs) + (k)))


			}
		}
	}	

	# trim first row of dataSet before saving file & assign column names
	dataSet  <-  dataSet[-1,]
	colnames(dataSet)  <-  c("h", "f", "C", "smInv", "aInv", "AInv", "sf", "smExt", "sm", "sfExt")

	# Export data as .csv to ./output/data
	filename <-  paste("./output/simData/dataPolySpaceFig", "_sMax", sMax, "_res", res, "_delta", delta, "_dj", delta_j, "_da", delta_a, "_dg", delta_gamma, ".csv", sep="")
	write.csv(dataSet, file=filename, row.names = FALSE)

}


##############################
#' For some reason, the titrations leave some aberrant
#' gaps in the sf_thresholds in the resulting data file.
#' This function fixes these gaps.
#' 
#' Parameters
#' df:	polySpace dataframe
cleanPolySpaceData  <-  function(data) {

    # Get plotting levels from df
    hLev  <-  unique(data$h)
    fLev  <-  unique(data$f)
    CLev  <-  unique(data$C)
    nHs   <-  length(hLev)
    nCs   <-  length(CLev)
    nFs   <-  length(fLev)
	newSfExt  <-  c()

	# loop through and fill NA gaps
	for(i in 1:nHs) {
        for(j in 1:nFs) {
			for(k in 1:nCs) {
				# Subset data by Dominance
				d  <-  data[data$h == hLev[i],]            
				# Subset data by Fertility Value
				d  <-  d[d$f == fLev[j],]     
				# Subset data by Selfing Rate
				d  <-  d[d$C == CLev[k],]     

				# loop over sfExt and identify NA gaps
				for(n in 1:length(d$sfExt)) {
					if(all(is.na(d$sfExt[1:n]))) {
						next
					}
					if(sum(!is.na(d$sfExt)) == 1) {
						next
					}
					
					if(is.na(d$sfExt[n])) {

						step  <- 1
						d$sfExt[c((n-1),(n+step))]
						while(any(is.na(d$sfExt[c((n-1),(n+step))])) && n+step <= nrow(d)) {
							step  <-  step + 1
						}
						gapEdges  <-  d$sfExt[c((n-1),(n+step))]
						gapFill  <-  rep((gapEdges[2] - gapEdges[1]), times=step)
						gapFill  <-  d$sfExt[n-1] + (gapFill / (seq(1:step) + 1))
						d$sfExt[n:(n+step-1)]  <-  gapFill
					}
				}

				newSfExt  <-  c(newSfExt, d$sfExt)
			}
		}
	}

	data$sfExt  <-  newSfExt
	return(data)
}



##############################
#' function to summarize polySpaceData to quantify 
#' demographically viable polymorphic parameter space
#' (i.e., hacky numerical integration of curves in 
#' SuppFig-polySpace.pdf file)
#' 
#' #' Parameters:
#' df:	data.frame object form dataPolySpaceFig_*.csv
quantPolySpace  <-  function(data, pars) {

    # Clean up aberrant smExt value
    data  <-  cleanPolySpaceData(data = data)
    data$smExt[90]  <-  mean(c(data$smExt[89], data$smExt[91]))

    # Extract plotting parameter values from df pars,
    # get levels for loops
    res      <-  pars$res
    sMax     <-  pars$sMax
    hLev     <-  unique(data$h)
    fLev     <-  unique(data$f)
    CLev     <-  unique(data$C)
    nHs      <-  length(hLev)
    nFs      <-  length(fLev)
    nCs      <-  length(CLev)
    results  <-  rep(0, times=9)

	# total area of sf x sm parameter space
	totSpace   <-  sum(rep(sMax*res, times=(sMax/res)))

    # loop through and subset data
	for(i in 1:nHs) {
		for(j in 1:nFs) {
			for(k in 1:nCs) {

				# Subset data by Dominance
				d  <-  data[data$h == hLev[i],]
				# Subset data by Fertility Value
				d  <-  d[d$f == fLev[j],]
				# Subset data by Selfing Rate
				d  <-  d[d$C == CLev[k],]

				# Calculate approx. total polymorphic space
				polySpace  <-  sum((d$aInv - d$AInv)*res)
				PrPoly     <-  polySpace/totSpace

				# Calculate approx. total extinction space
				if(length(d$smExt[!is.na(d$smExt)]) < 1) {
					extSpace   <-  0
					PrExt      <-  0
					viaPoly    <-  polySpace
					PrViaPoly  <-  PrPoly
				} else {
						extSm     <-  sum((sMax - d$smExt[!is.na(d$smExt)])*res)
						bottSm    <-  d$sf[!is.na(d$smExt)][1]
						extSf     <-  sum((bottSm - d$sfExt[!is.na(d$sfExt)])*res)
						extSpace  <-  extSm + extSf
						PrExt     <-  extSpace/totSpace

						# Calculate demographically viable polymorphic parameter space
						aFixExt    <-  (d$AInv - bottSm)*res
						aFixExt    <-  aFixExt[aFixExt > 0]
						aFixExt    <-  sum(aFixExt[!is.na(aFixExt)]) + extSf
						extPoly    <-  extSpace - aFixExt
						viaPoly    <-  polySpace - extPoly
						PrViaPoly  <-  viaPoly/totSpace
					}

				# append results
				results  <-  rbind(results, c(hLev[i], fLev[j], CLev[k], polySpace, PrPoly, extSpace, PrExt, viaPoly, PrViaPoly))

			}
		}
	}

	results  <-  results[-1,]
	colnames(results)  <-  c("h", "f", "C", "polySpace", "PrPoly", "extSpace", "PrExt", "viaPoly", "PrViaPoly")
	results  <-  as.data.frame(results, row.names=FALSE)

	return(results)
#	# Export data as .csv to ./output/data
#	filename <-  paste("./output/simData/dataPolySpaceFigTitrate", "_sMax", sMax, "_res", res, "_delta", pars$d, "_dj", pars$dj, "_da", pars$da, "_dg", pars$dg, ".csv", sep="")
#	write.csv(results, file=filename, row.names = FALSE)

}





##############################
#' function to make data to quantify demographically viable 
#' polymorphic parameter space for different dominance,
#' selfing rates, and inbreeding depression values.
makeDataDeltaPolyParamSpace  <-  function(sMax=0.15, res=0.0015, precision = 1e-4,
									 om = 2, g = 3, theta = c(0.6,0.6,0.05,6.5), theta_prime = 6.5, 
									 hVals= c(1/2, 1/4), dStar = 0.8, 
									 delta = 0, delta_j = 0, delta_a = 0, delta_gamma = 0,
									 tlimit = 10^5, eqThreshold = 1e-8) {


	# Storage for data
	dataSet  <-  rep(NA, times=11)

	print('Making Data For Fig.3')

	# Loop over dominance values
	for(i in 1:length(hVals)) {

		cat('\r', "hVal = ", hVals[i], '\n')

		Cs  <- c(seq(0, 0.5, by=0.025), 0.75, 0.9)
		
		# Define dStar ~ C
		deltaSeq  <-  predDelta(dStar=dStar, b=1/2, a=0.2, C=Cs)

		# Start progress bar
		print("Cs grad for delta")
		pb   <-  txtProgressBar(min=0, max=length(Cs), style=3)
		setTxtProgressBar(pb, 0)

		# loop over Cs for delta
		for(j in 1:length(Cs)){
			# Find Invasion Boundaries
			invData  <-  titrateInvBoundaries(sMax=sMax, res=res, precision=precision,
											  om = om, g = g, theta = theta, theta_prime = theta_prime, 
											  hf = hVals[i], hm = hVals[i], C = Cs[j], 
											  delta = deltaSeq[j], delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
											  tlimit = tlimit, eqThreshold = eqThreshold, verbose=FALSE, writeFile=FALSE)

				# Find Extinction Threshold
			extData  <-  extinctThreshTitrate(sMax=sMax, res=res, precision=precision,
											  om = om, g = g, theta = theta, theta_prime = theta_prime, 
											  hf = hVals[i], hm = hVals[i], C = Cs[j], 
											  delta = deltaSeq[j], delta_j = delta_j, delta_a = delta_a, delta_gamma = delta_gamma,
											  tlimit = tlimit, eqThreshold = eqThreshold, 
											  makePlots=FALSE, verbose = FALSE, writeFile=FALSE)

			# Append new data to dataSet
			dataSet  <-  rbind(dataSet, cbind(rep("delta", times=nrow(invData)),
											  rep(hVals[i], times=nrow(invData)),
											  rep(deltaSeq[j], times=nrow(invData)),
											  rep(Cs[j], times=nrow(invData)),
											  invData, 
											  extData))
			# Update Progress Bar
			setTxtProgressBar(pb, j)
		}

		# Start progress bar
		print("Cs grad for delta_j")
		pb   <-  txtProgressBar(min=0, max=length(Cs), style=3)
		setTxtProgressBar(pb, 0)

		# loop over Cs for delta_j
		for(j in 1:length(Cs)){
			# Find Invasion Boundaries
			invData  <-  titrateInvBoundaries(sMax=sMax, res=res, precision=precision,
											  om = om, g = g, theta = theta, theta_prime = theta_prime, 
											  hf = hVals[i], hm = hVals[i], C = Cs[j], 
											  delta = delta, delta_j = deltaSeq[j], delta_a = delta_a, delta_gamma = delta_gamma,
											  tlimit = tlimit, eqThreshold = eqThreshold, verbose=FALSE, writeFile=FALSE)

				# Find Extinction Threshold
			extData  <-  extinctThreshTitrate(sMax=sMax, res=res, precision=precision,
											  om = om, g = g, theta = theta, theta_prime = theta_prime, 
											  hf = hVals[i], hm = hVals[i], C = Cs[j], 
											  delta = delta, delta_j = deltaSeq[j], delta_a = delta_a, delta_gamma = delta_gamma,
											  tlimit = tlimit, eqThreshold = eqThreshold, verbose = FALSE, writeFile=FALSE)

			# Append new data to dataSet
			dataSet  <-  rbind(dataSet, cbind(rep("d_j", times=nrow(invData)),
											  rep(hVals[i], times=nrow(invData)),
											  rep(deltaSeq[j], times=nrow(invData)),
											  rep(Cs[j], times=nrow(invData)),
											  invData, 
											  extData))
			# Update Progress Bar
			setTxtProgressBar(pb, j)
		}

		# Start progress bar
		print("Cs grad for delta_a")
		pb   <-  txtProgressBar(min=0, max=length(Cs), style=3)
		setTxtProgressBar(pb, 0)

		# loop over Cs for delta_a
		for(j in 1:length(Cs)){
			# Find Invasion Boundaries
			invData  <-  titrateInvBoundaries(sMax=sMax, res=res, precision=precision,
											  om = om, g = g, theta = theta, theta_prime = theta_prime, 
											  hf = hVals[i], hm = hVals[i], C = Cs[j], 
											  delta = delta, delta_j = delta_j, delta_a = deltaSeq[j], delta_gamma = delta_gamma,
											  tlimit = tlimit, eqThreshold = eqThreshold, verbose=FALSE, writeFile=FALSE)

				# Find Extinction Threshold
			extData  <-  extinctThreshTitrate(sMax=sMax, res=res, precision=precision,
											  om = om, g = g, theta = theta, theta_prime = theta_prime, 
											  hf = hVals[i], hm = hVals[i], C = Cs[j], 
											  delta = delta, delta_j = delta_j, delta_a = deltaSeq[j], delta_gamma = delta_gamma,
											  tlimit = tlimit, eqThreshold = eqThreshold, verbose = FALSE, writeFile=FALSE)

			# Append new data to dataSet
			dataSet  <-  rbind(dataSet, cbind(rep("d_a", times=nrow(invData)),
											  rep(hVals[i], times=nrow(invData)),
											  rep(deltaSeq[j], times=nrow(invData)),
											  rep(Cs[j], times=nrow(invData)),
											  invData, 
											  extData))
			# Update Progress Bar
			setTxtProgressBar(pb, j)
		}

		# Start progress bar
		print("Cs grad for delta_g")
		pb   <-  txtProgressBar(min=0, max=length(Cs), style=3)
		setTxtProgressBar(pb, 0)

		# loop over Cs for delta_gamma
		for(j in 1:length(Cs)){
			# Find Invasion Boundaries
			invData  <-  titrateInvBoundaries(sMax=sMax, res=res, precision=precision,
											  om = om, g = g, theta = theta, theta_prime = theta_prime, 
											  hf = hVals[i], hm = hVals[i], C = Cs[j], 
											  delta = delta, delta_j = delta_j, delta_a = delta_a, delta_gamma = deltaSeq[j],
											  tlimit = tlimit, eqThreshold = eqThreshold, verbose=FALSE, writeFile=FALSE)

				# Find Extinction Threshold
			extData  <-  extinctThreshTitrate(sMax=sMax, res=res, precision=precision,
											  om = om, g = g, theta = theta, theta_prime = theta_prime, 
											  hf = hVals[i], hm = hVals[i], C = Cs[j], 
											  delta = delta, delta_j = delta_j, delta_a = delta_a, delta_gamma = deltaSeq[j],
											  tlimit = tlimit, eqThreshold = eqThreshold, verbose = FALSE, writeFile=FALSE)

			# Append new data to dataSet
			dataSet  <-  rbind(dataSet, cbind(rep("d_g", times=nrow(invData)),
											  rep(hVals[i], times=nrow(invData)),
											  rep(deltaSeq[j], times=nrow(invData)),
											  rep(Cs[j], times=nrow(invData)),
											  invData, 
											  extData))
			# Update Progress Bar
			setTxtProgressBar(pb, j)
		}

	}

	# trim first row of dataSet before saving file & assign column names
	dataSet  <-  dataSet[-1,]
	colnames(dataSet)  <-  c("Delta", "h", "deltaVal", "C", "smInv", "aInv", "AInv", "sf", "smExt", "sm", "sfExt")

	# Export data as .csv to ./output/data
	filename <-  paste("./output/simData/dataDeltaPolySpaceFig", "_sMax", sMax, "_res", res, "_dStar", dStar, "_f", theta[4], ".csv", sep="")
	write.csv(dataSet, file=filename, row.names = FALSE)

}




##############################
#' For some reason, the titrations leave some aberrant
#' gaps in the sf_thresholds in the resulting data file.
#' This function fixes these gaps.
#' 
#' Parameters
#' df:	dataDeltaPolySpace dataframe
cleanDeltaPolySpaceData  <-  function(data) {

    # Get plotting levels from df
    hLev  <-  unique(data$h)
    dLev  <-  unique(data$Delta)
    CLev  <-  unique(data$C)
    nHs   <-  length(hLev)
    nDs   <-  length(dLev)
    nCs   <-  length(CLev)
	newSfExt  <-  c()
	newaInv   <-  c()

	# loop through and fill NA gaps
	for(i in 1:nHs) {
        for(j in 1:nDs) {
			for(k in 1:nCs) {
				# Subset data by Dominance
				d  <-  data[data$h == hLev[i],]            
				# Subset data by Fertility Value
				d  <-  d[d$Delta == dLev[j],]     
				# Subset data by Selfing Rate
				d  <-  d[d$C == CLev[k],]     

				# loop over sfExt and identify NA gaps
				for(n in 1:length(d$sfExt)) {
					if(all(is.na(d$sfExt[1:n]))) {
						next
					}
					if(sum(!is.na(d$sfExt)) == 1) {
						next
					}
					
					if(is.na(d$sfExt[n])) {

						step  <- 1
						d$sfExt[c((n-1),(n+step))]
						while(any(is.na(d$sfExt[c((n-1),(n+step))])) && n+step <= nrow(d)) {
							step  <-  step + 1
						}
						gapEdges  <-  d$sfExt[c((n-1),(n+step))]
						gapFill  <-  rep((gapEdges[2] - gapEdges[1]), times=step)
						gapFill  <-  d$sfExt[n-1] + (gapFill / (seq(1:step) + 1))
						d$sfExt[n:(n+step-1)]  <-  gapFill
					}
				}
				newSfExt  <-  c(newSfExt, d$sfExt)

				# Check aInvade at C = 0
				if(j > 1){
					d$aInv[1]  <-  d$aInv[2]
				}
				newaInv  <-  c(newaInv, d$aInv)

			}
		}
	}

	data$sfExt  <-  newSfExt
	data$aInv   <-  newaInv
	return(data)
}






#   NEED TO FIX !!!!!
##############################
#' function to summarize polySpaceData to quantify 
#' demographically viable polymorphic parameter space
#' (i.e., hacky numerical integration of curves in 
#' SuppFig-polySpace.pdf file)
#' 
#' #' Parameters:
#' df:	data.frame object form dataDeltaPolySpaceFig_*.csv
quantDeltaPolySpace  <-  function(data, pars) {

    # Clean up aberrant smExt value
    data  <-  cleanDeltaPolySpaceData(data = data)
    data$smExt[6530:6544]  <-  0
    data$sfExt[6578:6579]  <-  0

    # Extract plotting parameter values from df pars,
    # get levels for loops
    res      <-  pars$res
    sMax     <-  pars$sMax
    hLev     <-  unique(data$h)
    dLev     <-  unique(data$Delta)
    CLev     <-  unique(data$C)
    nHs      <-  length(hLev)
    nDs      <-  length(dLev)
    nCs      <-  length(CLev)
    results  <-  rep(0, times=9)

	# total area of sf x sm parameter space
	totSpace   <-  sum(rep(sMax*res, times=(sMax/res)))

    # loop through and subset data
	for(i in 1:nHs) {
		for(j in 1:nDs) {
			for(k in 1:nCs) {

				# Subset data by Dominance
				d  <-  data[data$h == hLev[i],]
				# Subset data by Inbreeding Depression Parameter
				d  <-  d[d$Delta == dLev[j],]
				# Subset data by Selfing Rate
				d  <-  d[d$C == CLev[k],]

				# Calculate approx. total polymorphic space
				polySpace  <-  sum((d$aInv - d$AInv)*res)
				if(polySpace > totSpace) {
					polySpace  <-  totSpace
				}
				PrPoly     <-  polySpace/totSpace

				# Calculate approx. total extinction space
				if(length(d$smExt[!is.na(d$smExt)]) < 1) {
					extSpace   <-  0
					PrExt      <-  0
					viaPoly    <-  polySpace
					PrViaPoly  <-  PrPoly
					results    <-  rbind(results, c(hLev[i], dLev[j], CLev[k], polySpace, PrPoly, extSpace, PrExt, viaPoly, PrViaPoly))
					next
				} 
				if(all(d$smExt[!is.na(d$smExt)]  ==  0)) {
					extSpace   <-  1
					PrExt      <-  1
					viaPoly    <-  0
					PrViaPoly  <-  0
					results    <-  rbind(results, c(hLev[i], dLev[j], CLev[k], polySpace, PrPoly, extSpace, PrExt, viaPoly, PrViaPoly))
					next
				}
				extSm     <-  sum((sMax - d$smExt[!is.na(d$smExt)])*res)
				bottSm    <-  d$sf[!is.na(d$smExt)][1]
				extSf     <-  sum((bottSm - d$sfExt[!is.na(d$sfExt)])*res)
				extSpace  <-  extSm + extSf
				PrExt     <-  extSpace/totSpace

				# Calculate demographically viable polymorphic parameter space
				aFixExt    <-  (d$AInv - bottSm)*res
				aFixExt    <-  aFixExt[aFixExt > 0]
				aFixExt    <-  sum(aFixExt[!is.na(aFixExt)]) + extSf
				extPoly    <-  extSpace - aFixExt
				viaPoly    <-  polySpace - extPoly
				PrViaPoly  <-  viaPoly/totSpace

				# append results
				results  <-  rbind(results, c(hLev[i], dLev[j], CLev[k], polySpace, PrPoly, extSpace, PrExt, viaPoly, PrViaPoly))
			}
		}
	}

	results  <-  results[-1,]
	colnames(results)  <-  c("h", "Delta", "C", "polySpace", "PrPoly", "extSpace", "PrExt", "viaPoly", "PrViaPoly")
	results  <-  as.data.frame(results, row.names=FALSE)

	return(results)
#	# Export data as .csv to ./output/data
#	filename <-  paste("./output/simData/dataPolySpaceFigTitrate", "_sMax", sMax, "_res", res, "_delta", pars$d, "_dj", pars$dj, "_da", pars$da, "_dg", pars$dg, ".csv", sep="")
#	write.csv(results, file=filename, row.names = FALSE)

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
		lambda_i[i]  <-  max(Re(eigen(C*(1 - delta)*(USi[,,i]) + FSi[,,i] + (1 - C)*UXi[,,i] + FXi[,,i], symmetric=FALSE, only.values = TRUE)$values))
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
					"lambda_sim"   =  invadeEq$lambda_sim,
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
extinctThreshMimulus  <-  function(sMax = 0.15, res=0.01, precision=1e-4, 
								   datMat, theta.list, delta.list, useCompadre = FALSE,
								   hf = 1/2, hm = 1/2, C = 0, 
								   tlimit = 10^5, intInit = TRUE, Ainvade=FALSE, eqThreshold=1e-9, 
								   writeFile=FALSE, verbose=TRUE) {
	
	
	# vector of sf values to titrate along
	sfs          <-  seq(0, sMax, by=res)
	sms          <-  sfs
	extinctionThreshold  <-  cbind(rep(NA, times=2*length(sfs)),
								   rep(NA, times=2*length(sfs)))

	# Three values we need to keep track of
	titrateVals  <-  matrix(rbind(c(0, 0, NA),
								  c(0, 0, NA)), nrow=2,ncol=3)
	rownames(titrateVals)  <-  c("left", "right")
	colnames(titrateVals)  <-  c("sm", "sf", "lambda")

plot(NA, xlim=c(0,sMax), ylim=c(0,sMax))

	# loop over titration start values for upper-left triangle
	# of sf x sm space
	for(i in 1:(length(sfs))) {

		# Step through starting point on upper-left edges of sf x sm space
		titrateVals[1,c(1:2)]  <-  c(0, sfs[length(sfs)+1-i])
		titrateVals[2,c(1:2)]  <-  c(sms[i], sMax)
		titrateVals[,3]        <-  c(NA,NA)

points(	titrateVals[,2] ~ 	titrateVals[,1])

		# calculate population growth rate at these boundaries

		# start at boundaries of [0, sMax]
		L  <-  fwdSimMimulusDat(datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = TRUE,
								hf = hf, hm = hm, sf = titrateVals[1,2], sm = titrateVals[1,1], C = C,
								tlimit = tlimit, eqThreshold=eqThreshold, Ainvade = FALSE, intInit = intInit)
		R  <-  fwdSimMimulusDat(datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = TRUE,
								hf = hf, hm = hm, sf = titrateVals[2,2], sm = titrateVals[2,1], C = C,
								tlimit = tlimit, eqThreshold=eqThreshold, Ainvade = FALSE, intInit = intInit)
		
		# check to see if extinction outcomes differ. If not, go to next sf value
		if(all(L$lambda_sim > 1 && R$lambda_sim > 1)) {
			if(verbose) {
					cat('\r', paste("Titration", round(100*(i/(2*length(sfs)))), "% complete"))
					flush.console()
				}
			next
		} 
		if(!all(L$lambda_sim > 1 && R$lambda_sim > 1)) {

			# if they do, titrate to find threshold
			titrateVals[,3]  <-  c(L$lambda_sim, R$lambda_sim)
			titrateDelta     <-  eucDist(titrateVals[1,c(1:2)], titrateVals[2,c(1:2)])

			while(!all(titrateVals[,3] > 1) && titrateDelta > precision) {
				s_prop       <-  c(midPoint(titrateVals[1,1],titrateVals[2,1]),
								   midPoint(titrateVals[1,2],titrateVals[2,2]))

				titrateDelta  <-  eucDist(s_prop, titrateVals[1,c(1:2)])
				proposal      <-  fwdSimMimulusDat(datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = TRUE,
												   hf = hf, hm = hm, sf = s_prop[2], sm = s_prop[1], C = C,
												   tlimit = tlimit, eqThreshold=eqThreshold, Ainvade = FALSE, intInit = TRUE)

				# replace appropriate boundary value
				if(proposal$lambda_sim  < 1) {
					titrateVals[titrateVals[,3] < 1, ]  <-  c(s_prop, proposal$lambda_sim)
				}
				if(proposal$lambda_sim  > 1) {
					titrateVals[titrateVals[,3] > 1, ]  <-  c(s_prop, proposal$lambda_sim)
				}
				if(verbose) {
					cat('\r', paste("Titration", round(100*(i/(2*length(sfs)))), "% complete; Tit. ratio = ", round(titrateDelta/res, digits=2), "(done when < 1)"))
					flush.console()
				}

			}
		# Record sm threshold value
		extinctionThreshold[i,]  <-  s_prop
		points(extinctionThreshold[i,2] ~ extinctionThreshold[i,1], pch=3, cex=2)

		}

	}


	# loop over titration start values for lower-right triangle
	# of sf x sm space
	for(i in 1:(length(sfs))) {

		# Step through starting point on upper-left edges of sf x sm space
		titrateVals[1,c(1:2)]  <-  c(sfs[i], 0)
		titrateVals[2,c(1:2)]  <-  c(sMax, sms[length(sms)+1-i])
		titrateVals[,3]        <-  c(NA,NA)

points(	titrateVals[,2] ~ 	titrateVals[,1])

		# calculate population growth rate at these boundaries
		L  <-  fwdSimMimulusDat(datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = TRUE,
								hf = hf, hm = hm, sf = titrateVals[1,2], sm = titrateVals[1,1], C = C,
								tlimit = tlimit, eqThreshold=eqThreshold, Ainvade = FALSE, intInit = intInit)
		R  <-  fwdSimMimulusDat(datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = TRUE,
								hf = hf, hm = hm, sf = titrateVals[2,2], sm = titrateVals[2,1], C = C,
								tlimit = tlimit, eqThreshold=eqThreshold, Ainvade = FALSE, intInit = intInit)

		# check to see if extinction outcomes differ. If not, go to next sf value
		if(all(L$lambda_sim > 1 && R$lambda_sim > 1)) {
			if(verbose) {
					cat('\r', paste("Titration", round(100*((length(sfs)+i)/(2*length(sfs)))), "% complete"))
					flush.console()
				}
			next
		} 
		if(!all(L$lambda_sim > 1 && R$lambda_sim > 1)) {
			
			# if they do, titrate to find threshold
			titrateVals[,3]  <-  c(L$lambda_sim, R$lambda_sim)
			titrateDelta     <-  eucDist(titrateVals[1,c(1:2)], titrateVals[2,c(1:2)])

			while(!all(titrateVals[,3] > 1) && titrateDelta > precision) {
				s_prop       <-  c(midPoint(titrateVals[1,1],titrateVals[2,1]),
								   midPoint(titrateVals[1,2],titrateVals[2,2]))
				titrateDelta  <-  eucDist(s_prop, titrateVals[1,c(1:2)])
				proposal      <-  fwdSimMimulusDat(datMat=datMat, theta.list=theta.list, delta.list=delta.list, useCompadre = TRUE,
												   hf = hf, hm = hm, sf = s_prop[2], sm = s_prop[1], C = C,
												   tlimit = tlimit, eqThreshold=eqThreshold, Ainvade = FALSE, intInit = TRUE)
				# replace appropriate boundary value
				if(proposal$lambda_sim  < 1) {
					titrateVals[titrateVals[,3] < 1, ]  <-  c(s_prop, proposal$lambda_sim)
				}
				if(proposal$lambda_sim  > 1) {
					titrateVals[titrateVals[,3] > 1, ]  <-  c(s_prop, proposal$lambda_sim)
				}
				if(verbose) {
					cat('\r', paste("Titration", round(100*((length(sfs)+i)/(2*length(sfs)))), "% complete; Tit. ratio = ", round(titrateDelta/precision, digits=2), "(done when < 1)"))
					flush.console()
				}

			}

		# Record sm threshold value
		extinctionThreshold[(length(sfs) + i),]  <-  s_prop
		points(extinctionThreshold[(length(sfs) + i),2] ~ extinctionThreshold[(length(sfs) + i),1], pch=3, cex=2)

		}

	}

	# compile results as data frame
	results.df  <-  as.data.frame(extinctionThreshold)
	colnames(results.df)  <-  c("sm",
								"sf"
								)
	# Alter file name for if there is inbreeding depression
	if (all(delta.list == 0)){
		ID  <-  "FALSE"
	} else(ID  <-  "TRUE")

	# export data as .csv to ./output/data
	if(writeFile) {
			filename <-  paste("./output/simData/extThresholdMimulus_SfxSm", "_sMax", sMax, "_res", res, "_hf", hf, "_hm", hm, 
							"_C", C, "_useCompadre", useCompadre, "_ID", ID, ".csv", sep="")
			write.csv(results.df, file=filename, row.names = FALSE)
	} else{
			return(results.df)
	}
}








