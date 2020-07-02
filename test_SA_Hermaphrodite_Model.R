# First attempt to implement a Mendelian Matrix Population Model
# for a partially selfing population of hermaphrodites with selection
# through both sex functions.
#
#	For simplicity I'm just hacking Lotte's code here to see if
#	all the component matrices fit, and if the model does what 
#   I expect it to do.
#
# The model has:
# 	Two stages: juvenile and adult
# 	Three genotypes: AA, Aa, aa
#   A constant selfing rate: C_AA = C_Aa = C_aa = C, where 0 <= C <= 1
#   Constant inbreeding depression: 0 < delta < 1


rm(list=ls())
#####################
##  Dependencies
source('R/functions-MatModels.R')

Ainvade  <-  TRUE
intInit  <-  FALSE

# Simulation time"
tlimit <- 10^4

# dimensions
om  <-  2
g   <-  3

# Create identity and ones matrices
Ig   <-  diag(g)
Iom  <-  diag(om)
eg   <-  ones(c(g,1))
eom  <-  ones(c(om,1))
K    <-  vecperm(om,g)

# Population Selfing Rate and inbreeding depression
C      <-  1/4
delta  <-  0

#####################################################"

# PARAMETERS AND INITIAL CONDITIONS 

####################################################"

# BASELINE PARAMETERS for survival and 
# fertility through Female and Male function
# theta=[sigma_J, sigma_A, gamma, f_ii]
theta_temp  <-  c(0.6,0.6,0.05,6.8)
theta       <-  rep.col(theta_temp,3)

theta_primetemp <- c(0.6,0.6,0.05,6.8)
theta_prime <- rep.col(theta_primetemp,3)

## SELECTION PARAMETERS
hf  <-  1/2
hm  <-  1/2
sf  <-  0.1
sm  <-  0.12
fii        <-  c(1, 1 - hf*sf, 1 - sf)
fii_prime  <-  c(1 - sm, 1 - hm*sm, 1)

##SELECTION DIFFERENTIAL FEMALE-FUNCTION
##female function
theta[4,]  <-  theta[4,]*fii
#theta[4,]  <-  (1 - beta)*theta[4,] + beta*(theta[4,]*fii)

#SELECTION DIFFERENTIAL MALE-FUNCTION
# # f' differential
theta_prime[4,]  <-  theta_prime[4,]*fii_prime
#theta_prime[4,]  <-  (1 - beta)*theta_prime[4,] + beta*(theta_prime[4,]*fii_prime)

#Initial conditions
# nzero  <-  [juv_AA, ad_AA, juv_Aa, ad_Aa, juv_aa, ad_aa]
# All aa
if(Ainvade) {
	nzero  <- round(c(C*(1 - delta)*100*c(0,0,0,0,141.3,18.2),(1 - C)*100*c(0,0,0,0,141.3,18.2)))
} 
#All AA
if(!Ainvade) {
	nzero  <- round(c(C*(1 - delta)*100*c(141.3,18.2,0,0,0,0),(1 - C)*100*c(141.3,18.2,0,0,0,0)))
}
#All Aa
if(intInit) {
	nzero  <-  round(c(C*(1 - delta)*100*c(0,0,141.3,18.2,0,0),(1 - C)*100*c(0,0,141.3,18.2,0,0)))
}

#####################################################

# END OF PARAMETERS and INITIAL CONDITIONS

#####################################################

# Population and female fertility perameters
sigma_J    <-  theta[1,]
sigma_A    <-  theta[2,]
gamma      <-  theta[3,]
f          <-  theta[4,]
f_prime    <-  theta_prime[4,]
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

for (i in 1:3){
    USi[,,i]       <- rbind(c(sigma_J[i]*(1 - gamma[i]), 0         ),
                            c(sigma_J[i]*gamma[i],       sigma_A[i]))
    UXi[,,i]       <- rbind(c(sigma_J[i]*(1 - gamma[i]), 0         ),
                            c(sigma_J[i]*gamma[i],       sigma_A[i]))
    FSi[,,i]       <- rbind(c(0,C*f[i]*(1 - delta)),
					        c(0,0))
    FXi[,,i]       <- rbind(c(0,(1 - C)*f[i]),
					        c(0,0))
    FXi_prime[,,i] <- rbind(c(0,(1 - C)*f_prime[i]),
							c(0,0))
    
    # genotype-specific population growth rate
    eigs_temp[i,]  <-  eigen((FSi[,,i] + USi[,,i]) + (FXi[,,i] + UXi[,,i]), symmetric=FALSE, only.values = TRUE)$values
    lambda[i]      <-  max(eigs_temp)
    Nfun[,,i]      <-  solve((diag(2) - (C*(1 - delta)*USi[,,i] + (1 - C)*UXi[,,i])), diag(2))
    R[,,i]         <-  (FSi[,,i] + FXi[,,i]) %*% Nfun[,,i]
    #genotype-specific R0
    R_0[i] <- max(eigen(R[,,i],symmetric=FALSE, only.values = TRUE)$values)
}
# for i

Atilde_genotype  <-  zeros(c(2*om,2*om,g))
lambda_full      <-  zeros(c(1,g))
for (i in 1:3){
	Atilde_genotype[,,i]  <- rbind(cbind(USi[,,i] + FSi[,,i],zeros(c(2,2))), 
								   cbind(FXi[,,i], UXi[,,i])
								   )
    lambda_full[i]    <- max(eigen(Atilde_genotype[,,i],symmetric=FALSE, only.values = TRUE)$values)
}


d <- diag(g)

#CREATE BLOCK DIAGONAL MATRICES
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

###############################################################

# Simulating the Stage X genotype dynamics

###############################################################

for (i in 1:tlimit){
        nout[,i]  <-  n
    # Introduce new allele by introducing one heterozygote juvenile
    # produced by outcrossing
#if(i == 12) {
#	browser()
#}  
    if (i==10){
#      n[3]  <-  1
      n[9]  <-  1
    }

        #Creating the male gamete pool
        nX    <- n[1:6] + n[7:12] # all individuals can outcross
        ngam  <- ones(c(1,2)) %*% W_prime %*% blkFX_prime %*% nX
        
        #Equation 5 in the manuscript
        q_prime <- (W_prime%*%blkFX_prime%*%nX)/ngam[1] 

        UtildeS  <-  blkUS
        UtildeX  <-  blkUX

        HS <- rbind(c(1, 1/4, 0),
					c(0, 1/2, 0),
					c(0, 1/4, 1))

        HX <- zeros(c(g,g))
        for (j in 1:g){
            pi  <-  Ig[,j]
            qi  <-  W %*% pi #allele frequencies in oocytes of genotype i

            piprime  <-  Z %*% kronecker(qi,q_prime)
            # genotype frequencies in the offspring of mothers of genotype i
            # produced by outcrossing (equation 16 in de Vries and Caswell, 2018a (American Naturalist))

            HX[,j]  <-  piprime # the outcrossing parent-offspring matrix
        }

        blkHS  <-  kronecker(Iom,HS)
        blkHX  <-  kronecker(Iom,HX)

        FtildeS  <-  t(K) %*% blkHS %*% K %*% blkFS
        FtildeX  <-  t(K) %*% blkHX %*% K %*% blkFX
        Atilde   <-  rbind(cbind((UtildeS + FtildeS), FtildeS),
						   cbind(FtildeX            , (UtildeX + FtildeX)))

        nnext  <-  Atilde %*% n
#        nSXnext  <-  Atilde %*% c(n,n)

		if (sum(nnext) > 1e+200) {
			nnext  <-  nnext * 1e-10
			n <- nnext
		} else{
	        n <- nnext
       }

cat('\r', paste(i/tlimit,'%'))
}

Nout  <-  nout[1:6,] + nout[7:12,]
Nout[,1:30]
colSums(nout[,1:30])
# nOut   <-  nout[7:12,]
# nSelf  <-  nout[1:6,]
# nTot   <-  nOut + nSelf      
###############################################################
# Reminder of how output is organized
# n  <-  c[juv_AA, ad_AA, juv_Aa, ad_Aa, juv_aa, ad_aa]

par(mfrow=c(2,2))
plot( 1:tlimit, Nout[1,], log="y", ylim=c(10^-5,max(Nout)), type="l", col="red", lty=2)
lines(1:tlimit, Nout[2, ], col="red", type="l")
lines(1:tlimit, Nout[3, ], col="orange", type="l", lty=2)
lines(1:tlimit, Nout[4, ], col="orange", type="l")
lines(1:tlimit, Nout[5, ], col="green", type="l", lty=2)
lines(1:tlimit, Nout[6, ], col="green", type="l")

plot( 1:tlimit, colSums(Nout[c(2,4,6),]), log="y", type="l", col="red", lty=2)


#Genotype frequencies in females
#GENOTYPE FREQUENCIES, aggregate over all stages
temp         <-  kronecker(diag(g),ones(c(1,om))) %*% Nout
p_genotypes  <-  sweep(temp,2,colSums(temp),'/')
p_genotypes[,1:30]
plot(1:tlimit,p_genotypes[1,],col="red",type="l",ylim = c(0,1))
lines(1:tlimit,p_genotypes[2,],col="orange",type="l")
lines(1:tlimit,p_genotypes[3,],col="green",type="l")

