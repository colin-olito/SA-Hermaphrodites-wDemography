################################################################
#  Matrix Model Functions
#  Useful functions for implementing matrix models for stage-
#  structured populations
#
#  Author: Colin Olito, adapted from L. DeVries
#
#  NOTES:  
#		


#rm(list=ls())
#####################
##  Dependencies

#####################
##  Functions

# create an array of 0's with dimensions dims...
zeros  <-  function(dims) {
	array(0,dim=dims)
}

# create an array of 1's of with dimensions dims... 
ones  <-  function(dims) {
	array(1,dim=dims)
}

# emat(m,n,i,j) creates an m by n matrix with 1 in the i,j entry and zeros elsewhere
emat <- function(m,n,i,j){
	E <- matrix(0,nrow=m,ncol=n)
	E[i,j] <- 1
	return(E)
}

# vecperm
# function to calculate the vec permutation matrix of size m,n
# 4/9/03
# function p = vecperm(m,n)
vecperm <- function(m,n){
p <- as.vector(zeros(m*n))
a <- zeros(c(m,n))
for (i in 1:m){
  for (j in 1:n){
    e <- a
    e[i,j] <- 1
    p <- p + kronecker(e,t(e))
    }
}
return(p)
}


rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

squashIm <- function(x) {
	if (all(Im(z <- zapsmall(x))==0)) as.numeric(z) else x
}

eucDist <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2))
} 