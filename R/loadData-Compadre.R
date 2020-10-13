################################################################
#' Load and access demographic data from Compadre database
#' for use with our evolutionary demographic models
#'	
#'  Author: Colin Olito
#' 
#' 	COMPADRE Plant Matrix Database (2020). 
#' 	Available from: https://www.compadre-db.org
#'  [29 09 2020, Version 6.20.9.0 (release date Sep_10_2020)]

###############
# DEPENDENCIES
###############

#' install package 'remotes' if necessary
#' will already be installed if 'devtools' is installed
# install.packages("devtools")
# install.packages("remotes") 

#' install Rcompadre from current github build 
#' we set build_vignettes = FALSE to avoid a possible error. 
#' See: https://jonesor.github.io/Rcompadre/index.html
# remotes::install_github("jonesor/Rcompadre", build_vignettes = FALSE)

# load the package
library(Rcompadre)

# fetch the COMPADRE database
compadre <- cdb_fetch("compadre")


###########################
# Mimulus guttatus example
###########################

# Load data for Mimulus guttatus
cdb_Mg   <-  cdb_check_species(compadre, "Mimulus guttatus", return_db = TRUE)

# Author info and population names
cdb_Mg$Authors
cdb_Mg$MatrixPopulation
cdb_Mg$MatrixStartYear

# load class info for Eagle Meadows population from Petersen et al. (2016).
# this is a perennial montane population that was used as the 'local'
# population in this study. All other populations were grown in a common
# garden in the field at Eagle Meadows
classInfo        <-  matrixClass(cdb_Mg)
Mg_EM            <-  cdb_Mg[22]
Mg_EM_classInfo  <-  classInfo[[22]]

# Look at matrices
Mg_EM$mat
matList_Mg_EM  <-  list(
						"A"  =  matA(Mg_EM)[[1]],
						"U"  =  matU(Mg_EM)[[1]],
						"F"  =  matF(Mg_EM)[[1]],
						"C"  =  matC(Mg_EM)[[1]]
						)
for(i in 1:4) {
	rownames(matList_Mg_EM[[i]])  <-  classInfo[[22]][,2]
	colnames(matList_Mg_EM[[i]])  <-  classInfo[[22]][,2]
}
matList_Mg_EM

# Mistake in Data(?!?). Data from COMPADRE has NA's in U, 
# but according to Matrix 1 in Peterson et al. (2016) 
# and values reported in Table 1 from the Corregendum,
# these values should be as follows:
matList_Mg_EM$A[3,c(2,3)]  <-  0.179*8.71
matList_Mg_EM$U[3,c(2,3)]  <-  0.179*8.71

# Take a look at lambda for this population
popdemo::eigs(matList_Mg_EM$A, what="lambda")












