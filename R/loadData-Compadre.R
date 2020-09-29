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

# load class info for Eagle river population from Petersen et al. (2016).
# this is a perennial montane population that was used as the 'local'
# population in this study. All other populations were grown in a common
# garden in the field at Eagle River.
classInfo        <-  matrixClass(cdb_Mg)
Mg_ER            <-  cdb_Mg[22]
Mg_ER_classInfo  <-  classInfo[[22]]

# Look at matrices
Mg_ER$mat
matList_Mg_ER  <-  list(
						"A"  =  matA(Mg_ER)[[1]],
						"U"  =  matU(Mg_ER)[[1]],
						"F"  =  matF(Mg_ER)[[1]],
						"C"  =  matC(Mg_ER)[[1]]
						)
for(i in 1:4) {
	colnames(matList_Mg_ER[[i]])  <-  classInfo[[22]][,2]
}
matList_Mg_ER

# replace NA's with 0's in U matrix
# this ensures that all individuals 
# die after the 2nd season. No immortals allowed!
matList_Mg_ER$U[3,2:3]  <-  0

# Take a look at lambda for this population
popdemo::eigs(matList_Mg_ER$A, what="lambda")










