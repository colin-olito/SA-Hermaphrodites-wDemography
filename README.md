# The demographic costs of sexually antagonistic selection in partially selfing populations

## Overview

This is a GitHub repository for a theoretical evolutionary genetics research project that is now published under the title "*Demographic consequences of sexually antagonistic selection in partially selfing populations*" (doi: [XXX](https://doi.org/...)). Here you can find all of the necessary code to reproduce the simulations and figures from the published paper and appendices. Supplementary material for the paper is also available from the publisher [here](https://www.journals.uchicago.edu/toc/an/current), and on the Dryad/Zenodo digital repository [here](https://doi.org/10.5061/dryad.c2fqz619t), temporary reviewer link [here](https://datadryad.org/stash/share/81sAuXGEg8cSh-S9VVL0PfBCsl6YLkG1OIFBCvOefac).


## Abstract
When selection differs between the sexes, genes expressed by both males and females can experience sexually antagonistic (SA) selection, where beneficial alleles for one sex are deleterious for the other. Classic population genetics theory has been fundamental to understanding how and when SA genetic variation can be maintained by balancing selection, but these models have rarely considered the demographic consequences of coexisting alleles with deleterious fitness effects in each sex. In this paper we develop a stage-structured mendelian matrix model and jointly analyze the evolutionary and demographic consequences of SA selection in obligately outcrossing (i.e., dioecious/gonochorous) and partially selfing hermaphrodite populations. We focus on identifying when SA polymorphisms are maintained by balancing selection *and* the population growth rate remains positive. Additionally, we analyze the effects of inbreeding depression manifesting at different life-history stages and give an illustrative example of the potential for SA polymorphism in real populations using empirically estimated demographic rates for the hermaphroditic flowering plant *Mimulus guttatus*. Our results show that when population intrinsic growth rates approach one, extinction occurs across large swathes of parameter space favoring SA polymorphism or the fixation of male-beneficial alleles, and that inbreeding depression is a significant problem for maintaining SA polymorphism in partially selfing populations. Despite these demographic challenges, our example with *M. guttatus* appears to show that demographic rates observed in some real populations can sustain large regions of viable SA polymorphic space.


## Citing information
*Please cite the paper as*:

 Olito, C. and C. DeVries. 2022. Demographic consequences of sexually antagonistic selection in partially selfing populations. *American Naturalist* XX: XX--XX. doi: XXX

Full citing information will be provided when it is made [available through the publisher](https://www.journals.uchicago.edu/toc/an/current). You can also contact me directly if you would like a reprint. 


##  Instructions

This repository provides all code necessary to (1) rerun the simulations and (2) produce figures as .pdf's. To do so, please follow these basic steps:

1. Clone the repo using the following: `git clone https://github.com/colin-olito/SA-Hermaphrodites-wDemography`. Alternatively, on the project main page on GitHub, click on the green button `clone` or `download` and then click on `Download ZIP`.  
2. **Create the output directories** `./output/figs` and `./output/simData` on your local machine so that the simulation and figure files can be correctly saved.  
3. Check that you have a recent version of [`R`](https://www.r-project.org/) installed. 
4. Check that the following R package dependencies are correctly installed using `install.packages()`:  
	- `Rcompadre`
		- **NOTE:** At the time of acceptance one of the sub-dependencies (`popdemo`) only works with `R` v. >= 4.1.0. See [https://github.com/jonesor/Rcompadre](https://github.com/jonesor/Rcompadre) for more installation details.
		- The database can also be downloaded from [https://www.compadre-db.org/Data/CompadreDownload](https://www.compadre-db.org/Data/CompadreDownload). Save it to `./data/` in your local repository, and then run the line `compadre <- cdb_fetch("./data/COMPADRE_v.X.XX.X.X.RData.txt")` in `R` (replace `X`'s with the current version number of the database).
	- `extrafont`  
	- `plyr`  
	- `lattice`  
	- `latticeExtra`  
	- `wesanderson`  
	- `MASS`  
	- `raster`  
	- `akima`  
	- `foreach`  
	- `doParallel`  
	- `doSNOW`  
5. Run `run-Simulations.R` either interactively in R or in terminal, being sure to set the working directory to the root directory of the repo (e.g., `SA-Hermaphrodites-wDemography-master/`). The simulations will take some time to generate the output files. *We recommend doing this interactively and only running up to L.717*, which will avoid running many simulations contained in a coda to the main simulations.
6. Run `makeFigs.R` (up to L.108), which will read the simulation output files and generate the main figures in the paper and supplementary material.  


## Repostiory structure and contents 

The directories/files in this repostiory needed to reproduce the results for this study are as follows:  

- **`R`**   
	- `functions-Figs.R`  
	- `functions-MatModels.R`  
	- `functions-Simulations.R`  
	- `loadData-Compadre.R`  
	- `inv6-selection-coefficients.R`  
- **`data`**  
	- `Peterson_2016_Data.csv`  
	- `Willis_1993_Data.csv`  
	- `inv6_Figure_Measurements.csv`  
- **`output`***  
	- **`figs`**  
	- **`simData`**  
- `run-Simulations.R`  
- `makeFigs.R`  
- `LICENSE.txt`   

**Note:** * Output directories *must be created locally by the user* before running the simulations so that the results can be saved correctly.

### File & variable descriptions

Function & data processing files
- `functions-Figs.R`: collection of plotting functions called by `makeFigs.R`.   
- `functions-MatModels.R`: convenience functions for matrix models called by `functions-Simulations.R`.   
- `functions-Simulations.R`: workhorse simulation functions called by `run-Simulations.R`.  
- `loadData-Compadre.R`: loads the Compadre demographic database and extracts data for *M. guttatus.   
- `inv6-selection-coefficients.R`: Data processing file. Calculates selection coefficients for inv6 from raw figure measuremnts provided in `inv6_Figure_Measurements.csv`.  

Data files (variables in bullets)
- `Peterson_2016_Data.csv`: data from Peterson et al. (2016) used in this study.   
	- population: Population identifier  
	- D: seed bank survival rate  
	- G: seeed germination rate  
	- F: flower production  
	- O: ovules per flower  
	- A: seedling recruits  
	- S: overwinter survival  
	- R: rosette production  
- `Willis_1993_Data.csv`: data from Willis (1993) used in this study. Variables:
	- C: selfing rate  
	- delta_D: inbreeding depression term for seed bank survival rate  
	- delta_G: inbreeding depression term for seeed germination rate  
	- delta_F: inbreeding depression term for flower production  
	- delta_O: inbreeding depression term for ovules per flower  
	- delta_S: inbreeding depression term for overwinter survival  
- `inv6_Figure_Measurements.csv`: data extracted from figures in Lee. et al. (2016) used in this study.  
	- year: self explanatory  
	- no_inv6_flowers_inches: wild-type flower number, barplot length in inches.  
	- inv6_het_flowers_inches: inv6 heterozygote flower number, barplot length in inches  
	- no_inv6_fruits_inches: wild-type fruit number, barplot length in inches  
	- inv6_het_fruits_inches: inv6 heterozygote fruit number, barplot length in inches  
	- no_inv6_seeds_inches: wild-type seed number, barplot length in inches  
	- inv6_het_seeds_inches: inv6 heterozygote seed number, barplot length in inches  
	- no_inv6_PrViablePollen_inches: wild-type proportion viable pollen, barplot length in inches  
	- inv6_het_PrViablePollen_inches: inv6 heterozygote proportion viable pollen, barplot length in inches  
	- scaleUnits: Units of original plot y-axes.  
	- scaleAxis: Scale in original units for y-axes  
	- pptInches: Scale in inches  
- Files with the raw raw measurements of the figures in Lee et al. (2016) are available on Dryad/Zenodo under the name `inv6_measurements.pptx` (also as .odp, .pdf), or upon request from the authors.

Executables
- `run-Simulations.R`: Exectuable functions to run simulations; calls functions/objects defined in `functions-MatModels.R`, `functions-Simulations.R`, `loadData-Compadre.R`, and `inv6-selection-coefficients.R`).    
- `makeFigs.R`: executable plotting functions to generate .pdf figures using simulation results. Calls functions/object defined in `functions-Figs.R` and `inv6-selection-coefficients.R`.

License    
- `LICENSE.txt`: MIT license for this repository.  


### DRYAD
A copy of this repository was uploaded to Dryad at the time of manuscript acceptance [here](https://doi.org/10.5061/dryad.c2fqz619t), temporary Dryad [reviewer link](https://datadryad.org/stash/share/81sAuXGEg8cSh-S9VVL0PfBCsl6YLkG1OIFBCvOefac).



## An important note on data used in the study

All empirical data used to parameterize the simulations for the *Mimulus gutattus* case study were taken directly from tables and figures presented in previously published studies (i.e., no *new* data was generated for this study). Details of how these published data were used to parameterize our demographic model are presented in Appendix D to the main paper. The sources for the data can be found here:

- Demographic rates were taken directly from Table 1 in *Peterson et al. (2017) New Phytologist 216: 956–957*. **Note** that this is a corrected version of the table presented in the original publication *Peterson et al. (2016) New Phytologist (2016) 211: 345–356. doi: 10.1111/nph.13971*. The specific parameter values used in our simulations are provided in `./data/Peterson_2016_Data.csv`.

- Selfing and inbreeding depression estimates were calculated directly from Tables 1 & 2 in *Willis (1993) Heredity 71:145—154*. The specific parameter values used in our simulations are provided in `./data/Willis_1993_Data.csv`.

- inv6 selection coefficients were extracted from Figs. 2, 5, & 6 of *Lee et al. (2016) Genetics: 202, 1473–1484*. Raw measurements were performed in Powerpoint, and are provided as .pptx, .odp, and .pdf in the Dryad digital repository for this article [here](https://doi.org/10.5061/dryad.c2fqz619t), temporary Dryad [reviewer link](https://datadryad.org/stash/share/81sAuXGEg8cSh-S9VVL0PfBCsl6YLkG1OIFBCvOefac). The resulting measurements are provided here in `inv6_Figure_Measurements.csv` and back-calculation to selection coefficients is performed in `inv6-selection-coefficients.R`, which is called by `.R/functions-Figs.R` during plotting.

## Contact & bug reporting

Please report any bugs, problems, or issues by opening an issue on the inversionSize github [issues page](https://github.com/colin-olito/SA-Hermaphrodites-wDemography/issues), or inform me directly by sending a brief email detailing the problem you encountered to colin dot olito at biol dot lu dot se.

## Licence information

This repository is provided by the authors under the MIT License ([MIT](https://opensource.org/licenses/MIT)).

