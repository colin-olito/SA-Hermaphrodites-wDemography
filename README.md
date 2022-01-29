# The demographic costs of sexually antagonistic selection in partially selfing populations

## Overview

This is a GitHub repository for the development of a theoretical evolutionary genetics research project that is now published under the title "*Demographic consequences of sexually antagonistic selection in partially selfing populations*" (doi: [XXX](https://doi.org/...)). Here you can find all of the necessary code to reproduce the simulations presented in the published paper and appendices, and the LaTeX files used to compile the manuscript. Supplementary material for the paper is also available from the publisher [here](https://www.journals.uchicago.edu/toc/an/current), and on the Dryad digital repository [here](https://datadryad.org/stash/share/81sAuXGEg8cSh-S9VVL0PfBCsl6YLkG1OIFBCvOefac).


## Abstract
When selection differs between the sexes, genes expressed by both males and females can experience sexually antagonistic (SA) selection, where beneficial alleles for one sex are deleterious for the other. Classic population genetics theory has been fundamental to understanding how and when SA genetic variation can be maintained by balancing selection, but these models have rarely considered the demographic consequences of coexisting alleles with deleterious fitness effects in each sex. In this paper we develop a stage-structured mendelian matrix model and jointly analyze the evolutionary and demographic consequences of SA selection in obligately outcrossing (i.e., dioecious/gonochorous) and partially selfing hermaphrodite populations. We focus on identifying when SA polymorphisms are maintained by balancing selection *and* the population growth rate remains positive. Additionally, we analyze the effects of inbreeding depression manifesting at different life-history stages and give an illustrative example of the potential for SA polymorphism in real populations using empirically estimated demographic rates for the hermaphroditic flowering plant *Mimulus guttatus*. Our results show that when population intrinsic growth rates approach one, extinction occurs across large swathes of parameter space favoring SA polymorphism or the fixation of male-beneficial alleles, and that inbreeding depression is a significant problem for maintaining SA polymorphism in partially selfing populations. Despite these demographic challenges, our example with *M. guttatus* appears to show that demographic rates observed in some real populations can sustain large regions of viable SA polymorphic space.


## Citing information
*Please cite the paper as*:

 Olito, C. and C. DeVries. 2022. Demographic consequences of sexually antagonistic selection in partially selfing populations. *American Naturalist* XX: XX--XX. doi: XXX

Full citing information will be provided when it is made [available through the publisher](https://www.journals.uchicago.edu/toc/an/current). You can also contact me directly if you would like a reprint. 


## Structure of this repository & instructions for reproducing the results

The key directories in this repository that are needed to reproduce the results and manuscript are as follows:  

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
- **`doc`**  
	- `Refs2.bib`  
	- `SA-Hermaphrodites-wDemography.tex`  
	- **`Supplements`**  
		- `SupplementaryMaterial.tex`  
- **`output`**  
	- **`figs`**  
	- **`simData`**  
- `run-Simulations.R`  
- `makeFigs.R`    

**Note:** Output directories *must be created locally by the user* before running the simulations so that the results can be saved correctly.

In accordance with The American Naturalist's [guidelines for archiving Code with Data](http://comments.amnat.org/2021/12/guidelines-for-archiving-code-with-data.html), a clean version of this repository has been uploaded to Dryad [here](https://datadryad.org/stash/share/81sAuXGEg8cSh-S9VVL0PfBCsl6YLkG1OIFBCvOefac).


###  How to reproduce the results & manuscript

This repository provides all code necessary to (1) rerun the simulations, (2) produce figures as .pdf's, and (3) compile the LaTeX to produce the accepted manuscript, complete with embedded figures. To do so, please follow these basic steps:

1. Clone the repo, and create the output directories locally so that the simulation and figure files can be correctly saved and recalled later.  
2. Check that the following R package dependencies are correctly installed using `install.packages()`:  
	- `extrafont`  
	- `fontcm`
		- **A note on fonts**: We embeded CM Modern font in the .pdf figures produced for the paper. If you prefer to use a different font, instructions are provided on L.34 of `./R/functions-Figs.R`. If you prefer not to use the above font packages, comment out L.4-6 in `./R/functions-Figs.R`, choose a standard font as described on L.34, and comment out the `embed_fonts()` statements `in makeFigs.R`.
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
3. Run the file `run-Simulations.R` either interactively in R or in terminal using Rscript. The simulations will take some time to generate the output files. *We recommend doing this interactively and only running up to L.717*, which will avoid running many simulations contained in a coda to the main simulations.
4. Run the file `makeFigs.R` (up to L.108), which will read the simulation output files and generate the main figures in the paper and supplementary material.  
5. Compile the LaTeX file `SA-Hermaphrodites-wDemography.tex` to produce a .pdf version of the accepted manuscript.  
6. Compile the LaTeX file `SupplementaryMaterial.tex` to produce a .pdf version of the Supplementary Material.  

## An important note on data used in the study

All empirical data used to parameterize the simulations for the *Mimulus gutattus* case study were taken directly from tables presented in previously published studies (i.e., no new data sets were generated or analyzed for this study). Details of how these published data were used in our study is presented in Appendix D of the main paper. The sources for the data can be found here:

- Demographic rates were taken directly from Table 1 in *Peterson et al. (2017) New Phytologist 216: 956–957*. **Note** that this is a corrected version of the table presented in the original publication *Peterson et al. (2016) New Phytologist (2016) 211: 345–356. doi: 10.1111/nph.13971*. The specific parameters from Peterson et al. (2016) used in our simulations is provided in `./data/Peterson_2016_Data.csv`.

- Selfing and inbreeding depression estimates were calculated directly from Tables 1 & 2 in *Willis (1993) Heredity 71:145—154*. The specific parameters from Willis (1993) used in our simulations is provided in `./data/Willis_1993_Data.csv`.

- inv6 selection coefficients were back-calculated from data presented in Figs. 2, 5, & 6 of *Lee et al. (2016) Genetics: 202, 1473–1484*. Our raw measurements were done in Powerpoint, and are provided as .pptx, .odp, and .pdf in the Dryad digital repository for this article [here](https://datadryad.org/stash/share/81sAuXGEg8cSh-S9VVL0PfBCsl6YLkG1OIFBCvOefac). The resulting measurements are provided here in `inv6_Figure_Measurements.csv` and back-calculation to selection coefficients is performed in `inv6-selection-coefficients.R`, which is called by `.R/functions-Figs.R` during plotting.

## Contact & bug reporting

Please report any bugs, problems, or issues by opening an issue on the inversionSize github [issues page](https://github.com/colin-olito/SA-Hermaphrodites-wDemography/issues), or inform me directly by sending a brief email detailing the problem you encountered to colin.olito at biol dot lu dot se.