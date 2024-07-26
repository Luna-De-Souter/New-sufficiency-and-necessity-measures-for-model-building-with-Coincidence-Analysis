# New-sufficiency-and-necessity-measures-for-model-building-with-Coincidence-Analysis
[Unfinished version of README.]

This repository contains supplementary materials for the research paper "New sufficiency and necessity measures for model building with Coincidence Analysis". It includes materials to replicate all examples and simulation experiments presented in this paper. Furthermore, the script Instructions.R provides easy-to-follow instructions for CNA users, including those with limited R programming experience, to install and use the development version of the cna package that allows to run the cna algorithm with the new model-building measures presented in the paper.

The script "replication.R" ...

can be used to replicate the simulation experiments presented in the research paper. It includes code to install the development version of the cna R package that allows to run the algorithm with the new model-building measures from the file cna_3.5.3.4.tar.gz in this repository. It also sources the script Auxfuncs.R in this repository. 
It replicates the tables and CNA analyses used as illustrations in the paper. 
It then replicates the entire simulation experiment, including the data simulation and the running of CNA analyses.
The data generation (lines 62 to 109 in the Replication.R script) and especially running the CNA analyses (lines 110-163) are very computationally demanding. Therefore, readers whose main interest is having a closer look at the results of the experiment and/or recreating the plots presented in the paper are advised to skip this part of the script and instead study the results based on the file frscores.per.definition.RData and/or file scoreboards.per.defscore.RData, both of which are included in this repository, as done on lines 164-249 in the Replication.R script, and/or to recreate the plots by using file scoreboards.per.defscore.RData, which is included in this repository, as done on lines 250-318 in the Replication.R script.


[Check if the scoring after line 163 is computationally doable.]
