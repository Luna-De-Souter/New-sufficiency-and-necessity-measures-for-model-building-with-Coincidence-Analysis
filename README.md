# New-sufficiency-and-necessity-measures-for-model-building-with-Coincidence-Analysis

This repository contains supplementary materials for the research paper "New sufficiency and necessity measures for model building with Coincidence Analysis". 

The script "Instructions.R" provides easy-to-follow instructions for CNA users, including those with limited R programming experience, to install and use the development version of the cna package that allows to run the cna algorithm with the new model-building measures presented in the paper. 

The script "replication.R" contains code to replicate the examples and simulation experiments presented in the paper. Running this script requires the files "AuxFuncs.R", "cna_3.5.3.4.tar.gz", and "score_defs.R". As replicating the entire data simulation and data analysis of the experiment, which happens on lines 62-249 of the "replication.R" script, is very computationally demanding and time-consuming, readers whose main interest is to study the outcomes of the simulation experiments in detail are advised to skip these lines in the script. After downloading the file "scoreboards.per.defscore.RData", they can use the command 
"scoreboards.per.defscore <- readRDS("scoreboards.per.defscore.RData")" to access the results of the experiment and run the lines 250-319 in the "replication.R" script to regenerate the plots presented in the paper. The files "data_cs_1000.RData" and "frscores.per.definition.RData" contain the simulated datasets used in the experiment and the results of CNA analyses performed on those datasets in the experiment, respectively.
