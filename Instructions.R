# This script contains instructions to install a preliminary version of the cna
# package that can build cna models using weighted consistency, instead of 
# consistency, and weighted c-coverage, rather than coverage, as model-building
# measures. Furthermore, instructions on how to run a cna analysis with the 
# new model building measures are included.
# I thank Mathias Ambühl for implementing this preliminary version.
# Please email me at luna.souter@uib.no if you have questions.




### Instructions to install the new implementation ###

## Download the folder "cna_3.5.3.4.tar.gz" and make sure it is in the same 
## folder as this script.

# First, uninstall the current version of the cna package (If this gives 
# problems, restart your R session and try again. Hopefully, that solves
# any issues.)
remove.packages("cna")
# Set your R working directory to the folder containing this script and the 
# folder "cna_3.5.3.4.tar.gz".
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Now, install the version of the cna package with the new model-building
# measures.
install.packages("cna_3.5.3.4.tar.gz", repos = NULL, type ="source")
# Load cna.
library(cna)
# Test if the correct version is installed and loaded. (If no messages appear,
# everything is fine.)
stopifnot(packageVersion("cna") == "3.5.3.4")



### Instructions on how to use the new implementation ###


## Conduct cna analysis with weighted consistency and weighted c-coverage:

# First, define the new model-building measures of weighted consistency and weighted
# c-coverage, as follows. 
# (This code is unintelligible without knowing the codes
# we used in this preliminary implementation, and understanding them is not 
# needed for trying out these new measures.)
def_w_measures <- ccDef_ratio(list(
  cbind(c(1, 0, 0, 0, 0, 0),
        c(1, 0, -99-1i, 0, -99-1i, 0)),
  cbind(c(0, 1, 0, 0, 0, 0),
        c(0, 1, 0, -99, -99, 0))))
# Run the cna analysis with weighted consistency and weighted c-coverage.
cna(d.educate, con = 0.8, cov = 0.8, ccDef = def_w_measures)


# Conduct the cna analysis with regular consistency and coverage (these measures 
# are still the default in this preliminary implementation):
cna(d.educate, con = 0.8, cov = 0.8)


# Calculate weighted consistency and weighted c-coverage for model 
# (L <-> E)*(U + D <-> L) for dataset d.educate.
condTbl("(L <-> E)*(U + D <-> L)", d.educate, ccDef = def_w_measures)

# Calculate regular consistency and coverage for model (L <-> E)*(U + D <-> L) 
# for dataset d.educate.
condTbl("(L <-> E)*(U + D <-> L)", d.educate)

# After playing around with the new measures, detach and uninstall this 
# preliminary version of the package, and reinstall the official CRAN 
# version of the CNA package to go back to normal.
detach("package:cna", unload = TRUE)
remove.packages("cna")
install.packages("cna")



