## Code to check the accuracy of conclusions drawn using the MCMC samples produced when fitting SCR models that assume the state process can be represented by an inhomogeneous Poisson process

## Libraries we need
library(secr)
## Objects we need
load("../output/mona-inputs.RData")
## Functions we need
source("Functions.R")

## ---------------------------------------------------------------------------------------
# Checking the results for Figure 4
## ---------------------------------------------------------------------------------------

##### Loading in the MCMC samples #####

load("MCMC_Results/Figure4/InhomPP_18occ.RData")
load("MCMC_Results/Figure4/InhomPP_52occ.RData")
load("MCMC_Results/Figure4/InhomPP_111occ.RData")

## Discarding buring in (are discarding 1000 iterations as burn-in, see bayesian_code/Plots_Code.R for more info)
inhom.results.18occ <- inhom.results.18occ[-c(1:1000),]
inhom.results.52occ <- inhom.results.52occ[-c(1:1000),]
inhom.results.111occ <- inhom.results.111occ[-c(1:1000),]

##### Comparing results based on MCMC samples to results from secr.fit() #####

## Looking at the first plot in Row 2
check.inhom.mcmc(results=inhom.results.18occ, j=1, mask=mona_mask, array="3x3")

## Looking at the second plot in Row 2
check.inhom.mcmc(results=inhom.results.52occ, j=2, mask=mona_mask, array="3x3")

## Looking at the third plot in Row 2
check.inhom.mcmc(results=inhom.results.111occ, j=3, mask=mona_mask, array="3x3")

# There are some differences in both sets of results. However, things generally look good. The true value of lambda and sigma lies within the credible intervals in all three instances. The results from secr.fit() tend to align relatively closely to the results produced using the MCMC samples. So, it looks like our MCMC sampling here worked as desired. 

## ---------------------------------------------------------------------------------------
# Checking the results for Figure 5
## ---------------------------------------------------------------------------------------

##### Loading in the MCMC samples #####

load("MCMC_Results/Figure5/InhomPP_7occ.RData")
load("MCMC_Results/Figure5/InhomPP_25occ.RData")
load("MCMC_Results/Figure5/InhomPP_55occ.RData")

## Discarding buring in (are discarding 1000 iterations as burn-in, see bayesian_code/Plots_Code.R for more info)
inhom.results.7occ <- inhom.results.7occ[-c(1:1000),]
inhom.results.25occ <- inhom.results.25occ[-c(1:1000),]
inhom.results.55occ <- inhom.results.55occ[-c(1:1000),]

##### Comparing MCMC results to results from secr.fit() #####

## Looking at the first plot in Row 2
check.inhom.mcmc(results=inhom.results.7occ, j=1, mask=mona_mask, array="7x7")

## Looking at the second plot in Row 2
check.inhom.mcmc(results=inhom.results.25occ, j=2, mask=mona_mask, array="7x7")

## Looking at the third plot in Row 2
check.inhom.mcmc(results=inhom.results.55occ, j=3, mask=mona_mask, array="7x7")

# Overall, both sets of results align nicely, indicating our MCMC sampling worked as desired. However, true value of lambda0 and sigma does not fall in credible intervals for 7 occasions. But perhaps this is expected, due to a lack of information. 
