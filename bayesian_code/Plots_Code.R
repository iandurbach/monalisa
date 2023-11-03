## Code used to produce Bayesian versions of Figures 4 and 5
# Working directory should be the 'bayesian_code' folder

## Libraries we need
library(nimble)
library(coda)
library(nimbleSCR)
library(spatstat)
library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)
library(secr)
library(patchwork)
library(ggpubr)
## Objects we need
load("../output/mona-inputs.RData")
load("../output/mona-results.RData")
## Functions we need
source("Functions.R")

## ---------------------------------------------------------------------------------------
######################## Creating the objects needed for Figure 4 ########################
## ---------------------------------------------------------------------------------------

## In this figure, the first column uses simulated data w/ 18 sampling occasions. The second column uses simulated data w/ 52 sampling occasions, and the third column uses simulated data w/ 111 sampling occasions.

## ---------------------------------------------------------------------------------------
# Creating the data objects we need for all MCMC samples for Figure 4
## ---------------------------------------------------------------------------------------

## Function we will use
# The only argument we provide is 'j': this will be a value in {1, 2, 3} and represents the index of the objects we want to work with from the RData objects we have loaded in. If want objects generated using 18, 52 or 111 sampling occasions, j=1,2,3 respectively.
organise.data = function(j) {
  # Number of sampling occasions used for simulated data
  nocc  <- capthists_few_alloccs_3x3$noccasions[j]

  # Summing capture histories over all of the simulated sampling occasions
  all.dat <- capthists_few_alloccs_3x3$capthist[[j]]
  all.mat <- matrix(0, nrow=nrow(all.dat[,1,]), ncol=ncol(all.dat[,1,]))
  for (i in 1:nocc) {
    all.mat <- all.mat + all.dat[,i,]
  }

  # Trap locations
  trap.loc <- attributes(all.dat)$traps

4  # xlim, ylim (for our map area)
  xlim <- c(0.5, 50.5)
  ylim <- c(0.5, 50.5)

  # Creating the data object
  data <- list(encounter.data = all.mat, trap.loc = trap.loc, xlim = xlim, ylim = ylim, n.occasions = nocc)
  data
}

# Data object for 18 sampling occasions
data.18occ <- organise.data(1)

# Data object for 52 sampling occasions
data.52occ <- organise.data(2)

# Data object for 111 sampling occasions
data.111occ <- organise.data(3)

## ---------------------------------------------------------------------------------------
# Creating the objects we specifically need for the plots in Row 1 of Figure 4
## ---------------------------------------------------------------------------------------

## Row 1 consists of RACD maps. The SCR models that we fit to create these maps assume that the state process (the random process governing the distribution of the activity centres) is a homogeneous Poisson process.

##### Running the MCMC #####

## Uncomment the lines below if want to run the MCMC
## Running MCMC for simulated data from 18, 52 and 111 sampling occasions.
#results.18occ <- run.MCMC(data=data.18occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
#save(results.18occ, file="MCMC_Results/Figure4/HomPP_18occ.RData")

#results.52occ <- run.MCMC(data=data.52occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
#save(results.52occ, file="MCMC_Results/Figure4/HomPP_52occ.RData")

#results.111occ <- run.MCMC(data=data.111occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
#save(results.111occ, file="MCMC_Results/Figure4/HomPP_111occ.RData")

## Loading in the RData files containing the MCMC results
load("MCMC_Results/Figure4/HomPP_18occ.RData")
load("MCMC_Results/Figure4/HomPP_52occ.RData")
load("MCMC_Results/Figure4/HomPP_111occ.RData")

## Note that the burn-in iterations for these samples are discarded automatically, so we don't need to worry about this.
## Checking trace plots to make sure everything looks okay
# 18 sampling occasions -- looks good
check.trace.plots(results.18occ)
# 52 sampling occasions -- looks good
check.trace.plots(results.52occ)
# 111 sampling occasions -- looks good
check.trace.plots(results.111occ)

##### Creating the objects we need for the RACD plots #####

## Row 1 consists of RACD maps. So, we will create vectors that contain the posterior mean of the number of activity centres in each pixel -- these are the density values for each pixel in RACD maps that are based on MCMC results.
racd.18occ <- racd.density.vector(results=results.18occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
racd.52occ <-  racd.density.vector(results=results.52occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
racd.111occ <- racd.density.vector(results=results.111occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))

## ---------------------------------------------------------------------------------------
# Creating the objects we specifically need for the plots in Row 2 of Figure 4
## ---------------------------------------------------------------------------------------

## Row 2 consists of EACD maps. The SCR models that we fit to create these maps assume that the state process (the random process governing the distribution of the activity centres) is an inhomogeneous Poisson process.

##### 'pixel.info' object needed for MCMC #####

## Uncomment if want to run MCMC
#pixel.centres <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50)
#pixel.info <- cbind(pixel.centres, log.dblur)

##### Covariate value for each pixel #####

# Note, see 'Functions.R' for an explanation of this function. It basically just working with the data objects we have loaded into R by this point to extract the covariate values for each pixel.
log.dblur <- eacd.covariate()

##### Running the MCMC #####

## Uncomment the lines below if want to run the MCMC. Note that depending on the computer, all of the MCMC chains below can take below 5 or 10 minutes to run.
# 18 sampling occasions, saving the results
#inhom.results.18occ <- run.MCMC.inhom(data=data.18occ, pixel.info=pixel.info, M=300, inits.vec=c(10, 4, 0), n.iter=100000, n.burn=1000)
#save(inhom.results.18occ, file="MCMC_Results/Figure4/InhomPP_18occ.RData")

# 52 sampling occasions, saving the results
#inhom.results.52occ <- run.MCMC.inhom(data=data.52occ, pixel.info=pixel.info, M=300, inits.vec=c(10, 4, 0), n.iter=100000, n.burn=1000)
#save(inhom.results.52occ, file="MCMC_Results/Figure4/InhomPP_52occ.RData")

# 111 sampling occasions, saving the results
#inhom.results.111occ <- run.MCMC.inhom(data=data.111occ, pixel.info=pixel.info, M=300, inits.vec=c(10, 4, 0), n.iter=100000, n.burn=1000)
#save(inhom.results.111occ, file="MCMC_Results/Figure4/InhomPP_111occ.RData")

## Loading in the RData files containing the MCMC results
load("MCMC_Results/Figure4/InhomPP_18occ.RData")
load("MCMC_Results/Figure4/InhomPP_52occ.RData")
load("MCMC_Results/Figure4/InhomPP_111occ.RData")

## Note that with these MCMC samples, the burn-in isn't discarded automatically. So, each object contains data from 101,000 MCMC iterations
## Checking trace plots to decide how many iterations to discard as burn-in -- if we don't discard any iterations, we clearly see some burn-in on the trace plots. If we discard 1000 iterations as burn-in, the trace plots look good (see below)
# 18 sampling occasions
check.trace.plots(inhom.results.18occ[-c(1:1000),], inhom=T)
# 52 sampling occasions
check.trace.plots(inhom.results.52occ[-c(1:1000),], inhom=T)
# 111 sampling occasions
check.trace.plots(inhom.results.111occ[-c(1:1000),], inhom=T)

## So, discarding 1000 iterations as burn-in for our 3 MCMC samples
inhom.results.18occ <- inhom.results.18occ[-c(1:1000),]
inhom.results.52occ <- inhom.results.52occ[-c(1:1000),]
inhom.results.111occ <- inhom.results.111occ[-c(1:1000),]

##### Creating the objects we need for the EACD plots #####

## Creating vectors containing density values for each pixel when working with 18/52/111 sampling occasions
eacd.18occ <- eacd.density.vector(results=inhom.results.18occ, covariate=log.dblur, nPix=2500)
eacd.52occ <- eacd.density.vector(results=inhom.results.52occ, covariate=log.dblur, nPix=2500)
eacd.111occ <- eacd.density.vector(results=inhom.results.111occ, covariate=log.dblur, nPix=2500)

## ---------------------------------------------------------------------------------------
######################## Creating the objects needed for Figure 5 ########################
## ---------------------------------------------------------------------------------------

## In this figure, the first column uses simulated data w/ 7 sampling occasions. The second column uses simulated data w/ 25 sampling occasions, and the third column uses simulated data w/ 55 sampling occasions.

## ---------------------------------------------------------------------------------------
# Creating the data objects we need for all MCMC samples for Figure 5
## ---------------------------------------------------------------------------------------

## Function we will use
# The only argument we provide is 'j': this will be a value in {1, 2, 3} and represents the index of the objects we want to work with from the RData objects we have loaded in. If want objects generated using 7, 25 or 55 sampling occasions, j=1,2,3 respectively.
organise.data = function(j) {
  # Number of sampling occasions used for simulated data
  nocc  <- capthists_few_alloccs_7x7$noccasions[j]

  # Summing capture histories over all of the simulated sampling occasions
  all.dat <- capthists_few_alloccs_7x7$capthist[[j]]
  all.mat <- matrix(0, nrow=nrow(all.dat[,1,]), ncol=ncol(all.dat[,1,]))
  for (i in 1:nocc) {
    all.mat <- all.mat + all.dat[,i,]
  }

  # Trap locations
  trap.loc <- attributes(all.dat)$traps

  # xlim, ylim (for our map area)
  xlim <- c(0.5, 50.5)
  ylim <- c(0.5, 50.5)

  # Creating the data object
  data <- list(encounter.data = all.mat, trap.loc = trap.loc, xlim = xlim, ylim = ylim, n.occasions = nocc)
  data
}

# Data object for 7 sampling occasions
data.7occ <- organise.data(1)

# Data object for 25 sampling occasions
data.25occ <- organise.data(2)

# Data object for 55 sampling occasions
data.55occ <- organise.data(3)

## ---------------------------------------------------------------------------------------
# Creating the objects we specifically need for the plots in Row 1 of Figure 5
## ---------------------------------------------------------------------------------------

## Row 1 consists of RACD maps. The SCR models that we fit to create these maps assume that the state process (the random process governing the distribution of the activity centres) is a homogeneous Poisson process.

##### Running the MCMC #####

## Uncomment the lines below if want to run the MCMC
## Running MCMC for simulated data from 7, 25 and 55 sampling occasions.
#results.7occ <- run.MCMC(data=data.7occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
#save(results.7occ, file="MCMC_Results/Figure5/HomPP_7occ.RData")

#results.25occ <- run.MCMC(data=data.25occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
#save(results.25occ, file="MCMC_Results/Figure5/HomPP_25occ.RData")

#results.55occ <- run.MCMC(data=data.55occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
#save(results.55occ, file="MCMC_Results/Figure5/HomPP_55occ.RData")

## Loading in the RData files containing the MCMC results
load("MCMC_Results/Figure5/HomPP_7occ.RData")
load("MCMC_Results/Figure5/HomPP_25occ.RData")
load("MCMC_Results/Figure5/HomPP_55occ.RData")

## Note that the burn-in iterations for these samples are discarded automatically, so we don't need to worry about this.
## Checking trace plots to make sure everything looks okay
# 7 sampling occasions -- looks good
check.trace.plots(results.7occ)
# 25 sampling occasions
check.trace.plots(results.25occ)
# 55 sampling occasions
check.trace.plots(results.55occ)
# For 25 and 55 sampling occasions, we have very restricted mixing of N (N only ranges over a small range of values). As D is calculated directly from N, this means we also have restricted mixing of D. Note that for 25 and 55 sampling occasions, the lower limit of N is equal to the total number of observed animals, which seems sensible.
# Changing the value of M or the prior for sigma doesn't change this. Lambda0 and sigma seem to be mixing well, which seems to further reinforce that changing the priors for lambda0 and sigma would not help.
# It seems likely that in these cases, having many traps (49 traps) and many sampling occasions (25, 55 sampling occasions) means that we have enough information that we become fairly sure of the true value of N in the region of interest, resulting in the limited mixing we see. This seems to be supported by the fact that if we run: 'table(results.55occ[,"N"]); table(results.25occ[,"N"])', we can see that with 55 sampling occasions, the mixing for N is considerably more limited than for 25 occasions (we are even more sure about the value of N, due to the increase in sampling occasions). So, we will continue.

##### Creating the objects we need for the RACD plots #####

## Row 1 consists of RACD maps. So, we will create vectors that contain the posterior mean of the number of activity centres in each pixel
racd.7occ <- racd.density.vector(results=results.7occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
racd.25occ <-  racd.density.vector(results=results.25occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
racd.55occ <- racd.density.vector(results=results.55occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))

## ---------------------------------------------------------------------------------------
# Creating the objects we specifically need for the plots in Row 2 of Figure 5
## ---------------------------------------------------------------------------------------

## Row 2 consists of EACD maps. The SCR models that we fit to create these maps assume that the state process (the random process governing the distribution of the activity centres) is an inhomogeneous Poisson process.

##### 'pixel.info' object needed for MCMC #####

## Uncomment if want to run MCMC
#pixel.centres <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50)
#pixel.info <- cbind(pixel.centres, log.dblur)

##### Covariate value for each pixel #####

# Already found above when creating objects needed for row 2 of Figure 4 (same covariate is used!)
log.dblur

##### Running the MCMC #####

## Uncomment the lines below if want to run the MCMC. Note that all of the MCMC chains below take below 10 minutes to run.
# 7 sampling occasions, saving the results
#inhom.results.7occ <- run.MCMC.inhom(data=data.7occ, pixel.info=pixel.info, M=300, inits.vec=c(10, 4, 0), n.iter=100000, n.burn=1000)
#save(inhom.results.7occ, file="MCMC_Results/Figure5/InhomPP_7occ.RData")

# 25 sampling occasions, saving the results
#inhom.results.25occ <- run.MCMC.inhom(data=data.25occ, pixel.info=pixel.info, M=300, inits.vec=c(10, 4, 0), n.iter=100000, n.burn=1000)
#save(inhom.results.25occ, file="MCMC_Results/Figure5/InhomPP_25occ.RData")

# 55 sampling occasions, saving the results
#inhom.results.55occ <- run.MCMC.inhom(data=data.55occ, pixel.info=pixel.info, M=300, inits.vec=c(10, 4, 0), n.iter=100000, n.burn=1000)
#save(inhom.results.55occ, file="MCMC_Results/Figure5/InhomPP_55occ.RData")

## Loading in the RData files containing the MCMC results
load("MCMC_Results/Figure5/InhomPP_7occ.RData")
load("MCMC_Results/Figure5/InhomPP_25occ.RData")
load("MCMC_Results/Figure5/InhomPP_55occ.RData")

## Note that with these MCMC samples, the burn-in isn't discarded automatically. So, each object contains data from 101,000 MCMC iterations
## Checking trace plots to decide how many iterations to discard as burn-in -- if we don't discard any iterations, we clearly see some burn-in on the trace plots. If we discard 1000 iterations as burn-in, the trace plots look good (see below)
# 18 sampling occasions
check.trace.plots(inhom.results.7occ[-c(1:1000),], inhom=T)
# 52 sampling occasions
check.trace.plots(inhom.results.25occ[-c(1:1000),], inhom=T)
# 111 sampling occasions
check.trace.plots(inhom.results.55occ[-c(1:1000),], inhom=T)
# Overall, things look good. We once again see restricted mixing in the trace plots for N and D as we did above. Once again, we believe this is due to an increase in information meaning that we have increased confidence in the value of N (and therefore D), so we will continue.

## So, discarding 1000 iterations as burn-in for our 3 MCMC samples
inhom.results.7occ <- inhom.results.7occ[-c(1:1000),]
inhom.results.25occ <- inhom.results.25occ[-c(1:1000),]
inhom.results.55occ <- inhom.results.55occ[-c(1:1000),]

##### Creating the objects we need for the EACD plots #####

## Creating vectors containing density values for each pixel when working with 7/25/55 sampling occasions
eacd.7occ <- eacd.density.vector(results=inhom.results.7occ, covariate=log.dblur, nPix=2500)
eacd.25occ <- eacd.density.vector(results=inhom.results.25occ, covariate=log.dblur, nPix=2500)
eacd.55occ <- eacd.density.vector(results=inhom.results.55occ, covariate=log.dblur, nPix=2500)

## ---------------------------------------------------------------------------------------
########################### Objects needed for Figures 4 and 5  ##########################
## ---------------------------------------------------------------------------------------

##### Creating a data frame that contains all of the info required for both figures #####

## Creating a data frame, labelled 'predicted_densities_all', that summarises all of the information that we will use to create the plots included in Figures 4 and 5

## Function to summarise the data for the RACD maps
# The 'nocc' argument is the number of sampling occasions, and 'fig' is the number of the figure for which we wish to summarise data
racd.summary <- function(nocc, fig) {
  # Name of array we are working with -- if Figure 4, the name is '3x3' and if Figure 5, the name is '7x7'
  if (fig==4) {
    array <- "3x3"
  } else {
    if (fig==5) {
      array <- "7x7"
    }
  }
  # Obtaining the values for the necessary RACD map
  racd.vals <- get(paste0("racd.", nocc, "occ"))

  # Pixel centres that we are working with
  pixel.centres <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50)
  # Subtracting 0.5 so we are dealing with pixel edges
  pixel.edges <- pixel.centres - 0.5
  # We want to work with pixel edges rather than pixel centres, as this means that our resulting map will be coloured correctly (colours will extend to the edge of each pixel, based on density of that pixel). Otherwise, if we work with pixel centres, colours will go from centre to centre (so we won't be colouring each pixel correctly).
  # At the moment, 'pixel.edges' is missing pixel edges along the right and topmost edges of the map area. Once we create our data frame, we will do some manipulation to remedy this.

  # Data frame of information we want
  dat <- data.frame(x=pixel.edges[,1], y=pixel.edges[,2], covtype=rep("D~1", 2500), occasions=rep(nocc, 2500), array_size=rep(array, 2500), value=racd.vals)
  # Manipulating this data frame so it includes pixel edges along the right and topmost edges of the map
  dup1 <- dat[(dat$y==49.5 | dat$x==49.5),]
  save1 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save1$x=50.5
  save2 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save2$y=50.5 # If we don't run these two lines, then we'll miss these two sets of pixel edges in our data frame
  dup1$x[dup1$x==49.5] = 50.5; dup1$y[dup1$y==49.5] = 50.5 # Editing all of the entries in dup1, so that they represent pixel edges along the right and top edges of the map area
  dup1 <- rbind(dup1, save1, save2)
  # Putting everything together
  dat <- rbind(dat, dup1)

  # Returning this data frame
  dat
}

## Function to summarise the data for the EACD maps -- arguments are the same as the function above
eacd.summary <- function(nocc, fig) {
  if (fig==4) {
    array <- "3x3"
  } else {
    if (fig==5) {
      array <- "7x7"
    }
  }
  eacd.vals <- get(paste0("eacd.", nocc, "occ"))

  pixel.centres <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50)
  pixel.edges <- pixel.centres - 0.5

  dat <- data.frame(x=pixel.edges[,1], y=pixel.edges[,2], covtype=rep("D~log(Dblur)", 2500), occasions=rep(nocc, 2500), array_size=rep(array, 2500), value=eacd.vals)
  # Doing the same manipulation as above, so we work with pixel edges to colour our map
  dup1 <- dat[(dat$y==49.5 | dat$x==49.5),]
  save1 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save1$x=50.5
  save2 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save2$y=50.5
  dup1$x[dup1$x==49.5] = 50.5; dup1$y[dup1$y==49.5] = 50.5
  dup1 <- rbind(dup1, save1, save2)
  dat <- rbind(dat, dup1)

  dat
}

## Function to create the final 'predicted_densities_all' object
# Here, 'nocc' is the vector of sampling occasions we are working with; 'fig' is the corresponding figure number for each map; 'type' is the corresponding type of map we want to create (enter as 'RACD' or 'EACD')
overall.summary <- function(nocc, fig, type) {
  # Initialising data frame
  dat <- data.frame()
  for (i in 1:length(nocc)) {
    # If type="RACD", creating a data frame summarising the info for the RACD map corresponding to the given number of sampling occasions and given figure number
    if (type[i]=="RACD") {
      dat.add <- racd.summary(nocc=nocc[i], fig=fig[i])
    } else {
      # If type="EACD", creating a data frame summarising the info for the EACD map
      if (type[i]=="EACD") {
        dat.add <- eacd.summary(nocc=nocc[i], fig=fig[i])
      }
    }
    # Adding the data frame created by racd.summary() or eacd.summary() to our 'dat' data frame
    dat <- rbind(dat, dat.add)
  }
  # Returning the final data frame
  dat
}

## Creating the 'predicted_densities_all' data frame
nocc <- c(rep(c(18, 52, 111), 2), rep(c(7, 25, 55), 2))
fig <- c(rep(4, 6), rep(5, 6))
type <- c(rep(c(rep("RACD", 3), rep("EACD", 3)), 2))
predicted_densities_all <- overall.summary(nocc=nocc, fig=fig, type=type)

##### Object that contains information about the detectors in both figures #####

detectors_df_all <- res_acd %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()

# Saving objects we have created, for use in appendix.Rnw
#save(predicted_densities_all, file="fig4_fig5.RData")
#save(detectors_df_all, file="detectors.RData") # Same as 'detectors.RData' we save in 'UncertaintyPlots_Code.R'

##### Defining the max value of the colour scale #####

## We want this value to be the same for both Figures 4 and 5.

nn <- 3 # Number of simulated datasets we use in each figure
xx <- predicted_densities_all %>% filter(array_size == "3x3", occasions %in% capthists_few_alloccs_3x3$noccasions[1:nn])
maxval1 <- max(xx$value) # Max density value found across all plots in Figure 4
xx <- predicted_densities_all %>% filter(array_size == "7x7", occasions %in% capthists_few_alloccs_7x7$noccasions[1:nn])
maxval2 <- max(xx$value) # Max density value found across all plots in Figure 5
maxval <- max(maxval1, maxval2) # Using the higher of these two max values as the top value of the colour scale for plots in both figures

## ---------------------------------------------------------------------------------------
################################### Creating Figure 4 ####################################
## ---------------------------------------------------------------------------------------

##### Adding the column and row labels for Figure 4 to the 'predicted_densities_all' and 'detectors_df_all' objects #####

nn <- 3 # Number of different simulated datasets used in Figure 4
occ <- capthists_few_alloccs_3x3$noccasions # Number of sampling occasions for each dataset
asz <- c("3x3")
chs <- data.frame(do.call(rbind, lapply(capthists_few_alloccs_3x3$capthist, summary, terse = TRUE)))
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = chs$Detections, .y = chs$Animals, .f = paster) %>% unlist() # Column lables for Figure 4

## Adding the column labels for Figure 4 to the 'predicted_densities_all' and 'detectors_df_all' objects
predicted_densities_all$occasions2 <- factor(predicted_densities_all$occasions,
                                             levels = occ,
                                             labels = capthist_labels)
detectors_df_all$occasions2 <- factor(detectors_df_all$occasions,
                                      levels = occ,
                                      labels = capthist_labels)

## Adding the row labels for Figure 4 to the 'predicted_densities_all' and 'detectors_df_all' objects
predicted_densities_all$covtype2 <- factor(predicted_densities_all$covtype,
                                           levels = unique(predicted_densities_all$covtype),
                                           labels = c("Realised AC", "Expected AC"))
detectors_df_all$covtype2 <- factor(detectors_df_all$covtype,
                                    levels = unique(detectors_df_all$covtype),
                                    labels = c("Realised AC", "Expected AC"))

##### Creating and saving Figure 4 #####

fig4 <- predicted_densities_all %>%
  filter(occasions %in% occ[1:nn], array_size %in% asz) %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  scale_fill_distiller(limits=c(0, maxval)) +
  facet_grid(covtype2 ~ occasions2) +
  geom_point(data = detectors_df_all %>% filter(occasions %in% occ[1:nn], array_size %in% asz), inherit.aes = T,
             colour = "gray80", pch = 4, size = 2) +
  geom_point(data = simulated_points, inherit.aes = F, aes(x=x,y=y),
             colour = "darkorange", pch = 16, size = 1, alpha = 0.5) +
  coord_equal() +
  theme_classic(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.spacing=unit(-1, "lines"),
        strip.background = element_rect(fill=NA, colour = NA),
        legend.position="right", legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(1.3,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

fig4

## Saving Figure 4
ggsave("mona_3x3.png", fig4, width=8, height=6, dpi=600, bg="white")

## ---------------------------------------------------------------------------------------
################################### Creating Figure 5 ####################################
## ---------------------------------------------------------------------------------------

##### Adding the column and row labels for Figure 5 to the 'predicted_densities_all' and 'detectors_df_all' objects #####

nn <- 3
occ <-capthists_few_alloccs_7x7$noccasions
asz <- c("7x7")
chs <- data.frame(do.call(rbind, lapply(capthists_few_alloccs_7x7$capthist, summary, terse = TRUE)))
chs <- chs %>% dplyr::filter(Occasions %in% occ)
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = chs$Detections, .y = chs$Animals, .f = paster) %>% unlist()

## Adding the column labels for Figure 5 to the 'predicted_densities_all' and 'detectors_df_all' objects
predicted_densities_all$occasions2 <- factor(predicted_densities_all$occasions,
                                             levels = occ,
                                             labels = capthist_labels)
detectors_df_all$occasions2 <- factor(detectors_df_all$occasions,
                                      levels = occ,
                                      labels = capthist_labels)

## Adding the row labels for Figure 5 to the 'predicted_densities_all' and 'detectors_df_all' objects
predicted_densities_all$covtype2 <- factor(predicted_densities_all$covtype,
                                           levels = unique(predicted_densities_all$covtype),
                                           labels = c("Realised AC", "Expected AC"))
detectors_df_all$covtype2 <- factor(detectors_df_all$covtype,
                                    levels = unique(detectors_df_all$covtype),
                                    labels = c("Realised AC", "Expected AC"))

##### Creating and saving Figure 5 #####

fig5 <- predicted_densities_all %>%
  filter(occasions %in% occ[1:nn], array_size %in% asz) %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  scale_fill_distiller(limits=c(0, maxval)) +
  facet_grid(covtype2 ~ occasions2) +
  geom_point(data = detectors_df_all %>% filter(occasions %in% occ[1:nn], array_size %in% asz), inherit.aes = T,
             colour = "gray80", pch = 4, size = 2) +
  geom_point(data = simulated_points, inherit.aes = F, aes(x=x,y=y),
             colour = "darkorange", pch = 16, size = 1, alpha = 0.5) +
  theme_bw(base_size = 14) +
  coord_equal() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.spacing=unit(-1, "lines"),
        strip.background = element_rect(fill=NA, colour = NA),
        legend.position="right", legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(1.3,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

fig5

## Saving Figure 5
ggsave("mona_7x7.png", fig5, width=8, height=6, dpi=600, bg="white")
