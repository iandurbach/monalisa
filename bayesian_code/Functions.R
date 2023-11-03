## Functions needed to create the Bayesian plots (to run the MCMC, organise MCMC results, etc.)

# --------------------------------------------------------

## Function to generate MCMC samples when fitting a homogeneous density model

## The arguments we need to provide (in order) are:
# * 'data':  a data object created in a similar manner to the data objects from the beginning of 'bayesian_code/Plots_Code.R'. This is a list that should include elements labelled: 'encounter.data' (capture history matrix), 'trap.loc' (trap coordinates), 'xlim' (x-range of coordinates in map area), 'ylim' (y-range of coordinates in map area) and 'n.occasions' (number of sampling occasions)
# * 'M': the size of the superpopulation
# * 'n.iter': the number of MCMC iterations we want to run
# * 'n.adapt': the number of adaptation iterations we want (different from burn-in)
# * 'n.burn': the number of burn-in iterations that we want to use (these will be discarded for us)
# * 'lambda0.start': the initial value for lambda0 for the MCMC sampling
# * 'log_coeff.start': the initial value for 'log_coeff' for the MCMC sampling (see NIMBLE model below for an idea of what parameter this is)
# * 'thin': the value of the thinning parameter for the MCMC
# * 'parameters': the vector of parameters we want to monitor

run.MCMC <- function(data, M, n.iter=1000, n.adapt = 1000, n.burn = 100, lambda0.start = runif(1, 0, 50), log_coeff.start=runif(1, 0, 1), thin = 1, parameters) {
  ## Subsetting the data we'll use in the NIMBLE model:
  # Encounter matrix
  y <- data$encounter.data
  # Removing rows of zeroes in encounter data (in case they are there)
  y <- y[apply(y, 1, sum)>0, ]
  # Trap locations matrix
  traplocs <- data$trap.loc
  # Number of traps
  trap.no <- nrow(traplocs)
  # Number of sampling occasions over which data was collected
  n.occ <- data$n.occasions
  # xlim
  xlim <- data$xlim
  # ylim
  ylim <- data$ylim
  # Number of animals detected
  pop.size <- nrow(y)

  ## Data augmentation
  # Setting the size of the super-population
  M <- M
  # Adding all-0 rows (so that, in total, we have encounter data for M animals)
  y <- rbind(y, matrix(0, nrow=M-pop.size, ncol=ncol(y)))
  # Vector of 0's and 1's corresponding to our encounter data matrix - 1 if a 'real' individual (a detected individual), 0 for an 'added' individual (i.e. an individual that we are considering, but was never detected at a trap)
  z <- c(rep(1, pop.size), rep(0, M-pop.size))

  ## Setting the starting values for s (activity centres for each animal)
  ## For observed animals, we want these initialised activity centres to be at or near the traps at which individuals were captured. So, for these inviduals, we'll make the starting activity centre the 'mean' location of the traps at which they were caught. We do this below
  # First, generating 'random' activity centres for ALL individuals in the super-population
  sst <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
 # Now, for observed animals, making their initial activity centre the 'mean' of all of the traps at which they were detected
  for (i in 1:pop.size) {
    if (sum(y[i,])==0) {next}
    sst[i,1] <- mean(traplocs[y[i,]>0,1])
    sst[i,2] <- mean(traplocs[y[i,]>0,2])
  }

  ## NIMBLE model
  # We could monitor multiple things here: lambda0, psi, log_coeff, coeff, sigma, z-values and s-values (activity centres). Generally, we enter the argument: parameters = c("lambda0", "coeff", "sigma", "N", "D", "z", "s") -- that is, these are the parameters we tend to monitor!
  x <- nimbleCode({
    lambda0~dgamma(0.001, 0.001)
    psi ~ dunif(0,1)
    log_coeff ~ dunif(-10,10)
    coeff <- exp(log_coeff)
    sigma <- sqrt(1/(2*coeff))

    for (i in 1:M) {
      z[i] ~ dbern(psi)
      s[i,1] ~ dunif(xlim[1], xlim[2])
      s[i,2] ~ dunif(ylim[1], ylim[2])
      for (j in 1:trap.no) {
        d[i,j] <- sqrt((s[i,1] - traplocs[j,1])^2 + (s[i,2] - traplocs[j,2])^2)
        lambda[i,j] <- z[i] * lambda0 * exp(-coeff * d[i,j] * d[i,j])
        y[i,j] ~ dpois(lambda[i,j]*n.occ)
      }
    }

    N <-  sum(z[1:M])
    D <-  (N/((xlim[2] - xlim[1]) * (ylim[2] - ylim[1]))) * 10000
  })

  ## Data to provide to NIMBLE
  nim.data <- list(y=y, traplocs=traplocs)

  ## Constants to provide to NIMBLE
  constants <- list(n.occ=n.occ, M=M, trap.no=trap.no, xlim=xlim, ylim=ylim)

  ## Initial values for the MCMC
  inits <- list(lambda0=lambda0.start, log_coeff=log_coeff.start, s=sst, z=z)

  ## Parameters to monitor
  parameters <- parameters

  ## Running NIMBLE
  Rmodel <- nimbleModel(code=x, constants=constants, data=nim.data, inits=inits)
  conf <- configureMCMC(Rmodel, monitors=parameters, control = list(adaptInterval = n.adapt))
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

  ## Running the MCMC, generating and saving the final results
  results <- runMCMC(Cmcmc, niter=n.iter+n.burn, nburnin=n.burn, progressBar=TRUE, samplesAsCodaMCMC=T)
}

# --------------------------------------------------------

## Function to generate MCMC samples when fitting an inhomogeneous density model

## The arguments we need to provide (in order) are:
# * 'data': a data object, created so that it is a list containing the elements: 'encounter.data' (capture history matrix) and 'trap.loc' (trap coordinates) -- see 'bayesian_code/Plots_Code.R' for an example of the creation of such data objects
# * 'pixel.info': a data frame  with three columns. The first column gives the x-coordinates for pixel centres in the region of interest, the second column gives the y-coordinates of the pixel centres, and the third gives the associated covariate value for each pixel centre. NOTE we assume that: (1) these pixel centres are evenly-spaced, (2) the region of interest is square, so the number of pixels in the x- and y-directions is the square root of the number of rows in this data frame, (3) there are no areas in the region of interest where animals cannot go (so a mask of the study region will just be a matrix of 1's) and (4) the first two columns (of pixel centres) have been generated using the centres() function (which can be found below)
# * 'M': the size of the super-population
# * 'inits.vec': a vector containing the starting values for lambda0, sigma and beta1. The ordering is: c(lambda0, sigma, beta1). Later, we calculate sensible starting values for 'log_coeff' and 'beta0', so we won't provide them here
# * 'dmax': the dmax value to use for the getLocalObjects() function (see below)
# * 'n.iter': the number of iterations to run the MCMC for
# * 'n.burn': the number of burn-in iterations we want to use
# * 'parameters': a vector containing the labels of the parameters we want to monitor

run.MCMC.inhom <- function(data, pixel.info, M, inits.vec, dmax = 56, n.iter, n.burn, parameters=c("lambda0", "sigma", "N", "D", "beta0", "beta1")) {
  ## Therefore, subsetting the data we'll use in our NIMBLE model:
  # Encounter data
  y <- data$encounter.data
  # Checking if there are any rows of 0's -- if there are, returning an error because these capture histories are unobserveable
  if (any(apply(y,1,sum)==0)) stop("The data shouldn't include any all-0 capture histories")

  # Trap locations matrix
  traplocs <- data$trap.loc
  # Number of traps
  n.trap <- nrow(traplocs)
  # xlim, ylim (found based on centres in 'pixel.info' being evenly-spaced)
  xlim <- c(min(pixel.info[,1])-(0.5*abs(pixel.info[1,1] - pixel.info[2,1])), max(pixel.info[,1])+(0.5*abs(pixel.info[1,1] - pixel.info[2,1])))
  ylim <- c(min(pixel.info[,2])-(0.5*abs(pixel.info[1,1] - pixel.info[2,1])), max(pixel.info[,2])+(0.5*abs(pixel.info[1,1] - pixel.info[2,1])))
  # Number of animals detected
  n.observed <- nrow(y)

  # Number of pixels in the map region
  nPix <- nrow(pixel.info)

  # Pixel centres for the study region
  pixel.centres <- pixel.info[,1:2]
  # Creating a matrix that has an entry for each pixel centre in the survey region. Each entry will indicate the order in which that pixel occurs in the 'pixel.centres' object (which we assume is created by the centres() function found below)
  # Note that if we start at the bottom row and read each row from left to right, the indices increase in value
  pixel.centres.order <- matrix(1:nrow(pixel.info), ncol=sqrt(nPix), nrow=sqrt(nPix), byrow=T)
  pixel.centres.order <- pixel.centres.order[nrow(pixel.centres.order):1,]

  # Area of each pixel we are considering, calculated based on centres in 'pixel.info' being evenly spaced
  pixel.area <- abs(pixel.info[1,1] - pixel.info[2,1])^2


  ## Data augmentation
  # Setting the size of the super-population
  M <- M
  # Adding all-0 rows to our encounter data matrix (so that, in total, we have encounter data for M animals)
  y <- rbind(y, matrix(0, nrow=M-n.observed, ncol=ncol(y)))
  # Vector of 0's and 1's corresponding to our encounter data matrix - 1 if a 'real' individual (a detected individual), 0 for an 'added' individual (i.e. an individual that we are considering, but was never detected at a trap)
  z <- c(rep(1, n.observed), rep(0, M-n.observed))


  ## Setting the starting values for s (activity centres for each animal)
  ## For observed animals, we want these initialised activity centres to be at or near the traps at which individuals were captured. So, for these individuals, we'll make the starting activity centre the 'mean' location of the traps at which they were caught. We do this below
  # First, generating 'random' activity centres for ALL individuals in the super-population
  sst <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
  # Now, for observed animals, making their initial activity centre the 'mean' of all of the traps at which they were detected
  for (j in 1:n.observed) {
    if (sum(y[j,])==0) {next}
    sst[j,1] <- mean(traplocs[y[j,]>0,1])
    sst[j,2] <- mean(traplocs[y[j,]>0,2])
  }
  # And now, rounding these starting values so that they correspond to pixel centres
  sst <- round(sst)
  # Finding the corresponding pixel index for each row of 'sst' - i.e. for each pixel centre we have identified, we are finding the corresponding row index in the 'pixel.centres' object. This will be equal to the entry in 'pixel.centres.order' that represents the given pixel centre
  starting.pixel.indices <- vector("numeric", M)
  for (j in 1:M) {
    index <- which(sst[j,1]==pixel.centres[,1] & sst[j,2]==pixel.centres[,2])
    starting.pixel.indices[j] <- index
  }
  # And now, making sst equal to the 'starting.pixel.indices' object we have just created
  sst <- starting.pixel.indices
  # Finding sx and sy -- these are the row/column indices for each entry in 'pixel.centres.order' that is stored in 'sst'
  sx_sy_init <- matrix(0, ncol=2, nrow=length(sst))
  for (i in (1:length(sst))) {
    sx_sy_init[i,] <- which(pixel.centres.order==sst[i], arr.ind=T)
  }
  # Subsetting the initial values for sx
  sx_init <- sx_sy_init[,1]
  # Initial values for sy
  sy_init <- sx_sy_init[,2]
  # What we end up providing to our NIMBLE model in the way of initial activity centres is 'sst', 'sx' and 'sy'. So, we provide: the index of the pixel in which each animals' initial activity centre falls into (based on the 'pixel.centres' object), and the row/column indices that can be used to identify this pixel in the 'pixel.centres.order' object, respectively.

  # Covariate values we want to use in our model
  mona.densities <- pixel.info[,3]

  # As we want to use dpoisLocal_normal() below, we  need to scale the trap coordinates using the scaleCoordsToHabitatGrid() function
  # Before this, need to label the columns in 'pixel.centres' as 'x' and 'y'
  colnames(pixel.centres) <- c("x", "y")
  # Scaling the trap coordinates
  scaledtrapcoords <- scaleCoordsToHabitatGrid(coordsData = traplocs, coordsHabitatGridCenter = pixel.centres)
  scaledtrapcoords <- scaledtrapcoords$coordsDataScaled
  # Scaling the pixel centres, as well
  scaledpixelcentres <- scaleCoordsToHabitatGrid(coordsData = pixel.centres, coordsHabitatGridCenter = pixel.centres)
  scaledpixelcentres <- scaledpixelcentres$coordsDataScaled


  ## Defining the NIMBLE model
    code <- nimbleCode({
    # Priors -- will see that these are the same as for the homogeneous PP model (for the parameters that are common to both models)
    lambda0 ~ dgamma(0.001, 0.001)
    log_coeff ~ dunif(-10, 10)
    coeff  <- exp(log_coeff)
    sigma <- sqrt(1/(2*coeff))
    # (The parameters below are unique to the inhomogeneous PP model)
    beta0 ~ dunif(-50, 50)
    beta1 ~ dunif(-10, 10)

    # Specifying prior probabilities for each pixel
    DPix[1:nPix] <- exp(beta0 + beta1*(mona.densities[1:nPix]))
    mu[1:nPix] <- DPix[1:nPix] * pixel.area
    probs[1:nPix] <- mu[1:nPix]/EN

    EN <- sum(mu[1:nPix])  # Expected value of N, E(N)
    psi <- EN/M  # Data augmentation parameter

    for (k in 1:M) {
      z[k] ~ dbern(psi) # Whether or not the ith animal exists
      sx[k] ~ dunif(1, 51) # Prior for column of matrix that represents where in 'pixel.centres.order' the animal's activity centre lies
      sy[k] ~ dunif(1, 51) # Prior for row of matrix that represents where in 'pixel.centres.order' the animal's activity centre lies
      ind_x[k] <- trunc(sx[k]) # Finding the row of 'pixel.centres.order' that contains the sampled activity centre
      ind_y[k] <- trunc(sy[k]) # Finding the column of 'pixel.centres.order' that contains the sampled activity centre
      s[k] <- pixel.centres.order[ind_x[k], ind_y[k]] # Finding the corresponding entry in 'pixel.centres.order'. This is the row index of the pixel that contains the sampled activity centre, relative to the 'pixel.centres' object (so if we extracted the row of the 'pixel.centres' object that corresponds to the given index, we would have the centre of the pixel that the sampled activity centre falls into)
      ones[k] ~ dbern(probs[s[k]]) # Adds to the likelihood the probability of our sampled activity centre falling into the given pixel

      # Likelihood, based on the Wolverine NIMBLE example found at: https://nimble-dev.github.io/nimbleSCR/wolverine_example.html
      y[k, 1:nMaxDetectors] ~ dpoisLocal_normal(detNums = nbDetections[k],
                                                  detIndices = yDets[k, 1:nMaxDetectors],
                                                  lambda = lambda0,
                                                  s = scaledpixelcentres[s[k],1:2], # s[k] is the index of the pixel centre corresponding to the sampled activity centre -- so here, we are extracting the row of 'scaledpixelcentres' that gives us the scaled version of this pixel centre
                                                  sigma = sigma,
                                                  trapCoords = scaledtrapcoords[1:n.trap, 1:2],
                                                  localTrapsIndices = detectorIndex[1:n.cells,1:maxNBDets],
                                                  localTrapsNum = nDetectors[1:n.cells],
                                                  resizeFactor = ResizeFactor,
                                                  habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                  indicator = z[k])
    }
    N <- sum(z[1:M])  # Realised value of N
    D <- (N/((xlim[2] - xlim[1]) * (ylim[2] - ylim[1]))) * 10000
  }
  )

  ## Values that we want to provide to our NIMBLE model
  # Data
  data <- list(scaledtrapcoords = scaledtrapcoords,
               scaledpixelcentres = scaledpixelcentres,
               mona.densities = mona.densities,
               ones = rep(1, M),
               pixel.centres.order = pixel.centres.order)
  # Constants
  constants <- list(nPix=nPix, M=M, n.trap=n.trap, xlim=xlim, ylim=ylim,
                    pixel.area=pixel.area)
  # Initial values. We are setting the starting value for 'log_coeff' based on the starting value of sigma, and here we try to set a sensible starting value for beta0
  inits <- list (lambda0=inits.vec[1], log_coeff=log(1/(2*(inits.vec[2])^2)), z=z, beta0=log(n.observed/(nPix*pixel.area*10000)), beta1=inits.vec[3], s=sst, sx=sx_init, sy=sy_init)

  ## More constants that we want to provide to the NIMBLE model
  # And since we assume we have no unsuitable habitat, we are defining the 'habitatMask' argument as a matrix full of 1's. Note that we are using the largest value of dmax possible (see documentation of getLocalObjects() for more info)
  DetectorIndex <- getLocalObjects(habitatMask = matrix(1, ncol=50, nrow=50), coords = scaledtrapcoords, dmax = 52, resizeFactor = 1, plot.check=FALSE)
  # Generating the values of more constants that we want to provide to our NIMBLE model
  constants$y.maxDet <- dim(DetectorIndex$habitatGrid)[1]
  constants$x.maxDet <- dim(DetectorIndex$habitatGrid)[2]
  constants$ResizeFactor <- DetectorIndex$resizeFactor
  constants$n.cells <- nPix
  constants$maxNBDets <- DetectorIndex$numLocalIndicesMax
  data$detectorIndex <- DetectorIndex$localIndices
  data$nDetectors <- DetectorIndex$numLocalIndices
  data$habitatIDDet <- DetectorIndex$habitatGrid

  ## More data that we want to provide to the NIMBLE model
  # Generating a sparse representation of encounter data matrix
  ySparse <- getSparseY(x = y)
  data$y <- ySparse$y[,,1]
  data$yDets <- ySparse$detIndices[,,1]
  data$nbDetections <- ySparse$detNums[,1]
  constants$nMaxDetectors <- ySparse$maxDetNums

  ## Doing the clever NIMBLE stuff
  Rmodel <- nimbleModel(code, constants, data, inits, dimensions = list(pixel.centres.order = c(50,50)))
  # AF slice sampling for beta0 and beta1
  conf <- configureMCMC(Rmodel, monitors = parameters, print = FALSE)
  conf$removeSampler(c("beta0","beta1"))
  conf$addSampler(target = c("beta0","beta1"),
                  type = 'AF_slice',
                  control = list(adaptScaleOnly = TRUE),
                  silent = TRUE)
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  samples <- runMCMC(Cmcmc, niter = n.iter+n.burn, progressBar=TRUE)

  # Returning the MCMC samples
  samples
}

# --------------------------------------------------------

## Function to check the trace plots resulting from our MCMC samples
# The 'results' argument represents an MCMC sample where the parameters 'lambda0', 'sigma', 'N' and 'D' have been monitored
# If 'inhom'=TRUE, then we will also check trace plots for parameters labelled beta0 and beta1

check.trace.plots <-  function(results, inhom=FALSE) {
  if (inhom) {
    par(mfrow=c(3,2))
    plot(as.vector(results[,"lambda0"]), type='l', ylab=expression(lambda[0]))
    plot(as.vector(results[,"sigma"]), type='l', ylab=expression(sigma))
    plot(as.vector(results[,"N"]), type='l', ylab="N")
    plot(as.vector(results[,"D"]), type='l', ylab="D")
    plot(as.vector(results[,"beta0"]), type='l', ylab=expression(beta[0]))
    plot(as.vector(results[,"beta1"]), type='l', ylab=expression(beta[1]))
  } else {
    par(mfrow=c(2, 2))
    plot(as.vector(results[,"lambda0"]), type='l', ylab=expression(lambda[0]))
    plot(as.vector(results[,"sigma"]), type='l', ylab=expression(sigma))
    plot(as.vector(results[,"N"]), type='l', ylab="N")
    plot(as.vector(results[,"D"]), type='l', ylab="D")
  }
}

# --------------------------------------------------------

## Note that this function was written explicitly to work with the covariate used in the paper

## Function to generate the covariate values that we use when putting together the EACD plots. When running the function, the 'mona-inputs.RData' file must already be loaded into R.
eacd.covariate <- function() {
  # Subsetting the covariate values from the data loaded in from the 'mona-inputs.RData' file
  mona.densities <-  small_blurry_mona_df[,c("x", "y", "Dblur")]
  # Re-ordering 'mona.densities', so order of pixels matches order of pixels in 'pixel.centres' object (see above for creation of 'pixel.centres' object)
  split <-  split(mona.densities, mona.densities$y)
  mona.densities <-  do.call("rbind", split)
  rownames(mona.densities) = NULL
  # Now, subsetting "Dblur" vector only so is in corresponding order to centres in 'pixel.centres'
  dblur <-  mona.densities[,"Dblur"]
  # Logging the covariate, so we have the values of log(Dblur) (this is the covariate we will use to fit our SCR models)
  log.dblur <-  log(dblur)
  log.dblur
}

# --------------------------------------------------------

## Function to create the density values for an RACD map

# Here, 'xlim' and 'ylim' give the range of x- and y-coordinates for the map area. Also, 'results' refers to a set of MCMC samples generated using run.MCMC() and 'M' is the size of the super-population
# We are assuming that each pixel has an area of 1, so that 'xg' and 'yg' below contain the pixel centres. Then, we have that '(length(xg) - 1) * (length(yg) - 1)' below gives the number of pixels we are using in the map area

racd.density.vector <- function(xlim, ylim, results, M) {

  ## x- and y-range of coordinates in map area
  xg <- seq(xlim[1], xlim[2], by=1)
  yg <- seq(ylim[1], ylim[2], by=1)

  ## Extracting z-values
  # Names of variables that have been monitored
  names <- names(results[1,])
  # Extracting "z" values from MCMC results
  Z <- results[,grep("z", names)]

  ## Extracting activity centres
  # Extracting "s" values from MCMC results (i.e. extracting sampled activity centres)
  S <- results[,grep("s[^i]", names)]
  # x-coordinates of all activity centres
  Sx <- S[,1:M]
  # y-coordinates of all activity centres
  Sy <- S[,-(1:M)]

  ## For each MCMC iteration, storing the number of animals that 'exist' and have a sampled activity centre in each pixel -- so, we are building up the posterior distribution of the number of activity centres in each pixel
  # Number of pixel centres we are working with
  npix <-  (length(xg) - 1) * (length(yg) - 1)
  Dn.vals <-  matrix(0, nrow=nrow(results), ncol=npix)
  # For loop to find the posterior distributions mentioned above
  for (i in 1:nrow(results)) {
    if ((i %% 100) == 0) print(i) # Track progress
    Sxout <-  Sx[i,][Z[i,] == 1]
    Sxout <-  cut(Sxout, breaks=xg, include.lowest=TRUE)
    Syout <-  Sy[i,][Z[i,] == 1]
    Syout <-  cut(Syout, breaks=yg, include.lowest=TRUE)
    Dn.vals[i,] <-  as.vector(table(Sxout, Syout))
  }

  ## Posterior mean for number of activity centres in each pixel
  density.vector <- apply(Dn.vals, 2, mean)
  # Note that if our pixel area were not one, we would need to divide the posterior means here by the area of each pixel if we want an estimate of the number of activity centres per unit for each pixel
}

# --------------------------------------------------------

## Function to generate pixel centres across map area

# 'xlim' and 'ylim' are described as above. In addition, 'x.pixels' and 'y.pixels' give the number of pixels being used in the x- and y-direction.

library("spatstat")
centres <- function(xlim, ylim, x.pixels, y.pixels) {
  # Creating an object of class 'owin' representing our map area
  window.2 <- owin(xrange=xlim, yrange=ylim)
  # Generating our set of pixel centres
  points <- gridcentres(window.2, x.pixels, y.pixels)
  # Converting the result to a matrix
  centres <- as.matrix(cbind(points$x, points$y))
  # Printing the result
  centres
}

# --------------------------------------------------------

## See the paper for more information on the equation we use to calculate animal density in each pixel here

## Function to create vector of density values for the EACD maps we want to create. The arguments are:
# * 'results': an MCMC object, which contains samples for the parameters 'beta0' and 'beta1'
# * 'covariate': a vector containing the covariate values for each pixel (we are assuming the values are already logged)
# * 'nPix': the number of pixels that we are working with
eacd.density.vector <- function(results, covariate, nPix) {
  density <- numeric() # Initialising the object we will store density values in
  # For loop to calculate posterior mean density value for each pixel
  for (i in 1:nPix) {
    if(i%%100==0) print(i) # Track progress
    # Posterior distribution for the density of the ith pixel
    density.posterior <- exp(results[,'beta0'] + results[,'beta1'] * (covariate[i]))
    # Posterior mean of this density -- gives the density value for the ith pixel we will use in an EACD map
    density[i] <- mean(density.posterior)
  }
  # Returning the resulting vector
  density
}

# --------------------------------------------------------

## Note that this function was written explicitly to work with the simulated data used in the paper

## Function to compare the results we get from:
# (1) Fitting an SCR model with an inhomogeneous Poisson process as the state process, using run.MCMC.inhom(). We assume that the covariate used is log(Dblur) (see 'bayesian_code/Plots_Code.R' for info on this covariate)
# (2) Fitting the same model using secr.fit() from the 'secr' package
## The arugments are:
# * 'results': the MCMC samples produced using run.MCMC.inhom(). We assume that the parameters: 'beta0', 'beta1', 'lambda0' and 'sigma' have been monitored
# * 'mask': mask to use in secr.fit()
# * 'array': the name of the array we want to use. If working with Figure 4, array should be "3x3". If working with Figure 5, array should be "7x7".
# * 'j': this will be a value in {1, 2, 3} and represents the index of the objects we want to work with from the RData objects we need to have loaded in (to load in, use load("../output/mona-inputs.RData')). If array="3x3" and we want objects generated using 18, 52 or 111 sampling occasions, use j=1,2,3, respectively. If array="7x7" and we want objects generated using 7, 25 or 55 sampling occasions, use j=1,2,3, respectively.
# Note that the model we specify for secr.fit() is not an argument here -- we will be using "D~log(Dblur)" as the model every time we use this function. Same for the detection function -- for our purposes, we will always use 'detectfn=HHN'
check.inhom.mcmc <- function(results, j, mask, array) {
  if (array=="3x3") {
    # Number of sampling occasions we are working with
    nocc  <- capthists_few_alloccs_3x3$noccasions[j]
    # Capture history we want to work with
    capthist <- capthists_few_alloccs_3x3$capthist[[j]]
  } else {
    if (array=="7x7") {
      # Number of sampling occasions we are working with
      nocc  <- capthists_few_alloccs_7x7$noccasions[j]
      # Capture history we want to work with
      capthist <- capthists_few_alloccs_7x7$capthist[[j]]
    }
  }

  # Fitting the SCR model using secr.fit()
  fit <- secr.fit(capthist=capthist, model=D~log(Dblur), mask=mask, detectfn="HHN", trace=FALSE)

  # Comparing estimates of beta0
  mcmc.beta0 <- mean(results[,'beta0'])
  secr.fit.beta0 <- log(exp(summary(fit)$coef["D","beta"])/10000)

  # Comparing credible interval, confidence interval for beta1
  mcmc.beta1 <- quantile(results[,'beta1'], probs=c(0.025, 0.5, 0.975)) # 95% credible interval
  secr.fit.beta1 <- c(summary(fit)$coef["D.log(Dblur)","lcl"], summary(fit)$coef["D.log(Dblur)","ucl"]) # 95% confidence interval

  # Comparing credible interval, confidence interval for lambda0
  mcmc.lambda0 <- quantile(results[,'lambda0']/nocc, probs=c(0.025, 0.5, 0.975)) # 95% credible interval
  secr.fit.lambda0 <- c(summary(fit)$predicted["lambda0","lcl"], summary(fit)$predicted["lambda0","ucl"]) # 95% confidence interval

  # Comparing credible interval, confidence interval for sigma
  mcmc.sigma <- quantile(results[,'sigma'], probs=c(0.025, 0.5, 0.975)) # 95% credible interval
  secr.fit.sigma <- c(summary(fit)$predicted["sigma","lcl"], summary(fit)$predicted["sigma","ucl"]) # 95% confidence interval

  # Printing the results
  # Estimates of beta0
  cat("Estimates of beta0 \n", "MCMC estimate: ", mcmc.beta0, "\n", "secr.fit() estimate: ", secr.fit.beta0, "\n", sep="")
  # Intervals for beta1
  cat("\n", "Intervals for beta1", "\n", "95% credible interval: (", mcmc.beta1[1], ", ", mcmc.beta1[3], ") ", "\n", "95% confidence interval: (", secr.fit.beta1[1], ", ", secr.fit.beta1[2], ")", "\n", sep="")
  # Intervals for lambda0
  cat("\n", "Intervals for lambda0", "\n", "95% credible interval: (", mcmc.lambda0[1], ", ", mcmc.lambda0[3], ") ", "\n", "95% confidence interval: (", secr.fit.lambda0[1], ", ", secr.fit.lambda0[2], ")", "\n", "True value of lambda0: 0.1", "\n", sep="")
  # Intervals for sigma
  cat("\n", "Intervals for sigma", "\n", "95% credible interval: (", mcmc.sigma[1], ", ", mcmc.sigma[3], ") ", "\n", "95% confidence interval: (", secr.fit.sigma[1], ", ", secr.fit.sigma[2], ")", "\n", "True value of sigma=4", "\n", sep="")
}

# --------------------------------------------------------

## Note that the functions below were written explicitly to work with the MCMC results, covariate, etc. used in the paper

## Functions used to create the data objects that are used to put together the uncertainty figures in the Appendix

## Function to find values we will use in quantile plots
## Arguments are:
# * 'results': object containing MCMC results
# * 'covariate': vector containing covariate values for each pixel
# * 'nPix': number of pixels in survey area
# * 'quantile': scalar giving the quantile that we wish to calculate when looking at the posterior distn of the density for each pixel
eacd.quantile.plots <- function(results, covariate, nPix, quantile) {
  # Initialising vector we will store values in
  quantile.vals <- numeric()
  # For loop to calculate quantile value for posterior distribution for each pixel
  for (i in 1:nPix) {
    if(i%%100==0) print(i) # Track progress
    # Posterior distribution for the activity centre density in the ith pixel
    density.posterior <- exp(results[,'beta0'] + results[,'beta1'] * (covariate[i]))
    # Finding specified quantile for the posterior distribution for this pixel
    quantile.vals[i] <- quantile(density.posterior, probs=quantile)
  }
  # Returning this vector of values
  quantile.vals
}

## Function that creates data frames we use to create the data frames we require for our EACD uncertainty plots (below). The main purpose of this function is to manipulate the pixel centres we work with so that our final data frame for each individual plot contains pixel edges for the whole map area, so that the resulting plot is filled in correctly. (We want to colour pixels from pixel edge to pixel edge, rather than from pixel centre to pixel centre. The manipulation we do below achieves this -- if we don't do the manipulation, we'll apply colours from pixel centre to pixel centre, so each pixel won't be coloured correctly).
## Note that this function is specifically written to work with the pixel centres generated using the function call: centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50)
## Arguments are:
# * 'pixel.values': the values that we want to use to colour each pixel in our EACD uncertainty plot, we expect these will be produced using the eacd.quantile.plots() or eacd.density.vector() function (so the order of these values will match the order of the pixel centres we generate below)
# * 'nocc': a scalar indicating the number of sampling occasions we want to work with
eacd.quantile.df <- function(pixel.values, nocc) {
  # Array size we are working with
  if (nocc %in% c(18, 52, 111)) {
    array <- "3x3"
  } else {
    array <- "7x7"
  }
  # Pixel centres
  pixel.centres <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50)
  # Subtracting 0.5, so now contains pixel edges (excluding rightmost and topmost edges of map)
  pixel.edges <- pixel.centres - 0.5

  # Putting together a data frame containg information for our EACD uncertainty plot
  dat <- data.frame(x=pixel.edges[,1], y=pixel.edges[,2],
                    covtype=rep("D~log(Dblur)", 2500),
                    occasions=rep(nocc, 2500),
                    array_size=rep(array, 2500),
                    value=pixel.values)
  # Extracting and editing entries in this data frame, so that there are entries representing the topmost and rightmost edges of our map area (otherwise, these edges are left out).
  dup1 <- dat[(dat$y==49.5 | dat$x==49.5),]
  save1 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save1$x=50.5
  save2 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save2$y=50.5 # If we don't run these two lines, then we'll miss these two sets of pixel edges in our data frame
  dup1$x[dup1$x==49.5] = 50.5; dup1$y[dup1$y==49.5] = 50.5 # Editing all of the entries in dup1, so that they represent pixel edges along the right and top edges of the map area
  dup1 <- rbind(dup1, save1, save2)
  # Putting everything together
  dat <- rbind(dat, dup1)

  # Returning the data frame
  dat
}


## Function to calculate the 0.05 quantile values, EACD values and 0.95 quantile values for each pixel in our survey area, for a specified number of sampling occasions. The results will be collated and returned as a data frame, which will also include additional information we require for our uncertainty figure (including pixel edges, array size, number of sampling occasions, etc.)
## When running this function, we assume that the MCMC results for the given number of sampling occasions are already loaded into R (with burn-in already discarded).
## Arguments are:
# * 'results', 'covariate', 'nPix', 'nocc': same as above
eacd.quantile.info <- function(results, covariate, nPix, nocc) {
  # Retrieving object containing MCMC results we wish to use
  mcmc.results <- get(paste0("inhom.results.", nocc, "occ"))
  # Finding 5% quantile values for each pixel
  lowerq <- eacd.quantile.plots(results=mcmc.results, covariate=covariate, nPix=nPix, quantile=0.05)
  # EACD values for each pixel
  eacd <- eacd.density.vector(results=mcmc.results, covariate=covariate, nPix=nPix)
  # 95% quantile values for each pixel
  upperq <- eacd.quantile.plots(results=mcmc.results, covariate=covariate, nPix=nPix, quantile=0.95)

  # Putting together a data frame containing info for the 0.05 quantile plot
  dat.lowerq <- eacd.quantile.df(lowerq, nocc)
  # Doing the same for the EACD plot
  dat.eacd <- eacd.quantile.df(eacd, nocc)
  # And for the 0.95 quantile
  dat.upperq <- eacd.quantile.df(upperq, nocc)

  # Putting together and returning the final data frame
  dat <- rbind(dat.lowerq, dat.eacd, dat.upperq)
  dat
}

# --------------------------------------------------------

## Function we use to calculate the CV values in the CV plots shown in the Appendix OR standard deviation values in the standard deviation plots in the Appendix. Arguments are:
# * xlim, ylim: range of x- and y-coordinates in map area
# * results: MCMC results, created using run.MCMC()
# * M: superpopulation size
# * pixel.index: if 1 we are working with 1x1 pixels; if 2 we are working with 2x2 pixels, if 5 we are working with 5x5 pixels
# * cv: if TRUE, will return CV values (posterior standard deviation/posterior mean) for each pixel, if FALSE will return posterior standard deviation for each pixel
uncertainty.values.racd <- function(xlim, ylim, results, M, pixel.index, cv) {
  ## Range of x- and y-coordinates in map area
  xg <- seq(xlim[1], xlim[2], by=pixel.index)
  yg <- seq(ylim[1], ylim[2], by=pixel.index)

  ## Extracting z-values
  # Names of variables that have been monitored
  names <- names(results[1,])
  # Extracting "z" values from MCMC results
  Z <- results[,grep("z", names)]

  ## Extracting activity centres
  # Extracting "s" values from MCMC results (i.e. extracting sampled activity centres)
  S <- results[,grep("s[^i]", names)]
  # x-coordinates of all activity centres
  Sx <- S[,1:M]
  # y-coordinates of all activity centres
  Sy <- S[,-(1:M)]

  ## For each MCMC iteration, storing the number of animals alive and with their activity centres in each cell -- building up posterior distribution of number of activity centres in each pixel
  # Number of pixels we are working with
  npix <- 2500/(pixel.index^2)
  Dn.vals <-  matrix(0, nrow=nrow(results), ncol=npix)
  for (i in 1:nrow(results)) {
    if ((i %% 100) == 0) print(i) # Track progress
    Sxout <-  Sx[i,][Z[i,] == 1]
    Sxout <-  cut(Sxout, breaks=xg, include.lowest=TRUE)
    Syout <-  Sy[i,][Z[i,] == 1]
    Syout <-  cut(Syout, breaks=yg, include.lowest=TRUE)
    Dn.vals[i,] <-  as.vector(table(Sxout, Syout))
  }

  # Dividing total number of act cent in each aggreggated pixel by new pixel area, which is 4 -- this gives us the posterior we want to work with!
  Dn.vals <- Dn.vals/(pixel.index^2)

  ## Posterior mean for number of activity centres in each pixel
  posterior.mean <- apply(Dn.vals, 2, mean)
  ## Posterior standard deviation
  standard.deviation <- apply(Dn.vals, 2, sd)

  if (cv) {
    ## CV values
    cv.values <- (standard.deviation/posterior.mean) * 100
    # Returning the CV values
    cv.values
  } else {
    ## Returning posterior standard deviation values
    standard.deviation
  }
}

