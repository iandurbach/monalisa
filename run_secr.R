run_secr <- function(simulated_points,
                     secr.fitformula, 
                     mask,
                     nx, ny, # number of traps in x and y dim
                     dx, dy, # distance between traps in x and y dim
                     xorig, yorig, # origin of trap array
                     sigma, lambda0, # detection function
                     noccasions, 
                     capthist,
                     my.seed = sample(1:100000, 1)){
  
  print("next run")
  
  n_pts <- nrow(simulated_points) 
  
  # placeholders for inputs and results
  detectors_df <- data.frame(x = as.integer(), 
                             y = as.integer(), 
                             covtype = as.character(),
                             occasions = as.integer(),
                             array_size = as.character(),
                             array_spacing = as.character(),
                             array_origin = as.character(),
                             lambda0 = as.numeric(),
                             sigma = as.numeric(),
                             n_pts = as.integer())
  
  predicted_densities <- data.frame(x = as.integer(), 
                                    y = as.integer(), 
                                    prob_ac = as.numeric(),
                                    expnumber_ac = as.numeric(),
                                    covtype = as.character(),
                                    occasions = as.integer(),
                                    array_size = as.character(),
                                    array_spacing = as.character(),
                                    array_origin = as.character(),
                                    lambda0 = as.numeric(),
                                    sigma = as.numeric(),
                                    n_pts = as.integer())
  
  estimated_sigma <- data.frame(est_sigma = as.numeric(),
                                covtype = as.character(),
                                occasions = as.integer(),
                                array_size = as.character(),
                                array_spacing = as.character(),
                                array_origin = as.character(),
                                lambda0 = as.numeric(),
                                sigma = as.numeric(),
                                n_pts = as.integer())
  
  # make detector array
  detectors <- make.grid(nx = nx, ny = ny, spacex = dx, spacey = dy,
                         originxy = c(xorig, yorig), detector = "count")
  
  
  # fit model specified by secr.fitformula
  cfit <- secr.fit(capture_history, 
                   model = list(as.formula(secr.fitformula)),  
                   mask = mask, detectfn = "HHN", 
                   # starting value for D = n_pts * 4 as 50x50 is 0.25ha
                   start = list(D = n_pts * 4, lambda0 = lambda0, sigma = sigma),
                   ncores = 1, method = "Nelder-Mead", trace = TRUE)
  
  # extract predicted AC location surface from constant density models, and expected AC density surface
  # from covariate models
  # note: if there is any "1" in secr.fitformula then assume its a "no covariate" model
  if (str_detect(secr.fitformula, "1")) { 
    fxt <- fx.total(cfit)
    prob_ac <-  covariates(fxt)$D.sum
    # not used
    expnumber_ac <- prob_ac * region.N(cfit)["R.N","estimate"]
  } else {
    prob_ac <- covariates(predictDsurface(cfit))$D.0
    # not used
    expnumber_ac <- prob_ac * region.N(cfit)["R.N","estimate"]
  }
  
  # append results
  predicted_densities <- rbind(predicted_densities,
                               data.frame(x = mask$x, y = mask$y, 
                                          prob_ac = prob_ac,
                                          expnumber_ac = expnumber_ac,
                                          covtype = secr.fitformula,
                                          occasions = i,
                                          array_size = paste0(nx,"x",ny),
                                          array_spacing = paste0(dx,"_",dy),
                                          array_origin = paste0(xorig,"_",yorig),
                                          lambda0 = lambda0,
                                          sigma = sigma,
                                          n_pts = length(mask$x)))
  
  detectors_df <- rbind(detectors_df,
                        data.frame(x = detectors$x, y = detectors$y, 
                                   covtype = secr.fitformula,
                                   occasions = i,
                                   array_size = paste0(nx,"x",ny),
                                   array_spacing = paste0(dx,"_",dy),
                                   array_origin = paste0(xorig,"_",yorig),
                                   lambda0 = lambda0,
                                   sigma = sigma,
                                   n_pts = length(mask$x)))
  
  estimated_sigma <- rbind(estimated_sigma, 
                           data.frame(est_sigma = ifelse(str_detect(secr.fitformula,"1"),
                                                         exp(cfit$fit$par[3]),
                                                         exp(cfit$fit$par[4])),
                                      covtype = secr.fitformula,
                                      occasions = i,
                                      array_size = paste0(nx,"x",ny),
                                      array_spacing = paste0(dx,"_",dy),
                                      array_origin = paste0(xorig,"_",yorig),
                                      lambda0 = lambda0,
                                      sigma = sigma,
                                      n_pts = length(mask$x)))
  
  return(list(predicted_densities = predicted_densities, 
              detectors_df = detectors_df,
              estimated_sigma = estimated_sigma,
              capture_history = capture_history_max_occasions,
              mod = cfit))
  
}