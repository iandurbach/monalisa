## Code to create uncertainty plots for the EACD plots shown in Figures 4 and 5

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
library(viridis)
## Objects we need
load("../output/mona-inputs.RData")
load("../output/mona-results.RData")
## Functions we need
source("Functions.R")

#####################################################################################################

## For each EACD plot in Figures 4 and 5, we will create two uncertainty plots: one showing the lower 5% quantile for the posterior distribution of the density for each pixel, the other will show the upper 95% quantile.
## We will create two figures. One will show the uncertainty plots from Figure 4, along with the corresponding EACD plots. The other figure will show the same thing, but for Figure 5. We will refer to these figures as 'uncertainty figures' in the code below.

## First, we will create the objects needed to create both figures. We will then put together both figures at the end. This is because we want the same maximum value to be used when colouring both figures.

## ---------------------------------------------------------------------------------------
#################### Objects needed for uncertainty figure for Figure 4 ##################
## ---------------------------------------------------------------------------------------

## Loading in the MCMC results from Row 2 of Figure 4
load("MCMC_Results/Figure4/InhomPP_18occ.RData")
load("MCMC_Results/Figure4/InhomPP_52occ.RData")
load("MCMC_Results/Figure4/InhomPP_111occ.RData")
# Discarding burn-in
inhom.results.18occ <- inhom.results.18occ[-c(1:1000),]
inhom.results.52occ <- inhom.results.52occ[-c(1:1000),]
inhom.results.111occ <- inhom.results.111occ[-c(1:1000),]

## Covariate that we use for these EACD plots (see 'Functions.R' for an explanation of the function)
log.dblur <- eacd.covariate()

## Creating the data frames that will contain all of the info required to create the uncertainty figure for Figure 4
# 18 sampling occasions
info.18occ <- eacd.quantile.info(results=inhom.results.18occ, covariate=log.dblur, nPix=2500, nocc=18)
# 52 sampling occasions
info.52occ <- eacd.quantile.info(results=inhom.results.52occ, covariate=log.dblur, nPix=2500, nocc=52)
# 111 sampling occasions
info.111occ <- eacd.quantile.info(results=inhom.results.111occ, covariate=log.dblur, nPix=2500, nocc=111)

## ---------------------------------------------------------------------------------------
#################### Objects needed for uncertainty figure for Figure 5 ##################
## ---------------------------------------------------------------------------------------

## Loading in the MCMC results from Row 2 of Figure 5
load("MCMC_Results/Figure5/InhomPP_7occ.RData")
load("MCMC_Results/Figure5/InhomPP_25occ.RData")
load("MCMC_Results/Figure5/InhomPP_55occ.RData")
# Discarding burn-in
inhom.results.7occ <- inhom.results.7occ[-c(1:1000),]
inhom.results.25occ <- inhom.results.25occ[-c(1:1000),]
inhom.results.55occ <- inhom.results.55occ[-c(1:1000),]

## Covariate we use is the same as above
log.dblur <- eacd.covariate()

## Creating the data frames that will contain all of the info require to create the uncertainty figure for Figure 5
# 7 sampling occasions
info.7occ <- eacd.quantile.info(results=inhom.results.7occ, covariate=log.dblur, nPix=2500, nocc=7)
# 25 sampling occasions
info.25occ <- eacd.quantile.info(results=inhom.results.25occ, covariate=log.dblur, nPix=2500, nocc=25)
# 55 sampling occasions
info.55occ <- eacd.quantile.info(results=inhom.results.55occ, covariate=log.dblur, nPix=2500, nocc=55)

## ---------------------------------------------------------------------------------------
####################### Objects needed for both uncertainty figures ######################
## ---------------------------------------------------------------------------------------

## Collating all information for figures into one data frame
values_all <- rbind(info.18occ, info.52occ, info.111occ, info.7occ, info.25occ, info.55occ)

## Detectors
detectors_df_all <- res_acd %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()

# Saving the objects we have created, for us to use in appendix.Rnw -- if have already saved 'detectors.RData' from 'Plots_Code.R', then don't need to save 'detectors.RData' here, as both objects are the same!
#save(values_all, file="uncertainty_plots.RData")
#save(detectors_df_all, file="detectors.RData")

## Maximum value for colour scale for *all* plots
maxval <- max(values_all$value)

## ---------------------------------------------------------------------------------------
######################## Creating uncertainty figure for Figure 4 ########################
## ---------------------------------------------------------------------------------------

nn <- 3 # Number of different simulated datasets used in Figure 4
occ <- capthists_few_alloccs_3x3$noccasions # Number of sampling occasions for each dataset
asz <- c("3x3")
chs <- data.frame(do.call(rbind, lapply(capthists_few_alloccs_3x3$capthist, summary, terse = TRUE)))
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = chs$Detections, .y = chs$Animals, .f = paster) %>% unlist() # Row labels we'll use

## Adding faceting info to 'values_all' for rows that we will use in this uncertainty figure
values_all$covtype2 <- factor(values_all$covtype, levels=unique(values_all$covtype),
                              labels=c("Expected AC"))
values_all$quantile <- factor(rep(c(rep("0.05 quantile", 2601), rep("Mean", 2601), rep("0.95 quantile", 2601)), 3),
                              levels=c("0.05 quantile", "Mean", "0.95 quantile")) # (We use 2601 as the data frames we work with to create each individual plot contain 2601 entries)
values_all$occasions2 <- factor(values_all$occasions,
                                levels = occ,
                                labels = capthist_labels)

## Creating and saving the uncertainty figure
uncertainty.fig4 <- values_all %>%
  filter(occasions %in% occ[1:nn], array_size %in% asz) %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  scale_fill_distiller(limits=c(0, maxval)) +
  facet_grid(occasions2 ~ quantile) +
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

ggsave("mona_3x3_uncertainty.png", uncertainty.fig4, width=8, height=8, dpi=600, bg="white")

## ---------------------------------------------------------------------------------------
######################## Creating uncertainty figure for Figure 5 ########################
## ---------------------------------------------------------------------------------------

nn <- 3
occ <-capthists_few_alloccs_7x7$noccasions
asz <- c("7x7")
chs <- data.frame(do.call(rbind, lapply(capthists_few_alloccs_7x7$capthist, summary, terse = TRUE)))
chs <- chs %>% dplyr::filter(Occasions %in% occ)
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = chs$Detections, .y = chs$Animals, .f = paster) %>% unlist()

## Adding faceting info for rows that we will use to create this uncertainty figure
values_all$occasions2 <- factor(values_all$occasions,
                                levels = occ,
                                labels = capthist_labels)

## Creating and saving the uncertainty figure
uncertainty.fig5 <- values_all %>%
  filter(occasions %in% occ[1:nn], array_size %in% asz) %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  scale_fill_distiller(limits=c(0, maxval)) +
  facet_grid(occasions2 ~ quantile) +
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

ggsave("mona_7x7_uncertainty.png", uncertainty.fig5, width=8, height=8, dpi=600, bg="white")

#####################################################################################################

## Creating plots to illustrate that when we plot the CV for RACD plots, the result depends on pixel size. We will create 3 plots: one with 1x1 pixels (pixel size used for all maps in paper), one with 2x2 pixels (four pixels aggregated to be one pixel), one with 5x5 pixels (25 pixels aggregated to be one pixel)

## Data we'll be using
load("MCMC_Results/Figure4/HomPP_18occ.RData")

## Calculating the CV values
cv.1by1 <- uncertainty.values.racd(xlim=c(0.5, 50.5), ylim=c(0.5, 50.5), results=results.18occ, M=300, pixel.index=1, cv=TRUE)
cv.2by2 <- uncertainty.values.racd(xlim=c(0.5, 50.5), ylim=c(0.5, 50.5), results=results.18occ, M=300, pixel.index=2, cv=TRUE)
cv.5by5 <- uncertainty.values.racd(xlim=c(0.5, 50.5), ylim=c(0.5, 50.5), results=results.18occ, M=300, pixel.index=5, cv=TRUE)

## Creating a data frame containing all of the information we will use
# 1x1 pixels
pixels.1by1 <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50) - 0.5 # Editing pixel info so we are storing pixel edges instead of pixel centers. We want to colour pixels using pixel edges rather than pixel centres, so our maps will be filled in correctly.
cvinfo.1by1 <- data.frame(x=pixels.1by1[,1], y=pixels.1by1[,2], covtype=rep("D~1", 2500), occasions=rep(18, 2500), array_size=rep("3x3", 2500), value=cv.1by1, pixel.size=rep("1x1 pixels", 2500))
# 2x2 pixels
pixels.2by2 <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=25, y.pixels=25) - 1 # Storing pixel edges
cvinfo.2by2 <- data.frame(x=pixels.2by2[,1], y=pixels.2by2[,2], covtype=rep("D~1", 625), occasions=rep(18, 625), array_size=rep("3x3", 625), value=cv.2by2, pixel.size=rep("2x2 pixels", 625))
# 5x5 pixels
pixels.5by5 <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=10, y.pixels=10) - 2.5 # Storing pixel edges
cvinfo.5by5 <- data.frame(x=pixels.5by5[,1], y=pixels.5by5[,2], covtype=rep("D~1", 100), occasions=rep(18, 100), array_size=rep("3x3", 100), value=cv.5by5, pixel.size=rep("5x5 pixels", 100))
# Overall data frame
cv_values_all = rbind(cvinfo.1by1, cvinfo.2by2, cvinfo.5by5)
## Extracting and editing entries from this data frame, so it contains pixel edges for our whole map area (right now, are missing pixel edges for rightmost and topmost edge), along with their corresponding CV value. Doing this means that we will colour the pixels in our map area correctly.
# 1x1 pixels
dup1 <- cv_values_all[(cv_values_all$pixel.size=="1x1 pixels" & cv_values_all$y==49.5 | cv_values_all$pixel.size=="1x1 pixels" & cv_values_all$x==49.5),]
save1 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save1$x=50.5
save2 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save2$y=50.5 # If we don't run these two lines, then we'll miss these two sets of pixel edges in our data frame
dup1$x[dup1$x==49.5] = 50.5; dup1$y[dup1$y==49.5] = 50.5 # Editing all of the entries in dup1, so that they represent pixel edges along the right and top edges of the map area
dup1 <- rbind(dup1, save1, save2)
# 2x2 pixels
dup2 <- cv_values_all[(cv_values_all$pixel.size=="2x2 pixels" & cv_values_all$y==48.5 | cv_values_all$pixel.size=="2x2 pixels" & cv_values_all$x==48.5),]
save1 <- dup2[(dup2$x==48.5 & dup2$y==48.5),]; save1$x=50.5
save2 <- dup2[(dup2$x==48.5 & dup2$y==48.5),]; save2$y=50.5
dup2$x[dup2$x==48.5] = 50.5; dup2$y[dup2$y==48.5] = 50.5
dup2 <- rbind(dup2, save1, save2)
# 5x5 pixels
dup3 <- cv_values_all[(cv_values_all$pixel.size=="5x5 pixels" & cv_values_all$y==45.5 | cv_values_all$pixel.size=="5x5 pixels" & cv_values_all$x==45.5),]
save1 <- dup3[(dup3$x==45.5 & dup3$y==45.5),]; save1$x=50.5
save2 <- dup3[nrow(dup3),]; save2$y=50.5
dup3$x[dup3$x==45.5] = 50.5; dup3$y[dup3$y==45.5] = 50.5
dup3 <- rbind(dup3, save1, save2)
# Adding this extra information to 'cv_values_all'
cv_values_all  <- rbind(cv_values_all, dup1, dup2, dup3)
cv_values_all$valtype <- factor(rep("Realised AC CV", nrow(cv_values_all))) # Adding column for faceting in final plot

## Detectors
detectors_df_all <- res_acd %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()
detectors_df_all <- detectors_df_all %>% filter(covtype %in% c("D~1"), array_size %in% c("3x3"))

## Max CV value across all 3 plots
cvmaxval <- max(cv_values_all$value)

## Creating the figure!
nn <- 3
asz <- c("3x3")
occ <- 18
racd.cv.fig  <-  cv_values_all %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  scale_fill_distiller(limits=c(0,cvmaxval)) +
  facet_grid(valtype ~ pixel.size) +
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
        legend.position="right", legend.key.width = unit(0.56, "cm"),
        legend.key.height = unit(0.9,"cm"),
        panel.background=element_blank(), panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

racd.cv.fig

# Saving the objects that we used, for use in appendix.Rnw
#save(cv_values_all, file="cv_plots.RData")
#save(detectors_df_all, file="detectors_cv_plots.RData") # Will also use these detectors for the plots below :)

#ggsave(filename="cvfigure.png", plot=racd.cv.fig, bg="white")

# ------------------------------------------------------------

## Adding plots that show surfaces showing the posterior standard deviation also vary depending on pixel size:

load("MCMC_Results/Figure4/HomPP_18occ.RData") # Just in case anything changed with this object (v unlikely)

## Calculating the CV values
sd.1by1 <- uncertainty.values.racd(xlim=c(0.5, 50.5), ylim=c(0.5, 50.5), results=results.18occ, M=300, pixel.index=1, cv=FALSE)
sd.2by2 <- uncertainty.values.racd(xlim=c(0.5, 50.5), ylim=c(0.5, 50.5), results=results.18occ, M=300, pixel.index=2, cv=FALSE)
sd.5by5 <- uncertainty.values.racd(xlim=c(0.5, 50.5), ylim=c(0.5, 50.5), results=results.18occ, M=300, pixel.index=5, cv=FALSE)

## Creating a data frame containing all of the information we will use
# 1x1 pixels
pixels.1by1 <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50) - 0.5 # Editing pixel info so we are storing pixel edges instead of pixel centers
sdinfo.1by1 <- data.frame(x=pixels.1by1[,1], y=pixels.1by1[,2], covtype=rep("D~1", 2500), occasions=rep(18, 2500), array_size=rep("3x3", 2500), value=sd.1by1, pixel.size=rep("1x1 pixels", 2500))
# 2x2 pixels
pixels.2by2 <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=25, y.pixels=25) - 1 # Storing pixel edges
sdinfo.2by2 <- data.frame(x=pixels.2by2[,1], y=pixels.2by2[,2], covtype=rep("D~1", 625), occasions=rep(18, 625), array_size=rep("3x3", 625), value=sd.2by2, pixel.size=rep("2x2 pixels", 625))
# 5x5 pixels
pixels.5by5 <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=10, y.pixels=10) - 2.5 # Storing pixel edges
sdinfo.5by5 <- data.frame(x=pixels.5by5[,1], y=pixels.5by5[,2], covtype=rep("D~1", 100), occasions=rep(18, 100), array_size=rep("3x3", 100), value=sd.5by5, pixel.size=rep("5x5 pixels", 100))
# Overall data frame
sd_values_all = rbind(sdinfo.1by1, sdinfo.2by2, sdinfo.5by5)
## Extracting and editing entries from this data frame as we did above. Recall, we want to colour our map using pixel edges rather than pixel centres, so the pixels across the map will be filled in correctly. Our manipulation below allows for this.
dup1 <- sd_values_all[(sd_values_all$pixel.size=="1x1 pixels" & sd_values_all$y==49.5 | sd_values_all$pixel.size=="1x1 pixels" & sd_values_all$x==49.5),]
save1 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save1$x=50.5
save2 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save2$y=50.5
dup1$x[dup1$x==49.5] = 50.5; dup1$y[dup1$y==49.5] = 50.5
dup1 <- rbind(dup1, save1, save2)
# 2x2 pixels
dup2 <- sd_values_all[(sd_values_all$pixel.size=="2x2 pixels" & sd_values_all$y==48.5 | sd_values_all$pixel.size=="2x2 pixels" & sd_values_all$x==48.5),]
save1 <- dup2[(dup2$x==48.5 & dup2$y==48.5),]; save1$x=50.5
save2 <- dup2[(dup2$x==48.5 & dup2$y==48.5),]; save2$y=50.5
dup2$x[dup2$x==48.5] = 50.5; dup2$y[dup2$y==48.5] = 50.5
dup2 <- rbind(dup2, save1, save2)
# 5x5 pixels
dup3 <- sd_values_all[(sd_values_all$pixel.size=="5x5 pixels" & sd_values_all$y==45.5 | sd_values_all$pixel.size=="5x5 pixels" & sd_values_all$x==45.5),]
save1 <- dup3[(dup3$x==45.5 & dup3$y==45.5),]; save1$x=50.5
save2 <- dup3[nrow(dup3),]; save2$y=50.5
dup3$x[dup3$x==45.5] = 50.5; dup3$y[dup3$y==45.5] = 50.5
dup3 <- rbind(dup3, save1, save2)
# Adding this extra information to 'cv_values_all'
sd_values_all  <- rbind(sd_values_all, dup1, dup2, dup3)
sd_values_all$valtype <- factor(rep("Realised AC SD", nrow(sd_values_all))) # Adding column for faceting in final plot

## Saving object to use in Appendix.Rnw
#save(sd_values_all, file='sd_plots.RData')

## Max SD value across all 3 plots
sdmaxval <- max(sd_values_all$value)

## Creating the figure!
nn <- 3
asz <- c("3x3")
occ <- 18
racd.sd.fig  <-  sd_values_all %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  scale_fill_distiller(limits=c(0,sdmaxval), breaks=c(0, 0.1, 0.2), labels=c("0", "0.100", "0.200")) +
  facet_grid(valtype ~ pixel.size) +
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
        legend.position="right", legend.key.width = unit(0.56, "cm"),
        legend.key.height = unit(0.9,"cm"), legend.title = element_blank(),
        panel.background=element_blank(), panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

racd.sd.fig

## Combining CV and SD plots to build up figure

racd.cv.fig/racd.sd.fig

# ------------------------------------------------------------

## Want to add corresponding RACD plot onto left, generating this RACD plot and attaching it to create the final figure:

load("MCMC_Results/Figure4/HomPP_18occ.RData") # Just to be sure it is loaded as required

## Vector we need to create the map:
racd.18occ <- racd.density.vector(results=results.18occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))

## Summarising plot info. Note that once again, are manipulating data frame so it contains pixel edges rather than pixel centres, so the pixels in the map are coloured correctly
pixel.centres <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50)
pixel.edges <- pixel.centres - 0.5
dat <- data.frame(x=pixel.edges[,1], y=pixel.edges[,2], covtype=rep("D~1", 2500), occasions=rep(18, 2500), array_size=rep("3x3", 2500), value=racd.18occ)
dup1 <- dat[(dat$y==49.5 | dat$x==49.5),]
save1 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save1$x=50.5
save2 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save2$y=50.5
dup1$x[dup1$x==49.5] = 50.5; dup1$y[dup1$y==49.5] = 50.5
dup1 <- rbind(dup1, save1, save2)
dat <- rbind(dat, dup1)
predicted_densities_all <- dat
# Faceting variables we will use
chs <- data.frame(do.call(rbind, lapply(capthists_few_alloccs_3x3$capthist, summary, terse = TRUE)))
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = chs$Detections, .y = chs$Animals, .f = paster) %>% unlist()
predicted_densities_all$occasions2 <- rep(capthist_labels[1], 2601)
predicted_densities_all$covtype2 <- rep("Realised AC", 2601)

## Saving values to use in Appendix.Rnw
#save(predicted_densities_all, file="racd_18occ.RData")

## Max val for plot
maxval <- 0.2301 # Is 'maxval' in 'Plots_Code.R'. Are using it again here as want to use the same colour scale as in Figures 4 and 5

## Creating the RACD plot
nn <- 3
occ <- 18
asz <- c("3x3")
fig4.18occ <- predicted_densities_all %>%
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

fig4.18occ

## Okay, now trying to put all three pieces together to create the final figure :)
# One possible configuration
ggarrange(fig4.18occ, ggarrange(racd.cv.fig, racd.sd.fig, nrow=2), nrow=2, heights=c(1, 2))
# Another
#ggarrange(fig4.18occ, ggarrange(racd.cv.fig, racd.sd.fig, nrow=2), ncol=2, widths=c(1, 2))
ggsave("SD_CV_uncertainty.png", ggarrange(fig4.18occ, ggarrange(racd.cv.fig, racd.sd.fig, nrow=2), nrow=2, heights=c(1, 2)), width=8, height=10, bg="white")
