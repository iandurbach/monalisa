library(dplyr)
library(ggplot2)
library(imager)
library(secr)
library(stringr)

source("run_secr.R")

# user sets parameters
sigma <- 4
lambda0 <- 0.1

# 1) Make true expected AC density surface from the Mona Lisa image

## load hi-res image of part of mona lisa
mona <- load.image("hires_mona.jpg")

## make a 200x200 version, which is the true intensity surface used to generate ACs
big_mona <- resize(mona, 200, 200, interpolation_type = 2)
big_mona_df <- as.data.frame(big_mona) %>% 
  mutate(y = 51 - ((y - 1)*(49.75/199) + 0.625), x = (x - 1)*(49.75/199) + 0.625)

## make a 50x50 version, which is the intensity surface after splitting space into 1x1 unit cells
## note: this is the resolution used by the secr mask
small_mona <- resize(mona, 50, 50, interpolation_type = 2)
small_mona_df <- as.data.frame(small_mona) %>% 
  mutate(y = 51 - y)

## checks
summary(big_mona_df)
summary(small_mona_df)
small_mona_df %>% 
  ggplot(aes(x, y)) + geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = hcl.colors(16, palette = "Blues")) + coord_equal() 
big_mona_df %>% 
  ggplot(aes(x, y)) + geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = hcl.colors(16, palette = "Blues")) + coord_equal() 

# 2) Simulate activity centres

# secr assumes units are in meters
# true intensity surface has cell dim 0.25x0.25 (=1/16m2 area) and secr wants D in ACs / 10000m2
# number you want to end up with in each cell is D = expected number of ACs / 10000m2
# number you have is the pixel intensity P
# first turn that into Ej = expected number of ACs in that cell (i.e. per 1/16m2); this should sum to EN across all cells,
# where EN is the desired number of ACs in expectation (actual number is random)
# Ej = P / sum(P) * EN
# then turn that into secr D = expected number of ACs / 10000m2
# D = Ej * 16 * 10000
# putting together, D = P / sum(P) * EN * 16 * 10000

big_mona_mask <- read.mask(data = big_mona_df, spacing = 0.25)
mona_mask <- read.mask(data = small_mona_df, spacing = 1)

big_mona_df$D = big_mona_df$value/sum(big_mona_df$value) * (16*80*10000) 
small_mona_df$D = small_mona_df$value/sum(small_mona_df$value) * (1*80*10000) 

simulated_points <- sim.popn(D = big_mona_df$D, 
                             core = big_mona_mask, 
                             model2D = "IHP", 
                             nDist = "poisson",
                             seed = 112)

# check
big_mona_df %>% 
  ggplot(aes(x, y)) + geom_raster(aes(fill = value)) +
  geom_point(data = simulated_points, colour = "red") +
  scale_fill_gradientn(colours = hcl.colors(16, palette = "Blues")) + coord_equal() 

# 3) Make covariate

## a blurred version of the true surface
blurry_mona <- isoblur(big_mona,16) 

## resize back to 50x50, type 2 for averaging
small_blurry_mona <- resize(blurry_mona, 50, 50, interpolation_type = 2)

## small adjustments to x and y for plotting
blurry_mona_df <- as.data.frame(blurry_mona) %>% 
  mutate(y = 51 - ((y - 1)*(49.75/199) + 0.625), x = (x - 1)*(49.75/199) + 0.625) %>%
  mutate(Dblur = value/sum(value) * (16*80*10000))

small_blurry_mona_df <- as.data.frame(small_blurry_mona) %>% 
  mutate(y = 51 - y) %>%
  mutate(Dblur = value/sum(value) * (1*80*10000))

covariates(mona_mask)$Dblur <- small_blurry_mona_df$Dblur

## checks
blurry_mona_df %>% 
  ggplot(aes(x, y)) + geom_raster(aes(fill = Dblur)) +
  scale_fill_gradientn(colours = hcl.colors(16, palette = "Blues")) + coord_equal() 

small_blurry_mona_df %>% 
  ggplot(aes(x, y)) + geom_raster(aes(fill = Dblur)) +
  scale_fill_gradientn(colours = hcl.colors(16, palette = "Blues")) + coord_equal() 

# 4) Simulate capture histories

## make a 3x3 grid of detectors
detectors3x3 <- make.grid(nx = 3, ny = 3, spacex = 8, spacey = 8,
                          originxy = c(10, 16), detector = "count")

## make a 7x7 grid of detectors
detectors7x7 <- make.grid(nx = 7, ny = 7, spacex = 8, spacey = 8,
                          originxy = c(1.5, 1.5), detector = "count")

## make capthist for 3x3 grid

### make a capture history for some maximum number of occasions
### note: more than needed, but will subset
max_nocc <- 555
my.seed <- 1234
capture_history_max_occasions <- sim.capthist(detectors3x3, popn = simulated_points, detectfn = "HHN",
                                              detectpar = list(lambda0 = lambda0, sigma = sigma),
                                              noccasions = max_nocc,
                                              nsessions = 1,
                                              seed = my.seed)

### extract capture histories after a few specified numbers of occasions 
### note: these will be used as input for secr.fit
occs3x3 <- c(18, 52, 111)
cht <- list()
j <- 0
for(i in occs3x3){
  j <- j + 1
  cht[[j]] <- subset(capture_history_max_occasions, occasions = 1:i, dropunused = FALSE)
}

### data frame containing capthists
cht2 <- tibble(expand.grid(i=1:1, noccasions = occs3x3))
cht2 <- cht2 %>% mutate(capthist = cht)
rm(cht)

### merge with other parameters
capthists_few_alloccs_3x3 <- tibble(expand.grid(i=1,
                                            xorig = 10, yorig = 16,
                                            noccasions = occs3x3))
capthists_few_alloccs_3x3 <- left_join(capthists_few_alloccs_3x3, cht2, by = c("i", "noccasions"))
rm(cht2)

## make capthist for 7x7 grid

### make a capture history for some maximum number of occasions
### note: more than needed, but will subset
max_nocc <- 1201
my.seed = 123
capture_history_max_occasions <- sim.capthist(detectors7x7, popn = simulated_points, detectfn = "HHN",
                                              detectpar = list(lambda0 = lambda0, sigma = sigma),
                                              noccasions = max_nocc,
                                              nsessions = 1,
                                              seed = my.seed)

### extract capture histories after a few specified numbers of occasions 
### note: these will be used as input for secr.fit
occs7x7 <- c(7, 25, 55)
cht <- list()
j <- 0
for(i in occs7x7){
  j <- j + 1
  cht[[j]] <- subset(capture_history_max_occasions, occasions = 1:i, dropunused = FALSE)
}

### data frame containing capthists
cht2 <- tibble(expand.grid(i=1:1, noccasions = occs7x7))
cht2 <- cht2 %>% mutate(capthist = cht)
rm(cht)

### merge with other parameters
capthists_few_alloccs_7x7 <- tibble(expand.grid(i=1,
                                            xorig = 1.5, yorig = 1.5,
                                            noccasions = occs7x7))
capthists_few_alloccs_7x7 <- left_join(capthists_few_alloccs_7x7, cht2, by = c("i", "noccasions"))
rm(cht2)

## save
save(big_mona_df, small_mona_df, blurry_mona_df, small_blurry_mona_df,
     simulated_points, 
     big_mona_mask, mona_mask,
     capthists_few_alloccs_7x7, capthists_few_alloccs_3x3, 
     file = "output/mona-inputs.RData")

# 5) fit models

## set up data frames with all parameters, so can fit models by looping over rows

### runs for 3x3 array
occs3x3 <- capthists_few_alloccs_3x3$noccasions
parlist1 <- expand.grid(i=1:1,
                        secr.fitformula = c("D~1", "D~log(Dblur)"),
                        dx = 8, dy = 8, nx = 3, ny = 3, 
                        xorig = 10, yorig = 16, 
                        sigma = sigma, lambda0 = lambda0, 
                        noccasions = occs3x3,
                        stringsAsFactors = FALSE)

### merge with capture histories
parlist1 <- left_join(parlist1, capthists_few_alloccs_3x3,
                      by = c("i", "noccasions", "xorig", "yorig"))

### runs for 7x7 array
occs7x7 <- capthists_few_alloccs_7x7$noccasions
parlist2 <- expand.grid(i=1:1,
                        secr.fitformula = c("D~1","D~log(Dblur)"),
                        dx = 8, dy = 8, nx = 7, ny = 7, 
                        xorig = 1.5, yorig = 1.5, 
                        sigma = 4, lambda0 = 0.1, 
                        noccasions = occs7x7,
                        stringsAsFactors = FALSE)

### merge with capture histories
parlist2 <- left_join(parlist2, capthists_few_alloccs_7x7,
                      by = c("i", "noccasions", "xorig", "yorig"))

### combine
parlist3 <- rbind(parlist1, parlist2)
rm(parlist1, parlist2, capthists_few_alloccs_3x3, capthists_few_alloccs_7x7)

## fit models 
set.seed(123)
res_acd <- list()
for(i in 1:nrow(parlist3)){
  print(i)
  res_acd[[i]] <- run_secr(simulated_points = simulated_points,
                               mask = mona_mask,
                               secr.fitformula = parlist3$secr.fitformula[i], 
                               dx = parlist3$dx[i], dy = parlist3$dy[i], 
                               nx = parlist3$nx[i], ny = parlist3$ny[i], 
                               xorig = parlist3$xorig[i], yorig = parlist3$yorig[i], 
                               sigma = parlist3$sigma[i], lambda0 = parlist3$lambda0[i], 
                               noccasions = parlist3$noccasions[i], 
                               capthist = parlist3$capthist[[i]])
}

## save
save(res_acd, file = "output/mona-results.RData")
