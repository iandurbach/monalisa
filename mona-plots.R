library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)
library(secr)
library(patchwork)
library(ggpubr)

load("output/mona-inputs.RData")
load("output/mona-results.RData")

# process the outputs
predicted_densities_all <- res_acd %>% purrr::map_depth(1, "predicted_densities") %>% map_df(bind_rows)
detectors_df_all <- res_acd %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()

# small adjustment to deal with pixel edges instead of pixel centres so pixels are coloured correctly
predicted_densities_all$x <- predicted_densities_all$x - 0.5
predicted_densities_all$y <- predicted_densities_all$y - 0.5
dup1 <- predicted_densities_all[(predicted_densities_all$y==49.5 | predicted_densities_all$x==49.5),]
save1 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save1$x=50.5
save2 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save2$y=50.5 # if we don't run these two lines, then we'll miss these two sets of pixel edges in our data frame
dup1$x[dup1$x==49.5] = 50.5; dup1$y[dup1$y==49.5] = 50.5 # editing all of the entries in dup1, so that they represent pixel edges along the right and top edges of the map area
dup1 <- rbind(dup1, save1, save2)
predicted_densities_all <- rbind(predicted_densities_all, dup1)

# predicted densities per cell (1 cell = "1m2" = 1/10000ha)
predicted_densities_all$value <- predicted_densities_all$prob_ac / 10000

# check
predicted_densities_all %>% group_by(occasions, array_size, covtype) %>%
  summarize(mp = mean(prob_ac),
            me = mean(expnumber_ac),
            mv = mean(value),
            sumv = sum(value),
            vv = var(value))

# number of occasions
nn <- 3

# use same top of colour scale for all plots
xx <- predicted_densities_all %>% filter(array_size == "3x3", occasions %in% capthists_few_alloccs_3x3$noccasions[1:nn])
maxval1 <- max(xx$value)
xx <- predicted_densities_all %>% filter(array_size == "7x7", occasions %in% capthists_few_alloccs_7x7$noccasions[1:nn])
maxval2 <- max(xx$value)
maxval <- max(maxval1, maxval2)

# Plot 1: plot of input data 
# --------------------------

## Plot 1a: true surface
i1 <- big_mona_df %>%
  ggplot(aes(x, y)) + geom_raster(aes(fill = D/10000)) +
  scale_fill_distiller(limits = c(0,maxval)) + coord_equal() +
  theme_classic(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none", legend.key.width = unit(1, "cm"),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

## Plot 1x: true surface at mask resolution (not included in final plot)
i2 <- small_mona_df %>%
  ggplot(aes(x, y)) + geom_raster(aes(fill = D/10000)) +
  scale_fill_distiller(limits = c(0,maxval)) + coord_equal() +
  theme_classic(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none", legend.key.width = unit(1, "cm"),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

## Plot 1b: true AC locations
small_mona_df$x <- small_mona_df$x - 0.5
small_mona_df$y <- small_mona_df$y - 0.5
dup1 <- small_mona_df[(small_mona_df$y==49.5 | small_mona_df$x==49.5),]
save1 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save1$x=50.5
save2 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save2$y=50.5 # if we don't run these two lines, then we'll miss these two sets of pixel edges in our data frame
dup1$x[dup1$x==49.5] = 50.5; dup1$y[dup1$y==49.5] = 50.5 # editing all of the entries in dup1, so that they represent pixel edges along the right and top edges of the map area
dup1 <- rbind(dup1, save1, save2)
small_mona_df <- rbind(small_mona_df, dup1)

i3 <- small_mona_df %>%
  ggplot(aes(x, y)) +
  geom_raster(fill = "gray80") +
  geom_point(data = simulated_points, inherit.aes = F, aes(x=x,y=y),
             colour = "red", pch = 16, size = 1, alpha = 0.5) + coord_equal() +
  theme_classic(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none", legend.key.width = unit(1, "cm"),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

## Plot 1c: covariate
i4 <- small_blurry_mona_df %>%
  ggplot(aes(x, y)) + geom_raster(aes(fill = Dblur/10000)) +
  scale_fill_distiller(limits = c(0,maxval)) + coord_equal() +
  theme_classic(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none", legend.key.width = unit(1, "cm"),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

## combine
i5 <- ggarrange(i1,i3,i4,labels = c("(a)", "(b)", "(c)"), ncol = 3, nrow = 1,
                label.x = 0.5, label.y = 0, vjust = -0.2, hjust = 0.5,
                font.label = list(size = 14, color = "black", face = "plain", family = NULL))

i5

# Plot 2: density surfaces with 3x3 grid 
# --------------------------------------

## plotting details 

occ <- capthists_few_alloccs_3x3$noccasions
asz <- c("3x3")

### make nice facet labels 
chs <- data.frame(do.call(rbind, lapply(capthists_few_alloccs_3x3$capthist, summary, terse = TRUE)))
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = chs$Detections, .y = chs$Animals, .f = paster) %>% unlist()

### rename occasions as factor with desired facet labels
predicted_densities_all$occasions2 <- factor(predicted_densities_all$occasions,
                                             levels = occ,
                                             labels = capthist_labels)

detectors_df_all$occasions2 <- factor(detectors_df_all$occasions,
                                      levels = occ,
                                      labels = capthist_labels)

predicted_densities_all$covtype2 <- factor(predicted_densities_all$covtype,
                                           levels = unique(predicted_densities_all$covtype),
                                           labels = c("Predicted AC", "Expected AC"))

### rename surface type as factor with desired facet labels
detectors_df_all$covtype2 <- factor(detectors_df_all$covtype,
                                    levels = unique(detectors_df_all$covtype),
                                    labels = c("Predicted AC", "Expected AC"))

## final plot
p2a <- predicted_densities_all %>%
  filter(occasions %in% occ[1:nn], array_size %in% asz) %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  scale_fill_distiller(limits = c(0,maxval)) +
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
        legend.position="none", legend.key.width = unit(3, "cm"),
        legend.key.height = unit(0.7,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p2a

# Plot 3: density surfaces with 7x7 grid 
# --------------------------------------

## plotting details

occ <-capthists_few_alloccs_7x7$noccasions
asz <- c("7x7")

### make nice facet labels 
chs <- data.frame(do.call(rbind, lapply(capthists_few_alloccs_7x7$capthist, summary, terse = TRUE)))
chs <- chs %>% dplyr::filter(Occasions %in% occ)
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = chs$Detections, .y = chs$Animals, .f = paster) %>% unlist()

### rename factors
predicted_densities_all$occasions2 <- factor(predicted_densities_all$occasions,
                                             levels = occ,
                                             labels = capthist_labels)

detectors_df_all$occasions2 <- factor(detectors_df_all$occasions,
                                      levels = occ,
                                      labels = capthist_labels)

predicted_densities_all$covtype2 <- factor(predicted_densities_all$covtype,
                                           levels = unique(predicted_densities_all$covtype),
                                           labels = c("Predicted AC", "Expected AC"))

detectors_df_all$covtype2 <- factor(detectors_df_all$covtype,
                                    levels = unique(detectors_df_all$covtype),
                                    labels = c("Predicted AC", "Expected AC"))

## final plot
p2b <- predicted_densities_all %>%
  filter(occasions %in% occ[1:nn], array_size %in% asz) %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  scale_fill_distiller(limits = c(0,maxval)) +
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
        legend.position="none", legend.key.height = unit(0.7,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p2b

# save
ggsave("paper/figures/mona_inputdata.jpg", i5, width=7.5, height=2.5, dpi = 600)
ggsave("paper/figures/mona_3x3.jpg", p2a, width=8, height=6, dpi = 600)
ggsave("paper/figures/mona_7x7.jpg", p2b, width=8, height=6, dpi = 600)


