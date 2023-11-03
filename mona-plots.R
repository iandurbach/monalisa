library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)
library(secr)
library(patchwork)
library(ggpubr)

load("output/revision/mona-inputs.RData")
load("output/revision/mona-results.RData")

# # if not plotting all capthists
# capthists_few_alloccs_3x3 <- c(capthists_few_alloccs_3x3[c(1:3,22),])
# capthists_few_alloccs_7x7 <- c(capthists_few_alloccs_7x7[c(1:3,24),])

# process the outputs
predicted_densities_all <- res_acd %>% purrr::map_depth(1, "predicted_densities") %>% map_df(bind_rows)

# dealing with pixel edges instead of pixel centres, so pixels are coloured correctly
predicted_densities_all$x <- predicted_densities_all$x - 0.5
predicted_densities_all$y <- predicted_densities_all$y - 0.5
dup1 <- predicted_densities_all[(predicted_densities_all$y==49.5 | predicted_densities_all$x==49.5),]
save1 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save1$x=50.5
save2 <- dup1[(dup1$x==49.5 & dup1$y==49.5),]; save2$y=50.5 # if we don't run these two lines, then we'll miss these two sets of pixel edges in our data frame
dup1$x[dup1$x==49.5] = 50.5; dup1$y[dup1$y==49.5] = 50.5 # editing all of the entries in dup1, so that they represent pixel edges along the right and top edges of the map area
dup1 <- rbind(dup1, save1, save2)
predicted_densities_all <- rbind(predicted_densities_all, dup1)

# predicted densities per cell (1 cell = "1m2" = 1/10000ha)
predicted_densities_all$value <- predicted_densities_all$prob_ac/ 10000

# standardised densities (not used)
predicted_densities_all <- predicted_densities_all %>% group_by(array_size, occasions) %>%
  mutate(stdvalue = (value - min(value)) / (max(value) - min(value)))

# check that its worked
predicted_densities_all %>% group_by(occasions, array_size, covtype) %>%
  summarize(mp = mean(prob_ac),
            me = mean(expnumber_ac),
            mv = mean(value),
            sumv = sum(value),
            vv = var(value))

detectors_df_all <- res_acd %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()

nn <- 3
occ <- capthists_few_alloccs_3x3$noccasions
asz <- c("3x3")

chs <- data.frame(do.call(rbind, lapply(capthists_few_alloccs_3x3$capthist, summary, terse = TRUE)))
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = chs$Detections, .y = chs$Animals, .f = paster) %>% unlist()

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

# use same top of colour scale for all plots
xx <- predicted_densities_all %>% filter(array_size == "3x3", occasions %in% capthists_few_alloccs_3x3$noccasions[1:nn])
maxval1 <- max(xx$value)
xx <- predicted_densities_all %>% filter(array_size == "7x7", occasions %in% capthists_few_alloccs_7x7$noccasions[1:nn])
maxval2 <- max(xx$value)
maxval <- max(maxval1, maxval2)

### Plot of input data

# plots of inputs

i1 <- big_mona_df %>%
  ggplot(aes(x, y)) + geom_raster(aes(fill = D/10000)) +
  #scale_fill_gradientn(colours = hcl.colors(16, palette = "Light grays"), limits = c(0,maxval)) + coord_equal() +
  scale_fill_distiller(limits = c(0,maxval)) + coord_equal() +
  theme_classic(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none", legend.key.width = unit(1, "cm"),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

i2 <- small_mona_df %>%
  ggplot(aes(x, y)) + geom_raster(aes(fill = D/10000)) +
  #scale_fill_gradientn(colours = hcl.colors(16, palette = "Light grays"), limits = c(0,maxval)) + coord_equal() +
  scale_fill_distiller(limits = c(0,maxval)) + coord_equal() +
  theme_classic(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none", legend.key.width = unit(1, "cm"),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

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

i4 <- small_blurry_mona_df %>%
  ggplot(aes(x, y)) + geom_raster(aes(fill = Dblur/10000)) +
  #scale_fill_gradientn(colours = hcl.colors(16, palette = "Light grays"), limits = c(0,maxval)) + coord_equal() +
  scale_fill_distiller(limits = c(0,maxval)) + coord_equal() +
  theme_classic(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none", legend.key.width = unit(1, "cm"),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

i5 <- ggarrange(i1,i3,i4,labels = c("(a)", "(b)", "(c)"), ncol = 3, nrow = 1,
                label.x = 0.5, label.y = 0, vjust = -0.2, hjust = 0.5,
                font.label = list(size = 14, color = "black", face = "plain", family = NULL))

i5

ggsave("paper/mona_inputdata.png", i5, width=7.5, height=2.5, dpi = 600)

### Density surfaces with 3x3 grid

p2a <- predicted_densities_all %>%
  filter(occasions %in% occ[1:nn], array_size %in% asz) %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  #scale_fill_gradientn(colours = hcl.colors(16, palette = "Light grays"), limits = c(0,maxval)) +
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

### Density surfaces with 7x7 grid

occ <-capthists_few_alloccs_7x7$noccasions
asz <- c("7x7")

chs <- data.frame(do.call(rbind, lapply(capthists_few_alloccs_7x7$capthist, summary, terse = TRUE)))
chs <- chs %>% dplyr::filter(Occasions %in% occ)
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = chs$Detections, .y = chs$Animals, .f = paster) %>% unlist()

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

p2b <- predicted_densities_all %>%
  filter(occasions %in% occ[1:nn], array_size %in% asz) %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  #scale_fill_gradientn(colours = hcl.colors(16, palette = "Light grays"), limits = c(0,maxval)) +
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

ggsave("output/mona_3x3.png", p2a, width=8, height=6, dpi = 600)
ggsave("output/mona_7x7.png", p2b, width=8, height=6, dpi = 600)


