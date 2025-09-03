
library(phyloseq)
rarefied_physeq <- readRDS("rarefied_physeq.rds")
coral<-subset_samples(rarefied_physeq, Treatment =="Sclero")
coral<-subset_samples(rarefied_physeq, Treatment =="Lobo")

wild_ids  <- c("wild_0","wild_6","wild_19")
tank_ids  <- c(paste0("tank_",1:12),"tank__10","tank_19")
types     <- c(wild_ids, tank_ids)
ps_sel    <- subset_samples(coral, All_type %in% types)
bc_dist   <- phyloseq::distance(ps_sel, method = "bray",binary=F)
bc_mat    <- as.matrix(bc_dist)

wild_samps <- sample_names(ps_sel)[ sample_data(ps_sel)$All_type %in% wild_ids ]
tank_samps <- sample_names(ps_sel)[ sample_data(ps_sel)$All_type %in% tank_ids ]

dist_mat   <- bc_mat[wild_samps, tank_samps, drop = FALSE]


library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
dist_df <- as.data.frame(dist_mat) %>%
  rownames_to_column(var = "WildRep") %>%
  pivot_longer(
    cols      = -WildRep,
    names_to  = "TankRep",
    values_to = "BrayDist"
  )

meta_df <- sample_data(coral) %>%
  as("data.frame") %>%
  rownames_to_column(var = "Sample")   

dist_df$TankType <- meta_df$Time[ match(dist_df$TankRep, meta_df$Sample) ]

tank_levels <- c(paste0(1:12),19)
dist_df$TankType <- factor(dist_df$TankType, levels = tank_levels)

library(ggplot2)




library(ggplot2)
library(ggplot2)
library(scales) 
my_blue <-"#4b7c8a"
ggplot(dist_df, aes(x = factor(TankType), y = BrayDist, fill = factor(TankType))) +
  geom_boxplot(
    fill          = alpha(my_blue, 0.3), 
    color         = my_blue,
    outlier.shape = 21,
    outlier.color = my_blue,
    width         = 0.8
  )  +
  geom_jitter(
    aes(color = factor(TankType)),
    size   = 0.8,
    width  = 0.15,
    alpha  = 0.3,
    show.legend = FALSE
  ) +
  stat_summary(
    fun        = median,
    geom       = "line",
    aes(group   = 1),
    color      = "#4b7c8a",
    size       = 0.8
  ) +
  stat_summary(
    fun        = median,
    geom       = "point",
    aes(group   = 1),
    shape      = 21,
    size       = 2.5,
    fill       = "#4b7c8a",
    color      = "#4b7c8a"
  ) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  scale_color_grey(start = 0.3, end = 0.5) +
  scale_x_discrete(
    name  = "Week",
    labels = c("1","2","3","4","5","6","7","8","9","10","11","12","19")
  ) +
  scale_y_continuous(
    name   = "Bray–Curtis Distance",
    breaks = seq(0, 1, by = 0.2)
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.border      = element_rect(color = "black", fill = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.text.x       = element_text(size = 12, angle = 0, hjust = 1),      
    axis.text.y       = element_text(size = 12),                            
    axis.title.x      = element_text(size = 16, face = "bold", vjust = -0.5), 
    axis.title.y      = element_text(size = 16, face = "bold", vjust =  1.5)  
  )+ theme(
    axis.title.x = element_blank()
  )+
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2) 
  )
