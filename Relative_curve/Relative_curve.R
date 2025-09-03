library(phyloseq)
library(dplyr)
library(pheatmap)
library(rlang)
library(ggthemes)
library(phyloseq)
library(RColorBrewer)
library(ggplot2)
library(vegan)
library(ggpubr)
library(dplyr)
library(MicrobiomeStat)
#library(rpl)
library(Rarefy)
library(ggh4x)
library(randomcoloR)
library(forcats)
#library(MicEco)
library(eulerr)

rm(list=ls()) 
rarefied_physeq <- readRDS("/Users/linyanzhi/Desktop/EXP-珊瑚代謝體及微生物組成實驗/bacteria/rarefied_physeq.rds")

lobo_ps <- subset_samples(rarefied_physeq, Treatment == "Sclero")

physeq_genus   <- tax_glom(lobo_ps, taxrank = "Genus")
physeq_genus_rel <- transform_sample_counts(
  physeq_genus, 
  function(x) 100 * x / sum(x)
)


sig_genera <- read.csv("sinu_lefse.csv", stringsAsFactors=FALSE)$Genus

phy_sig <- subset_taxa(physeq_genus_rel, Genus %in% sig_genera)
df <- psmelt(phy_sig)
df$environment2 <- factor(df$environment2, levels = c("wild","tank"))

df_sum <- df %>%
  group_by(Genus, environment2, Time) %>%
  summarise(
    MeanAbund = mean(Abundance),
    SD        = sd(Abundance),
    N         = n(),
    SE        = SD / sqrt(N),
    .groups   = "drop"
  )

library(phyloseq); library(dplyr); library(ggplot2); library(scales); library(ggforce)

ggplot(df_sum, aes(x = Time, y = MeanAbund, color = environment2, group = environment2)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap_paginate(~ Genus, ncol = 5, nrow = 5, scales = "free_y", page = 4) +
  scale_y_continuous(
    labels = number_format(accuracy = 0.1),
    expand = expansion(c(0, 0.02))
  ) +
  scale_color_manual(values = c(wild = "#ff7f0e", tank = "#1f77b4")) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text       = element_text(size = 8, face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y      = element_text(size = 6),
    panel.spacing    = unit(0.5, "lines"),
    legend.position  = "bottom"
  ) +
  labs(
    x     = NULL,
    y     = "Relative abundance (%)",
    color = "Environment"
  )+
  geom_errorbar(aes(ymin=MeanAbund-SD, ymax=MeanAbund+SD), 
                width=0.2, position=position_dodge(0.2))+  # 白底
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    strip.background = element_rect(fill = "grey95", color = "black"),   
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 12, face = "bold")
  )


