setwd("/Users/linyanzhi/Desktop/PAPER-Soft coral environemntal acclimation_insight frome microbiome and metabolome dynamics/data")

rm(list=ls()) 
#library####
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
library(ggh4x)
library(randomcoloR)
library(forcats)
#library(MicEco)
library(eulerr)
library(pheatmap)
rarefied_physeq <- readRDS("rarefied_physeq.rds")
lobophyseq<-subset_samples(rarefied_physeq,Treatment=="Lobo")
sinuphyseq<-subset_samples(rarefied_physeq,Treatment=="Sclero")

#bacteria_top50 genus####
physeq_genus <- tax_glom(lobophyseq, "Genus")
lobo_top <- prune_taxa(names(sort(taxa_sums(physeq_genus),TRUE)[1:50]), physeq_genus)

lobo_topz <- microbiome::transform(lobo_top, "Z")

genus_table <- as.data.frame(otu_table(lobo_topz))
sample_table <- as.data.frame(sample_data(lobo_topz))
tax_table <- as.data.frame(phyloseq::tax_table(lobo_topz))
row.names(genus_table) <- tax_table$Genus
row.names(genus_table) <- paste(tax_table$Family, tax_table$Genus, sep = " | ")
colnames(genus_table) <- sample_table$All_type

tank10 <- genus_table[, 1:3]
tank11 <- genus_table[, 4:6]
tank12 <- genus_table[, 7:9]
tank19 <- genus_table[, 10:12]
tank1 <- genus_table[, 13:15]
tank2 <- genus_table[, 16:18]
tank3 <- genus_table[, 19:21]
tank4 <- genus_table[, 22:24]
tank5 <- genus_table[, 25:27]
tank6 <- genus_table[, 28:30]
tank7 <- genus_table[, 31:33]
tank8 <- genus_table[, 34:36]
tank9 <- genus_table[, 37:39]
wild0<- genus_table[, 40:42]
wild6 <- genus_table[, 43:45]
wild19 <- genus_table[, 46:48]

heatmap_data <- data.frame(wild0 = rowMeans(wild0),tank1 = rowMeans(tank1),tank2 = rowMeans(tank2),tank3 = rowMeans(tank3),tank4 = rowMeans(tank4),tank5 = rowMeans(tank5),wild6 = rowMeans(wild6),tank6 = rowMeans(tank6),tank7 = rowMeans(tank7),tank8 = rowMeans(tank8),tank9 = rowMeans(tank9),tank10 = rowMeans(tank10), tank11 = rowMeans(tank11),tank12 = rowMeans(tank12),wild19 = rowMeans(wild19),tank19 = rowMeans(tank19))

write.csv(heatmap_data,"lobo_bac_heat.csv")


library(tidyverse)
library(cowplot)
library(pheatmap)
library(dplyr)
rm(list=ls())
# readCSV
laa_df    <- read.csv("lobo_aa_z_heat.csv",     row.names=1, check.names=FALSE)
llip_df   <- read.csv("lobo_lipid_heat_log.csv",  row.names=1, check.names=FALSE)
saa_df    <- read.csv("sinu_aa_z_heat.csv",     row.names=1, check.names=FALSE)
slip_df   <- read.csv("sinu_lipid_heat_log.csv",  row.names=1, check.names=FALSE)
lbac_df    <- read.csv("lobo_bac_heat.csv",     row.names=1, check.names=FALSE)
sbac_df    <- read.csv("sinu_bac_heat.csv",     row.names=1, check.names=FALSE)
# order
stage_order <- c(paste0("tank",c(1,2,3,4,5,6,7,8,9,10,11,12,19)), paste0("wild", c(0,6,19)))

#add_missing_col  
add_missing_cols <- function(df, cols){
  missing <- setdiff(cols, colnames(df))
  df[missing] <- NA
  df <- df[, cols]
  return(df)
}

laa2  <- add_missing_cols(laa_df,  stage_order)
llip2 <- add_missing_cols(llip_df, stage_order)
saa2  <- add_missing_cols(saa_df,  stage_order)
slip2 <- add_missing_cols(slip_df, stage_order)
lbac_df2  <- add_missing_cols(lbac_df,  stage_order)
sbac_df2 <- add_missing_cols(sbac_df, stage_order)


tank_cols <- grep("^tank", colnames(slip2), value=TRUE)

pt <- pheatmap(
  slip2[, tank_cols],
  cluster_rows   = TRUE,
  cluster_cols   = FALSE,        
  silent         = TRUE,
  cellwidth      = 15,
  cellheight     = 10,
  fontsize_row   = 12,
  fontsize_col   = 12,
  color          = colorRampPalette(c("lightblue","white","darkred"))(50),
  na_col         = "white"
)
pt

#row dendrogram
row_hclust <- pt$tree_row

#visulize
pheatmap(
  slip2,
  cluster_rows   = row_hclust,
  cluster_cols   = FALSE,     
  cellwidth      = 15,
  cellheight     = 10,
  fontsize_row   = 12,
  fontsize_col   = 12,
  color          = colorRampPalette(c("lightblue","white","darkred"))(50),
  na_col         = "white"
)
