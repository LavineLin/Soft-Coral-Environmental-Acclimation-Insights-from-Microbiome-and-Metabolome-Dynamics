setwd("/Users/linyanzhi/Desktop/PAPER-Soft coral environemntal acclimation_insight frome microbiome and metabolome dynamics/data")


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
library(Rarefy)
library(ggh4x)
library(randomcoloR)
library(forcats)
library(eulerr)
rm(list=ls())
#readfile####
feature_table_bacteria <- read.csv(file="feature_table.csv", header = T, check.names = TRUE, row.names = 1)
tax<-read.csv(file="taxa_table.csv",header=T,check.names = T, row.names = 1)
metadata_bacteria <- read.csv(file = "metadata.csv", header = T, check.names = T, row.names = 1)


#Cleansing taxonomy table
tax[tax$Class=="",]$Class <- tax[tax$Class=="",]$Phylum 

replace_class = "uncultured" 
tax[grep(pattern = replace_class,tax$Class),]$Class <- 
  tax[grep(pattern = replace_class,tax$Class),]$Phylum

tax[tax$Order=="",]$Order <- tax[tax$Order=="",]$Class

replace_order = "uncultured"
tax[grep(pattern = replace_order,tax$Order),]$Order <-
  tax[grep(pattern = replace_order,tax$Order),]$Class

tax[tax$Family=="",]$Family <- tax[tax$Family=="",]$Order

replace_family = "uncultured"
tax[grep(pattern = replace_family,tax$Family),]$Family <- 
  tax[grep(pattern =  replace_family,tax$Family),]$Order

tax[tax$Genus =="",]$Genus <- tax[tax$Genus =="",]$Family

replace_genus = "uncultured"  
tax[grep(pattern = replace_genus,tax$Genus),]$Genus <-
  tax[grep(pattern = replace_genus,tax$Genus),]$Family


tax[tax$Species =="",]$Species <- tax[tax$Species =="",]$Genus

replace_species = "uncultured"
tax[grep(pattern = replace_species,tax$Species),]$Species <-  
  tax[grep(pattern = replace_species,tax$Species),]$Genus

empty_rows<- which(tax$Phylum == ""|tax$Class == ""|tax$Order == ""|tax$Family == ""|tax$Genus == ""|tax$Species == "")

tax$Phylum[empty_rows] <- tax$Kingdom[empty_rows]
tax$Class[empty_rows] <- tax$Kingdom[empty_rows]
tax$Order[empty_rows] <- tax$Kingdom[empty_rows]
tax$Family[empty_rows] <- tax$Kingdom[empty_rows]
tax$Genus[empty_rows] <- tax$Kingdom[empty_rows]
tax$Species[empty_rows] <- tax$Kingdom[empty_rows]

#Taxonomic filtering
tax_bacteria <- as.matrix(tax)
TAX_bacteria <- tax_table(tax_bacteria)
meta_bacteria <- as.data.frame(metadata_bacteria)
META_bacteria <- sample_data(meta_bacteria)
OTU_bacteria <- otu_table(feature_table_bacteria, taxa_are_rows = TRUE)

physeq_bacteria <- phyloseq(OTU_bacteria, TAX_bacteria, META_bacteria)
physeq_bacteria <- subset_taxa(physeq_bacteria, Kingdom != "d__Eukaryota")
physeq_bacteria <- subset_taxa(physeq_bacteria, Kingdom != "Unassigned")
physeq_bacteria <- subset_taxa(physeq_bacteria, Order != "o__Chloroplast")
physeq_bacteria <- subset_taxa(physeq_bacteria, Family != "f__Mitochondria")
physeq_bacteria

#draw rarefaction curve
col <- c("darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
raremax = min(sample_sums(physeq_bacteria)) #get the minimum point

tab <- otu_table(physeq_bacteria)
class(tab) <- "matrix"
tab <- t(tab)
rarecurve(tab, step=50, cex=0.5, xlab="Sequence sample size", ylab="Species richness", label=T, col=col)+
  geom_text(title("Rarefaction Curve"))+
  abline(v=raremax) #add vertical line with the minimum point
raremax
rarefied_physeq = rarefy_even_depth(physeq_bacteria, rngseed = 1, sample.size = raremax, replace = F)

#save_physeq
saveRDS(rarefied_physeq, file = "rarefied_physeq.rds")


