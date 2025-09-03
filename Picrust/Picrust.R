#picrust_data_preprocessing
library(dplyr)
library(tidyr)
df <- read.table("pred_metagenome_contrib_sinu.tsv",
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# calculatetaxon_rel_function_abundance sum
df_sum <- df %>%
  group_by(sample, function.) %>%
  summarise(total_rel_abun = sum(taxon_rel_function_abun, na.rm = TRUE)) %>%
  ungroup()

df_wide <- df_sum %>%
  pivot_wider(
    names_from  = sample,
    values_from = total_rel_abun,
    values_fill = 0
  )

write.table(df_wide, "ko_abundance_rel_sinu.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
library(MARco)
if (!requireNamespace("BiocManager", quietly =TRUE))                 
  install.packages("BiocManager")                 

pkgs <-c("gage", "vegan", "fossil", "vioplot", "igraph",                 
         "SpiecEasi", "DESeq2", "pheatmap", "MASS",                 
         "gage", "cluster", "clusterSim")                 

for (pkg in pkgs) {                 
  if (!requireNamespace(pkg, quietly =TRUE))                 
    BiocManager::install(pkg)                 
}
library(gage)
abundance_file <- read.table("ko_abundance_rel_sinu.tsv", header = TRUE, sep = "\t", check.names = FALSE,row.names = 1)

#Convert to metabolic pathways
gs.table <- KO2path(predtable = abundance_file,GeneSet = "pathway") 
head(gs.table)
#gs.table.p <- KO2path(predtable = abundance_file,GeneSet = "metabolic",prop = TRUE) # present with proportion
#head(gs.table.p)

#Convert to KEGG modules
md.table <- KO2path(predtable = abundance_file,GeneSet = "module")
head(md.table)
#Annotate with upper levels
lv.table <- path.lvs(gs.table)
head(lv.table)

write.csv(lv.table, file = "picrust_table_rel_sinu.csv", row.names = TRUE)


#readfile####
feature_table_bacteria <- read.csv(file="picrust_table_rel_sinu.csv", header = T, check.names = F, row.names = 1)
tax<-read.csv(file="lv_tax.csv",header=T,check.names = T, row.names = 1)
metadata_bacteria <- read.csv(file = "metadata.csv", header = T, check.names = T, row.names = 1)


tax_bacteria <- as.matrix(tax)
TAX_bacteria <- tax_table(tax_bacteria)
meta_bacteria <- as.data.frame(metadata_bacteria)
META_bacteria <- sample_data(meta_bacteria)
#colnames(feature_table_bacteria) <- row.names(meta_bacteria)
OTU_bacteria <- otu_table(feature_table_bacteria, taxa_are_rows = TRUE)

physeq_bacteria <- phyloseq(OTU_bacteria, TAX_bacteria, META_bacteria)

physeq_bacteria


#匯出檔案
saveRDS(physeq_bacteria, file = "picrust_rel_sinu.rds")
#heatmap####
rm(list=ls())
library(pheatmap)
picrust <- readRDS("picrust_rel_lobo.rds")

picrust_level<-subset_samples(picrust, Treatment =="Lobo")


lipid_metabolism <- subset_taxa(picrust_level,Level_2 %in% c("1.3 Lipid metabolism","1.5 Amino acid metabolism","1.6 Metabolism of other amino acids"))
lipid_metabolism <- subset_taxa(picrust_level,Level_3 %in% c("Biosynthesis of secondary metabolites"))

# zscore 
lobo_topz <- microbiome::transform(lipid_metabolism, "Z")
phyloseq::otu_table(lipid_metabolism) <- phyloseq::otu_table(lobo_topz, taxa_are_rows = TRUE)

genus_table <- as.data.frame(otu_table(lobo_topz))
sample_table <- as.data.frame(sample_data(lobo_topz))
sample_table$SampleID <- rownames(sample_table)
tax_table <- as.data.frame(phyloseq::tax_table(lobo_topz))
row.names(genus_table) <- tax_table$Level_3
colnames(genus_table) <- sample_table$Env_Time

x_axis_order<-c("wild_0","tank_1","tank_2","tank_3","tank_4","tank_5","wild_6","tank_6","tank_7","tank_8","tank_9","tank_10","tank_11","tank_12","tank_19","wild_19")
x_axis_order <- as.character(x_axis_order)
#lobo
#mean
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
wild0 <- genus_table[, 40:42]
wild6 <- genus_table[, 43:45]
wild19 <- genus_table[, 46:48]


#sinu
tank10 <- genus_table[, 1:3]
tank11 <- genus_table[, 4:6]
tank12 <- genus_table[, 7:9]
tank19 <- genus_table[, 10:12]
tank1 <- genus_table[, 13:15]
tank2 <- genus_table[, 16:17]
tank3 <- genus_table[, 18:20]
tank4 <- genus_table[, 21:23]
tank5 <- genus_table[, 24:25]
tank6 <- genus_table[, 26:28]
tank7 <- genus_table[, 29:31]
tank8 <- genus_table[, 32:34]
tank9 <- genus_table[, 35:37]
wild0<- genus_table[, 38:40]
wild6 <- genus_table[, 41:43]
wild19 <- genus_table[, 44:46]



heatmap_data <- data.frame(
  wild0 = rowMeans(wild0), tank1 = rowMeans(tank1), tank2 = rowMeans(tank2),
  tank3 = rowMeans(tank3), tank4 = rowMeans(tank4), tank5 = rowMeans(tank5),
  wild6 = rowMeans(wild6), tank6 = rowMeans(tank6), tank7 = rowMeans(tank7),
  tank8 = rowMeans(tank8), tank9 = rowMeans(tank9), tank10 = rowMeans(tank10),
  tank11 = rowMeans(tank11), tank12 = rowMeans(tank12), wild19 = rowMeans(wild19),
  tank19 = rowMeans(tank19)
)

#order
stage_order <- c(paste0("tank",c(1,2,3,4,5,6,7,8,9,10,11,12,19)), paste0("wild", c(0,6,19)))

#add)missing
add_missing_cols <- function(df, cols){
  missing <- setdiff(cols, colnames(df))
  df[missing] <- NA
  df <- df[, cols]
  return(df)
}

heatmap_data2  <- add_missing_cols(heatmap_data,  stage_order)

p1<-pheatmap(
  heatmap_data2,
  cluster_cols = F,cluster_rows = F,
  treeheight_row = 30,
  height = 7,
  width = 12,   # 整體寬度
  cellwidth      = 15,
  cellheight     = 10,
  fontsize_row   = 12,
  fontsize_col   = 12,
  color          = colorRampPalette(c("lightblue","white","darkred"))(50),
  na_col         = "white"
)
p1



tank_cols <- grep("^tank", colnames(heatmap_data2), value=TRUE)
pt <- pheatmap(
  heatmap_data2[, tank_cols],
  cluster_rows   = TRUE,
  cluster_cols   = FALSE,         # 如果你要切群可以保留
  silent         = TRUE,
  cellwidth      = 15,
  cellheight     = 10,
  fontsize_row   = 12,
  fontsize_col   = 12,
  color          = colorRampPalette(c("lightblue","white","darkred"))(50),
  na_col         = "white"
)
pt

row_hclust <- pt$tree_row

pheatmap(
  heatmap_data2,
  cluster_rows   = row_hclust,
  cluster_cols   = FALSE,       # or remove if 不要切群
  cellwidth      = 15,
  cellheight     = 10,
  fontsize_row   = 12,
  fontsize_col   = 12,
  color          = colorRampPalette(c("lightblue","white","darkred"))(50),
  na_col         = "white"
)

