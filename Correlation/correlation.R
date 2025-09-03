
library(mixOmics)
library(caret)
library(dplyr)
library(tidyverse) 
library(microeco) 
library(magrittr)
library(phyloseq)
library(vegan)
library(randomForest)
library(ggplot2)
library(dplyr)
library(tibble)
library(corrplot)
library(Hmisc)
library(RColorBrewer)
library(ggpubr)
library(ggplotify)
library(svglite)
library(tidyr)
library(reshape2)

X <- read.csv("bacteria_t.csv", row.names = 1, check.names = F)
M <- read.csv("metabolite_t.csv", row.names = 1)  
metadata_bacteria <- read.csv(file = "metadata.csv", header = T, check.names = T, row.names = 1)

taxa <- read.csv("tax_table.csv", row.names = 1)
#lobo_deseq
deseq_taxa<-c("Tyr","Glu","Thr")
#sinu_deseq
deseq_taxa<- c(
  "Tau",
  "PC.O.32.2",
  "X4.8.12.15.18.eicosapentaenoic.acid",
  "X28.6.10Z.13Z.16Z.19Z.22Z.25Z.",
  "GABA",
  "LPC.O.18.0",
  "Glu",
  "X6Z.9Z.12Z.15Z.18Z.tetracosapentaenoic.acid",
  "Docosahexaenoic.acid",
  "LPE.O.18.1",
  "DGCC.16.0_22.6",
  "Arachidonic.acid","ST.29.2.O5.3",
  "SL.10.0.O.16.0.O"
)

#sinu
lefse_taxa<- c(
  "f__Microtrichaceae",
  "g__Enhydrobacter",
  "g__Endozoicomonas",
  "g__Pantoea",
  "g__Mycoplasma",
  "g__Grimontia",
  "g__Subgroup_22",
  "g__Sphingomonas",
  "g__Mycobacterium",
  "g__Pelagibius",
  "g__37-13",
  "g__Marinoscillum",
  "c__Gammaproteobacteria",
  "g__Haliangium",
  "g__Ilumatobacter",
  "f__Kiloniellaceae",
  "f__Arenicellaceae",
  "o__Rhizobiales",
  "g__Rheinheimera",
  "g__Candidatus_Endoecteinascidia",
  "g__UBA10353_marine_group",
  "g__BD7-8",
  "f__Unknown_Family",
  "g__LD1-PA32",
  "g__Pir4_lineage",
  "f__Sandaracinaceae",
  "g__Lentisphaera",
  "g__Corynebacterium",
  "o__Actinomarinales",
  "g__Synechococcus_CC9902",
  "g__SAR116_clade",
  "g__NB1-j",
  "g__Schlegelella",
  "g__Subgroup_10",
  "f__Cyclobacteriaceae",
  "g__Ekhidna",
  "g__Thiomicrorhabdus",
  "g__OM60(NOR5)_clade",
  "o__Defluviicoccales",
  "o__Ardenticatenales",
  "g__Halioglobus",
  "f__Thiotrichaceae",
  "g__Woeseia",
  "g__Peredibacter",
  "g__Psychrobacter",
  "g__Cutibacterium",
  "g__Sulfurimonas",
  "f__Spirochaetaceae",
  "g__Candidatus_Hepatoplasma",
  "c__Bacilli",
  "f__Terasakiellaceae",
  "f__Methyloligellaceae",
  "g__Clade_III",
  "g__Aurantivirga",
  "f__Hyphomonadaceae",
  "d__Bacteria",
  "g__Erythrobacter",
  "g__HIMB11",
  "g__Flavobacterium",
  "g__Salinirepens",
  "g__Pseudoalteromonas",
  "f__Arcobacteraceae",
  "g__Alteromonas",
  "g__Candidatus_Bacilloplasma",
  "g__AEGEAN-169_marine_group",
  "g__NS4_marine_group",
  "g__Thiothrix",
  "f__Rhodobacteraceae",
  "g__SAR324_clade(Marine_group_B)",
  "g__PB19",
  "f__Cellvibrionaceae",
  "g__Shewanella",
  "g__Marinicaulis",
  "g__Ruegeria",
  "g__Nitrincolaceae",
  "g__Dasania",
  "g__Algicola",
  "g__Persicirhabdus",
  "g__Cognatishimia",
  "g__BD1-7_clade",
  "g__Tropicibacter",
  "g__Rhodobacteraceae",
  "g__Shimia",
  "g__Mesoflavibacter",
  "g__C1-B045",
  "g__Pseudobacteriovorax",
  "g__Pseudoteredinibacter",
  "g__Aestuariibacter",
  "g__028H05-P-BN-P5",
  "g__Pseudophaeobacter",
  "g__[Caedibacter]_taeniospiralis_group",
  "g__Labrenzia",
  "g__SM1A02",
  "g__Marinobacterium",
  "g__Owenweeksia"
)

#lobo


lefse_taxa<- c(
  "g__Brevibacterium",
  "g__Altererythrobacter",
  "g__SAR116_clade",
  "g__Candidatus_Hepatoplasma",
  "g__Sulfurimonas",
  "g__BD7-8",
  "g__Sulfurovum",
  "f__Spirochaetaceae",
  "g__B2M28",
  "g__Thiomicrorhabdus",
  "g__Pantoea",
  "f__Terasakiellaceae",
  "g__Endozoicomonas",
  "g__Lactobacillus",
  "c__Bacilli",
  "g__Comamonas",
  "g__Corynebacterium",
  "g__Saccharimonadales",
  "g__Sphingomonas",
  "g__Alteromonas",
  "g__Nautella",
  "g__Tenacibaculum",
  "f__Rhodobacteraceae",
  "g__Neptuniibacter",
  "g__Kordiimonas",
  "g__Acanthopleuribacter",
  "g__Tropicibacter",
  "g__Reichenbachiella",
  "g__Vibrio",
  "g__Blastocatella",
  "g__Crocinitomix",
  "g__Cellvibrio",
  "f__Hyphomonadaceae",
  "g__Cognatishimia",
  "g__Thalassobius",
  "o__Micavibrionales",
  "g__OM182_clade",
  "g__Vicingus",
  "g__Candidatus_Kaiserbacteria",
  "g__Francisella",
  "f__Micavibrionaceae",
  "g__Halioxenophilus",
  "g__Marinobacterium",
  "g__Sulfitobacter",
  "g__Marinobacter",
  "g__Aestuariicella",
  "f__Methyloligellaceae",
  "g__Seonamhaeicola",
  "g__SM2D12",
  "g__Algicola",
  "g__Pseudobacteriovorax",
  "g__Gilvimarinus",
  "f__Simkaniaceae",
  "g__Legionella",
  "f__Chlamydiaceae"
)

rarefied_physeq <- readRDS("rarefied_physeq.rds")
Coral_subset <- subset_samples(rarefied_physeq, 
                               Treatment == "Lobo" & environment2 == "tank" & Time %in% c(1, 6, 10,19))

Coral_genus <- tax_glom(Coral_subset, taxrank = "Genus")

tax_tab <- as.data.frame(tax_table(Coral_genus))


otus <- rownames(tax_tab)[ tax_tab$Genus %in% lefse_taxa ]
matched_otus <- rownames(tax_tab)[tax_tab$Genus %in% lefse_taxa]
top30_physeq <- prune_taxa(matched_otus, Coral_genus)
top30_physeq_rel <- transform_sample_counts(top30_physeq, function(x) x / sum(x))
top30_abund <- otu_table(top30_physeq_rel) %>% as.data.frame()


abund_t <- as.data.frame(t(top30_abund))  
abund_t$SampleID <- rownames(abund_t)


otus <- rownames(taxa)[ taxa$Genus %in% lefse_taxa ]

X_lefse <- X[ , colnames(X) %in% otus, drop = FALSE ]
M_selected <- M[, colnames(M) %in% deseq_taxa]

common_samples <- Reduce(intersect, list(rownames(X_lefse), rownames(M_selected), rownames(metadata_bacteria)))



M_selected <- M[, colnames(M) %in% deseq_taxa]
M_sub <- M_selected[common_samples, , drop=FALSE]
meta_sub <- metadata_bacteria[common_samples, ]
meta_lobo <- meta_sub %>% filter(Treatment == "Lobo")
meta_lobo_tank <- meta_lobo %>% filter(Environment == "wild")
samples_lobo <- rownames(meta_lobo_tank)

M_lobo <- M_sub[samples_lobo, , drop=FALSE]
M_lobo$SampleID <- rownames(M_lobo)
env_data_merged <- merge(abund_t, M_lobo, by = "SampleID")  

#lobo
env_vars<-c("Tyr","Glu","Thr")
#sinu
env_vars<-c(
  "Tau",
  "PC.O.32.2",
  "X4.8.12.15.18.eicosapentaenoic.acid",
  "X28.6.10Z.13Z.16Z.19Z.22Z.25Z.",
  "GABA",
  "LPC.O.18.0",
  "Glu",
  "X6Z.9Z.12Z.15Z.18Z.tetracosapentaenoic.acid",
  "Docosahexaenoic.acid",
  "LPE.O.18.1",
  "DGCC.16.0_22.6",
  "Arachidonic.acid",
  "ST.29.2.O5.3",
  "SL.10.0.O.16.0.O"
)


genus_vars <- rownames(top30_abund) 

env_mat <- as.matrix(env_data_merged[, env_vars])
genus_mat <- as.matrix(env_data_merged[, genus_vars])
genus_mat <- genus_mat[, apply(genus_mat, 2, var) > 0]

cor_result <- rcorr(genus_mat, env_mat, type = "spearman")

cor_matrix <- cor_result$r   # correlation value
p_matrix <- cor_result$P     # p-values

write.csv(cor_matrix, file = "cor_matrix_sinu_tank.csv")
#去excel整理
cor <- read.csv("cor_lobo_wild.csv", row.names = 1, check.names = FALSE)
p <- read.csv("p_lobo_wild.csv", row.names = 1, check.names = FALSE)


cor <- cor %>%
  rename_with(~ recode(.x,
                       "X4.8.12.15.18.eicosapentaenoic.acid"         = "Eicosapentaenoic.acid",
                       "X28.6.10Z.13Z.16Z.19Z.22Z.25Z."              = "28:6 FA",
                       "X6Z.9Z.12Z.15Z.18Z.tetracosapentaenoic.acid" = "Tetracosapentaenoic.acid"
  ))
p <- p %>%
  rename_with(~ recode(.x,
                       "X4.8.12.15.18.eicosapentaenoic.acid"         = "Eicosapentaenoic.acid",
                       "X28.6.10Z.13Z.16Z.19Z.22Z.25Z."              = "28:6 FA",
                       "X6Z.9Z.12Z.15Z.18Z.tetracosapentaenoic.acid" = "Tetracosapentaenoic.acid"
  ))

cor_df <- melt(as.matrix(cor),
               varnames = c("Genus","Metabolite"),
               value.name = "rho")

p_df   <- melt(as.matrix(p),
               varnames = c("Genus","Metabolite"),
               value.name = "P")


tax <- read.csv("tax_table.csv", row.names = 1)
cor_df$Genus_name <- tax[as.character(cor_df$Genus), "Genus"]
p_df$Genus_name <- tax[as.character(p_df$Genus), "Genus"]


deseq_taxa<- c("ST.29.2.O5.3",
  "28:6 FA",
  "Eicosapentaenoic.acid",
  "PC.O.32.2",
  "Tau",
  "Arachidonic.acid",
  "Tetracosapentaenoic.acid",
  "Glu",
  "GABA",
  "SL.10.0.O.16.0.O",
  "LPE.O.18.1",
  "DGCC.16.0_22.6"
)

cor_df$Genus_name <- factor(cor_df$Genus_name, levels = rev(lefse_taxa))
p_df$Genus_name <- factor(p_df$Genus_name, levels = rev(lefse_taxa))

cor_df$Metabolite <- factor(cor_df$Metabolite, levels = rev(deseq_taxa))
p_df$Metabolite <- factor(p_df$Metabolite, levels = rev(deseq_taxa))
# 轉成寬格式矩陣，行為 Genus_name，列為 Env_var，值是 rho
cor_mat <- cor_df %>%
  select(Genus_name, Metabolite, rho) %>%
  pivot_wider(names_from = Metabolite, values_from = rho) %>%
  column_to_rownames(var = "Genus_name") %>%
  as.matrix()


cor_mat <- cor_df %>%
  select(Genus_name, Metabolite, rho) %>%
  pivot_wider(
    names_from = Metabolite,
    values_from = rho,
    values_fn = mean  # 或其他聚合函數，如 max, min, first
  ) %>%
  column_to_rownames(var = "Genus_name") %>%
  as.matrix()

cor_mat <- cor_df %>%
  select(Genus_name, Metabolite, rho) %>%
  pivot_wider(
    names_from = Metabolite,
    values_from = rho,
    values_fn = mean
  ) %>%
  # 把 Genus_name 當成欄位
  as.data.frame() %>%
  mutate(Genus_name = factor(Genus_name, levels = lefse_taxa)) %>%
  arrange(Genus_name) %>%
  # 設成 rownames，並移除 Genus_name 欄
  column_to_rownames(var = "Genus_name") %>%
  as.matrix()

#lobo
cor_mat<- cor_df %>%
  filter(Genus_name %in% lefse_taxa) %>%
  mutate(Genus_name = factor(Genus_name, levels = lefse_taxa)) %>%
  arrange(Genus_name) %>%
  pivot_wider(names_from = Metabolite, values_from = rho) %>%
  column_to_rownames(var = "Genus_name") %>%
  as.matrix() %>%
  .[, -1] %>%
  { matrix(as.numeric(.), nrow = nrow(.), ncol = ncol(.), dimnames = dimnames(.)) }

# 同理轉換 p-value 矩陣
p_mat <- p_df %>%
  select(Genus_name, Metabolite, P) %>%
  na.omit() %>% 
  pivot_wider(names_from = Metabolite, values_from = P,
              values_fn = mean) %>%
  column_to_rownames(var = "Genus_name") %>%  # 移除含 NA 的列
  as.matrix()

p_mat<- p_df %>%
  filter(Genus_name %in% lefse_taxa) %>%
  mutate(Genus_name = factor(Genus_name, levels = lefse_taxa)) %>%
  arrange(Genus_name) %>%
  pivot_wider(names_from = Metabolite, values_from = P) %>%
  column_to_rownames(var = "Genus_name") %>%
  as.matrix() %>%
  .[, -1] %>%
  { matrix(as.numeric(.), nrow = nrow(.), ncol = ncol(.), dimnames = dimnames(.)) }



library(corrplot)
library(RColorBrewer)
library(fields)  # for image.plot

# 你的顏色調色盤
col_pal <- rev(colorRampPalette(brewer.pal(n = 7, "RdBu"))(200))

# 1. 畫主圖，不顯示色標 (cl.pos = "n")
corrplot(cor_mat, p.mat = p_mat, method = "circle",
         col = col_pal,
         tl.cex = 0.5,
         tl.col = "black", tl.srt = 90,
         type = "full",
         sig.level = c(0.001, 0.01, 0.05),
         insig = "label_sig",
         pch.cex = 0.8, pch.col = "black",
         cl.pos = "n",      # 不顯示色標
         is.corr = TRUE)


# 2. 單獨畫色標（scale bar）
# 設定繪圖區域大小跟位置（可調整）
par(mar=c(5,4,4,2)+0.1)  # 調整邊界

image.plot(legend.only = TRUE, zlim = c(-1, 1), col = col_pal,
           legend.lab = "Correlation",
           legend.line = 2.5,
           legend.width = 0.7)


library(corrplot)
library(RColorBrewer)

# 調色盤
col_pal <- rev(colorRampPalette(brewer.pal(7, "RdBu"))(200))

# 1) 畫底圖：所有相關係數的彩色實心圓
corrplot(
  cor_mat,
  p.mat    = p_mat,
  method   = "circle",
  col      = col_pal,
  tl.cex   = 0.5,
  tl.col   = "black",
  tl.srt   = 90,
  type     = "full",
  cl.pos   = "n",
  is.corr  = TRUE,
  insig    = "blank"    # 完全不要標示不顯著
)


met_order <- c("ST.29.2.O5.3",
               "28:6 FA",
               "Eicosapentaenoic.acid",
               "PC.O.32.2",
               "Tau",
               "Arachidonic.acid",
               "Tetracosapentaenoic.acid",
               "Glu",
               "GABA",
               "SL.10.0.O.16.0.O",
               "LPE.O.18.1",
               "DGCC.16.0_22.6"
)

# 1) 先把 cor_mat 和 p_mat 的列都按这个顺序重排
cor_mat_ord <- cor_mat[, met_order]
p_mat_ord   <- p_mat[,   met_order]

corrplot(
  cor_mat,
  method   = "circle",
  col      = col_pal,
  tl.cex   = 0.5,
  tl.col   = "black",
  tl.srt   = 90,
  type     = "full",
  cl.pos   = "n",     # 不顯示色標
  is.corr  = TRUE
)

corrplot(
  cor_mat,
  method   = "circle",
  col      = col_pal,
  tl.cex   = 0.5,
  tl.col   = "black",
  tl.srt   = 90,
  type     = "full",
  cl.pos   = "n",     # 不顯示色標
  is.corr  = TRUE
)

# 2) 疊加：只對 p<0.05 的格子，用空心方框畫出來
sig_idx <- which(p_mat < 0.05, arr.ind = TRUE)
n        <- nrow(cor_mat)

# corrplot 的座標系：x = 列號 (Metabolite)，y = row 反向 (Genus)
for(k in seq_len(nrow(sig_idx))) {
  i <- sig_idx[k, 1]      # row index in cor_mat
  j <- sig_idx[k, 2]      # col index in cor_mat
  # draw 空心方框 (pch=0)，cex 大小微調到蓋滿格子
  points(
    x    = j,
    y    = n - i + 1,
    pch  = 0,
    cex  = 1.3,
    lwd  = 1.5,
    col  = "black"
  )
}

for(k in seq_len(nrow(sig_idx))) {
  i <- sig_idx[k, 1]      # row index in cor_mat
  j <- sig_idx[k, 2]      # col index in cor_mat
  # draw 空心方框 (pch=0)，cex 大小微調到蓋滿格子
  points(
    x    = j,
    y    = n - i + 1,
    pch  = 0,
    cex  = 2,
    lwd  = 1.5,
    col  = "black"
  )
}
