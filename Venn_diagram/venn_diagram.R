
rm(list=ls())

rarefied_physeq <- readRDS("rarefied_physeq.rds")
rarefied_physeq_tank<-subset_samples(rarefied_physeq, environment2 =="tank")
water<-subset_samples(rarefied_physeq_tank, Treatment =="Sea water")
lobo<-subset_samples(rarefied_physeq_tank, Treatment =="Lobo")
sclero<-subset_samples(rarefied_physeq_tank, Treatment =="Sclero")
tank<-merge_phyloseq(water,lobo,sclero)


library(phyloseq)
library(ggvenn)

#  OTU table
otu_table_water <- otu_table(water)
otu_table_lobo <- otu_table(lobo)
otu_table_sclero <- otu_table(sclero)

if (taxa_are_rows(rarefied_physeq)) {
  otu_table_water <- t(otu_table_water)
  otu_table_lobo <- t(otu_table_lobo)
  otu_table_sclero <- t(otu_table_sclero)
}

# calculate presence OTUs
otus_water <- colnames(otu_table_water)[colSums(otu_table_water) > 0]
otus_lobo <- colnames(otu_table_lobo)[colSums(otu_table_lobo) > 0]
otus_sclero <- colnames(otu_table_sclero)[colSums(otu_table_sclero) > 0]


venn_list <- list(
  "Sea water" = otus_water,
  "Sclerophytum" = otus_sclero,
  "Lobophytum"=otus_lobo
)

venn_list <- list(
  "Sea water" = otus_water,
  "Lobophytum"=otus_lobo
)

venn_list <- list(
  "Sea water" = otus_water,
  "Sclerophytum" = otus_sclero
)
# visulize
ggvenn(
  venn_list,
  fill        = c("#A9CDEB", "#F6D49B", "#B4E4B2"),  # 低飽和色
  fill_alpha  = 0.5,
  stroke_size = 0.4,   # 畫圈邊線稍加粗
  set_name_size = 5,   # 調大集合名稱
  text_size   = 6      # 調大交集數字
) +
  theme(text = element_text(size = 14)) 

ggvenn(
  venn_list,
  fill        = c("#A9CDEB", "#F6D49B"),  # 低飽和色
  fill_alpha  = 0.5,
  stroke_size = 0.4,   # 畫圈邊線稍加粗
  set_name_size = 5,   # 調大集合名稱
  text_size   = 6      # 調大交集數字
) +
  theme(text = element_text(size = 14)) 

ggvenn(
  venn_list,
  fill        = c("#A9CDEB", "#B4E4B2"),  # 低飽和色
  fill_alpha  = 0.5,
  stroke_size = 0.4,   # 畫圈邊線稍加粗
  set_name_size = 5,   # 調大集合名稱
  text_size   = 6      # 調大交集數字
) +
  theme(text = element_text(size = 14)) 


