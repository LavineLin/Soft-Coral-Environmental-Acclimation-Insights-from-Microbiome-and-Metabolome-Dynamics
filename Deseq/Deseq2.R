
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
library(Rarefy)
library(ggh4x)
library(randomcoloR)
library(forcats)
#library(MicEco)
library(eulerr)
rm(list=ls())
Amino <- readRDS("Amino_acid.rds")
L_amino<-subset_samples(Amino, Treatment =="Lobo")
S_amino<-subset_samples(Amino, Treatment =="Sclero")
L_lipid<-readRDS("L_lipidomic.rds")
S_lipid<-readRDS("S_lipidomic.rds")


L_metabolite<-merge_phyloseq(L_amino,L_lipid)
S_metabolite<-merge_phyloseq(S_amino,S_lipid)

L_metabolite<-merge_phyloseq(L_amino,L_lipid)
L_metabolite_1<-subset_samples(L_metabolite,  Group == "week1" | Group == "week0")
S_metabolite_1<-subset_samples(S_metabolite,  Group == "week1" | Group == "week0")

L_metabolite_6<-subset_samples(L_metabolite,  Time == "19" )
S_metabolite_6<-subset_samples(S_metabolite,  Time == "6")
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(ggrepel)
sample_data(L_metabolite_6)$environment2 <- factor(
  sample_data(L_metabolite_6)$environment2,
  levels = c("tank", "wild")
)

dds <- phyloseq_to_deseq2(L_metabolite_6, ~ environment2)

dds <- DESeq(dds, fitType = "parametric")

res <- results(
  dds,
  contrast = c("environment2", "wild", "tank"),
  alpha    = 0.05         
)

df <- as.data.frame(res)
df$FeatureID <- rownames(df)
df$minusLog10P <- -log10(df$pvalue)

df$Significance <- "NotSig"
df$Significance[df$pvalue < 0.05 & df$log2FoldChange >=  1] <- "Up"    # wild > tank
df$Significance[df$pvalue < 0.05 & df$log2FoldChange <= -1] <- "Down"  # wild < tank


old_ids <- rownames(df)

new_ids <- recode(old_ids,
                  "X4.8.12.15.18.eicosapentaenoic.acid"         = "Eicosapentaenoic.acid",
                  "X28.6.10Z.13Z.16Z.19Z.22Z.25Z."              = "28:6 FA",
                  "X6Z.9Z.12Z.15Z.18Z.tetracosapentaenoic.acid" = "Tetracosapentaenoic.acid"
                  
)
rownames(df) <- new_ids

df <- df %>%
  mutate(
    FeatureID = recode(FeatureID,
                       "X4.8.12.15.18.eicosapentaenoic.acid"         = "Eicosapentaenoic.acid",
                       "X28.6.10Z.13Z.16Z.19Z.22Z.25Z."              = "28:6 FA",
                       "X6Z.9Z.12Z.15Z.18Z.tetracosapentaenoic.acid" = "Tetracosapentaenoic.acid"
    )
  )
# color
cols <- c(
  "Up"     = "#ff7f0e",  # (wild > tank)
  "Down"   = "#1f77b4",  #  (wild < tank)
  "NotSig" = "grey70"
)
p_volcano <- ggplot(df, aes(x = log2FoldChange, y = minusLog10P)) +
  geom_point(aes(color = Significance),
             size = 8, alpha = 0.9) +
  scale_color_manual(values = cols) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed",
             color = "grey50") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "grey50") +
  theme_minimal(base_size = 14) +
  labs(
    x     = "Log Fold Change ",
    y     = "-Log(p-value)",
    color = ""
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )


sig_df <- subset(df, Significance%in% c("Up", "Down"))

p_volcano_labeled <- p_volcano +
  geom_text_repel(
    data          = sig_df,
    aes(label     = FeatureID),
    size          = 8,           
    max.overlaps  = 15,         
    box.padding   = 0.3,         
    point.padding = 0.8,          
    segment.size  = 0.3,          
    segment.color = "grey50"      
  )+
  theme_bw(base_size = 14) +
  theme(
    axis.title.x = element_text(size =22), 
    axis.title.y = element_text(size = 22),  
    axis.text.x  = element_text(size = 18),   
    axis.text.y  = element_text(size = 18),   
    panel.border     = element_rect(color = "black", fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(hjust = 0.5)
  )

print(p_volcano_labeled)

