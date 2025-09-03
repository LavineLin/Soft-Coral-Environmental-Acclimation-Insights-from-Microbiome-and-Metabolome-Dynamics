rm(list=ls())
library(tidyverse)
library(microeco)
library(magrittr) 


feature_tables <- read.csv('featuretable_sinu.csv', row.names = 1)
sample_tables <- read.csv('metadata_sinu.csv',row.names = 1)
feature_tablel <- read.csv('featuretable_lobo.csv', row.names = 1)
sample_tablel <- read.csv('metadata_lobo.csv',row.names = 1)


feature_tablel <- feature_tablel %>%
  mutate(across(everything(),
                ~ as.numeric(as.character(.))))

feature_tables <- feature_tables %>%
  mutate(across(everything(),
                ~ as.numeric(as.character(.))))


feature_tablel[is.na(feature_tablel)] <- 0
keep_cols <- colSums(feature_tablel) > 0
feature_tablel     <- feature_tablel[, keep_cols]

feature_tables[is.na(feature_tables)] <- 0
keep_cols <- colSums(feature_tables) > 0
feature_tables     <- feature_tables[, keep_cols]


subsetlobo_metadata <- subset(sample_tablel,  Group == "week1" | Group == "week0")
subsetsinu_metadata <- subset(sample_tables,  Group == "week1" | Group == "week0")

subsetlobo_metadata2 <- subset(sample_tablel,  Group == "week6" )
subsetsinu_metadata2 <- subset(sample_tables,  Group == "week6" )

subsetlobo_metadata3 <- subset(sample_tablel,  Group == "week19" )
subsetsinu_metadata3 <- subset(sample_tables,  Group == "week19" )



datasets1 <- microtable$new(sample_table = subsetsinu_metadata,
                           otu_table = feature_tables, 
                           tax_table = tax_table)

datasets2 <- microtable$new(sample_table = subsetsinu_metadata2,
                           otu_table = feature_tables, 
                           tax_table = tax_table)

datasets3 <- microtable$new(sample_table = subsetsinu_metadata3,
                           otu_table = feature_tables, 
                           tax_table = tax_table)

lefsel3 <- trans_diff$new(dataset = datasets3, 
                         method = "lefse",
                         group = "Group3", 
                         alpha = 0.05, 
                         taxa_level = "Genus", 
                         lefse_min_subsam = 8,
                         p_adjust_method = "none"
)

lefsel2 <- trans_diff$new(dataset = datasets2, 
                         method = "lefse", 
                         group = "Group2",
                         alpha = 0.05,
                         taxa_level = "Genus", 
                         lefse_min_subsam = 8,
                         p_adjust_method = "none"
)

lefsel1 <- trans_diff$new(dataset = datasets1, 
                         method = "lefse",
                         group = "Group1", 
                         alpha = 0.05,
                         taxa_level = "Genus",
                         lefse_min_subsam = 8,
                         p_adjust_method = "none"
)

library(dplyr)
library(tidyr)
library(ggplot2)

#LDA>2.5
df1 <- lefsel1$res_diff%>%
  filter(LDA > 2.5) %>% 
  select(Taxa, LDA, Group) %>% 
  rename(
    LDA_lefses1   = LDA,
    Group_lefses1 = Group
  )

df2 <- lefsel2$res_diff %>%
  filter(LDA > 2.5)%>% 
  select(Taxa, LDA, Group) %>% 
  rename(
    LDA_lefses2   = LDA,
    Group_lefses2 = Group
  )

df3 <- lefsel3$res_diff%>%
  filter(LDA > 2.5) %>% 
  select(Taxa, LDA, Group) %>% 
  rename(
    LDA_lefses3   = LDA,
    Group_lefses3 = Group
  )

df_merged <- df1 %>%
  full_join(df2, by = "Taxa") %>%
  full_join(df3, by = "Taxa")

df_merged2 <- df_merged %>%
  replace_na(list(
    LDA_lefses1   = 0,
    LDA_lefses2   = 0,
    LDA_lefses3   = 0,
    Group_lefses1 = "NotSig",
    Group_lefses2 = "NotSig",
    Group_lefses3 = "NotSig"
  )) %>%
  separate(
    col   = Taxa,
    into  = c("Kingdom","Phylum","Class","Order","Family","Genus"),
    sep   = "\\|",     
    fill  = "right",  
    remove = TRUE       
  )

write_csv(df_merged2,"sinu_lefse.csv")



library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

df <- read.csv("sinu_lefse.csv", check.names=FALSE, stringsAsFactors=FALSE)

df2 <- df %>% 
  select(
    Genus,
    Week0 = Week0,
    Week6 = Week6,
    Week19 = Week19,
    Group_lefses1,   
    Group_lefses2, 
    Group_lefses3    
  )

df_long <- df2 %>%
  pivot_longer(
    cols      = c(Week0, Week6, Week19),
    names_to  = "Time",
    values_to = "LDA_score"
  ) %>%
  mutate(
    Time  = factor(Time, levels = c("Week0","Week6","Week19")),
    Panel = case_when(
      Time == "Week0"  ~ Group_lefses1,
      Time == "Week6"  ~ Group_lefses2,
      Time == "Week19" ~ Group_lefses3
    )
  ) %>%
  filter(!is.na(Panel))%>%
  filter(Panel != "NotSig") 




library(ggplot2)

df_wild <- df_long %>% filter(Panel == "Wild")
df_tank <- df_long %>% filter(Panel == "Tank")

df_wild2 <- df_wild %>%
  group_by(Genus) %>%
  mutate(maxLDA = max(LDA_score)) %>%
  ungroup() %>%
  mutate(Genus = fct_reorder(Genus, maxLDA))

df_tank2 <- df_tank %>%
  group_by(Genus) %>%
  mutate(maxLDA = max(LDA_score)) %>%
  ungroup() %>%
  mutate(Genus = fct_reorder(Genus, maxLDA))




df_tank2_sorted <- df_tank2 %>%
  mutate(Time = factor(Time, levels = c("Week0", "Week6", "Week19"))) %>%
  arrange(Time, LDA_score) %>%
  mutate(Genus = fct_inorder(Genus))

#  Tank
p_tank_sorted <- ggplot(df_tank2_sorted, aes(x = LDA_score, y = Genus, fill = Time)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey50") +
  scale_fill_manual(
    name   = "Time",
    values = c("Week0"="#4C72B0", "Week6"="#55A868", "Week19"="#C44E52")
  ) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "LDA score", y = NULL) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 12),
    legend.position    = "top"
  )

p_tank_sorted


df_wild2_sorted <- df_wild2 %>%
  mutate(Time = factor(Time, levels = c("Week0", "Week6", "Week19"))) %>%
  arrange(Time, desc(LDA_score)) %>%
  mutate(Genus = fct_inorder(Genus))

#wild
p_wild_sorted <- ggplot(df_wild2_sorted, aes(x = LDA_score, y = Genus, fill = Time)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey50") +
  scale_fill_manual(
    name   = "Time",
    values = c("Week0"="#4C72B0", "Week6"="#55A868", "Week19"="#C44E52")
  )+
  scale_x_continuous(expand = c(0,0)) +
  labs( x = "LDA score", y = NULL) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 12),
    legend.position    = "top"
  )

library(patchwork)
combined_plot <- p_wild_sorted + p_tank_sorted + 
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = 'top')

combined_plot

