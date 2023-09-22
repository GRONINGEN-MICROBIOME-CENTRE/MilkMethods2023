###############################################################################
##### Plotting of shared genera (at least 1 pair) in different sample comparisons
### Author(s): Asier Fernández / Johanne E. Spreckels 
### Last updated: 21st September, 2023
################################################################################

#****************
# Load libraries
#****************
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)

# Set working directory
setwd("~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/")

# Read tables with the number of pairs sharing genera (minimum shared in 1 pair)
data_milk_feces_mother <- read_excel("Shared_genera_barplots/Shared_genera_1_Pair.xlsx", sheet = "Milk-Feces Mother")
data_milk_oral_baby <- read_excel("Shared_genera_barplots/Shared_genera_1_Pair.xlsx", sheet = "Milk Mother-Oral Baby")
data_milk_feces_baby <- read_excel("Shared_genera_barplots/Shared_genera_1_Pair.xlsx", sheet = "Milk-Feces Baby")

#Rename F0428 to Lachnospiraceae_F0428
data_milk_oral_baby$Genus[7] <- "Lachnospiraceae_F0428"

# Keep the data tables in a list
data_tables <- list(data_milk_feces_mother, data_milk_oral_baby, data_milk_feces_baby)
names(data_tables) <- c("Milk-Feces Mother", "Milk-Oral Baby", "Milk-Feces Baby")


###########################
# Data processing
###########################
#Reformat each data table to 
# (A) Convert it to long format including the counts and the total number of pairs per genus
DFs_long_format <- lapply(data_tables, function(df) {
  df %>%
    pivot_longer(cols = c(`Pairs sharing`, `Pairs not sharing`), names_to = "Sharing", values_to = "count")  %>% 
    group_by(Genus, `Sequencing method`) %>% 
    reframe(Sharing= Sharing,
            count = count,
            total = sum(count))  %>%
    mutate(`Sequencing method` = as.factor(`Sequencing method`))
})


#*****************
#A) Milk-Oral Baby
#*****************
genus_order_milk_oral_16S23S <- c("Streptococcus", "Rothia", "Gemella", "Pauljensenia", "Veillonella",
                                  "Alloprevotella", "Granulicatella", "Leptotrichia", "Pseudoleptotrichia","Staphylococcus")
genus_order_milk_oral_16S <- c("Streptococcus", "Veillonella", "Haemophilus", "Pauljensenia", "Gemella", "Rothia", 
                               "Prevotella", "Staphylococcus", "Cutibacterium", "Lancefieldella", "Actinomyces",
                               "Alloprevotella", "Corynebacterium", "Neisseria", "Streptobacillus", "Bifidobacterium",
                               "Campylobacter", "Lachnospiraceae_F0428", "Leptotrichia", "Porphyromonas", 
                               "Saccharimonas")
                         
                               
Milk_oral_16S23S <- DFs_long_format$`Milk-Oral Baby`[DFs_long_format$`Milk-Oral Baby`$`Sequencing method` == "16S23S", ]
Milk_oral_16S <- DFs_long_format$`Milk-Oral Baby`[DFs_long_format$`Milk-Oral Baby`$`Sequencing method` == "V3V4", ]

Milk_oral_16S23S_ordered <- Milk_oral_16S23S[order(match(Milk_oral_16S23S$Genus, genus_order_milk_oral_16S23S)),]
Milk_oral_16S_ordered <- Milk_oral_16S[order(match(Milk_oral_16S$Genus, genus_order_milk_oral_16S)),]

DFs_long_format$`Milk-Oral Baby` <- rbind(Milk_oral_16S23S_ordered, Milk_oral_16S_ordered)


#******************
#B) Milk-Feces Baby
#******************
genus_order_Milk_feces_baby_16S23S <- c("Streptococcus", "Staphylococcus", "Cutibacterium", "Lancefieldella", "Veillonella")
genus_order_Milk_feces_baby_16S <- c("Streptococcus", "Veillonella", "Staphylococcus", "Haemophilus","Rothia", "Bifidobacterium",
                                     "Pauljensenia","Prevotella", "Clostridium", "Corynebacterium", "Cutibacterium", "Gemella",
                                     "Klebsiella", "Lancefieldella")

Milk_feces_baby_16S23S <- DFs_long_format$`Milk-Feces Baby`[DFs_long_format$`Milk-Feces Baby`$`Sequencing method` == "16S23S", ]
Milk_feces_baby_16S <- DFs_long_format$`Milk-Feces Baby`[DFs_long_format$`Milk-Feces Baby`$`Sequencing method` == "V3V4", ]

Milk_feces_baby_16S23S_ordered <- Milk_feces_baby_16S23S[order(match(Milk_feces_baby_16S23S$Genus, genus_order_Milk_feces_baby_16S23S)),]
Milk_feces_baby_16S_ordered <- Milk_feces_baby_16S[order(match(Milk_feces_baby_16S$Genus, genus_order_Milk_feces_baby_16S)),]

DFs_long_format$`Milk-Feces Baby` <- rbind(Milk_feces_baby_16S23S_ordered, Milk_feces_baby_16S_ordered)

#********************
#C) Milk-Feces Mother
#********************
DFs_long_format$`Milk-Feces Mother` <- DFs_long_format$`Milk-Feces Mother`[c(9,10,11,12,1,2,3,4,5,6,7,8),]

# Add new factor levels to reorder
genus_order_milk_oral_baby <- unique(DFs_long_format$`Milk-Oral Baby`$Genus)
DFs_long_format$`Milk-Oral Baby`$Genus <- factor(DFs_long_format$`Milk-Oral Baby`$Genus, levels = genus_order_milk_oral_baby)

genus_order_milk_feces_baby <- unique(DFs_long_format$`Milk-Feces Baby`$Genus)
DFs_long_format$`Milk-Feces Baby`$Genus <- factor(DFs_long_format$`Milk-Feces Baby`$Genus, levels = genus_order_milk_feces_baby)

genus_order_milk_feces_mother <- unique(DFs_long_format$`Milk-Feces Mother`$Genus)
DFs_long_format$`Milk-Feces Mother`$Genus <- factor(DFs_long_format$`Milk-Feces Mother`$Genus, levels = genus_order_milk_feces_mother)

###########################
# Generate barplots
###########################

#******************
#A) Milk-Oral Infant
#******************
DF_Milk_Oral_Baby_V3V4_plot <- DFs_long_format$`Milk-Oral Baby`[DFs_long_format$`Milk-Oral Baby`$`Sequencing method`=="V3V4",]
DF_Milk_Oral_Baby_V3V4_plot$Genus <- factor(DF_Milk_Oral_Baby_V3V4_plot$Genus,
                                                 levels = unique(DF_Milk_Oral_Baby_V3V4_plot$Genus))

pdf('Shared_genera_barplots/Plots/Milk_Oral_Baby_V3V4_barplot.pdf', width=7, height=4.5)
ggplot(DF_Milk_Oral_Baby_V3V4_plot, aes(x =Genus, y = count, fill = Sharing)) +
  geom_bar(stat = "identity", position = position_stack()) +
  facet_grid_paginate(~ `Sequencing method`, scales = "free_x", space = "free_x", switch = "y") +
  labs(y = "Number of pairs", fill = "Sharing") +
  scale_fill_manual(values = c("#CCCCCC","#0072B2"), labels = c("No", "Yes")) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14)) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        strip.text.x = element_text(size = 13),
        strip.background = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(size =13, angle =90, hjust = 1),
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=13), 
        axis.title.x = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(y = "Number of pairs", fill = "Present")
dev.off()

pdf('Shared_genera_barplots/Plots/Milk_Oral_Baby_16S23S_barplot.pdf', width=5, height=4.5)
ggplot(DFs_long_format$`Milk-Oral Baby`[DFs_long_format$`Milk-Oral Baby`$`Sequencing method`=="16S23S",], aes(x =Genus, y = count, fill = Sharing)) +
  geom_bar(stat = "identity", position = position_stack()) +
  facet_grid_paginate(~ `Sequencing method`, scales = "free_x", space = "free_x", switch = "y") +
  labs(y = "Number of pairs", fill = "Sharing") +
  scale_fill_manual(values = c("#CCCCCC","#0072B2"), labels = c("No", "Yes")) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14)) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        strip.text.x = element_text(size = 13),
        strip.background = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(size =13, angle =90, hjust = 1),
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=13), 
        axis.title.x = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(y = "Number of pairs", fill = "Present")
dev.off()

#******************
#B) Milk-Feces Baby
#******************
DF_Milk_Feces_Baby_V3V4_plot <- DFs_long_format$`Milk-Feces Baby`[DFs_long_format$`Milk-Feces Baby`$`Sequencing method`=="V3V4",]
DF_Milk_Feces_Baby_V3V4_plot$Genus <- factor(DF_Milk_Feces_Baby_V3V4_plot$Genus,
                                            levels = unique(DF_Milk_Feces_Baby_V3V4_plot$Genus))

pdf('Shared_genera_barplots/Plots/Milk_Feces_Baby_V3V4_barplot.pdf', width=5.5, height=4.5)
ggplot(DF_Milk_Feces_Baby_V3V4_plot, aes(x =Genus, y = count, fill = Sharing)) +
  geom_bar(stat = "identity", position = position_stack()) +
  facet_grid_paginate(~ `Sequencing method`, scales = "free_x", space = "free_x", switch = "y") +
  labs(y = "Number of pairs", fill = "Sharing") +
  scale_fill_manual(values = c("#CCCCCC","#0072B2"), labels = c("No", "Yes")) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14)) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        strip.text.x = element_text(size = 13),
        strip.background = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(size =13, angle =90, hjust = 1),
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=13), 
        axis.title.x = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(y = "Number of pairs", fill = "Present")
dev.off()

pdf('Shared_genera_barplots/Plots/Milk_Feces_Baby_16S23S_barplot.pdf', width=3.5, height=4.5)
ggplot(DFs_long_format$`Milk-Feces Baby`[DFs_long_format$`Milk-Feces Baby`$`Sequencing method`=="16S23S",], aes(x =Genus, y = count, fill = Sharing)) +
  geom_bar(stat = "identity", position = position_stack()) +
  facet_grid_paginate(~ `Sequencing method`, scales = "free_x", space = "free_x", switch = "y") +
  labs(y = "Number of pairs", fill = "Sharing") +
  scale_fill_manual(values = c("#CCCCCC","#0072B2"), labels = c("No", "Yes")) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14)) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        strip.text.x = element_text(size = 13),
        strip.background = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(size =13, angle =90, hjust = 1),
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=13), 
        axis.title.x = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(y = "Number of pairs", fill = "Present")
dev.off()

#********************
#C) Milk-Fecal Mother
#********************

pdf('Shared_genera_barplots/Plots/Milk_Feces_Mother_V3V4_barplot.pdf', width=3.5, height=4.5)
ggplot(DFs_long_format$`Milk-Feces Mother`[DFs_long_format$`Milk-Feces Mother`$`Sequencing method`=="V3V4",], aes(x =Genus, y = count, fill = Sharing)) +
  geom_bar(stat = "identity", position = position_stack()) +
  facet_grid_paginate(~ `Sequencing method`, scales = "free_x", space = "free_x", switch = "y") +
  labs(y = "Number of pairs", fill = "Sharing") +
  scale_fill_manual(values = c("#CCCCCC","#0072B2"), labels = c("No", "Yes")) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14)) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        strip.text.x = element_text(size = 13),
        strip.background = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(size =13, angle =90, hjust = 1),
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=13), 
        axis.title.x = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(y = "Number of pairs", fill = "Present")
dev.off()

pdf('Shared_genera_barplots/Plots/Milk_Feces_Mother_16S23S_barplot.pdf', width=2.5, height=4.5)
ggplot(DFs_long_format$`Milk-Feces Mother`[DFs_long_format$`Milk-Feces Mother`$`Sequencing method`=="16S23S",], aes(x =Genus, y = count, fill = Sharing)) +
  geom_bar(stat = "identity", position = position_stack()) +
  facet_grid_paginate(~ `Sequencing method`, scales = "free_x", space = "free_x", switch = "y") +
  labs(y = "Number of pairs", fill = "Sharing") +
  scale_fill_manual(values = c("#CCCCCC","#0072B2"), labels = c("No", "Yes")) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14)) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        strip.text.x = element_text(size = 13),
        strip.background = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(size =13, angle =90, hjust = 1),
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=13), 
        axis.title.x = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(y = "Number of pairs", fill = "Present")
dev.off()
