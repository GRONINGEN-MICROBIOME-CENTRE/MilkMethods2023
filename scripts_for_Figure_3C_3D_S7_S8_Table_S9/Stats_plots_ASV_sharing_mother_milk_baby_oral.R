################################################################################
##### ASVs transmission between maternal milk and infant oral samples:  statistical analysis and plotting
### Author(s): Asier Fernández / Johanne E. Spreckels
### Last updated: 19th September, 2023
################################################################################

#****************
# Load libraries
#****************
library(dplyr)
library(tidyr)
library(ggplot2)

# Set working directory
setwd("~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/")

# Read results files
final_PSE_data <- read.delim("Result_tables/PSE_table_maternal_milk_infant_feces.txt", header = T, check.names = F) 
cont_tables_V3V4 <- readRDS("16S/Mother_Milk-Baby_Feces/Contingency_tables_16S.rds")
cont_tables_16S23S <- readRDS("16S23S/Mother_Milk-Baby_Feces/Contingency_tables_16S23S.rds")

###############################################################
# Generate PSE results table and test statistical significance 
###############################################################

#******************************
# Test PSE differences in related vs unrelated pairs
#******************************
# Reformat the data table 
final_PSE_data$Pair_ID <- c(rownames(final_PSE_data))
final_PSE_data_long <- final_PSE_data %>%
  pivot_longer(cols = -c(Pair_ID, Relatedness), names_to = "Genus_Tech", values_to = "Value") %>%
  separate(Genus_Tech, into = c("Genus", "Technology"), sep = "_")

final_PSE_data_summary <- final_PSE_data_long %>%
  dplyr::group_by(Genus, Technology, Value, Relatedness) %>%
  dplyr::summarize(count = n()) %>%
  tidyr::pivot_wider(names_from = Value, values_from = count, values_fill = 0)

final_PSE_data_summary$all_no_PSE <- rowSums(data.frame(final_PSE_data_summary$no_PSE, final_PSE_data_summary$`NA`))
final_PSE_data_summary$Genus_Tech <- paste(final_PSE_data_summary$Genus, final_PSE_data_summary$Technology, sep = "_")


# Generate contingency tables and estimate p-values for each genus
Fischer_tests <- data.frame(matrix(nrow=length(unique(final_PSE_data_summary$Genus_Tech)), ncol=1))
colnames(Fischer_tests) <- "p_value"
rownames(Fischer_tests) <- unique(final_PSE_data_summary$Genus_Tech)

for (i in 1:length(rownames(Fischer_tests))) {
  genus_data <- final_PSE_data_summary[final_PSE_data_summary$Genus_Tech==rownames(Fischer_tests)[i],]
  cont_table <- data.frame(rbind(c(genus_data[genus_data$Relatedness=="Related",7], #7
                                   genus_data[genus_data$Relatedness=="Related",4]),
                                 c(genus_data[genus_data$Relatedness=="Unrelated",7], #7
                                   genus_data[genus_data$Relatedness=="Unrelated",4])))
  cont_table <- matrix(unlist(cont_table), 2)
  colnames(cont_table) <- c("no_PSE", "PSE")
  rownames(cont_table) <- c("Related", "Unrelated")
  test <- fisher.test(cont_table)
  Fischer_tests[i,1] <- test$p.value  
}

#Estimate the FDR per sequencing method
Fischer_tests$FDR <- NA

for (i in unique(final_PSE_data_summary$Technology)) {
  technology <- i #select the technology
  row_numbers <- grep(technology,rownames(Fischer_tests)) #select the row numbers in the dataframe
  p_values <- Fischer_tests[row_numbers,"p_value"]
  FDR_values <- p.adjust(p_values, method = "BH")
  Fischer_tests$FDR[row_numbers] <- FDR_values
}

final_PSE_data_summary <- merge(final_PSE_data_summary, Fischer_tests, by.x="Genus_Tech", by.y=0)

# Reformat results table in long format
final_PSE_data_summary_long <- final_PSE_data_summary %>%
  pivot_longer(cols = c("PSE", "no_PSE", "NA"), names_to = "Value", values_to = "count")

# Add the total counts, percentages and propotion for each group to the dataframe (Relatedness, Genus, Technology)
total_count <- final_PSE_data_summary_long %>% 
  dplyr::group_by(Relatedness, Genus, Technology) %>% 
  dplyr::summarize(total = sum(count))

final_PSE_data_summary_long <- dplyr::left_join(final_PSE_data_summary_long, total_count)
final_PSE_data_summary_long$percent <- paste0(round(final_PSE_data_summary_long$count / final_PSE_data_summary_long$total * 100), "%")
final_PSE_data_summary_long$proportion <- final_PSE_data_summary_long$count / final_PSE_data_summary_long$total * 100
final_PSE_data_summary_long$all_no_PSE <- NULL

#******************************
# Test PSE differences in unrelated pairs between sequencing methods
#******************************
# Subset the data frame with previous results per technology
results_V3V4 <- subset(final_PSE_data_summary , Technology == "V3V4")
results_16S23S <- subset(final_PSE_data_summary , Technology == "16S23S")

# Create an empty dataframe to store results and list the genera
genus_list <- unique(final_PSE_data_summary$Genus)
sequencing_PSE_results <- data.frame(Genus= character(), p_value = numeric(), Relatedness=character(), stringsAsFactors = FALSE)

# Generate the contingency tables and estimate the p-values (Fischer test)
for (Genus in genus_list) {
  # create contingency table for each Genus
  Cont_table_related <- rbind(results_V3V4[results_V3V4$Genus == Genus & results_V3V4$Relatedness == "Related", c("PSE", "all_no_PSE")],
                              results_16S23S[results_16S23S$Genus == Genus & results_16S23S$Relatedness == "Related", c("PSE", "all_no_PSE")])
  Cont_table_unrelated <- rbind(results_V3V4[results_V3V4$Genus == Genus & results_V3V4$Relatedness == "Unrelated", c("PSE", "all_no_PSE")],
                                results_16S23S[results_16S23S$Genus == Genus & results_16S23S$Relatedness == "Unrelated", c("PSE", "all_no_PSE")])
  rownames(Cont_table_related) <- c("V3V4", "16S23S")
  rownames(Cont_table_unrelated) <- c("V3V4", "16S23S")
  # Perform Fisher's exact test
  fisher_test_related <- fisher.test(Cont_table_related)
  fisher_test_unrelated <- fisher.test(Cont_table_unrelated)
  # Add results to the data frame
  sequencing_PSE_results <- rbind(sequencing_PSE_results, 
                                  data.frame(Genus= Genus, 
                                             p_value = fisher_test_related$p.value, 
                                             Relatedness = "Related"))
  sequencing_PSE_results <- rbind(sequencing_PSE_results, 
                                  data.frame(Genus= Genus, 
                                             p_value = fisher_test_unrelated$p.value, 
                                             Relatedness = "Unrelated"))
}

# Select only "Unrelated" comparisons and estimate FDR (BH) correction
sequencing_PSE_results <- sequencing_PSE_results[sequencing_PSE_results$Relatedness == "Unrelated",]
sequencing_PSE_results$FDR <- p.adjust(sequencing_PSE_results$p_value, method = "BH")
sequencing_PSE_results <- sequencing_PSE_results[,c("Genus", "Relatedness", "p_value", "FDR")]

#******************************
# Generate the plots
#******************************
# Generate the stacked bar plot for A) V3-V4 and B) 16S23S
#A) V3-V4
final_PSE_data_summary_long_V3V4 <- final_PSE_data_summary_long[final_PSE_data_summary_long$Technology == "V3V4",] #alternative dataframe for plotting without V1V9
final_PSE_data_summary_long_V3V4$Value <- factor(final_PSE_data_summary_long_V3V4$Value, levels = c("NA", "no_PSE", "PSE"))

pdf('Plots/PSE_relatedness_V3V4_mother_milk_baby_oral.pdf', width=9.2, height=4)
ggplot(final_PSE_data_summary_long_V3V4, aes(x = Relatedness, y = proportion, fill = Value)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ Genus, switch = "y") +
  labs(y = "Number of pairs", fill = "Sharing") +
  scale_fill_manual(values = c("#CCCCCC","#E69F00","#0072B2"), labels = c("Genus not present","No", "Yes")) +
  scale_y_continuous(breaks=c(0,25, 50, 75, 100), expand=c(0.07,0),limits=c(0,107)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        strip.text.x = element_text(size = 13),
        strip.background = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(size =13, angle =90, hjust = 1),
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=13), 
        axis.title.x = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  # Add percentage labels for the PSE group
  geom_label(aes(label = percent),
             fill = "#FFFFFF",
             data = subset(final_PSE_data_summary_long_V3V4, Value == "PSE"),
             position = position_fill(vjust = 50), size = 4.7) +
  labs(y = "Percentage of pairs with PSEs", fill = "Sharing")
dev.off()

#B) 16S-23S
final_PSE_data_summary_long_16S23S <- final_PSE_data_summary_long[final_PSE_data_summary_long$Technology == "16S23S",] #alternative dataframe for plotting without V1V9
final_PSE_data_summary_long_16S23S$Value <- factor(final_PSE_data_summary_long_16S23S$Value, levels = c("NA", "no_PSE", "PSE"))

pdf('Plots/PSE_relatedness_16S23S_mother_milk_baby_oral.pdf', width=9.2, height=4)
ggplot(final_PSE_data_summary_long_16S23S, aes(x = Relatedness, y = proportion, fill = Value)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ Genus, switch = "y") +
  labs(y = "Number of pairs", fill = "Sharing") +
  scale_fill_manual(values = c("#CCCCCC","#E69F00","#0072B2"), labels = c("Genus not present","No", "Yes")) +
  scale_y_continuous(breaks=c(0,25, 50, 75, 100), expand=c(0.07,0),limits=c(0,107)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        strip.text.x = element_text(size = 13),
        strip.background = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(size =13, angle =90, hjust = 1),
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=13), 
        axis.title.x = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  # Add percentage labels for the PSE group
  geom_label(aes(label = percent),
             fill = "#FFFFFF",
             data = subset(final_PSE_data_summary_long_16S23S, Value == "PSE"),
             position = position_fill(vjust = 50), size = 4.7) +
  labs(y = "Percentage of pairs with PSEs", fill = "Sharing")
dev.off()

########################
# Save results
########################
write.table(final_PSE_data_summary_long, file ="~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/Result_tables/Mother_milk_Baby_Oral_final_PSE_results.txt", sep="\t")
write.table(sequencing_PSE_results, file ="~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/Result_tables/Mother_milk_Baby_Oral_final_PSE_results_Seq_Technology.txt", sep="\t")

