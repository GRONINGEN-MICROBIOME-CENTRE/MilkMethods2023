### ========================== SCRIPT 3 - COMPARISON OF BACTERIAL ENRICHMENT METHODS ========================== ###

### SCRIPT:          SCRIPT 3 - COMPARISON OF BACTERIAL ENRICHMENT METHODS
##
## DESCRIPTION:      Script to generate/run
##                    - Table S5
##
## AUTHORS:          Johanne Spreckels
##
## NOTES:            To use this script, set correct paths to the following folders/files:
##                    - MY_PATH_TO_ENRICHMENT_METHOD_COMPARISON_DATA/Data_enrichment_method_comparison_n39.txt
##                    - MY_PATH_TO_ENRICHMENT_METHOD_COMPARISON_RESULTS/
##
## DATE OF CREATION: Script generated for publication in September 2023

### CONTENTS OF THIS FILE
## 0. IMPORT DATA
## 1. SUBSET SEQUENCING DATA BY SAMPLE TYPE
## 2. SUMMARY STATISTICS TABLE FOR MILK SAMPLES
## 3. ADD STATISITCS TO SUMMARY STATISTICS TABLE FOR MILK SAMPLES
## 4. SAVE TABLE S5


### ===== 0. IMPORT DATA ===== ###

## import Data_enrichment_method_comparison_n39.txt
data <- read.table("MY_PATH_TO_ENRICHMENT_METHOD_COMPARISON_DATA/Data_enrichment_method_comparison_n39.txt", header=T, sep="\t", stringsAsFactors=T)


### ===== 1. SUBSET SEQUENCING DATA BY SAMPLE TYPE ===== ###

## set sample_ids as rownames
rownames(data) <- data$sample_ID

## subset by sample type
mock <- data[data$sample_type=="Mock",] #3 samples
milk <- data[data$sample_type=="Milk",] #27 samples
negctrl <- data[data$sample_type=="Negative_control",] #9 samples


### ===== 2. SUMMARY STATISTICS TABLE FOR MILK SAMPLES ===== ###

## split milk data frame by enrichment method
milk_NE <- milk[milk$enrichment_method=="No_enrichment",]
milk_benz <- milk[milk$enrichment_method=="Osmotic_lysis_and_benzonase",]
milk_PMA <- milk[milk$enrichment_method=="Osmotic_lysis_and_PMA",]

## function to retrieve median (range)
getSumStats_milk <- function(inputdata){
  SumStats_milk <- c(nrow(inputdata),
                     paste0(median(inputdata$DNA_yield_ng), " (", min(inputdata$DNA_yield_ng), " - ", max(inputdata$DNA_yield_ng), ")"),
                     paste0(nrow(inputdata[inputdata$mgx_seq_successful=="yes" & !is.na(inputdata$mgx_seq_successful),]), "/", nrow(inputdata)),
                     paste0(median(inputdata[inputdata$mgx_seq_successful=="yes" & !is.na(inputdata$mgx_seq_successful),"n_reads_total"]), " (", min(inputdata[inputdata$mgx_seq_successful=="yes" & !is.na(inputdata$mgx_seq_successful),"n_reads_total"]), " - ", max(inputdata[inputdata$mgx_seq_successful=="yes" & !is.na(inputdata$mgx_seq_successful),"n_reads_total"]), ")"),
                     paste0(median(inputdata[inputdata$mgx_seq_successful=="yes" & !is.na(inputdata$mgx_seq_successful),"n_reads_microbial"]), " (", min(inputdata[inputdata$mgx_seq_successful=="yes" & !is.na(inputdata$mgx_seq_successful),"n_reads_microbial"]), " - ", max(inputdata[inputdata$mgx_seq_successful=="yes" & !is.na(inputdata$mgx_seq_successful),"n_reads_microbial"]), ")"),
                     paste0(median(inputdata[inputdata$mgx_seq_successful=="yes" & !is.na(inputdata$mgx_seq_successful),"perc_reads_microbial"])*100, " (", min(inputdata[inputdata$mgx_seq_successful=="yes" & !is.na(inputdata$mgx_seq_successful),"perc_reads_microbial"])*100, " - ", max(inputdata[inputdata$mgx_seq_successful=="yes" & !is.na(inputdata$mgx_seq_successful),"perc_reads_microbial"])*100, ")"),
                     paste0(median(inputdata[inputdata$mgx_seq_successful=="yes" & !is.na(inputdata$mgx_seq_successful),"n_reads_human"]), " (", min(inputdata[inputdata$mgx_seq_successful=="yes" & !is.na(inputdata$mgx_seq_successful),"n_reads_human"]), " - ", max(inputdata[inputdata$mgx_seq_successful=="yes" & !is.na(inputdata$mgx_seq_successful),"n_reads_human"]), ")"),
                     paste0(median(inputdata[inputdata$mgx_seq_successful=="yes" & !is.na(inputdata$mgx_seq_successful),"perc_reads_human"])*100, " (", min(inputdata[inputdata$mgx_seq_successful=="yes" & !is.na(inputdata$mgx_seq_successful),"perc_reads_human"])*100, " - ", max(inputdata[inputdata$mgx_seq_successful=="yes" & !is.na(inputdata$mgx_seq_successful),"perc_reads_human"])*100, ")"))
  
  names(SumStats_milk) <- c("n (samples)", "DNA yield (ng)", "n/N (successfully sequenced/total number of samples)", "n (total reads)", "n (microbial reads)", "% (microbial reads)", "n (human reads)", "% (human reads)")
  
  return(SumStats_milk)
}

SumStats_milk <- as.data.frame(cbind(getSumStats_milk(milk_NE),
                                     getSumStats_milk(milk_benz),
                                     getSumStats_milk(milk_PMA)))
colnames(SumStats_milk) <- c("No_enrichment", "Osmotic_lysis_and_benzonase", "Osmotic_lysis_and_PMA")
SumStats_milk$Info <- rownames(SumStats_milk)
rownames(SumStats_milk) <- 1:nrow(SumStats_milk)
SumStats_milk <- SumStats_milk[,c(4,1:3)]
SumStats_milk


### ===== 3. ADD STATISITCS TO SUMMARY STATISTICS TABLE FOR MILK SAMPLES ===== ###

## 4.1 Compare DNA_yield_ng between all 3 enrichment methods (Kruskal-Wallis test)
kw <- kruskal.test(DNA_yield_ng ~ enrichment_method, data=milk) #p-value = 6.075e-06


## 4.2 Compare number of samples successfully sequenced, only comparing No_enrichment and Osmotic_lysis_and_PMA (Fisher's test)
# create contingency table for use in Fisher's test
milk_matrix <- data.frame("mgx_successful" = c(9, 5),
                          "mgx_failed" = c(0, 4),
                          row.names = c("No_enrichment", "Osmotic_lysis_and_PMA"))
ft <- fisher.test(milk_matrix) #p-value = 0.08235


## 4.3 Compare number and percentage of microbial/human reads, only comparing No_enrichment and Osmotic_lysis_and_PMA (Wilcoxon test)
performWilcoxonTextNEvsPMA <- function(inputdata, columns_of_interest){
  name_x <- c()
  name_y <- c()
  n_NE <- c()
  n_PMA <- c()
  p <- c()
  statistic <- c()
  
  for (i in columns_of_interest){
    #perform paired Wilcoxon test
    WT <- wilcox.test(x=as.numeric(inputdata[inputdata$enrichment_method=="No_enrichment" & !is.na(inputdata[,i]), i]),
                      y=as.numeric(inputdata[inputdata$enrichment_method=="Osmotic_lysis_and_PMA" & !is.na(inputdata[,i]), i]),
                      alternative="two.sided", paired=F, exact=T)
    
    #save data
    name_x <- c(name_x, "enrichment_method")
    name_y <- c(name_y, colnames(inputdata)[i])
    n_NE <- c(n_NE, nrow(inputdata[inputdata$enrichment_method=="No_enrichment" & !is.na(inputdata[,i]),]))
    n_PMA <- c(n_PMA, nrow(inputdata[inputdata$enrichment_method=="Osmotic_lysis_and_PMA" & !is.na(inputdata[,i]),]))
    statistic <- c(statistic, "wilcoxon_test")
    p <- c(p, WT$p.value)
  }
  
  #save all results in 1 data frame to export
  results <- data.frame(name_x=name_x, name_y=name_y, n_NE=n_NE, n_PMA=n_PMA, statistic=statistic, p_value=p)
  
  #calculate FDR
  # results$BH_adjusted_p_value <- p.adjust(results$p_value, method="BH")
  
  return(results)
}

wt <- performWilcoxonTextNEvsPMA(milk, columns_of_interest = c(18,20,22,19,21))


## 4.4 Add all statistic results to SumStats_milk table
SumStats_milk$p <- c(NA, kw$p.value, ft$p.value, wt$p_value)


## 4.5 Add multiple testing correction (Benjamini-Hochberg)
SumStats_milk$FDR <- p.adjust(SumStats_milk$p, method="BH")
SumStats_milk


### ===== 4. SAVE TABLE S5 ===== ###

write.table(SumStats_milk, file="MY_PATH_TO_ENRICHMENT_METHOD_COMPARISON_RESULTS/Table_S5.txt", sep="\t", row.names=F, quote=F)




