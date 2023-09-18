### ========================== SCRIPT 6 - COMPARISON OF SEQUENCING METHODS - PILOT SAMPLES ========================== ###

### SCRIPT:          SCRIPT 6 - COMPARISON OF SEQUENCING METHODS - PILOT SAMPLES
##
## DESCRIPTION:      Script to generate/run
##                    - Figure 3B
##                    - Figure S4
##                    - Figure S5
##                    - Figure S6
##                    - Table S8
##
## AUTHORS:          Johanne Spreckels
##
## NOTES:            To use this script, set correct paths to the following folders/files:
##                    - MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_PILOT_DATA/
##                    - MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_PILOT_RESULTS/
##
## DATE OF CREATION: Script generated for publication in September 2023

### CONTENTS OF THIS FILE
## 0. IMPORT DATA
## 1. SELECT SAMPLES OF INTEREST 
## 2. PLOT TOTAL BACTERIAL LOAD BY SAMPLE TYPE (FIGURE S4)
## 3. CALCULATE RELATIVE ABUNDANCES (BASED ON ASVs)
## 4. SUBSET SEQUENCING DATA BY TAXONOMIC LEVEL
## 5. MERGE RELATIVE ABUNDANCES BACK WITH METADATA
## 6. AITCHISON DISTANCES, ADONIS AND PCoA PLOTS (FIGURE 3B)
##     6.1 CLR-TRANSFORM RELATIVE ABUNDANCES
##     6.2 CALCULATE EUCLIDIAN DISTANCES, RUN ADONIS AND MAKES PCoAs (FIGURE 3B)
## 7. RELATIVE ABUNDANCE PLOTS BY SAMPLE TYPE AND SEQUENCING METHOD (FIGURE S5)
##     7.1 PREPARE DATA FOR PLOTTING, INCL. FILTERING ON MOST PREVALENT AND MOST ABUNDANT BACTERIA
##     7.2 CREATE RELATIVE ABUNDANCE PLOT ON GENUS LEVEL (FIGURE S5)
## 8. CHECK FOR MISSING MOCK BACTERIA IN LIFELINES NEXT SAMPLES (TABLE S8)
## 9. CHECK GENUS RICHNESS BY SAMPLE TYPE AND SEQUENCING METHOD (FIGURE S6)
##     9.1 PREPARE DATA FRAME
##     9.2 CALCULATE CLASSIFIED GENERA RICHNESS PER SAMPLE TYPE AND SEQUENCING METHOD
##     9.3 PLOT GENERA RICHNESS PER SAMPLE TYPE AND SEQUENCING METHOD
##     9.4 STATISTICS FOR GENERA RICHNESS PER SAMPLE TYPE AND SEQUENCING METHOD


### ===== 0. IMPORT DATA ===== ###

## set working directory
setwd("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_PILOT_DATA/")

## import Data_sequencing_method_comparison_n136.txt
df <- read.table("Data_sequencing_method_comparison_n136.txt", header=T, sep="\t")


### ===== 1. SELECT SAMPLES OF INTEREST ===== ###

## set sample_ids as rownames
rownames(df) <- df$alias

## select NEXT pilot samples (n=112)
data <- df[df$sample_set=="pilot_sample",]


### ===== 2. PLOT TOTAL BACTERIAL LOAD BY SAMPLE TYPE (FIGURE S4) ===== ###

library(ggplot2)

## select qpcr data
qpcrdata <- data[,c(4,16,19)]
qpcrdata <- qpcrdata[!duplicated(qpcrdata),] #ensure every sample is only included 1x in the dataset

## sort sample_origin_type levels before plotting
qpcrdata$sample_origin_type <- factor(qpcrdata$sample_origin_type, levels = c("Human_milk", "Infant_oral_cavity", "Infant_faeces", "Maternal_faeces"))

## plot log10-transformed 16S rRNA gene copies per 1 ng sample DNA
qpcrplot <- ggplot(qpcrdata, aes(x=sample_origin_type, y=qPCR_16S_rRNA_gene_copies_per_1ng_DNA_log10))+
  geom_jitter(aes(color=sample_origin_type, size=2), alpha=0.5)+
  geom_boxplot(alpha=0)+
  theme_bw()+
  theme(legend.position="none",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())+
  scale_y_continuous(limits = c(2.5,12.5))+
  scale_color_manual(values = c("#E69F00", "#CC79A7", "#009E73", "#0072B2"))+
  labs(x="", y="log10 (16S rRNA gene copies / 1 ng DNA)")
qpcrplot
# ggsave("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_PILOT_RESULTS/Figure_S4.pdf", qpcrplot, dpi=500, width=15, height=15, units="cm", useDingbats=F)


## Kruskal-Wallis test
KT <- kruskal.test(qPCR_16S_rRNA_gene_copies_per_1ng_DNA_log10 ~ sample_origin_type, data=qpcrdata)
#p-value = 1.841e-09

## run Dunn's posthoc test
library(FSA)

# Dunn's posthoc test for DNA yield
DT1 <- dunnTest(qPCR_16S_rRNA_gene_copies_per_1ng_DNA_log10 ~ sample_origin_type, data=qpcrdata, method="bh")
dunnTestRes <- as.data.frame(DT1$res)
# write.table(dunnTestRes, file="MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Figure_S4_Dunns_test_results_total_bact_load_by_sample_type.txt", sep="\t", row.names=F, quote=F)


### ===== 3. CALCULATE RELATIVE ABUNDANCES (BASED ON ASVs) ===== ###

## save only columns with absolute bacterial abundances
bact <- data[,22:ncol(data)]

## calculate relative bacterial abundances
data_relab <- (bact/rowSums(bact))
table(rowSums(data_relab), useNA="ifany") #1 sample has rowSums=NA (the maternal feces sample from pair5 failed 16S sequencing)

## for the negative control that failed sequencing, set relative abundances to 0
data_relab[rownames(data_relab)=="22Jun148-DL052",] <- 0
table(rowSums(data_relab), useNA="ifany") #now the 1 sample has rowSums=0, all others have rowSums=1 (100%)

## merge relative abundances back with metadata
data2 <- merge(data[,1:21], data_relab, by="row.names")
data2 <- data2[,-1] #remove column Row.names


### ===== 4. SUBSET SEQUENCING DATA BY TAXONOMIC LEVEL ===== ###

## function to subset data frame by taxonomic level (genus levels)
subset.by.taxlevel <- function(inputdata, pattern){ 
  colnames(inputdata) <- gsub(pattern, "", colnames(inputdata)) #shorten colnames after selected taxonomic level
  taxon <- t(rowsum(t(inputdata), group = colnames(inputdata), na.rm = T)) #combine reads from columns with the same name in one column
  taxon <- taxon[,colSums(taxon)>0] #remove columns where colSums==0
  return(taxon)
}

## subset data by taxonomic level
data_genus <- subset.by.taxlevel(inputdata = data2[,22:ncol(data2)], pattern = "_s__.*") #112x280 genera


### ===== 5. MERGE RELATIVE ABUNDANCES BACK WITH METADATA ===== ###

# genus
df_genus <- merge(data2[,c(4,15,16,20:21)], data_genus, by="row.names")
df_genus <- df_genus[,-1]


### ===== 6. AITCHISON DISTANCES, ADONIS AND PCoA PLOTS (FIGURE 3B) ===== ###

## exclude bacteria that are unclassified on genus level
df_genus2 <- df_genus[,-grep("g__NA", colnames(df_genus))]

## set sample names as rownames
rownames(df_genus2) <- c(paste0(df_genus2$sequencing_method, "_", df_genus2$subjectId))


### 6.1 CLR-TRANSFORM RELATIVE ABUNDANCES

# function for clr transformation
do_clr_externalWeighting = function(interest_matrix, core_matrix) {
  if(any(interest_matrix==0)) interest_matrix = interest_matrix + min(interest_matrix[interest_matrix>0])/2
  if(any(core_matrix==0)) core_matrix = core_matrix + min(core_matrix[core_matrix>0])/2
  #estimate weighting parameter
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  
  Gmean_core = apply(core_matrix, 1, gm_mean)
  
  #do transformation
  data_prepared = cbind(Gmean_core,interest_matrix)
  data_transformed = t(apply(data_prepared,1,function(x){
    log(x / x[1])[-1]
  }))
  colnames(data_transformed) = colnames(data_transformed)
  rownames(data_transformed) = rownames(data_transformed)
  data_transformed
}

df_genus2_clr <- df_genus2
df_genus2_clr[,6:ncol(df_genus2_clr)] <- as.data.frame(do_clr_externalWeighting(df_genus2_clr[,6:ncol(df_genus2_clr)], df_genus2_clr[,6:ncol(df_genus2_clr)]))


### 6.2 CALCULATE EUCLIDIAN DISTANCES, RUN ADONIS AND MAKES PCoAs (FIGURE 3B)

library(vegan)

## function for adonis and pcoa plots with Euclidian distances - all by sequencing method, sample type and pair ID
run.adonis.and.plot.pcoas.Ait <- function(inputdata){
  
  ## AITCHISON DISTANCE MATRIX / + METADATA DATA FRAME
  print("Creating Aitchison distance matrix and data frame")
  ait <- vegdist(inputdata, method="euclidian")
  aitm <- cmdscale(ait, k=3) #create distance matrix for first 3 components
  aitd <- data.frame(aitm) #convert matrix to data frame
  
  print("Calculating variance explained")
  beta.cmd <- cmdscale(as.matrix(ait), k=2, eig=T)
  PCoA1 <- beta.cmd$eig[1]/sum(beta.cmd$eig)*100
  PCoA2 <- beta.cmd$eig[2]/sum(beta.cmd$eig)*100
  print(paste0("PCoA1 explains ", PCoA1, "% of the variance"))
  print(paste0("PCoA2 explains ", PCoA2, "% of the variance"))
  
  ## add metadata columns to Aitchison distance data frame (ensure that rownames indicate sample info)
  # add column for sequencing method
  aitd[grep("16S_V3V4", rownames(aitd)),4] <- "16S"
  aitd[grep("16S-ITS-23S", rownames(aitd)),4] <- "16S-ITS-23S"
  
  # add column for sample_type
  aitd[grep("Human_milk", rownames(aitd)),5] <- "A_Human_milk"
  aitd[grep("Maternal_faeces", rownames(aitd)),5] <- "D_Maternal_faeces"
  aitd[grep("Infant_oral_cavity", rownames(aitd)),5] <- "B_Infant_oral_cavity"
  aitd[grep("Infant_faeces", rownames(aitd)),5] <- "C_Infant_faeces"
  
  # add column for pair_ID
  aitd[grep("Pair1_", rownames(aitd)),6] <- "Pair1"
  aitd[grep("Pair2_", rownames(aitd)),6] <- "Pair2"
  aitd[grep("Pair3_", rownames(aitd)),6] <- "Pair3"
  aitd[grep("Pair4_", rownames(aitd)),6] <- "Pair4"
  aitd[grep("Pair5_", rownames(aitd)),6] <- "Pair5"
  aitd[grep("Pair6_", rownames(aitd)),6] <- "Pair6"
  aitd[grep("Pair7_", rownames(aitd)),6] <- "Pair7"
  aitd[grep("Pair8_", rownames(aitd)),6] <- "Pair8"
  aitd[grep("Pair9_", rownames(aitd)),6] <- "Pair9"
  aitd[grep("Pair10_", rownames(aitd)),6] <- "Pair10"
  aitd[grep("Pair11_", rownames(aitd)),6] <- "Pair11"
  aitd[grep("Pair12_", rownames(aitd)),6] <- "Pair12"
  aitd[grep("Pair13_", rownames(aitd)),6] <- "Pair13"
  aitd[grep("Pair14_", rownames(aitd)),6] <- "Pair14"

  ## fix colnames and data structure of Aitchison distance + metadata data frame
  colnames(aitd)[4:6] <- c("sequencing_method", "sample_type", "pair_ID")
  for (i in 4:6){aitd[,i] <- as.factor(as.character(aitd[,i]))}
  
  
  ## RUN ADONIS
  print("Running Adonis")
  print("Running Adonis for sequencing method")
  adonisMethod <- adonis(ait ~ aitd$sequencing_method)
  adonisMethodRes <- adonisMethod$aov.tab
  
  print("Running Adonis for sample type")
  adonisSampleType <- adonis(ait ~ aitd$sample_type)
  adonisSampleTypeRes <- adonisSampleType$aov.tab
  
  print("Running Adonis for pair ID")
  adonisPair <- adonis(ait ~ aitd$pair_ID)
  adonisPairRes <- adonisPair$aov.tab
  
  adonisResults <- data.frame(x = c("sequencing_method", "sample_type", "pair_ID"),
                              R2 = c(adonisMethodRes$R2[1], adonisSampleTypeRes$R2[1], adonisPairRes$R2[1]),
                              p = c(adonisMethodRes$`Pr(>F)`[1], adonisSampleTypeRes$`Pr(>F)`[1], adonisPairRes$`Pr(>F)`[1]))
  adonisResults$BH_adj_p <- p.adjust(adonisResults$p, method = "BH")
  
  # write.table(adonisResults, paste0("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_PILOT_RESULTS/NEXT_pilot_all_samples_Aitchison_adonis_results_genus.txt"), sep="\t", row.names=F, quote=F)

  
  ## PCoA PLOTS
  print("Creating PCoA plots")
  
  # PCoA plot for X1 (PC1) and X2 (PC2), colored by sample_type and shaped by sequencing_method
  pcoa.sample.type.PC1.PC2 <- ggplot(aitd, aes(x=X1, y=X2, color=sample_type, shape=sequencing_method))+
    geom_point(alpha=0.5, aes(size=2))+
    stat_ellipse(aes(group=sample_type))+
    labs(x=paste0("PC1 (", signif(PCoA1,5), "%)"), y=paste0("PC2 (", signif(PCoA2,5), "%)"), color="Sample type", shape="Sequencing method", title="Aitchison PCoA PC1-PC2 - Genus level")+
    theme_bw()+
    theme(panel.grid=element_blank())+
    scale_color_manual(values=c("#E69F00","#CC79A7","#009E73","#0072B2"))
  pcoa.sample.type.PC1.PC2
  # ggsave(paste0("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_PILOT_RESULTS/Figure_3B.pdf"), pcoa.sample.type.PC1.PC2, dpi=500, width=15, height=10, units="cm", useDingbats=F)

}

# all samples (exclude one maternal sample which had read count=0)
run.adonis.and.plot.pcoas.Ait(df_genus2_clr[df_genus2_clr$n_reads_total>0, 6:ncol(df_genus2_clr)])


## function for adonis and pcoa plots for one sample type at a time
run.adonis.and.plot.pcoas.single <- function(inputdata, samples){
  
  ## AITCHISON DISTANCE MATRIX / + METADATA DATA FRAME
  print("Creating Aitchison distance matrix and data frame")
  ait <- vegdist(inputdata, method="euclidian")
  aitm <- cmdscale(ait, k=3) #create distance matrix for first 3 components
  aitd <- data.frame(aitm) #convert matrix to data frame
  
  ## add metadata columns to Aitchison distance data frame (ensure that rownames indicate sample info)
  # add column for sequencing method
  aitd[grep("16S_V3V4", rownames(aitd)),4] <- "16S"
  aitd[grep("16S-ITS-23S", rownames(aitd)),4] <- "16S-ITS-23S"
  
  # add column for sample_type
  aitd[grep("Human_milk", rownames(aitd)),5] <- "Human_milk"
  aitd[grep("Matenal_faeces", rownames(aitd)),5] <- "Maternal_faeces"
  aitd[grep("Infant_oral_cavity", rownames(aitd)),5] <- "Infant_oral_cavity"
  aitd[grep("Infant_faeces", rownames(aitd)),5] <- "Infant_faeces"
  
  # add column for Pair_ID
  aitd[grep("Pair1_", rownames(aitd)),6] <- "Pair1"
  aitd[grep("Pair2_", rownames(aitd)),6] <- "Pair2"
  aitd[grep("Pair3_", rownames(aitd)),6] <- "Pair3"
  aitd[grep("Pair4_", rownames(aitd)),6] <- "Pair4"
  aitd[grep("Pair5_", rownames(aitd)),6] <- "Pair5"
  aitd[grep("Pair6_", rownames(aitd)),6] <- "Pair6"
  aitd[grep("Pair7_", rownames(aitd)),6] <- "Pair7"
  aitd[grep("Pair8_", rownames(aitd)),6] <- "Pair8"
  aitd[grep("Pair9_", rownames(aitd)),6] <- "Pair9"
  aitd[grep("Pair10_", rownames(aitd)),6] <- "Pair10"
  aitd[grep("Pair11_", rownames(aitd)),6] <- "Pair11"
  aitd[grep("Pair12_", rownames(aitd)),6] <- "Pair12"
  aitd[grep("Pair13_", rownames(aitd)),6] <- "Pair13"
  aitd[grep("Pair14_", rownames(aitd)),6] <- "Pair14"
  
  ## fix colnames and data structure of Aitchison distance + metadata data frame
  colnames(aitd)[4:6] <- c("sequencing_method", "sample_type", "pair_ID")
  for (i in 4:6){aitd[,i] <- as.factor(as.character(aitd[,i]))}
  
  
  ## RUN ADONIS
  print("Running Adonis")
  print("Running Adonis for sequencing method")
  adonisMethod <- adonis(ait ~ aitd$sequencing_method)
  adonisMethodRes <- adonisMethod$aov.tab
  
  print("Running Adonis for pair ID")
  adonisPair <- adonis(ait ~ aitd$pair_ID)
  adonisPairRes <- adonisPair$aov.tab
  
  adonisResults <- data.frame(x = c("sequencing_method", "pair_ID"),
                              R2 = c(adonisMethodRes$R2[1], adonisPairRes$R2[1]),
                              p = c(adonisMethodRes$`Pr(>F)`[1],  adonisPairRes$`Pr(>F)`[1]))
  adonisResults$BH_adj_p <- p.adjust(adonisResults$p, method = "BH")
  
  write.table(adonisResults, paste0("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_PILOT_RESULTS/", samples, "_Aitchison_adonis_results_genus.txt"), sep="\t", row.names=F, quote=F)


  # PCoA plot for X1 (PC1) and X2 (PC2), colored by sequencing_method
  print("Plot PCoA")
  pcoa.seq.method.PC1.PC2 <- ggplot(aitd, aes(x=X1, y=X2, color=sequencing_method))+
    geom_point(alpha=0.5)+
    stat_ellipse(aes(group=sequencing_method))+
    labs(x="PC1", y="PC2", color="Sequencing method", title="Aitchison PCoA PC1-PC2 - Genus level")+
    theme_bw()+
    theme(panel.grid=element_blank())+
    scale_color_manual(values=c("#56B4E9","#CC79A7","#E69F00"))
  pcoa.seq.method.PC1.PC2
  ggsave(paste0("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_PILOT_RESULTS/", samples, "_Aitchison_PCoA_plot_X1_X2_sequencing_method_genus.pdf"), pcoa.seq.method.PC1.PC2, dpi=500, width=15, height=10, units="cm", useDingbats=F)
  
}

# only milk samples
run.adonis.and.plot.pcoas.single(df_genus2_clr[df_genus2_clr$n_reads_total>0 & df_genus2_clr$sample_origin_type=="Human_milk", 6:ncol(df_genus2_clr)], samples = "NEXT_pilot_milk_16S_and_16S-ITS-23S_samples_")

# only oral samples
run.adonis.and.plot.pcoas.single(df_genus2_clr[df_genus2_clr$n_reads_total>0 & df_genus2_clr$sample_origin_type=="Infant_oral_cavity", 6:ncol(df_genus2_clr)], samples = "NEXT_pilot_oral_16S_and_16S-ITS-23S_samples_")

# only maternal faeces samples
run.adonis.and.plot.pcoas.single(df_genus2_clr[df_genus2_clr$n_reads_total>0 & df_genus2_clr$sample_origin_type=="Maternal_faeces", 6:ncol(df_genus2_clr)], samples = "NEXT_pilot_mother_feces_16S_and_16S-ITS-23S_samples_")

# only infant faeces samples
run.adonis.and.plot.pcoas.single(df_genus2_clr[df_genus2_clr$n_reads_total>0 & df_genus2_clr$sample_origin_type=="Infant_faeces", 6:ncol(df_genus2_clr)], samples = "NEXT_pilot_infant_feces_16S_and_16S-ITS-23S_samples_")

# adjust p values from seq method for each sample type (order: human milk, infant oral cavity, maternal faeces, infant faeces)
p.adjust(c(0.158, 0.001, 0.001, 0.587), method="BH")
#FDR: 0.2106667 0.0020000 0.0020000 0.5870000


### ===== 7. RELATIVE ABUNDANCE PLOTS BY SAMPLE TYPE AND SEQUENCING METHOD (FIGURE S5) ===== ###

## only use data on genus level

### 7.1 PREPARE DATA FOR PLOTTING, INCL. FILTERING ON MOST PREVALENT AND MOST ABUNDANT BACTERIA

## exclude bacteria that are unclassified on phylum level
df_genus3 <- df_genus[,-grep("g__NA", colnames(df_genus))] #the relative abundance of unclassified genera is going into 'Other'

## shorten colnames
# colnames(df_genus3)[c(10,36,48,61,70,71,73,74,85,95,109,111,116,118,120,128,135,139,148,149,177,179:181,197,199,203)] <- gsub(".*p__", "", colnames(df_genus3)[grep("c__NA", colnames(df_genus3))])
colnames(df_genus3) <- gsub(".*g__", "", colnames(df_genus3))

## sort bacterial columns alphabetically
df_genus3 <- df_genus3[,c(1:5, order(colnames(df_genus3)[6:ncol(df_genus3)])+5)]
                                                       

## prepare data for filtering
rownames(df_genus3) <- c(paste0(df_genus3$sequencing_method, "_", df_genus3$subjectId))
bact_genus <- df_genus3[,6:ncol(df_genus3)]

## abundance filtering: filter on bacteria with at least 2% relative abundance in any sample
maxvector <- c()
for (i in 1:ncol(bact_genus)){maxvector[i] <- max(bact_genus[,i])}
bact_genus[113,] <- maxvector
bact_genus_filt <- bact_genus[,bact_genus[113,]>0.02] #103 genera kept, 136 removed
bact_genus <- bact_genus[-113,]
bact_genus_filt <- bact_genus_filt[-113,]

# prevalence filtering: filter on bacteria present in at least 14 samples
bact_genus_filt2 <- bact_genus_filt[,colSums(bact_genus_filt>0)>=14] #27 genera kept, 76 removed

# sum up excluded genera in 'Other'
bact_genus_filt2$Other <- 1-rowSums(bact_genus_filt2)
df_genus3_filt <- merge(df_genus3[,c(1:5)], bact_genus_filt2, by="row.names") #merge back with metadata
df_genus3_filt <- df_genus3_filt[,-1] #remove Row.names column
df_genus3_filt[df_genus3_filt$n_reads_total==0,33] <- 0 #set 'Other' to 0 for the one sample that failed sequencing


## convert data frame from wide to long format
library(reshape)
genus_long <- melt(df_genus3_filt, id.vars=colnames(df_genus3_filt)[c(1:5)], variable_name="Genus")

## ensure factors show in the desired order in the plot
genus_long$pair_ID <- factor(genus_long$pair_ID,
                             levels = c("Pair1", "Pair2", "Pair3", "Pair4", "Pair5",
                                        "Pair6", "Pair7", "Pair8", "Pair9", "Pair10",
                                        "Pair11", "Pair12", "Pair13", "Pair14"))

genus_long$sample_origin_type <- factor(genus_long$sample_origin_type,
                                        levels = c("Human_milk", "Infant_oral_cavity", "Infant_faeces", "Maternal_faeces"))

genus_long$sequencing_method <- factor(genus_long$sequencing_method,
                                       levels = c("16S-ITS-23S", "16S_V3V4"))


### 7.2 CREATE RELATIVE ABUNDANCE PLOT ON GENUS LEVEL (FIGURE S5)

## create stacked barplot by sample_type and sequencing method
GenusPlot <- ggplot(genus_long, aes(x=pair_ID, y=value, fill=Genus))+
  geom_bar(stat="identity", position="fill")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90),
        text = element_text(size=20),
        legend.position = "bottom")+
  labs(x="", y="Relative abundance", title="Genera relative abundances by sample type and sequencing method")+
  scale_y_continuous(limits=c(0,1),
                     breaks=seq(from=0,to=1,by=0.2),
                     labels=c("0%","20%","40%","60%","80%","100%"))+
  facet_grid(sample_origin_type~sequencing_method, scales="free")+
  scale_fill_manual(values=c("maroon3", "gold3", "darkorange1", "rosybrown", "slategray3",
                             "palegreen", "cadetblue2", "firebrick1", "wheat", "lavender",
                             "blue", "#56B4E9", "mediumpurple", "goldenrod4", "deeppink3",
                             "plum1", "pink4", "#F0E442", "seagreen1", "plum3",
                             "cyan1", "brown4", "linen", "darkgreen", "#E69F00",
                             "navy", "darkseagreen2", "grey"))
GenusPlot

## save plot
# ggsave("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_PILOT_RESULTS/Figure_S5.pdf", GenusPlot, dpi=500, width=40, height=50, units="cm", useDingbats=F)


### ===== 8. CHECK FOR MISSING MOCK BACTERIA IN LIFELINES NEXT SAMPLES (TABLE S8) ===== ###

# Escherichia, Limosilactobacillus and Pseudomonas were missing in the mock communty with 16S sequencing
# -> check if they were detected in biological samples with 16S sequencing

## genus level
rownames(df_genus) <- c(paste0(df_genus$sequencing_method, "_", df_genus$subjectId))
df_genus_16S23S <- df_genus[df_genus$sequencing_method=="16S-ITS-23S",]
df_genus_16S23S_mock_bacteria <- df_genus_16S23S[,grep("Escherichia|Limosilactobacillus|Pseudomonas", colnames(df_genus_16S23S))]
df_genus_16S23S_mock_bacteria_2 <- merge(df_genus_16S23S[,1:5], df_genus_16S23S_mock_bacteria, by="row.names")  
df_genus_16S23S_mock_bacteria_2 <- df_genus_16S23S_mock_bacteria_2[,-1]  
df_genus_16S23S_mock_bacteria_3 <- df_genus_16S23S_mock_bacteria_2[,c(2,3,5,7,6,8,9)]
colnames(df_genus_16S23S_mock_bacteria_3) <- gsub(".*_g__", "", colnames(df_genus_16S23S_mock_bacteria_3))
df_genus_16S23S_mock_bacteria_3 <- df_genus_16S23S_mock_bacteria_3[order(df_genus_16S23S_mock_bacteria_3$sample_origin_type),]
# write.table(df_genus_16S23S_mock_bacteria_3, "MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_PILOT_RESULTS/Table_S8.txt", sep="\t", row.names=F, quote=F)


### ===== 9. CHECK GENUS RICHNESS BY SAMPLE TYPE AND SEQUENCING METHOD (FIGURE S6) ===== ###

### 9.1 PREPARE DATA FRAME

## ensure sample info is shown ni rownames
rownames(data) <- c(paste0(data$sequencing_method, "_", data$pair_ID, "__", data$sample_origin_type))

## subset data by taxonomic level (select only genus level)
data_genus <- subset.by.taxlevel(inputdata = data[,22:ncol(data)], pattern = "_s__.*") #112x280 genera

## only select classified bacterial genera
data_genus_cf <- data_genus[,-grep("g__NA", colnames(data_genus))]

## shorten colnames
colnames(data_genus_cf) <- gsub(".*g__", "", colnames(data_genus_cf))

## split data frame by sequencing method
data_genus_cf_16S <- data_genus_cf[grep("16S_V3V4", rownames(data_genus_cf)),]
data_genus_cf_16S <- data_genus_cf_16S[rownames(data_genus_cf_16S)!="16S_V3V4_Pair5__Maternal_faeces",]
data_genus_cf_16S23S <- data_genus_cf[grep("16S-ITS-23S", rownames(data_genus_cf)),]


### 9.2 CALCULATE CLASSIFIED GENERA RICHNESS PER SAMPLE TYPE AND SEQUENCING METHOD

## 16S data
milk_genus_richness_16S <- as.data.frame(specnumber(data_genus_cf_16S[grep("Human_milk", rownames(data_genus_cf_16S)),]))
inforal_genus_richness_16S <- as.data.frame(specnumber(data_genus_cf_16S[grep("Infant_oral_cavity", rownames(data_genus_cf_16S)),]))
inffec_genus_richness_16S <- as.data.frame(specnumber(data_genus_cf_16S[grep("Infant_faeces", rownames(data_genus_cf_16S)),]))
matfec_genus_richness_16S <- as.data.frame(specnumber(data_genus_cf_16S[grep("Maternal_faeces", rownames(data_genus_cf_16S)),]))

## 16S-ITS-23S data
milk_genus_richness_16S23S <- as.data.frame(specnumber(data_genus_cf_16S23S[grep("Human_milk", rownames(data_genus_cf_16S23S)),]))
inforal_genus_richness_16S23S <- as.data.frame(specnumber(data_genus_cf_16S23S[grep("Infant_oral_cavity", rownames(data_genus_cf_16S23S)),]))
inffec_genus_richness_16S23S <- as.data.frame(specnumber(data_genus_cf_16S23S[grep("Infant_faeces", rownames(data_genus_cf_16S23S)),]))
matfec_genus_richness_16S23S <- as.data.frame(specnumber(data_genus_cf_16S23S[grep("Maternal_faeces", rownames(data_genus_cf_16S23S)),]))

## fix colnames
colnames(milk_genus_richness_16S) <- "genera_richness"
colnames(inforal_genus_richness_16S) <- "genera_richness"
colnames(inffec_genus_richness_16S) <- "genera_richness"
colnames(matfec_genus_richness_16S) <- "genera_richness"

colnames(milk_genus_richness_16S23S) <- "genera_richness"
colnames(inforal_genus_richness_16S23S) <- "genera_richness"
colnames(inffec_genus_richness_16S23S) <- "genera_richness"
colnames(matfec_genus_richness_16S23S) <- "genera_richness"

## combine results in 1 data frame
genera_richness <- as.data.frame(rbind(milk_genus_richness_16S, inforal_genus_richness_16S, inffec_genus_richness_16S, matfec_genus_richness_16S,
                                       milk_genus_richness_16S23S, inforal_genus_richness_16S23S, inffec_genus_richness_16S23S, matfec_genus_richness_16S23S))
genera_richness$sequencing_method <- gsub("_Pair.*", "", rownames(genera_richness))
genera_richness$sample_origin_type <- gsub(".*__", "", rownames(genera_richness))
genera_richness$sequencing_sample_type <- c(paste0(genera_richness$sequencing_method, "_", genera_richness$sample_origin_type))
genera_richness

## ensure factors are ordered as desired for plot
genera_richness$sample_origin_type <- factor(genera_richness$sample_origin_type,
                                             levels = c("Human_milk", "Infant_oral_cavity", "Infant_faeces", "Maternal_faeces"))
genera_richness$sequencing_method <- factor(genera_richness$sequencing_method,
                                            levels = c("16S-ITS-23S", "16S_V3V4"))


### 9.3 PLOT GENERA RICHNESS PER SAMPLE TYPE AND SEQUENCING METHOD

grplot <- ggplot(genera_richness, aes(x=sequencing_method, y=genera_richness))+
  geom_jitter(aes(color=sample_origin_type, size=2), alpha=0.5)+
  geom_boxplot(alpha=0)+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())+
  scale_color_manual(values = c("#E69F00", "#CC79A7", "#009E73", "#0072B2"))+
  labs(x="", y="Genera richness")+
  facet_grid(.~sample_origin_type, scale="free_x")
grplot
# ggsave("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_PILOT_RESULTS/Figure_S6.pdf", grplot, dpi=500, width=15, height=15, units="cm", useDingbats=F)


### 9.4 STATISTICS FOR GENERA RICHNESS PER SAMPLE TYPE AND SEQUENCING METHOD

## Wilcoxon tests
WT_milk <- wilcox.test(genera_richness[genera_richness$sequencing_method=="16S-ITS-23S" & genera_richness$sample_origin_type=="Human_milk","genera_richness"],
                  genera_richness[genera_richness$sequencing_method=="16S_V3V4" & genera_richness$sample_origin_type=="Human_milk","genera_richness"],
                  alternative="two.sided", paired=F, exact=T)

WT_inforal <- wilcox.test(genera_richness[genera_richness$sequencing_method=="16S-ITS-23S" & genera_richness$sample_origin_type=="Infant_oral_cavity","genera_richness"],
                      genera_richness[genera_richness$sequencing_method=="16S_V3V4" & genera_richness$sample_origin_type=="Infant_oral_cavity","genera_richness"],
                      alternative="two.sided", paired=F, exact=T)

WT_inffec <- wilcox.test(genera_richness[genera_richness$sequencing_method=="16S-ITS-23S" & genera_richness$sample_origin_type=="Infant_faeces","genera_richness"],
                      genera_richness[genera_richness$sequencing_method=="16S_V3V4" & genera_richness$sample_origin_type=="Infant_faeces","genera_richness"],
                      alternative="two.sided", paired=F, exact=T)

WT_matfec <- wilcox.test(genera_richness[genera_richness$sequencing_method=="16S-ITS-23S" & genera_richness$sample_origin_type=="Maternal_faeces","genera_richness"],
                      genera_richness[genera_richness$sequencing_method=="16S_V3V4" & genera_richness$sample_origin_type=="Maternal_faeces","genera_richness"],
                      alternative="two.sided", paired=F, exact=T)

## combine results from Wilcoxon tests in 1 data frame
WT_results <- data.frame(sample_origin_type=c("Human_milk", "Infant_oral_cavity", "Infant_faeces", "Maternal_faeces"),
                         p=c(WT_milk$p.value, WT_inforal$p.value, WT_inffec$p.value, WT_matfec$p.value))

## add multiple testing correction
WT_results$FDR <- p.adjust(WT_results$p, method="BH")
WT_results

## save results table
# write.table(WT_results, file="MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Wilcoxon_test_results_genera_richness_by_sequencing_method_and_sample_type.txt", sep="\t", row.names=F, quote=F)

