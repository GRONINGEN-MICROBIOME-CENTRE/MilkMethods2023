### ========================== SCRIPT 5 - COMPARISON OF SEQUENCING METHODS ========================== ###

### SCRIPT:          SCRIPT 5 - COMPARISON OF SEQUENCING METHODS
##
## DESCRIPTION:      Script to generate/run
##                    - Figure 2
##                    - Figure S3
##                    - Table S7
##
## AUTHORS:          Johanne Spreckels
##
## NOTES:            To use this script, set correct paths to the following folders/files:
##                    - MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_DATA/Data_sequencing_method_comparison_n136.txt
##                    - MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/
##
## DATE OF CREATION: Script generated for publication in September 2023

### CONTENTS OF THIS FILE
## 0. IMPORT DATA
## 1. CALCULATE RELATIVE ABUNDANCES (BASED ON ASVs)
## 2. SELECT SAMPLES OF INTEREST
## 3. SUBSET SEQUENCING DATA BY TAXONOMIC LEVEL
## 4. MERGE RELATIVE BACTERIAL ABUNDANCES BACK WITH METADATA
## 5. MOCK COMMUNITY COMPOSITION BY SEQUENCING METHOD (FIGURE 2B AND FIGURE S3)
##     5.1 ADD THEORETICAL MOCK COMPOSITION AND PREPARE DATA FOR PLOTTING 
##     5.2 CLR-TRANSFORM RELATIVE ABUNDANCES
##     5.3 CALCULATE EUCLIDIAN DISTANCES FROM CLR-TRANSFORMED ABUNDANCES (= AITCHISON DISTANCES)
##     5.4 GENUS LEVEL: CREATE DENDROGRAMS AND RELATIVE ABUNDANCE PLOTS (FIGURE 2B)
##     5.5 SPECIES LEVEL: CREATE DENDROGRAMS AND RELATIVE ABUNDANCE PLOTS (FIGURE S3)
## 6. MILK RELATIVE ABUNDANCES BY DNA ISOLATION METHOD (FIGURE 2C)
##     6.1 PREPARE DATA FOR PLOTTING, INCL. FILTERING ON MOST PREVALENT AND MOST ABUNDANT BACTERIA
##     6.2 GENUS LEVEL: CREATE RELATIVE ABUNDANCE PLOT (FIGURE 2C)
## 7. PREPARE DATA FROM NEGATIVE CONTROLS FOR ADONIS AND PCoA PLOTS LATER
## 8. AITCHISON DISTANCES, ADONIS AND PCoA PLOTS (FIGURE 2D)
## 9. COMPARISON OF MILK BACTERIAL RELATIVE ABUNDANCES OF MOST PREVALENT GENERA BETWEEN SEQUENCING METHODS (TABLE S7)
## 10. CHECK FOR THE MOST ABUNDANT GENERA IN MILK BY DNA ISOLATION METHOD


### ===== 0. IMPORT DATA ===== ###

## set working directory
setwd("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_DATA/")

## import Data_sequencing_method_comparison_n136.txt
df <- read.table("Data_sequencing_method_comparison_n136.txt", header=T, sep="\t")


### ===== 1. CALCULATE RELATIVE ABUNDANCES (BASED ON ASVs) ===== ###

## set alias as rownames
rownames(df) <- df$alias

## save only columns with absolute bacterial abundances
bact <- df[,22:ncol(df)]

## calculate relative bacterial abundances
data_relab <- (bact/rowSums(bact))
table(rowSums(data_relab), useNA="ifany") #2 sample have rowSums=NA (1 isolation negative control and the maternal faecal sample from pair5 failed 16S sequencing)

## for the negative control that failed sequencing, set relative abundances to 0
data_relab[is.na(rowSums(data_relab)),] <- 0
table(rowSums(data_relab), useNA="ifany") #now the 2 samples have rowSums=0, all others have rowSums=1 (100%)

## merge relative abundances back with metadata
data <- merge(df[,c(1:21)], data_relab, by="row.names")
data <- data[,-1] #remove column Row.names


### ===== 2. SELECT SAMPLES OF INTEREST ===== ###

## set alias as rownames
rownames(data) <- data$alias

## select by sample type
mock    <- data[data$sample_set=="method_test_sample" & data$sample_type=="Mock",]             #2 samples
negctrl <- data[data$sample_set=="method_test_sample" & data$sample_type=="Negative_control",] #7 samples
milk    <- data[data$sample_set=="method_test_sample" & data$sample_type=="Milk",]             #15 samples


### ===== 3. SUBSET SEQUENCING DATA BY TAXONOMIC LEVEL ===== ###

## function to subset data frame by taxonomic level (genus and species levels)
subset.by.taxlevel <- function(inputdata, pattern){ 
  colnames(inputdata) <- gsub(pattern, "", colnames(inputdata)) #shorten colnames after selected taxonomic level
  taxon <- t(rowsum(t(inputdata), group = colnames(inputdata), na.rm = T)) #combine reads from columns with the same name in one column
  taxon <- taxon[,colSums(taxon)>0] #remove columns where colSums==0
  return(taxon)
}

## subset mock data by taxonomic level
mock_genus <- subset.by.taxlevel(inputdata = mock[,22:ncol(mock)], pattern = "_s__.*") #2x8 genera
mock_species <- subset.by.taxlevel(inputdata = mock[,22:ncol(mock)], pattern = "_ASV__.*") #2x10 species

## subset negctrl data by taxonomic level
negctrl_genus <- subset.by.taxlevel(inputdata = negctrl[,22:ncol(negctrl)], pattern = "_s__.*") #7x21 genera

## subset milk data by taxonomic level
milk_genus <- subset.by.taxlevel(inputdata = milk[,22:ncol(milk)], pattern = "_s__.*") #15x65 genera


### ===== 4. MERGE RELATIVE BACTERIAL ABUNDANCES BACK WITH METADATA ===== ###
  
## Mock samples
# genus
df_genus_mock <- merge(mock[,c(2,4,20,21)], mock_genus, by="row.names")
df_genus_mock <- df_genus_mock[,-1]

# species
df_species_mock <- merge(mock[,c(2,4,20,21)], mock_species, by="row.names")
df_species_mock <- df_species_mock[,-1]


## Negative controls
# genus
df_genus_negctrl <- merge(negctrl[,c(2,4,13,20,21)], negctrl_genus, by="row.names")
df_genus_negctrl <- df_genus_negctrl[,-1]


## Milk samples
# genus
df_genus_milk <- merge(milk[,c(2,4,13,20,21)], milk_genus, by="row.names")
df_genus_milk <- df_genus_milk[,-1]


### ===== 5. MOCK COMMUNITY COMPOSITION BY SEQUENCING METHOD (FIGURE 2B AND FIGURE S3) ===== ###

### 5.1 ADD THEORETICAL MOCK COMPOSITION AND PREPARE DATA FOR PLOTTING

## Theoretical mock composition in 16S sequencing:
#    % Bacterial species         Note
#  4.2 Pseudomonas aeruginosa    
# 10.1 Escherichia coli          
# 10.4 Salmonella enterica       
# 18.4 Lactobacillus fermentum   now called Limosilactobacillus fermentum
#  9.9 Enterococcus faecalis     
# 15.5 Staphylococcus aureus     
# 14.1 Listeria monocytogenes    
# 17.4 Bacillus subtilis    

## genus
for (i in c(1:3)){df_genus_mock[,i] <- as.character(df_genus_mock[,i])}
df_genus_mock[3,] <- c("Theor-Mock", "Theor-Mock", "Theor-Mock", NA, 0.174, 0.099, 0.184, 0.141, 0.155, 0.101, 0.104, 0.042)
for (i in c(1:3)){df_genus_mock[,i] <- as.factor(as.character(df_genus_mock[,i]))}
for (i in c(4:12)){df_genus_mock[,i] <- as.numeric(as.character(df_genus_mock[,i]))}
colnames(df_genus_mock) <- gsub(".*g__", "", colnames(df_genus_mock))

## species
for (i in c(1:3)){df_species_mock[,i] <- as.character(df_species_mock[,i])}
df_species_mock[3,] <- c("Theor-Mock", "Theor-Mock", "Theor-Mock", NA, 0, 0, 0.099, 0.184, 0.141, 0, 0.155, 0, 0.104, 0.042)
colnames(df_species_mock) <- gsub("_s__NA", "_unclassified", colnames(df_species_mock))
colnames(df_species_mock) <- gsub(".*s__", "", colnames(df_species_mock))
colnames(df_species_mock) <- gsub(".*g__", "", colnames(df_species_mock))
colnames(df_species_mock) <- gsub("\\.", "_", colnames(df_species_mock))

# add columns for all species that were not identified with 16S
df_species_mock$Bacillus_subtilis <- c(rep(0,2), 0.174)
df_species_mock$Escherichia_coli <- c(rep(0,2), 0.101)

for (i in c(1:3)){df_species_mock[,i] <- as.factor(as.character(df_species_mock[,i]))}
for (i in c(4:16)){df_species_mock[,i] <- as.numeric(as.character(df_species_mock[,i]))}

# sort bacterial columns alphabetically
df_genus_mock <- df_genus_mock[,c(1:4, order(colnames(df_genus_mock)[5:12])+4)]
df_species_mock <- df_species_mock[,c(1:4, order(colnames(df_species_mock)[5:16])+4)]


### 5.2 CLR-TRANSFORM RELATIVE ABUNDANCES

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

# ensure that the subjectIds show as rownames
rownames(df_genus_mock) <- df_genus_mock$subjectId
rownames(df_species_mock) <- df_species_mock$subjectId

## mocks
df_genus_mock_clr <- df_genus_mock[,5:12] #save bacterial columns separately
df_genus_mock_clr <- as.data.frame(do_clr_externalWeighting(df_genus_mock_clr, df_genus_mock_clr)) #clr-transform
df2_genus_mock_clr <- merge(df_genus_mock[,1:4], df_genus_mock_clr, by="row.names") #merge back with metadata
df2_genus_mock_clr <- df2_genus_mock_clr[,-1] #remove rownames column

df_species_mock_clr <- df_species_mock[,5:16] #save bacterial columns separately
df_species_mock_clr <- as.data.frame(do_clr_externalWeighting(df_species_mock_clr, df_species_mock_clr)) #clr-transform
df2_species_mock_clr <- merge(df_species_mock[,1:4], df_species_mock_clr, by="row.names") #merge back with metadata
df2_species_mock_clr <- df2_species_mock_clr[,-1] #remove rownames column


### 5.3 CALCULATE EUCLIDIAN DISTANCES FROM CLR-TRANSFORMED ABUNDANCES (= AITCHISON DISTANCES)

library(vegan)
library(ggplot2)
library(ggdendro)
library(reshape)
library(ggpubr)
library(cowplot)

# ensure that the subjectId show as rownames
rownames(df2_genus_mock_clr)  <- df2_genus_mock_clr$subjectId
rownames(df2_species_mock_clr) <- df2_species_mock_clr$subjectId

# calculate Aitchison distances for each taxonomic level (Aitchison = Euclidian distance on clr-transformed relative abundances)
ait_genus   <- vegdist(df2_genus_mock_clr[,5:ncol(df2_genus_mock_clr)], method="euclidian")
ait_species <- vegdist(df2_species_mock_clr[,5:ncol(df2_species_mock_clr)], method="euclidian") #note: not all bacteria are classified on the species level!!


### 5.4 GENUS LEVEL: CREATE DENDROGRAMS AND RELATIVE ABUNDANCE PLOTS (FIGURE 2B)

## genus
# create dendrogram, extract dendrogram data and plot ggplot dendrogram
dendrogenus <- hclust(ait_genus, method="complete")
dendrodatagenus <- dendro_data(as.dendrogram(dendrogenus), type="rectangle")
dendrogenusPlot <- ggplot(dendrodatagenus$segments)+
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))+
  # geom_text(data=dendrodatagenus$labels, aes(x, y, label=label),
  # hjust=1, angle=90, size=3)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  scale_y_continuous(limits=c(0,5))+
  labs(y="Aitchison\ndistance", title="Mock communities")
dendrogenusPlot

## convert data frame from wide to long format
genus_mock_long <- melt(df_genus_mock, id.vars=colnames(df_genus_mock)[c(1:4)], variable_name="Genus")

## sort samples as in dendrograms
# -> order: "16S-ITS-23S-Mock", "16S-Mock", "Theor-Mock"
genus_mock_long$subjectId <- factor(genus_mock_long$subjectId,
                                    levels=c("16S-ITS-23S-Mock", "16S-Mock", "Theor-Mock"))

## create stacked barplot by DNA isolation method
GenusPlotMock <- ggplot(genus_mock_long, aes(x=subjectId, y=value, fill=Genus))+
  geom_bar(stat="identity", position="fill")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90))+
  labs(x="", y="Relative abundance")+
  scale_y_continuous(breaks=seq(from=0,to=1,by=0.2),
                     labels=c("0%","20%","40%","60%","80%","100%"))+
  scale_fill_manual(values=c("#CC79A7","#009E73","#56B4E9","#F0E442","#D55E00","#0072B2","darkred","#E69F00"))
GenusPlotMock

## save legend separately
genuslegend <- get_legend(GenusPlotMock)
GenusLegend <- as_ggplot(genuslegend)

## save plot without legend
GenusPlotMockPure <- GenusPlotMock+theme(legend.position = "none")

## combine dendrograms and stacked barplots
GenusPlotsMock <- ggdraw(plot_grid(dendrogenusPlot, GenusPlotMockPure, nrow=2, align="v", rel_heights=c(1/3, 1/1)))

## add legend
GenusPlotsMock <- plot_grid(GenusPlotsMock, GenusLegend, nrow=1)
GenusPlotsMock

## save plot
# ggsave("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/Figure_2B.pdf", GenusPlotsMock, dpi=500, width=20, height=12, units="cm", useDingbats=F)


### 5.5 SPECIES LEVEL: CREATE DENDROGRAMS AND RELATIVE ABUNDANCE PLOTS (FIGURE S3)

## species
# create dendrogram, extract dendrogram data and plot ggplot dendrogram
dendrospecies <- hclust(ait_species, method="complete")
dendrodataspecies <- dendro_data(as.dendrogram(dendrospecies), type="rectangle")
dendrospeciesPlot <- ggplot(dendrodataspecies$segments)+
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))+
  # geom_text(data=dendrodataspecies$labels, aes(x, y, label=label),
  # hjust=1, angle=90, size=3)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  scale_y_continuous(limits=c(0,8))+
  labs(y="Aitchison\ndistance", title="Species level")
dendrospeciesPlot

## convert data frame from wide to long format
species_mock_long <- melt(df_species_mock, id.vars=colnames(df_species_mock)[c(1:4)], variable_name="Species")

## sort samples as in dendrograms
# -> order: "16S-ITS-23S-Mock", "16S-Mock", "Theor-Mock"
species_mock_long$subjectId <- factor(species_mock_long$subjectId,
                                      levels=c("Theor-Mock", "16S-ITS-23S-Mock", "16S-Mock"))

## create stacked barplot by DNA isolation method
SpeciesPlotMock <- ggplot(species_mock_long, aes(x=subjectId, y=value, fill=Species))+
  geom_bar(stat="identity", position="fill")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90))+
  labs(x="", y="Relative abundance")+
  scale_y_continuous(breaks=seq(from=0,to=1,by=0.2),
                     labels=c("0%","20%","40%","60%","80%","100%"))+
  scale_fill_manual(values=c("pink","mistyrose1","#CC79A7","#009E73","#56B4E9","lightskyblue1","#F0E442","#D55E00","#0072B2","darkred","#E69F00","goldenrod1"))
SpeciesPlotMock

## save legend separately
specieslegend <- get_legend(SpeciesPlotMock)
SpeciesLegend <- as_ggplot(specieslegend)

## save plot without legend
SpeciesPlotMockPure <- SpeciesPlotMock+theme(legend.position = "none")

## combine dendrograms and stacked barplots
SpeciesPlots <- ggdraw(plot_grid(dendrospeciesPlot, SpeciesPlotMockPure, nrow=2, align="v", rel_heights=c(1/3, 1/1)))

## add legend
SpeciesPlots <- plot_grid(SpeciesPlots, SpeciesLegend, nrow=1)
SpeciesPlots

## save plot
# ggsave("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/Figure_S3.pdf", SpeciesPlots, dpi=500, width=20, height=12, units="cm", useDingbats=F)


### ===== 6. MILK RELATIVE ABUNDANCES BY DNA ISOLATION METHOD (FIGURE 2C) ===== ###

## only use data on genus level

### 6.1 PREPARE DATA FOR PLOTTING, INCL. FILTERING ON MOST PREVALENT AND MOST ABUNDANT BACTERIA

## exclude bacteria that are unclassified on genus level
df_genus_milk_cf <- df_genus_milk[,-grep("g__NA", colnames(df_genus_milk))] #the relative abundance of excluded bacteria is going into 'Other'

## shorten colnames
colnames(df_genus_milk_cf)[c(25)] <- "Coleofasciculaceae_SIO2C1"
colnames(df_genus_milk_cf) <- gsub(".*g__", "", colnames(df_genus_milk_cf))

## sort bacterial columns alphabetically
df_genus_milk_cf <- df_genus_milk_cf[,c(1:5, order(colnames(df_genus_milk_cf)[6:ncol(df_genus_milk_cf)])+5)]

## add column showing replicate number
df_genus_milk_cf$replicate_number <- as.factor(as.character(gsub(".*_", "", df_genus_milk_cf$subjectId)))

## add column showing sequencing methods and replicate numbers
df_genus_milk_cf$method_replicate <- as.factor(as.character(c(paste0(df_genus_milk_cf$sequencing_method, "_", df_genus_milk_cf$replicate_number))))
df_genus_milk_cf <- df_genus_milk_cf[,c(1:3,62,63,4:5,6:61)]

## prepare data for filtering
rownames(df_genus_milk_cf) <- df_genus_milk_cf$alias
bact_genus_milk <- df_genus_milk_cf[,8:63]

## abundance filtering: filter on bacteria with at least 2% relative abundance in any milk sample
maxvector <- c()
for (i in 1:ncol(bact_genus_milk)){maxvector[i] <- max(bact_genus_milk[,i])}
bact_genus_milk[16,] <- maxvector
bact_genus_milk_filt <- bact_genus_milk[,bact_genus_milk[16,]>0.02] #18 genera kept, 38 removed
bact_genus_milk <- bact_genus_milk[-16,]
bact_genus_milk_filt <- bact_genus_milk_filt[-16,]

# prevalence filtering: filter on bacteria present in at least 2 milk samples
bact_genus_milk_filt2 <- bact_genus_milk_filt[,colSums(bact_genus_milk_filt>0)>=2] #18 genera kept, 0 removed

# sum up excluded genera in 'Other'
bact_genus_milk_filt2$Other <- 1-rowSums(bact_genus_milk_filt2)
df_genus_milk_filt <- merge(df_genus_milk_cf[,c(1:7)], bact_genus_milk_filt2, by="row.names") #merge back with metadata
df_genus_milk_filt <- df_genus_milk_filt[,-1] #remove Row.names column

## convert data frame from wide to long format
genus_milk_long <- melt(df_genus_milk_filt, id.vars=colnames(df_genus_milk_filt)[c(1:7)], variable_name="Genus")

## resort method_replicates sequencing method (show 16S-ITS-23S before 16S)
genus_milk_long$method_replicate <- factor(genus_milk_long$method_replicate,
                                           levels = c("16S-ITS-23S_1", "16S-ITS-23S_2",
                                                      "16S_V3V4_1", "16S_V3V4_2", "16S_V3V4_3"))

### 6.2 GENUS LEVEL: CREATE RELATIVE ABUNDANCE PLOT (FIGURE 2C)

## create stacked barplot by sequencing method showing all milk samples next to each other
GenusPlotMilk <- ggplot(genus_milk_long, aes(x=method_replicate, y=value, fill=Genus))+
  geom_bar(stat="identity", position="fill")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90))+
  labs(x="", y="Relative abundance", title="Milk samples")+
  scale_y_continuous(limits=c(0,1.15),
                     breaks=seq(from=0,to=1,by=0.2),
                     labels=c("0%","20%","40%","60%","80%","100%"))+
  facet_grid(.~sample_ID, scales="free")+
  scale_fill_manual(values=c("#CC79A7","lavender","blue","gold","powderblue",
                             "deeppink3","plum1","pink4","darksalmon","darkcyan",
                             "seagreen1","orangered","mistyrose","#0072B2","darkgreen",
                             "#E69F00","navy","darkseagreen2","grey"))
GenusPlotMilk

## save plot
# ggsave("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/Figure_2C.pdf", GenusPlotMilk, dpi=500, width=29, height=14, units="cm", useDingbats=F)


### ===== 7. PREPARE DATA FROM NEGATIVE CONTROLS FOR ADONIS AND PCoA PLOTS LATER ===== ###

## exclude archaea and bacteria that are unclassified on genus level
df_genus_negctrl_cf <- df_genus_negctrl[,-grep("g__NA", colnames(df_genus_negctrl))] #the relative abundance of excluded bacteria is going into 'Other'

## shorten colnames
colnames(df_genus_negctrl_cf)[c(8,13)] <- c("Dysgonomonadaceae_UBA4179","Polyangiaceae_JADJCA01")
colnames(df_genus_negctrl_cf) <- gsub(".*g__", "", colnames(df_genus_negctrl_cf))

## sort bacterial columns alphabetically
df_genus_negctrl_cf <- df_genus_negctrl_cf[,c(1:5, order(colnames(df_genus_negctrl_cf)[6:ncol(df_genus_negctrl_cf)])+5)]

# ## sum up excluded reads (NA on phylum level) in 'Other'
# df_genus_negctrl_cf$Other <- 1-rowSums(df_genus_negctrl_cf[,6:ncol(df_genus_negctrl_cf)])
# df_genus_negctrl_cf$Other[5] <- 0 #ensure that the neg ctrl with 0 reads does not have rel abundances showing in the graph


### ===== 8. AITCHISON DISTANCES, ADONIS AND PCoA PLOTS (FIGURE 2D) ===== ###

## merge data for negctrls and milk samples for genus level

#genus, only including bacteria that were classified on the genus level
df_genus <- as.data.frame(merge(df_genus_negctrl_cf, df_genus_milk_cf[,-c(4:5)], all=T))
for (i in 6:ncol(df_genus)){df_genus[is.na(df_genus[,i]),i] <- 0}
df_genus$method_sample <- c(paste0(df_genus$sequencing_method, "_", df_genus$subjectId))
df_genus <- df_genus[,c(1:4,74,5:73)]
rownames(df_genus) <- df_genus$method_sample

## clr-transform relative abundances
df_genus_clr <- df_genus[,7:74] #save bacterial columns separately
df_genus_clr <- as.data.frame(do_clr_externalWeighting(df_genus_clr, df_genus_clr)) #clr-transform
df2_genus_clr <- merge(df_genus[,1:6], df_genus_clr, by="row.names") #merge back with metadata
df2_genus_clr <- df2_genus_clr[,-1] #remove rownames column
rownames(df2_genus_clr) <- df2_genus_clr$method_sample

#note that 1 negative control has no seq data (0 reads): "16S_V3V4_Isol-NC-1" and ensure that it is included from adonis and PCoAs (choose only milk samples or samples with >0 reads)


## function for adonis and pcoa plots to compare negctrl and milk samples
run.adonis.and.plot.pcoas <- function(inputdata){
  
  ## AITCHISON DISTANCE MATRIX / + METADATA DATA FRAME
  print(paste0("Preparing Aitchison distance matrix"))
  ait <- vegdist(inputdata, method="euclidian")
  aitm <- cmdscale(ait, k=3) #create distance matrix for first 3 components
  aitd <- data.frame(aitm) #convert matrix to data frame
  
  print("Calculating variance explained")
  beta.cmd <- cmdscale(as.matrix(ait), k=2, eig=T)
  PCoA1 <- beta.cmd$eig[1]/sum(beta.cmd$eig)*100
  PCoA2 <- beta.cmd$eig[2]/sum(beta.cmd$eig)*100
  print(paste0("PCoA1 explains ", PCoA1, "% of the variance"))
  print(paste0("PCoA2 explains ", PCoA2, "% of the variance"))
  
  ## add metadata columns to Aitchison distance data frame (ensure that rownames indicate sample IDs)
  print(paste0("Adding metadata to Aitchison distance matrix"))
  # add column for sequencing method
  aitd[grep("16S_V3V4", rownames(aitd)),4] <- "16S_V3V4"
  aitd[grep("16S-ITS-23S", rownames(aitd)),4] <- "16S-ITS-23S"
  
  # add column for sample_type
  aitd[grep("Milk", rownames(aitd)),5] <- "Milk"
  aitd[grep("NC", rownames(aitd)),5] <- "Negative_control"
  
  # add column for sample_name_for_plots
  aitd[grep("Milk-1", rownames(aitd)),6] <- "Milk-1"
  aitd[grep("Milk-2", rownames(aitd)),6] <- "Milk-2"
  aitd[grep("Milk-3", rownames(aitd)),6] <- "Milk-3"
  aitd[grep("NC", rownames(aitd)),6] <- "Negative_control"
  
  ## fix colnames and data structure of Aitchison distance + metadata data frame
  print(paste0("Fixing structure of created data frame"))
  colnames(aitd)[4:6] <- c("sequencing_method", "sample_type", "sample_name_for_plots")
  for (i in 4:6){aitd[,i] <- as.factor(as.character(aitd[,i]))}
  
  
  ## RUN ADONIS
  print(paste0("Running Adonis"))
  print(paste0("Running Adonis for DNA isolation method"))
  adonisMethod <- adonis(ait ~ aitd$sequencing_method)
  adonisMethodRes <- adonisMethod$aov.tab
  
  print(paste0("Running Adonis for sample type"))
  adonisSampleType <- adonis(ait ~ aitd$sample_type)
  adonisSampleTypeRes <- adonisSampleType$aov.tab
  
  print(paste0("Running Adonis for sample"))
  adonisSample <- adonis(ait ~ aitd$sample_name_for_plots)
  adonisSampleRes <- adonisSample$aov.tab
  
  print(paste0("Combining Adonis results in data frame"))
  adonisResults <- data.frame(x = c("sequencing_method", "sample_type", "sample"),
                              R2 = c(adonisMethodRes$R2[1], adonisSampleTypeRes$R2[1], adonisSampleRes$R2[1]),
                              p = c(adonisMethodRes$`Pr(>F)`[1], adonisSampleTypeRes$`Pr(>F)`[1], adonisSampleRes$`Pr(>F)`[1]))
  adonisResults$BH_adj_p <- p.adjust(adonisResults$p, method = "BH")
  
  write.table(adonisResults, paste0("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/negctrl_and_milk_Aitchison_adonis_results_genus.txt"), sep="\t", row.names=F, quote=F)
  
  
  ## PCoA PLOTS
  print(paste0("Start plotting PCoAs"))
  
  # PCoA plot for X1 (PC1) and X2 (PC2), colored by sample_type and shaped by sequencing_method
  pcoa.sample.type <- ggplot(aitd, aes(x=X1, y=X2, color=sample_type, shape=sequencing_method))+
    geom_point(alpha=0.7, aes(size=2))+
    stat_ellipse(aes(group=sample_type))+
    labs(x=paste0("PC1 (", signif(PCoA1,5), "%)"), y=paste0("PC2 (", signif(PCoA2,5), "%)"), color="Sample type", shape="Sequencing method", title="PCoA")+
    theme_bw()+
    theme(panel.grid=element_blank())+
    scale_color_manual(values=c("#E69F00","#000000"))
  pcoa.sample.type
  ggsave(paste0("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/Figure_2D_part1.pdf"), pcoa.sample.type, dpi=500, width=15, height=10, units="cm", useDingbats=F)
  
  # PCoA plot for X1 (PC1) and X2 (PC2), colored by sample and shaped by sequencing_method
  pcoa.sample <- ggplot(aitd, aes(x=X1, y=X2, color=sample_name_for_plots, shape=sequencing_method))+
    geom_point(alpha=0.7, aes(size=2))+
    stat_ellipse(aes(group=sample_type))+
    labs(x=paste0("PC1 (", signif(PCoA1,5), "%)"), y=paste0("PC2 (", signif(PCoA2,5), "%)"), color="Sample", shape="Sequencing method")+
    ggtitle(paste0("Genus level"))+
    theme_bw()+
    theme(panel.grid=element_blank())+
    scale_color_manual(values=c("violetred4","sienna3","goldenrod2","#000000"))
  pcoa.sample
  ggsave(paste0("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/Figure_2D_part2.pdf"), pcoa.sample, dpi=500, width=15, height=10, units="cm", useDingbats=F)
  
}

# negctrl and milk samples
run.adonis.and.plot.pcoas(df2_genus_clr[df2_genus_clr$n_reads_total>0, 7:ncol(df2_genus_clr)])


## function for adonis and pcoa plots for milk samples only
run.adonis.and.plot.pcoas.milk <- function(inputdata){
  
  ## AITCHISON DISTANCE MATRIX / + METADATA DATA FRAME
  print(paste0("Preparing Aitchison distance matrix"))
  ait <- vegdist(inputdata, method="euclidian")
  aitm <- cmdscale(ait, k=3) #create distance matrix for first 3 components
  aitd <- data.frame(aitm) #convert matrix to data frame
  
  ## add metadata columns to Aitchison distance data frame (ensure that rownames indicate sample IDs)
  print(paste0("Adding metadata to Aitchison distance matrix"))
  # add column for sequencing method
  aitd[grep("16S_V3V4", rownames(aitd)),4] <- "16S_V3V4"
  aitd[grep("16S-ITS-23S", rownames(aitd)),4] <- "16S-ITS-23S"
  
  # add column for sample_name_for_plots
  aitd[grep("Milk-1", rownames(aitd)),5] <- "Milk-1"
  aitd[grep("Milk-2", rownames(aitd)),5] <- "Milk-2"
  aitd[grep("Milk-3", rownames(aitd)),5] <- "Milk-3"
  
  ## fix colnames and data structure of Aitchison distance + metadata data frame
  print(paste0("Fixing structure of created data frame"))
  colnames(aitd)[4:5] <- c("sequencing_method", "sample_name_for_plots")
  for (i in 4:5){aitd[,i] <- as.factor(as.character(aitd[,i]))}
  
  
  ## RUN ADONIS
  print(paste0("Running Adonis"))
  print(paste0("Running Adonis for DNA isolation method"))
  adonisMethod <- adonis(ait ~ aitd$sequencing_method)
  adonisMethodRes <- adonisMethod$aov.tab
  
  print(paste0("Running Adonis for sample"))
  adonisSample <- adonis(ait ~ aitd$sample_name_for_plots)
  adonisSampleRes <- adonisSample$aov.tab
  
  print(paste0("Combining Adonis results in data frame"))
  adonisResults <- data.frame(x = c("sequencing_method", "sample"),
                              R2 = c(adonisMethodRes$R2[1], adonisSampleRes$R2[1]),
                              p = c(adonisMethodRes$`Pr(>F)`[1], adonisSampleRes$`Pr(>F)`[1]))
  adonisResults$BH_adj_p <- p.adjust(adonisResults$p, method = "BH")
  
  write.table(adonisResults, paste0("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/milk_Aitchison_adonis_results_genus.txt"), sep="\t", row.names=F, quote=F)
  
  
  ## PCoA plot for X1 (PC1) and X2 (PC2), colored by sample and shaped by DNA_isolation_method
  print(paste0("Start plotting PCoAs"))
  
  pcoa.sample <- ggplot(aitd, aes(x=X1, y=X2, color=sample_name_for_plots, shape=sequencing_method))+
    geom_point(alpha=0.7)+
    labs(x="PC1", y="PC2", color="Sample", shape="Sequencing method")+
    theme_bw()+
    theme(panel.grid=element_blank())+
    scale_color_manual(values=c("#CC79A7","#E69F00","#0072B2"))
  pcoa.sample
  ggsave(paste0("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/milk_Aitchison_PCoA_plot_X1_X2_sample_and_sequencing_method_genus.pdf"), pcoa.sample, dpi=500, width=15, height=10, units="cm", useDingbats=F)
  
}

# milk samples
run.adonis.and.plot.pcoas.milk(df2_genus_clr[grep("Milk", df2_genus_clr$sample_ID),7:ncol(df2_genus_clr)])


### ===== 9. COMPARISON OF MILK BACTERIAL RELATIVE ABUNDANCES OF MOST PREVALENT GENERA BETWEEN SEQUENCING METHODS (TABLE S7) ===== ###

## find most prevalent classified bacterial genera
## 15 milk samples, 70% of 15 -> 11 samples
rownames(df_genus_milk_cf) <- as.factor(as.character(paste0(df_genus_milk_cf$sequencing_method, "_", df_genus_milk_cf$subjectId)))
bact_genus_milk2 <- df_genus_milk_cf[,8:63]
bact_genus_milk_prev <- bact_genus_milk2[,colSums(bact_genus_milk2>0)>=11] #7 genera are present in at least 70% of the milk samples
colnames(bact_genus_milk_prev)
# "Corynebacterium" "Cutibacterium"   "Gemella"         "Rothia"          "Staphylococcus"  "Streptococcus"   "Veillonella"  
# compared to the 16S-SILVA DNA isolation comparison, now "Acinetobacter" (present in 10 milk samples) and "Enhydrobacter" (detected in  0 milk samples with FANGORN database) are missing
# -> compare the relative abundances of the most prevalent bacterial genera between sequencing methods

## add metadata columns
# add sequencing_method column to data frame
bact_genus_milk_prev[grep("16S_V3V4", rownames(bact_genus_milk_prev)),8] <- "16S_V3V4"
bact_genus_milk_prev[grep("16S-ITS-23S", rownames(bact_genus_milk_prev)),8] <- "16S-ITS-23S"

# add sample column to data frame
bact_genus_milk_prev[grep("Milk-1", rownames(bact_genus_milk_prev)),9] <- "Milk-1"
bact_genus_milk_prev[grep("Milk-2", rownames(bact_genus_milk_prev)),9] <- "Milk-2"
bact_genus_milk_prev[grep("Milk-3", rownames(bact_genus_milk_prev)),9] <- "Milk-3"

colnames(bact_genus_milk_prev)[8:9] <- c("sequencing_method", "sample")

## clr-transform relative abundances
bact_genus_milk_prev_clr <- bact_genus_milk_prev
bact_genus_milk_prev_clr[,1:7] <- as.data.frame(do_clr_externalWeighting(bact_genus_milk_prev_clr[,1:7], bact_genus_milk_prev_clr[,1:7]))

## histograms before/after clr-transformation of relative abundances
# before transformation
pdf("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/histograms_milk_rel_ab_raw.pdf")
for (i in c(1:7)){
  hmohist <- ggplot(bact_genus_milk_prev,
                    aes(x=as.numeric(bact_genus_milk_prev[,i]))) + geom_histogram() + labs(x=colnames(bact_genus_milk_prev)[i])
  print(hmohist)
}
dev.off()

# after transformation
pdf("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/histograms_milk_rel_ab_clr-transformed.pdf")
for (i in c(1:7)){
  hmohist <- ggplot(bact_genus_milk_prev_clr,
                    aes(x=as.numeric(bact_genus_milk_prev_clr[,i]))) + geom_histogram() + labs(x=colnames(bact_genus_milk_prev_clr)[i])
  print(hmohist)
}
dev.off()


## Linear model for comparison of clr-trnsformed relative bacterial abundances
# function for lm with correction for milk sample
run.lm.cor.sample <- function(inputdata, ycolumn){
  R2 <- c()
  p_models <- c()
  
  y <- c()
  n <- c()
  stat <- c()
  
  for (j in ycolumn){
    
    #first model without the phenotype of interest (i)
    m0 <- lm(inputdata[,j] ~ inputdata$sample, data=inputdata)
    sum0 <- summary(m0)
    
    #second model with the phenotype of interest (i)
    m1 <- lm(inputdata[,j] ~ inputdata$sample + inputdata$sequencing_method, data=inputdata)
    sum1 <- summary(m1)
    
    #compare models and save R2 diff and p value from model comparison
    an1 <- anova(m0, m1)
    R2 <- c(R2, sum1$r.squared-sum0$r.squared)
    p_models <- c(p_models, an1[2,6])
    
    stat <- c(stat, paste0("lm_correction_sample"))
    
    x <- rep("sequencing_method", length(ycolumn))
    y <- c(y, colnames(inputdata)[j])
    n <- c(n, nrow(inputdata))
  }
  
  #save results in data frame
  res <- data.frame(statistic=stat, x=x, bacteria=y, n_total=n, R2=R2, p_models=p_models)
  
  #multiple testing correction
  res$BH_adj_p_value <- p.adjust(res$p_models, method="BH")
  
  #order rows by BH_adj_p_value and p_value and resort columns
  res <- res[order(res$BH_adj_p_value, res$p_models),]
  return(res)
}

# run model
lmResCorSample <- run.lm.cor.sample(inputdata = bact_genus_milk_prev_clr, ycolumn = c(1:7))
lmResCorSample


## function to retrieve median (range) for untransformed relative bacterial abundances of each highly prevalent bacterium for each sequencing method
getSumStatsRA <- function(inputdata, seq_method){
  RA <- c()
  
  # calculate median and range
  for (i in 1:7){
    RA <- c(RA, paste0(median(inputdata[inputdata$sequencing_method==seq_method,i])*100,
                       " (", min(inputdata[inputdata$sequencing_method==seq_method,i])*100, " - ",
                       max(inputdata[inputdata$sequencing_method==seq_method,i])*100, ")"))
  }
  
  return(RA)
}

SumStatsRA <- as.data.frame(cbind(getSumStatsRA(bact_genus_milk_prev, seq_method="16S_V3V4"),
                                  getSumStatsRA(bact_genus_milk_prev, seq_method="16S-ITS-23S")))
colnames(SumStatsRA) <- c("16S_V3V4","16S-ITS-23S")
SumStatsRA$bacteria <- colnames(bact_genus_milk_prev)[1:7]
library(dplyr)
SumStatsRA <- left_join(SumStatsRA, lmResCorSample[,c(3,6:7)], by="bacteria")
SumStatsRA <- SumStatsRA[,c(3,2,1,4:5)]
SumStatsRA
# write.table(SumStatsRA, "MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/Table_S7.txt", sep="\t", row.names=F, quote=F)


### ===== 10. CHECK FOR THE MOST ABUNDANT GENERA IN MILK BY DNA ISOLATION METHOD ===== ###

## filter on classified bacteria with at least 15% abundance in any milk sample
maxvector <- c()
for (i in 1:ncol(bact_genus_milk2)){maxvector[i] <- max(bact_genus_milk2[,i])}
bact_genus_milk2[16,] <- maxvector
bact_genus_milk_abun <- bact_genus_milk2[,bact_genus_milk2[16,]>=0.15] #5 genera kept, 58 removed
bact_genus_milk2 <- bact_genus_milk2[-16,]
bact_genus_milk_abun <- bact_genus_milk_abun[-16,]

colnames(bact_genus_milk_abun)
# "Corynebacterium" "Cutibacterium"   "Gemella"         "Staphylococcus"  "Streptococcus"  
# compared to the 16S-SILVA DNA isolation comparison, now "Exiguobacterium" (present in 4 samples, 5% RA) is missing


## add metadata columns
# add sequencing_method column to data frame
bact_genus_milk_abun[grep("16S_V3V4", rownames(bact_genus_milk_abun)),6] <- "16S_V3V4"
bact_genus_milk_abun[grep("16S-ITS-23S_", rownames(bact_genus_milk_abun)),6] <- "16S-ITS-23S"

# add sample column to data frame
bact_genus_milk_abun[grep("Milk-1", rownames(bact_genus_milk_abun)),7] <- "Milk-1"
bact_genus_milk_abun[grep("Milk-2", rownames(bact_genus_milk_abun)),7] <- "Milk-2"
bact_genus_milk_abun[grep("Milk-3", rownames(bact_genus_milk_abun)),7] <- "Milk-3"

colnames(bact_genus_milk_abun)[6:7] <- c("sequencing_method", "sample")


## resort and save table showing bacteria with at least 15% relative abundance in any milk sample
bact_genus_milk_abun <- bact_genus_milk_abun[,c(6:7,1:5)]
bact_genus_milk_abun
# write.table(bact_genus_milk_abun, "MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/milk_most_abundant_genera_by_sequencing_method.txt", sep="\t", row.names=F, quote=F)



