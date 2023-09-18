### ========================== SCRIPT 2 - COMPARISON OF DNA ISOLATION METHODS ========================== ###

### SCRIPT:          SCRIPT 2 - RESULTS BY DNA ISOLATION METHOD
##
## DESCRIPTION:      Script to generate/run
##                    - Figure 1B-D
##                    - Figure S1
##                    - Figure S2
##                    - Table S3
##                    - check the most abundant genera in milk by DNA isolation method
##
## AUTHORS:          Johanne Spreckels
##
## NOTES:            To use this script, set correct paths to the following folders/files:
##                    - MY_PATH_TO_DNA_ISOLATION_COMPARISON_DATA/Data_DNA_isolation_comparison_n42.txt
##                    - MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/
##
## DATE OF CREATION: Script generated for publication in September 2023


### CONTENTS OF THIS FILE
## 0. IMPORT DATA
## 1. CALCULATE RELATIVE ABUNDANCES (BASED ON ASVs)
## 2. SUBSET SEQUENCING DATA BY SAMPLE TYPE
## 3. SUBSET SEQUENCING DATA BY TAXONOMIC LEVEL
## 4. MERGE RELATIVE BACTERIAL ABUNDANCES BACK WITH METADATA
## 5. MOCK COMMUNITY COMPOSITION BY DNA ISOLATION METHOD (FIGURE 1B AND FIGURE S1)
##     5.1 ADD THEORETICAL MOCK COMPOSITION AND PREPARE DATA FOR PLOTTING
##     5.2 CLR-TRANSFORM RELATIVE ABUNDANCES
##     5.3 CALCULATE EUCLIDIAN DISTANCES FROM CLR-TRANSFORMED ABUNDANCES (= AITCHISON DISTANCES)
##     5.4 GENUS LEVEL: CREATE DENDROGRAMS AND RELATIVE ABUNDANCE PLOTS (FIGURE 1B)
##     5.5 SPECIES LEVEL: CREATE DENDROGRAMS AND RELATIVE ABUNDANCE PLOTS (FIGURE S1)
## 6. NEGATIVE CONTROLS RELATIVE ABUNDANCES BY DNA ISOLATION METHOD (FIGURE S2)
##     6.1 PREPARE DATA FOR PLOTTING, INCL. FILTERING ON MOST PREVALENT AND MOST ABUNDANT BACTERIA
##     6.2 GENUS LEVEL: CREATE RELATIVE ABUNDANCE PLOT (FIGURE S2)
## 7. MILK RELATIVE ABUNDANCES BY DNA ISOLATION METHOD (FIGURE 1C)
##     7.1 PREPARE DATA FOR PLOTTING, INCL. FILTERING ON MOST PREVALENT AND MOST ABUNDANT BACTERIA
##     7.2 GENUS LEVEL: CREATE RELATIVE ABUNDANCE PLOT (FIGURE 1C)
## 8. AITCHISON DISTANCES, ADONIS AND PCoA PLOTS (FIGURE 1D)
## 9. COMPARISON OF MILK BACTERIAL RELATIVE ABUNDANCES OF MOST PREVALENT GENERA BETWEEN DNA ISOLATION METHODS (TABLE S3)
## 10. CHECKING MOST ABUNDANT GENERA IN MILK BY DNA ISOLATION METHOD



### ===== 0. IMPORT DATA ===== ###

## import Data_DNA_isolation_comparison_n42.txt data file
# note that this data file contains ASVs, which are not assigned on certain taxonomic levels, incl. on kingdom and phylum level
data <- read.table("MY_PATH_TO_DNA_ISOLATION_COMPARISON_DATA/Data_DNA_isolation_comparison_n42.txt", header=T, sep="\t", stringsAsFactors=T)


### ===== 1. CALCULATE RELATIVE ABUNDANCES (BASED ON ASVs) ===== ###

## set sample_ids as rownames
rownames(data) <- data$sample_ID

## save only columns with absolute bacterial abundances
bact <- data[,19:ncol(data)]

## calculate relative bacterial abundances
data_relab <- (bact/rowSums(bact))
table(rowSums(data_relab), useNA="ifany") #1 sample has rowSums=NA (the BM_N1_B_NC_1 negative control sample failed sequencing)

## for the negative control that failed sequencing, set relative abundances to 0
data_relab[12,] <- 0
table(rowSums(data_relab), useNA="ifany") #now the 1 sample has rowSums=0, all others have rowSums=1 (100%)

## merge relative abundances back with metadata
data2 <- merge(data[,1:18], data_relab, by="row.names")
data2 <- data2[,-1] #remove column Row.names


### ===== 2. SUBSET SEQUENCING DATA BY SAMPLE TYPE ===== ###

## set sample_ids as rownames
rownames(data2) <- data2$sample_ID

## subset by sample type
mock <- data2[data2$sample_type=="Mock",] #5 samples
milk <- data2[data2$sample_type=="Milk",] #27 samples
negctrl <- data2[data2$sample_type=="Negative_control",]  #10 samples; note that 1 PSK-isolated negative control (BM_N1_B_NC_1) has 0 reads!


### ===== 3. SUBSET SEQUENCING DATA BY TAXONOMIC LEVEL ===== ###

## function to subset data frame by taxonomic level
subset.by.taxlevel <- function(inputdata, pattern){ 
  colnames(inputdata) <- gsub(pattern, "", colnames(inputdata)) #shorten colnames after selected taxonomic level
  taxon <- t(rowsum(t(inputdata), group = colnames(inputdata), na.rm = T)) #combine reads from columns with the same name in one column
  taxon <- taxon[,colSums(taxon)>0] #remove columns where colSums==0
  return(taxon)
}

## subset mock data by taxonomic level
mock_genus <- subset.by.taxlevel(inputdata = mock[,19:ncol(mock)], pattern = "_s__.*") #5x9 genera
mock_species <- subset.by.taxlevel(inputdata = mock[,19:ncol(mock)], pattern = "_ASV__.*") #5x11 species

## subset negctrl data by taxonomic level
negctrl_genus <- subset.by.taxlevel(inputdata = negctrl[,19:ncol(negctrl)], pattern = "_s__.*") #10x75 genera

## subset milk data by taxonomic level
milk_genus <- subset.by.taxlevel(inputdata = milk[,19:ncol(milk)], pattern = "_s__.*") #27x167 genera


### ===== 4. MERGE RELATIVE BACTERIAL ABUNDANCES BACK WITH METADATA ===== ###

## Mock samples

# genus
df_genus_mock <- merge(mock[,c(13,14,15,4,18)], mock_genus, by="row.names")
df_genus_mock <- df_genus_mock[,-1]

# species
df_species_mock <- merge(mock[,c(13,14,15,4,18)], mock_species, by="row.names")
df_species_mock <- df_species_mock[,-1]


## Negative controls

# genus
df_genus_negctrl <- merge(negctrl[,c(13,14,15,4,18)], negctrl_genus, by="row.names")
df_genus_negctrl <- df_genus_negctrl[,-1]


## Milk samples

# genus
df_genus_milk <- merge(milk[,c(13,14,15,4,18)], milk_genus, by="row.names")
df_genus_milk <- df_genus_milk[,-1]


### ===== 5. MOCK COMMUNITY COMPOSITION BY DNA ISOLATION METHOD (FIGURE 1B AND FIGURE S1) ===== ###

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
for (i in 1:4){df_genus_mock[,i] <- as.character(df_genus_mock[,i])}
df_genus_mock[6,] <- c("Theoretical_Mock", "Theoretical_Mock", "Theoretical_Mock", "Theoretical_Mock", NA, 0, 0.174, 0.099, 0.184, 0.141, 0.155, 0.101, 0.104, 0.042)
df_genus_mock$DNA_isolation_method[5] <- "DNA_Mock" # ensure the DNA isolation method column also shows names for the DNA mock for plotting later
for (i in 1:4){df_genus_mock[,i] <- as.factor(as.character(df_genus_mock[,i]))}
for (i in 5:14){df_genus_mock[,i] <- as.numeric(as.character(df_genus_mock[,i]))}
colnames(df_genus_mock) <- gsub(".*g__", "", colnames(df_genus_mock))

## species
for (i in 1:4){df_species_mock[,i] <- as.character(df_species_mock[,i])}
df_species_mock[6,] <- c("Theoretical_Mock", "Theoretical_Mock", "Theoretical_Mock", "Theoretical_Mock", NA, 0, 0, 0, 0, 0, 0.155, 0, 0, 0.104, 0, 0)
df_species_mock$DNA_isolation_method[5] <- "DNA_Mock" # ensure the DNA isolation method column also shows names for the DNA mock for plotting later
colnames(df_species_mock) <- gsub(".*g__", "", colnames(df_species_mock))
colnames(df_species_mock) <- gsub("_s__NA", "_unclassified", colnames(df_species_mock))
colnames(df_species_mock) <- gsub("s__", "", colnames(df_species_mock))

# add columns for all species that were not identified with 16S
df_species_mock$Bacillus_subtilis <- c(rep(0,5), 0.174)
df_species_mock$Enterococcus_faecalis <- c(rep(0,5), 0.099)
df_species_mock$Limosilactobacillus_fermentum <- c(rep(0,5), 0.184)
df_species_mock$Listeria_monocytogenes <- c(rep(0,5), 0.141)
df_species_mock$Escherichia_coli <- c(rep(0,5), 0.101)
df_species_mock$Pseudomonas_aeruginosa <- c(rep(0,5), 0.042)

for (i in 1:4){df_species_mock[,i] <- as.factor(as.character(df_species_mock[,i]))}
for (i in 5:16){df_species_mock[,i] <- as.numeric(as.character(df_species_mock[,i]))}

# sort bacterial columns alphabetically
df_genus_mock <- df_genus_mock[,c(1:5, order(colnames(df_genus_mock)[6:14])+5)]
df_species_mock <- df_species_mock[,c(1:5, order(colnames(df_species_mock)[6:22])+5)]


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

## ensure sample_ID is shown as rownames
rownames(df_genus_mock) <- df_genus_mock$sample_ID
rownames(df_species_mock) <- df_species_mock$sample_ID

## mocks
df_genus_mock_clr <- df_genus_mock[,6:14] #save bacterial columns separately
df_genus_mock_clr <- as.data.frame(do_clr_externalWeighting(df_genus_mock_clr, df_genus_mock_clr)) #clr-transform
df2_genus_mock_clr <- merge(df_genus_mock[,1:5], df_genus_mock_clr, by="row.names") #merge back with metadata
df2_genus_mock_clr <- df2_genus_mock_clr[,-1] #remove rownames column

df_species_mock_clr <- df_species_mock[,6:22] #save bacterial columns separately
df_species_mock_clr <- as.data.frame(do_clr_externalWeighting(df_species_mock_clr, df_species_mock_clr)) #clr-transform
df2_species_mock_clr <- merge(df_species_mock[,1:5], df_species_mock_clr, by="row.names") #merge back with metadata
df2_species_mock_clr <- df2_species_mock_clr[,-1] #remove rownames column


### 5.3 CALCULATE EUCLIDIAN DISTANCES FROM CLR-TRANSFORMED ABUNDANCES (= AITCHISON DISTANCES)

library(vegan)
library(ggplot2)
library(ggdendro)
library(reshape)
library(ggpubr)
library(cowplot)

# ensure that the sample_IDs show as rownames
rownames(df2_genus_mock_clr)   <- df2_genus_mock_clr$sample_ID
rownames(df2_species_mock_clr) <- df2_species_mock_clr$sample_ID

# calculate Aitchison distances for each taxonomic level (Aitchison = Euclidian distance on clr-transformed relative abundances)
ait_genus   <- vegdist(df2_genus_mock_clr[,6:ncol(df2_genus_mock_clr)], method="euclidian")
ait_species <- vegdist(df2_species_mock_clr[,6:ncol(df2_species_mock_clr)], method="euclidian") #note: not all bacteria are classified on the species level!!


### 5.4 GENUS LEVEL: CREATE DENDROGRAMS AND RELATIVE ABUNDANCE PLOTS (FIGURE 1B)

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
  scale_y_continuous(limits=c(0,4))+
  labs(y="Aitchison\ndistance", title="Mock communities")
dendrogenusPlot

## convert data frame from wide to long format
genus_mock_long <- melt(df_genus_mock, id.vars=colnames(df_genus_mock)[c(1:5)], variable_name="Genus")

## sort DNA isolation methods as in dendrograms
# -> order: "FastStoolMiniKit", "MilkDNAKit", "DNA_Mock", "Theoretical_Mock", "PowerSoilProKit", "MagMAXKit"
genus_mock_long$DNA_isolation_method <- factor(genus_mock_long$DNA_isolation_method,
                                                levels=c("FastStoolMiniKit", "MilkDNAKit", "DNA_Mock", "Theoretical_Mock", "PowerSoilProKit", "MagMAXKit"))

## create stacked barplot by DNA isolation method
GenusPlotMock <- ggplot(genus_mock_long, aes(x=DNA_isolation_method, y=value, fill=Genus))+
  geom_bar(stat="identity", position="fill")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90))+
  labs(x="", y="Relative abundance")+
  scale_y_continuous(breaks=seq(from=0,to=1,by=0.2),
                     labels=c("0%","20%","40%","60%","80%","100%"))+
  scale_fill_manual(values=c("#CC79A7", "palegreen","#009E73","#56B4E9","#F0E442","#D55E00","#0072B2","darkred","#E69F00"))
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
# ggsave("MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Figure_1B.pdf", GenusPlotsMock, dpi=500, width=20, height=12, units="cm", useDingbats=F)


### 5.5 SPECIES LEVEL: CREATE DENDROGRAMS AND RELATIVE ABUNDANCE PLOTS (FIGURE S1)

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
  scale_y_continuous(limits=c(0,30))+
  labs(y="Aitchison\ndistance", title="Species level")
dendrospeciesPlot

## convert data frame from wide to long format
species_mock_long <- melt(df_species_mock, id.vars=colnames(df_species_mock)[c(1:5)], variable_name="Species")

## sort DNA isolation methods as in dendrograms
# -> order: "Theoretical_Mock", "FastStoolMiniKit", "MilkDNAKit", "DNA_Mock", "PowerSoilProKit", "MagMAXKit"
species_mock_long$DNA_isolation_method <- factor(species_mock_long$DNA_isolation_method,
                                               levels=c("Theoretical_Mock", "FastStoolMiniKit", "MilkDNAKit", "DNA_Mock", "PowerSoilProKit", "MagMAXKit"))

## create stacked barplot by DNA isolation method
SpeciesPlotMock <- ggplot(species_mock_long, aes(x=DNA_isolation_method, y=value, fill=Species))+
  geom_bar(stat="identity", position="fill")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90))+
  labs(x="", y="Relative abundance")+
  scale_y_continuous(breaks=seq(from=0,to=1,by=0.2),
                     labels=c("0%","20%","40%","60%","80%","100%"))+
  scale_fill_manual(values=c("#CC79A7","pink","palegreen","#009E73","yellowgreen","#56B4E9","lightskyblue1","#F0E442","palegoldenrod","#D55E00","lightsalmon","#0072B2","skyblue3","darkred","indianred3","#E69F00","goldenrod1"))
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
# ggsave("MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Figure_S1.pdf", SpeciesPlots, dpi=500, width=20, height=12, units="cm", useDingbats=F)


### ===== 6. NEGATIVE CONTROLS RELATIVE ABUNDANCES BY DNA ISOLATION METHOD (FIGURE S2) ===== ###

## only use data on genus level

### 6.1 PREPARE DATA FOR PLOTTING, INCL. FILTERING ON MOST PREVALENT AND MOST ABUNDANT BACTERIA

## ensure library preparation negative control is shown in DNA isolation methods column
df_genus_negctrl$DNA_isolation_method <- as.character(as.factor(df_genus_negctrl$DNA_isolation_method))
df_genus_negctrl$DNA_isolation_method[10] <- "Library_preparation"

## exclude bacteria that are unclassified on genus level
df_genus_negctrl_cf <- df_genus_negctrl[,-grep("g__NA", colnames(df_genus_negctrl))] #the relative abundance of excluded bacteria is going into 'Other'

## shorten colnames
colnames(df_genus_negctrl_cf) <- gsub(".*g__", "", colnames(df_genus_negctrl_cf))

## sort bacterial columns alphabetically
df_genus_negctrl_cf <- df_genus_negctrl_cf[,c(1:5, order(colnames(df_genus_negctrl_cf)[6:ncol(df_genus_negctrl_cf)])+5)]

## add column showing isolation number
df_genus_negctrl_cf$DNA_isolation_number <- df_genus_negctrl_cf$sample_ID
df_genus_negctrl_cf$DNA_isolation_number <- gsub(".*_NC_", "NC-", df_genus_negctrl_cf$DNA_isolation_number)
df_genus_negctrl_cf$DNA_isolation_number[10] <- "NC-1"
df_genus_negctrl_cf <- df_genus_negctrl_cf[,c(1:4,63,5:62)]

## prepare data for filtering
rownames(df_genus_negctrl_cf) <- df_genus_negctrl_cf$sample_ID
bact_genus_negctrl <- df_genus_negctrl_cf[,7:ncol(df_genus_negctrl_cf)]

## abundance filtering: filter on bacteria with at least 2% relative abundance in any negative control
maxvector <- c()
for (i in 1:ncol(bact_genus_negctrl)){maxvector[i] <- max(bact_genus_negctrl[,i])}
bact_genus_negctrl[11,] <- maxvector
bact_genus_negctrl_filt <- bact_genus_negctrl[,bact_genus_negctrl[11,]>0.02] #49 genera kept, 8 removed
bact_genus_negctrl <- bact_genus_negctrl[-11,]
bact_genus_negctrl_filt <- bact_genus_negctrl_filt[-11,]

# prevalence filtering: filter on bacteria present in at least 2 negative controls
bact_genus_negctrl_filt2 <- bact_genus_negctrl_filt[,colSums(bact_genus_negctrl_filt>0)>=2] #22 genera kept, 27 removed

# sum up excluded genera in 'Other'
bact_genus_negctrl_filt2$Other <- 1-rowSums(bact_genus_negctrl_filt2)
bact_genus_negctrl_filt2$Other[3] <- 0 #ensure that the neg ctrl with 0 reads does not have rel abundances showing in the graph
df_genus_negctrl_cf_filt <- merge(df_genus_negctrl_cf[,1:6], bact_genus_negctrl_filt2, by="row.names") #merge back with metadata
df_genus_negctrl_cf_filt <- df_genus_negctrl_cf_filt[,-1] #remove Row.names column

## convert data frame from wide to long format
genus_negctrl_long <- melt(df_genus_negctrl_cf_filt, id.vars=colnames(df_genus_negctrl_cf_filt)[c(1:6)], variable_name="Genus")

## sort DNA isolation methods as for mocks
genus_negctrl_long$DNA_isolation_method <- factor(genus_negctrl_long$DNA_isolation_method,
                                               levels=c("FastStoolMiniKit", "MilkDNAKit", "PowerSoilProKit", "MagMAXKit", "Library_preparation"))
levels(genus_negctrl_long$DNA_isolation_method) <- c("FSK", "MDK", "PSK", "MXK", "LibPrep")


### 6.2 GENUS LEVEL: CREATE RELATIVE ABUNDANCE PLOT (FIGURE S2)

## create stacked barplot by DNA isolation method
GenusPlotNegctrl <- ggplot(genus_negctrl_long, aes(x=DNA_isolation_number, y=value, fill=Genus))+
  geom_bar(stat="identity", position="fill")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90))+
  labs(x="", y="Relative abundance", title="Negative controls")+
  scale_y_continuous(limits=c(0,1.15),
                     breaks=seq(from=0,to=1,by=0.2),
                     labels=c("0%","20%","40%","60%","80%","100%"))+
  facet_grid(.~DNA_isolation_method, scales="free")+
  scale_fill_manual(values=c("seagreen3","lightgoldenrod1","#CC79A7","mediumpurple4","lavender","orangered2","blue","turquoise","#56B4E9","goldenrod4","sienna1","#F0E442",
                             "mistyrose","deeppink4","#0072B2","orchid","cornsilk1","#E69F00","navy","palegreen4","darkorchid4","darkseagreen1","grey"))
GenusPlotNegctrl

## save plot
# ggsave("MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Figure_S2.pdf", GenusPlotNegctrl, dpi=500, width=29.1, height=12, units="cm", useDingbats=F)


### ===== 7. MILK RELATIVE ABUNDANCES BY DNA ISOLATION METHOD (FIGURE 1C) ===== ###

## only use data on genus level

### 7.1 PREPARE DATA FOR PLOTTING, INCL. FILTERING ON MOST PREVALENT AND MOST ABUNDANT BACTERIA

## exclude bacteria that are unclassified on genus level
df_genus_milk_cf <- df_genus_milk[,-grep("g__NA", colnames(df_genus_milk))] #the relative abundance of excluded bacteria is going into 'Other'

## shorten colnames
colnames(df_genus_milk_cf)[c(8,64,71,77,79,80,85)] <- c("Actinomycetaceae_F0332", "Ruminococcaceae_CAG.352", "Peptostreptococcales.Tissierellales_Family.XI_W5053", "Myxococcaceae_P3OB.42", "Saccharimonadaceae_TM7a", "Phycisphaeraceae_SM1A02", "Caulobacteraceae_PMMR1")
colnames(df_genus_milk_cf) <- gsub(".*g__", "", colnames(df_genus_milk_cf))

## sort bacterial columns alphabetically
df_genus_milk_cf <- df_genus_milk_cf[,c(1:5, order(colnames(df_genus_milk_cf)[6:ncol(df_genus_milk_cf)])+5)]

## add column showing isolation number
df_genus_milk_cf$DNA_isolation_number <- df_genus_milk_cf$sample_ID
df_genus_milk_cf$DNA_isolation_number <- gsub(".*_", "", df_genus_milk_cf$DNA_isolation_number)

## add column showing DNA isolation method and DNA isolation number
df_genus_milk_cf$method_isolation <- c(paste0(df_genus_milk_cf$DNA_isolation_method, "_", df_genus_milk_cf$DNA_isolation_number))
df_genus_milk_cf <- df_genus_milk_cf[,c(1:4,130:131,5:129)]

## prepare data for filtering
rownames(df_genus_milk_cf) <- df_genus_milk_cf$sample_ID
bact_genus_milk <- df_genus_milk_cf[,8:ncol(df_genus_milk_cf)]

## abundance filtering: filter on bacteria with at least 2% relative abundance in any milk sample
maxvector <- c()
for (i in 1:ncol(bact_genus_milk)){maxvector[i] <- max(bact_genus_milk[,i])}
bact_genus_milk[28,] <- maxvector
bact_genus_milk_filt <- bact_genus_milk[,bact_genus_milk[28,]>0.02] #31 genera kept, 93 removed
bact_genus_milk <- bact_genus_milk[-28,]
bact_genus_milk_filt <- bact_genus_milk_filt[-28,]

# prevalence filtering: filter on bacteria present in at least 2 milk samples
bact_genus_milk_filt2 <- bact_genus_milk_filt[,colSums(bact_genus_milk_filt>0)>=2] #24 genera kept, 7 removed

# sum up excluded genera in 'Other'
bact_genus_milk_filt2$Other <- 1-rowSums(bact_genus_milk_filt2)
df_genus_milk_cf_filt <- merge(df_genus_milk_cf[,1:7], bact_genus_milk_filt2, by="row.names") #merge back with metadata
df_genus_milk_cf_filt <- df_genus_milk_cf_filt[,-1] #remove Row.names column

## convert data frame from wide to long format
genus_milk_long <- melt(df_genus_milk_cf_filt, id.vars=colnames(df_genus_milk_cf_filt)[c(1:7)], variable_name="Genus")

## sort DNA isolation methods as for mocks
genus_milk_long$method_isolation <- factor(genus_milk_long$method_isolation,
                                           levels=c("FastStoolMiniKit_1", "FastStoolMiniKit_2",
                                                    "MilkDNAKit_1", "MilkDNAKit_2",
                                                    "PowerSoilProKit_1", "PowerSoilProKit_2", "PowerSoilProKit_3",
                                                    "MagMAXKit_1", "MagMAXKit_2"))


### 7.2 GENUS LEVEL: CREATE RELATIVE ABUNDANCE PLOT (FIGURE 1C)

## create stacked barplot by milk sample showing all DNA isolation methods next to each other
GenusPlotMilk <- ggplot(genus_milk_long, aes(x=method_isolation, y=value, fill=Genus))+
  geom_bar(stat="identity", position="fill")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90))+
  labs(x="", y="Relative abundance", title="Milk samples")+
  scale_y_continuous(limits=c(0,1.15),
                     breaks=seq(from=0,to=1,by=0.2),
                     labels=c("0%","20%","40%","60%","80%","100%"))+
  facet_grid(.~subjectId, scales="free")+
  scale_fill_manual(values=c("seagreen3","lightgoldenrod1","#CC79A7","lavender","blue","turquoise","#56B4E9","gold","sienna1","thistle4","deeppink3","plum1","paleturquoise2",
                             "darksalmon","chocolate4","firebrick3","darkcyan","burlywood","mistyrose","#0072B2","darkgreen","#E69F00","navy","darkseagreen2","grey"))
GenusPlotMilk

## save plot
# ggsave("MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Figure_1C.pdf", GenusPlotMilk, dpi=500, width=30, height=12, units="cm", useDingbats=F)


### ===== 8. AITCHISON DISTANCES, ADONIS AND PCoA PLOTS (FIGURE 1D) ===== ###

#genus, only including bacteria that were classified on the genus level
df_genus <- as.data.frame(merge(df_genus_negctrl_cf, df_genus_milk_cf, all=T))
df_genus <- df_genus[,c(1:5,64,6,7:63,65:147)]
for (i in 8:ncol(df_genus)){df_genus[is.na(df_genus[,i]),i] <- 0}
rownames(df_genus) <- df_genus$sample_ID

## clr-transform relative abundances
df_genus_clr <- df_genus[,8:147] #save bacterial columns separately
df_genus_clr <- as.data.frame(do_clr_externalWeighting(df_genus_clr, df_genus_clr)) #clr-transform
df2_genus_clr <- merge(df_genus[,1:7], df_genus_clr, by="row.names") #merge back with metadata
df2_genus_clr <- df2_genus_clr[,-1] #remove rownames column
rownames(df2_genus_clr) <- df2_genus_clr$sample_ID


## function for adonis and pcoa plots to compare negctrl and milk samples (PCoA plots are used to generate Figure 1D)
run.adonis.and.plot.pcoas <- function(inputdata, samples, taxlevel){
  
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
  # add column for DNA isolation method
  aitd[grep("_A_", rownames(aitd)),4] <- "MilkDNAKit"
  aitd[grep("_B_", rownames(aitd)),4] <- "PowerSoilProKit"
  aitd[grep("_D_", rownames(aitd)),4] <- "MagMAXKit"
  aitd[grep("_E_", rownames(aitd)),4] <- "FastStoolMiniKit"
  aitd[grep("_NTC", rownames(aitd)),4] <- "Library_preparation"
  
  # add column for sample_type
  aitd[grep("_M", rownames(aitd)),5] <- "Milk"
  aitd[grep("_NC", rownames(aitd)),5] <- "Negative_control"
  aitd[grep("_NTC", rownames(aitd)),5] <- "Negative_control"
  
  # add column for sample_name_for_plots
  aitd[grep("_M1", rownames(aitd)),6] <- "Milk-3"
  aitd[grep("_M2", rownames(aitd)),6] <- "Milk-2"
  aitd[grep("_M3", rownames(aitd)),6] <- "Milk-1"
  aitd[grep("_NC", rownames(aitd)),6] <- "Negative_control"
  aitd[grep("_NTC", rownames(aitd)),6] <- "Negative_control"
  
  ## fix colnames and data structure of Aitchison distance + metadata data frame
  print(paste0("Fixing structure of created data frame"))
  colnames(aitd)[4:6] <- c("DNA_isolation_method", "sample_type", "sample_name_for_plots")
  for (i in 4:6){aitd[,i] <- as.factor(as.character(aitd[,i]))}
  
  ## sort levels for DNA isolation method
  aitd$DNA_isolation_method <- factor(aitd$DNA_isolation_method,
                                     levels=c("FastStoolMiniKit", "MilkDNAKit", "PowerSoilProKit", "MagMAXKit", "Library_preparation"))

  
  ## RUN ADONIS
  print(paste0("Running Adonis"))
  print(paste0("Running Adonis for DNA isolation method"))
  adonisMethod <- adonis(ait ~ aitd$DNA_isolation_method)
  adonisMethodRes <- adonisMethod$aov.tab

  print(paste0("Running Adonis for sample type"))
  adonisSampleType <- adonis(ait ~ aitd$sample_type)
  adonisSampleTypeRes <- adonisSampleType$aov.tab

  print(paste0("Running Adonis for sample"))
  adonisSample <- adonis(ait ~ aitd$sample_name_for_plots)
  adonisSampleRes <- adonisSample$aov.tab

  print(paste0("Combining Adonis results in data frame"))
  adonisResults <- data.frame(x = c("DNA_isolation_method", "sample_type", "sample"),
                              R2 = c(adonisMethodRes$R2[1], adonisSampleTypeRes$R2[1], adonisSampleRes$R2[1]),
                              p = c(adonisMethodRes$`Pr(>F)`[1], adonisSampleTypeRes$`Pr(>F)`[1], adonisSampleRes$`Pr(>F)`[1]))
  adonisResults$BH_adj_p <- p.adjust(adonisResults$p, method = "BH")

  write.table(adonisResults, paste0("MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/", samples, "_Aitchison_adonis_results.", taxlevel, ".txt"), sep="\t", row.names=F, quote=F)


  ## PCoA PLOTS
  print(paste0("Start plotting PCoAs"))
  
  # PCoA plot for X1 (PC1) and X2 (PC2), colored by sample_type and shaped by DNA_isolation_method
  pcoa.sample.type <- ggplot(aitd, aes(x=X1, y=X2, color=sample_type, shape=DNA_isolation_method))+
    geom_point(alpha=0.7, aes(size=2))+
    stat_ellipse(aes(group=sample_type))+
    labs(x=paste0("PC1 (", signif(PCoA1,5), "%)"), y=paste0("PC2 (", signif(PCoA2,5), "%)"), color="Sample type", shape="DNA isolation method", title="PCoA")+
    theme_bw()+
    theme(panel.grid=element_blank())+
    scale_color_manual(values=c("#E69F00","#000000"))
  pcoa.sample.type
  ggsave(paste0("MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/", samples, "_Aitchison_PCoA_plot_X1_X2_sample_type_and_DNA_isolation_method.", taxlevel, ".pdf"), pcoa.sample.type, dpi=500, width=15, height=10, units="cm", useDingbats=F)

  # PCoA plot for X1 (PC1) and X2 (PC2), colored by sample and shaped by DNA_isolation_method
  pcoa.sample <- ggplot(aitd, aes(x=X1, y=X2, color=sample_name_for_plots, shape=DNA_isolation_method))+
    geom_point(alpha=0.7, aes(size=2))+
    stat_ellipse(aes(group=sample_type))+
    labs(x=paste0("PC1 (", signif(PCoA1,5), "%)"), y=paste0("PC2 (", signif(PCoA2,5), "%)"), color="Sample", shape="DNA isolation method")+
    ggtitle(paste0(taxlevel, " level"))+
    theme_bw()+
    theme(panel.grid=element_blank())+
    scale_color_manual(values=c("violetred4","sienna3","goldenrod2","#000000"))
  pcoa.sample
  ggsave(paste0("MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/", samples, "_Aitchison_PCoA_plot_X1_X2_sample_and_DNA_isolation_method.", taxlevel, ".pdf"), pcoa.sample, dpi=500, width=15, height=10, units="cm", useDingbats=F)

  # return(pcoa.sample.type)
  # return(pcoa.sample)
}

# negctrl and milk samples
run.adonis.and.plot.pcoas(df2_genus_clr[df2_genus_clr$n_reads_dada2_final>0,8:ncol(df2_genus_clr)], samples = "negctrl_and_milk", taxlevel = "genus")


## function for adonis and pcoa plots for milk samples only
run.adonis.and.plot.pcoas.milk <- function(inputdata, taxlevel){
  
  ## AITCHISON DISTANCE MATRIX / + METADATA DATA FRAME
  print(paste0("Preparing Aitchison distance matrix"))
  ait <- vegdist(inputdata, method="euclidian")
  aitm <- cmdscale(ait, k=3) #create distance matrix for first 3 components
  aitd <- data.frame(aitm) #convert matrix to data frame
  
  ## add metadata columns to Aitchison distance data frame (ensure that rownames indicate sample IDs)
  print(paste0("Adding metadata to Aitchison distance matrix"))
  # add column for DNA_isolation_method
  aitd[grep("_A_", rownames(aitd)),4] <- "MilkDNAKit"
  aitd[grep("_B_", rownames(aitd)),4] <- "PowerSoilProKit"
  aitd[grep("_D_", rownames(aitd)),4] <- "MagMAXKit"
  aitd[grep("_E_", rownames(aitd)),4] <- "FastStoolMiniKit"
  
  # add column for sample_name_for_plots
  aitd[grep("_M1", rownames(aitd)),5] <- "Milk-3"
  aitd[grep("_M2", rownames(aitd)),5] <- "Milk-2"
  aitd[grep("_M3", rownames(aitd)),5] <- "Milk-1"
  
  ## fix colnames and data structure of Aitchison distance + metadata data frame
  print(paste0("Fixing structure of created data frame"))
  colnames(aitd)[4:5] <- c("DNA_isolation_method", "sample_name_for_plots")
  for (i in 4:5){aitd[,i] <- as.factor(as.character(aitd[,i]))}
  
  ## sort levels for DNA_isolation_method
  aitd$DNA_isolation_method <- factor(aitd$DNA_isolation_method,
                                     levels=c("FastStoolMiniKit", "MilkDNAKit", "PowerSoilProKit", "MagMAXKit"))
  
  ## sort levels for sample_name_for_plots
  aitd$sample_name_for_plots <- factor(aitd$sample_name_for_plots,
                                      levels=c("Milk-1", "Milk-2", "Milk-3"))
  
  
  ## RUN ADONIS
  print(paste0("Running Adonis"))
  print(paste0("Running Adonis for DNA isolation method"))
  adonisMethod <- adonis(ait ~ aitd$DNA_isolation_method)
  adonisMethodRes <- adonisMethod$aov.tab
  
  print(paste0("Running Adonis for sample"))
  adonisSample <- adonis(ait ~ aitd$sample_name_for_plots)
  adonisSampleRes <- adonisSample$aov.tab
  
  adonisResults <- data.frame(x = c("DNA_isolation_method", "sample"),
                              R2 = c(adonisMethodRes$R2[1], adonisSampleRes$R2[1]),
                              p = c(adonisMethodRes$`Pr(>F)`[1], adonisSampleRes$`Pr(>F)`[1]))
  adonisResults$BH_adj_p <- p.adjust(adonisResults$p, method = "BH")
  
  write.table(adonisResults, paste0("MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/milk_Aitchison_adonis_results.", taxlevel, ".txt"), sep="\t", row.names=F, quote=F)
  
  
  ## PCoA plot for X1 (PC1) and X2 (PC2), colored by sample and shaped by DNA_isolation_method
  print(paste0("Start plotting PCoAs"))
  pcoa.sample <- ggplot(aitd, aes(x=X1, y=X2, color=sample_name_for_plots, shape=DNA_isolation_method))+
    geom_point(alpha=0.7)+
    labs(x="PC1", y="PC2", color="Sample", shape="DNA isolation method")+
    theme_bw()+
    theme(panel.grid=element_blank())+
    scale_color_manual(values=c("#CC79A7","#E69F00","#0072B2"))
  pcoa.sample
  ggsave(paste0("MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/milk_Aitchison_PCoA_plot_X1_X2_sample_and_DNA_isolation_method.", taxlevel, ".pdf"), pcoa.sample, dpi=500, width=15, height=10, units="cm", useDingbats=F)
}

# only milk samples
run.adonis.and.plot.pcoas.milk(df2_genus_clr[df2_genus_clr$sample_type=="Milk",8:ncol(df2_genus_clr)], taxlevel = "genus")


### ===== 9. COMPARISON OF MILK BACTERIAL RELATIVE ABUNDANCES OF MOST PREVALENT GENERA BETWEEN DNA ISOLATION METHODS (TABLE S3) ===== ###

## find most prevalent classified bacterial genera
## 27 milk samples, 70% of 27 -> 19 samples
bact_genus_milk_prev <- bact_genus_milk[,colSums(bact_genus_milk>0)>=19] #9 genera are present in at least 70% of the milk samples
colnames(bact_genus_milk_prev)
# "Acinetobacter"   "Corynebacterium" "Cutibacterium"   "Enhydrobacter"   "Gemella"         "Rothia"          "Staphylococcus"  "Streptococcus"   "Veillonella"
# -> compare the relative abundances of the most prevalent bacterial genera between DNA isolation methods

## add metadata columns
# add DNA_isolation_method column to data frame
bact_genus_milk_prev[grep("A_", rownames(bact_genus_milk_prev)),10] <- "MilkDNAKit"
bact_genus_milk_prev[grep("B_", rownames(bact_genus_milk_prev)),10] <- "PowerSoilProKit"
bact_genus_milk_prev[grep("D_", rownames(bact_genus_milk_prev)),10] <- "MagMAXKit"
bact_genus_milk_prev[grep("E_", rownames(bact_genus_milk_prev)),10] <- "FastStoolMiniKit"

# add sample column to data frame
bact_genus_milk_prev[grep("M1", rownames(bact_genus_milk_prev)),11] <- "Milk-3"
bact_genus_milk_prev[grep("M2", rownames(bact_genus_milk_prev)),11] <- "Milk-2"
bact_genus_milk_prev[grep("M3", rownames(bact_genus_milk_prev)),11] <- "Milk-1"

colnames(bact_genus_milk_prev)[10:11] <- c("DNA_isolation_method", "sample")

## clr-transform relative abundances
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

bact_genus_milk_prev_clr <- bact_genus_milk_prev
bact_genus_milk_prev_clr[,1:9] <- as.data.frame(do_clr_externalWeighting(bact_genus_milk_prev_clr[,1:9], bact_genus_milk_prev_clr[,1:9]))

## histograms before/after clr-transformation of relative abundances
# before transformation
pdf("MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/histograms_milk_rel_ab_raw.pdf")
for (i in c(1:9)){
  hmohist <- ggplot(bact_genus_milk_prev,
                    aes(x=as.numeric(bact_genus_milk_prev[,i]))) + geom_histogram() + labs(x=colnames(bact_genus_milk_prev)[i])
  print(hmohist)
}
dev.off()

# after transformation
pdf("MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/histograms_milk_rel_ab_clr-transformed.pdf")
for (i in c(1:9)){
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
    m1 <- lm(inputdata[,j] ~ inputdata$sample + inputdata$DNA_isolation_method, data=inputdata)
    sum1 <- summary(m1)
    
    #compare models and save R2 diff and p value from model comparison
    an1 <- anova(m0, m1)
    R2 <- c(R2, sum1$r.squared-sum0$r.squared)
    p_models <- c(p_models, an1[2,6])
    
    stat <- c(stat, paste0("lm_correction_sample"))
    
    x <- rep("DNA_isolation_method", length(ycolumn))
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
lmResCorSample <- run.lm.cor.sample(inputdata = bact_genus_milk_prev_clr, ycolumn = c(1:9))
lmResCorSample


## function to retrieve median (range) for untransformed relative bacterial abundances of each highly prevalent bacterium for each DNA isolation method
getSumStatsRA <- function(inputdata, isolation_method){
  RA <- c()
  
  # calculate median and range
  for (i in 1:9){
    RA <- c(RA, paste0(median(inputdata[inputdata$DNA_isolation_method==isolation_method,i])*100,
            " (", min(inputdata[inputdata$DNA_isolation_method==isolation_method,i])*100, " - ",
            max(inputdata[inputdata$DNA_isolation_method==isolation_method,i])*100, ")"))
  }
  
  return(RA)
}

SumStatsRA <- as.data.frame(cbind(getSumStatsRA(bact_genus_milk_prev, isolation_method="FastStoolMiniKit"),
                                  getSumStatsRA(bact_genus_milk_prev, isolation_method="MilkDNAKit"),
                                  getSumStatsRA(bact_genus_milk_prev, isolation_method="PowerSoilProKit"),
                                  getSumStatsRA(bact_genus_milk_prev, isolation_method="MagMAXKit")))
colnames(SumStatsRA) <- c("FSK", "MDK", "PSK", "MXK")
SumStatsRA$bacteria <- colnames(bact_genus_milk_prev)[1:9]
library(dplyr)
SumStatsRA <- left_join(SumStatsRA, lmResCorSample[,c(3,6:7)], by="bacteria")
SumStatsRA <- SumStatsRA[,c(5,1:4,6:7)]
SumStatsRA
# write.table(SumStatsRA, "MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Table_S3.txt", sep="\t", row.names=F, quote=F)


### ===== 10. CHECKING MOST ABUNDANT GENERA IN MILK BY DNA ISOLATION METHOD ===== ###

## filter on classified bacteria with at least 15% abundance in any milk sample
maxvector <- c()
for (i in 1:ncol(bact_genus_milk)){maxvector[i] <- max(bact_genus_milk[,i])}
bact_genus_milk[28,] <- maxvector
bact_genus_milk_abun <- bact_genus_milk[,bact_genus_milk[28,]>=0.15]
bact_genus_milk <- bact_genus_milk[-28,]
bact_genus_milk_abun <- bact_genus_milk_abun[-28,]

colnames(bact_genus_milk_abun)
# "Corynebacterium" "Cutibacterium"   "Exiguobacterium" "Gemella"         "Staphylococcus"  "Streptococcus" 


## add metadata columns
# add DNA_isolation_method column to data frame
bact_genus_milk_abun[grep("A_", rownames(bact_genus_milk_abun)),7] <- "MilkDNAKit"
bact_genus_milk_abun[grep("B_", rownames(bact_genus_milk_abun)),7] <- "PowerSoilProKit"
bact_genus_milk_abun[grep("D_", rownames(bact_genus_milk_abun)),7] <- "MagMAXKit"
bact_genus_milk_abun[grep("E_", rownames(bact_genus_milk_abun)),7] <- "FastStoolMiniKit"

# add sample column to data frame
bact_genus_milk_abun[grep("M1", rownames(bact_genus_milk_abun)),8] <- "Milk-3"
bact_genus_milk_abun[grep("M2", rownames(bact_genus_milk_abun)),8] <- "Milk-2"
bact_genus_milk_abun[grep("M3", rownames(bact_genus_milk_abun)),8] <- "Milk-1"

colnames(bact_genus_milk_abun)[7:8] <- c("DNA_isolation_method", "sample")


## resort and save table showing bacteria with at least 15% relative abundance in any milk sample
bact_genus_milk_abun <- bact_genus_milk_abun[,c(7:8,1:6)]
bact_genus_milk_abun
# write.table(bact_genus_milk_abun, "MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/milk_most_abundant_genera_by_DNA_isolation_method.txt", sep="\t", row.names=F, quote=F)





