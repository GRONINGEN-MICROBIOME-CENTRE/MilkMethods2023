################################################################################
##### Phylogenetic analysis of Streptococcus ASVs in mothers and infants
### Author(s): Asier Fernández / Johanne E. Spreckels 
### Last updated: 20th September, 2023
################################################################################

#****************
# Load libraries
#****************
library(stringr)
library(msa)
library(seqinr)
library(ape)
library(phangorn)
library(treeio)
library(ggtree)


#********************
# Defining functions
#********************

#A. MSA_result_genus: generate a MSA of the ASVs of selected genera using ClustalW 
# It outputs a FASTA file with ASVs sequences, a FASTA file with the alignment and a PDF with the alignment (including logos)
# output_dir : full path to directory where one folder will be stored for each genus with its results
# description: information about sequencing technology and pairs of samples to compare
MSA_result_genus <- function(ASVs, genus_name, output_dir, technology) {
  genus <- genus_name #name of the genus
  dir <- paste0(output_dir ,genus)
  dir.create(dir) #Create a directory for the genus
  setwd (dir)  #Move to the directory
  genus_ASVs <- ASVs[grep (genus, names(ASVs))] #select ASVs of the specific genus
  filename <- paste0(paste(genus, sep = "_"), ".fa")
  write.fasta(genus_ASVs, names(genus_ASVs), filename) # generate a FASTA formatted file
  MSA <-  msa(filename, method = "ClustalW", type = "DNA", order = "input") # do the MSA
  ASV_names_for_plot <- gsub("\\..*", "", names(genus_ASVs))# get names of the ASVs sequences to be displayed later in the alignment and the tree
  ASV_names_for_plot <- gsub(paste0(".*", technology, "_"), "", ASV_names_for_plot)
  names(MSA@unmasked) <- ASV_names_for_plot
  try({msaPrettyPrint(MSA, output="pdf", file = paste(genus,  "MSA.pdf", sep= "_"), #get a PDF with the alignment
                      showNames ="left", shadingMode="similar",
                      shadingColors="blues", showLogo="top", logoColors="rasmol",
                      showLegend=T, askForOverwrite=FALSE)})
  msa_align <- msaConvert(MSA, type="seqinr::alignment") #convert the alignment for later processing with seqinr
  msa_align$nam <- ASV_names_for_plot #add names for the plot
  filename_aln <- paste0(paste(genus, sep = "_"), "_aln.fa")
  write.fasta(as.list(msa_align$seq), names(genus_ASVs), filename_aln) #write alignment fasta (as ALN file is empty)
}

#B. Phylo_tree_genus: generate a plot of the phylogenetic tree generated by RaxML-NG
# It outputs a PDF with the phylogenetic tree
# directory : full path to directory where one the output PDF will be created (folder in which subfolders of different genera are present)
# description: information about sequencing technology and genus
Phylo_tree_genus <- function(genus_name, directory) {
  setwd(directory)
  genus <- genus_name
  setwd (genus)  #Move to the directory
  try({ treeML <- read.tree(list.files("./", "support$")) # read tree with bootstrap support
  treeML <- midpoint(treeML) # add midpoint root
  treeML$edge.length <- log(treeML$edge.length*(1/min(treeML$edge.length))) # branch distances to log scale
  treeML$root.edge <- 0.01 # set tree root lenght
  #Create a dataframe with pair_ID, sample type and sample origin information to use it in the plot
  pair_ID <- factor(extracted_numbers <- as.integer(gsub("^.*Pair(\\d+)\\..*$", "\\1",treeML$tip.label)))
  Sample_type <- str_extract(treeML$tip.label, "(Milk|Faeces|Oral_swab)")
  Sample_origin <- str_extract(treeML$tip.label, "(Mother|Infant)")
  data <- data.frame(label=treeML$tip.label, pair_ID=pair_ID, Sample_type=Sample_type, Sample_origin=Sample_origin) # create metadata table
  tree_with_data <- full_join(treeML, data, by='label') #S4 class #important not to do this in the tree used in ggtree (S3) as the root would disappear
  colors <- c("#88aee1", "#29403b", "#35afa6", "#a83250", "#62c049", "#8711ac", 
              "#92a654", "#53348e", "#d48f4d", "#572f3c", "#fe7446", "#1d6d1f", "#e771dd", "#276cb6")
  colors <- colors[1:length(unique(pair_ID))]
  tree <- ggtree(treeML, size =1) + #generate the tree
    geom_rootedge(linewidth = 1) +
    geom_tippoint(aes(colour = tree_with_data@data$pair_ID, shape = interaction(tree_with_data@data$Sample_type, tree_with_data@data$Sample_origin)), size=10) + 
    scale_color_manual(values=colors)  + 
    scale_shape_manual(values = c("Faeces.Infant" = 16, "Milk.Mother" = 17, "Faeces.Mother" = 15, "Oral_swab.Infant" = 18),
                       labels = c("Infant Feces", "Infant Oral", "Maternal Feces", "Maternal Milk")) +
    geom_treescale(y = -5, fontsize = 12, offset = 3.5, linesize =2) +
    labs(colour = "Pair", shape ="Sample Type and Origin") + 
    theme(text=element_text(size = 40), legend.key.size = unit(3.5, 'lines')) +
    coord_flip()
  })
  pdf(file = paste0(genus, ".pdf"),   # Generate PDF with the phylogenetic tree
      width = 35, # The width of the plot in inches
      height = 10)
  print(tree)
  dev.off() 
}

# Set working directory
setwd("~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/")

# Read metadata
metadata <- read.delim("Data_sequencing_method_comparison_n136.txt") 

#########################################################################
# Phylogenetic tree of Streptococcus ASVs in maternal and infant samples
#########################################################################

#################################################
#A. 16S - Phylogenetic tree
#################################################

#*****************
# A1.Get metadata
#*****************
# A1.Get metadata (all mother and infant samples)
metadata_16S <- metadata[metadata$sample_set == "pilot_sample" &
                           metadata$sequencing_method == "16S_V3V4",]

#***********************************************************************
# A2.Get all Streptococcus ASVs in these samples
#***********************************************************************
# For each sample select those ASVs that are present
# Concatenate all sequences in a list
sample_ASVs_16S <- metadata_16S[,22:ncol(metadata_16S)]
rownames(sample_ASVs_16S) <- paste(metadata_16S$sequencing_method,metadata_16S$sample_origin, metadata_16S$sample_type, metadata_16S$pair_ID, sep = "_")
ASV_final_list_16S <- list()

for (i in 1:nrow(sample_ASVs_16S)) {
  sample_name <- rownames(sample_ASVs_16S)[i]
  ASV <- gsub(".*_", "",  colnames(sample_ASVs_16S)[which(sample_ASVs_16S[i,] != 0)]) #sequence
  name_full <- gsub("_ASV.*", "\\1", colnames(sample_ASVs_16S)[which(sample_ASVs_16S[i,] != 0)]) #full-name
  name <- gsub("_s__.*", "", name_full) #remove everything after genus level
  name <- gsub(".*g__", "g__", name) # remove everything before genus level (get genus name)
  indexes_to_update <- which(str_sub(name, - 2, - 2) == "_") #get indexes of genera names that need to be updated (extra undescore in the end)
  name[indexes_to_update] <- sub("_[^_]+$", "", name[indexes_to_update])  #update names    
  if (any(grepl("__NA",name_full))) { #if genus is missing get the lowest-level taxonomy
    NAs <- which(name == "g__NA")
    name [NAs] <-  gsub("..__NA.*", "", name_full[which(name == "g__NA")])
    name [NAs] <-  gsub(".*_", "", name [NAs])
  }
  ASV_list <- ASV
  names(ASV_list) <- paste(rep(sample_name, length(name)), name, sep = ".")
  ASV_final_list_16S <- append(ASV_final_list_16S, ASV_list) #final list of ASVs with names and sequences
}

# Check ASV names in case there are genera that need to be renamed: names(ASV_final_list_16S)
# Get all ASVs of Streptococcus
names(ASV_final_list_16S) <- make.unique(names(ASV_final_list_16S), sep = '_') #distinguish between the names of different ASVs with same taxonomy (adding _1, _2, etc)
Streptococcus_ASVs_16S <- ASV_final_list_16S[grep("Streptococcus", names(ASV_final_list_16S))]

#***********************************************************************
# A3.Get the MSA for Streptococcus
#***********************************************************************
MSA_result_genus(Streptococcus_ASVs_16S, "Streptococcus", "~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/Streptococcus_phylogeny/16S/", "V3V4")
# Run RaxML-NG using the output aln.fa files 

#***********************************************************************
# A4.Plot phylogenetic tree
#***********************************************************************
Phylo_tree_genus( "Streptococcus", "~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/Streptococcus_phylogeny/16S")

#################################################
#B. 16S-23S - Phylogenetic tree
#################################################

#*****************
# A1.Get metadata
#*****************
# A1.Get metadata (all mother and infant samples)
metadata_16S23S <- metadata[metadata$sample_set == "pilot_sample" &
                           metadata$sequencing_method == "16S-ITS-23S",]

#***********************************************************************
# A2.Get all Streptococcus ASVs in these samples
#***********************************************************************
# For each sample select those ASVs that are present
# Concatenate all sequences in a list
sample_ASVs_16S23S <- metadata_16S23S[,22:ncol(metadata_16S23S)]
rownames(sample_ASVs_16S23S) <- paste(metadata_16S23S$sequencing_method,metadata_16S23S$sample_origin, metadata_16S23S$sample_type, metadata_16S23S$pair_ID, sep = "_")
ASV_final_list_16S23S <- list()

for (i in 1:nrow(sample_ASVs_16S23S)) {
  sample_name <- rownames(sample_ASVs_16S23S)[i]
  ASV <- gsub(".*_", "",  colnames(sample_ASVs_16S23S)[which(sample_ASVs_16S23S[i,] != 0)]) #sequence
  name_full <- gsub("_ASV.*", "\\1", colnames(sample_ASVs_16S23S)[which(sample_ASVs_16S23S[i,] != 0)]) #full-name
  name <- gsub("_s__.*", "", name_full) #remove everything after genus level
  name <- gsub(".*g__", "g__", name) # remove everything before genus level (get genus name)
  indexes_to_update <- which(str_sub(name, - 2, - 2) == "_") #get indexes of genera names that need to be updated (extra undescore in the end)
  name[indexes_to_update] <- sub("_[^_]+$", "", name[indexes_to_update])  #update names    
  if (any(grepl("__NA",name_full))) { #if genus is missing get the lowest-level taxonomy
    NAs <- which(name == "g__NA")
    name [NAs] <-  gsub("..__NA.*", "", name_full[which(name == "g__NA")])
    name [NAs] <-  gsub(".*_", "", name [NAs])
  }
  ASV_list <- ASV
  names(ASV_list) <- paste(rep(sample_name, length(name)), name, sep = ".")
  ASV_final_list_16S23S <- append(ASV_final_list_16S23S, ASV_list) #final list of ASVs with names and sequences
}

# Check ASV names in case there are genera that need to be renamed: names(ASV_final_list_16S)
# Get all ASVs of Streptococcus
names(ASV_final_list_16S23S) <- make.unique(names(ASV_final_list_16S23S), sep = '_') #distinguish between the names of different ASVs with same taxonomy (adding _1, _2, etc)
Streptococcus_ASVs_16S23S <- ASV_final_list_16S23S[grep("Streptococcus", names(ASV_final_list_16S23S))]

#***********************************************************************
# A3.Get the MSA for Streptococcus
#***********************************************************************
MSA_result_genus(Streptococcus_ASVs_16S23S, "Streptococcus", "~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/Streptococcus_phylogeny/16S23S/", "16S23S")
# Run RaxML-NG using the output aln.fa files 

#***********************************************************************
# A4.Plot phylogenetic tree
#***********************************************************************
Phylo_tree_genus( "Streptococcus", "~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/Streptococcus_phylogeny/16S23S")
