# STEP 1. load packages and datasets
#load required packages
#BiocManager::install("NPARC")
# Load necessary packages
library(dplyr)
library(purrr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(broom)
library(knitr)
library(NPARC)
#data("stauro_TPP_data_tidy")
#Before applying any further transformations and filters we create a copy of the imported data.
#df <- stauro_TPP_data_tidy

#expression set contains proteofrom melting data
pg_proteoforms<-readRDS("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/inputs/pg_proteoforms_vsn_filtered.RDS")

expression_data <- as.data.frame(exprs(pg_proteoforms) )#contains quantitive data. relative abundance
phenotype_data  <- pData(pg_proteoforms)  #contains temperature
feature_data<- featureData(pg_proteoforms)@data #contains Novel (T/F)
expression_data_t<-as.data.frame(t(expression_data))
#add proteoform_id as a row in expression data
#expression_data$proteoform_id<-rownames(expression_data)

#STEP 2, generate df same as NPARC example

#check if all rowname match
all(rownames(expression_data_t) == rownames(phenotype_data))

#test
proteoform_name<-'A2M_1'
merged_data <- data.frame(proteoform_id = proteoform_name, sample_id = phenotype_data$sample_id, temperature = phenotype_data$temperature,relAbundance=expression_data_t[[proteoform_name]])
head(merged_data)


# Generate the full dataset using purrr::map_dfr for efficient row-binding
full_data <- map_dfr(feature_data$proteoform_id, 
                     ~ phenotype_data %>%
                       mutate(proteoform_id = .x, 
                              relAbundance = expression_data_t[[.x]]))

# Check the dimensions of the full data
dim(full_data)
#should have 20716*176=3046016 rows

#keep the desired columns
full_data <- full_data[, c("proteoform_id", "sample_id", "temperature", "relAbundance")]


###STEP 3. Add other informations

# Filter feature_data to keep only the specified columns using base R
feature_data <- feature_data[, c("proteoform_id", "ioi", "id", "protein_ids", 
                                 "membership", "num_ambiguous_peptides", 
                                 "ambiguous_peptides_only", "ambiguity_ratio")]
head(feature_data)

#merge
full_data_2 <- merge(full_data, feature_data, by = "proteoform_id")
head(full_data_2)
write.table(full_data_2,'/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/outputs/NPARC_data_2.txt',sep = "\t")
#######


###STEP 4. check if these proteoforms contain at least 1 novel peptides
#reload
pf_anno<-read.delim('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 2-RSS analysis/Input/ProteoformAnnotation.txt')
full_data_2<-read.delim('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/outputs/NPARC_data_2.txt')
#load data that contains Type and peptide-membership info
type_info<-read.delim('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/outputs/peptide_info.txt')
pep_pf_info<-read.delim('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/inputs/Pf_Annotation_new.txt')
table(type_info$info)
dim(type_info)
dim(pep_pf_info)
head(type_info,n=3)
head(pep_pf_info,n=3)
length(unique(pep_pf_info$proteoform_id))#[1] 20751
head(unique(pep_pf_info$proteoform_id))
#remove all the numbers and + in peptide
# Remove all numbers and + signs in the 'peptide' column
pep_pf_info$peptide <- gsub("[0-9+\\.]", "", pep_pf_info$peptide)
# View the first few rows of the updated pep_pf_info to confirm the changes
head(pep_pf_info,n=3)

# Assign row names from type_info to pep_pf_info based on matching Peptide.sequence and protein_ids
rownames(pep_pf_info) <- rownames(type_info)[match(paste(pep_pf_info$peptide, pep_pf_info$protein_ids),
                                                   paste(type_info$Peptide.sequence, type_info$protein_ids))]

# View the first few rows of pep_pf_info to confirm the changes
head(pep_pf_info)
#sort out proteoforms that contains novel peptide (detected)
novel_peptides <- type_info$Peptide.sequence[type_info$info == "Undetected"]
head(novel_peptides)
novel_pf <- unique(pep_pf_info$proteoform_id[pep_pf_info$peptide %in% novel_peptides])

