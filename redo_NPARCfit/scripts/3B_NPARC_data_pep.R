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

#expression set contains proteofrom melting data
pg_proteoforms<-readRDS("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/inputs/pg_proteoforms_vsn_filtered.RDS")
rss_filter<-read.delim('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/inputs/rss_data_SampleFiltered.txt')
#load datasets
load("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 2-RSS analysis/neoantigen_objects.RData")
expression_data <- as.data.frame(exprs(peptides_novel_normalized) )#contains quantitive data. relative abundance
phenotype_data  <- pData(peptides_novel_normalized)  #contains temperature
feature_data<- featureData(peptides_novel_normalized)@data #contains Novel (T/F)
expression_data_t<-as.data.frame(t(expression_data))

#STEP 2. filter peptides
table(rss_filter$Type)
peptide_list <- rownames(rss_filter)
#conv_list <- rownames(rss_filter[rss_filter$Type == 'Conventional', ])

# Filter phenotype_data
feature_data <- feature_data[rownames(feature_data) %in% peptide_list, ]
feature_data$peptide_id<-rownames(feature_data)
# Filter expression_data
expression_data <- expression_data[rownames(expression_data) %in% peptide_list, ]
expression_data_t<-as.data.frame(t(expression_data))

#STEP 3. generate df same as NPARC example
#check if all rowname match
all(rownames(expression_data_t) == rownames(phenotype_data))
#test
peptide_name<-'peptide145480'
merged_data <- data.frame(peptide_id =peptide_name, sample_id = phenotype_data$sample_id, temperature = phenotype_data$temperature,relAbundance=expression_data_t[[peptide_name]])
head(merged_data)

## Assuming feature_data, phenotype_data, and expression_data_t are already loaded in your environment
full_data <- map_dfr(feature_data$peptide_id, 
                     ~ phenotype_data %>%
                       mutate(peptide_id = .x, 
                              relAbundance = expression_data_t[[.x]]))

# Display the first few rows of the full_data to confirm the changes
head(full_data)
dim(full_data)
#should have 51777*176=9112572 rows

#keep the desired columns
full_data <- full_data[, c("peptide_id", "sample_id", "temperature", "relAbundance","sample_name",'subtype')]
head(full_data)

###STEP 3. Add other informations
# Filter feature_data to keep only the specified columns using base R
feature_data_2 <- feature_data[, c("id", "peptide", "protein_ids", 
                                  'Peptide_AA','Novel',"peptide_id")]
head(feature_data_2)

# add cell line count
feature_data_2$cell_line_count <- rss_filter$cell_line_count[match(rownames(feature_data_2), rownames(rss_filter))]
head(feature_data_2,n=2)

#merge
full_data_2 <- merge(full_data, feature_data_2, by = "peptide_id")
colnames(full_data_2)
head(full_data_2)
dim(full_data_2)
table(full_data_2$Novel)


write.table(full_data_2,'/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/inputs/NPARC_data_PepFilter.txt',sep = "\t")
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

