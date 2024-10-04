# STEP 1. load packages and datasets
library(dplyr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(broom)
library(knitr)
library(NPARC)
library(MASS)
library(car)
#load datasets
load("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 2-RSS analysis/neoantigen_objects.RData")
expression_data <- as.data.frame(exprs(peptides_novel_normalized) )#contains quantitive data. relative abundance
phenotype_data  <- pData(peptides_novel_normalized)  #contains temperature
feature_data<- featureData(peptides_novel_normalized)@data #contains Novel (T/F)

#read Husen meta file and peptide
husen_metadata<-read.delim('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 3B-bio interpretation/inputs/all_cell_lines_sample_meta.txt')
husen_peptide<-read.delim('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 3B-bio interpretation/inputs/novel_filteredBlastSpecAI_normalized_final_peptidetable_canalt.tsv')

#read rss data
rss_data<-read.delim("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 3B-bio interpretation/inputs/rss_data_full.txt")

#subset by type
rss_data_novel<-subset(rss_data,Type=="Novel")
rss_data_conventional<-subset(rss_data,Type=="Conventional")
head(rss_data_novel)

# STEP 2. Data pre processing. select novel peptides, map cell line information from husen's data to ALL

# 1.select novel peptide dataset
novel_peptides<-rownames(rss_data_novel)
conventional_peptides<-rownames(rss_data_conventional)
pep_list<-rownames(rss_data)

#extract novel data
#feature_novel <- feature_data[rownames(feature_data) %in% novel_peptides, ]
#expression_novel<-expression_data[rownames(feature_data) %in% novel_peptides, ]

#remove unnecessary datasets
#rm(expression_data,feature_data,peptides_novel,peptides_novel_vsn,psms_eset,psm_husen)

#2. intersect cell lines that bother present in Husen and ALL
# Modify sample names
phenotype_data$sample_name <- gsub("_BR2", "", phenotype_data$sample_name)
husen_metadata$cell_line <- gsub(" Replicate 2", "", husen_metadata$cell_line)

# Filter datasets
common_samples <- intersect(phenotype_data$sample_name, husen_metadata$cell_line)
common_samples
phenotype_data_filtered <- phenotype_data[phenotype_data$sample_name %in% common_samples, ]
husen_metadata_filtered <- husen_metadata[husen_metadata$cell_line %in% common_samples, ]

# Get sample IDs from the filtered metadata
sample_ids_husen <- husen_metadata_filtered$sample_id
sample_ids_phenotype <- phenotype_data_filtered$sample_id
length(sample_ids_phenotype)
head(sample_ids_phenotype)

# Discard the rows which has valid data in any of the columns other
#than husen_metadata_filtered$sample_id. Do the same to feature_novel. 
#Filter out the rows which has valid data in any of the columns other 
#than phenotype_data_filtered$sample_id. 

# Filter expression_novel
expression_filtered <- expression_data[apply(expression_data[, -which(names(expression_data) %in% sample_ids_phenotype)], 1, function(row) all(is.na(row))), ]
dim(expression_filtered)
#[1] 51777   176
#513->143 rows
detected_peptides<-rownames(expression_filtered)
head(detected_peptides)
length(detected_peptides)

# Assign a new column 'HusenDetect' in rss_data_novel
rss_data$HusenDetect <- rownames(rss_data) %in% detected_peptides

#discard all columns with FALSE Husen detect
rss_data_filtered <- rss_data %>% filter(HusenDetect == TRUE)

# View the modified rss_data_novel
head(rss_data)
table(rss_data$HusenDetect)# T=138/143
head(rss_data_filtered)
table(rss_data_filtered$Type)# T=138 novel

#add log RSS
rss_data_filtered$RSS_log<-log(rss_data_filtered$RSS)
write.table(rss_data_filtered,'/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/inputs/rss_data_SampleFiltered.txt',sep = "\t", row.names = TRUE)
###########
####extra, get peptides_info_type.######
rss_data <- rss_data %>%
  mutate(info = case_when(
    Type == "Conventional" ~ "Conventional",
    Type == "Novel" & HusenDetect == TRUE ~ "Novel",
    Type == "Novel" & HusenDetect == FALSE ~ "Undetected",
    TRUE ~ NA_character_  # This handles any other cases or NAs
  ))
dim(rss_data)
table(rss_data$info)
#Conventional        Novel   Undetected 
#222905          143          370 
#write.table(rss_data,'/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/outputs/peptide_info.txt',sep = "\t")
#add annotation file that contains peptide-proteofom relationship
pf_anno<-read.delim('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 2-RSS analysis/Input/ProteoformAnnotation.txt')
dim(pf_anno)
head(pf_anno,n=2)
head(rss_data,n=2)


#############RSS novel husen preprocessing finished###########

#*choose groups of interest & ioi for the 'alter' group
colnames(rss_data_novel) #you can check here
group_interest<-"HusenDetect" 
unique(rss_data_novel[[group_interest]])
alter_locations <-TRUE

rss_data_conventional_info<-rss_data_novel

#assign alter and null group
novel_detected <- rss_data_conventional_info %>%
  filter(get(group_interest) %in% alter_locations, !is.na(RSS_log), !is.na(get(group_interest)))
novel_undetected <- rss_data_conventional_info %>%
  filter(!get(group_interest) %in% alter_locations, !is.na(RSS_log), !is.na(get(group_interest)))


#######plot and table
table(rss_data_filtered$Type)
generate_boxplot(rss_data_filtered, "Type", "RSS_log")
generate_density_plot_fill(rss_data_filtered,"Type", "RSS_log")

##stat test

# Split the data into two groupsy
group_novel <- rss_data_filtered$RSS[rss_data_filtered$Type == 'Novel']
group_conventional <- rss_data_filtered$RSS[rss_data_filtered$Type == 'Conventional']

# Perform Wilcoxon rank-sum test
wilcox.test(group_novel, group_conventional)

# Perform Kolmogorov-Smirnov test
ks.test(group_novel, group_conventional)

# Perform Brown-Forsythe test
leveneTest(RSS_log ~ Type, data = rss_data_filtered, center = median)

#save the dataset:
#write.table(rss_data_conventional_info,'/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 3B-bio interpretation/outputs/rss_data_novel_detect.txt', sep = "\t", row.names = TRUE)



########
#read rss data
rss_data<-read.delim("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 3B-bio interpretation/inputs/rss_data_full.txt")
rss_data<-rss_data%>%filter(!is.na(RSS))
rss_data$RSS_log<-log(rss_data$RSS)
table(rss_data$Type)
rss_data<- rss_data[
  rss_data$Type == "Conventional" | (rss_data$Type == "Novel" & rownames(rss_data) %in% detected_peptides),
]

table(rss_data$Type)
generate_boxplot(rss_data, "Type", "RSS_log")
generate_density_plot_fill(rss_data,"Type", "RSS_log")

summary(rss_data$cell_line_count)

# Define quantile ranges
quantiles <- c(2, 6, 14, 20, Inf)

# Create a new column `cl_quantile` and assign quantile categories based on `cell_lin_count`
rss_data$cl_quantile <- cut(as.numeric(rss_data$cell_line_count), breaks = quantiles, labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)

# View the updated data
head(rss_data)
generate_density_plot_line(rss_data, "cl_quantile", "RSS_log")

Q1_data <- rss_data %>% filter(cl_quantile == "Q1")
Q2_data <- rss_data %>% filter(cl_quantile == "Q2")
table(Q1_data$Type)
table(Q2_data$Type)

generate_boxplot(Q1_data, "Type", "RSS_log")
generate_density_plot_fill(Q1_data, "Type", "RSS_log")

generate_boxplot(Q2_data, "Type", "RSS_log")
generate_density_plot_fill(Q2_data, "Type", "RSS_log")

# Split the data into two group
group_novel <- Q1_data$RSS[Q1_data$Type == 'Novel']
group_conventional <-Q1_data$RSS[Q1_data$Type == 'Conventional']

group_novel <- Q2_data$RSS[Q2_data$Type == 'Novel']
group_conventional <-Q2_data$RSS[Q2_data$Type == 'Conventional']

# Perform Wilcoxon rank-sum test
wilcox.test(group_novel, group_conventional)
ks.test(group_novel, group_conventional)

#write.table('')
