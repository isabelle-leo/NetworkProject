#load datasets
load("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 2-RSS analysis/neoantigen_objects.RData")

#load required packages
#BiocManager::install("NPARC")
library(ggplot2)
library(cowplot)
library(NPARC)

#create datasets for different data
expression_data <- as.data.frame(exprs(peptides_novel_normalized) )#contains quantitive data. relative abundance
phenotype_data  <- pData(peptides_novel_normalized)  #contains temperature
feature_data<- featureData(peptides_novel_normalized)@data #contains Novel (T/F)
expression_data_t<-as.data.frame(t(expression_data))

# Extract novel and conventional peptides
novel_peptides <- rownames(feature_data)[feature_data$Novel]
head(novel_peptides)
conventional_peptides <- rownames(feature_data)[!feature_data$Novel]
head(conventional_peptides)


#Cmd+opt+E run from here
# Assign peptide to be plotted
peptide_name<-'peptide7777'
#novel peptides: peptide1542" "peptide1857" "peptide1983" "peptide2013" "peptide2702" "peptide2742....

# Redo analysis, create dataframe that includes necessary data
temperature_data <- data.frame(sample_id = rownames(phenotype_data), 
                               temperature = phenotype_data$temperature,
                               sample_name = phenotype_data$sample_name)  # Include sample_name
relAbundance_data <- data.frame(sample_id = rownames(expression_data_t), 
                                relAbundance = expression_data_t[[peptide_name]])
merged_data_sample <- merge(temperature_data, relAbundance_data, by = "sample_id")

merged_data_sample<-na.omit(merged_data_sample)

# Create the main plot
sample_plot_orig <- ggplot(merged_data_sample, aes(x = temperature, y = relAbundance)) +
  geom_point(color = 'red', alpha = 0.5) +
  theme_bw() +
  ggtitle(peptide_name)

print(sample_plot_orig)

#fit model 
nullFit_sample <- NPARC:::fitSingleSigmoid(x = merged_data_sample$temperature, y = merged_data_sample$relAbundance)
nullPredictions <- broom::augment(nullFit_sample)
merged_data_sample$nullPrediction <- nullPredictions$.fitted
merged_data_sample$nullResiduals <- nullPredictions$.resid
#get RSS
rss_sample <- sum(merged_data_sample$nullResiduals^2)

#add fitted model line and RSS value to the plot
sample_plot <- sample_plot_orig + 
  geom_line(data = merged_data_sample, aes(y = nullPrediction))+
  annotate("text", x = Inf, y = Inf, label = paste("RSS =", round(rss_sample, 2)), hjust = 1.1, vjust = 1.1)
print(sample_plot)

# Determine the type based on the Novel column
peptide_type <- ifelse(feature_data[peptide_name, "Novel"], "Novel", "Conventional")

# Create a text box with annotations
annotation_text <- paste(
  " Peptide sequence:", feature_data[peptide_name, "Peptide_AA"], "\n",
  "Gene ID:", feature_data[peptide_name, "id"], "\n",
  "Type:", peptide_type, "\n",
  "Cell lines count:", length(unique(merged_data_sample$sample_name))
)

# Add annotation text to the plot
sample_plot_anno <- sample_plot +
  labs(subtitle = annotation_text) 

# Print the annotated plot
print(sample_plot_anno)
print(paste("Cell lines:", peptide_name,paste(unique(merged_data_sample$sample_name), collapse = ", ")))


      