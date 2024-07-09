# STEP 1. load packages and datasets
#load datasets
load("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 2-RSS analysis/neoantigen_objects.RData")
#load required packages
BiocManager::install("NPARC")
# Load necessary packages
library(dplyr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(broom)
library(knitr)
library(NPARC)

# SETP 2. Data preprocessing and exploration
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

# 3. Function to fit null model and calculate RSS
calculate_rss <- function(peptide_name, phenotype_data, expression_data_t) {
  # Extract relevant data
  temperature_data <- data.frame(sample_id = rownames(phenotype_data), temperature = phenotype_data$temperature)
  relAbundance_data <- data.frame(sample_id = rownames(expression_data_t), relAbundance = expression_data_t[[peptide_name]])
  
  # Merge data frames by sample_id
  merged_data <- merge(temperature_data, relAbundance_data, by = "sample_id")
  
  # Remove missing values
  merged_data <- na.omit(merged_data)
  
  # Fit null model
  nullFit <- try(NPARC:::fitSingleSigmoid(x = merged_data$temperature, y = merged_data$relAbundance), silent = TRUE)
  
  if (inherits(nullFit, "try-error")) {
    return(NA)
  } else {
    # Calculate RSS
    nullPredictions <- broom::augment(nullFit)
    merged_data$nullResiduals <- nullPredictions$.resid
    rss <- sum(merged_data$nullResiduals^2)
    return(rss)
  }
}

# 4. Calculate RSS for novel and conventional peptides
novel_rss <- sapply(novel_peptides, calculate_rss, phenotype_data = phenotype_data, expression_data_t = expression_data_t)
conventional_rss <- sapply(conventional_peptides, calculate_rss, phenotype_data = phenotype_data, expression_data_t = expression_data_t)
#check how many peptides are there
length(novel_rss)
length(conventional_rss)

# 5. Create dataframe for RSS value and save
rss_data <- data.frame(
  Type = c(rep("Novel", length(novel_rss)), rep("Conventional", length(conventional_rss))),
  RSS = c(novel_rss, conventional_rss)
)

#save to local
#save(rss_data, file = "/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 2/Outputs/rss_data.RData")
write.table(rss_data, file = "/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 2/Outputs/rss_data_vsn2.txt", sep = "\t", row.names = TRUE)

#####extra: plotting melting curve######

# Assign peptide to be plotted
peptide_name<-'peptide2013'
#novel peptides: peptide1542" "peptide1857" "peptide1983" "peptide2013" "peptide2702" "peptide2742....


# Redo analysis, create dataframe that includes necessary data
temperature_data <- data.frame(sample_id = rownames(phenotype_data), 
                               temperature = phenotype_data$temperature,
                               sample_name = phenotype_data$sample_name)  # Include sample_name
relAbundance_data <- data.frame(sample_id = rownames(expression_data_t), 
                                relAbundance = expression_data_t[[peptide_name]])
merged_data_sample <- merge(temperature_data, relAbundance_data, by = "sample_id")

# Create the main plot
sample_plot_orig <- ggplot(merged_data_sample, aes(x = temperature, y = relAbundance)) +
  geom_point(color = 'red', alpha = 0.5) +
  theme_bw() +
  ggtitle(peptide_name, subtitle = "(novel)")

print(sample_plot_orig)

# Create a text box with annotations
annotation_text <- paste(
  "Peptide sequence:", feature_data[peptide_name, "Peptide_AA"], "\n",
  "Gene ID:", feature_data[peptide_name, "id"], "\n",
  "Cell lines:", paste(unique(merged_data_sample$sample_name), collapse = ", ")
)

annotation_plot <- ggdraw() + draw_label(annotation_text, x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5, size = 7)

# Combine the main plot with the text box below it
sample_plot_orig <- plot_grid(main_plot, annotation_plot, ncol = 1, rel_heights = c(3, 1))

# Print the plot
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

############################


##part II. validation for noise
#load saved data
rss_data<-read.delim("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 2/Outputs/rss_data.txt")
novel_rss<-subset(rss_data,Type=="Novel")$RSS
conventional_rss<-subset(rss_data,Type=="Conventional")$RSS

## plot the distribution
# Identify the MS1 area columns
ms1_cols <- grep('ms1_area', colnames(feature_data), value = TRUE)

# Create a data frame with the MS1 values
ms1_data <- feature_data[, ms1_cols]

# Apply log10 transformation to each column in ms1_data
log10_ms1_data <- log10(ms1_data)

library(reshape2)

# Create a long format data frame for ggplot
log10_ms1_data_long <- melt(log10_ms1_data, variable.name = "Set", value.name = "log10_MS1_Area")

# Plot distribution across sets using ggplot2
ggplot(log10_ms1_data_long, aes(x = log10_MS1_Area, color = Set)) +
  geom_density() +
  labs(title = "Density Plot of log10 MS1 Area for Each Set", x = "log10 MS1 Area", y = "Density") +
  theme_minimal()

##plot mean values
# Calculate mean values for each column in ms1_data
mean_values <- sapply(ms1_cols, function(col) mean(feature_data[[col]], na.rm = TRUE))

# Create a data frame for plotting
mean_values_df <- data.frame(Set = ms1_cols, MeanValue = mean_values)

# Plot using ggplot2
ggplot(mean_values_df, aes(x = Set, y = MeanValue, fill = Set)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Mean MS1 Area Values Across Sets", x = "Set", y = "Mean MS1 Area") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##################

### Part III. re-run analysis using low intensity conv peptides. save the result 
# Calculate the mean value for each peptide across all sets

mean_values <- rowMeans(ms1_data, na.rm = TRUE)

# Subset peptides with mean value < 10^6
# Calculate the mean value for each peptide across all sets
mean_values <- rowMeans(ms1_data, na.rm = TRUE)

# Subset peptides with mean value < 10^6
peptides_low_int <- rownames(ms1_data)[mean_values < 10^6]

#check
length(peptides_low_int)#1492

# Find the intersection of peptides_low_int and conventional_peptides
conventional_peptides_low <- intersect(peptides_low_int, conventional_peptides)
length(conventional_peptides_low) #1490

#re-run analysis
#calculate rss for new peptide list
conventional_rss_low <- sapply(conventional_peptides_low, calculate_rss, phenotype_data = phenotype_data, expression_data_t = expression_data_t)
rss_data_low <- data.frame(
  Type = c(rep("Novel", length(novel_rss)), rep("Conventional", length(conventional_rss_low))),
  RSS = c(novel_rss, conventional_rss_low)
)
#save to local
write.table(rss_data_low, file = "/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 2/Outputs/rss_data_low.txt", sep = "\t", row.names = TRUE)




##################

### Part 4. rerun analysis, using ....
peptides_table<-read.delim("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 2/Input/peptides_table.txt")
dim(peptides_table)
#[1] 223418    403
colnames(peptides_table)
