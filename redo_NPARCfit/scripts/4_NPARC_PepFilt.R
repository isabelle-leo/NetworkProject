#1. Load # Load necessary packages
library(dplyr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(broom)
library(knitr)
library(NPARC)
library(data.table)

##load data
peptide_data<-read.delim('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/inputs/NPARC_data_PepFilter.txt')
#peptide_data<-fread('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/inputs/NPARC_data_PepFilter.txt')

# 2. pre process data

# Rename the column from Novel to Type
peptide_data <- peptide_data %>%
  rename(Type = Novel)
# Replace TRUE with "Novel" and FALSE with "Conventional"
peptide_data$Type <- ifelse(peptide_data$Type, "Novel", "Conventional")
head(peptide_data,n=3)

peptide_data %>% 
  mutate(Type = factor(Type), 
         sample_name = factor(sample_name), 
         cell_line_count = factor(cell_line_count)) 
#######
df<-peptide_data
colnames(df)
unique(df$sample_name)
#df %<>% filter(cell_line_count>= 5) #novel=37, cov=14692
table(df$Type)

#add replicate sample name 
df$sample_name2<- gsub("_BR2", "", df$sample_name)

#remove columns with NA in abundance
df <- df[!is.na(df$relAbundance), ]

# Count full curves per peptide
df %<>%
  group_by( Peptide_AA) %>%
  mutate(n = n()) %>%
  group_by(sample_name) %>%
  mutate(max_n = max(n)) %>% 
  ungroup
table(distinct(df, Peptide_AA, n)$n)


# Subset df to df_novel and df_conv based on Type using dplyr
df_novel <- df %>%
  filter(Type == "Novel")

df_conv <- df %>%
  filter(Type == "Conventional")


df_conv_subset <- df_conv[df_conv$peptide_id %in% sample(unique(df_conv$peptide_id), 1500), ]
df_conv_subset<-df_conv

# 3. perform NPARC fitting

#novel peptides
BPPARAM <- BiocParallel::SerialParam(progressbar = FALSE)
#BPPARAM <- BiocParallel::SerialParam(progressbar = TRUE)


fits_novel <- NPARCfit(x = df_novel$temperature, 
                 y = df_novel$relAbundance, 
                 id = df_novel$peptide_id, 
                 groupsNull = NULL,
                 groupsAlt = df_novel$sample_name, 
                 BPPARAM = BPPARAM,
                 returnModels = FALSE)
saveRDS(fits_novel,'/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/outputs/fits_novel.RDS')

#conventional
fits_conv <- NPARCfit(x = df_conv_subset$temperature, 
                       y = df_conv_subset$relAbundance, 
                       id = df_conv_subset$peptide_id, 
                       groupsNull = NULL, 
                       groupsAlt = df_conv_subset$sample_name, 
                       BPPARAM = BPPARAM,
                       returnModels = FALSE)
#saveRDS(fits_conv,'/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/outputs/fits_conv.RDS')
################# Model Fitting Finished ###############


#load fitted model;

metrics_novel<-fits_novel$metrics %>% 
  mutate(modelType = factor(modelType), nCoeffs = factor(nCoeffs), nFitted = factor(nFitted), group = factor((group)))

predict_novel<-fits_novel$predictions %>% 
  mutate(modelType = factor(modelType), group = factor((group))) 

metrics_conv<-fits_conv$metrics %>% 
  mutate(modelType = factor(modelType), nCoeffs = factor(nCoeffs), nFitted = factor(nFitted), group = factor((group)))  
  
predict_conv<-fits_conv$predictions %>% 
  mutate(modelType = factor(modelType), group = factor((group))) 

#####calculate rss
rss_novel <- metrics_novel %>%
  filter(modelType == 'alternative') %>%
  group_by(id) %>%
  summarize(RSS = mean(rss,na.rm = TRUE))

rss_conv <- metrics_conv %>%
  filter(modelType == 'alternative') %>%
  group_by(id) %>%
  summarize(RSS = mean(rss,na.rm = TRUE))

rss_novel$RSS_log<-log(rss_novel$RSS)
rss_conv$RSS_log<-log(rss_conv$RSS)

# Display the resulting table
print(rss_novel)
print(rss_conv)

#4. process result data

# Combine the data frames and add a distinguishing column in one step
combined_metrics <- rbind(
  cbind(rss_novel, group = "Novel"),
  cbind(rss_conv, group = "Conventional")
)

#na.omit(combined_metrics)
head(combined_metrics)


#add necessary data, make new rss table
rss_old<-read.delim("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 3B-bio interpretation/inputs/rss_data_full.txt")
head(combined_metrics)
head(rss_old)
dim(combined_metrics)# 51777     4

# Filter rss_old based on rownames being present in combined_metrics$id
rss_filtered <- rss_old[rownames(rss_old) %in% combined_metrics$id, ]
# Remove the RSS column and add RSS and RSS_log from combined_metrics
rss_filtered$RSS <- NULL
# Ensure that the row names in rss_filtered match the ids in combined_metrics
rss_filtered <- rss_filtered %>%
  mutate(RSS = combined_metrics$RSS[match(rownames(rss_filtered), combined_metrics$id)],
         RSS_log = combined_metrics$RSS_log[match(rownames(rss_filtered), combined_metrics$id)])
head(rss_filtered)
head(combined_metrics)
table(rss_filtered$Type)

## investigate the effect of cell line count
#plot cell line count
generate_boxplot(rss_filtered,'cell_line_count','RSS_log')
generate_density_plot_line(rss_filtered,'cell_line_count','RSS_log')
#generate_boxplot(rss_filtered[rss_filtered$Type == 'Novel', ], 'cell_line_count', 'RSS_log')
#generate_boxplot(rss_filtered[rss_filtered$Type == 'Conventional', ], 'cell_line_count', 'RSS_log')

cor.test(rss_filtered$RSS, rss_filtered$cell_line_count, method = "spearman")

#write.table(rss_filtered,'/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4B-redo bio interpretation/inputs/rss_data_mean.txt',sep = "\t",row.names = TRUE)
#saveRDS(combined_metrics,'/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/outputs/model_result_4')


# 5. Plotting
# Generate the combined density plot
ggplot(rss_filtered, aes(x = RSS_log, fill = Type)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of RSS for Alternative Model Type",
       x = "RSS (log scale)",
       y = "Density") +
  theme_linedraw() 


ggplot(combined_metrics, aes(y = RSS_log, fill = group)) +
  geom_boxplot(alpha = 0.5) +
  labs(title = "Density Plot of RSS for Alternative Model Type",
       y = "RSS (log scale)",
       x = "Density") +
  theme_linedraw() 

#plotting novel peptides with more cell line count
ggplot(rss_filtered %>% 
         filter((Type == "Novel" & cell_line_count >= 4) | Type == "Conventional"), 
       aes(x = RSS_log, fill = Type)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of RSS for Alternative Model Type",
       x = "RSS (log scale)",
       y = "Density") +
  theme_linedraw()

ggplot(rss_filtered %>% 
         filter((Type == "Novel" & cell_line_count >= 4) | Type == "Conventional"), 
       aes(y = RSS_log, fill = Type)) +
  geom_boxplot(alpha = 0.5) +
  labs(title = "BoxPlot of RSS for Alternative Model Type",
       x = "RSS (log scale)",
       y = "Density") +
  theme_linedraw()
#6. stat test

#summary 
summary (rss_novel$RSS_log)
summary(rss_conv$RSS_log)

wilcox.test(RSS ~ group, data = combined_metrics)
print(ks.test(combined_metrics$RSS_log[combined_metrics$group == "Novel"], 
              combined_metrics$RSS_log[combined_metrics$group == "Conventional"]))
kruskal.test(RSS ~ group, data = combined_metrics)
#leveneTest(RSS ~ group, data = combined_metrics)
print(fligner.test(RSS ~ group, data = combined_metrics))


# 7. example curve
#plotting functions
plot_peptide <- function(peptide_name) {
  # Filter the data for the selected peptide
  merged_data_sample <- df_novel[df_novel$peptide_id == peptide_name, ]
  
  # Ensure sample_name is a factor
  merged_data_sample$sample_name <- as.factor(merged_data_sample$sample_name)
  
  # Create the main plot
  sample_plot_orig_2 <- ggplot(merged_data_sample, aes(x = temperature, y = relAbundance)) +
    geom_point(aes(color = sample_name)) +
    theme_bw() +
    ggtitle(paste(peptide_name, "(Type=Novel)")) +
    scale_color_discrete(name = "Sample Name")
  
  # Filter predictions and metrics for the selected peptide
  stk4Predictions <- filter(fits_novel$predictions, id == peptide_name)
  stk4Metrics <- filter(fits_novel$metrics, id == peptide_name)
  
  # Calculate RSS values
  rssNull <- filter(stk4Metrics, modelType == "null")$rss
  rssAlt <- mean(filter(stk4Metrics, modelType == "alternative")$rss) 
  
  # Add fitted lines and RSS annotations to the plot
  sample_plot_final <- sample_plot_orig_2 +
    geom_line(data = filter(stk4Predictions, modelType == "alternative"), 
              aes(x = x, y = .fitted, color = factor(group))) +
    geom_line(data = filter(stk4Predictions, modelType == "null"), 
              aes(x = x, y = .fitted)) +
    annotate("text", x = Inf, y = Inf, label = paste("RSS null =", signif(rssNull, 3)), hjust = 1.1, vjust = 1.1) +
    annotate("text", x = Inf, y = Inf, label = paste("RSS alter =", signif(rssAlt, 3)), hjust = 1.1, vjust = 3)
  # Print the plot
  print(sample_plot_final)
}

# Example usage:
plot_peptide("peptide146384")
plot_peptide("peptide144295")
plot_peptide("peptide145732")
plot_peptide("peptide146601")
plot_peptide("peptide147758")

