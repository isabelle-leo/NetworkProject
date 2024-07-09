# STEP 0. load file and package

# Load necessary libraries
library(ggplot2)
library(MASS)
library(car)
library(moments)
library(dplyr)

#load annotation file
info_data<-read.delim("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 3B-bio interpretation/inputs/CellAtlasAnnotationsVsSubCellBarCode.txt")

#read rss data
rss_data<-read.delim("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 3B-bio interpretation/inputs/rss_data_anno_new.txt")
#do not use na.omit!!!

#subset by type
rss_data_novel<-subset(rss_data,Type=="Novel")
rss_data_conventional<-subset(rss_data,Type=="Conventional")
head(rss_data_conventional)

### STEP 1. Load desired info data to conv dataset

#check info data
colnames(info_data)
#define columns to be assigned
columns_assign<-c( "Gene",
                   "Uniprot",
                  "nrLocations",
                  "Type",
                  "HPAbin",
                  "LOCATION",
                  "SCneigh",
                  "SC_ClassCount",
                  "Reliability")
 
# Perform an inner join to merge and keep only matching rows
rss_data_conventional_info <- rss_data_conventional %>%
  inner_join(info_data[columns_assign], by = c("Gene.ID" = "Gene")) %>%
  filter(!is.na(RSS)) #remove rows without RSS

head(rss_data_conventional_info)
dim(rss_data_conventional_info)

#log transform RSS
rss_data_conventional_info$RSS_log<-log(rss_data_conventional_info$RSS)
rss_data_novel$RSS_log<-log(rss_data_novel$RSS)
head(rss_data_conventional_info)
head(rss_data_novel)

### STEP 2 assign feature/group that you like to investigate
#set feature of interest

colnames(rss_data_conventional_info)
group_interest<-"nrLocations"  
category_counts <- table(rss_data_conventional_info[[group_interest]])
print(category_counts)

### STEP 3 summarize and visualize distribution

# Create multiple QQ plots for peptides in different location
qqplots <- generate_qqplot(rss_data_conventional_info, group = group_interest,feature='5')
print(qqplots)

# Density plots
density_plot_fill <- generate_density_plot_fill(rss_data_conventional_info, group_interest, "RSS_log")
print(density_plot_fill)

density_plot_line <-generate_density_plot_line(rss_data_conventional_info, group_interest, "RSS_log")
print(density_plot_line)

#for novel
ggplot() +
  geom_density(data = rss_data_novel %>% filter(!is.na(RSS_log)), aes(x = RSS_log, fill = "Novel"), alpha = 0.5) +
  labs(title = "RSS Novel", x = "RSS_log", y = "Density") +
  theme_light() #novel peak~2.7

# Creating a faceted density plot
ggplot(rss_data_conventional_info, aes(x = RSS_log, fill = group_interest)) +
  geom_density(alpha = 0.6) +  # Adjust transparency with alpha
  facet_wrap(~ LOCATION, scales = "fixed", ncol = 5) +  # Adjust the number of columns based on your display
  labs(title = "Faceted Density Plots of RSS_log by Subcellular Location",
       x = "RSS_log",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "none")+  # Hide the legend
  xlim(-5, 5)  # Set x-axis limits

# Boxplot
boxplot <- generate_boxplot(rss_data_conventional_info, group_interest, "RSS_log")
print(boxplot)

# Prepare the subsets
#novel group
novel_data <- filter(rss_data_novel, !is.na(RSS))
# Alternative group of interest, specifically filtering for '4_Mitochondria'
conventional_alter <- filter(rss_data_conventional_info, get(group_interest) == '5', !is.na(RSS_log))
# Null control group, including any data not categorized as '4_Mitochondria'
conventional_null <- filter(rss_data_conventional_info, get(group_interest) != '5', !is.na(RSS_log))

# 3-in-1 density plot
ggplot() +
  geom_density(data = conventional_null, aes(x = RSS_log, fill = "null"), alpha = 0.5) +
  geom_density(data = conventional_alter, aes(x = RSS_log, fill = "alter"), alpha = 0.5) +
  geom_density(data = novel_data, aes(x = RSS_log, fill = "Novel"), alpha = 0.5) +
  labs(title = paste("RSS_log Distribution by", group_interest),
       x = "RSS_log", y = "Density") +
  theme_light() +
  theme(legend.position = "top")


# Generate density plots
ggplot() +
  geom_density(data = rss_data_novel %>% filter(!is.na(RSS_log)), aes(x = RSS_log, fill = "Novel"), alpha = 0.5) +
  geom_density(data = rss_data_conventional_info %>% filter(nrLocations == 5, !is.na(RSS_log)), aes(x = RSS_log, fill = "Conventional,nrLocations=5"), alpha = 0.5) +
  labs(title = "Distribution of RSS_log, Novel v.s nrLocation=5", x = "RSS_log", y = "Density") +
  theme_light()+
  theme(legend.position = "top")

##QQ plots for Novel, test and control group
q# QQ plot for novel_data
qqnorm(novel_data$RSS_log, main = "QQ Plot for novel_data",)
qqline(novel_data$RSS_log, col = 'red')

# QQ plot for conventional_alter
qqnorm(conventional_alter$RSS_log, main = "QQ Plot for group of interest (nr=5)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(conventional_alter$RSS_log, col = 'red')

# QQ plot for conventional_null
qqnorm(conventional_null$RSS_log, main = "QQ Plot for nr<5", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(conventional_null$RSS_log, col = 'red')
### STEP 3. Statistical comparison/tests

# 3.1 Kruskal-Wallis rank sum test for more than 3 group
kruskal_test_result <- kruskal.test(as.formula(paste("RSS_log ~", group_interest)), data = rss_data_conventional_info)
print(kruskal_test_result)

# dunn test
library(dunn.test)
# Perform Dunn's test with/without Bonferroni correction
dunn_test_results <- dunn.test(x = rss_data_conventional_info$RSS, 
                               g = rss_data_conventional_info[[group_interest]], 
                               method = "bonferroni")
print(dunn_test_results)

# 3.2

# Wilcox test novel v.s. nrLocation=5
wilcox_results <- wilcox.test(novel_data$RSS_log, conventional_alter$RSS_log, alternative = "two.sided")
print(wilcox_results)

# Wilcox test  nrLocation=5 v.s. <5
wilcox_results <- wilcox.test(conventional_null$RSS_log, conventional_alter$RSS_log)
print(wilcox_results)

#Wilcox test Novel v.s. nrLocation<5
wilcox_results <- wilcox.test(novel_data$RSS_log, conventional_null$RSS_log)
print(wilcox_results)

## 3.3 Kolmogorov-Smirnov Tests of distribution shape

# KS test Novel vs. nrLocation=5
ks_results <- ks.test(novel_data$RSS_log, conventional_alter$RSS_log)
print(ks_results)

# KS test nrLocation=5 vs. nrLocation<5
ks_results <- ks.test(conventional_alter$RSS_log, conventional_null$RSS_log)
print(ks_results)

# KS test Novel vs. nrLocation<5
ks_results <- ks.test(novel_data$RSS_log, conventional_null$RSS_log)
print(ks_results)


#4.4 correlation
cor.test(rss_data_conventional_info$nrLocations,rss_data_conventional_info$RSS,method='spearman')
