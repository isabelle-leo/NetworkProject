# STEP 0. load file and package
# Load necessary libraries
library(ggplot2)
library(MASS)
library(car)
library(moments)
library(dplyr)
library(data.table)
#load annotation file

info_data<-read.delim('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 3B-bio interpretation/inputs/DisorderAnnotation.txt')

#read rss data
rss_data<-read.delim("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 3B-bio interpretation/inputs/rss_data_full.txt")
#do not use na.omit!!!

#subset by type
rss_data_novel<-subset(rss_data,Type=="Novel")
rss_data_conventional<-subset(rss_data,Type=="Conventional")
head(rss_data_conventional)

### STEP 1. Load desired info data to conv dataset &assign disorder lv
rss_data_conventional_info  <- merge(rss_data_conventional, info_data, by.x = "Gene.ID", by.y = "hgnc_symbol")
#filter out na in RSS
rss_data_conventional_info <- rss_data_conventional_info [!is.na(rss_data_conventional_info $RSS), ]

#assign disorder level based on value

# Assign 'disorder_lv' based on 'prediction_disorder_alphafold' values
rss_data_conventional_info$disorder_lv_alphafold <- with(rss_data_conventional_info, 
                                                         ifelse(prediction_disorder_alphafold < 0.1, "structured", 
                                                                ifelse(prediction_disorder_alphafold >= 0.5, "disordered", "medium")))

# Assign 'disorder_lv' based on 'prediction_disorder_mobidb_lite' values
rss_data_conventional_info$disorder_lv_mobidb <- with(rss_data_conventional_info, 
                                                      ifelse(prediction_disorder_mobidb_lite < 0.1, "structured", 
                                                             ifelse(prediction_disorder_mobidb_lite >= 0.5, "disordered", "medium")))

# Print the head of the updated data frame
print(head(rss_data_conventional_info))

#log transform RSS
rss_data_conventional_info$RSS_log<-log(rss_data_conventional_info$RSS)
rss_data_novel$RSS_log<-log(rss_data_novel$RSS)
head(rss_data_conventional_info)
head(rss_data_novel)    

#make quantile bin
# Create quantile bins based on 'prediction_disorder_mobidb_lite'
rss_data_conventional_info$disorder_quantile <- cut(
  rss_data_conventional_info$prediction_disorder_mobidb_lite, 
  breaks = quantile(rss_data_conventional_info$prediction_disorder_mobidb_lite, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("Q1", "Q2", "Q3", "Q4")
)

### STEP 2 assign feature/group that you like to investigate
#set feature of interest

colnames(rss_data_conventional_info)
group_interest<-"disorder_quantile"  
category_counts <- table(rss_data_conventional_info[[group_interest]])
print(category_counts)


### STEP 3 summarize and visualize distribution

#histogram of distribution of disorder rate 
hist(rss_data_conventional_info$prediction_disorder_mobidb_lite,breaks = 50,main = "Histogram of MobiDB-lite disorder",xlab = 'MobiDB-lite disorder')
hist(rss_data_conventional_info$prediction_disorder_alphafold,breaks = 50,,main = "Histogram of AlphaFold-disorder",ylab = 'MobiDB-lite disorder')

# Create multiple QQ plots for peptides in different location
#qqplots <- generate_qqplot(rss_data_conventional_info, group = group_interest,feature='5')
#print(qqplots)

# Density plots
density_plot_fill <- generate_density_plot_fill(rss_data_conventional_info, group_interest, "RSS_log")
print(density_plot_fill)

density_plot_line <-generate_density_plot_line(rss_data_conventional_info, group_interest, "RSS_log")
print(density_plot_line)

#for novel
ggplot() +
  geom_density(data = rss_data_novel %>% filter(!is.na(RSS_log)), aes(x = RSS_log, fill = "Novel"), alpha = 0.5) +
  labs(title = "RSS Novel", x = "RSS_log", y = "Density") +
  theme_light() #novel peak~0.27, -0.8-1

#boxplot
boxplot <- generate_boxplot(rss_data_conventional_info, group_interest, "RSS_log")
print(boxplot)


####Prepare the subsets!
colnames(rss_data_conventional_info)
#novel group
novel_data <- filter(rss_data_novel, !is.na(RSS))
# Alternative group of interest, specifically filtering for '4_Mitochondria'
conventional_alter <- filter(rss_data_conventional_info, get(group_interest) =='Q4', !is.na(RSS_log))
# Null control group, including any data not categorized as '4_Mitochondria'
conventional_null <- filter(rss_data_conventional_info, get(group_interest) =='Q1', !is.na(RSS_log))

# 3-in-1 density plot
ggplot() +
  geom_density(data = conventional_null, aes(x = RSS_log, fill = "null"), alpha = 0.5) +
  geom_density(data = conventional_alter, aes(x = RSS_log, fill = "alter"), alpha = 0.5) +
  geom_density(data = novel_data, aes(x = RSS_log, fill = "Novel"), alpha = 0.5) +
  labs(title = paste("RSS_log Distribution by", group_interest),
       x = "RSS_log", y = "Density") +
  theme_light() +
  theme(legend.position = "top")

### STEP 4. Statistical comparison/tests

# 4.1 Kruskal-Wallis rank sum test for more than 3 group
kruskal_test_result <- kruskal.test(as.formula(paste("RSS ~", group_interest)), data = rss_data_conventional_info)
print(kruskal_test_result)

# dunn test
library(dunn.test)
# Perform Dunn's test with/without Bonferroni correction
#dunn_test_results <- dunn.test(x = rss_data_conventional_info$RSS, 
                              # g = rss_data_conventional_info[[group_interest]], 
                              # method = "bonferroni")
#print(dunn_test_results)

# 4.2

# Wilcox test novel v.s. alter
wilcox.test(novel_data$RSS_log, conventional_alter$RSS_log, alternative = "two.sided")

# Wilcox test  alter v.s. null
wilcox.test(conventional_null$RSS_log, conventional_alter$RSS_log)

#Wilcox test Novel v.s. null
wilcox.test(novel_data$RSS_log, conventional_null$RSS_log)

## 4.3 Kolmogorov-Smirnov Tests of distribution shape

# KS test Novel vs. alter
ks_results <- ks.test(novel_data$RSS_log, conventional_alter$RSS_log)

# KS test alter v.s. null
ks_results <- ks.test(conventional_alter$RSS_log, conventional_null$RSS_log)

# KS test Novel vs. null
ks_results <- ks.test(novel_data$RSS_log, conventional_null$RSS_log)

##4.4 Extra tests:

#correlation analysis (spearman non-parametric)
cor.test(rss_data_conventional_info$RSS, rss_data_conventional_info$,method='spearman')
