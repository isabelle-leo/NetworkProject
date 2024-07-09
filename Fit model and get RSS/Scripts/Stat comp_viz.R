### STEP 0. Load data and package
#read rss data
rss_data<-read.delim("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 2/Outputs/rss_data.txt")
#remove na
rss_data<-na.omit(rss_data)
#assign rss for novel/conv peptides to list
novel_rss<-subset(rss_data,Type=="Novel")$RSS
conventional_rss<-subset(rss_data,Type=="Conventional")$RSS

#check the data
summary(rss_data)

# Load necessary libraries
library(ggplot2)
library(MASS)
library(car)
library(moments)

### STEP 1. log transform data and and assign to new dataset 
# Log transformation
log_novel_rss <- log(novel_rss)
log_conventional_rss <- log(conventional_rss)

#assign transformed value to new data frames
rss_data_log <- data.frame(
  RSS = c(log_novel_rss, log_conventional_rss),
  Type = factor(c(rep("Novel", length(log_novel_rss)), rep("Conventional", length(log_conventional_rss))))
)

### STEP 2 summarize and visualize distribution
#summary

# Get summary statistics for Novel and Conventional peptides
summary(rss_data[rss_data$Type == "Novel", "RSS"])
summary(rss_data[rss_data$Type == "Conventional", "RSS"])

#density plot for both type
ggplot(rss_data_log, aes(x = RSS, fill = Type)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Log-transformed RSS Values",
       x = "Log-transformed RSS",
       y = "Density")+
  theme_minimal()

# boxplots
ggplot(rss_data_log, aes(x = Type, y = RSS, fill = Type)) +
  geom_boxplot() +
  labs(title = "Boxplot of Log-transformed RSS Values",
       x = "Peptide Type",
       y = "Log-transformed RSS")+
  theme_minimal()

### STEP 3. Statistical comparison/tests

## 3.1 Check stat summary & normality
summary(novel_rss)
length(novel_rss)
summary(conventional_rss)
length(conventional_rss)

##check normarlity (can be skipped)

#shapiro.test(log_novel_rss)
#shapiro.test(log_conventional_rss)#can't because sample size too large

# Q-Q plot for log_conventional_rss
#qqnorm(log_novel_rss,main = "Normal Q-Q plot: log_novel_rss")
#qqline(log_novel_rss, col = "red")


# Q-Q plot for log_conventional_rss
#qqnorm(log_conventional_rss, main = "Normal Q-Q plot: log_conv_rss")
#qqline(log_conventional_rss, col = "red")

## 3.2 Statistic tests

#If sample normally distributed, use t-test (it's not)
t.test(novel_rss,conventional_rss)#p=0.268

#If the samples are not normally distributed, use wilcox.test (transformation does not affect result)
wilcox.test(novel_rss, conventional_rss, alternative = "two.sided")
#p=0.1237

#Test variance 
leveneTest(RSS ~ Type, data = rss_data)#p=0.9432

##Calculate kurtosis for each type
kurtosis_novel <- kurtosis(novel_rss, na.rm = TRUE)
kurtosis_conventional <- kurtosis(conventional_rss, na.rm = TRUE)
# Display the kurtosis values
cat("Kurtosis for Novel RSS: ", kurtosis_novel, "\n")
cat("Kurtosis for Conventional RSS: ", kurtosis_conventional, "\n")

## Kolmogorov-Smirnov Tests of distribution shape
ks.test(novel_rss, conventional_rss)


############## END ##################

### extra:
#try different transform method and visulize
# Remove NAs from novel_rss
novel_rss <- novel_rss[!is.na(novel_rss)]

# Log transformation
log_novel_rss <- log(novel_rss)

# Square root transformation
sqrt_novel_rss <- sqrt(novel_rss)

# Reciprocal transformation
reciprocal_novel_rss <- 1 / (novel_rss + 1)

# Create a dataframe for plotting
transformed_data <- data.frame(
  RSS = c(novel_rss, log_novel_rss, sqrt_novel_rss, reciprocal_novel_rss),
  Transformation = rep(c("Original", "Log", "Square Root", "Reciprocal"), each = length(novel_rss))
)

# Plotting the distributions
ggplot(transformed_data, aes(x = RSS, fill = Transformation)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Transformation, scales = "free") +
  labs(title = "Density Plots for Different Transformations of Novel RSS (vsn)",
       x = "Transformed RSS",
       y = "Density") +
  theme_minimal()


##plot conventional
# Remove NAs from novel_rss
conventional_rss <- conventional_rss[!is.na(conventional_rss)]

# Log transformation
log_conventional_rss <- log(conventional_rss)

# Square root transformation
sqrt_coventional_rss <- sqrt(conventional_rss)

# Reciprocal transformation
reciprocal_conventional_rss <- 1 / (conventional_rss + 1)

# Create a dataframe for plotting
transformed_data_conv <- data.frame(
  RSS = c(conventional_rss, log_conventional_rss, sqrt_coventional_rss, reciprocal_conventional_rss),
  Transformation = rep(c("Original", "Log", "Square Root", "Reciprocal"), each = length(conventional_rss))
)

# Plotting the distributions
ggplot(transformed_data_conv, aes(x = RSS, fill = Transformation)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Transformation, scales = "free") +
  labs(title = "Density Plots for Different Transformations of Conv RSS (vsn) ",
       x = "Transformed RSS",
       y = "Density") +
  theme_minimal()




