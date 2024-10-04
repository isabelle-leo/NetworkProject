library(dplyr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(broom)
library(knitr)
library(NPARC)
library(data.table)
library(car)

#load data
fits_novel<-readRDS('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/outputs/fits_novel.RDS')
fits_conv<-readRDS('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/outputs/fits_conv.RDS')

#load fitted model;
metrics_novel<-fits_novel$metrics %>% 
  mutate(modelType = factor(modelType), nCoeffs = factor(nCoeffs), nFitted = factor(nFitted), group = factor((group)))

predict_novel<-fits_novel$predictions %>% 
  mutate(modelType = factor(modelType), group = factor((group))) 

metrics_conv<-fits_conv$metrics %>% 
  mutate(modelType = factor(modelType), nCoeffs = factor(nCoeffs), nFitted = factor(nFitted), group = factor((group)))  

predict_conv<-fits_conv$predictions %>% 
  mutate(modelType = factor(modelType), group = factor((group))) 

#extract aumc
aumc_novel <- metrics_novel %>%
  filter(modelType == 'alternative') %>%
  group_by(id) %>%
  summarize(tm = mean(tm,na.rm = TRUE))

aumc_conv <- metrics_conv %>%
  filter(modelType == 'alternative') %>%
  group_by(id) %>%
  summarize(tm = mean(tm,na.rm = TRUE))

# Combine the data frames and add a distinguishing column in one step
combined_metrics <- rbind(
  cbind(aumc_novel, Type = "Novel"),
  cbind(aumc_conv, Type = "Conventional")
)


#na.omit(combined_metrics)
head(combined_metrics)

################# Initial Visualization & Statistical Test##########=

###plotting aumc
ggplot(combined_metrics, aes(x = tm, fill = Type)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of AUMC for Alternative Model Type",
       x = "AUMC",
       y = "Density") +
  theme_linedraw() 

# Remove top and bottom 5% extreme values
filtered_data <- combined_metrics %>%
  filter(between(tm, quantile(tm, 0.01, na.rm = TRUE), 
                 quantile(tm, 0.99, na.rm = TRUE)))

# Plot the density plot after filtering
ggplot(filtered_data, aes(x = tm, fill = Type)) +
  geom_density(alpha = 0.5) +
  labs(title = "Melting temperature distribution",
       x = "Melting Temperature",
       y = "Density") +
  theme_linedraw()

ggplot(combined_metrics, aes(y = log_tm, fill = Type)) +
  geom_boxplot(alpha = 0.5) +
  labs(title = "BoxPlot of AUMC for Alternative Model Type",
       x = "AUMC",
       y = "Density") +
  theme_linedraw() 

#summary
summary (aumc_novel$AUMC)
summary(aumc_conv$AUMC)

#normality check
shapiro.test(aumc_novel$AUMC)
#shapiro.test(aumc_conv$AUMC)

# QQ plot for visual inspection of normality
par(mfrow = c(1, 2)) # Create two plots side by side
qqnorm(aumc_novel$AUMC, main = "QQ Plot for Novel AUMC")
qqline(aumc_novel$AUMC, col = "red")
qqnorm(aumc_conv$AUMC, main = "QQ Plot for Conventional AUMC")
qqline(aumc_conv$AUMC, col = "red")

install.packages("nortest")
library(nortest)

# Anderson-Darling test for normality
ad.test(aumc_novel$AUMC) #normally distributed
ad.test(aumc_conv$AUMC) # not normally distributed

## statistical tests
wilcox.test(AUMC ~ Type, data = combined_metrics)#test median
leveneTest(AUMC ~ Type, data =  combined_metrics, center = median) #test variances
ks.test(AUMC ~ Type, data =  combined_metrics) #test distribution shape (if two group come from same distribution)

################ Add more information to the data, combine rss and aumc #############
aumc_data<-combined_metrics
rss_data<-read.delim("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4B-redo bio interpretation/inputs/aumc_rss_data.txt")
dim(rss_data)
dim(aumc_data)

aumc_rss_data <- rss_data %>%
  mutate(tm = aumc_data$tm[match(rss_data$id, aumc_data$id)])
 
 write.table(aumc_rss_data,'/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4B-redo bio interpretation/inputs/aumc_rss_data.txt',sep = "\t",row.names = TRUE)
 
 
 
 
