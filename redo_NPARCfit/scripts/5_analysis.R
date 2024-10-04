library(dplyr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(broom)
library(knitr)
library(NPARC)
library(data.table)

#1. load data
combined_metrics<-read.delim('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4B-redo bio interpretation/inputs/rss_data_mean.txt')
#load mdoel
fits_conv<-readRDS('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/outputs/fits_conv.RDS')
rss_filtered<-combined_metrics

#2. plotting
# Generate the combined density plot
ggplot(combined_metrics, aes(x = RSS_log, fill = Type)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of RSS for Alternative Model Type",
       x = "RSS (log scale)",
       y = "Density") +
  theme_light() 


ggplot(combined_metrics, aes(y = RSS_log, fill = Type)) +
  geom_boxplot(alpha = 0.5) +
  labs(title = "BoxPlot of RSS for Alternative Model Type",
       x = "RSS (log scale)",
       y = "Density") +
  theme_light() 

#plotting novel peptides with more cell line count
ggplot(rss_filtered %>% 
         filter((Type == "Novel" & cell_line_count >= 6) | Type == "Conventional"), 
       aes(x = RSS_log, fill = Type)) +
  geom_density(alpha = 0.5) +
  labs(title = "RSS distribution, Novel cell line count>=6",
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

#3. stat test

#summary 
table(rss_filtered$Type)
summary (rss_filtered$RSS[combined_metrics$Type == "Novel"])
summary (rss_filtered$RSS[combined_metrics$Type == "Conventional"])

#median/mean
wilcox.test(RSS ~ Type, data = combined_metrics)
#shape
print(ks.test(combined_metrics$RSS[combined_metrics$Type == "Novel"], 
              combined_metrics$RSS[combined_metrics$Type == "Conventional"]))

#variance
car::leveneTest(RSS_log ~ Type, data = combined_metrics)
#kruskal.test(RSS ~ Type, data = combined_metrics)
print(fligner.test(RSS_log ~ Type, data = combined_metrics))

#### subset data
summary(rss_filtered$RSS[rss_filtered$Type == "Novel" & rss_filtered$cell_line_count >= 6])

wilcox.test(
  rss_filtered %>% filter(Type == "Novel" & cell_line_count >= 4) %>% pull(RSS_log),
  rss_filtered %>% filter(Type == "Conventional") %>% pull(RSS_log),
  alternative = "two.sided"
)

ks.test(
  rss_filtered %>% filter(Type == "Novel" & cell_line_count >= 4) %>% pull(RSS_log),
  rss_filtered %>% filter(Type == "Conventional") %>% pull(RSS_log),
  alternative = "two.sided"
)

#variance
car::leveneTest(RSS_log ~ Type, data = rss_filtered %>% filter((Type == "Novel" & cell_line_count >= 4) | Type == "Conventional"))

#kruskal.test(RSS ~ Type, data = combined_metrics)
print(fligner.test(RSS_log ~ Type, data = combined_metrics))



library(ggplot2)

# Remove outliers using the 1st and 99th percentiles, ignoring NA values
filtered_data <- rss_filtered %>%
  filter(between(RSS_log, quantile(RSS_log, 0.01, na.rm = TRUE), 
                          quantile(RSS_log, 0.99, na.rm = TRUE)))

# Create density plot
ggplot(filtered_data, aes(x = RSS_log, fill = Type)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot (Excluding Extreme Values 0.01)",
       x = "RSS (log scale)",
       y = "Density") +
  theme_linedraw()

ggplot(filtered_data %>% 
         filter((Type == "Novel" & cell_line_count >= 4) | Type == "Conventional"), 
       aes(x = RSS_log, fill = Type)) +
  geom_density(alpha = 0.5) +
  labs(title = "RSS distribution, Novel cell line count>=4(outlier removed)",
       subtitle = "p=0.068",
       x = "RSS (log scale)",
       y = "Density") +
  theme_linedraw()
library(ggplot2)

# Remove outliers using IQR method
# Remove outliers using IQR and quantiles, ignoring NA values
filtered_data <- rss_filtered %>%
  filter(between(RSS_log, 
                 quantile(RSS_log, 0.25, na.rm = TRUE) - 1.5 * IQR(RSS_log, na.rm = TRUE),
                 quantile(RSS_log, 0.75, na.rm = TRUE) + 1.5 * IQR(RSS_log, na.rm = TRUE)))
# Create density plot
ggplot(filtered_data, aes(y = RSS_log, fill = Type)) +
  geom_boxplot(alpha = 0.5) +
  labs(title = "Density Plot (Excluding Extreme Values)",
       x = "RSS (log scale)",
       y = "Density") +
  theme_linedraw()

