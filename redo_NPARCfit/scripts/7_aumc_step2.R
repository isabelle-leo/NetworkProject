library(dplyr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(broom)
library(knitr)
library(NPARC)
library(data.table)

#load data
aumc_rss_data<-read.delim('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4B-redo bio interpretation/inputs/aumc_rss_data.txt')
head(aumc_rss_data,n=3)

#1. plotting novel peptides with more cell line count
table(aumc_rss_data$Type)
ggplot(aumc_rss_data, aes(x = AUMC, fill = Type)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of AUMC for Alternative Model Type",
       x = "AUMC",
       y = "Density") +
  theme_bw() 

#n=72 for cell line>=4, n=37 for cell line>=6, n=21 for cell line >=8
ggplot(aumc_rss_data %>% 
         filter((Type == "Novel" & cell_line_count >= 4) | Type == "Conventional"), 
       aes(x = AUMC, fill = Type)) +
  geom_density(alpha = 0.5) +
  labs(title = "AUMC Distribution, Novel cell line count >= 4",
       x = "AUMC",
       y = "Density") +
  theme_bw()

ggplot(aumc_rss_data %>% 
         filter((Type == "Novel" & cell_line_count >= 6) | Type == "Conventional"), 
       aes(x = AUMC, fill = Type)) +
  geom_density(alpha = 0.5) +
  labs(title = "AUMC Distribution, Novel cell line count >= 6",
       x = "AUMC",
       y = "Density") +
  theme_bw()

#calculating sample size
aumc_rss_data %>%
  filter(Type == "Novel" & cell_line_count >= 8) %>%
  nrow()

#All-in- one
ggplot(aumc_rss_data %>%
         filter(Type == "Novel") %>%
         mutate(group = case_when(cell_line_count >= 6 ~ "Novel, cell line >= 6",
                                  cell_line_count >= 4 ~ "Novel, cell line >= 4",
                                  cell_line_count >= 0 ~ "Novel, all")),
       aes(x = AUMC,  color= group)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of AUMC for Novel Peptides by Cell Line Count",
       x = "AUMC",
       y = "Density") +
  theme_bw()

#stat test
# Prepare the data with the group variable
kruskal_data <- aumc_rss_data %>%
  filter(Type == "Novel") %>%
  mutate(group = case_when(
    cell_line_count >= 6 ~ "Novel, cell line >= 6",
    cell_line_count >= 4 ~ "Novel, cell line >= 4",
    TRUE ~ "All Novel"
  ))

# Perform the Kruskal-Wallis test
kruskal.test(AUMC ~ group, data = kruskal_data)


#summary
summary(aumc_rss_data %>%
          filter(Type == "Novel" & cell_line_count >= 6) %>%
          select(AUMC, RSS_log))
summary(aumc_rss_data %>%
          filter(Type == "Conventional" & cell_line_count >= 0) %>%
          select(AUMC, RSS_log))


#2. effect of cell line count on aumc
generate_boxplot(aumc_rss_data,'cell_line_count','AUMC')

# 3. Relationship between RSS and AUMC
ggplot(aumc_rss_data, aes(y = RSS_log, x = AUMC)) +
  geom_point(alpha = 0.2,color='purple') +
  geom_smooth(method = 'glm', se=TRUE, color = "black") +
  labs(title = "Relationship between AUMC and RSS (log scale)",
       y = "RSS (log scale)",
       x = "AUMC") +
  theme_bw()

# Load necessary package
library(gridExtra)
# Novel plot with rho, axes reversed
plot_novel <- ggplot(aumc_rss_data %>% filter(Type == "Novel"), aes(y = RSS_log, x = AUMC)) +
  geom_point(alpha = 0.3, color = "deepskyblue3") +
  geom_smooth(method = 'glm', se = TRUE, color = "black") +
  labs(title = "RSS (log) vs. AUMC - Novel",
       x = "AUMC",
       y = "RSS (log scale)") +
  annotate("text", x = Inf, y = Inf, label = "rho = 0.3168", hjust = 1.1, vjust = 1.2, size = 4) +
  theme_bw()

# Conventional plot with rho, axes reversed
plot_conv <- ggplot(aumc_rss_data %>% filter(Type == "Conventional"), aes(y = RSS_log, x = AUMC)) +
  geom_point(alpha = 0.3, color = "coral2") +
  geom_smooth(method = 'glm', se = TRUE, color = "black") +
  labs(title = "RSS (log) vs. AUMC - Conventional",
       x = "AUMC",
       y = "RSS (log scale)") +
  annotate("text", x = Inf, y = Inf, label = "rho = 0.5008", hjust = 1.1, vjust = 1.2, size = 4) +
  theme_bw()
# Display side by side
grid.arrange(plot_novel, plot_conv, ncol = 2)


#correlation test
cor.test(aumc_rss_data$AUMC,aumc_rss_data$RSS_log, method = "spearman")#for all,0.5003485 ,0.480798 
#for novel
cor.test(aumc_rss_data %>% filter(Type == "Novel") %>% pull(RSS_log),
         aumc_rss_data %>% filter(Type == "Novel") %>% pull(AUMC), 
         method = "spearman",na.rm=TRUE) #rho=0.3168193 ,p-value = 0.000125
sum(!is.na(aumc_rss_data %>% filter(Type == "Novel") %>% pull(RSS_log)) & !is.na(aumc_rss_data %>% filter(Type == "Novel") %>% pull(AUMC)))
#for conv
cor.test(aumc_rss_data %>% filter(Type == "Conventional") %>% pull(RSS_log),
         aumc_rss_data %>% filter(Type == "Conventional") %>% pull(AUMC), 
         method = "spearman",na.rm=TRUE) #rho=0.5007541 
sum(!is.na(aumc_rss_data %>% filter(Type == "Conventional") %>% pull(RSS_log)) & !is.na(aumc_rss_data %>% filter(Type == "Conventional") %>% pull(AUMC)))

#comparing coefficient using z difference
zdifference <- function(r1, r2, n1, n2){
  zd <- (atanh(r1)-atanh(r2))/sqrt(1/(n1-3)+1/(n2-3))
  p <- 1 - pnorm(abs(zd))
  print(paste("Z Difference: ", zd))
  print(paste("One-Tailed P-Value: ", p))
  print(paste("Two-Tailed P-Value: ", (2*p)))
}

zdifference(0.3168193 , 0.5007541 , 143, 51632)

