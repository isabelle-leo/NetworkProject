
aumc_rss_data<-read.delim('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4B-redo bio interpretation/inputs/aumc_rss_data.txt')
head(aumc_rss_data,n=3)

# Remove top and bottom 5% extreme values
aumc_rss_data <- aumc_rss_data %>%
  filter(between(tm, quantile(tm, 0.01, na.rm = TRUE), 
                 quantile(tm, 0.99, na.rm = TRUE)))

ggplot(aumc_rss_data, aes(y = AUMC, x = tm)) +
  geom_point(alpha = 0.2,color='purple') +
  geom_smooth(method = , se=TRUE, color = "black") +
  labs(title = "Relationship between AUMC and tm ",
       x = "tm",
       y = "AUMC") +
  theme_bw()

ggplot(aumc_rss_data, aes(y = RSS_log, x = tm)) +
  geom_point(alpha = 0.2,color='purple') +
  geom_smooth(method = , se=TRUE, color = "black") +
  labs(title = "Relationship between AUMC and tm ",
       x = "tm",
       y = "RSS_log") +
  theme_bw()

cor.test(aumc_rss_data$AUMC,aumc_rss_data$tm, method = "spearman")
sum(!is.na(aumc_rss_data$AUMC) & !is.na(aumc_rss_data$tm))

library(gridExtra)
plot_novel <- ggplot(aumc_rss_data %>% filter(Type == "Novel"), aes(x = tm, y = AUMC)) +
  geom_point(alpha = 0.3, color = "deepskyblue3") +
  geom_smooth(method = 'glm', se = TRUE, color = "black") +
  labs(title = "AUMC v.s. tm-Novel",
       x = "tm",
       y = "AUMC") +
  annotate("text", x = Inf, y = Inf, label = "rho = 0.924", hjust = 1.1, vjust = 1.2, size = 4) +
  theme_bw()

# Conventional plot with rho, axes reversed
plot_conv <- ggplot(aumc_rss_data %>% filter(Type == "Conventional"), aes(x = tm, y = AUMC)) +
  geom_point(alpha = 0.3, color = "coral2") +
  geom_smooth(method = 'glm', se = TRUE, color = "black") +
  labs(title = "AUMC v.s. tm- Conventional",
       x = "tm",
       y = "AUMC") +
  annotate("text", x = Inf, y = Inf, label = "rho = 0.952", hjust = 1.1, vjust = 1.2, size = 4) +
  theme_bw()
# Display side by side
grid.arrange(plot_novel, plot_conv, ncol = 2)

#for novel
cor.test(aumc_rss_data %>% filter(Type == "Novel") %>% pull(tm),
         aumc_rss_data %>% filter(Type == "Novel") %>% pull(AUMC), 
         method = "spearman",na.rm=TRUE) #rho=0.923899
sum(!is.na(aumc_rss_data %>% filter(Type == "Novel") %>% pull(tm)) & !is.na(aumc_rss_data %>% filter(Type == "Novel") %>% pull(AUMC)))#check sample size
#for conv
cor.test(aumc_rss_data %>% filter(Type == "Conventional") %>% pull(tm),
         aumc_rss_data %>% filter(Type == "Conventional") %>% pull(AUMC), 
         method = "spearman",na.rm=TRUE) #rho=0.9519375 
sum(!is.na(aumc_rss_data %>% filter(Type == "Conventional") %>% pull(tm)) & !is.na(aumc_rss_data %>% filter(Type == "Conventional") %>% pull(AUMC)))#check sample size

#comparing coefficient using z difference
zdifference <- function(r1, r2, n1, n2){
  zd <- (atanh(r1)-atanh(r2))/sqrt(1/(n1-3)+1/(n2-3))
  p <- 1 - pnorm(abs(zd))
  print(paste("Z Difference: ", zd))
  print(paste("One-Tailed P-Value: ", p))
  print(paste("Two-Tailed P-Value: ", (2*p)))
}

zdifference(0.923899  ,0.9519375  , 131, 49937) 




cor.test(aumc_rss_data$RSS_log,aumc_rss_data$tm, method = "spearman")
sum(!is.na(aumc_rss_data$AUMC) & !is.na(aumc_rss_data$tm))

library(gridExtra)
plot_novel <- ggplot(aumc_rss_data %>% filter(Type == "Novel"), aes(x = tm, y = RSS_log)) +
  geom_point(alpha = 0.3, color = "deepskyblue3") +
  geom_smooth(method = 'glm', se = TRUE, color = "black") +
  labs(title = "RSS_log v.s. tm-Novel",
       x = "tm",
       y = "RSS_log") +
  annotate("text", x = Inf, y = Inf, label = "rho = 0.277", hjust = 1.1, vjust = 1.2, size = 4) +
  theme_bw()

# Conventional plot with rho, axes reversed
plot_conv <- ggplot(aumc_rss_data %>% filter(Type == "Conventional"), aes(x = tm, y = RSS_log)) +
  geom_point(alpha = 0.3, color = "coral2") +
  geom_smooth(method = 'glm', se = TRUE, color = "black") +
  labs(title = "RSS_log v.s. tm- Conventional",
       x = "tm",
       y = "RSS_log") +
  annotate("text", x = Inf, y = Inf, label = "rho = 0.419", hjust = 1.1, vjust = 1.2, size = 4) +
  theme_bw()
# Display side by side
grid.arrange(plot_novel, plot_conv, ncol = 2)

#for novel
cor.test(aumc_rss_data %>% filter(Type == "Novel") %>% pull(tm),
         aumc_rss_data %>% filter(Type == "Novel") %>% pull(RSS_log), 
         method = "spearman",na.rm=TRUE) #rho=0.923899
sum(!is.na(aumc_rss_data %>% filter(Type == "Novel") %>% pull(tm)) & !is.na(aumc_rss_data %>% filter(Type == "Novel") %>% pull(RSS_log)))#check sample size
#for conv
cor.test(aumc_rss_data %>% filter(Type == "Conventional") %>% pull(tm),
         aumc_rss_data %>% filter(Type == "Conventional") %>% pull(RSS_log), 
         method = "spearman",na.rm=TRUE) #rho=0.9519375 
sum(!is.na(aumc_rss_data %>% filter(Type == "Conventional") %>% pull(tm)) & !is.na(aumc_rss_data %>% filter(Type == "Conventional") %>% pull(RSS_log)))#check sample size
zdifference(0.276795   ,0.4191812 , 131, 49937) 
