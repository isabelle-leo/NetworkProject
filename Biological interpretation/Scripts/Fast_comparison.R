generate_p_value_table <- function(novel_data, conventional_alter, conventional_null, group_interest, p_thre = 0.05) {
  # Function to format p-value and add asterisk if below threshold
  format_p_value <- function(p_value, threshold) {
    p_value <- as.numeric(p_value)
    formatted_p_value <- format(p_value, digits = 4)
    if (p_value < threshold) {
      formatted_p_value <- paste0(formatted_p_value, "*")
    }
    return(formatted_p_value)
  }
  
  # Perform Wilcox tests
  wilcox_novel_alter <- format_p_value(wilcox.test(novel_data$RSS_log, conventional_alter$RSS_log, alternative = "two.sided")$p.value, p_thre)
  wilcox_null_alter <- format_p_value(wilcox.test(conventional_null$RSS_log, conventional_alter$RSS_log)$p.value, p_thre)
  wilcox_novel_null <- format_p_value(wilcox.test(novel_data$RSS_log, conventional_null$RSS_log)$p.value, p_thre)
  
  # Perform KS tests
  ks_novel_alter <- format_p_value(ks.test(novel_data$RSS_log, conventional_alter$RSS_log)$p.value, p_thre)
  ks_null_alter <- format_p_value(ks.test(conventional_alter$RSS_log, conventional_null$RSS_log)$p.value, p_thre)
  ks_novel_null <- format_p_value(ks.test(novel_data$RSS_log, conventional_null$RSS_log)$p.value, p_thre)
  
  # Create data frame for Wilcox test results
  wilcox_results <- data.frame(
    Novel = c(wilcox_novel_null, wilcox_novel_alter),
    Alter = c(wilcox_null_alter, "-")
  )
  rownames(wilcox_results) <- c("Null", "Alter")
  
  # Create data frame for KS test results
  ks_results <- data.frame(
    Novel = c(ks_novel_null, ks_novel_alter),
    Alter = c(ks_null_alter, "-")
  )
  rownames(ks_results) <- c("Null", "Alter")
  
  # Get unique values of group_interest
  unique_alter <- unique(conventional_alter[[group_interest]])
  unique_null <- unique(conventional_null[[group_interest]])
  
  # Combine results into a list
  results <- list(
    "Wilcox test" = wilcox_results,
    "Kolmogorov-Smirnov test" = ks_results,
    "Unique_Alter_Group" = paste(unique_alter, collapse = ", "),
    "Unique_Null_Group" = paste(unique_null, collapse = ", ")
  )
  
  return(results)
}


## Assign groups for stat comparision

#novel group
novel_data <- filter(rss_data_novel, !is.na(RSS))

#*choose groups of interest & ioi for the 'alter' group
colnames(rss_data_conventional_info) #you can check here
unique(rss_data_conventional_info[[group_interest]])

####assign groups&ioi####!!

group_interest<-"disorder_quantile" 
alter_locations <-c('Q4') #shift-cmd-enter (run all)

#assign alter and null group
conventional_alter <- rss_data_conventional_info %>%
  filter(get(group_interest) %in% alter_locations, !is.na(RSS_log), !is.na(get(group_interest)))
conventional_null <- rss_data_conventional_info %>%
  filter(!get(group_interest) %in% alter_locations, !is.na(RSS_log), !is.na(get(group_interest)))

# generate 3-in-1 density plot
ggplot() +
  geom_density(data = conventional_null, aes(x = RSS_log, fill = "null"), alpha = 0.5) +
  geom_density(data = conventional_alter, aes(x = RSS_log, fill = "alter"), alpha = 0.5) +
  geom_density(data = novel_data, aes(x = RSS_log, fill = "Novel"), alpha = 0.5) +
  labs(title = paste("RSS_log Distribution by", group_interest),
       x = "RSS_log", y = "Density") +
  theme_light() +
  theme(legend.position = "top")

results <- generate_p_value_table(novel_data, conventional_alter, conventional_null, group_interest,p_thre=0.01)



##OPTIONAL: curtosis value##

# Calculate kurtosis for each type
#kurtosis_novel <- kurtosis(novel_data$RSS_log, na.rm = TRUE)
#kurtosis_alter <- kurtosis(conventional_alter$RSS_log, na.rm = TRUE)
#kurtosis_null <- kurtosis(conventional_null$RSS_log, na.rm = TRUE)

#print results
kruskal.test(as.formula(paste("RSS_log ~", group_interest)), data = rss_data_conventional_info)
print(results)

#OPTIONAL: Display the kurtosis values
#cat("Kurtosis for Novel RSS: ", kurtosis_novel, "\n")
#cat("Kurtosis for Conventional Alter RSS: ", kurtosis_alter, "\n")
#cat("Kurtosis for Conventional Null RSS: ", kurtosis_null, "\n")

