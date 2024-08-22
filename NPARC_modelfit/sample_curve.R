#plotting functions
plot_peptide_n <- function(peptide_name) {
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
plot_peptide_c <- function(peptide_name) {
  # Filter the data for the selected peptide
  merged_data_sample <- df_conv_subset[df_conv_subset$peptide_id == peptide_name, ]
  
  # Ensure sample_name is a factor
  merged_data_sample$sample_name <- as.factor(merged_data_sample$sample_name)
  
  # Create the main plot
  sample_plot_orig_2 <- ggplot(merged_data_sample, aes(x = temperature, y = relAbundance)) +
    geom_point(aes(color = sample_name)) +
    theme_bw() +
    ggtitle(paste(peptide_name, "(Type=Conventional)")) +
    scale_color_discrete(name = "Sample Name")
  
  # Filter predictions and metrics for the selected peptide
  stk4Predictions <- filter(fits_conv$predictions, id == peptide_name)
  stk4Metrics <- filter(fits_conv$metrics, id == peptide_name)
  
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
plot_peptide_n("peptide146384")
plot_peptide_n("peptide144295")
plot_peptide_n("peptide145732")
plot_peptide_c("peptide144448")
plot_peptide_c("peptide145869")
