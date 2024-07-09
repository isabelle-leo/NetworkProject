generate_qqplot <- function(data, group, feature = NULL) {
  plot_list <- list()
  features <- if (is.null(feature)) unique(data[[group]]) else feature
  for (f in features) {
    subset_data <- data[data[[group]] == f, ]
    p <- ggplot(subset_data, aes(sample = RSS_log)) +
      stat_qq(alpha = 0.5, color='blue') +
      stat_qq_line(color = "red") +
      ggtitle(paste("QQ Plot of Log-transformed RSS for", f))+
      theme_minimal()
    plot_list[[f]] <- p
  }
  plot_list
}

generate_density_plot_line <- function(data, group, value, novel_data = NULL, novel = TRUE) {
  data <- data[!is.na(data[[value]]) & !is.na(data[[group]]), ]
  # Convert the group variable to a factor
  data[[group]] <- as.factor(data[[group]])
  
  # Start building the plot with the main data
  p <- ggplot(data, aes_string(x = value, color = group, group = group)) +
    geom_density() +
    labs(title = paste("Density Plot of", value, "by", group),
         x = value, y = "Density", color = group) +
    theme_light()
  
  # Conditionally add the density plot for novel data if 'novel' is TRUE
  if (novel && !is.null(novel_data)) {
    p <- p + geom_density(data = novel_data, aes(x = value, fill = "Novel"), alpha = 0.5)
    p <- p + scale_fill_manual(values = c("Novel" = "blue"))
  }
  
  return(p)
}

generate_density_plot_fill <- function(data, group, value) {
  # Filter out NA values from the specified value and group columns
  data <- data[!is.na(data[[value]]) & !is.na(data[[group]]), ]
  # Convert the group variable to a factor
  data[[group]] <- as.factor(data[[group]])
  # Create the density plot
  ggplot(data, aes_string(x = value, fill = group)) +
    geom_density(alpha=0.5) +
    labs(title = paste("Density Plot of", value, "by", group),
         x = value,
         y = "Density",
         fill = group) +
    theme_light()
}

generate_boxplot <- function(data, group, value) {
  # Filter out NA values from the specified value and group columns
  data <- data[!is.na(data[[value]]) & !is.na(data[[group]]), ]
  # Convert the group variable to a factor
  data[[group]] <- as.factor(data[[group]])
  ggplot(data, aes_string(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.alpha = 0.3) +
    labs(title = paste("Boxplot of", value, "Values"),
         x = group,
         y = value) +
    theme_light()
}


