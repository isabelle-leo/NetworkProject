library(dplyr)
library(magrittr)
library(ggplot2)
library(tibble)
rss_data_novel<-read.delim('/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/inputs/rss_data_novel_detect.txt')
rss_data<-read.delim("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 3B-bio interpretation/inputs/rss_data_full.txt")
table(rss_data$Type)
table(rss_data_novel$HusenDetect)


# Convert rownames to column in both datasets
rss_data$pep_id <- rownames(rss_data)
rss_data_novel$pep_id <- rownames(rss_data_novel)

rss_data_new<-merge(x=rss_data,y=rss_data_novel,by='pep_id',all=TRUE)

#drop useless columns
rss_data_new <- rss_data_new[, !names(rss_data_new) %in% c("Type.y", "RSS.y", "Gene.ID.y", "Peptide.sequence.y", 
                                                           "cell_line_count.y", "protein_ids.y", "RSS_log", "pep_id")]

# View the column names to confirm removal
colnames(rss_data_new)
colnames(rss_data_new) <- gsub("\\.x$", "", colnames(rss_data_new))
# View the column names to confirm the changes
colnames(rss_data_new)


rss_data_new$HusenDetect <- ifelse(is.na(rss_data_new$HusenDetect), 'conventional', 
                                   ifelse(rss_data_new$HusenDetect, 'detected', 'undetected'))
rss_data_new$RSS_log<-log(rss_data_new$RSS)
# View the updated HusenDetect column to confirm the changes
head(rss_data_new$HusenDetect)
table(rss_data_new$HusenDetect)
head(rss_data_new)

generate_density_plot_fill(rss_data_new,'Type','RSS_log')
generate_density_plot_fill(rss_data_new,'HusenDetect','RSS_log')
generate_density_plot_line(rss_data_new,'HusenDetect','RSS_log')

# Generate boxplot for cell_line_count for the three groups in HusenDetect
ggplot(rss_data_new, aes(fill = HusenDetect, y = cell_line_count)) +
  geom_boxplot() +
  labs(title = "Cell Line Count by HusenDetect Group",
       x = "HusenDetect Group",
       y = "Cell Line Count") +
  theme_light()

generate_density_plot_line(rss_data_new,'HusenDetect','cell_line_count')

# Generate histograms for cell_line_count for each group in HusenDetect
ggplot(rss_data_new, aes(x = cell_line_count,fill = HusenDetect)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~ HusenDetect, scales = "free_y") +
  labs(title = "Histogram of Cell Line Count by HusenDetect Group",
       x = "Cell Line Count",
       y = "Frequency") +
  theme_minimal()

write.table(rss_data_new,'/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/Step 4-redo analysis/inputs/rss_data_4.txt',sep = "\t", row.names = TRUE)
