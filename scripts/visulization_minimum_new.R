# Load required libraries
library(dplyr)
library(ggplot2)
library(igraph)

###define functions
# Function to process graphs object and create hits_everything (data of intrest)
create_hits_everything <- function(graphs_comms) {
  id_vector <- c()
  gene_vector <- c()
  membership_vector <- c()
  
  for (i in names(graphs_comms)) {
    graph <- graphs_comms[[i]]
    unique_memberships <- unique(vertex_attr(graph, "membership"))
    id_vector <- c(id_vector, unique(paste0(graph$ioi, "_", unique_memberships)))
    gene_vector <- c(gene_vector, rep(graph$ioi, length(unique_memberships)))
    membership_vector <- c(membership_vector, unique_memberships)
  }
  
  hits_everything <- data.frame(
    id = id_vector, 
    hgnc_symbol = gene_vector,
    membership = membership_vector,
    pAdj = "N.D", 
    test = "meltome"
  )
  
  hits_everything
}


merge_func <- function(hits_everything, ALL_data) {
  merged_data <- left_join(hits_everything, ALL_data, by = "id")
  merged_data <- merged_data[, c("id", "hgnc_symbol", "membership")]
  return(merged_data)
}


gene_proteoform_count <- function(merged_data) {
  # Count the number of unique proteoforms for each gene
  gene_counts <- merged_data %>%
    group_by(hgnc_symbol) %>%
    summarise(proteoform_count = n_distinct(membership)) %>%
    arrange(desc(proteoform_count))

  return(gene_counts)
}

#function to create barplot
plot_func <- function(gene_counts) {
  # Count the number of genes for each proteoform count
  count_per_proteoform <- gene_counts %>%
    group_by(proteoform_count) %>%
    summarise(gene_count = n())
  
  ggplot(count_per_proteoform, aes(x = factor(proteoform_count), y = gene_count)) +
    geom_bar(stat = "identity", fill = 'steelblue') +  
    geom_text(aes(label = gene_count), vjust = -0.3, color = "black") +                     
    labs(x = "Number of Proteoforms", y = "Number of Genes", title = "Number of Genes my membership.no") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)
    )
}





###load data
#load ALL data
ALL_data <- readRDS("/Users/labrat/Downloads/ALL_new.RDS")

#load graph object
graphs_comms <- readRDS("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/PSM-proteoform pipeline/inputs/graphs_comms 1.RDS")



# process the graphs object to find data of interest
hits_everything <- create_hits_everything(graphs_comms)

# Merge data
merged_data <- merge_func(hits_everything, ALL_data)

# Count proteoforms
gene_counts <- gene_proteoform_count(merged_data)

# Print the gene proteoform counts
gene_counts_df<-as.data.frame(gene_counts)
head(gene_counts_df)
# Plot the data
plot_func(gene_counts)


