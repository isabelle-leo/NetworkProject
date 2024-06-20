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

# Function to merge data and create the bar plot
merge_and_plot <- function(hits_everything, specific_data, data_type) {
  merged_data <- left_join(hits_everything, specific_data, by = "id")
  
  gene_counts <- merged_data %>%
    group_by(membership) %>%
    tally(name = "n")
  
  ggplot(gene_counts, aes(x = factor(membership), y = n, fill = factor(membership))) +
    geom_bar(stat = "identity") +  
    geom_text(aes(label = n), vjust = -0.3, color = "black") +                     
    labs(x = "Membership", y = "Number of Gene Symbols", title = paste("ALL meltome_", data_type, sep = "")) +
    theme_bw()+
    theme(legend.position = "none")
}
#test
plot <- merge_and_plot(hits_everything, ALL_data, data_type)
print(plot)

###load data
#load ALL data
ALL_data <- readRDS("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/PSM-proteoform pipeline/outputs/all_from_CLL_default_normcenter.RDS")

#load graph object
graphs_comms <- readRDS("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/PSM-proteoform pipeline/inputs/graphs_comms 1.RDS")

# process the graphs object to find data of interest
hits_everything <- create_hits_everything(graphs_comms)

# Generate and display the plot
data_type <- "defalt"  # Define the data type, for example default
plot <- merge_and_plot(hits_everything, ALL_data, data_type)
print(plot)
#save plot
ggsave("plot.png", plot = plot, width = 10, height = 6, dpi = 300)

