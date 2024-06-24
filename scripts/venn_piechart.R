# Load required libraries
library(dplyr)
library(igraph)
library(ggplot2)
library(ggvenn)

### Load data
# Load ALL data
ALL_output <- readRDS("/Users/labrat/Downloads/ALL_new.RDS")

# Load graph object
graphs_comm <- readRDS("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/PSM-proteoform pipeline/inputs/graphs_comms 1.RDS")

# Load psm table
psm_table <- read.delim("/Users/labrat/Library/CloudStorage/OneDrive-共享的库-Onedrive/Sweden/MTLS/summer intern/💻Research project/PSM-proteoform pipeline/inputs/target_psmtable.txt")

# Function to remove mass modifications from peptide sequences
normalize_peptide <- function(peptide) {
  gsub("\\+\\d+\\.\\d+", "", peptide)
}

# Function to process graphs object and create hits_everything (data of interest)
create_hits_everything <- function(graphs_comms) {
  data_list <- lapply(names(graphs_comms), function(i) {
    graph <- graphs_comms[[i]]
    data.frame(
      name = vertex_attr(graph, "name"),
      gene_symbol = vertex_attr(graph, "id"),
      protein_ids = vertex_attr(graph, "protein_ids"),
      first_protein_id = vertex_attr(graph, "first_protein_id"),
      membership = vertex_attr(graph, "membership"),
      Peptide = vertex_attr(graph, "name"),  # Assuming 'name' contains the peptide sequence
      NormalizedPeptide = sapply(vertex_attr(graph, "name"), normalize_peptide),
      pAdj = "N.D",
      test = "meltome"
    )
  })
  
  hits_everything <- do.call(rbind, data_list)
  return(hits_everything)
}

# Create hits_everything dataframe from graphs_comms
hits_everything_graph <- create_hits_everything(graphs_comm)

# Print the first few rows of hits_everything
print(head(hits_everything_graph, n = 3))

# Normalize peptide sequences in psm_table and add as a new column
psm_table <- psm_table %>%
  mutate(NormalizedPeptide = sapply(Peptide, normalize_peptide))

# Normalize the 'name' column in hits_everything_graph and store in a new column 'NormalizedPeptide'
hits_everything_graph$NormalizedPeptide <- sapply(hits_everything_graph$name, normalize_peptide)

# Print the first few rows of hits_everything_graph to verify
print(head(hits_everything_graph, n = 2))


###analysis and plotting
# Extract normalized peptides from hits_everything_graph and psm_table
peptides_graph <- unique(hits_everything_graph$NormalizedPeptide)
peptides_psm <- unique(psm_table$NormalizedPeptide)

# Find overlapping peptides
overlapping_peptides <- intersect(peptides_graph, peptides_psm)

# check the summary of peptides
cat( 
  "Dimensions of hits_everything_graph: ", dim(hits_everything_graph), "\n",
  "Dimensions of psm_table: ", dim(psm_table), "\n",
  "Unique Peptides number in hits_everything_graph: ", length(unique(hits_everything_graph$NormalizedPeptide)), "\n",
  "Unique Peptides number in psm_table: ", length(unique(psm_table$NormalizedPeptide)), "\n",
  "Number of overlapping peptides: ", length(intersect(unique(hits_everything_graph$NormalizedPeptide), unique(psm_table$NormalizedPeptide))), "\n"
)
print(head(overlapping_peptides))

# Prepare data for Venn diagram
venn_data <- list(
  Graphs_com = peptides_graph,
  PSM_Table = peptides_psm
)

# Create Venn diagram
ggvenn(venn_data,  stroke_size = 0.5, set_name_size = 4)+
  ggtitle("Peptides in ALL meltome vs. proteome ")

#pie chart
# Calculate counts for unique and overlapping peptides
counts <- c(
  length(setdiff(psm_table$NormalizedPeptide, hits_everything_graph$NormalizedPeptide)), 
  length(setdiff(hits_everything_graph$NormalizedPeptide, psm_table$NormalizedPeptide)), 
  length(intersect(psm_table$NormalizedPeptide, hits_everything_graph$NormalizedPeptide))
)
labels <- c("PSM Table Only", "Graphs Only", "Overlap")

# Create labels with counts and percentages
total_peptides <- sum(counts)
percentage <- round(counts / total_peptides * 100, 1)
pie_labels <- paste0(labels, "\n", counts, " (", percentage, "%)")

# Plot pie chart
pie(counts, labels = pie_labels, main = "Peptides found in PSM Table and Graphs", col = c("#FF9999", "#99CCFF", "#CCCCFF"))
#check if data is correct
# Find the number of unique values in the name column
 length(unique(hits_everything$name))
 length(unique(psm_table$Peptide))

##note: the graph only contains the unique, normalised peptide seq.which means all duplicated peptide seq after normalizing is removed.


###looking into proteoforms

# Create a new id column using gene_symbol and membership
hits_everything_graph <- hits_everything_graph %>%
  mutate(id = paste(gene_symbol, membership, sep = "_"))

# Display the head of the modified dataframe
head(hits_everything_graph)


# Calculate counts for unique and overlapping proteoforms
counts <- c(
  length(setdiff(ALL_data$id, hits_everything_graph$id)), 
  length(setdiff(hits_everything_graph$id, ALL_data$id)), 
  length(intersect(ALL_data$id, hits_everything_graph$id))
)
labels <- c("ALL Data Only", "Hits Only", "Overlap")

# Create labels with counts and percentages
total_proteoforms <- sum(counts)
percentage <- round(counts / total_proteoforms * 100, 1)
pie_labels <- paste0(labels, "\n", counts, " (", percentage, "%)")

# Print data to check if it is correct
print(data.frame(Category = labels, Count = counts, Percentage = percentage))

# Plot pie chart
pie(counts, labels = pie_labels, main = "Proteoforms found in ALL Data and Graphs_comm", col = c("#FF9999", "#99CCFF", "#CCCCFF"))
#check data
length(unique(ALL_output$id))#number of different proteomes found in two datasets
length(unique(hits_everything_graph$id))




#extra: look into gene symbol

# Extract gene symbols from both datasets
genes_graph <- unique(hits_everything_graph$gene_symbol)
genes_psm <- unique(psm_table$Gene.Symbol)

# Find overlapping gene symbols
overlapping_genes <- intersect(genes_graph, genes_psm)


# Prepare data for Venn diagram
venn_data <- list(
  Graphs_comm = genes_graph,
  PSM_Table = genes_psm
)

# Create Venn diagram
ggvenn(venn_data, fill_color = c("#FF9999", "#99CCFF"), stroke_size = 0.5, set_name_size = 4)+
  ggtitle("Genes in ALL meltome vs. proteome ")

#create pie chart 
## Calculate counts for unique and overlapping gene symbols
counts <- c(
  length(setdiff(genes_graph, genes_psm)), 
  length(setdiff(genes_psm, genes_graph)), 
  length(intersect(genes_graph, genes_psm))
)
labels <- c("Graphs Only", "PSM Table Only", "Overlap")

# Create labels with counts and percentages
total_genes <- sum(counts)
percentage <- round(counts / total_genes * 100, 1)
pie_labels <- paste0(labels, "\n", counts, " (", percentage, "%)")

# Plot pie chart
pie(counts, labels = pie_labels, main = "Genes found in ALL proteome and meltome", col = c("#FF9999", "#99CCFF", "#CCCCFF"))

