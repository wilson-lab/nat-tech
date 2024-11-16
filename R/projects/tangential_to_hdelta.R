#################################################
### Hemibrain Connectivity Analysis and Visualization
### 
### This script analyzes and visualizes connectivity data from the Hemibrain dataset,
### focusing on tangential cell inputs to hDelta neurons in the Drosophila brain.
#################################################

### Load required libraries
# Data manipulation and visualization
library(tidyverse)
library(ggraph)
library(igraph)
library(tidygraph)
library(plotly)
# Database connection and manipulation
library(DBI)
library(RSQLite)
# Heatmap generation and color palettes
library(pheatmap)
library(RColorBrewer)

### Data source paths
# Paths to various SQLite databases containing connectivity data
# Uncomment and use the appropriate path based on your system and requirements
hemibrain.sql <- "/Volumes/neurobio/wilsonlab/banc/connectivity/hemibrain_v.1.2.1_data.sqlite"
# hemibrain.sql <- "/n/data1/hms/neurobio/wilson/banc/connectivity/hemibrain_v.1.2.1_data.sqlite" # On O2
# fafb.sql <- "/Volumes/neurobio/wilsonlab/banc/connectivity/flywire_783_data.sqlite"
# manc.sql <- "/Volumes/neurobio/wilsonlab/banc/connectivity/manc_1.2.1_data.sqlite"
# banc.sql <- "/Volumes/neurobio/wilsonlab/banc/connectivity/banc_data.sqlite"

### Output directory for generated images
hemibrain.images <- "inst/images/hemibrain/"

### Data Retrieval and Preprocessing
# Connect to the Hemibrain SQLite database
con <- dbConnect(RSQLite::SQLite(), hemibrain.sql)

# List and print all tables in the database
tables <- dbListTables(con)
print(tables)

# Helper function to get the mode (most frequent value) of a vector
get_mode <- function(x) {
  tbl <- table(x)
  names(tbl)[which.max(tbl)]
}

# Retrieve and process metadata
hb.meta <- dplyr::tbl(con, "meta") %>%
  # ... [rest of the code for processing metadata]
  
  # Lazy load the edgelist
  hb.elist <- dplyr::tbl(con, "edgelist")

### Analysis Setup
# Filter for tangential cells and hDelta neurons
hb.tangential <- hb.meta %>%
  dplyr::filter(grepl("^FB",cell_type))
hb.hdelta <- hb.meta %>%
  dplyr::filter(grepl("hDelta",cell_type))

# Collect IDs for tangential cells and hDelta neurons
tangential_ids <- hb.tangential$bodyid
hdelta_ids <- hb.hdelta$bodyid

#############################################
### collect data from large SQL data base ###
#############################################

### NOTE: if this is slow, it could be because you are reading from remote, copy and store locally for faster performance ###

# Now, use these IDs in the query
# Faster if counr filter is stricter
# hb.elist.th <- hb.elist %>%
#   dplyr::filter(count >= 20) %>%
#   dplyr::filter(pre %in% !!tangential_ids) %>%
#   dplyr::filter(post %in% !!hdelta_ids) %>%
#   dplyr::collect()

# Construct the SQL query
query <- paste0("
  SELECT *
  FROM edgelist
  WHERE count >= 5
    AND pre IN ('", paste(tangential_ids, collapse = "','"), "')
    AND post IN ('", paste(hdelta_ids, collapse = "','"), "')
")

# Execute the query and collect the results
hb.elist.th <- dbGetQuery(con, query)

# First, let's merge the edgelist with meta data
hb.elist.merged <- hb.elist.th %>%
  dplyr::select(-pre_top_nt, -post_top_nt) %>%
  dplyr::left_join(hb.meta %>% select(bodyid, cell_type, top_nt), by = c("pre" = "bodyid")) %>%
  dplyr::rename(pre_cell_type = cell_type, pre_top_nt = top_nt) %>%
  dplyr::left_join(hb.meta %>% select(bodyid, cell_type, top_nt), by = c("post" = "bodyid")) %>%
  dplyr::rename(post_cell_type = cell_type, post_top_nt = top_nt) %>%
  dplyr::filter(post_label %in% c("axon","dendrite"),
                pre_label %in% c("axon","dendrite"))

# Reshape and collapse the data by cell_type
hb.elist.collapsed <- hb.elist.merged %>%
  dplyr::mutate(connection_type = paste(pre_label, "to", post_label)) %>%
  dplyr::group_by(pre_cell_type, post_cell_type, connection_type) %>%
  dplyr::summarise(count = sum(count), 
                   pre_top_nt = first(pre_top_nt),
                   post_top_nt = first(post_top_nt),
                   .groups = 'drop') %>%
  dplyr::arrange(dplyr::desc(count))

############
### plot ###
############

# Step 1: Aggregate the data
aggregated_data <- hb.elist.collapsed %>%
  filter(connection_type %in% c("axon to dendrite", "axon to axon")) %>%
  group_by(pre_cell_type, post_cell_type, connection_type) %>%
  summarise(total_count = sum(count), .groups = "drop")

# Step 2: Reshape the aggregated data
heatmap_data <- aggregated_data %>%
  mutate(connection_type = factor(connection_type, levels = c("axon to dendrite", "axon to axon"))) %>%
  unite("post_cell_connection", post_cell_type, connection_type, sep = "_") %>%
  pivot_wider(names_from = post_cell_connection, values_from = total_count, values_fill = 0) %>%
  arrange(pre_cell_type) %>%  # Sort rows alphabetically
  column_to_rownames("pre_cell_type")

# Step 3: Create column annotation for axon vs dendrite
column_annotation <- data.frame(
  connection_type = ifelse(grepl("axon to axon$", colnames(heatmap_data)), "axon", "dendrite"),
  row.names = colnames(heatmap_data)
)

# Step 4: Create row annotation for pre_top_nt
row_annotation <- hb.elist.collapsed %>%
  distinct(pre_cell_type, pre_top_nt) %>%
  column_to_rownames("pre_cell_type")

# Ensure row_annotation matches heatmap_data rows
row_annotation <- row_annotation[rownames(heatmap_data), , drop = FALSE]

# Step 5: Create column order
column_order <- names(heatmap_data) %>%
  str_extract("^[^_]+") %>%
  unique() %>%
  sort() %>%
  sapply(function(x) {
    connections <- c(paste(x, "axon to dendrite", sep = "_"),
                     paste(x, "axon to axon", sep = "_"))
    connections[connections %in% names(heatmap_data)]
  }) %>%
  unlist()

# Ensure column_annotation only includes columns present in heatmap_data
column_annotation <- column_annotation[column_order, , drop = FALSE]

# Step 6: Define colors for axon vs dendrite and pre_top_nt
pre_top_nt_colors <- c(
  acetylcholine = "#EF7C12",
  glutamate = "#8FDA04",
  gaba = "#1BB6AF",
  serotonin = "#FBBB48",
  dopamine = "#F4E3C7",
  octopamine = "#C70E7B",
  unknown = "grey"
)
ann_colors <- list(
  connection_type = c(axon = "#EE4244", dendrite = "#007BC3"),
  pre_top_nt = pre_top_nt_colors
)

# Step 7: Create the heatmap
pheatmap(heatmap_data[, column_order],
         cluster_rows = FALSE,  # Disable row clustering
         cluster_cols = FALSE,
         annotation_col = column_annotation,
         annotation_row = row_annotation,
         annotation_colors = ann_colors,
         show_colnames = TRUE,
         show_rownames = TRUE,
         labels_col = str_extract(column_order, "^[^_]+"),
         fontsize_row = 6,  # Adjust as needed
         fontsize_col = 8,  # Adjust as needed
         angle_col = 90,
         annotation_names_col = FALSE,
         annotation_names_row = FALSE,
         annotation_legend = TRUE)

# Step 8: Save
ggsave(filename=file.path(hemibrain.images,"hemibrain_tangential_neuron_to_hdelta_plot.png"), 
       heatmap_plot, width = 12, height = 10, units = "in")


