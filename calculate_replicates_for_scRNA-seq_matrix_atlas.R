

#Piece of code in R to estimate replicates from single cell original dataframes: 


#Stimating replicates:
# Get the column names
col_names <- colnames(pn1.int)

# Identify unique patterns (e.g., B_cells, Basophils, Endo) in column names
# Assuming patterns are unique names without trailing numbers
patterns <- unique(gsub(".\\d+", "", col_names))  # removes trailing numbers after point if present

# Initialize an empty list to store the resulting mean columns
mean_columns <- list()

for (pattern in patterns) {
  # Find columns matching the current pattern
  matching_cols <- grep(paste0("^", pattern), col_names, value = TRUE)
  
  # Check if there are columns that match the pattern
  if (length(matching_cols) > 0) {
    # Calculate group size to obtain exactly 10 mean columns for each pattern
    group_size <- ceiling(length(matching_cols) / 10)
    
    # Split the columns into 10 groups (or as close as possible) by calculated group size
    groupings <- split(matching_cols, ceiling(seq_along(matching_cols) / group_size))
    
    # Calculate mean for each group and add to mean_columns list
    for (i in seq_len(min(10, length(groupings)))) {  # limit to 5 groups
      group <- groupings[[i]]
      mean_columns[[paste0(pattern, "_mean_", i)]] <- rowMeans(pn1.int[group], na.rm = TRUE)
    }
  }
}

# Combine mean columns into a new data frame
pn1_mean_df <- as.data.frame(mean_columns)

pn1_mean_names <- colnames(pn1_mean_df)
pn1_mean_names <- gsub("_mean_[0-9]*$","",pn1_mean_names)
colnames(pn1_mean_df) <- pn1_mean_names

#Reorder columns:
pn1_mean_ord_df <-  pn1_mean_df[, order(names(pn1_mean_df))]