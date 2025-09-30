# Title: Function for calculating mean infiltration of large flow panels
# Description: Function to analyze large flow .csv files and generate a heatmap from proportion data
# Date: 2025.04.22
# Author: Bell Wu

library(pheatmap)
library(dplyr)

# input_df = rows are samples, columns represent variable. Contains sample name as row names
# rm_cols = columns to remove

ConvertFCSMeans <- function (input_df, sampleID , treatment, rm_cols = "base_ID", scale_fun = scale) {
  # if else statement to test whether or not df has row names  
  if (!all(row.names(input_df) == as.character(seq_len(nrow(input_df))))) {
    input_df$base_ID <- sub("_.*", "", row.names(input_df)) # add base_ID to the data frame if row names exist
  } else {
    input_df$base_ID <- sub("_.*", "", sampleID)
  }
  input_df <- input_df[ , !(is.na(names(input_df)) | names(input_df) == "")] # remove all columns with no names
  input_df$Treatment <- treatment
  vars <- setdiff(colnames(input_df), "base_ID") # sets the variables to pull from input df
  vars <- setdiff(vars, rm_cols)
  # so instead of aggregating based of the SampleID, I want to add another element to this function which is the treatment aspect of it
  input_df <- aggregate(input_df[vars],
                        by = list(base_ID = input_df$base_ID,
                                  treatment = input_df$Treatment),
                        FUN = mean) # calculate means
  # row.names(input_df) = input_df$base_ID
  input_col <- paste(input_df$base_ID, input_df$treatment, sep = "-")
  input_df <- input_df %>%
    select(-c("base_ID", "treatment"))
  input_df <- as.matrix(t(input_df))
  scaled_matrix <- t(apply(input_df, 1, scale_fun))
  colnames(scaled_matrix) = input_col
  return(scaled_matrix)
}

# RANDOM TESTING -------------------------------------------------------------

