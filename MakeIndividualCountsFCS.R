# Title: Function for creating a counts matrix to plot individual HM from .fcs files
# Description: Function to analyze large flow .csv files and generate a heatmap from proportion data
# Date: 2025.04.22
# Author: Bell Wu

# input_df should be a counts matrix with samples as row names
# rm_rows = data rows that you want to remove from the final output

MakeIndividualCountsFCS <- function(input_df, sampleID, rm_rows = "", scale_fun = scale) { 
  input_df <- input_df[ , !(is.na(names(input_df)) | names(input_df) == "")]
  input_df <- input_df %>%
    select(-contains(rm_rows))
  input_df <- t(input_df)
  df_col <- colnames(input_df)
  counts_matrix <- t(apply(input_df, 1, scale_fun))
  colnames(counts_matrix) = sampleID
  return(counts_matrix)
}

# RANDOM TESTING -------------------------------------------------------------

# MakeIndividualCountsFCS <- function(input_df, rm_rows = "", scale_fun = scale) { 
#   input_df <- input_df[ , !(is.na(names(input_df)) | names(input_df) == "")]
#   vars <- setdiff(colnames(CD8_activation), vars)
#   test <- CD8_activation %>%
#     select(-contains(vars))
#   input_df <- t(input_df)
#   df_col <- colnames(input_df)
#   counts_matrix <- t(apply(input_df, 1, scale_fun))
#   colnames(counts_matrix) = df_col
#   return(counts_matrix)
# }
# vars
# head(CD8_activation)
# head(test)
