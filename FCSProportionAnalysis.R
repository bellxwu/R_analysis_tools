# Title: Proportion analyses for large flow panels
# Description: Analyze large flow .csv files and generate a proportion data frame
# Date: 2025.04.23
# Author: Bell Wu

# input_df = .csv file from flowjo
# rm_vars = variables you want to remove
# population = cell population count you want to divide all variables in data frame by. Input with ""

GetProportions <- function(input_df, sampleID,  rm_vars, population) {
  sample_ID <- sampleID
  vars <- colnames(input_df)
  denom <- input_df[[population]]
  final_vars <- setdiff(vars, rm_vars)
  select_df <- input_df %>%
    select(any_of(final_vars))
  output_df <- select_df[ , FALSE]
  for (i in seq_along(select_df)) {
    varname <- final_vars[i]
    output_df[[varname]] <- (select_df[[varname]] / denom) * 100
  }
  output_df$SampleID <- sample_ID
  return(output_df)
}
