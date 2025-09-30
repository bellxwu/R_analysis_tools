# Title: GetAbsoluteCell
# Description: Function to get the absolute cell count from a .csv file from flowjo
# Date: 2025.06.27
# Author: Bell Wu

GetAbsoluteCell <- function(input_df, Treatment, SampleID, bead_vol, cell_vol, bead_conc) {
# remove characters  
  cell_counts <- input_df |> 
    select(!where(is.character))
# function for calculating abs cell counts  
  abs_count <- function(cell_counts, bead_vol, ebead, cell_vol, bead_conc) {
    (cell_counts * bead_vol) / (ebead * cell_vol) * bead_conc * (bead_vol + cell_vol)
  }
# calculate the absolute count
  if (any(colnames(input_df) == "eBeads")) {
    abs_count <- abs_count(cell_counts = cell_counts,
                         bead_vol = bead_vol,
                         ebead = cell_counts$eBeads,
                         cell_vol = cell_vol,
                         bead_conc = bead_conc)
  } else {
      print("No column named eBeads")
  }
  output_df <- abs_count
  output_df$Treatment <- Treatment
  output_df$SampleID <- SampleID
  output_df <- as_tibble(output_df)
  return(output_df)
}
# Testing ---------
