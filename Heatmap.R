# Title: Heatmap analyses for large flow panels
# Description: Analyze large flow .csv files and generate a heatmap from proportion data
# Date: 2025.04.17
# Author: Bell Wu

# set working directory and read .csv
setwd("~/desktop/Lok Lab/*HuMice Project - Model generation/*Cytokine panel, cell-line analyses/2025.04.02 Spleen Hu-Cell-line/")
CD4_cytokines <- read.csv("2025.04.02 Hu-Cell-line CD4 Spleen analyses (CD4 cytokine).csv")
head(CD4_cytokines)

# 1.1 Setting up the data into a counts matrix --------------------------------

# transpose data such that samples are columns
CD4_cytokines <- t(CD4_cytokines)
# have only the counts matrix with sample name embedded into column name
sample_ID <- CD4_cytokines["Sample.ID", ]
colnames(CD4_cytokines) = sample_ID 
CD4_counts <- as.matrix(CD4_cytokines[-c(1:2), ]) # remove non-numeric row

str(CD4_counts)
dim(CD4_counts)
# set matrix to numeric 
class(CD4_counts) <- "numeric" 

# 1.2 Plotting as heatmap of individual values --------------------------------

# need to scale data accordingly
cal_z_score <- function(x) {
  (x-mean(x)) / sd(x)
}

# calculate the z-score (ie. standardize the counts)
CD4_z <- t(apply(CD4_counts, 1, cal_z_score))
CD4_s <-t(apply(CD4_counts, 1, scale)) # scale and calculating z-score achieves the same
colnames(CD4_s) = sample_ID 

# remove CD3/live
CD4_z <- CD4_z[-1, ]

# set colour for heatmap
library(pheatmap)
library(RColorBrewer)

# set colour for heatmap.
h_colour <- colorRampPalette(brewer.pal(8, "RdBu"))(100)
# # to centre around 0 if needed
# range_val <- max(abs(CD4_z))
# breaks <- seq(-range_val, range_val, length.out = 101)

# build annotation column 
sample_subtype <- data.frame()
sample_subtype <- data.frame (samples = c("SCLC-N", "SCLC-N", "SCLC-N", # H1694
                                          "SCLC-N", "SCLC-N", "SCLC-N", # H446
                                          "SCLC-A", "SCLC-A", "SCLC-A", "SCLC-A", "SCLC-A", # H69
                                          "SCLC-I", "SCLC-I","SCLC-I","SCLC-I","SCLC-I", # H841
                                          "SCLC-I","SCLC-I","SCLC-I", # SBC5
                                          "SCLC-A", "SCLC-A", "SCLC-A", "SCLC-A", "SCLC-A")) # SHP77
row.names(sample_subtype) = sample_ID


library(pheatmap)
GMA <- pheatmap(CD4_z, scale = "none", fontsize_row = 10, cellheight = 10, cellwidth = 20,
                width = 400, height = 1000, cluster_rows = FALSE, color = rev(h_colour),
                annotation_col = sample_subtype,
                border_color = "BLACK")

pdf(file = "Hu-Cell-line pheatmap [CD4 activation].pdf", width = 15, height = 5)
grid::grid.newpage()
grid::grid.draw(GMA$gtable)
dev.off()

# 1.3 Plotting averages instead of individual samples -------------------------

library(dplyr)
library(stringr)
# want to have samples for rows. 
head(CD4_cytokines)

# remove .fcs file names to create a base ID.
base_ID <- sub("_.*", "", CD4_cytokines$Sample.ID)
CD4_cytokines$base_ID <- base_ID

# add base IDs to the data frame
CD4_cytokines$base_ID <- base_ID
# Select the necessary variables
CD4_vars <- setdiff(colnames(CD4_cytokines), c("Sample.ID", "Treatment", "base_ID"))
# use the by interface of the aggregate function to create new data frame of means
means_CD4 <- aggregate(
  CD4_cytokines[CD4_vars],
  by = list(base_ID = CD4_cytokines$base_ID),
  FUN = mean
)
# create the counts matrix with the samples as columns
row.names(means_CD4) = means_CD4$base_ID
CD4_mcounts <- t(means_CD4)
CD4_mcounts <- CD4_mcounts[row.names(CD4_mcounts) != "base_ID", ]
class(CD4_mcounts) <- "numeric"
# normalize with scale function
sample_ID <- colnames(CD4_mcounts)
CD4_s <-t(apply(CD4_mcounts, 1, scale))
colnames(CD4_s) = sample_ID

# 1.4 Sourcing R script from previous function (CD4) ----------------------------
library(dplyr)
source("~/R_programming/R_analyses_scripts/ConvertFCSMeans.R")
setwd("~/desktop/Lok Lab/*HuMice Project - Model generation/*Cytokine panel, cell-line analyses/2025.04.02 Spleen Hu-Cell-line/")
CD4_cytokines <- read.csv("2025.04.02 Hu-Cell-line CD4 Spleen analyses (CD4 cytokine).csv")
head(CD4_cytokines)
CD4_means <- ConvertFCSMeans(CD4_cytokines, rm_cols = c("Sample.ID", "Treatment", "CD3.Live", "CD4.CD3",
                                                        "CD25.CD4"))
head(CD4_means)

# build heatmap
library(pheatmap)
library(RColorBrewer)
# set colour for heatmap.

h_colour <- colorRampPalette(brewer.pal(8, "RdBu"))(100)
GMA <- pheatmap(CD4_means, scale = "none", fontsize_row = 10, cellheight = 10, cellwidth = 20,
                width = 400, height = 1000, cluster_rows = FALSE, color = rev(h_colour),
                border_color = "BLACK"
                )

pdf(file = "Hu-Cell-line pheatmap means [CD4 activation].pdf", width = 15, height = 5)
grid::grid.newpage()
grid::grid.draw(GMA$gtable)
dev.off()

# 1.5 HM of individual samples for CD8 ----------------------------------------
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# preparing the counts matrix
setwd("~/Desktop/Lok Lab/*HuMice Project - Model generation/*HuMice Baseline TINF analyses/2025.04.08 Hu-Cell-line Spleen [aPDL1 activation]/")
CD8_activation <- read.csv("2025.04.24 Hu-cell-line spleen aPDL1 activation, unbiased.csv")
row.names(CD8_activation) = CD8_activation$Sample.ID..1M.H1694.
head(CD8_activation)

source("~/R_programming/R_analyses_scripts/MakeIndividualCountsFCS.R")

vars <- c("Live", "CD4.CD3", "CD8.CD3", "Sample", "Treat")
CD8_z <- MakeIndividualCountsFCS(CD8_activation, rm_rows = vars, scale_fun = cal_z_score)
head(CD8_z)

# preparing annotation column
sample_subtype <- data.frame (Subtype = c("SCLC-N", "SCLC-N", "SCLC-N", # H1694
                                          "SCLC-N", "SCLC-N", "SCLC-N", # H446
                                          "SCLC-A", "SCLC-A", "SCLC-A", "SCLC-A", "SCLC-A", # H69
                                          "SCLC-I", "SCLC-I","SCLC-I","SCLC-I","SCLC-I", # H841
                                          "SCLC-I","SCLC-I","SCLC-I", # SBC5
                                          "SCLC-A", "SCLC-A", "SCLC-A", "SCLC-A", "SCLC-A")) # SHP77
row.names(sample_subtype) = CD8_colnames

# building heatmap
h_colour <- colorRampPalette(c("#E31A1C", "#FFFFFF" ,"#1F78B4"))(100)
range_val <- max(abs(CD8_z))
breaks <- seq(-range_val, range_val, length.out = 101)

HM_CD8 <- pheatmap(CD8_z, scale = "none", fontsize_row = 10, cellheight = 10, cellwidth = 20,
                   width = 400, height = 1000, color = rev(h_colour), breaks = breaks,
                   annotation_col = sample_subtype,
                   border_color = "BLACK")

pdf(file = "Hu-Cell-line pheatmap individual [PDL1 activation] _test.pdf", width = 15, height = 5)
grid::grid.newpage()
grid::grid.draw(HM_CD8$gtable)
dev.off()


# 1.5 HM of mean of samples for CD8 -------------------------------------------

setwd("~/Desktop/Lok Lab/*HuMice Project - Model generation/*HuMice Baseline TINF analyses/2025.04.08 Hu-Cell-line Spleen [aPDL1 activation]/")
CD8_activation <- read.csv("2025.04.24 Hu-cell-line spleen aPDL1 activation, unbiased.csv")
head(CD8_activation)
CD8_cols <- colnames(CD8_activation)

CD8_means <- ConvertFCSMeans(CD8_activation, rm_cols = c("Sample.ID..1M.H1694.", "Treatment",
                                                        "CD3.Live", "CD8.CD3", "CD4.CD3"
                                                        ))
head(CD8_means)
annot_col <- data.frame (Subtype = c("NE-high", "NE-high", "NE-high", "NE-low", "NE-low", "NE-high"))
row.names(annot_col) <- colnames(CD8_means)
ann_colour <- list(
  Subtype = c("NE-high" = "#984ea3", "NE-low" = "#ff7f00")) # "SCLC-A" = "#4daf4a",

h_colour <- colorRampPalette(brewer.pal(8, "RdBu"))(100)
GMA <- pheatmap(CD8_means, scale = "none", fontsize_row = 10, cellheight = 10, cellwidth = 20,
                width = 400, height = 1000, cluster_rows = TRUE, color = rev(h_colour),
                border_color = "BLACK", annotation_col = annot_col, annotation_colors = ann_colour
)

pdf(file = "Hu-Cell-line pheatmap means, NE-high vs -low [PDL1 activation].pdf", width = 15, height = 5)
grid::grid.newpage()
grid::grid.draw(GMA$gtable)
dev.off()
