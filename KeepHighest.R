# Description: KeepHighest value if there are duplicates from a data frame.
# Author: Bell Wu
# Date created: 2025.09.30
# counts = numerical count matrix
# sym = column for gene names

keepHighest = function(counts, sym){
  counts = as.data.frame(counts)
  #Add the gene symbols to a column named 'tmp'
  counts$tmp = sym
  #create object with removed tmp list do use rowMeans function
  object = subset(counts, select = -tmp)
  #Order the counts based on overall expression (highest to lowest)
  counts = counts[order(rowMeans(object), decreasing = T),]
  #Remove any rows that don't have a gene symbol in 'tmp'
  counts = counts[complete.cases(counts$tmp), ]
  #Identify the duplicate gene symbols
  dups = unique(counts$tmp[duplicated(counts$tmp)])
  #If there are duplicates...
  if (length(dups) >= 1) {
    print("Duplicate gene symbols detected. Keeping symbols with highest overall ecountspression...")
    #Use this to store the row numbers to remove
    to_remove = c()
    for (i in 1:length(dups)) {
      #Identify the row numbers of each gene with duplicate symbols
      ind = which(counts$tmp %in% dups[i], arr.ind = T)
      #Store all row numbers except the first (the one with the highest expression) in to_remove
      to_remove = c(to_remove, ind[-1])
    }
    #Remove all the lower expressed rows
    counts = counts[-to_remove, ]
  }
  #Set the rownames to the gene symbols
  rownames(counts) <- counts$tmp
  #Remove the 'tmp' column
  counts <- counts[, -(which(colnames(counts) %in% "tmp"))]
  #Return the counts
  return(counts)
}