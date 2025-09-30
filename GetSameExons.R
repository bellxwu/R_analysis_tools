# Description: Function to directly select exons from ensembl API.
# Used to identify common sequences within exons of multiple ensembl transcripts.
# For final design of sgRNA to KO selected exon for sequence.
# Date created: 2025.09.14
# Author: Bell Wu

# Required packages: httr2, jsonlite, tidyverse, rlang
# Requires ENSG ID input. Thus need to run get_ensembl to get ensembl ID 

# Changelog:
# 2025.09.29 Combine get_ensembl with GetSameExons (so all one function and not need to run both)

# 1.0 Combines get wrapper to select exons: -------------------------------------------------
GetSameExons = function(gene_id) {

  # find ensembl gene ID for gene of choice:    
  path = paste0("xrefs/symbol/homo_sapiens/", gene_id)
  get_ensembl = function(path, query = list()) {
    request("http://rest.ensembl.org") |> # set the main URL to request
      req_url_path(path) |> # use path to alter the path of the main URL
      req_headers("Accept" = "application/json") |> # accept only JSON files
      req_url_query(!!!query) |> # modify query components for GET URL
      req_perform() |>  # perform request
      resp_body_json(simplifyVector = TRUE)
  }
  # get path for gene_id
  gene_id = get_ensembl(path)
  gene_id = gene_id$id
  
  # get all known exons for gene ensembl ID   
  get_ensembl_exon = function(gene_id) {
    query = list("feature" = "exon")
    get_ensembl(paste0("overlap/id/", gene_id),
                query)
  }
  exon_df = get_ensembl_exon(gene_id = gene_id)
  tx_exons  <- split(exon_df$id, exon_df$Parent) # split by transcripts
  common_id <- Reduce(intersect, tx_exons) # find the exon present in all IDs 
  if (length(common_id) == 0) {
    # unlist and table how exon transcript appearance
    t = table(unlist(tx_exons)) |> 
      data.frame() 
    
    # select exon transcripts in top 10% freq
    exon_t = t |> 
      dplyr::filter(Freq >= quantile(t$Freq, 0.9)) # find above 90% quantile
    # create exon_ids
    exon_ids = exon_t$Var1
    colnames(exon_t)[colnames(exon_t) == "Var1"] = "exon_id" # rename column 
    # pull exon data from ensembl REST
    exons = request("http://rest.ensembl.org") |> # set the main URL to request
      req_url_path("lookup/id") |> # use path to alter the path of the main URL
      req_headers("Accept" = "application/json") |> # accept only JSON files
      req_body_json(list("ids" = exon_ids)) |> # modify query components for GET URL
      req_perform() |>  # perform request
      resp_body_json(simplifyVector = TRUE)
    
    # find exon starts to determine early exons
    exon_starts = lapply(exons, function(x) x[["start"]])
    early_exon = unlist(exon_starts)
    # add exon strand information to table
    # create a df for the exon id with strands
    strand_ids = data.frame(exon_id = exon_df$exon_id,
                            strands = exon_df$strand)
    strand_ids = strand_ids[!duplicated(strand_ids$exon_id), ] # remove duplicates
    # match the ids to the tabled results
    exon_t = dplyr::inner_join(exon_t, strand_ids, by = "exon_id") # join the dfs 
    # add exon starts to the table
    exon_t$Start = early_exon[match(exon_t$exon_id, names(early_exon))]
    # if there are two strand indicators, remove the transcripts of fewest strands
    count_strands = exon_t |> 
      dplyr::count(strands) |> 
      dplyr::filter(n == max(n))
    # remove fewest
    exon_t = exon_t |> 
      dplyr::filter(strands == count_strands$strands)
    
    # create if-else statement for if -1 or +1 strand
    if (unique(exon_t$strands == -1)) { # if reverse strand
      exon_t = exon_t[order(exon_t$Start, decreasing = TRUE), ] # order by ascending order of exons
    } else {
      exon_t = exon_t[order(exon_t$Start, decreasing = FALSE), ]
    }
    return(exon_t)
    
  } else { 
    return(common_id)
  }
}