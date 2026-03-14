suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(tidyverse)
})

de <- read_csv(
  "data/intermediate/differential_gene_expression_results.csv",
  show_col_types = FALSE
  )

# Rank the gene lists by the t value of limma

ranked_gene_lists <- de |> 
  split(paste(de$cell_type, de$contrast, sep = ";")) |> 
  map(~{
    .x |> 
      drop_na(t) |> 
      arrange(desc(t)) |> 
      pull(t, name = "symbol")
  })

output_file <- "data/intermediate/gene_set_enrichment_results.csv"

gsea_header <- c("cell_type", "contrast", "ID", "Description", "setSize", 
                 "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", 
                 "rank", "leading_edge", "core_enrichment")  

write_lines(paste(gsea_header, collapse = ","), output_file)

# GSEA against Gene Ontology: Biological Processes (BP)

ranked_gene_lists |> 
  imap(~{
    gsea <- gseGO(
    .x, 
    OrgDb = org.Hs.eg.db, 
    ont = "BP",
    keyType = "SYMBOL",
    minGSSize    = 10,    
    maxGSSize    = 500, 
    pvalueCutoff = 1, 
    verbose = FALSE
  )
  # Write results to file
    
  gsea@result |> 
    mutate(
      cell_type = str_split_i(.y, ";", 1),
      contrast =  str_split_i(.y, ";", 2)
      ) |> 
    select(all_of(gsea_header)) |> 
    tibble() |> 
    write_csv(file = output_file, append = TRUE)
    
  }, .progress = "GSEA")
