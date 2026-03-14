suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(SummarizedExperiment)
})

all_genes_sc_file <- "data/intermediate/sc_reference_all_genes.tsv"

all_genes_sc <- read_tsv(
  all_genes_sc_file,
  name_repair = "minimal",
  show_col_types = FALSE
  ) 

colnames(all_genes_sc)[1] <- "rowname"

all_genes_sc <- all_genes_sc |> column_to_rownames() |> as.matrix()

# Bulk matrix from Hershberger et al. 

load("data/raw/rnaseq_hersh.Rdata")

counts <- assay(se_hersh, "counts")

# Convert to CPMs

all_genes_bulk <- counts |> 
  asplit(2) |> 
  map2(colSums(counts), ~(.x / .y) * 1e6)  |>
  map(~signif(.x, 2)) %>%
  do.call(cbind, .)

# Find intersection of genes

common_genes <- sort(intersect(rownames(all_genes_sc), rownames(all_genes_bulk)))

# Subset sc matrix and write 

sc <- all_genes_sc[common_genes, ]

header_sc <- paste0( rep("\t", ncol(sc)), colnames(sc), collapse = "")

sc_file <- file.path("data", "intermediate", "sc_reference_subsetted.txt")

write_lines(header_sc, sc_file)

write.table(
  sc,
  file = sc_file,
  sep = "\t",
  append = TRUE,
  quote = FALSE,
  col.names = FALSE
) 

glue("Single cell reference matrix written to {sc_file}")

# Write bulk file

bulk <- all_genes_bulk[common_genes, ]

bulk_file <- file.path("data", "intermediate", "bulk_mixture.txt")

header_bulk <-  paste0( rep("\t", ncol(bulk)), colnames(bulk), collapse = "")

write_lines(header_bulk, bulk_file)

write.table(
  bulk,
  file = bulk_file,
  sep = "\t",
  append = TRUE, 
  quote = FALSE, 
  col.names = FALSE
) 

glue("Bulk mixture matrix written to {bulk_file}")


