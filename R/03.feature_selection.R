suppressPackageStartupMessages({
  library(tidyverse)
  library(SummarizedExperiment)
  library(limma)
})

# Cibersortx High Resolution Mode only allows for 1000 genes. 
# We will select the one that show most variation in our subgroups 
# of interests in contrast to a baseline group

# Intersection of genes between the sc reference and the bulk, 
# which are used in the downstream deconvolution

gene_list <- read_tsv("data/intermediate/bulk_mixture.txt", col_select = 1) |> 
  pull("...1")

load("data/raw/rnaseq_hersh.Rdata")

subgroups <- levels(se_hersh$sub_group)

# The previously defined molecular subgroups are discarded

old_molecular_groups <- c("del5q-IB", "SF3B1-IB", "Complex")

keep_samples <- !(se_hersh$sub_group %in% old_molecular_groups) & !is.na(se_hersh$sub_group)

se_hersh <- se_hersh[, keep_samples]

# The morphological ones are merged into one to act as a baseline

subgroup_baseline <- se_hersh$sub_group |> 
  fct_recode(
    "morphological" = "MDS-LB",
    "morphological" = "MDS-IB1",
    "morphological" = "MDS-IB2"
  ) |>
  fct_relevel("morphological") |> 
  fct_drop()

counts <- assay(se_hersh, "counts")

lib_sizes <- colSums(counts)

# Transform to log2 CPMs, remove no-variance genes

log2_cpm <- log2((t(t(counts) / lib_sizes) * 1e6) + 1)

log2_cpm <- log2_cpm[apply(log2_cpm, 1, \(row) var(row) != 0 ), ]

# multiclass modeling with limma 

design <- model.matrix( ~ subgroup_baseline, colData(se_hersh))

fit <- lmFit(log2_cpm, design) |> eBayes(trend = TRUE)

# We select the top 1000 genes by bayes F statistic

table <- topTable(fit, number = Inf, sort.by = "F") 

table |> 
  rownames_to_column() |> 
  filter(rowname %in% gene_list) |> 
  head(1000) |> 
  pull(rowname) |> 
  write_lines("data/intermediate/gene_subset.txt")



