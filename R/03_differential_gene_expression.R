suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(SummarizedExperiment)
  library(limma)
})

load("data/raw/rnaseq_hersh.Rdata")

sample_metadata <- colData(se_hersh) |> 
  as.data.frame() |> 
  select(sub_group) |> 
  drop_na()

# Sanitize subgroup names, will be reverted back at in the final results

level_name_fix_map <- levels(sample_metadata$sub_group) |> 
  set_names() |> 
  map_chr(~str_to_snake(make.names(.x)))

sample_metadata <- mutate(sample_metadata, sub_group = level_name_fix_map[sub_group])

glue("Reading fraction expression matrices")

# Genes with only NAs or 0 variance are removed. 
# Matrices are log_2 transformed

expression_matrices <- Sys.glob("data/cibersortx_results/CIBERSORTx_Job3_output/CIBERSORTxHiRes_Job3_*_Window20.txt") |> 
  set_names(~str_extract(.x, "CIBERSORTxHiRes_Job3_(.+?)_Window20.txt", group = 1)) |> 
  map(~read_tsv(
    .x,
    col_select = c("GeneSymbol", rownames(sample_metadata)),
    col_types = paste(c("c",rep("d", nrow(sample_metadata))), collapse = "")
    )) |> 
  map(~column_to_rownames(.x, "GeneSymbol")) |> 
  map(as.matrix) |> 
  map(~.x[apply(.x, 1, \(row) sum(!is.na(row)) != 0),]) |> 
  map(~.x[apply(.x, 1, \(row) var(row) != 0),]) |> 
  map(~log2(.x)) |> 
  iwalk(~print(glue("{.y}: dimensions {dim(.x)[1]} x {dim(.x)[2]}")))

coldata <- map_dfc(set_names(level_name_fix_map), ~sample_metadata$sub_group == .x)

diff_exp_grid <- expand_grid(
  cell_type = names(expression_matrices),
  subgroup = level_name_fix_map
  )

glue("Calculating differential expression across {length(expression_matrices)} cell type matrices ", 
     "and {length(level_name_fix_map)} contrasts. ({nrow(diff_exp_grid)} experiments)")

limma_results <- diff_exp_grid |> 
  asplit(1) |>
  head(15) |> 
  map(~{
    subgroup <- .x[['subgroup']]
    cell_type <- .x[['cell_type']]
    formula <- reformulate(subgroup)
    design <- model.matrix(formula, data = coldata)
    fit <- lmFit(expression_matrices[[cell_type]], design) 
    fit <- eBayes(fit, trend = TRUE)
    topTable(fit, number = Inf, coef = 2) |> 
      rownames_to_column("symbol") |> 
      mutate(
        cell_type = cell_type,
        subgroup = names(level_name_fix_map[level_name_fix_map ==  subgroup]), 
        contrast = glue("{subgroup} vs rest")
        )
    })  |> 
  bind_rows() 

limma_results |> 
  select(cell_type, contrast, symbol, logFC, t, P.Value, adj.P.Val) |> 
  mutate(across(where(is.double), ~signif(.x, 3))) |> 
  write_csv("results/differential_gene_expression_results.csv")

