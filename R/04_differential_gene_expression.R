suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(SummarizedExperiment)
  library(limma)
})

load("data/raw/rnaseq_hersh.Rdata")

# Sanitize subgroup names, 
# will be reverted back to original in the final results

subgroups_pretty <- levels(colData(se_hersh)$sub_group)
subgroups_sanitized <- str_to_snake(make.names(subgroups_pretty))
subgroups_map <- set_names(subgroups_sanitized, subgroups_pretty)

sample_metadata <- colData(se_hersh) |> 
  as.data.frame() |> 
  select(sub_group) |> 
  drop_na() |> 
  mutate(sub_group = fct_relabel(sub_group, ~subgroups_map[.x])) 

glue("Reading fraction expression matrices")

expression_matrices <- Sys.glob(
  "data/cibersortx_results/CIBERSORTx_Job3_output/CIBERSORTxHiRes_Job3_*_Window20.txt"
) |> 
  set_names(~str_extract(.x, "CIBERSORTxHiRes_Job3_(.+?)_Window20.txt", group = 1)) |> 
  map(~read_tsv(
    .x,
    col_select = c("GeneSymbol", rownames(sample_metadata)),
    col_types = paste(c("c",rep("d", nrow(sample_metadata))), collapse = "")
  )) |> 
  map(~column_to_rownames(.x, "GeneSymbol")) |> 
  map(as.matrix) |> 
  # Genes with only NAs or 0 variance are removed. 
  map(~.x[apply(.x, 1, \(row) sum(!is.na(row)) != 0),]) |> 
  map(~.x[apply(.x, 1, \(row) var(row) != 0),]) |> 
  # Matrices are log_2 transformed
  map(~log2(.x)) |> 
  iwalk(~print(glue("{.y}: dimensions {dim(.x)[1]} x {dim(.x)[2]}")))

# ~ 0 + ... creates a design matrix without an intercept
design <- model.matrix(~ 0 + sub_group, data = sample_metadata)

colnames(design) <- str_remove(colnames(design), "sub_group")

# "One-vs-All" comparison
# Each sub_group is compared against the unweighted average of all other sub_groups
# Ej: 
# ~ ezh_2 - (tet_2_bi + x_7 + stag_2 + del_5_q_ib + sf_3_b_1_ib ... ) / 9

contrast_strings <- subgroups_map |> 
  imap_chr(~ {
    rest <- discard_at(subgroups_map, .y) |> paste(collapse = " + ")
    n_minus_1 <- length(subgroups_sanitized) - 1
    glue("{.x} - ({rest})/{n_minus_1}")
    }
  ) |> 
  set_names(~glue("{.x} vs All"))

contrast_matrix <- do.call(
  makeContrasts,
  c(as.list(contrast_strings), list(levels = design))
)

# Fit all contrasts to all cell type expression matrices

fits <- expression_matrices |> 
  map(~lmFit(.x, design)) |> 
  map(~contrasts.fit(.x, contrast_matrix)) |> 
  map(eBayes)

# Save fitted models to a grid

grid <- expand_grid(
   fit = fits,
   contrast = names(contrast_strings)
 )

# Pairwise (pw) "Subgroup vs baseline" comparisons
# Here we contrast four subgroups of interest against a baseline group

sample_metadata_pw <- sample_metadata |> 
  mutate(sub_group = fct_recode(
    sub_group, 
    "morphological" = "mds_lb", 
    "morphological" = "mds_ib_1", 
    "morphological" = "mds_ib_2"   
    )) 

# All subgroups are included in the model, no intercept

design_pw <- model.matrix(~ 0 + sub_group, data = sample_metadata_pw)
colnames(design_pw) <- str_remove(colnames(design_pw), "sub_group")

# Only contrasts of newly defined subgroups vs morphological will be collected

denovo_subgroups <- c("EZH2", "TET2-bi", "-7", "STAG2") |> 
  set_names() |> 
  map_chr(~subgroups_map[.x])

# These designs are like:
# ~ ezh_2 - morphological (mds_ib1 U mds_ib2 U mds_lb)

contrast_strings_pw <-  denovo_subgroups |> 
  set_names(~glue("{.x} vs morphological")) |> 
  map_chr(~glue("{.x} - morphological "))

contrast_matrix_pw <- do.call(
  makeContrasts,
  c(as.list(contrast_strings_pw), list(levels = design_pw))
)

# Fit all contrasts to all cell type expression matrices

fits_pw <- expression_matrices |> 
  map(~lmFit(.x, design_pw)) |> 
  map(~contrasts.fit(.x, contrast_matrix_pw)) |> 
  map(eBayes)

# Save fitted models to a grid

grid_pw <- expand_grid(
  fit = fits_pw, 
  contrast = names(contrast_strings_pw)
) 

# For each run of the experiment, extract the diff expression 
# values for each subgroup, format and write to file

bind_rows(grid, grid_pw) |> 
  mutate(cell_type = names(fit)) |> 
  mutate(stats = map2(fit, contrast,
    ~rownames_to_column(topTable(.x, coef = .y, number = Inf), "symbol")
    )) |>
  select(-fit) |>
  unnest(stats) |> 
  mutate(across(where(is.double), ~signif(.x, 3))) |> 
  write_csv("data/intermediate/differential_gene_expression_results.csv")
