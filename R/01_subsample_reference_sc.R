suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(Seurat)
  library(GEOquery)
})

set.seed(2026)

# We use the annotated bone marrow single cell atlas GSE245108 as a reference 

geo_accession <- "GSE245108"
sc_data_dir <- file.path("data", "raw", geo_accession)
supp_files <- getGEOSuppFiles(GEO = geo_accession, fetch_files = FALSE, baseDir = sc_data_dir)

# We only use cells from healthy patients
# We collapse their top level of annotation to a higher one, 
# consisting of only 5 classes 

glue("Subsampling dataset")

annotation <- supp_files |> 
  filter(str_detect(fname, "cell-annotation.txt.gz")) |> 
  pull(url) |> 
  read_tsv(show_col_types = FALSE) |> 
  rename_with(str_to_snake) |> 
  mutate(across(contains("level"), str_to_snake)) |> 
  filter(!level_1 %in% c("ba_ma_eo", "stroma")) |> 
  mutate(top_level = case_when(
    level_1 %in% c("megakaryocytic", "erythroid") ~ "megakaryocytic_erythroid", 
    level_1 %in% c("dendritic", "monocytic", "gmp") ~ "dc_monocytic_gmp", 
    level_1 %in% c("early_lymphoid", "b_cell") ~ "early_lymphoid_bcell", 
    TRUE ~ level_1
  )) |> 
  separate(uid, into = c("barcode", "sample"), sep = "\\.", remove = FALSE) 

write_csv(annotation, file.path("data", "intermediate", "sc_annotation.csv"))

# Cibersortx has limited space and memory resources, so
# we apply a stratified subsampling approach to the dataset, to collect a 
# balance and representative sample that can be used by cibersortx

subsampled_df <- annotation |>
  group_by(sample, top_level) |>
  slice_sample(n = 30) |>
  ungroup() |>
  group_by(top_level) |>
  slice_sample(n = 300) |>
  mutate(
    match = str_match(sample, "(\\w\\w\\d\\d)_\\d+_(.+)$"),
    short_sample_name = if_else(
      is.na(match[, 2]),
      str_remove(sample, "_CITE_GEX"),
      paste(match[, 2], match[, 3], sep = "-")
    ), 
    short_uid = paste(barcode, short_sample_name, sep = ".")
  ) |> 
  select(-match)

write_csv(subsampled_df, file.path("data", "intermediate", "sc_annotation_subsampled.csv"))

# Download GEO files

glue("Downloading GEO files")

selected_samples <- unique(subsampled_df$short_sample_name)

supp_files |> 
  mutate(short_sample_name = str_match(
    basename(fname),
    "GSE245108_(.+)_filtered_feature_bc_matrix"
    )[,2]
    ) |> 
  filter(short_sample_name %in% selected_samples) |> 
  pull(url) |> 
  walk( ~ system(glue("wget -q -P {sc_data_dir} -c {.x}")))

# Pipeline to read and preprocess the single cell files

glue("Preprocessing raw files")

seurat_list <- selected_samples |>
  set_names() |> 
  map_chr(~glue("GSE245108_{.x}_filtered_feature_bc_matrix.h5")) |> 
  map_chr(~file.path(sc_data_dir, .x)) |>
  # Read h5 files
  map(quietly( ~ Read10X_h5(.x)), .progress = "read")  |>
  map( ~ pluck(.x, "result")) |>
  map("Gene Expression") |> 
  # Convert to Seurat
  imap(quietly( ~ CreateSeuratObject(counts = .x, project = .y)), .progress = "convert") |>
  map( ~ pluck(.x, "result")) |>
  map( ~ subset(.x, features = rownames(.x)[str_detect(rownames(.x), "Gene Expression")])) |> 
  # Subset to sampled cells
  imap( ~.x[, paste(colnames(.x), .y, sep = ".") %in% subsampled_df$short_uid]) |>
  # Filter by n Features and Mitochondrial percentage
  map(
    ~ subset(
      .x,
      subset =
        nFeature_RNA > 500 &
        nFeature_RNA < 10000 &
        PercentageFeatureSet(.x, pattern = "^MT-") < 10
    ),
    .progress = "filter"
  )

# Merge, normalize to CPMs

glue("Merging data")

merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list)
) |>
  # Relative Count normalization with a scaling factor of 1e6, effectively CPMs
  NormalizeData(
    normalization.method = "RC",
    scale.factor = 1e6,
    verbose = FALSE
  ) |>
  JoinLayers()

# Annotate with cell type

merged_seurat$cell_type <- colnames(merged_seurat) |> 
  map_chr( ~ subsampled_df$top_level[paste(
    subsampled_df$short_sample_name, subsampled_df$barcode, sep = "_"
  ) == .x])

Idents(merged_seurat) <- "cell_type"

# Transform to matrixx

matrix <- GetAssayData(merged_seurat, assay = "RNA", layer = "data") |>
  as.data.frame() |>
  mutate(across(everything(), ~ signif(.x, 2))) |> 
  as.matrix()

header <- paste0(
  rep("\t", length(merged_seurat$cell_type)),
  merged_seurat$cell_type,
  collapse = ""
)  

sc_file <- file.path("data", "intermediate", "sc_reference_all_genes.tsv")

write_lines(header, sc_file)

write.table(
  matrix,
  file = sc_file,
  sep = "\t",
  append = TRUE,
  quote = FALSE,
  col.names = FALSE
) 

glue("Subsampled matrix written to {sc_file}")