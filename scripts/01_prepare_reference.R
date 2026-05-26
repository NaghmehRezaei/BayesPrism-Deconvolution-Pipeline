# ============================================================
# 01_prepare_reference.R
# Prepare bulk RNA-seq matrix and single-cell reference matrix
# for BayesPrism deconvolution
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# -----------------------------
# User settings
# -----------------------------

bulk_count_file <- "data/bulk_counts.csv"
reference_file  <- "data/single_cell_reference.txt"

sample_metadata <- data.frame(
  original_sample = c(
    "Control_Rep1", "Control_Rep2", "Control_Rep3",
    "Treatment_Rep1", "Treatment_Rep2", "Treatment_Rep3"
  ),
  public_sample = c(
    "Control_Rep1", "Control_Rep2", "Control_Rep3",
    "Treatment_Rep1", "Treatment_Rep2", "Treatment_Rep3"
  ),
  group = c(
    rep("Control", 3),
    rep("Treatment", 3)
  ),
  stringsAsFactors = FALSE
)

bulk_gene_column <- "name"
reference_gene_column <- "GeneSymbol"

# -----------------------------
# Load bulk RNA-seq counts
# -----------------------------

bulk_raw <- read.csv(
  bulk_count_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(bulk_gene_column %in% colnames(bulk_raw))
stopifnot(all(sample_metadata$original_sample %in% colnames(bulk_raw)))

bulk <- bulk_raw[, c(bulk_gene_column, sample_metadata$original_sample)]

colnames(bulk) <- c("GeneSymbol", sample_metadata$public_sample)

bulk <- bulk |>
  filter(!is.na(GeneSymbol), GeneSymbol != "") |>
  group_by(GeneSymbol) |>
  summarise(
    across(
      all_of(sample_metadata$public_sample),
      ~ sum(as.numeric(.x), na.rm = TRUE)
    ),
    .groups = "drop"
  )

bulk_mat <- as.data.frame(bulk)
rownames(bulk_mat) <- bulk_mat$GeneSymbol
bulk_mat$GeneSymbol <- NULL
bulk_mat <- as.matrix(bulk_mat)
mode(bulk_mat) <- "numeric"

# -----------------------------
# Load single-cell reference
# -----------------------------

ref_raw <- read.delim(
  reference_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(reference_gene_column %in% colnames(ref_raw))

ref_genes <- ref_raw[[reference_gene_column]]

ref_mat <- ref_raw[, setdiff(colnames(ref_raw), reference_gene_column)]
ref_mat <- as.matrix(ref_mat)
mode(ref_mat) <- "numeric"
rownames(ref_mat) <- ref_genes

# -----------------------------
# Extract cell-type labels
# -----------------------------
# Assumption:
# Reference columns are named like:
# Hepatocytes.1, Hepatocytes.2, Endothelial_cells.1, etc.

cell_type_labels <- gsub("\\.[0-9]+$", "", colnames(ref_mat))

# -----------------------------
# Match genes between bulk and reference
# -----------------------------

common_genes <- intersect(rownames(bulk_mat), rownames(ref_mat))

message("Number of common genes: ", length(common_genes))

bulk_use <- bulk_mat[common_genes, , drop = FALSE]
ref_use  <- ref_mat[common_genes, , drop = FALSE]

keep_genes <- rowSums(bulk_use) > 0 & rowSums(ref_use) > 0

bulk_use <- bulk_use[keep_genes, , drop = FALSE]
ref_use  <- ref_use[keep_genes, , drop = FALSE]

message("Number of genes after filtering: ", nrow(bulk_use))

# -----------------------------
# Format for BayesPrism
# -----------------------------
# BayesPrism expects:
# mixture/reference as samples or cells x genes

bulk_bp <- t(bulk_use)
sc_bp   <- t(ref_use)

stopifnot(all(colnames(bulk_bp) == colnames(sc_bp)))
stopifnot(length(cell_type_labels) == nrow(sc_bp))

# -----------------------------
# Save processed objects
# -----------------------------

saveRDS(
  bulk_bp,
  file = "results/objects/bulk_for_bayesprism.rds"
)

saveRDS(
  sc_bp,
  file = "results/objects/reference_for_bayesprism.rds"
)

saveRDS(
  cell_type_labels,
  file = "results/objects/cell_type_labels.rds"
)

write.csv(
  sample_metadata,
  file = "results/objects/sample_metadata.csv",
  row.names = FALSE
)

message("Step 01 complete.")
message("Saved processed BayesPrism input objects in results/objects/")
