# ============================================================
# 02_run_bayesprism.R
# Run BayesPrism deconvolution using processed bulk/reference
# ============================================================

suppressPackageStartupMessages({
  library(BayesPrism)
})

# -----------------------------
# Create output folders
# -----------------------------

dir.create("results", showWarnings = FALSE)
dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load processed objects from Step 01
# -----------------------------

bulk_bp <- readRDS("results/objects/bulk_for_bayesprism.rds")
sc_bp <- readRDS("results/objects/reference_for_bayesprism.rds")
cell_type_labels <- readRDS("results/objects/cell_type_labels.rds")
sample_metadata <- read.csv("results/objects/sample_metadata.csv")

# -----------------------------
# Sanity checks
# -----------------------------

message("Bulk matrix dimensions:")
print(dim(bulk_bp))

message("Reference matrix dimensions:")
print(dim(sc_bp))

message("Number of cell labels:")
print(length(cell_type_labels))

message("Cell type counts:")
print(table(cell_type_labels))

stopifnot(all(colnames(bulk_bp) == colnames(sc_bp)))
stopifnot(length(cell_type_labels) == nrow(sc_bp))

# -----------------------------
# Create BayesPrism object
# -----------------------------

message("Creating BayesPrism object...")

myPrism <- new.prism(
  reference = sc_bp,
  mixture = bulk_bp,
  input.type = "count.matrix",
  cell.type.labels = cell_type_labels,
  cell.state.labels = cell_type_labels,
  key = NULL,
  outlier.cut = 0.01,
  outlier.fraction = 0.1
)

print(myPrism)

# -----------------------------
# Run BayesPrism
# -----------------------------

message("Running BayesPrism...")

bp_res <- run.prism(
  prism = myPrism,
  n.cores = 4
)

# -----------------------------
# Extract final cell-type fractions
# -----------------------------

theta_final <- get.fraction(
  bp_res,
  which.theta = "final",
  state.or.type = "type"
)

message("Final cell-type fraction matrix:")
print(dim(theta_final))
print(head(theta_final))

# -----------------------------
# Save results
# -----------------------------

saveRDS(
  bp_res,
  file = "results/objects/bayesprism_result.rds"
)

saveRDS(
  theta_final,
  file = "results/objects/celltype_fractions_final.rds"
)

write.csv(
  theta_final,
  file = "results/tables/celltype_fractions_final.csv",
  row.names = TRUE
)

# -----------------------------
# Group mean summary
# -----------------------------

theta_df <- as.data.frame(theta_final)
theta_df$Sample <- rownames(theta_df)

theta_df <- merge(
  sample_metadata,
  theta_df,
  by.x = "public_sample",
  by.y = "Sample"
)

fraction_summary <- aggregate(
  theta_df[, colnames(theta_final)],
  by = list(Group = theta_df$group),
  FUN = mean
)

write.csv(
  fraction_summary,
  file = "results/tables/celltype_fraction_group_means.csv",
  row.names = FALSE
)

message("Step 02 complete.")
message("Saved BayesPrism result and cell-type fractions.")