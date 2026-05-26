# ============================================================
# 05_differential_expression.R
# Cell-type-specific differential expression using limma
# ============================================================

suppressPackageStartupMessages({
  library(BayesPrism)
  library(limma)
  library(dplyr)
  library(tibble)
})

# -----------------------------
# Output folders
# -----------------------------

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/ranked_genes", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load BayesPrism result and metadata
# -----------------------------

bp_res <- readRDS("results/objects/bayesprism_result.rds")
sample_metadata <- read.csv("results/objects/sample_metadata.csv")

# -----------------------------
# Extract cell-type-specific expression
# -----------------------------
# This object contains inferred expression values per cell type.

Z_list <- get.exp(
  bp_res,
  state.or.type = "type"
)

# -----------------------------
# Make sure group order is correct
# -----------------------------

sample_metadata$group <- factor(
  sample_metadata$group,
  levels = c("Control", "Treatment")
)

# -----------------------------
# Differential expression function
# -----------------------------

run_celltype_de <- function(cell_type_name) {
  
  message("Running limma DE for: ", cell_type_name)
  
  Z <- Z_list[[cell_type_name]]
  
  # Z should be samples x genes
  sample_ids <- rownames(Z)
  
  meta_use <- sample_metadata %>%
    filter(public_sample %in% sample_ids)
  
  meta_use <- meta_use[match(sample_ids, meta_use$public_sample), ]
  
  stopifnot(all(meta_use$public_sample == sample_ids))
  
  group <- factor(
    meta_use$group,
    levels = c("Control", "Treatment")
  )
  
  logZ <- log2(Z + 1)
  
  design <- model.matrix(~ group)
  colnames(design) <- c("Intercept", "Treatment_vs_Control")
  
  fit <- lmFit(t(logZ), design)
  fit <- eBayes(fit)
  
  de <- topTable(
    fit,
    coef = "Treatment_vs_Control",
    number = Inf,
    adjust.method = "BH"
  ) %>%
    rownames_to_column("Gene") %>%
    arrange(adj.P.Val) %>%
    mutate(CellType = cell_type_name)
  
  # Save full DE table
  write.csv(
    de,
    file = file.path(
      "results/tables",
      paste0("DE_", cell_type_name, "_Treatment_vs_Control_limma.csv")
    ),
    row.names = FALSE
  )
  
  # Save ranked gene file for GSEA
  ranked <- de %>%
    filter(!is.na(t)) %>%
    select(Gene, t) %>%
    arrange(desc(t))
  
  write.table(
    ranked,
    file = file.path(
      "results/ranked_genes",
      paste0(cell_type_name, "_Treatment_vs_Control_limma_tstat.rnk")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  return(de)
}

# -----------------------------
# Run DE for all cell types
# -----------------------------

DE_list <- lapply(
  names(Z_list),
  run_celltype_de
)

names(DE_list) <- names(Z_list)

DE_all <- bind_rows(DE_list)

write.csv(
  DE_all,
  file = "results/tables/DE_all_celltypes_Treatment_vs_Control_limma.csv",
  row.names = FALSE
)

# -----------------------------
# Summary table
# -----------------------------

de_summary <- DE_all %>%
  group_by(CellType) %>%
  summarise(
    tested_genes = n(),
    nominal_p_0.05 = sum(P.Value < 0.05, na.rm = TRUE),
    FDR_0.25 = sum(adj.P.Val < 0.25, na.rm = TRUE),
    FDR_0.05 = sum(adj.P.Val < 0.05, na.rm = TRUE),
    top_gene = Gene[which.min(adj.P.Val)],
    top_FDR = min(adj.P.Val, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  de_summary,
  file = "results/tables/DE_summary_all_celltypes.csv",
  row.names = FALSE
)

message("Step 05 complete.")
message("Differential expression tables and ranked files saved.")