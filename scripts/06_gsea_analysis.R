# ============================================================
# 06_gsea_analysis.R
# Gene Set Enrichment Analysis using limma-ranked genes
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

# -----------------------------
# Output folders
# -----------------------------

dir.create("results/gsea", recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# User settings
# -----------------------------

species_name <- "Mus musculus"

gene_set_collections <- list(
  Hallmark = list(category = "H", subcategory = NULL),
  Reactome = list(category = "C2", subcategory = "CP:REACTOME"),
  GO_BP = list(category = "C5", subcategory = "GO:BP")
)

# -----------------------------
# Load ranked files
# -----------------------------

ranked_files <- list.files(
  "results/ranked_genes",
  pattern = "\\.rnk$",
  full.names = TRUE
)

stopifnot(length(ranked_files) > 0)

# -----------------------------
# Function to load gene sets
# -----------------------------

get_gene_sets <- function(collection_name) {
  
  collection <- gene_set_collections[[collection_name]]
  
  if (is.null(collection$subcategory)) {
    m_df <- msigdbr(
      species = species_name,
      category = collection$category
    )
  } else {
    m_df <- msigdbr(
      species = species_name,
      category = collection$category,
      subcategory = collection$subcategory
    )
  }
  
  gene_sets <- split(
    m_df$gene_symbol,
    m_df$gs_name
  )
  
  return(gene_sets)
}

# -----------------------------
# Run GSEA for one ranked file
# -----------------------------

run_gsea_one <- function(ranked_file, collection_name) {
  
  cell_type <- basename(ranked_file)
  cell_type <- gsub("_Treatment_vs_Control_limma_tstat.rnk", "", cell_type)
  
  message("Running ", collection_name, " GSEA for: ", cell_type)
  
  ranked_df <- read.table(
    ranked_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  colnames(ranked_df) <- c("Gene", "Score")
  
  ranks <- ranked_df$Score
  names(ranks) <- ranked_df$Gene
  
  ranks <- ranks[!is.na(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)
  
  gene_sets <- get_gene_sets(collection_name)
  
  fgsea_res <- fgsea(
    pathways = gene_sets,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    nperm = 10000
  )
  
  fgsea_res <- fgsea_res %>%
    arrange(padj) %>%
    mutate(
      CellType = cell_type,
      Collection = collection_name
    )
  
  out_file <- file.path(
    "results/gsea",
    paste0(cell_type, "_", collection_name, "_GSEA.csv")
  )
  
  write.csv(
    fgsea_res,
    out_file,
    row.names = FALSE
  )
  
  return(fgsea_res)
}

# -----------------------------
# Run all GSEA analyses
# -----------------------------

all_gsea <- list()

for (collection_name in names(gene_set_collections)) {
  
  for (ranked_file in ranked_files) {
    
    res <- run_gsea_one(
      ranked_file = ranked_file,
      collection_name = collection_name
    )
    
    all_gsea[[paste(collection_name, ranked_file, sep = "_")]] <- res
  }
}

gsea_all <- bind_rows(all_gsea)

write.csv(
  gsea_all,
  "results/gsea/GSEA_all_celltypes_all_collections.csv",
  row.names = FALSE
)

# -----------------------------
# Summary plot
# -----------------------------

plot_df <- gsea_all %>%
  filter(padj < 0.25) %>%
  group_by(Collection, CellType) %>%
  slice_min(order_by = padj, n = 10, with_ties = FALSE) %>%
  ungroup()

if (nrow(plot_df) > 0) {
  
  p <- ggplot(
    plot_df,
    aes(
      x = NES,
      y = reorder(pathway, NES),
      color = padj,
      size = abs(NES)
    )
  ) +
    geom_point() +
    facet_grid(
      Collection ~ CellType,
      scales = "free_y",
      space = "free_y"
    ) +
    theme_bw() +
    labs(
      title = "Top Enriched Pathways Across Cell Types",
      x = "Normalized Enrichment Score",
      y = NULL,
      color = "FDR",
      size = "|NES|"
    )
  
  ggsave(
    "results/plots/gsea_summary_dotplot.png",
    p,
    width = 16,
    height = 12,
    dpi = 300
  )
}

message("Step 06 complete.")
message("GSEA results saved in results/gsea/")