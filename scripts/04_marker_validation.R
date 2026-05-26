# ============================================================
# 04_marker_validation.R
# Cell-type marker validation
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# -----------------------------
# Output folders
# -----------------------------

dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load BayesPrism result
# -----------------------------

bp_res <- readRDS(
  "results/objects/bayesprism_result.rds"
)

# -----------------------------
# Extract posterior expression
# -----------------------------

expr <- get.exp(bp_res)

# genes x celltypes matrix
expr_mat <- expr$mu

# -----------------------------
# Example marker list
# Replace with your preferred markers
# -----------------------------

marker_genes <- list(
  Hepatocytes = c("Alb","Apoa1","Ttr"),
  Endothelial = c("Pecam1","Kdr","Vwf"),
  Macrophages = c("Lyz2","Adgre1","Csf1r"),
  Fibrogenic = c("Col1a1","Col1a2","Acta2")
)

marker_df <- do.call(
  rbind,
  lapply(names(marker_genes), function(ct){
    
    data.frame(
      Gene = marker_genes[[ct]],
      ExpectedCellType = ct,
      stringsAsFactors = FALSE
    )
  })
)

# -----------------------------
# Keep genes present
# -----------------------------

marker_df <- marker_df[
  marker_df$Gene %in% rownames(expr_mat),
]

# -----------------------------
# Build plotting table
# -----------------------------

plot_df <- do.call(
  rbind,
  lapply(marker_df$Gene, function(g){
    
    vals <- expr_mat[g, ]
    
    data.frame(
      Gene = g,
      CellType = names(vals),
      Expression = as.numeric(vals)
    )
  })
)

plot_df <- left_join(
  plot_df,
  marker_df,
  by = "Gene"
)

# -----------------------------
# Z-score per gene
# -----------------------------

plot_df <- plot_df %>%
  group_by(Gene) %>%
  mutate(
    Zscore = scale(Expression)[,1]
  ) %>%
  ungroup()

# -----------------------------
# Dotplot
# -----------------------------

p <- ggplot(
  plot_df,
  aes(
    x = CellType,
    y = Gene
  )
) +
  geom_point(
    aes(
      size = Expression,
      color = Zscore
    )
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    title = "Marker Validation",
    x = NULL,
    y = NULL
  )

ggsave(
  "results/plots/marker_validation_dotplot.png",
  p,
  width = 10,
  height = 6,
  dpi = 300
)

# -----------------------------
# Save table
# -----------------------------

write.csv(
  plot_df,
  "results/tables/marker_validation_table.csv",
  row.names = FALSE
)

message("Step 04 complete.")