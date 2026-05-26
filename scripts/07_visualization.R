# ============================================================
# 07_visualization.R
# Final publication-quality visualizations
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# -----------------------------
# Output folder
# -----------------------------

dir.create("results/final_figures",
           recursive = TRUE,
           showWarnings = FALSE)

# ============================================================
# Figure 1
# Cell-type fraction comparison
# ============================================================

fraction_summary <- read.csv(
  "results/tables/fraction_summary.csv"
)

p1 <- ggplot(
  fraction_summary,
  aes(
    x = CellType,
    y = Mean,
    fill = group
  )
) +
  geom_col(
    position = position_dodge()
  ) +
  geom_errorbar(
    aes(
      ymin = Mean - SD,
      ymax = Mean + SD
    ),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Cell-Type Composition",
    x = NULL,
    y = "Estimated Fraction"
  )

ggsave(
  "results/final_figures/Figure1_CellTypeComposition.png",
  p1,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# Figure 2
# Significant pathway counts
# ============================================================

gsea_all <- read.csv(
  "results/gsea/GSEA_all_celltypes_all_collections.csv"
)

pathway_summary <- gsea_all %>%
  filter(padj < 0.25) %>%
  group_by(CellType, Collection) %>%
  summarise(
    Significant_Pathways = n(),
    .groups = "drop"
  )

p2 <- ggplot(
  pathway_summary,
  aes(
    x = CellType,
    y = Significant_Pathways,
    fill = Collection
  )
) +
  geom_col() +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    title = "Significant Pathway Counts",
    y = "Number of Pathways"
  )

ggsave(
  "results/final_figures/Figure2_PathwayCounts.png",
  p2,
  width = 10,
  height = 6,
  dpi = 300
)

# ============================================================
# Figure 3
# Top Hallmark pathways
# ============================================================

hallmark_top <- gsea_all %>%
  filter(
    Collection == "Hallmark",
    padj < 0.25
  ) %>%
  group_by(CellType) %>%
  slice_max(
    order_by = abs(NES),
    n = 10
  ) %>%
  ungroup()

if (nrow(hallmark_top) > 0) {
  
  p3 <- ggplot(
    hallmark_top,
    aes(
      x = NES,
      y = reorder(pathway, NES),
      color = padj,
      size = abs(NES)
    )
  ) +
    geom_point() +
    facet_wrap(
      ~ CellType,
      scales = "free_y"
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = "Top Hallmark Pathways",
      x = "Normalized Enrichment Score",
      y = NULL
    )
  
  ggsave(
    "results/final_figures/Figure3_HallmarkDotplot.png",
    p3,
    width = 14,
    height = 10,
    dpi = 300
  )
}

# ============================================================
# Figure 4
# Differential expression summary
# ============================================================

de_summary <- read.csv(
  "results/tables/DE_summary_all_celltypes.csv"
)

p4 <- ggplot(
  de_summary,
  aes(
    x = reorder(CellType, FDR_0.25),
    y = FDR_0.25
  )
) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Differentially Expressed Genes",
    x = NULL,
    y = "Genes with FDR < 0.25"
  )

ggsave(
  "results/final_figures/Figure4_DEGeneCounts.png",
  p4,
  width = 8,
  height = 6,
  dpi = 300
)

message("Step 07 complete.")
message("Figures saved to results/final_figures/")