# ============================================================
# 03_fraction_analysis.R
# Cell-type fraction statistics and visualization
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# -----------------------------
# Create output folders
# -----------------------------

dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load inputs
# -----------------------------

theta_final <- readRDS(
  "results/objects/celltype_fractions_final.rds"
)

sample_metadata <- read.csv(
  "results/objects/sample_metadata.csv"
)

# -----------------------------
# Long-format table
# -----------------------------

frac_df <- as.data.frame(theta_final)

frac_df$Sample <- rownames(frac_df)

frac_long <- frac_df %>%
  pivot_longer(
    cols = -Sample,
    names_to = "CellType",
    values_to = "Fraction"
  )

frac_long <- left_join(
  frac_long,
  sample_metadata,
  by = c("Sample" = "public_sample")
)

# -----------------------------
# Summary statistics
# -----------------------------

frac_summary <- frac_long %>%
  group_by(group, CellType) %>%
  summarise(
    Mean = mean(Fraction),
    SD = sd(Fraction),
    N = n(),
    .groups = "drop"
  )

write.csv(
  frac_summary,
  "results/tables/fraction_summary.csv",
  row.names = FALSE
)

# -----------------------------
# Statistical comparison
# -----------------------------

cell_types <- unique(frac_long$CellType)

stats_list <- lapply(cell_types, function(ct){
  
  tmp <- frac_long %>%
    filter(CellType == ct)
  
  pval <- tryCatch(
    t.test(Fraction ~ group, data = tmp)$p.value,
    error = function(e) NA
  )
  
  data.frame(
    CellType = ct,
    PValue = pval
  )
})

fraction_stats <- bind_rows(stats_list)

fraction_stats$FDR <- p.adjust(
  fraction_stats$PValue,
  method = "BH"
)

write.csv(
  fraction_stats,
  "results/tables/fraction_statistics.csv",
  row.names = FALSE
)

# -----------------------------
# Boxplot
# -----------------------------

p1 <- ggplot(
  frac_long,
  aes(
    x = group,
    y = Fraction,
    fill = group
  )
) +
  geom_boxplot(
    alpha = 0.7,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.1,
    size = 2
  ) +
  facet_wrap(
    ~ CellType,
    scales = "free_y"
  ) +
  theme_bw() +
  labs(
    title = "Cell-Type Fraction Comparison",
    x = NULL,
    y = "Estimated Fraction"
  )

ggsave(
  "results/plots/celltype_fraction_boxplots.png",
  p1,
  width = 12,
  height = 8,
  dpi = 300
)

# -----------------------------
# Mean fraction plot
# -----------------------------

mean_df <- frac_long %>%
  group_by(group, CellType) %>%
  summarise(
    Mean = mean(Fraction),
    SD = sd(Fraction),
    .groups = "drop"
  )

p2 <- ggplot(
  mean_df,
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
  theme_bw() +
  labs(
    title = "Mean Cell-Type Fractions",
    y = "Fraction"
  )

ggsave(
  "results/plots/celltype_fraction_barplot.png",
  p2,
  width = 8,
  height = 6,
  dpi = 300
)

message("Step 03 complete.")
message("Fraction tables and plots saved.")