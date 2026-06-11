#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

project_root <- "/mnt/data/serinharmanci/preps_revision_052726"
analysis_dir <- file.path(project_root, "analysis_fine_grained_mouse_cortical_subclasses")
pipe_dir <- file.path(analysis_dir, "mouse_model_pipeline")
results_dir <- file.path(pipe_dir, "results", "pca_baseline_regenerated")
figures_dir <- file.path(pipe_dir, "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

family_order <- c("IT", "ET", "CT", "Vip", "Sncg", "Lamp5", "Pvalb", "Sst")

assign_true_family <- function(subclass_clean) {
  case_when(
    grepl("(^| )IT( |-|$)", subclass_clean) ~ "IT",
    grepl("(^| )ET( |-|$)", subclass_clean) ~ "ET",
    grepl("(^| )CT( |-|$)|L6b", subclass_clean) ~ "CT",
    grepl("Vip", subclass_clean) ~ "Vip",
    grepl("Sncg", subclass_clean) ~ "Sncg",
    grepl("Lamp5", subclass_clean) ~ "Lamp5",
    grepl("Pvalb", subclass_clean) ~ "Pvalb",
    grepl("Sst", subclass_clean) ~ "Sst",
    TRUE ~ NA_character_
  )
}

make_concordance_summary <- function(pca_df, cutoffs = c(0.0, 0.5)) {
  cutoff_labels <- c(
    `0` = "p >= 0.0\nall cells",
    `0.5` = "p >= 0.5\nfiltered"
  )

  bind_rows(lapply(cutoffs, function(cutoff) {
    pca_df %>%
      filter(prob >= cutoff) %>%
      group_by(true_family) %>%
      summarise(
        threshold = cutoff_labels[as.character(cutoff)],
        n_cells = n(),
        n_correct = sum(pred_rna_family == true_family),
        concordance = 100 * n_correct / n_cells,
        .groups = "drop"
      )
  })) %>%
    mutate(
      family = if_else(as.character(true_family) == "ET", "ET / PT", as.character(true_family)),
      family = factor(family, levels = rev(c("IT", "ET / PT", "CT", "Vip", "Sncg", "Lamp5", "Pvalb", "Sst"))),
      threshold = factor(threshold, levels = unname(cutoff_labels[as.character(cutoffs)])),
      label = paste0(sprintf("%.1f%%", concordance), "\n", "n=", format(n_cells, big.mark = ","))
    )
}

make_heatmap <- function(heat_df) {
  ggplot(heat_df, aes(x = threshold, y = family, fill = concordance)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = label), color = "white", size = 6.5, lineheight = 0.95) +
    scale_fill_gradientn(
      colors = c("#b10026", "#f46d43", "#ffffbf", "#66bd63", "#006837"),
      limits = c(0, 100),
      name = "Concordance (%)"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = "PCA baseline WMB concordance",
      subtitle = "CTX M1-family matched cells; train/eval fine-tuning cells removed"
    ) +
    theme_minimal(base_size = 19) +
    theme(
      plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
      plot.subtitle = element_text(size = 18, hjust = 0.5),
      axis.text.x = element_text(size = 18, face = "bold"),
      axis.text.y = element_text(size = 18),
      panel.grid = element_blank(),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15),
      legend.key.height = unit(1.3, "cm"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

pca_path <- file.path(analysis_dir, "results", "pca_baseline", "pca_wmb_combined_predictions.csv")
# full 35,493-cell fine-tune+eval exclusion (matches the corrected main Figure 5)
exclude_path <- file.path(pipe_dir, "results", "wmb_finetune_exclusion", "wmb_subclass_finetune_eval_cell_ids_35493.csv")

if (!file.exists(pca_path)) stop("Missing PCA predictions: ", pca_path)
if (!file.exists(exclude_path)) stop("Missing fine-tuning exclusion list: ", exclude_path)

exclude_ids <- read_csv(exclude_path, show_col_types = FALSE) %>%
  pull(cell_id) %>%
  unique()

pca_df <- read_csv(pca_path, show_col_types = FALSE) %>%
  transmute(
    cell_id = cell_label,
    pred_rna_family,
    prob,
    subclass_clean,
    class,
    true_family = assign_true_family(subclass_clean)
  ) %>%
  filter(!cell_id %in% exclude_ids, !is.na(true_family)) %>%
  mutate(true_family = factor(true_family, levels = family_order))

heat_df <- make_concordance_summary(pca_df)
write_csv(heat_df, file.path(results_dir, "WMB_PCA_family_concordance_excluding_finetune_cells_regenerated.csv"))

all_cells <- heat_df %>%
  filter(threshold == "p >= 0.0\nall cells") %>%
  transmute(true_family, n_cells, n_correct, concordance_pct = concordance)
write_csv(all_cells, file.path(results_dir, "WMB_PCA_family_concordance_all_cells_only_regenerated.csv"))

plot <- make_heatmap(heat_df)
out_prefix <- file.path(figures_dir, "Figure_WMB_PCA_baseline_concordance")
ggsave(paste0(out_prefix, ".png"), plot, width = 11.5, height = 7.5, dpi = 300, bg = "white")
ggsave(paste0(out_prefix, ".pdf"), plot, width = 11.5, height = 7.5, bg = "white")

message("[saved] ", paste0(out_prefix, ".png"))
