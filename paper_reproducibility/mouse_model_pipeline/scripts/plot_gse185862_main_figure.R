#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tidyr)
})

project_root <- "/mnt/data/serinharmanci/preps_revision_052726"
analysis_dir <- file.path(project_root, "analysis_fine_grained_mouse_cortical_subclasses")
figures_dir <- file.path(analysis_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

feature_order <- c(
  "Input resistance (MOhm)",
  "Latency (ms)",
  "AP amplitude (mV)",
  "Rheobase (pA)",
  "Sag ratio",
  "Membrane time constant (ms)",
  "AP threshold (mV)",
  "Upstroke-to-downstroke ratio",
  "ISI adaptation index",
  "AP width (ms)"
)

family_order <- c("IT", "ET", "CT", "Vip", "Sncg", "Lamp5", "Pvalb", "Sst")
display_family <- c(IT = "IT", ET = "ET / PT", CT = "CT", Vip = "Vip", Sncg = "Sncg",
                    Lamp5 = "Lamp5", Pvalb = "Pvalb", Sst = "Sst")
heatmap_order <- c("IT", "ET / PT", "L6 CT CTX", "Vip", "Sncg", "Lamp5", "Pvalb", "Sst")
family_colors <- c(
  IT = "#4C78A8",
  ET = "#F58518",
  CT = "#54A24B",
  Vip = "#E45756",
  Sncg = "#B279A2",
  Lamp5 = "#FF9DA6",
  Pvalb = "#9D755D",
  Sst = "#72B7B2"
)

scale01_by_feature <- function(long_df) {
  long_df %>%
    group_by(feature) %>%
    mutate(
      v_min = min(value, na.rm = TRUE),
      v_max = max(value, na.rm = TRUE),
      scaled = if_else(v_max > v_min, (value - v_min) / (v_max - v_min), 0.5)
    ) %>%
    ungroup() %>%
    select(-v_min, -v_max)
}

matrix_to_long <- function(df, n_cells_df = NULL) {
  out <- df %>%
    select(group, all_of(feature_order)) %>%
    mutate(group = factor(as.character(group), levels = family_order)) %>%
    arrange(group) %>%
    pivot_longer(cols = all_of(feature_order), names_to = "feature", values_to = "value") %>%
    scale01_by_feature() %>%
    mutate(feature = factor(feature, levels = feature_order))

  if (!is.null(n_cells_df)) {
    names(n_cells_df)[1:2] <- c("group", "n_cells")
    out <- out %>% left_join(n_cells_df, by = "group")
  } else {
    out$n_cells <- NA_real_
  }
  out
}

read_m1_matrix <- function() {
  df <- read_csv(
    file.path(analysis_dir, "tables", "DotPlot_M1_patchseq_real_ephys_CTX_M1_family_matched.medians.csv"),
    show_col_types = FALSE
  )
  names(df)[1] <- "group"
  df %>% select(group, all_of(feature_order))
}

read_m1_counts <- function() {
  df <- read_csv(
    file.path(analysis_dir, "tables", "DotPlot_M1_patchseq_real_ephys_CTX_M1_family_matched.ncells.csv"),
    show_col_types = FALSE
  )
  names(df)[1:2] <- c("group", "n_cells")
  df
}

read_previous_gse_ephys <- function() {
  gse_raw <- read_csv(
    file.path(project_root, "results/hierarchical_broad_then_fine/DotPlot_GSE185862_expression_matrix_SSv4_predicted_ephys_by_CTX_neuronal_mapped_labels_mean_ephys_by_label.csv"),
    show_col_types = FALSE
  ) %>%
    filter(mapped_training_RNA_family %in% family_order)

  gse_raw %>%
    group_by(group = mapped_training_RNA_family) %>%
    summarise(
      across(all_of(feature_order), ~ weighted.mean(.x, w = .data$n_cells, na.rm = TRUE)),
      n_cells = sum(n_cells, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(group = factor(group, levels = family_order)) %>%
    arrange(group)
}

make_previous_concordance_heatmap <- function() {
  heat_df <- read_csv(
    file.path(analysis_dir, "results/reviewer_1_comment1_results/GSE185862_CTX_M1_family_matched_L6CT_only_concordance_no_PCA_heatmap_values.csv"),
    show_col_types = FALSE
  ) %>%
    mutate(
      row_label = factor(row_label, levels = rev(heatmap_order)),
      column_label = factor(column_label, levels = c("PREPS\np >= 0.0\nall cells", "PREPS\np >= 0.5\nfiltered")),
      label = paste0(sprintf("%.1f%%", concordance_pct), "\n", "n=", format(n_cells, big.mark = ","))
    )

  ggplot(heat_df, aes(x = column_label, y = row_label, fill = concordance_pct)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = label), color = "black", size = 6.1, lineheight = 0.95) +
    scale_fill_gradientn(
      colors = c("#b10026", "#f46d43", "#ffffbf", "#66bd63", "#006837"),
      limits = c(0, 100),
      name = "Concordance (%)"
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 19) +
    theme(
      axis.text.x = element_text(face = "bold", size = 15),
      axis.text.y = element_text(size = 18),
      panel.grid = element_blank(),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15),
      legend.key.height = unit(1.3, "cm"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

make_dotplot <- function(long_df, title) {
  label_df <- long_df %>%
    distinct(group, n_cells) %>%
    mutate(
      group_chr = as.character(group),
      group_label = paste0(display_family[group_chr], " (n=", format(round(n_cells), big.mark = ","), ")")
    )

  plot_df <- long_df %>%
    mutate(group_chr = as.character(group)) %>%
    left_join(label_df %>% select(group_chr, group_label), by = "group_chr") %>%
    mutate(group_label = factor(group_label, levels = rev(label_df$group_label)))

  ggplot(plot_df, aes(x = feature, y = group_label)) +
    geom_point(aes(size = scaled, fill = scaled), shape = 21, color = "grey25", stroke = 0.2) +
    scale_fill_viridis_c(option = "magma", begin = 0.15, end = 0.85, limits = c(0, 1), name = "Average Score") +
    scale_size(range = c(2.8, 12.0), limits = c(0, 1), guide = "none") +
    labs(x = NULL, y = NULL, title = title) +
    theme_bw(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", size = 21, hjust = 0),
      axis.text.x = element_text(angle = 42, hjust = 1, vjust = 1, size = 15),
      axis.text.y = element_text(size = 16),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      legend.position = "right"
    )
}

make_scatter_panel <- function(m1_df, gse_df) {
  m1_long <- matrix_to_long(m1_df) %>% select(group, feature, m1_scaled = scaled)
  gse_long <- matrix_to_long(gse_df) %>% select(group, feature, gse_scaled = scaled)
  scatter_df <- inner_join(m1_long, gse_long, by = c("group", "feature"))
  pear <- suppressWarnings(cor.test(scatter_df$m1_scaled, scatter_df$gse_scaled, method = "pearson"))
  spear <- suppressWarnings(cor.test(scatter_df$m1_scaled, scatter_df$gse_scaled, method = "spearman"))
  label <- sprintf("Pearson r = %.2f\nSpearman rho = %.2f\nn = %d", unname(pear$estimate), unname(spear$estimate), nrow(scatter_df))
  write_csv(
    tibble(
      comparison = "display_scaled01_flattened",
      n = nrow(scatter_df),
      pearson_r = unname(pear$estimate),
      pearson_p = pear$p.value,
      spearman_rho = unname(spear$estimate),
      spearman_p = spear$p.value
    ),
    file.path(analysis_dir, "results/gse185862_current_16k_model/M1_vs_GSE185862_previous_models_ephys_correlation_summary.csv")
  )

  ggplot(scatter_df, aes(x = m1_scaled, y = gse_scaled, color = group)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.55) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.75) +
    geom_point(size = 3.8, alpha = 0.9) +
    scale_color_manual(values = family_colors, labels = display_family, drop = FALSE) +
    annotate("label", x = 0.03, y = 0.97, label = label, hjust = 0, vjust = 1, size = 5.1) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = TRUE) +
    labs(
      x = "M1 measured ephys (scaled 0-1)",
      y = "GSE185862 predicted ephys (scaled 0-1)",
      color = "Family",
      title = "M1 vs GSE185862 family-profile correlation"
    ) +
    theme_bw(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", size = 21, hjust = 0),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 17),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15),
      panel.grid.minor = element_blank()
    )
}

m1_df <- read_m1_matrix()
m1_n <- read_m1_counts()
gse_all <- read_previous_gse_ephys()
gse_df <- gse_all %>% select(group, all_of(feature_order))
gse_n <- gse_all %>% select(group, n_cells)

p_concordance <- make_previous_concordance_heatmap()
p_scatter <- make_scatter_panel(m1_df, gse_df)
p_m1 <- make_dotplot(matrix_to_long(m1_df, m1_n), "M1 Patch-seq measured ephys")
p_gse <- make_dotplot(matrix_to_long(gse_df, gse_n), "GSE185862 PREPS-predicted ephys")

final_plot <- (p_concordance + p_scatter) / (p_m1 + p_gse) +
  plot_layout(widths = c(1.25, 1.15), heights = c(1.05, 1.05)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 26))

for (name in c(
  "Figure_GSE185862_previous_models_2x2_concordance_ephys",
  "Figure_GSE185862_current_16k_2x2_concordance_ephys",
  "R1MP1_supp_GSE185862_ephys_dotplots_R",
  "R1MP1_GSE185862_2x2_concordance_ephys_R"
)) {
  out_prefix <- file.path(figures_dir, name)
  ggsave(paste0(out_prefix, ".png"), final_plot, width = 21.6, height = 14.8, dpi = 300, bg = "white")
  ggsave(paste0(out_prefix, ".pdf"), final_plot, width = 21.6, height = 14.8, bg = "white")
  message("[saved] ", paste0(out_prefix, ".png"))
}
