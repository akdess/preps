#!/usr/bin/env Rscript
# CANONICAL paper Figure 5 (WMB) generator.
# Same colors, panels, and layout as the original, but the validation set EXCLUDES the
# FULL 35,493 fine-tune+eval cells (24-subclass Mouse-Geneformer fine-tune), not the
# incomplete 16,000-cell list. Output is the official Figure5_WMB_concordance_ephys.{png,pdf}.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork); library(readr); library(tidyr)
})

project_root <- "/mnt/data/serinharmanci/preps_revision_052726"
comment_dir <- file.path(project_root, "reviewer_1_comment1")
clean_dir <- file.path(project_root, "reviewer_1_major_point1_clean")
pipe_dir <- file.path(project_root, "analysis_fine_grained_mouse_cortical_subclasses/mouse_model_pipeline")
figures_dir <- file.path(pipe_dir, "figures")
clean_figures_dir <- file.path(clean_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# *** only change vs original: full 35,493-cell exclusion list ***
finetune_ids_path <- file.path(pipe_dir, "results/wmb_finetune_exclusion/wmb_subclass_finetune_eval_cell_ids_35493.csv")

feature_order <- c(
  "Input resistance (MOhm)", "Latency (ms)", "AP amplitude (mV)", "Rheobase (pA)",
  "Sag ratio", "Membrane time constant (ms)", "AP threshold (mV)",
  "Upstroke-to-downstroke ratio", "ISI adaptation index", "AP width (ms)"
)
family_order <- c("IT", "ET", "CT", "Vip", "Sncg", "Lamp5", "Pvalb", "Sst")
family_colors <- c(IT = "#4C78A8", ET = "#F58518", CT = "#54A24B", Vip = "#E45756",
                   Sncg = "#B279A2", Lamp5 = "#FF9DA6", Pvalb = "#9D755D", Sst = "#72B7B2")

assign_wmb_family <- function(subclass) {
  case_when(
    grepl("^006 L4/5 IT CTX Glut|^007 L2/3 IT CTX Glut|^004 L6 IT CTX Glut|^005 L5 IT CTX Glut", subclass) ~ "IT",
    grepl("^022 L5 ET CTX Glut", subclass) ~ "ET",
    grepl("^029 L6b CTX Glut|^030 L6 CT CTX Glut", subclass) ~ "CT",
    grepl("^046 Vip Gaba", subclass) ~ "Vip",
    grepl("^047 Sncg Gaba", subclass) ~ "Sncg",
    grepl("^049 Lamp5 Gaba|^050 Lamp5 Lhx6 Gaba", subclass) ~ "Lamp5",
    grepl("^051 Pvalb chandelier Gaba|^052 Pvalb Gaba", subclass) ~ "Pvalb",
    grepl("^053 Sst Gaba", subclass) ~ "Sst",
    TRUE ~ NA_character_
  )
}

finetune_ids <- read_csv(finetune_ids_path, show_col_types = FALSE) %>% pull(cell_id) %>% unique()

load_filtered_wmb <- function(pred_subdirs) {
  pred_paths <- file.path(comment_dir, "preps_workflow_ctx_m1_family", pred_subdirs, "all_predictions.csv")
  wmb_pred <- bind_rows(lapply(pred_paths, read_csv, show_col_types = FALSE)) %>%
    filter(!cell_id %in% finetune_ids)
  meta_path <- file.path(project_root,
    "results/preps_ephys_mouse_geneformer__ft_wmb_isocortex_subclass_2000perCellType__hidden+logits+probs",
    "wmb_predictions/wmb_isocortex_combined_predicted_ephys_embs.csv")
  wmb_meta <- read_csv(meta_path, show_col_types = FALSE) %>%
    select(cell_id, subclass) %>%
    mutate(true_family = assign_wmb_family(subclass)) %>%
    filter(!is.na(true_family)) %>%
    distinct(cell_id, true_family)
  wmb_pred %>% inner_join(wmb_meta, by = "cell_id") %>%
    mutate(true_family = factor(true_family, levels = family_order))
}

scale01_by_feature <- function(long_df) {
  long_df %>% group_by(feature) %>%
    mutate(v_min = min(value, na.rm = TRUE), v_max = max(value, na.rm = TRUE),
           scaled = if_else(v_max > v_min, (value - v_min) / (v_max - v_min), 0.5)) %>%
    ungroup() %>% select(-v_min, -v_max)
}

read_matrix <- function(path, family_col_candidates = c("family", "RNA family", "subclass_label")) {
  df <- read_csv(path, show_col_types = FALSE)
  family_col <- family_col_candidates[family_col_candidates %in% names(df)][1]
  df %>% rename(group = all_of(family_col))
}

matrix_to_long <- function(df, group_order = NULL, n_cells_df = NULL) {
  out <- df %>% select(group, all_of(feature_order)) %>%
    pivot_longer(cols = all_of(feature_order), names_to = "feature", values_to = "value") %>%
    scale01_by_feature() %>%
    mutate(feature = factor(feature, levels = feature_order),
           group = if (is.null(group_order)) factor(group, levels = unique(group)) else factor(group, levels = group_order))
  if (!is.null(n_cells_df)) { names(n_cells_df)[1:2] <- c("group", "n_cells"); out <- out %>% left_join(n_cells_df, by = "group") }
  else { out$n_cells <- NA_real_ }
  out
}

make_concordance_summary <- function(wmb_filtered, cutoffs = c(0.0)) {
  summarize_cutoff <- function(cutoff, label) {
    wmb_filtered %>% filter(RNA_family_prob >= cutoff) %>%
      group_by(true_family) %>%
      summarise(threshold = label, n_cells = n(), n_correct = sum(RNA_family == true_family),
                concordance = 100 * n_correct / n_cells, .groups = "drop")
  }
  cutoff_labels <- c(`0` = "p >= 0.0\nall cells", `0.5` = "p >= 0.5\nfiltered")
  bind_rows(lapply(cutoffs, function(cutoff) summarize_cutoff(cutoff, cutoff_labels[as.character(cutoff)]))) %>%
    mutate(family = if_else(as.character(true_family) == "ET", "ET / PT", as.character(true_family)),
           family = factor(family, levels = rev(c("IT", "ET / PT", "CT", "Vip", "Sncg", "Lamp5", "Pvalb", "Sst"))),
           threshold = factor(threshold, levels = unname(cutoff_labels[as.character(cutoffs)])),
           label = paste0(sprintf("%.1f%%", concordance), "\n", "n=", format(n_cells, big.mark = ",")))
}

make_concordance_heatmap <- function(heat_df, title, subtitle) {
  ggplot(heat_df, aes(x = threshold, y = family, fill = concordance)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = label), color = "white", size = 6.5, lineheight = 0.95) +
    scale_fill_gradientn(colors = c("#b10026", "#f46d43", "#ffffbf", "#66bd63", "#006837"),
                         limits = c(0, 100), name = "Concordance (%)") +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_minimal(base_size = 19) +
    theme(plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
          plot.subtitle = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 18),
          panel.grid = element_blank(), legend.title = element_text(size = 16), legend.text = element_text(size = 15),
          legend.key.height = unit(1.3, "cm"), plot.margin = margin(5, 5, 5, 5))
}

make_dotplot <- function(long_df, title, row_order = family_order) {
  label_df <- long_df %>% distinct(group, n_cells) %>% mutate(group_chr = as.character(group)) %>%
    filter(group_chr %in% row_order) %>% mutate(group_chr = factor(group_chr, levels = row_order)) %>%
    arrange(group_chr) %>% mutate(group_label = paste0(as.character(group_chr), " (n=", format(round(n_cells), big.mark = ","), ")"))
  plot_df <- long_df %>% mutate(group_chr = as.character(group)) %>%
    left_join(label_df %>% select(group_chr, group_label), by = "group_chr") %>%
    filter(!is.na(group_label)) %>% mutate(group_label = factor(group_label, levels = rev(label_df$group_label)))
  ggplot(plot_df, aes(x = feature, y = group_label)) +
    geom_point(aes(size = scaled, fill = scaled), shape = 21, color = "grey25", stroke = 0.2) +
    scale_fill_viridis_c(option = "magma", begin = 0.15, end = 0.85, limits = c(0, 1), name = "Average Score") +
    scale_size(range = c(2.8, 12.0), limits = c(0, 1), guide = "none") +
    labs(x = NULL, y = NULL, title = title) + theme_bw(base_size = 18) +
    theme(plot.title = element_text(face = "bold", size = 21, hjust = 0),
          axis.text.x = element_text(angle = 42, hjust = 1, vjust = 1, size = 15), axis.text.y = element_text(size = 16),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.25), panel.grid.minor = element_blank(),
          legend.title = element_text(size = 15), legend.text = element_text(size = 14), legend.position = "right")
}

make_scatter_panel <- function(m1_df, wmb_df, summary_df) {
  m1_long <- matrix_to_long(m1_df, group_order = family_order) %>% select(group, feature, m1_scaled = scaled)
  wmb_long <- matrix_to_long(wmb_df, group_order = family_order) %>% select(group, feature, wmb_scaled = scaled)
  scatter_df <- inner_join(m1_long, wmb_long, by = c("group", "feature"))
  s <- summary_df %>% filter(comparison == "display_scaled01_flattened") %>% slice(1)
  label <- sprintf("Spearman rho = %.2f\nn = %d", s$spearman_rho, s$n)
  ggplot(scatter_df, aes(x = m1_scaled, y = wmb_scaled, color = group)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.55) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.75) +
    geom_point(size = 3.8, alpha = 0.9) +
    scale_color_manual(values = family_colors, drop = FALSE) +
    annotate("label", x = 0.03, y = 0.97, label = label, hjust = 0, vjust = 1, size = 5.1) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = TRUE) +
    labs(x = "M1 measured ephys (scaled 0-1)", y = "WMB predicted ephys (scaled 0-1)",
         color = "Subclass", title = "M1 vs WMB family-profile correlation") +
    theme_bw(base_size = 18) +
    theme(plot.title = element_text(face = "bold", size = 21, hjust = 0), axis.text = element_text(size = 16),
          axis.title = element_text(size = 17), legend.title = element_text(size = 16), legend.text = element_text(size = 15),
          panel.grid.minor = element_blank())
}

m1_path <- file.path(clean_figures_dir, "DotPlot_M1_patchseq_real_ephys_CTX_M1_family_matched.medians.csv")
m1_n_path <- file.path(clean_figures_dir, "DotPlot_M1_patchseq_real_ephys_CTX_M1_family_matched.ncells.csv")
m1_df <- read_matrix(m1_path)
m1_n <- read_csv(m1_n_path, show_col_types = FALSE)

wmb_celltype_filtered <- load_filtered_wmb(c("wmb_isocortex_1_m1_celltype_exons_preds_exact", "wmb_isocortex_2_m1_celltype_exons_preds_exact"))
wmb_ephys_filtered <- load_filtered_wmb(c("wmb_isocortex_1_m1_patchseq_predictions", "wmb_isocortex_2_m1_patchseq_predictions"))

heat_df_all_cells <- make_concordance_summary(wmb_celltype_filtered, cutoffs = c(0.0))

wmb_summary <- wmb_ephys_filtered %>% group_by(group = true_family) %>%
  summarise(across(all_of(feature_order), ~ median(.x, na.rm = TRUE)), n_cells = n(), .groups = "drop") %>%
  mutate(group = factor(as.character(group), levels = family_order)) %>% arrange(group)
wmb_df <- wmb_summary %>% select(group, all_of(feature_order))
wmb_n <- wmb_summary %>% select(group, n_cells)

m1_long_for_corr <- matrix_to_long(m1_df, group_order = family_order) %>% select(group, feature, m1_scaled = scaled)
wmb_long_for_corr <- matrix_to_long(wmb_df, group_order = family_order) %>% select(group, feature, wmb_scaled = scaled)
corr_df <- inner_join(m1_long_for_corr, wmb_long_for_corr, by = c("group", "feature"))
pear <- suppressWarnings(cor.test(corr_df$m1_scaled, corr_df$wmb_scaled, method = "pearson"))
spear <- suppressWarnings(cor.test(corr_df$m1_scaled, corr_df$wmb_scaled, method = "spearman"))
summary_df <- tibble(comparison = "display_scaled01_flattened", n = nrow(corr_df),
                     pearson_r = unname(pear$estimate), spearman_rho = unname(spear$estimate))

p_concordance_all <- make_concordance_heatmap(heat_df_all_cells,
  title = "WMB concordance",
  subtitle = NULL)
p_scatter <- make_scatter_panel(m1_df, wmb_df, summary_df)
p_m1 <- make_dotplot(matrix_to_long(m1_df, group_order = family_order, n_cells_df = m1_n), "M1 Patch-seq measured ephys")
p_wmb <- make_dotplot(matrix_to_long(wmb_df, group_order = family_order, n_cells_df = wmb_n), "WMB predicted ephys")

final_plot_no_p05 <- (p_concordance_all + p_scatter) / (p_m1 + p_wmb) +
  plot_layout(widths = c(1.0, 1.15), heights = c(1.05, 1.05)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 26))

out_prefix <- file.path(figures_dir, "Figure5_WMB_concordance_ephys")
ggsave(paste0(out_prefix, ".png"), final_plot_no_p05, width = 20.8, height = 14.8, dpi = 300, bg = "white")
ggsave(paste0(out_prefix, ".pdf"), final_plot_no_p05, width = 20.8, height = 14.8, bg = "white")
message("[saved] ", paste0(out_prefix, ".png"))
message(sprintf("Pearson r = %.3f  Spearman rho = %.3f  n = %d", summary_df$pearson_r, summary_df$spearman_rho, summary_df$n))
