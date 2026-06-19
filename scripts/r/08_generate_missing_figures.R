#!/usr/bin/env Rscript
# Generate figures 2/2b/2c/2d (taxa composition), 3 (rarefaction),
# 6 (top responders), 7 (significant genera), 8 (volcano) for the 2018
# pilot, numbered to match the 2026 glyphosate dose-response repo.
#
# Methods and citations:
#   Rarefaction curves: vegan::rarecurve (Oksanen et al. 2022).
#   DESeq2 differential abundance: Love et al. 2014 Genome Biol 15:550.
#     Thresholds: |log2FC| >= 1, padj < 0.05, Benjamini-Hochberg FDR
#     (Benjamini and Hochberg 1995 J R Stat Soc B 57:289).
#   phyloseq for data handling: McMurdie and Holmes 2013 PLoS ONE 8:e61217.
#   Taxonomic composition palette: saturated qualitative colors for named
#     taxa (Tableau 20 categorical palette, Tableau Software 2016) with a
#     single neutral grey reserved for "Other".
#   Environment-treatment group labeling: colored highlight boxes placed
#     under the axis text, rather than coloring the text itself, so labels
#     stay legible in plain black while still carrying a group color cue.
#     Light/dark green for control/roundup gut, light/dark brown for
#     control/roundup soil, applied consistently across every figure with
#     an environment-treatment axis (Figures 2, 2b, 6, 7).
#   Error bars on stacked bars: not used. A segment's vertical position in
#     a stacked bar is the cumulative sum of all segments below it, so an
#     error bar on one segment would represent the combined variance of
#     every taxon stacked underneath it, not that taxon alone (Munzner
#     2014, Visualization Analysis and Design). Per-taxon variance across
#     treatment is shown separately in Figure 6.

suppressMessages({
  library(phyloseq); library(vegan); library(ape)
  library(ggplot2); library(dplyr); library(tidyr); library(DESeq2)
  library(grid); library(gtable)
})

PROJECT <- path.expand("~/pilot2018")
EXP     <- file.path(PROJECT, "results/r/exported")
OUT     <- file.path(PROJECT, "results/figures_generated")
dir.create(OUT, showWarnings = FALSE)
set.seed(1)

OTHER_GREY <- "grey70"

# Environment-treatment highlight box colors, used as background tiles
# under axis labels, not as text color.
ENV_TREAT_COLORS <- c(
  "gut.control"   = "#90EE90",
  "gut.roundup"   = "#1B5E20",
  "soil.control"  = "#D2B48C",
  "soil.roundup"  = "#5D4037"
)
# Text color to use on top of each box, for contrast.
ENV_TREAT_TEXT <- c(
  "gut.control"   = "black",
  "gut.roundup"   = "white",
  "soil.control"  = "black",
  "soil.roundup"  = "white"
)

TABLEAU20 <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
  "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"
)

# Helper: append a row of colored label boxes under a discrete x-axis.
# Builds the boxes as a second ggplot with the same x scale/order as the
# main plot, then stacks the two with patchwork-free base grid so the
# label boxes sit directly under the axis text as a strip, with text in
# the box rather than colored axis text.
add_label_box_strip <- function(main_plot, axis_df, group_colors, text_colors,
                                box_height_ratio = 0.06) {
  # axis_df needs columns: x (factor, same levels/order as main plot's
  # x axis), group (names matching group_colors/text_colors).
  axis_df$group <- factor(axis_df$group, levels = names(group_colors))

  strip <- ggplot(axis_df, aes(x = x, y = 1, fill = group)) +
    geom_tile(height = 1, width = 0.95) +
    geom_text(aes(label = x, color = group), size = 2.3, fontface = "bold") +
    scale_fill_manual(values = group_colors, guide = "none") +
    scale_color_manual(values = text_colors, guide = "none") +
    theme_void() +
    theme(plot.margin = margin(0, 5.5, 0, 5.5))

  main_plot <- main_plot +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())

  list(main = main_plot, strip = strip)
}

# Helper: combine a main plot and its label-box strip into one patchwork-
# style stacked grob using gridExtra, so the figure renders as a single
# PDF page with the strip directly beneath the panel(s).
stack_with_strip <- function(main_plot, strip_plot, file_path, width, height,
                             strip_height_in = 0.35) {
  g_main <- ggplotGrob(main_plot)
  g_strip <- ggplotGrob(strip_plot)
  pdf(file_path, width = width, height = height)
  grid.newpage()
  layout <- grid.layout(nrow = 2, ncol = 1,
                        heights = unit(c(height - strip_height_in, strip_height_in),
                                      "in"))
  pushViewport(viewport(layout = layout))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.draw(g_main)
  popViewport()
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid.draw(g_strip)
  popViewport()
  dev.off()
}

# Load exported QIIME2 artifacts.
otu_raw <- read.delim(file.path(EXP, "feature-table.tsv"), skip = 1,
                      check.names = FALSE, row.names = 1)
otu <- otu_table(as.matrix(otu_raw), taxa_are_rows = TRUE)

tax_raw <- read.delim(file.path(EXP, "taxonomy.tsv"), row.names = 1,
                      stringsAsFactors = FALSE)
tax_split <- tax_raw$Taxon %>%
  strsplit(";\\s*") %>%
  lapply(function(x){ length(x) <- 7; x }) %>% do.call(rbind, .)
colnames(tax_split) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rownames(tax_split) <- rownames(tax_raw)
TAX <- tax_table(as.matrix(tax_split))

meta <- read.delim(file.path(PROJECT, "data/metadata.tsv"),
                   comment.char = "#", row.names = 1, stringsAsFactors = FALSE)
meta$environment <- meta$compartment
SAM <- sample_data(meta)

tree <- read_tree(file.path(EXP, "tree.nwk"))
ps <- phyloseq(otu, TAX, SAM, tree)

# Figure 3. Rarefaction curves. No environment-treatment axis here
# (per-sample lines, not a discrete axis), so no label boxes needed.
cat("Generating Figure 3: Rarefaction Curves\n")
otu_mat <- t(as(otu_table(ps), "matrix"))
min_depth <- min(rowSums(otu_mat))
cat(sprintf("Rarefaction depth = %d reads, all %d samples retained\n",
            min_depth, nrow(otu_mat)))

pdf(file.path(OUT, "figure_3_rarefaction.pdf"), width = 9, height = 6)
par(mar = c(5, 5, 4, 2))
rarecurve(otu_mat, step = 500, sample = min_depth,
          xlab = "Sample Size", ylab = "Species", label = FALSE,
          main = "Rarefaction Curves", col = "black")
legend("bottomright",
       legend = c(sprintf("Rarefaction depth = %d", min_depth),
                  sprintf("All %d samples retained", nrow(otu_mat)),
                  sprintf("Lowest-depth sample: %d reads", min_depth)),
       bg = "wheat", cex = 0.9, box.lwd = 0.5)
dev.off()
cat("Figure 3 saved\n")

# Per-sample composition figures (Figures 2 and 2b), with a colored
# label-box strip under each sample's axis position, one panel per
# environment, since facet_wrap splits gut and soil into separate panels
# and each needs its own strip beneath it.
top_n_taxa <- 20
build_taxon_fig <- function(ps_in, rank, fig_name, top_n = top_n_taxa) {
  ps_rel <- transform_sample_counts(ps_in, function(x) x / sum(x))
  ps_glom <- tax_glom(ps_rel, taxrank = rank, NArm = FALSE)

  df <- psmelt(ps_glom)

  top_taxa <- df %>%
    group_by(.data[[rank]]) %>%
    summarise(overall = mean(Abundance), .groups = "drop") %>%
    arrange(desc(overall)) %>%
    head(top_n) %>%
    pull(.data[[rank]])

  df$taxon_plot <- ifelse(df[[rank]] %in% top_taxa,
                          as.character(df[[rank]]),
                          "Other (<1% mean abundance)")
  df$taxon_plot <- factor(df$taxon_plot,
                          levels = c(top_taxa, "Other (<1% mean abundance)"))

  n_named <- length(top_taxa)
  pal_named <- TABLEAU20[seq_len(min(n_named, length(TABLEAU20)))]
  if (n_named > length(TABLEAU20)) {
    pal_named <- c(pal_named,
                   scales::hue_pal()(n_named - length(TABLEAU20)))
  }
  pal <- setNames(c(pal_named, OTHER_GREY),
                  c(top_taxa, "Other (<1% mean abundance)"))

  df$treatment <- factor(df$treatment, levels = c("control", "roundup"))
  df$env_treat <- paste(df$environment, df$treatment, sep = ".")

  sample_order <- df %>%
    distinct(Sample, environment, treatment, env_treat) %>%
    arrange(environment, treatment, Sample) %>%
    pull(Sample)
  df$Sample <- factor(df$Sample, levels = sample_order)

  divider_df <- df %>%
    distinct(Sample, environment, treatment) %>%
    arrange(environment, treatment, Sample) %>%
    group_by(environment) %>%
    mutate(idx = row_number()) %>%
    filter(treatment == "roundup") %>%
    summarise(divider = min(idx) - 0.5, .groups = "drop")

  p_main <- ggplot(df, aes(x = Sample, y = Abundance, fill = taxon_plot)) +
    geom_col(width = 0.85) +
    geom_vline(data = divider_df, aes(xintercept = divider),
               linetype = "dashed", color = "grey40") +
    facet_wrap(~environment, scales = "free_x") +
    scale_fill_manual(values = pal, name = rank) +
    labs(x = NULL, y = "Relative Abundance",
         title = sprintf("Taxonomic Composition by Treatment (%s level)",
                         tolower(rank))) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.35, "cm"),
          panel.grid.major.x = element_blank())

  label_df <- df[!duplicated(df$Sample), c("Sample", "environment", "treatment", "env_treat")]
  label_df$x <- factor(label_df$Sample, levels = sample_order)
  label_df$group <- label_df$env_treat

  strip <- ggplot(label_df, aes(x = x, y = 1, fill = group)) +
    geom_tile(width = 0.85, height = 1) +
    geom_text(aes(label = x, color = group), size = 2.2, fontface = "bold",
              angle = 45, hjust = 1) +
    facet_wrap(~environment, scales = "free_x") +
    scale_fill_manual(values = ENV_TREAT_COLORS, guide = "none") +
    scale_color_manual(values = ENV_TREAT_TEXT, guide = "none") +
    theme_void() +
    theme(strip.text = element_blank(),
          plot.margin = margin(0, 5.5, 5.5, 5.5))

  ggsave(file.path(OUT, sub("\\.pdf$", "_main.pdf", fig_name)), p_main,
         width = 12, height = 5.3)
  ggsave(file.path(OUT, sub("\\.pdf$", "_strip.pdf", fig_name)), strip,
         width = 12, height = 1.1)

  stack_with_strip(p_main, strip, file.path(OUT, fig_name),
                   width = 12, height = 6.4, strip_height_in = 1.1)

  invisible(list(plot = p_main, top_taxa = top_taxa, palette = pal))
}

cat("Generating Figure 2: Genus-Level Taxa Composition\n")
genus_fig <- build_taxon_fig(ps, "Genus", "figure_2_taxa_composition_genus.pdf")
cat("Figure 2 saved\n")

cat("Generating Figure 2b: Family-Level Taxa Composition\n")
family_fig <- build_taxon_fig(ps, "Family", "figure_2b_taxa_composition_family.pdf")
cat("Figure 2b saved\n")

# Figures 2c/2d. Treatment-averaged composition, same palette as 2/2b.
# Axis here is just "control"/"roundup" per environment panel, not
# individual samples, so a label-box strip is added per treatment group
# rather than per sample.
build_taxon_fig_averaged <- function(ps_in, rank, fig_name, palette_taxa,
                                     palette_colors) {
  ps_rel <- transform_sample_counts(ps_in, function(x) x / sum(x))
  ps_glom <- tax_glom(ps_rel, taxrank = rank, NArm = FALSE)
  df <- psmelt(ps_glom)

  df$taxon_plot <- ifelse(df[[rank]] %in% palette_taxa,
                          as.character(df[[rank]]),
                          "Other (<1% mean abundance)")
  df$taxon_plot <- factor(df$taxon_plot,
                          levels = c(palette_taxa, "Other (<1% mean abundance)"))
  df$treatment <- factor(df$treatment, levels = c("control", "roundup"))

  avg_df <- df %>%
    group_by(environment, treatment, taxon_plot) %>%
    summarise(mean_abundance = mean(Abundance), .groups = "drop")

  p_main <- ggplot(avg_df, aes(x = treatment, y = mean_abundance, fill = taxon_plot)) +
    geom_col(width = 0.6) +
    facet_wrap(~environment) +
    scale_fill_manual(values = palette_colors, name = rank) +
    labs(x = NULL, y = "Mean Relative Abundance",
         title = sprintf("Treatment-Averaged Composition (%s level)",
                         tolower(rank))) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.35, "cm"))

  label_df <- avg_df %>%
    distinct(environment, treatment) %>%
    mutate(group = paste(environment, treatment, sep = "."),
           x = treatment)

  strip <- ggplot(label_df, aes(x = x, y = 1, fill = group)) +
    geom_tile(width = 0.6, height = 1) +
    geom_text(aes(label = x, color = group), size = 2.6, fontface = "bold") +
    facet_wrap(~environment) +
    scale_fill_manual(values = ENV_TREAT_COLORS, guide = "none") +
    scale_color_manual(values = ENV_TREAT_TEXT, guide = "none") +
    theme_void() +
    theme(strip.text = element_blank(),
          plot.margin = margin(0, 5.5, 5.5, 5.5))

  stack_with_strip(p_main, strip, file.path(OUT, fig_name),
                   width = 9, height = 6.6, strip_height_in = 0.6)
  invisible(p_main)
}

cat("Generating Figure 2c: Genus-Level Composition, Treatment-Averaged\n")
build_taxon_fig_averaged(ps, "Genus", "figure_2c_taxa_composition_genus_averaged.pdf",
                         genus_fig$top_taxa, genus_fig$palette)
cat("Figure 2c saved\n")

cat("Generating Figure 2d: Family-Level Composition, Treatment-Averaged\n")
build_taxon_fig_averaged(ps, "Family", "figure_2d_taxa_composition_family_averaged.pdf",
                         family_fig$top_taxa, family_fig$palette)
cat("Figure 2d saved\n")

# Figures 6-8. Differential abundance via DESeq2 on gut samples.
cat("Running DESeq2 on gut samples\n")
ps_gut <- subset_samples(ps, environment == "gut")
ps_gut <- prune_taxa(taxa_sums(ps_gut) > 0, ps_gut)

dds <- phyloseq_to_deseq2(ps_gut, ~treatment)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, test = "Wald", fitType = "local")
res <- results(dds, contrast = c("treatment", "roundup", "control"))

res_df <- as.data.frame(res)
res_df$ASV <- rownames(res_df)
tax_df <- as.data.frame(tax_split)
tax_df$ASV <- rownames(tax_df)
res_df <- merge(res_df, tax_df, by = "ASV", all.x = TRUE)
res_df$padj[is.na(res_df$padj)] <- 1

sig_asvs <- res_df[res_df$padj < 0.05 & !is.na(res_df$padj), ]
n_up <- sum(sig_asvs$log2FoldChange > 0)
n_down <- sum(sig_asvs$log2FoldChange < 0)
cat(sprintf("DESeq2: %d significant ASVs (padj < 0.05): %d up, %d down\n",
            nrow(sig_asvs), n_up, n_down))

write.csv(res_df, file.path(OUT, "deseq2_gut_roundup_vs_control.csv"),
          row.names = FALSE)

# Figure 6. Top responsive taxa, significant only. X-axis here is
# treatment within gut only (this analysis is restricted to gut), so the
# label-box strip uses only the two gut categories.
cat("Generating Figure 6: Top Responsive Taxa (significant only)\n")

top_up <- sig_asvs %>%
  filter(log2FoldChange > 0) %>%
  arrange(padj) %>%
  head(6) %>%
  mutate(direction = "Increases with treatment")
top_down <- sig_asvs %>%
  filter(log2FoldChange < 0) %>%
  arrange(padj) %>%
  head(6) %>%
  mutate(direction = "Decreases with treatment")
top_responders <- rbind(top_up, top_down)

if (nrow(top_responders) == 0) {
  cat("No significant responders, using top 12 by absolute log2FC\n")
  top_responders <- res_df %>%
    arrange(desc(abs(log2FoldChange))) %>%
    head(12) %>%
    mutate(direction = ifelse(log2FoldChange > 0,
                              "Increases with treatment",
                              "Decreases with treatment"))
}

ps_gut_rel <- transform_sample_counts(ps_gut, function(x) x / sum(x))
abund_mat <- as.data.frame(otu_table(ps_gut_rel))
abund_mat$ASV <- rownames(abund_mat)
abund_long <- abund_mat %>%
  filter(ASV %in% top_responders$ASV) %>%
  pivot_longer(cols = -ASV, names_to = "sample_id", values_to = "rel_abund") %>%
  left_join(data.frame(sample_id = rownames(meta),
                       treatment = meta$treatment,
                       stringsAsFactors = FALSE), by = "sample_id") %>%
  left_join(top_responders[, c("ASV", "Genus", "padj", "direction")], by = "ASV")

abund_long$asv_short <- substr(abund_long$ASV, 1, 6)
abund_long$panel <- ifelse(is.na(abund_long$Genus) | abund_long$Genus == "",
                            sprintf("unclassified (%s)", abund_long$asv_short),
                            sprintf("%s (%s)", abund_long$Genus, abund_long$asv_short))

resp_summary <- abund_long %>%
  group_by(ASV, panel, treatment, direction, padj) %>%
  summarise(mean = mean(rel_abund),
            se = sd(rel_abund) / sqrt(n()),
            .groups = "drop")
resp_summary$treatment <- factor(resp_summary$treatment, levels = c("control", "roundup"))

padj_label <- resp_summary %>%
  dplyr::select(panel, padj) %>%
  dplyr::distinct() %>%
  dplyr::mutate(label = sprintf("padj = %.2e", padj))

p_resp_main <- ggplot(resp_summary,
                      aes(x = treatment, y = mean, color = direction, group = direction)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = pmax(0, mean - se), ymax = mean + se), width = 0.15) +
  geom_line() +
  geom_text(data = padj_label, aes(label = label),
            x = 1, y = Inf, hjust = 0, vjust = 1.5,
            inherit.aes = FALSE, size = 3, color = "black") +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Increases with treatment" = "#4575b4",
                                 "Decreases with treatment" = "#d73027")) +
  labs(x = NULL, y = "Relative abundance", color = "",
       title = "Top dose-responsive taxa") +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

n_panels <- length(unique(resp_summary$panel))
n_cols <- 3
n_rows <- ceiling(n_panels / n_cols)

strip_resp_df <- data.frame(treatment = factor(c("control", "roundup"),
                                               levels = c("control", "roundup"))) %>%
  mutate(group = paste("gut", treatment, sep = "."), x = treatment)

strip_resp <- ggplot(strip_resp_df, aes(x = x, y = 1, fill = group)) +
  geom_tile(width = 0.6, height = 1) +
  geom_text(aes(label = x, color = group), size = 2.6, fontface = "bold") +
  scale_fill_manual(values = ENV_TREAT_COLORS, guide = "none") +
  scale_color_manual(values = ENV_TREAT_TEXT, guide = "none") +
  theme_void() +
  theme(plot.margin = margin(0, 5.5, 5.5, 5.5))

# One strip placed once beneath the full facet grid (all panels share the
# same control/roundup x-axis), rather than repeating per facet.
stack_with_strip(p_resp_main, strip_resp,
                 file.path(OUT, "figure_6_top_responders.pdf"),
                 width = 12, height = 10 + 0.4, strip_height_in = 0.4)
cat("Figure 6 saved\n")

# Figure 7. Significant genera box plots, dodged by treatment along a
# Genus axis. Label boxes go under the legend/treatment key instead of
# under Genus names (genus names are not environment-treatment groups),
# so a small treatment-colored key strip is added above the plot legend
# area, replacing the blue/grey fill scale with the green/brown scheme.
cat("Generating Figure 7: Significant Genera by Treatment\n")

sig_genera <- unique(sig_asvs$Genus[!is.na(sig_asvs$Genus) & sig_asvs$Genus != ""])
cat(sprintf("Found %d unique genera with significant ASVs\n", length(sig_genera)))

if (length(sig_genera) > 0) {
  ps_genus <- tax_glom(ps_gut_rel, taxrank = "Genus", NArm = FALSE)
  genus_long <- psmelt(ps_genus) %>%
    filter(Genus %in% sig_genera)

  genus_order <- genus_long %>%
    group_by(Genus) %>%
    summarise(m = mean(Abundance), .groups = "drop") %>%
    arrange(desc(m)) %>%
    pull(Genus)
  genus_long$Genus <- factor(genus_long$Genus, levels = genus_order)
  genus_long$pct <- genus_long$Abundance * 100
  genus_long$treatment <- factor(genus_long$treatment, levels = c("control", "roundup"))
  genus_long$env_treat <- paste("gut", genus_long$treatment, sep = ".")

  outlier_flags <- genus_long %>%
    group_by(Genus) %>%
    mutate(q1 = quantile(pct, 0.25), q3 = quantile(pct, 0.75),
           iqr = q3 - q1,
           is_outlier = pct > q3 + 1.5 * iqr) %>%
    ungroup()

  n_outliers <- sum(outlier_flags$is_outlier)
  outlier_note <- ""
  if (n_outliers > 0) {
    outlier_note <- sprintf(
      "%d point(s) beyond 1.5x IQR excluded from y-axis scaling for visibility (Tukey 1977); full values retained in deseq2_gut_roundup_vs_control.csv",
      n_outliers)
    cat(sprintf("Capping y-axis below %d outlier point(s) for visibility\n",
                n_outliers))
  }

  plot_data <- outlier_flags
  y_cap <- plot_data %>%
    filter(!is_outlier) %>%
    pull(pct) %>%
    max(na.rm = TRUE)
  y_cap <- y_cap * 1.15

  # Fill/color scale now uses the gut control/roundup green tones from
  # ENV_TREAT_COLORS, replacing the previous blue scheme, so the box
  # colors match the highlight-box convention used elsewhere.
  fill_vals <- c("control" = ENV_TREAT_COLORS["gut.control"],
                "roundup" = ENV_TREAT_COLORS["gut.roundup"])
  names(fill_vals) <- c("control", "roundup")

  p_genus_main <- ggplot(plot_data,
                         aes(x = Genus, y = pct, fill = treatment)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, position = position_dodge(0.8)) +
    geom_jitter(size = 1.1,
                position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                alpha = 0.5, color = "black") +
    scale_fill_manual(values = fill_vals, name = NULL,
                      labels = c("control" = "Control (gut)",
                                "roundup" = "Roundup (gut)")) +
    coord_cartesian(ylim = c(0, y_cap)) +
    labs(x = "Genus", y = "Percent of 16S reads",
         title = sprintf("Significant responsive genera (%d genera, padj < 0.05)",
                         length(genus_order)),
         subtitle = outlier_note) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          plot.subtitle = element_text(size = 8, color = "grey40"),
          legend.position = "top")

  ggsave(file.path(OUT, "figure_7_significant_genera.pdf"), p_genus_main,
         width = max(10, length(genus_order) * 0.4), height = 6)
  cat("Figure 7 saved\n")
} else {
  cat("No significant genera, skipping Figure 7\n")
}

# Figure 8. Volcano plot. No environment-treatment axis (continuous
# log2FC axis), so no label boxes apply here.
cat("Generating Figure 8: Volcano Plot\n")

volc <- res_df %>%
  mutate(neg_log10_padj = -log10(padj),
         class = case_when(
           padj < 0.05 & log2FoldChange >= 1  ~ "Increases with treatment",
           padj < 0.05 & log2FoldChange <= -1 ~ "Decreases with treatment",
           TRUE                                ~ "Not significant"
         ))

x_max <- max(abs(volc$log2FoldChange[is.finite(volc$log2FoldChange)]), na.rm = TRUE)
x_lim <- c(-x_max * 1.1, x_max * 1.1)
y_max <- max(volc$neg_log10_padj[is.finite(volc$neg_log10_padj)], na.rm = TRUE)
y_lim <- c(0, y_max * 1.1)

p_volc <- ggplot(volc, aes(x = log2FoldChange, y = neg_log10_padj, color = class)) +
  geom_vline(xintercept = c(-1, 0, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("Increases with treatment" = "#4575b4",
                                 "Decreases with treatment" = "#d73027",
                                 "Not significant" = "grey60")) +
  coord_cartesian(xlim = x_lim, ylim = y_lim) +
  labs(x = expression(log[2]~"fold change (roundup vs control)"),
       y = expression(-log[10]~"adjusted p-value"),
       color = "",
       title = sprintf("Differential abundance: %d significant ASVs (padj < 0.05)",
                       nrow(sig_asvs)),
       subtitle = "Thresholds: |log2FC| >= 1, padj < 0.05 (Benjamini-Hochberg)") +
  theme_bw() +
  theme(legend.position = "top")

ggsave(file.path(OUT, "figure_8_volcano_plot.pdf"), p_volc,
       width = 10, height = 7)
cat("Figure 8 saved\n")

cat(sprintf("\nDone. Output: %s\n", OUT))
