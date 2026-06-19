#!/usr/bin/env Rscript
# Generate figures 3, 8, 9, 10, plus family-level taxa composition for 2018 pilot.
#
# Methods and citations:
#   Rarefaction curves: vegan::rarecurve (Oksanen et al. 2022).
#   DESeq2 differential abundance: Love et al. 2014 Genome Biol 15:550.
#     Standard thresholds: |log2FC| >= 1, padj < 0.05 with Benjamini-Hochberg
#     FDR correction (Benjamini and Hochberg 1995 J R Stat Soc B 57:289).
#   phyloseq for data handling: McMurdie and Holmes 2013 PLoS ONE 8:e61217.

suppressMessages({
  library(phyloseq); library(vegan); library(ape)
  library(ggplot2); library(dplyr); library(tidyr); library(DESeq2)
})

PROJECT <- path.expand("~/pilot2018")
EXP     <- file.path(PROJECT, "results/r/exported")
OUT     <- file.path(PROJECT, "results/figures_generated")
dir.create(OUT, showWarnings = FALSE)
set.seed(1)

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

# Figure 3. Rarefaction curves using vegan::rarecurve.
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

# Figure 2b. Treatment-averaged taxonomic composition at family level.
cat("Generating Figure 2b: Family-Level Taxa Composition\n")
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps_fam <- tax_glom(ps_rel, taxrank = "Family", NArm = FALSE)

fam_df <- psmelt(ps_fam) %>%
  group_by(environment, treatment, Family) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop")

top_fam <- fam_df %>%
  group_by(Family) %>%
  summarise(overall = mean(mean_abundance), .groups = "drop") %>%
  arrange(desc(overall)) %>%
  head(15) %>%
  pull(Family)

fam_df$Family_plot <- ifelse(fam_df$Family %in% top_fam,
                              as.character(fam_df$Family),
                              "Other (<1% mean abundance)")
fam_df$Family_plot <- factor(fam_df$Family_plot,
                              levels = c(top_fam, "Other (<1% mean abundance)"))

p_fam <- ggplot(fam_df, aes(x = treatment, y = mean_abundance, fill = Family_plot)) +
  geom_col() +
  facet_wrap(~environment) +
  labs(x = "Treatment", y = "Mean relative abundance", fill = "Family",
       title = "Treatment-averaged taxonomic composition (family level)") +
  theme_bw() +
  theme(legend.text = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"))

ggsave(file.path(OUT, "figure_2b_taxa_composition_family.pdf"), p_fam,
       width = 10, height = 6)
cat("Figure 2b saved\n")

# Figures 8-10. Differential abundance via DESeq2 on gut samples.
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

# Figure 8. Top responsive taxa.
cat("Generating Figure 8: Top Responsive Taxa (significant only)\n")

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

# Use dplyr::summarise explicitly and pull padj per panel directly.
padj_label <- resp_summary %>%
  dplyr::select(panel, padj) %>%
  dplyr::distinct() %>%
  dplyr::mutate(label = sprintf("padj = %.2e", padj))

p_resp <- ggplot(resp_summary,
                 aes(x = treatment, y = mean, color = direction, group = direction)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = pmax(0, mean - se), ymax = mean + se), width = 0.15) +
  geom_line() +
  geom_text(data = padj_label, aes(label = label),
            x = 1, y = Inf, hjust = 0, vjust = 1.5,
            inherit.aes = FALSE, size = 3,
            color = "black") +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Increases with treatment" = "#4575b4",
                                 "Decreases with treatment" = "#d73027")) +
  labs(x = "Treatment", y = "Relative abundance",
       color = "",
       title = "Top dose-responsive taxa") +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 0))

ggsave(file.path(OUT, "figure_8_top_responders.pdf"), p_resp,
       width = 12, height = 10)
cat("Figure 8 saved\n")

# Figure 9. Top genera among significant responders.
cat("Generating Figure 9: Significant Genera by Treatment\n")

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

  p_genus <- ggplot(genus_long,
                    aes(x = Genus, y = pct, fill = treatment)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(0.8)) +
    geom_jitter(aes(color = treatment), size = 1.2,
                position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                alpha = 0.6) +
    scale_fill_manual(values = c("control" = "#bdd7e7", "roundup" = "#2171b5")) +
    scale_color_manual(values = c("control" = "#08519c", "roundup" = "#08306b")) +
    labs(x = "Genus", y = "Percent of 16S reads",
         title = sprintf("Significant responsive genera (%d genera, padj < 0.05)",
                         length(genus_order))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  ggsave(file.path(OUT, "figure_9_significant_genera.pdf"), p_genus,
         width = max(10, length(genus_order) * 0.4), height = 6)
  cat("Figure 9 saved\n")
} else {
  cat("No significant genera, skipping Figure 9\n")
}

# Figure 10. Volcano plot.
cat("Generating Figure 10: Volcano Plot\n")

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

ggsave(file.path(OUT, "figure_10_volcano_plot.pdf"), p_volc,
       width = 10, height = 7)
cat("Figure 10 saved\n")

cat(sprintf("\nDone. Output: %s\n", OUT))
