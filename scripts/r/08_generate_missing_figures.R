#!/usr/bin/env Rscript
# 08_generate_missing_figures.R -- Generate missing figures for 2018 pilot
# Figures 3 (rarefaction), 8-10 (differential abundance)

suppressMessages({
  library(phyloseq); library(vegan); library(ape)
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
})

PROJECT <- path.expand("~/pilot2018")
EXP     <- file.path(PROJECT, "results/r/exported")
OUT     <- file.path(PROJECT, "results/figures_generated")
dir.create(OUT, showWarnings = FALSE)
set.seed(1)

# Load exported QIIME2 artifacts
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

tree_path <- file.path(EXP, "tree.nwk")
tree <- read_tree(tree_path)

ps <- phyloseq(otu, TAX, SAM, tree)

# FIGURE 3: RAREFACTION CURVES
cat("Generating Figure 3: Rarefaction Curves...\n")

min_depth <- min(sample_sums(ps))
cat(sprintf("Rarefying to depth = %d reads\n", min_depth))

rarefaction_list <- lapply(1:5, function(i) {
  ps_rar <- rarefy_even_depth(ps, sample.size = min_depth, 
                              rngseed = i, replace = FALSE, verbose = FALSE)
  alpha_i <- estimate_richness(ps_rar, measures = "Observed")
  alpha_i$sample_name <- rownames(alpha_i)
  alpha_i$iteration <- i
  return(alpha_i)
})

rarefaction_data <- do.call(rbind, rarefaction_list)
rarefaction_data$depth <- min_depth
rarefaction_data <- merge(rarefaction_data, 
                          meta[c("environment", "treatment")], 
                          by = "row.names")
rownames(rarefaction_data) <- rarefaction_data$Row.names
rarefaction_data$Row.names <- NULL

p_rarefaction <- ggplot(rarefaction_data, aes(x = depth, y = Observed, 
                                               group = sample_name, 
                                               color = environment)) +
  geom_line(alpha = 0.6, linewidth = 1) +
  facet_wrap(~treatment) +
  labs(x = "Rarefaction Depth (reads)", y = "Observed Richness", 
       title = "Rarefaction Curves") +
  theme_minimal()

ggsave(file.path(OUT, "figure_3_rarefaction.pdf"), p_rarefaction, 
       width = 10, height = 6)
cat("Figure 3 saved\n")

# FIGURES 8-10: DIFFERENTIAL ABUNDANCE (Gut only)
cat("Generating Figures 8-10: Differential Abundance Analysis...\n")

ps_gut <- subset_samples(ps, environment == "gut")
ps_gut_rel <- transform_sample_counts(ps_gut, function(x) x / sum(x))

# Get abundance data
abund_df <- as.data.frame(otu_table(ps_gut_rel))
abund_df$ASV <- rownames(abund_df)

abund_long <- reshape2::melt(abund_df, id.vars = "ASV", 
                              variable.name = "sample_name",
                              value.name = "rel_abundance")

# Add taxonomy
tax_df <- as.data.frame(tax_split)
tax_df$ASV <- rownames(tax_df)
abund_long <- merge(abund_long, tax_df, by = "ASV", all.x = TRUE)

# Add metadata
meta$sample_name <- rownames(meta)
abund_long <- merge(abund_long, meta[c("environment", "treatment", "sample_name")], 
                    by = "sample_name", all.x = TRUE)

# Calculate mean abundance by ASV and treatment
asv_means <- aggregate(rel_abundance ~ ASV + treatment + Genus, 
                       data = abund_long, 
                       FUN = function(x) c(mean = mean(x), 
                                          sd = sd(x), 
                                          n = length(x)))
asv_means <- cbind(asv_means[,1:3], asv_means[,4])
colnames(asv_means)[4:6] <- c("mean_abundance", "sd_abundance", "n")
asv_means$se <- asv_means$sd_abundance / sqrt(asv_means$n)

# Fold change: Roundup vs Control
control_means <- asv_means[asv_means$treatment == "control", c("ASV", "mean_abundance")]
colnames(control_means) <- c("ASV", "Control")
roundup_means <- asv_means[asv_means$treatment == "roundup", c("ASV", "mean_abundance")]
colnames(roundup_means) <- c("ASV", "Roundup")

fc_data <- merge(control_means, roundup_means, by = "ASV", all = TRUE)
fc_data[is.na(fc_data)] <- 1e-6
fc_data$log2fc <- log2(fc_data$Roundup / fc_data$Control)

# Add genus info
fc_data <- merge(fc_data, 
                 unique(abund_long[c("ASV", "Genus")]), 
                 by = "ASV", all.x = TRUE)
fc_data <- fc_data[order(abs(fc_data$log2fc), decreasing = TRUE), ]

# FIGURE 8: TOP RESPONDERS
cat("Generating Figure 8: Top Responders...\n")

top_up <- fc_data[fc_data$log2fc > 0, "ASV"][1:6]
top_down <- fc_data[fc_data$log2fc < 0, "ASV"][1:6]
top_asv <- c(na.omit(top_up), na.omit(top_down))

responder_data <- abund_long[abund_long$ASV %in% top_asv, ]
responder_summary <- aggregate(rel_abundance ~ ASV + treatment + Genus, 
                               data = responder_data, 
                               FUN = function(x) c(mean = mean(x), 
                                                  se = sd(x) / sqrt(length(x))))
responder_summary <- cbind(responder_summary[,1:3], responder_summary[,4])
colnames(responder_summary)[4:5] <- c("mean", "se")
responder_summary$label <- ifelse(!is.na(responder_summary$Genus), 
                                  responder_summary$Genus, 
                                  responder_summary$ASV)

p_responders <- ggplot(responder_summary, aes(x = treatment, y = mean, fill = treatment)) +
  geom_col() +
  geom_errorbar(aes(ymin = pmax(0, mean - se), ymax = mean + se), 
                width = 0.2) +
  facet_wrap(~label, scales = "free_y", nrow = 4) +
  labs(x = "Treatment", y = "Relative Abundance", 
       title = "Top Dose-Responsive Taxa (Gut)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT, "figure_8_top_responders.pdf"), p_responders, 
       width = 12, height = 10)
cat("Figure 8 saved\n")

# FIGURE 9: TOP 32 GENERA BOX PLOTS
cat("Generating Figure 9: Top 32 Genera...\n")

top_genera <- aggregate(rel_abundance ~ Genus, 
                        data = abund_long, 
                        FUN = mean)
top_genera <- top_genera[order(top_genera$rel_abundance, decreasing = TRUE), ]
top_32 <- head(top_genera$Genus, 32)

genera_data <- abund_long[abund_long$Genus %in% top_32, ]
genera_data$Genus <- factor(genera_data$Genus, levels = top_32)

p_genera <- ggplot(genera_data, aes(x = treatment, y = rel_abundance, fill = treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  facet_wrap(~Genus, scales = "free_y", nrow = 8) +
  labs(x = "Treatment", y = "Relative Abundance", 
       title = "Top 32 Genera by Mean Abundance (Gut)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(OUT, "figure_9_top_32_genera.pdf"), p_genera, 
       width = 14, height = 16)
cat("Figure 9 saved\n")

# FIGURE 10: VOLCANO-STYLE PLOT
cat("Generating Figure 10: Volcano Plot...\n")

mean_abund <- aggregate(rel_abundance ~ ASV, 
                        data = abund_long, 
                        FUN = mean)
colnames(mean_abund)[2] <- "mean_abund"

volcano_data <- merge(fc_data, mean_abund, by = "ASV")
volcano_data$class <- ifelse(abs(volcano_data$log2fc) > 1,
                             ifelse(volcano_data$log2fc > 0, "Enriched", "Depleted"),
                             "NS")

p_volcano <- ggplot(volcano_data, aes(x = log2fc, y = mean_abund)) +
  geom_point(aes(color = class), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Enriched" = "blue", "Depleted" = "red", "NS" = "gray")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray", alpha = 0.3) +
  labs(x = "log2(Fold Change: Roundup/Control)", 
       y = "Mean Relative Abundance",
       title = "Differential Abundance: Control vs Roundup (Gut)",
       color = "Classification") +
  theme_minimal()

ggsave(file.path(OUT, "figure_10_volcano_plot.pdf"), p_volcano, 
       width = 10, height = 7)
cat("Figure 10 saved\n")

cat("\nAll figures generated successfully!\n")
cat(sprintf("Saved to: %s\n", OUT))
