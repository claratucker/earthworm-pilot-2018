#!/usr/bin/env Rscript
# 02_denoising_stats_figure.R -- Generate Figure 1: Denoising stats visualization

suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

PROJECT <- path.expand("~/pilot2018")
OUT     <- file.path(PROJECT, "results/figures")
dir.create(OUT, showWarnings = FALSE)
set.seed(1)

# Load denoising stats - skip BOTH header (row 1) and types (row 2)
denoise_file <- file.path(PROJECT, "results/qiime2/denoising-stats-export/stats.tsv")
denoise <- read.delim(denoise_file, skip = 2, stringsAsFactors = FALSE)

# Set proper column names
colnames(denoise) <- c("sample_id", "input", "filtered", "pct_filtered", 
                       "denoised", "non_chimeric", "pct_nonchimeric")

# Load metadata
meta <- read.delim(file.path(PROJECT, "data/metadata.tsv"),
                   comment.char = "#", row.names = 1, stringsAsFactors = FALSE)
meta$sample_id <- rownames(meta)

# Merge
denoise <- merge(denoise, meta[c("sample_id", "compartment", "treatment")], 
                 by = "sample_id", all.x = TRUE)
denoise$environment <- denoise$compartment

# Create figure: retention by sample
p_denoise <- ggplot(denoise, aes(x = reorder(sample_id, -pct_nonchimeric), 
                                   y = pct_nonchimeric, 
                                   fill = environment)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", pct_nonchimeric)), 
            vjust = -0.3, size = 2.5) +
  facet_wrap(~treatment, scales = "free_x") +
  labs(x = "Sample", y = "Percent Retained (%)", 
       title = "Read Retention After DADA2 and Chimera Removal",
       subtitle = "Percentage of input reads retained through full pipeline") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom")

ggsave(file.path(OUT, "figure_1_denoising_stats.pdf"), p_denoise, 
       width = 11, height = 6)

cat("Figure 1 saved to: ", file.path(OUT, "figure_1_denoising_stats.pdf\n"))

# Print summary stats
cat("\n=== DENOISING SUMMARY ===\n")
cat(sprintf("Input reads: min=%d, max=%d, mean=%.0f\n",
            min(denoise$input), max(denoise$input), mean(denoise$input)))
cat(sprintf("Non-chimeric retention: min=%.1f%%, max=%.1f%%, mean=%.1f%%\n",
            min(denoise$pct_nonchimeric), max(denoise$pct_nonchimeric), 
            mean(denoise$pct_nonchimeric)))
