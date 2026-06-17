#!/usr/bin/env Rscript
# Alpha and beta diversity visualizations with Pochron's feedback
# Requires: alpha_diversity.csv and Bray-Curtis distance matrix

library(ggplot2)
library(vegan)

# Read alpha diversity
alpha <- read.csv("results/r/alpha_diversity.csv", stringsAsFactors=FALSE)

# Color gradient: control light, roundup dark
colors <- c("control" = "#CCCCCC", "roundup" = "#333333")

# Alpha diversity plots
pdf("results/r/fig_alpha_diversity.pdf", width=12, height=5)

# Shannon diversity by treatment
p1 <- ggplot(alpha, aes(x=treatment, y=Shannon, color=treatment)) +
  geom_boxplot(alpha=0.5) +
  geom_jitter(width=0.2, size=2) +
  scale_color_manual(values=colors) +
  facet_wrap(~compartment) +
  theme_minimal() +
  labs(title="Shannon Diversity Index", x="Treatment", y="Shannon Index")

# Observed ASVs
p2 <- ggplot(alpha, aes(x=treatment, y=Observed, color=treatment)) +
  geom_boxplot(alpha=0.5) +
  geom_jitter(width=0.2, size=2) +
  scale_color_manual(values=colors) +
  facet_wrap(~compartment) +
  theme_minimal() +
  labs(title="Observed ASVs", x="Treatment", y="Number of ASVs")

print(p1)
print(p2)
dev.off()

cat("\nWrote: results/r/fig_alpha_diversity.pdf\n")

# Summary statistics
cat("\n=== Alpha Diversity Summary ===\n")
for (comp in c("gut", "soil")) {
  cat("\n", comp, ":\n")
  subset_data <- alpha[alpha$compartment == comp, ]
  
  for (treat in c("control", "roundup")) {
    subset_treat <- subset_data[subset_data$treatment == treat, ]
    cat(treat, "- Shannon mean:", round(mean(subset_treat$Shannon), 2),
        "sd:", round(sd(subset_treat$Shannon), 2), "\n")
  }
}

