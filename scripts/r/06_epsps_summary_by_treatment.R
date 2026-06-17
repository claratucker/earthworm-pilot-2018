#!/usr/bin/env Rscript
# Summarize EPSPS classes

# Read the EPSPS classifications
epsps <- read.csv("results/epsps/epsps_classified.tsv", 
                  sep="\t", header=FALSE, stringsAsFactors=FALSE)

colnames(epsps) <- c("sequence_id", "primary_class", "classes", "sensitivity",
                     "identity_I", "identity_II", "identity_III", "identity_IV",
                     "is_unclassified", "is_too_divergent", "notes")

# Convert identity columns to numeric
epsps$identity_I <- as.numeric(epsps$identity_I)
epsps$identity_II <- as.numeric(epsps$identity_II)
epsps$identity_III <- as.numeric(epsps$identity_III)
epsps$identity_IV <- as.numeric(epsps$identity_IV)

cat("\n=== EPSPS Classification Summary (164 Genera) ===\n\n")

class_counts <- table(epsps$primary_class)
for (cls in names(class_counts)) {
  pct <- class_counts[cls] / nrow(epsps) * 100
  cat(sprintf("%s: %d (%.1f%%)\n", cls, class_counts[cls], pct))
}

cat("\n=== Confidence Gap (Identity_I - Identity_II) ===\n")
gap <- epsps$identity_I - epsps$identity_II
cat("Mean:  ", round(mean(gap, na.rm=T), 2), "\n")
cat("Median:", round(median(gap, na.rm=T), 2), "\n")
cat("Min:   ", round(min(gap, na.rm=T), 2), "\n")
cat("Max:   ", round(max(gap, na.rm=T), 2), "\n")

cat("\n=== Breakdown by Class ===\n")
for (cls in sort(unique(epsps$primary_class))) {
  n <- sum(epsps$primary_class == cls)
  cat(sprintf("%-20s: %3d\n", cls, n))
}

cat("\n=== Sample sizes ===\n")
cat("Total genera:", nrow(epsps), "\n")
cat("Sensitive (Class I):   ", sum(epsps$primary_class == "I"), "\n")
cat("Resistant (II+III+IV): ", sum(epsps$primary_class %in% c("II", "III", "IV")), "\n")

