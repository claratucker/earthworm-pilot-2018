#!/usr/bin/env Rscript
# 05_epsps_overlay.R -- join PICRUSt2-predicted aroA carriage and EPSPSClass
# sensitivity classes onto the community taxonomy and treatment.
#
# EXPLORATORY ONLY. Output is a hypothesis about community-level EPSPS class
# composition under glyphosate, built on a prediction (PICRUSt2) classified via
# reference proteins (EPSPSClass). NOT a measurement of EPSPS in the samples.
# See 04_epsps_exploratory.sh header for the full caveat and citations.

suppressMessages({
  library(phyloseq); library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

PROJECT <- path.expand("~/pilot2018")
EXP  <- file.path(PROJECT, "results/r/exported")
EOUT <- file.path(PROJECT, "results/epsps")
OUT  <- file.path(PROJECT, "results/r"); dir.create(OUT, showWarnings = FALSE)

# ---- Rebuild the phyloseq object (taxonomy + table + metadata) --------------
otu_raw <- read.delim(file.path(EXP,"feature-table.tsv"), skip=1, check.names=FALSE, row.names=1)
otu <- otu_table(as.matrix(otu_raw), taxa_are_rows=TRUE)
tax_raw <- read.delim(file.path(EXP,"taxonomy.tsv"), row.names=1, stringsAsFactors=FALSE)
tax_split <- tax_raw$Taxon %>% strsplit(";\\s*") %>%
  lapply(function(x){length(x)<-7; x}) %>% do.call(rbind,.)
colnames(tax_split) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rownames(tax_split) <- rownames(tax_raw)
TAX <- tax_table(as.matrix(tax_split))
meta <- read.delim(file.path(PROJECT,"data/metadata.tsv"), comment.char="#",
                   row.names=1, stringsAsFactors=FALSE)
ps <- phyloseq(otu, TAX, sample_data(meta))

# ---- 1. Which genera are predicted to carry aroA? (the list to curate refs) -
aroA <- read_tsv(file.path(EOUT,"asv_aroA_predicted.tsv"), show_col_types=FALSE)
aroA_pos <- aroA %>% filter(K00800_predicted_copies > 0) %>% pull(asv)
tax_df <- as.data.frame(as(tax_table(ps),"matrix"))
tax_df$asv <- rownames(tax_df)
genera_needed <- tax_df %>% filter(asv %in% aroA_pos) %>%
  distinct(Genus) %>% filter(!is.na(Genus), Genus!="") %>% arrange(Genus)
write_csv(genera_needed, file.path(OUT,"epsps_genera_to_curate.csv"))
cat("Genera predicted to carry aroA (curate one reference EPSPS protein each):\n")
print(as.data.frame(genera_needed))
cat(sprintf("\n%d genera total.\n", nrow(genera_needed)))
cat("Build refs/epsps/reference_epsps.faa, header format: >Genus|accession\n\n")

# ---- 2. If classification exists, overlay class onto genera/treatment -------
cls_path <- file.path(EOUT,"epsps_classified.tsv")
if (file.exists(cls_path)) {
  cls <- read_tsv(cls_path, show_col_types=FALSE) %>%
    mutate(Genus = sub("\\|.*$","", query_id))   # header was Genus|accession

  # Relative abundance of each genus per sample, joined to predicted class
  ps_g  <- tax_glom(ps, "Genus", NArm=TRUE)
  rel   <- transform_sample_counts(ps_g, function(x) x/sum(x))
  long  <- psmelt(rel) %>%
    left_join(cls %>% select(Genus, primary_class, sensitivity), by="Genus") %>%
    mutate(sensitivity = ifelse(is.na(sensitivity),"Unknown/Unpredicted",sensitivity))

  # Community-weighted predicted sensitivity fraction per sample
  sens_by_sample <- long %>%
    group_by(Sample, compartment, treatment, sensitivity) %>%
    summarise(rel_abund = sum(Abundance), .groups="drop")
  write_csv(sens_by_sample, file.path(OUT,"epsps_sensitivity_by_sample.csv"))

  p <- ggplot(sens_by_sample,
              aes(treatment, rel_abund, fill=sensitivity)) +
    geom_boxplot(outlier.shape=NA, alpha=.6, position=position_dodge()) +
    geom_point(position=position_dodge(width=.75), size=1.5) +
    facet_wrap(~compartment) +
    labs(title="PREDICTED EPSPS sensitivity composition (EXPLORATORY)",
         subtitle="PICRUSt2 aroA prediction x EPSPSClass on reference proteins; not measured",
         y="Community-weighted relative abundance") +
    theme_bw()
  ggsave(file.path(OUT,"fig_epsps_predicted_sensitivity.png"), p, width=8, height=4.5, dpi=150)
  cat("Wrote epsps_sensitivity_by_sample.csv and figure.\n")

  # Descriptive only: predicted resistant fraction, gut, control vs roundup
  res_gut <- sens_by_sample %>% filter(compartment=="gut", sensitivity=="Resistant")
  cat("\nPredicted RESISTANT fraction in gut (descriptive means):\n")
  print(res_gut %>% group_by(treatment) %>% summarise(mean_resistant=mean(rel_abund)))
  cat("\nNOTE: descriptive only. Do not attach a p-value to a prediction-on-prediction.\n")
} else {
  cat("No epsps_classified.tsv yet. Curate refs/epsps/reference_epsps.faa for the\n")
  cat("genera listed above, run 04_epsps_exploratory.sh, then re-run this script.\n")
}
