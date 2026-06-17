#!/usr/bin/env Rscript
# 03_pilot_analysis.R -- descriptive + pilot analysis of the 2018 earthworm
# gut/soil glyphosate dataset.
#
# DESIGN (from the mapping file):
#   compartment x treatment, UNBALANCED
#     gut  (.In)  : 8 control + 8 roundup = 16
#     soil (.Soil): 3 control + 3 roundup = 6
#
# STATISTICAL STANCE (deliberately conservative; this is a PILOT):
#   * Gut treatment contrast (control gut n=8 vs roundup gut n=8): TESTABLE.
#     PERMANOVA at n=8/group is defensible for a large/moderate effect
#     (Anderson 2001 Austral Ecol 26:32; the n>=8 rationale for soil community
#     heterogeneity follows Prosser 2010 Environ Microbiol 12:1806).
#   * Soil treatment contrast (n=3 vs 3) and gut-vs-soil compartment contrast:
#     DESCRIPTIVE ONLY. Reported as ordination + composition, no p-values used
#     as evidence (permutation space too small; would mislead).
#   * NO DESeq2/ANCOMBC differential abundance: dispersion is not estimable at
#     this n and would produce false precision.
#   * Primary deliverable is an EFFECT-SIZE estimate to power the real study.
#
# Packages: phyloseq (McMurdie & Holmes 2013 PLoS ONE 8:e61217),
#   vegan (Oksanen et al. 2022), ape, ggplot2.

suppressMessages({
  library(phyloseq); library(vegan); library(ape)
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
})

PROJECT <- path.expand("~/pilot2018")
EXP     <- file.path(PROJECT, "results/r/exported")
OUT     <- file.path(PROJECT, "results/r"); dir.create(OUT, showWarnings = FALSE)
set.seed(1)

# ---- Load exported QIIME artifacts ------------------------------------------
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
SAM <- sample_data(meta)

# The tree is the ROOTED tree exported from rooted-tree.qza; `qiime tools export`
# writes it as tree.nwk. (If this file is missing, you exported the unrooted tree
# or skipped the export -- re-run the export of rooted-tree.qza in step 02.)
tree_path <- file.path(EXP, "tree.nwk")
if (!file.exists(tree_path)) {
  stop("tree.nwk not found at ", tree_path,
       "\n  -> export rooted-tree.qza in 02_qiime2_upstream.sh: ",
       "qiime tools export --input-path rooted-tree.qza --output-path ", EXP)
}
tree <- read_tree(tree_path)

ps <- phyloseq(otu, TAX, SAM, tree)
cat("Loaded phyloseq object:\n"); print(ps)

# ---- Library-size sanity (gut host-removal effect shows up here) ------------
libsizes <- data.frame(sample = sample_names(ps),
                       reads = sample_sums(ps),
                       meta[sample_names(ps), c("compartment","treatment")])
write_csv(libsizes, file.path(OUT, "library_sizes.csv"))
cat("\nLibrary sizes after host/chloroplast removal:\n"); print(libsizes)

# ---- Rarefaction for diversity (single even depth) --------------------------
# Pilot: rarefy to the min sample depth for comparability. Report the depth.
min_depth <- min(sample_sums(ps))
cat(sprintf("\nRarefying to even depth = %d reads.\n", min_depth))
ps_rar <- rarefy_even_depth(ps, sample.size = min_depth,
                            rngseed = 1, replace = FALSE, verbose = FALSE)

# =============================================================================
# 1. ALPHA DIVERSITY (descriptive)
# =============================================================================
alpha <- estimate_richness(ps_rar, measures = c("Observed","Shannon"))
alpha$pd <- picante::pd(t(as(otu_table(ps_rar),"matrix")),
                        phy_tree(ps_rar), include.root = FALSE)$PD
alpha <- cbind(alpha, meta[rownames(alpha), c("compartment","treatment")])
write_csv(cbind(sample = rownames(alpha), alpha), file.path(OUT, "alpha_diversity.csv"))

p_alpha <- ggplot(alpha, aes(treatment, Shannon, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = .6) +
  geom_jitter(width = .15, size = 2) +
  facet_wrap(~compartment) +
  labs(title = "Shannon diversity (pilot; descriptive)",
       subtitle = "Gut n=8/treatment, Soil n=3/treatment") +
  theme_bw()
ggsave(file.path(OUT, "fig_alpha_shannon.png"), p_alpha, width = 7, height = 4, dpi = 150)

# =============================================================================
# 2. BETA DIVERSITY / ORDINATION (descriptive for all; test only gut)
# =============================================================================
bray <- phyloseq::distance(ps_rar, method = "bray")
wuf  <- phyloseq::distance(ps_rar, method = "wunifrac")

ord <- ordinate(ps_rar, method = "PCoA", distance = "bray")
p_ord <- plot_ordination(ps_rar, ord, color = "treatment", shape = "compartment") +
  geom_point(size = 4) +
  stat_ellipse(aes(group = treatment_compartment), type = "norm",
               linetype = 2, show.legend = FALSE) +
  labs(title = "PCoA (Bray-Curtis) -- points are the data, ellipses illustrative",
       subtitle = "No significance claimed for soil (n=3) or compartment contrast") +
  theme_bw()
# treatment_compartment may not survive sample_data subsetting; guard it:
if (!"treatment_compartment" %in% colnames(meta)) p_ord$layers <- p_ord$layers[-2]
ggsave(file.path(OUT, "fig_pcoa_bray.png"), p_ord, width = 7, height = 5, dpi = 150)

# =============================================================================
# 3. THE ONE TESTABLE CONTRAST: glyphosate effect WITHIN the gut (n=8 vs 8)
# =============================================================================
ps_gut <- subset_samples(ps_rar, compartment == "gut")
gut_meta <- data.frame(sample_data(ps_gut))
bray_gut <- phyloseq::distance(ps_gut, method = "bray")
wuf_gut  <- phyloseq::distance(ps_gut, method = "wunifrac")

set.seed(1)
adon_bray <- adonis2(bray_gut ~ treatment, data = gut_meta, permutations = 9999)
adon_wuf  <- adonis2(wuf_gut  ~ treatment, data = gut_meta, permutations = 9999)

# Dispersion check (PERMANOVA assumes comparable dispersion; betadisper guards
# against a location result that is really a spread difference). Anderson 2006
# Biometrics 62:245.
disp <- betadisper(bray_gut, gut_meta$treatment)
disp_test <- permutest(disp, permutations = 9999)

sink(file.path(OUT, "gut_treatment_permanova.txt"))
cat("=== GUT-ONLY glyphosate contrast (control gut n=8 vs roundup gut n=8) ===\n\n")
cat(">> PERMANOVA, Bray-Curtis:\n");   print(adon_bray)
cat("\n>> PERMANOVA, weighted UniFrac:\n"); print(adon_wuf)
cat("\n>> Betadisper (dispersion) test, Bray-Curtis:\n"); print(disp_test)
cat("\nThe Bray-Curtis R2 above is the PILOT EFFECT SIZE for powering the\n")
cat("full study (feed it into the micropower simulation in section 4).\n")
sink()
cat("\nGut PERMANOVA written. Bray R2 =", round(adon_bray$R2[1], 3),
    " p =", adon_bray$`Pr(>F)`[1], "\n")

# =============================================================================
# 4. POWER PROJECTION for the real study, using the pilot effect size
# =============================================================================
# Use the observed gut Bray-Curtis R2 as the anchor. Report what replication the
# full design would need to detect an effect of this magnitude with power ~0.8.
# If the micropower package is available, run the Kelly et al. 2015 simulation
# (Bioinformatics 31:2461); otherwise emit the R2 + a guidance note.
r2_obs <- adon_bray$R2[1]
sink(file.path(OUT, "power_projection.txt"))
cat("Observed pilot effect size (gut, Bray-Curtis PERMANOVA R2):", round(r2_obs,3), "\n\n")
if (requireNamespace("micropower", quietly = TRUE)) {
  cat("micropower available: run Kelly-style simulation across n per group,\n")
  cat("targeting the observed R2. (Kelly et al. 2015 Bioinformatics 31:2461.)\n")
} else {
  cat("micropower not installed. Heuristic: an R2 of this size is a 'large'\n")
  cat("community effect; n=8/group typically gives >0.8 power for R2>~0.15,\n")
  cat("while moderate effects (R2 0.05-0.10) need n~12-20/group. Install\n")
  cat("micropower for a defensible proposal figure (Kelly et al. 2015).\n")
}
sink()

# =============================================================================
# 5. COMPOSITION: gut vs soil, and which soil taxa appear enriched in the gut
# =============================================================================
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_rel   <- transform_sample_counts(ps_genus, function(x) x / sum(x))
top_genera <- names(sort(taxa_sums(ps_rel), decreasing = TRUE))[1:20]
df_bar <- psmelt(prune_taxa(top_genera, ps_rel))
p_bar <- ggplot(df_bar, aes(Sample, Abundance, fill = Genus)) +
  geom_col() + facet_grid(~compartment + treatment, scales = "free_x", space = "free_x") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
                     legend.position = "right", legend.text = element_text(size = 6)) +
  labs(title = "Top-20 genera (relative abundance)",
       subtitle = "Gut-vs-soil filter comparison; descriptive")
ggsave(file.path(OUT, "fig_genus_barplot.png"), p_bar, width = 11, height = 6, dpi = 150)

# Shared vs unique ASVs: which soil taxa survive gut transit (the filter question)
soil_asvs <- taxa_names(prune_taxa(taxa_sums(subset_samples(ps, compartment=="soil"))>0, ps))
gut_asvs  <- taxa_names(prune_taxa(taxa_sums(subset_samples(ps, compartment=="gut" ))>0, ps))
venn <- data.frame(
  category = c("soil_only","gut_only","shared"),
  n_asvs   = c(length(setdiff(soil_asvs, gut_asvs)),
               length(setdiff(gut_asvs, soil_asvs)),
               length(intersect(soil_asvs, gut_asvs))))
write_csv(venn, file.path(OUT, "shared_asvs_soil_gut.csv"))
cat("\nShared/unique ASV counts (gut-as-filter):\n"); print(venn)

cat("\nAll outputs written to", OUT, "\n")
cat("Remember to: aws s3 sync ~/pilot2018/results/ s3://earthworm-pilot-2018/results/\n")
