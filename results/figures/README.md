# Figures: 2018 Earthworm-Glyphosate Pilot Study

Numbering matches the 2026 glyphosate dose-response repo where the analysis is equivalent, for direct side-by-side comparison.

## Figure 1: Sequencing and Denoising
figure_1_denoising_stats.pdf
Read retention through DADA2 denoising and chimera removal. Input 20,528 to 161,535 reads per sample. Non-chimeric retention 44.6 to 90.4 percent, mean 57.4 percent.

## Figure 2: Community Composition
figure_2_taxa_composition_genus.pdf
Per-sample taxonomic composition at genus level, stacked bars, one bar per sample. Sample labels colored by environment and treatment: light green control gut, dark green roundup gut, light brown control soil, dark brown roundup soil. Dashed line separates control from roundup within each environment panel.

figure_2b_taxa_composition_family.pdf
Same as above at family level.

figure_2c_taxa_composition_genus_averaged.pdf
Treatment-averaged version of Figure 2, replicates collapsed to one bar per treatment per environment. Same taxon color palette as Figure 2. No error bars: a stacked bar's segment position is the cumulative sum of all segments below it, so a per-segment error bar would represent the combined variance of every taxon underneath it, not that taxon alone. Per-taxon variance across treatment is shown in Figure 6.

figure_2d_taxa_composition_family_averaged.pdf
Same as above at family level.

## Figure 3: Rarefaction
figure_3_rarefaction.pdf
Rarefaction curves, all 22 samples, depth 10,903 reads (minimum sample depth, all samples retained).

## Figure 4: Alpha Diversity
figure_4_alpha_diversity_by_environment.pdf
Shannon diversity and observed richness by environment and treatment. Gut significantly more diverse than soil (Kruskal-Wallis p = 0.0402).

## Figure 5: Beta Diversity
figure_5a_nmds_all_samples.pdf
NMDS, all 22 samples, Bray-Curtis. PERMANOVA environment effect R-squared 0.0791, p = 0.0001.

figure_5b_beta_diversity_gut.pdf
NMDS, gut only (n=16). Treatment effect R-squared 0.0724, p = 0.2690, not significant.

figure_5c_beta_diversity_soil.pdf
NMDS, soil only (n=6, descriptive). Treatment effect R-squared 0.1944, p = 0.4990, not significant.

## Figure 6: Top Dose-Responsive Taxa
figure_6_top_responders.pdf
Top 6 enriched and 6 depleted significant ASVs (DESeq2, padj < 0.05) in gut, control vs roundup. Mean relative abundance with standard error by treatment.

## Figure 7: Significant Genera
figure_7_significant_genera.pdf
Box plots, the 17 genera containing at least one significant ASV (padj < 0.05). Points beyond 1.5x IQR excluded from y-axis scaling for visibility (Tukey 1977); full values in deseq2_gut_roundup_vs_control.csv.

## Figure 8: Volcano Plot
figure_8_volcano_plot.pdf
Log2 fold change (roundup vs control) versus -log10 adjusted p-value, gut ASVs. Thresholds |log2FC| >= 1, padj < 0.05 (Benjamini-Hochberg).

## Figure 9: EPSPS Sensitivity
figure_9_epsps_by_treatment.pdf
EPSPS Class I (sensitive) and Class II/III (resistant) relative abundance by environment and treatment. No significant treatment effect (all Wilcoxon p > 0.79).

## Generation scripts
scripts/r/02_denoising_stats_figure.R generates Figure 1.
scripts/r/08_generate_missing_figures.R generates Figures 2, 2b, 2c, 2d, 3, 6, 7, 8.
Figures 4, 5a, 5b, 5c, 9 come from the original QIIME2/R pilot pipeline (scripts/r/03_pilot_analysis.R, scripts/r/05_epsps_overlay.R, scripts/r/06_epsps_summary_by_treatment.R, scripts/r/07_diversity_figures.R).
