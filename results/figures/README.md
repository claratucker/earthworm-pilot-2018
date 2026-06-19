# Figures: 2018 Earthworm-Glyphosate Pilot Study

Numbering matches the 2026 glyphosate dose-response repo where the analysis is equivalent, for direct side-by-side comparison.

## Figure 1: Sequencing and Denoising
figure_1_denoising_stats.pdf
Read retention through DADA2 denoising and chimera removal. Input 20,528 to 161,535 reads per sample. Non-chimeric retention 44.6 to 90.4 percent, mean 57.4 percent.

## Figure 2: Community Composition (Genus Level)
figure_2_taxa_composition_genus.pdf
Per-sample stacked bars, genus level, one bar per sample, faceted by environment (gut, soil). Dashed line separates control from roundup within each environment panel. Top 20 genera by mean relative abundance shown individually; all remaining genera summed into a single "Other" category.

figure_2c_taxa_composition_genus_averaged.pdf
Treatment-averaged version of Figure 2: replicates collapsed to one bar per treatment per environment. Same taxon color palette as Figure 2. No error bars: a stacked bar's segment position is the cumulative sum of all segments below it, so a per-segment error bar would represent the combined variance of every taxon stacked underneath it, not that taxon alone (Munzner 2014, Visualization Analysis and Design). Per-taxon variance across treatment is shown in Figure 6 instead.

## Figure 2b: Community Composition (Family Level)
figure_2b_taxa_composition_family.pdf
Same as Figure 2, at family level.

figure_2d_taxa_composition_family_averaged.pdf
Same as Figure 2c, at family level.

## Figure 3: Rarefaction
figure_3_rarefaction.pdf
Rarefaction curves, all 22 samples (vegan::rarecurve, Oksanen et al. 2022). Depth set to 10,903 reads, the minimum sample depth; all 22 samples retained.

## Figure 4: Alpha Diversity
figure_4_alpha_diversity_by_environment.pdf
Shannon diversity and observed richness by environment and treatment. Gut significantly more diverse than soil (Kruskal-Wallis Shannon p = 0.0402, Richness p = 0.0045). No significant treatment effect within either environment.

## Figure 5: Beta Diversity
figure_5a_nmds_all_samples.pdf
NMDS, all 22 samples, Bray-Curtis distance. PERMANOVA environment effect R-squared = 0.0791, p = 0.0001. PERMDISP p = 0.9881 (homogeneous dispersion). NMDS stress = 0.468.

figure_5b_beta_diversity_gut.pdf
NMDS, gut only (n=16, testable contrast). Treatment effect R-squared = 0.0724, p = 0.2690, not significant. NMDS stress = 0.445.

figure_5c_beta_diversity_soil.pdf
NMDS, soil only (n=6, descriptive only, see Methods note on sample size). Treatment effect R-squared = 0.1944, p = 0.4990, not significant. NMDS stress = 0.373.

## Figure 6: Top Dose-Responsive Taxa (Gut Only)
figure_6_top_responders.pdf
Top 6 enriched and 6 depleted significant ASVs (DESeq2, Love et al. 2014, padj < 0.05) in gut, roundup vs control. Mean relative abundance with standard error by treatment. Restricted to gut (n=8/group); soil (n=3/group) lacks sufficient replication for reliable DESeq2 dispersion estimation and is not tested for differential abundance (see Methods note).

## Figure 7: Significant Genera (Gut Only)
figure_7_significant_genera.pdf
Box plots, the 17 genera containing at least one significant ASV (padj < 0.05) in the gut contrast. Box and point fill colors use the gut control/roundup palette (light green control, dark green roundup). Points beyond 1.5x IQR excluded from y-axis scaling for visibility (Tukey 1977); full values retained in deseq2_gut_roundup_vs_control.csv.

## Figure 8: Volcano Plot (Gut Only)
figure_8_volcano_plot.pdf
Log2 fold change (roundup vs control) versus -log10 adjusted p-value, all tested gut ASVs. Thresholds: |log2FC| >= 1, padj < 0.05 (Benjamini-Hochberg). 33 significant ASVs: 10 enriched, 23 depleted under roundup.

## Figure 9: EPSPS Sensitivity
figure_9_epsps_by_treatment.pdf
EPSPS Class I (glyphosate-sensitive) and Class II/III (glyphosate-resistant) relative abundance by environment and treatment. No significant treatment effect in either environment (all Wilcoxon p > 0.79).

## Methods note: differential abundance restricted to gut
DESeq2 (Love et al. 2014, Genome Biol 15:550) requires adequate replication per group to estimate per-gene dispersion reliably. At n=3 per group, soil dispersion estimates are not reliable; a single outlier sample can dominate the estimated variance for any given taxon. This mirrors the original pilot design rationale (see scripts/r/03_pilot_analysis.R header): the gut contrast (n=8/group) is treated as the one statistically testable comparison in this pilot, while soil (n=3/group) and the gut-versus-soil comparison are reported descriptively (composition and ordination only, Figures 2, 2b, 5a, 5c), with no p-values assigned. Figures 6, 7, and 8 (differential abundance) are gut-only for this reason.

## Methods note: error bars on stacked composition plots
There is no standard convention for placing error bars directly on segments within a single stacked bar, since a segment's vertical position is the cumulative sum of every segment below it; an error bar on one segment would represent the combined variance of everything stacked underneath it, not that taxon alone. Variance for individual taxa is shown separately, as small multiples with point and SE per treatment, in Figure 6.

## Generation scripts
scripts/r/02_denoising_stats_figure.R generates Figure 1.
scripts/r/08_generate_missing_figures.R generates Figures 2, 2b, 2c, 2d, 3, 6, 7, 8.
Figures 4, 5a, 5b, 5c, 9 come from the original pilot pipeline (scripts/r/03_pilot_analysis.R, scripts/r/05_epsps_overlay.R, scripts/r/06_epsps_summary_by_treatment.R, scripts/r/07_diversity_figures.R).

## Note on raw QIIME2 outputs

The large QIIME2 pipeline binaries (.qza, .qzv, raw and trimmed fastq.gz, ~983 MB total) are not stored in this git repository; they exceed GitHub's 100 MB per-file limit and are pipeline intermediates rather than independent source data. They are stored at:

s3://earthworm-pilot-2018/results/qiime2/

This includes demux.qza, demux-untrimmed.qza, trimmed.qza, forward.fastq.gz, table.qza, taxonomy.qza, rooted-tree.qza, and related intermediates. The small derived exports needed for downstream analysis (denoising stats, feature table, taxonomy, tree, all as plain text) are committed to this repo under results/qiime2_denoising_stats_export/ and results/r/exported/.
