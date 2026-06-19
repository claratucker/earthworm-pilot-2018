# Figures - 2018 Earthworm-Glyphosate Pilot Study

All analysis figures are stored here as PDFs, organized by analysis type.

## Sequencing & Quality Control

**figure_1_denoising_stats.pdf**
Read retention through DADA2 denoising and chimera removal pipeline. Input: 20,528-161,535 reads/sample. Non-chimeric retention: 44.6-90.4% (mean 57.4%).

## Community Composition & Taxonomy

**figure_5_taxa_composition.pdf**
Treatment-averaged relative abundance at genus level. Gut dominated by Staphylococcus and Pseudomonas; soil shows Allivibrio and Colwellia as major contributors.

## Diversity Analysis

**figure_3_rarefaction.pdf**
Rarefaction curves showing sequencing depth saturation. Rarefaction depth: 10,903 reads. All samples plateau, indicating adequate sequencing depth.

**figure_4_alpha_diversity_by_environment.pdf**
Shannon diversity and observed richness by environment and treatment. Gut significantly more diverse than soil (Kruskal-Wallis p=0.0402).

## Beta Diversity / Ordination

**figure_6_nmds_all_samples_gut_vs_soil.pdf**
NMDS ordination of all 22 samples (gut and soil, control and Roundup). PERMANOVA environment effect R²=0.0791, p=0.0001. NMDS stress=0.468.

**figure_1_beta_diversity_gut.pdf**
NMDS gut samples only (n=16, control and Roundup). Treatment effect R²=0.0724, p=0.2690 (not significant). NMDS stress=0.445.

**figure_2_beta_diversity_soil.pdf**
NMDS soil samples only (n=6, control and Roundup, descriptive). Treatment effect R²=0.1944, p=0.4990 (not significant). NMDS stress=0.373.

## Differential Abundance (Gut, Control vs Roundup)

**figure_8_top_responders.pdf**
Top 6 enriched and 6 depleted taxa under Roundup in gut samples. Mean relative abundance with error bars.

**figure_9_top_32_genera.pdf**
Distribution of 32 most abundant genera in gut samples by treatment. Box plots with individual replicates.

**figure_10_volcano_plot.pdf**
Log2 fold change (Roundup/Control) versus mean relative abundance for all gut ASVs. Blue=enriched in Roundup (log2FC>1), red=depleted (log2FC<-1), gray=minimal change.

## Functional Analysis

**figure_3_epsps_by_treatment.pdf**
EPSPS sensitivity class composition (Class I sensitive vs Class II/III resistant) by environment and treatment. No significant treatment effects (all Wilcoxon p>0.79).

## Generation Scripts

Figure generation scripts are in `scripts/r/`:
- `02_denoising_stats_figure.R` - Generates Figure 1 (denoising stats)
- `08_generate_missing_figures.R` - Generates Figures 3, 8, 9, 10 (rarefaction, responders, genera, volcano)

Generated: June 19, 2026
