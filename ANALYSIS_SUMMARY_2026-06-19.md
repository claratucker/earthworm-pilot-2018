# 2018 Earthworm Glyphosate Pilot: Analysis Summary (Update)

Follow-up to ANALYSIS_SUMMARY_2026-06-17.md. That document's "Next Priority" asked whether glyphosate treatment shifts the sensitive/resistant EPSPS balance, and whether soil differs from gut microbiome. Both are now answered below, along with a PERMANOVA bug fix and new differential abundance analysis.

## Project Overview

Unchanged from 2026-06-17: glyphosate-exposed earthworms (Lumbricus terrestris) and associated soil microbiome. Treatments: Control (C) vs Roundup-treated (R). Environments: earthworm gut (In) and surrounding soil (Soil). 22 samples total: 8 Control Gut, 8 Roundup Gut, 3 Control Soil, 3 Roundup Soil.

## Completed Since 2026-06-17

### 1. PERMANOVA Bug Fixed

The original environment-effect PERMANOVA (Gut vs Soil) reported R²=0.0039, p=0.422, no significant separation, contradicting alpha diversity, taxonomic composition, and visual inspection. Root cause: the original implementation approximated the test using raw group-mean pairwise distances instead of the correct sum-of-squares decomposition (McArdle and Anderson 2001), and never computed a real permutation p-value. Fixed and re-run with 999 permutations.

Corrected results:
- Environment effect (Gut vs Soil, n=22): R²=0.0791, pseudo-F=1.717, p=0.003 (significant)
- Treatment effect within Gut (n=16): R²=0.0724, pseudo-F=1.093, p=0.269 (NS, unchanged conclusion)
- Treatment effect within Soil (n=6): R²=0.1944, pseudo-F=0.965, p=0.499 (NS, unchanged conclusion)

Treatment conclusions were unchanged by the fix; the environment-effect conclusion flipped from non-significant (wrong) to significant (correct). PERMDISP (Anderson 2006) added as a companion test: F=0.0002, p=0.988, no significant dispersion difference between environments.

### 2. Treatment-Level Statistics (answers the "Next Priority" from 2026-06-17)

Does glyphosate shift the sensitive/resistant EPSPS balance? No.
- Gut: Class I mean Control 0.254 vs Roundup 0.184 (p=0.693, NS). Class II mean Control 0.495 vs Roundup 0.539 (p=1.000, NS).
- Soil: Class I mean Control 0.081 vs Roundup 0.100 (p=1.000, NS). Class II mean Control 0.353 vs Roundup 0.277 (p=1.000, NS).

Does soil differ from gut microbiome? Yes, substantially.
- Alpha diversity: Gut Shannon 3.96 ± 0.51 vs Soil 3.16 ± 0.84 (p=0.016). Gut Observed richness 124.9 ± 45.2 vs Soil 75.5 ± 16.0 (p=0.0045).
- Beta diversity: environment effect R²=0.0791, p=0.003 (see PERMANOVA fix above).
- Taxonomic composition: visibly distinct between environments; Aliivibrio notably more abundant in Soil, concentrated in 2 of 6 soil samples (one per treatment), consistent with microhabitat heterogeneity rather than a treatment or processing effect.

### 3. Differential Abundance Analysis (New)

Not part of the original pipeline. The analyses above (PERMANOVA, Mann-Whitney U) all test community-level or per-sample summary effects. DESeq2 (Love et al. 2014) was added to test individual ASVs directly, restricted to gut (n=8/group), since soil's n=3/group does not support reliable dispersion estimation.

Results: 33 significant ASVs at padj < 0.05 (10 enriched, 23 depleted under Roundup). PERMANOVA (p=0.269) is an omnibus test of overall composition, while DESeq2 finds individual taxa that move significantly even when the community as a whole does not shift detectably. Cross-checked against EPSPS class: the 33 significant ASVs are not stratified by Class I/II, so the taxon-level response is not organized around glyphosate target-site sensitivity.

Output: deseq2_gut_roundup_vs_control.csv (all tested ASVs). Figures: top responders (point + SE), significant-genera box plots (17 genera with at least one significant ASV), volcano plot (|log2FC| >= 1, padj < 0.05 thresholds).

## Updated Key Results

EPSPS classification (unchanged from 2026-06-17, v1.0.7): Class I 32.1%, Class II 46.1%, Class III 3.0%, Class IV 10.9%, Mixed ~14.6%, Unclassified 7.3%, across 164 genera.

Community-level treatment effect: none detected, in either environment, by PERMANOVA, alpha diversity, or EPSPS class balance.

Individual-taxon treatment effect (gut only): 33 of tested ASVs significant by DESeq2, not stratified by EPSPS class.

Environment effect: significant and consistent across every analysis (alpha diversity, beta diversity, composition).

## Key Files

PERMANOVA correction outputs: results/beta_diversity_by_environment/, results/beta_diversity_combined/ (regenerated with corrected statistics).

EPSPS treatment-level results (unchanged paths from 2026-06-17): ~/pilot2018/results/epsps/epsps_classified.tsv, ~/pilot2018/results/epsps/epsps_summary.txt.

DESeq2 output (new): results/r/deseq2_gut_roundup_vs_control.csv.

Figures (renumbered to match 2026 dose-response repo): results/figures/, see results/figures/README.md for the full figure-by-figure breakdown.

Results document: RESULTS.md (community-level pipeline, with PERMANOVA correction note); separate results summary document reconciled against RESULTS.md this session with DESeq2 section added on top (see WORKING_PIPELINE_LOG.md, 2026-06-19 entries).

## Remaining Open Items

Functional gene analysis beyond EPSPS (nifH, narG/nirK, ppx, chitinase) — not started, as of 2026-06-17, still not started.

Active-site residue classifier (independent EPSPS validation) — not started, as of 2026-06-17, still not started.

Comparison section in the results summary document, contextualizing 2018 pilot findings against the 2026 dose-response study — not yet written.
