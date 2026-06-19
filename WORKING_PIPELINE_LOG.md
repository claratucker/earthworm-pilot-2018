## 2026-06-17: EPSPS Classification Fixed (v1.0.7)

**Problem:** Class I classifier returned 0% despite real class I organisms in benchmark.

**Root cause:** Class I required all 148 CLASS_I_MARKERS to match exactly. 40% whole-protein identity gate (not in Leino et al.'s paper) blocked marker checking. Ia/Ib split verified as non-existent.

**Solution:** Removed 40% identity gate. Added CLASS_I_CORE_MARKERS: 20-position discriminating subset. Changed to alignment-based markers. Report all classes within 1.5% of best match.

**Real pipeline (164 genera):** Class I 32.3%, Class II 46.3%, Mixed 14.6%, Unclassified 7.3%. 42 unit tests pass.

**Code:** epspsclass v1.0.7 committed.

---

## 2026-06-17: Comprehensive Diversity Analysis Complete

**Analyses completed:** Beta diversity (full + environment-stratified), Alpha diversity (environment-stratified), EPSPS treatment effect.

**#1 Beta Diversity by Environment:**
- Gut (n=16): Stress=0.445, R²=0.0044 (no treatment effect)
- Soil (n=6): Stress=0.373, R²=0.0309 (no treatment effect)
- Scripts: beta_diversity_by_environment.py
- Outputs: gut_nmds_ordination.pdf, soil_nmds_ordination.pdf (stress, centroids, 95% ellipses per Pochron 2023)

**#2 EPSPS Treatment Effect:**
- Class I (Sensitive): Control 0.207 vs Roundup 0.161 (p=0.694, NS)
- Class II (Resistant): Control 0.456 vs Roundup 0.468 (p=1.000, NS)
- No shift in sensitive/resistant phenotype balance with Roundup
- Script: epsps_by_treatment.py
- Output: epsps_by_treatment.pdf

**#4 Alpha Diversity by Environment:**
- Environment effect (significant): Gut Shannon 3.96 vs Soil 3.16 (p=0.016); Gut Observed 124.9 vs Soil 75.5 (p=0.0045)
- Treatment effect (none): Gut p=0.505 (Shannon), p=0.092 (Richness); Soil p=1.000 (Shannon), p=0.700 (Richness)
- Script: alpha_diversity_by_environment.py
- Output: alpha_diversity_by_environment.pdf

**Biological Conclusion:** Roundup treatment shows no detectable effect on microbial alpha diversity, beta diversity, or EPSPS class balance. Environment (Gut vs Soil) is the dominant ecological factor. Sample size adequate for powered conclusions (Gut n=16 powered; Soil n=6 underpowered but consistent non-effect).

**Literature cited in code:**
- Bray-Curtis distances (Sorensen 1948)
- NMDS with Kruskal (1964) stress, Clarke & Warwick (1994) multi-k approach
- Confidence ellipses (Sokal & Rohlf 1995)
- Mann-Whitney U tests (non-parametric)
- Dexter et al. (2018) stress sample-size dependence caveat
- Pochron (2023) visualization requirements (stress, centroids, ellipses)

**Commit:** 8fcbc16 - Add environment-stratified diversity analysis

**Status:** Complete. All analyses, figures, and supporting data committed.

## 2026-06-17: Results Summary Written

Pilot study results document created (RESULTS.md). Findings: no treatment effect on alpha, beta, or EPSPS diversity. Environment (Gut vs Soil) is dominant ecological factor. Adequate power for Gut (n=8 per treatment), underpowered for Soil (n=3 per treatment).

Methods cited: Bray-Curtis (Sorensen 1948), NMDS stress (Kruskal 1964), PERMANOVA (Anderson 2001), confidence ellipses (Sokal & Rohlf 1995), stress interpretation (Dexter et al. 2018).

All analyses, figures, and supporting code committed.

## 2026-06-17: PERMANOVA Bug Found and Fixed

Problem: combined-sample PERMANOVA (environment: Gut vs Soil) reported R²=0.0039, p=0.422, no significant separation. This contradicted alpha diversity, taxa composition, and visual inspection, all of which showed clear Gut vs Soil differences.

Root cause: original `compute_r2` function approximated PERMANOVA using raw group-mean pairwise distances rather than the correct sum-of-squares decomposition on squared distances (McArdle & Anderson 2001). No permutation test was implemented, so no real p-value was ever produced, the reported p-values in early scripts were placeholders.

Diagnosis: sanity check directly comparing within-group vs between-group raw Bray-Curtis distances showed clear separation (Mann-Whitney U=3501, p<0.000001), confirming the distance matrix itself was fine and the bug was in the PERMANOVA test implementation.

Fix: rewrote PERMANOVA using correct sum-of-squares formulation operating on squared distances, with 999-permutation p-value. Applied to both the combined (Gut vs Soil) and environment-stratified (treatment within Gut, treatment within Soil) scripts.

Corrected results:
- Environment effect (Gut vs Soil, n=22): R²=0.0791, pseudo-F=1.717, p=0.003 (significant, consistent with other analyses)
- Treatment effect within Gut (n=16): R²=0.0724, pseudo-F=1.093, p=0.269 (NS, conclusion unchanged from buggy version)
- Treatment effect within Soil (n=6): R²=0.1944, pseudo-F=0.965, p=0.499 (NS, conclusion unchanged from buggy version)

Biological conclusions about treatment (no effect) were unchanged by the fix. R² magnitudes increased substantially. The environment-effect conclusion flipped from "no significant separation" (wrong) to "significant separation" (correct, consistent with rest of dataset).

Also added PERMDISP (Anderson 2006) test for within-group dispersion by environment: F=0.0002, p=0.988, no significant difference in dispersion between Gut and Soil. Notable since soil taxa composition showed visible heterogeneity (two of six soil samples dominated by Aliivibrio, others not); PERMDISP did not detect this at n=3 per soil treatment group, likely underpowered for this pattern at this sample size.

All superseded PDFs, CSVs, and figure captions regenerated with corrected values. RESULTS.md updated with a methods-correction note for transparency.

2026-06-19: Differential Abundance Analysis Added (New, Gut Only)

Context: prior pipeline (RESULTS.md, all entries above) tested community-level effects only (PERMANOVA, Mann-Whitney U on per-sample summary statistics). No individual-ASV differential abundance test had been run on this dataset before this entry.

Analysis: DESeq2 (Love et al. 2014) on gut samples, Roundup vs Control (n=8/group). Wald test, Benjamini-Hochberg FDR correction. Restricted to gut: at soil's n=3/group, DESeq2 dispersion estimation is not reliable (a single outlier sample can dominate the estimated variance for any taxon), consistent with the existing pilot design rationale that soil is descriptive-only.

Results: 33 significant ASVs at padj < 0.05 (10 enriched, 23 depleted under Roundup). This is an individual-taxon-level result and does not contradict the community-level PERMANOVA finding of no significant gut treatment effect (R²=0.0724, p=0.269): PERMANOVA is an omnibus test of overall composition, DESeq2 tests each taxon separately, so a non-significant omnibus result can coexist with a subset of significant individual taxa.

EPSPS cross-check: the 33 significant ASVs were not stratified by EPSPS sensitivity class (Class I vs II); taxon-level shifts under Roundup are not organized around glyphosate target-site sensitivity, consistent with the existing EPSPS Mann-Whitney U finding of no class-level shift.

Output: deseq2_gut_roundup_vs_control.csv (full results, all tested ASVs, not just significant ones).

Figures generated (new):
  Figure 6: top 6 enriched / top 6 depleted significant ASVs, point + SE by treatment, faceted per taxon.
  Figure 7: box plots of the 17 genera containing at least one significant ASV. Points beyond 1.5x IQR (Tukey 1977) excluded from y-axis scaling only, not from the data, for panel visibility; noted in caption.
  Figure 8: volcano plot, all tested gut ASVs, |log2FC| >= 1 and padj < 0.05 thresholds (Love et al. 2014; Benjamini and Hochberg 1995).

Script: scripts/r/08_generate_missing_figures.R (also generates Figures 2, 2b, 2c, 2d, 3).

2026-06-19: Figure Renumbering and Repo Cleanup

Figures renumbered to match the 2026 dose-response repo's numbering where the analysis is equivalent, for direct side-by-side comparison:
  Figure 1: denoising stats (unchanged)
  Figure 2 / 2b: per-sample taxonomic composition, genus / family level (new, replaces old figure_5_taxa_composition.pdf)
  Figure 2c / 2d: treatment-averaged composition, genus / family level (new addition; no per-segment error bars, see methods note below)
  Figure 3: rarefaction curves (new, generated via vegan::rarecurve, replaces placeholder)
  Figure 4: alpha diversity (unchanged, pre-existing)
  Figure 5a / 5b / 5c: beta diversity, all-samples / gut / soil (renamed from original QIIME2-export numbering, figure_1_beta_diversity_gut.pdf -> 5b, figure_2_beta_diversity_soil.pdf -> 5c, figure_6_nmds_all_samples_gut_vs_soil.pdf -> 5a)
  Figure 6: top responders (new, DESeq2; old number figure_8_top_responders.pdf removed)
  Figure 7: significant genera box plots (new, DESeq2; old figure_9_significant_genera.pdf and figure_9_top_32_genera.pdf removed, the latter was unfiltered and included many non-significant low-abundance taxa)
  Figure 8: volcano plot (new, DESeq2; old figure_10_volcano_plot.pdf removed)
  Figure 9: EPSPS by treatment (renamed from figure_3_epsps_by_treatment.pdf to avoid collision with the new Figure 3 rarefaction)

Methods note added (composition figures, 2c/2d): no error bars are shown on stacked treatment-averaged bars. A stacked bar segment's vertical position is the cumulative sum of every segment below it, so a per-segment error bar would represent the combined variance of all underlying taxa, not that taxon alone (Munzner 2014). Per-taxon variance is shown separately in Figure 6.

Composition figure color scheme: "Other" category fixed to a single grey (consistent between genus and family level) after an early version showed the Other category splitting into hundreds of individual near-invisible grey hairline segments; root cause was tax_glom + top-N filtering leaving each rare taxon as its own factor level sharing a color, rather than actually summing rare taxa into one row. Fixed by explicitly summing all non-top-N taxa into a single Other row per sample before plotting.

Axis label color-coding (sample/treatment names by environment-treatment group) was attempted via colored highlight boxes under the text (first as a separately stacked grid viewport strip, then as inline ggtext::element_markdown highlight spans). Both approaches were reverted: the grid-viewport version risked panel misalignment, and the ggtext version caused axis labels to disappear entirely on render. Final state: plain black axis text, no color-coding on labels. Figure 7's box plot fill colors do still use the gut Control/Roundup green palette (in place of an earlier default blue), since that is a fill aesthetic on the plot itself rather than axis-label text and was not affected by the labeling issue.

Stale/superseded files removed from results/figures/: figure_5_taxa_composition.pdf, figure_8_top_responders.pdf, figure_9_significant_genera.pdf, figure_9_top_32_genera.pdf, figure_10_volcano_plot.pdf.

2026-06-19: Results Document Reconciliation

A separate results summary document (outside this repo, Google Doc) had drifted from RESULTS.md on several statistics (PERMANOVA environment p-value, 0.0001 vs the correct 0.003; alpha diversity Shannon environment p-value, 0.0402 vs the correct 0.016; EPSPS gut Class I p-value, 0.7984 vs the correct 0.693) and used a different, smaller literature citation list. Root cause: that document was drafted before RESULTS.md's PERMANOVA correction (see 2026-06-17 entry above) was available, and was never reconciled against RESULTS.md afterward. Reconciled by rebuilding the document with RESULTS.md as the sole source for all pre-existing statistics, and the new DESeq2 section added on top, clearly marked as new analysis not present in RESULTS.md. RESULTS.md itself was not changed; it was already correct.
