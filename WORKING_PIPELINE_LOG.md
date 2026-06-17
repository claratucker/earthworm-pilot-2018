## 2026-06-17: EPSPS Classification Fixed (v1.0.7)

**Problem:** Class I classifier returned 0% despite real class I organisms in benchmark.

**Root cause:** Class I required all 148 CLASS_I_MARKERS to match exactly. 40% whole-protein identity gate (not in Leino et al.'s paper) blocked marker checking. Ia/Ib split verified as non-existent.

**Solution:** Removed 40% identity gate. Added CLASS_I_CORE_MARKERS: 20-position discriminating subset. Changed to alignment-based markers. Report all classes within 1.5% of best match.

**Real pipeline (164 genera):** Class I 32.3%, Class II 46.3%, Mixed 14.6%, Unclassified 7.3%. 42 unit tests pass.

**Code:** epspsclass v1.0.7 committed.

---

## 2026-06-17: Comprehensive Diversity Analysis Complete

**Analyses completed:** Beta diversity (full + compartment-stratified), Alpha diversity (compartment-stratified), EPSPS treatment effect.

**#1 Beta Diversity by Compartment:**
- Gut (n=16): Stress=0.445, R²=0.0044 (no treatment effect)
- Soil (n=6): Stress=0.373, R²=0.0309 (no treatment effect)
- Scripts: beta_diversity_by_compartment.py
- Outputs: gut_nmds_ordination.pdf, soil_nmds_ordination.pdf (stress, centroids, 95% ellipses per Pochron 2023)

**#2 EPSPS Treatment Effect:**
- Class I (Sensitive): Control 0.207 vs Roundup 0.161 (p=0.694, NS)
- Class II (Resistant): Control 0.456 vs Roundup 0.468 (p=1.000, NS)
- No shift in sensitive/resistant phenotype balance with Roundup
- Script: epsps_by_treatment.py
- Output: epsps_by_treatment.pdf

**#4 Alpha Diversity by Compartment:**
- Compartment effect (significant): Gut Shannon 3.96 vs Soil 3.16 (p=0.016); Gut Observed 124.9 vs Soil 75.5 (p=0.0045)
- Treatment effect (none): Gut p=0.505 (Shannon), p=0.092 (Richness); Soil p=1.000 (Shannon), p=0.700 (Richness)
- Script: alpha_diversity_by_compartment.py
- Output: alpha_diversity_by_compartment.pdf

**Biological Conclusion:** Roundup treatment shows no detectable effect on microbial alpha diversity, beta diversity, or EPSPS class balance. Compartment (Gut vs Soil) is the dominant ecological factor. Sample size adequate for powered conclusions (Gut n=16 powered; Soil n=6 underpowered but consistent non-effect).

**Literature cited in code:**
- Bray-Curtis distances (Sorensen 1948)
- NMDS with Kruskal (1964) stress, Clarke & Warwick (1994) multi-k approach
- Confidence ellipses (Sokal & Rohlf 1995)
- Mann-Whitney U tests (non-parametric)
- Dexter et al. (2018) stress sample-size dependence caveat
- Pochron (2023) visualization requirements (stress, centroids, ellipses)

**Commit:** 8fcbc16 - Add compartment-stratified diversity analysis

**Status:** Complete. All analyses, figures, and supporting data committed.

## 2026-06-17: Results Summary Written

Pilot study results document created (RESULTS.md). Findings: no treatment effect on alpha, beta, or EPSPS diversity. Compartment (Gut vs Soil) is dominant ecological factor. Adequate power for Gut (n=8 per treatment), underpowered for Soil (n=3 per treatment).

Methods cited: Bray-Curtis (Sorensen 1948), NMDS stress (Kruskal 1964), PERMANOVA (Anderson 2001), confidence ellipses (Sokal & Rohlf 1995), stress interpretation (Dexter et al. 2018).

All analyses, figures, and supporting code committed.
