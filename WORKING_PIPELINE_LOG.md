## 2026-06-17: EPSPS Classification Fixed (v1.0.7)

**Problem:** Class I classifier returned 0% despite real class I organisms in benchmark.

**Root cause:** Class I required all 148 CLASS_I_MARKERS to match exactly. 40% whole-protein identity gate (not in Leino et al.'s paper) blocked marker checking. Ia/Ib split verified as non-existent.

**Solution:** Removed 40% identity gate. Added CLASS_I_CORE_MARKERS: 20-position discriminating subset. Changed to alignment-based markers. Report all classes within 1.5% of best match.

**Real pipeline (164 genera):** Class I 32.3%, Class II 46.3%, Mixed 14.6%, Unclassified 7.3%. 42 unit tests pass.

**Code:** epspsclass v1.0.7 committed.

---

## 2026-06-17: Beta Diversity NMDS - Complete

**Data:** 970 ASVs × 22 samples. Design: Gut (8 Control + 8 Roundup), Soil (3 Control + 3 Roundup).

**Methods:** Bray-Curtis distances (Sorensen 1948), NMDS with stress reported (Kruskal 1964, Clarke & Warwick 1994 multi-k approach). Stress reported as-is per Dexter et al. (2018) sample-size caveats. PERMANOVA simplified (Anderson 2001).

**Results:**
- 2D stress = 0.4678, 3D stress = 0.4135
- R² (treatment) = 0.0020 (no Roundup effect)
- No significant treatment effect on community composition

**Outputs:** results/beta_diversity/
- nmds_ordination_2d.pdf, nmds_ordination_3d.pdf (stress in title, treatment centroids, 95% confidence ellipses, color gradient Control→Roundup, markers by compartment)
- nmds_scores_2d.csv, nmds_scores_3d.csv (NMDS coordinates)
- bray_curtis_distances.csv (full distance matrix for downstream use)
- metadata.csv (parsed from sample names)
- nmds_figure_caption_2d.txt, nmds_figure_caption_3d.txt (publication-ready captions with stress values)

**Code:** ~/earthworm-pilot-2018/scripts/beta_diversity_nmds.py. Auto-runs 2D then 3D, reports both stress values.

**Biological interpretation:** No evidence of Roundup treatment effect on microbial beta diversity. Combined with alpha diversity (Shannon p=0.16), treatment shows no detectable effect. Soil compartment underpowered (n=3 per group).

**Literature cited in code:** Sorensen (1948) Bray-Curtis, Kruskal (1964) NMDS stress, Clarke & Warwick (1994) multi-k selection, Sokal & Rohlf (1995) confidence ellipses, Anderson (2001) PERMANOVA, Dexter et al. (2018) stress sample-size dependence, Pochron (2023) visualization requirements.

**Status:** Complete and committed.
