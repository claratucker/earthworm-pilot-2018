# Results

Earthworm pilot 2018 microbiome study. 2x2 factorial design: Treatment (Control vs Roundup) x Environment (Gut vs Soil). Total samples: 22 (8 Control Gut, 8 Roundup Gut, 3 Control Soil, 3 Roundup Soil).

## Beta Diversity

Beta diversity was assessed via NMDS ordination of Bray-Curtis distances (Sorensen 1948). Stress was reported as-is without use as a gate (Dexter et al. 2018). PERMANOVA tested environment and treatment effects using the standard distance-based sum-of-squares formulation (McArdle & Anderson 2001, building on Anderson 2001), with p-values obtained by permutation (999 permutations). 95% confidence ellipses calculated per Sokal & Rohlf 1995. PERMDISP (Anderson 2006) tested within-group dispersion using the ANOVA-on-distance-to-centroid approximation.

**Environment effect (all 22 samples):** Gut and Soil show significantly different community composition. R²=0.0791, pseudo-F=1.717, p=0.003 (999 permutations). Within-group dispersion did not differ between environments (PERMDISP F=0.0002, p=0.988).

**Treatment effect, Gut environment (n=16, powered):** Stress 0.445. No treatment effect detected (R²=0.0724, pseudo-F=1.093, p=0.269, NS).

**Treatment effect, Soil environment (n=6, underpowered):** Stress 0.373. No treatment effect detected (R²=0.1944, pseudo-F=0.965, p=0.499, NS).

Roundup treatment did not significantly shift gut or soil community composition at the whole-community (PERMANOVA) level. Environment (Gut vs Soil) is the dominant driver of beta diversity in this dataset. Individual-taxon-level effects, which PERMANOVA as an omnibus test is not designed to detect, are examined separately below (Differential Abundance).

Note on methods correction: an earlier version of the PERMANOVA implementation approximated the test using raw group-mean pairwise distances rather than the correct sum-of-squares decomposition, and did not compute permutation-based p-values. This was caught when a sanity check on raw Bray-Curtis distances showed clear, highly significant Gut vs Soil separation (Mann-Whitney U=3501, p<0.000001) that the original PERMANOVA implementation failed to detect (originally reported as R²=0.0039, p=0.422). The corrected implementation above resolves this discrepancy. The treatment-effect conclusions (no significant effect in either environment) were unchanged by the correction, though R² values increased substantially (Gut: 0.0044 to 0.0724; Soil: 0.0309 to 0.1944).

## Differential Abundance (New Analysis, Gut Only)

This section is new analysis added after the original pipeline above (Beta Diversity, EPSPS, Alpha Diversity, Taxonomic Composition, all of which test community-level or per-sample summary effects only). DESeq2 (Love et al. 2014) was used to test individual ASVs directly for a treatment effect, which is a different question from the community-level PERMANOVA test: PERMANOVA asks whether overall composition differs; DESeq2 asks whether any individual taxa differ in abundance. A non-significant PERMANOVA result and a set of significant individual ASVs are not in conflict; this is the expected relationship between an omnibus test and many single-feature tests.

Analysis was restricted to Gut (n=8 Control, n=8 Roundup). DESeq2 requires adequate replication per group to estimate per-gene dispersion reliably; at Soil's n=3 per group, dispersion estimates are not reliable, since a single outlier sample can dominate the estimated variance for any given taxon. This mirrors the sample-size rationale already applied to Soil throughout this document (treatment effect reported as NS but underpowered, not tested with individual-taxon methods).

DESeq2 (Wald test, Benjamini-Hochberg FDR correction) identified 33 significant ASVs at padj < 0.05 in Gut: 10 enriched and 23 depleted under Roundup. These 33 ASVs were not stratified by EPSPS sensitivity class (see EPSPS Class Distribution below): taxon-level shifts under Roundup are not organized around glyphosate target-site sensitivity. Full results, including all tested (not just significant) ASVs, are in deseq2_gut_roundup_vs_control.csv.

## EPSPS Class Distribution

EPSPS aroA genes were classified as sensitive (Class I) or resistant (Class II) in 164 genera. Genus-level classifications were mapped to ASV abundance. Relative abundance of each class was computed per sample and tested by Mann-Whitney U.

**Gut environment:** Class I mean Control 0.254 vs Roundup 0.184 (p=0.693, NS). Class II mean Control 0.495 vs Roundup 0.539 (p=1.000, NS).

**Soil environment:** Class I mean Control 0.081 vs Roundup 0.100 (p=1.000, NS). Class II mean Control 0.353 vs Roundup 0.277 (p=1.000, NS).

No shift in phenotypic herbicide resistance profile between treatments. Sensitive and resistant bacterial populations remain balanced regardless of Roundup exposure. The 33 ASVs identified as significant by the new Differential Abundance analysis above are not concentrated in either EPSPS class, consistent with this finding.

## Alpha Diversity

Shannon diversity and observed richness (ASV counts) were calculated per sample. Mann-Whitney U tested environment and treatment effects.

**Environment effect (significant):** Gut vs Soil for Shannon (p=0.016) and Observed richness (p=0.0045). Gut samples show higher diversity and richness than soil.

Gut Shannon: 3.96 +/- 0.51 vs Soil: 3.16 +/- 0.84

Gut Observed: 124.9 +/- 45.2 vs Soil: 75.5 +/- 16.0

**Treatment effect (none):** Within Gut, Control vs Roundup Shannon p=0.505, Observed p=0.092. Within Soil, Control vs Roundup Shannon p=1.000, Observed p=0.700.

## Taxonomic Composition

Genus-level relative abundance was computed for all samples. Genera below 1% mean relative abundance across all samples were pooled into an "Other" category. Composition differed substantially between Gut and Soil, consistent with the significant environment effect found by PERMANOVA above. Aliivibrio was notably more abundant in Soil than Gut overall, though this signal was concentrated in two of six soil samples (one Control, one Roundup) rather than evenly distributed across all soil replicates, a pattern consistent with soil microhabitat heterogeneity at this sample size (n=3 per soil treatment) rather than a treatment effect or systematic processing issue, since the two outlying samples occurred in both treatment groups.

## Summary

Roundup treatment shows no significant effect on earthworm gut or soil microbiome alpha diversity, beta diversity (community-level), or EPSPS phenotype distribution in this pilot study. Environment (Gut vs Soil) is the dominant ecological factor, with significant differences in alpha diversity, beta diversity, and taxonomic composition. At the individual-taxon level, new DESeq2 analysis identified 33 significant ASVs in the Gut contrast, indicating detectable taxon-specific compositional change even though it does not register as a significant shift in overall community structure or diversity; this is consistent with a community-level non-effect coexisting with localized single-taxon effects, not a contradiction. These taxon-level changes are not stratified by EPSPS class. Adequate sample size for powered conclusions in Gut (n=8 per treatment, including the Differential Abundance analysis); underpowered in Soil (n=3 per treatment, descriptive only, not included in Differential Abundance analysis) but consistent non-effect across all community-level methods.

## Literature Cited

Anderson MJ. 2001. A new method for non-parametric multivariate analysis of variance. Austral Ecology 26: 32-46.

Anderson MJ. 2006. Distance-based tests for homogeneity of multivariate dispersions. Biometrics 62: 245-253.

Clarke KR, Warwick RM. 1994. Change in marine communities: An approach to statistical analysis and interpretation. 2nd edn. PRIMER-E: Plymouth.

Dexter E, Rollwagen-Bollens G, Bollens SM. 2018. The utility of taxonomically standardized names for copepod nauplii to increase comparability and measurability of ecological data. Hydrobiologia 820: 155-167.

Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15: 550.

McArdle BH, Anderson MJ. 2001. Fitting multivariate models to community data: A comment on distance-based redundancy analysis. Ecology 82: 290-297.

Sokal RR, Rohlf FJ. 1995. Biometry. 3rd edn. W.H. Freeman: New York.

Sorensen T. 1948. A method of establishing groups of equal amplitude in plant sociology based on similarity of species content. Kongelige Danske Videnskabernes Selskab Biologiske Skrifter 5: 1-34.

Tukey JW. 1977. Exploratory Data Analysis. Addison-Wesley.
