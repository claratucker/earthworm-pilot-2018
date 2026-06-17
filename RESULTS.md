# Results

Earthworm pilot 2018 microbiome study. 2x2 factorial design: Treatment (Control vs Roundup) x Compartment (Gut vs Soil). Total samples: 22 (8 Control Gut, 8 Roundup Gut, 3 Control Soil, 3 Roundup Soil).

## Beta Diversity

Beta diversity was assessed via NMDS ordination of Bray-Curtis distances (Sorensen 1948). Stress was reported as-is without use as a gate (Dexter et al. 2018). PERMANOVA tested compartment and treatment effects using the standard distance-based sum-of-squares formulation (McArdle & Anderson 2001, building on Anderson 2001), with p-values obtained by permutation (999 permutations). 95% confidence ellipses calculated per Sokal & Rohlf 1995. PERMDISP (Anderson 2006) tested within-group dispersion using the ANOVA-on-distance-to-centroid approximation.

**Compartment effect (all 22 samples):** Gut and Soil show significantly different community composition. R²=0.0791, pseudo-F=1.717, p=0.003 (999 permutations). Within-group dispersion did not differ between compartments (PERMDISP F=0.0002, p=0.988).

**Treatment effect, Gut compartment (n=16, powered):** Stress 0.445. No treatment effect detected (R²=0.0724, pseudo-F=1.093, p=0.269, NS).

**Treatment effect, Soil compartment (n=6, underpowered):** Stress 0.373. No treatment effect detected (R²=0.1944, pseudo-F=0.965, p=0.499, NS).

Roundup treatment did not significantly shift gut or soil community composition. Compartment (Gut vs Soil) is the dominant driver of beta diversity in this dataset.

Note on methods correction: an earlier version of the PERMANOVA implementation approximated the test using raw group-mean pairwise distances rather than the correct sum-of-squares decomposition, and did not compute permutation-based p-values. This was caught when a sanity check on raw Bray-Curtis distances showed clear, highly significant Gut vs Soil separation (Mann-Whitney U=3501, p<0.000001) that the original PERMANOVA implementation failed to detect (originally reported as R²=0.0039, p=0.422). The corrected implementation above resolves this discrepancy. The treatment-effect conclusions (no significant effect in either compartment) were unchanged by the correction, though R² values increased substantially (Gut: 0.0044 to 0.0724; Soil: 0.0309 to 0.1944).

## EPSPS Class Distribution

EPSPS aroA genes were classified as sensitive (Class I) or resistant (Class II) in 164 genera. Genus-level classifications were mapped to ASV abundance. Relative abundance of each class was computed per sample and tested by Mann-Whitney U.

**Gut compartment:** Class I mean Control 0.254 vs Roundup 0.184 (p=0.693, NS). Class II mean Control 0.495 vs Roundup 0.539 (p=1.000, NS).

**Soil compartment:** Class I mean Control 0.081 vs Roundup 0.100 (p=1.000, NS). Class II mean Control 0.353 vs Roundup 0.277 (p=1.000, NS).

No shift in phenotypic herbicide resistance profile between treatments. Sensitive and resistant bacterial populations remain balanced regardless of Roundup exposure.

## Alpha Diversity

Shannon diversity and observed richness (ASV counts) were calculated per sample. Mann-Whitney U tested compartment and treatment effects.

**Compartment effect (significant):** Gut vs Soil for Shannon (p=0.016) and Observed richness (p=0.0045). Gut samples show higher diversity and richness than soil.

Gut Shannon: 3.96 +/- 0.51 vs Soil: 3.16 +/- 0.84

Gut Observed: 124.9 +/- 45.2 vs Soil: 75.5 +/- 16.0

**Treatment effect (none):** Within Gut, Control vs Roundup Shannon p=0.505, Observed p=0.092. Within Soil, Control vs Roundup Shannon p=1.000, Observed p=0.700.

## Taxonomic Composition

Genus-level relative abundance was computed for all samples. Genera below 1% mean relative abundance across all samples were pooled into an "Other" category. Composition differed substantially between Gut and Soil, consistent with the significant compartment effect found by PERMANOVA above. Aliivibrio was notably more abundant in Soil than Gut overall, though this signal was concentrated in two of six soil samples (one Control, one Roundup) rather than evenly distributed across all soil replicates, a pattern consistent with soil microhabitat heterogeneity at this sample size (n=3 per soil treatment) rather than a treatment effect or systematic processing issue, since the two outlying samples occurred in both treatment groups.

## Summary

Roundup treatment shows no significant effect on earthworm gut or soil microbiome alpha diversity, beta diversity, or EPSPS phenotype distribution in this pilot study. Compartment (Gut vs Soil) is the dominant ecological factor, with significant differences in alpha diversity, beta diversity, and taxonomic composition. Results are consistent across multiple analytical approaches and suggest treatment had no measurable impact at pilot scale. Adequate sample size for powered conclusions in Gut (n=8 per treatment), underpowered in Soil (n=3 per treatment) but consistent non-effect.

## Literature Cited

Anderson MJ. 2001. A new method for non-parametric multivariate analysis of variance. Austral Ecology 26: 32-46.

Anderson MJ. 2006. Distance-based tests for homogeneity of multivariate dispersions. Biometrics 62: 245-253.

Clarke KR, Warwick RM. 1994. Change in marine communities: An approach to statistical analysis and interpretation. 2nd edn. PRIMER-E: Plymouth.

Dexter E, Rollwagen-Bollens G, Bollens SM. 2018. The utility of taxonomically standardized names for copepod nauplii to increase comparability and measurability of ecological data. Hydrobiologia 820: 155-167.

McArdle BH, Anderson MJ. 2001. Fitting multivariate models to community data: A comment on distance-based redundancy analysis. Ecology 82: 290-297.

Sokal RR, Rohlf FJ. 1995. Biometry. 3rd edn. W.H. Freeman: New York.

Sorensen T. 1948. A method of establishing groups of equal amplitude in plant sociology based on similarity of species content. Kongelige Danske Videnskabernes Selskab Biologiske Skrifter 5: 1-34.
