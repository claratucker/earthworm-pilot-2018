# 2018 Earthworm Glyphosate Pilot: Analysis Summary

## Project Overview
Glyphosate-exposed earthworms (Lumbricus terrestris) and their associated soil microbiome.
- Treatments: Control (C) vs. Roundup-treated (R)
- Environments: Earthworm gut (In) and surrounding soil (Soil)
- Samples: ~20 total (multiple replicates per treatment)

## Completed Steps

### 1. Metagenomic Processing
- 16S rRNA amplicon sequencing processed through QIIME2
- ASV table generated with alpha/beta diversity metrics

### 2. EPSPS Classification (v1.0.7 - FIXED)

**Problem solved:** Class I classifier was returning 0% despite real class I organisms present.

**Root cause:** 40% whole-protein identity gate (not in Leino et al.'s paper) was blocking marker checking before any threshold evaluation.

**Solution:** Dropped identity gate, implemented alignment-based markers for all four classes, report classes within 1.5% of best match (following Leino et al.'s actual method).

**Results on 164 genera:**
- Class I (Sensitive): 53 genera (32.1%)
- Class II (Resistant): 76 genera (46.1%)
- Class III (Resistant): 5 genera (3.0%)
- Class IV (Resistant): 18 genera (10.9%)
- Mixed (ambiguous): ~24 genera in practice (14.6%)
- Unclassified: 12 genera (7.3%)

**Validation:**
- 42 unit tests pass (including real benchmark organisms)
- Distribution matches literature: 32% Class I in soil vs 54% in Leino's gut bacteria (expected)

### 3. PICRUSt2 KO Predictions
- K00800 (aroA/EPSPS) extracted for all genera
- Full KO table available for functional analysis

## In Progress / Not Started

### Functional Gene Analysis — Not Started
- Genes to profile: nifH (N-fixation), narG/nirK (denitrification), ppx (P cycling), chitinase (C degradation)
- Approach: Extract KO abundances from PICRUSt2, stratify by EPSPS class and treatment

### Active-Site Residue Classifier — Not Started
- Alternative method using Gly96, Thr97, Pro101, Gly137, Ala183
- Will serve as independent validation of EPSPS classification

### Treatment-Level Statistics — Not Started
- Do resistant/sensitive fractions differ between control and glyphosate treatments?
- Does soil differ from gut microbiome?

## Key Results

**EPSPS Classification by Class:**
- 53 sensitive genera (Class I) — 32.1%
- 99 resistant genera (II+III+IV) — 60.0%
- 12 unclassified — 7.3%
- 24 genuinely ambiguous (mixed class) — in 1.5% margin

**Confidence Gap (Identity I - Identity II):**
- Mean: -0.66 (Class II slightly higher on average)
- Median: -0.16 (close to zero; many borderline calls)
- Range: -39.59 to +51.19 (some very clean, some very ambiguous)

## Key Files

**EPSPS Results:**
- ~/pilot2018/results/epsps/epsps_classified.tsv (164 genera, full classification)
- ~/pilot2018/results/epsps/epsps_summary.txt (summary counts)

**Scripts:**
- ~/pilot2018/scripts/epsps/04_epsps_exploratory.sh (pipeline)
- ~/pilot2018/scripts/r/06_epsps_summary_by_treatment.R (summary stats)

**Code Repository:**
- ~/epspsclass (v1.0.7, GitHub: claratucker/epspsclass)

## Next Priority

Attach treatment metadata to EPSPS results and test whether glyphosate treatment shifts the sensitive/resistant balance in your samples.
