#!/usr/bin/env python3
"""Sanity check: do raw Bray-Curtis distances actually separate Gut from Soil?

Compares within-group distances (Gut-Gut, Soil-Soil) against between-group
distances (Gut-Soil) directly from the distance matrix, with no embedding,
no R-squared formula, just the raw numbers. If environment is a real
driver of composition, between-group distances should be clearly larger
than within-group distances.
"""

import pandas as pd
import numpy as np

distance_matrix = pd.read_csv('results/beta_diversity_combined/bray_curtis_distances.csv', index_col=0)
metadata = pd.read_csv('results/beta_diversity_combined/metadata.csv', index_col=0)

samples = distance_matrix.index.tolist()
gut_samples = metadata[metadata['environment'] == 'Gut'].index.tolist()
soil_samples = metadata[metadata['environment'] == 'Soil'].index.tolist()

within_gut = []
within_soil = []
between = []

for i in range(len(samples)):
    for j in range(i + 1, len(samples)):
        s1, s2 = samples[i], samples[j]
        d = distance_matrix.loc[s1, s2]
        if s1 in gut_samples and s2 in gut_samples:
            within_gut.append(d)
        elif s1 in soil_samples and s2 in soil_samples:
            within_soil.append(d)
        else:
            between.append(d)

within_gut = np.array(within_gut)
within_soil = np.array(within_soil)
between = np.array(between)

print("Raw Bray-Curtis distance comparison")
print("=" * 50)
print(f"Within-Gut  (n={len(within_gut)} pairs):  mean={within_gut.mean():.4f}, sd={within_gut.std():.4f}")
print(f"Within-Soil (n={len(within_soil)} pairs):  mean={within_soil.mean():.4f}, sd={within_soil.std():.4f}")
print(f"Between Gut-Soil (n={len(between)} pairs): mean={between.mean():.4f}, sd={between.std():.4f}")
print()

from scipy.stats import mannwhitneyu
u_stat, p_val = mannwhitneyu(np.concatenate([within_gut, within_soil]), between)
print(f"Within (Gut+Soil combined) vs Between: Mann-Whitney U={u_stat:.1f}, p={p_val:.6f}")
