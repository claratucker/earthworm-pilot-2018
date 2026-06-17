#!/usr/bin/env python3
"""Combined NMDS (all 22 samples) testing Gut vs Soil separation, plus
PERMDISP test for within-group dispersion (heterogeneity) by compartment.

PERMANOVA implemented per McArdle & Anderson (2001), the standard
distance-based formulation building on Anderson (2001). Total sum of
squares is computed from squared pairwise distances directly
(SS_total = (1/n) * sum of all squared distances), and within-group sum
of squares is computed the same way restricted to each group, avoiding
the simple group-mean-distance approximation, which does not correctly
reflect centroid separation in distance space.

PERMDISP (Anderson 2006) reported using the ANOVA-on-distance-to-centroid
approximation (not full permutation PERMDISP), computed on an MDS
embedding of the same distance matrix.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS
from scipy.stats import f_oneway
import warnings
warnings.filterwarnings('ignore')

from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms


def parse_sample_metadata(sample_names):
    metadata_list = []
    for sample in sample_names:
        treatment = 'Control' if sample[0] == 'C' else 'Roundup'
        compartment = 'Gut' if '.In' in sample else 'Soil'
        metadata_list.append({'sample': sample, 'treatment': treatment, 'compartment': compartment})
    return pd.DataFrame(metadata_list).set_index('sample')


def load_feature_table(filepath):
    return pd.read_csv(filepath, sep='\t', index_col=0, skiprows=1)


def bray_curtis_distance(sample1, sample2):
    numerator = np.sum(np.abs(sample1 - sample2))
    denominator = np.sum(sample1 + sample2)
    return numerator / denominator if denominator > 0 else 0.0


def compute_distance_matrix(feature_table):
    samples = feature_table.columns.tolist()
    n_samples = len(samples)
    distance_matrix = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            dist = bray_curtis_distance(feature_table.iloc[:, i].values, feature_table.iloc[:, j].values)
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist
    return pd.DataFrame(distance_matrix, index=samples, columns=samples)


def nmds_ordination(distance_matrix, n_dimensions=2, random_state=42):
    mds = MDS(n_components=n_dimensions, dissimilarity='precomputed', random_state=random_state)
    scores = mds.fit_transform(distance_matrix)
    ordination_distances = pdist(scores)
    original_distances = pdist(distance_matrix, metric='euclidean')
    stress = np.sqrt(np.sum((original_distances - ordination_distances) ** 2) / np.sum(original_distances ** 2))
    nmds_df = pd.DataFrame(scores, index=distance_matrix.index, columns=[f'NMDS{i+1}' for i in range(n_dimensions)])
    return nmds_df, stress


def permanova_sum_of_squares(dist_sq_matrix, group_labels):
    """Compute PERMANOVA pseudo-F components per McArdle & Anderson (2001).

    dist_sq_matrix: square matrix of squared distances (n x n)
    group_labels: array of group assignments, length n
    Returns SS_within (sum across groups) and SS_total.
    """
    n = dist_sq_matrix.shape[0]
    ss_total = dist_sq_matrix.sum() / (2 * n)

    ss_within = 0
    unique_groups = np.unique(group_labels)
    for group in unique_groups:
        idx = np.where(group_labels == group)[0]
        n_i = len(idx)
        if n_i < 2:
            continue
        group_sq_dists = dist_sq_matrix[np.ix_(idx, idx)]
        ss_within += group_sq_dists.sum() / (2 * n_i)

    return ss_within, ss_total


def permanova_test(distance_matrix, metadata, group_column, n_permutations=999, random_state=42):
    """PERMANOVA (McArdle & Anderson 2001) with permutation-based p-value."""
    groups = metadata[group_column].values
    dist_sq_matrix = distance_matrix.values ** 2
    n = len(groups)
    unique_groups = np.unique(groups)
    a = len(unique_groups)

    ss_within_obs, ss_total = permanova_sum_of_squares(dist_sq_matrix, groups)
    ss_between_obs = ss_total - ss_within_obs
    r_squared = ss_between_obs / ss_total

    df_between = a - 1
    df_within = n - a
    pseudo_f_obs = (ss_between_obs / df_between) / (ss_within_obs / df_within)

    rng = np.random.RandomState(random_state)
    permuted_f = np.zeros(n_permutations)
    for i in range(n_permutations):
        permuted_labels = rng.permutation(groups)
        ss_within_perm, _ = permanova_sum_of_squares(dist_sq_matrix, permuted_labels)
        ss_between_perm = ss_total - ss_within_perm
        permuted_f[i] = (ss_between_perm / df_between) / (ss_within_perm / df_within)

    p_value = (np.sum(permuted_f >= pseudo_f_obs) + 1) / (n_permutations + 1)

    return {'r_squared': r_squared, 'pseudo_f': pseudo_f_obs, 'p_value': p_value, 'n_permutations': n_permutations}


def permdisp_test(distance_matrix, metadata, group_column):
    """PERMDISP (Anderson 2006), ANOVA-on-distance-to-centroid approximation."""
    groups = metadata[group_column].values
    unique_groups = np.unique(groups)

    mds = MDS(n_components=min(len(distance_matrix) - 1, 10),
              dissimilarity='precomputed', random_state=42)
    coords = mds.fit_transform(distance_matrix)

    distances_to_centroid = {}
    for group in unique_groups:
        idx = np.where(groups == group)[0]
        group_coords = coords[idx]
        centroid = group_coords.mean(axis=0)
        dists = np.sqrt(((group_coords - centroid) ** 2).sum(axis=1))
        distances_to_centroid[group] = dists

    f_stat, p_value = f_oneway(*distances_to_centroid.values())

    summary = {
        group: {'mean_dist_to_centroid': dists.mean(), 'sd': dists.std(), 'n': len(dists)}
        for group, dists in distances_to_centroid.items()
    }

    return {'f_statistic': f_stat, 'p_value': p_value, 'group_summary': summary}


def confidence_ellipse(x, y, n_std=2.0):
    cov = np.cov(x, y)
    if len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        return 0.1, 0.1, np.mean(x), np.mean(y), 0
    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    return scale_x, scale_y, mean_x, mean_y, 0


def plot_combined_nmds(nmds_df, metadata, stress_value, permanova_results, output_file):
    plot_data = nmds_df.reset_index().merge(metadata.reset_index(), left_on='index', right_on='sample')
    plot_data = plot_data.rename(columns={'index': 'sample'})

    fig, ax = plt.subplots(figsize=(12, 10))

    colors = {'Gut': '#2E8B57', 'Soil': '#8B5A3C'}

    for compartment in ['Gut', 'Soil']:
        subset = plot_data[plot_data['compartment'] == compartment]
        ax.scatter(subset['NMDS1'], subset['NMDS2'], c=colors[compartment], s=150,
                  alpha=0.6, edgecolors='black', linewidth=0.8, label=compartment)

        centroid_x, centroid_y = subset['NMDS1'].mean(), subset['NMDS2'].mean()
        ax.scatter(centroid_x, centroid_y, c=colors[compartment], marker='+', s=400,
                  linewidth=3, edgecolors='black', zorder=5)

        scale_x, scale_y, mean_x, mean_y, _ = confidence_ellipse(subset['NMDS1'].values, subset['NMDS2'].values)
        if scale_x > 0 and scale_y > 0:
            ellipse = Ellipse((0, 0), width=scale_x*2, height=scale_y*2,
                            facecolor=colors[compartment], alpha=0.15,
                            edgecolor=colors[compartment], linewidth=2, linestyle='--')
            transf = transforms.Affine2D().scale(scale_x, scale_y).translate(mean_x, mean_y)
            ellipse.set_transform(transf + ax.transData)
            ax.add_patch(ellipse)

    ax.set_xlabel('NMDS1', fontsize=12)
    ax.set_ylabel('NMDS2', fontsize=12)

    if stress_value < 0.1:
        stress_interp = "Good fit"
    elif stress_value < 0.2:
        stress_interp = "Acceptable"
    else:
        stress_interp = "Moderate distortion"

    ax.set_title(f'NMDS: All Samples, Bray-Curtis Distances\n(Stress = {stress_value:.3f}, {stress_interp})',
                fontsize=13, fontweight='bold')
    ax.legend(loc='best', fontsize=11, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle=':')

    stats_text = (f'All samples (n={len(plot_data)})\n'
                 f'R² (compartment) = {permanova_results["r_squared"]:.4f}\n'
                 f'pseudo-F = {permanova_results["pseudo_f"]:.4f}\n'
                 f'p = {permanova_results["p_value"]:.4f} '
                 f'({permanova_results["n_permutations"]} permutations)')
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, verticalalignment='top',
           fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")


if __name__ == '__main__':
    import os
    os.makedirs('results/beta_diversity_combined', exist_ok=True)

    print("Loading feature table...")
    feature_table = load_feature_table('r/exported/feature-table.tsv')
    metadata = parse_sample_metadata(feature_table.columns)

    print("Computing Bray-Curtis distances (all 22 samples)...")
    distance_matrix = compute_distance_matrix(feature_table)

    print("Running NMDS (2D)...")
    nmds_scores, stress = nmds_ordination(distance_matrix, n_dimensions=2)
    print(f"Stress: {stress:.4f}")
    print()

    print("Running PERMANOVA (compartment: Gut vs Soil)...")
    permanova_results = permanova_test(distance_matrix, metadata, 'compartment', n_permutations=999)
    print(f"R² = {permanova_results['r_squared']:.4f}")
    print(f"pseudo-F = {permanova_results['pseudo_f']:.4f}")
    print(f"p = {permanova_results['p_value']:.4f} ({permanova_results['n_permutations']} permutations)")
    print()

    print("Running PERMDISP (within-group dispersion by compartment)...")
    permdisp_results = permdisp_test(distance_matrix, metadata, 'compartment')
    print(f"F = {permdisp_results['f_statistic']:.4f}")
    print(f"p = {permdisp_results['p_value']:.4f}")
    print()
    print("Mean distance to group centroid (higher = more internal variability):")
    for group, stats in permdisp_results['group_summary'].items():
        print(f"  {group}: {stats['mean_dist_to_centroid']:.4f} ± {stats['sd']:.4f} (n={stats['n']})")
    print()

    print("Creating plot...")
    plot_combined_nmds(nmds_scores, metadata, stress, permanova_results,
                       'results/beta_diversity_combined/nmds_all_samples_compartment.pdf')
    plot_combined_nmds(nmds_scores, metadata, stress, permanova_results,
                       'results/figures/figure_6_nmds_all_samples_gut_vs_soil.pdf')

    nmds_scores.to_csv('results/beta_diversity_combined/nmds_scores.csv')
    distance_matrix.to_csv('results/beta_diversity_combined/bray_curtis_distances.csv')
    metadata.to_csv('results/beta_diversity_combined/metadata.csv')

    with open('results/beta_diversity_combined/permdisp_results.txt', 'w') as f:
        f.write("PERMDISP: Within-group dispersion by compartment\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"F statistic: {permdisp_results['f_statistic']:.4f}\n")
        f.write(f"p-value: {permdisp_results['p_value']:.4f}\n\n")
        f.write("Mean distance to group centroid:\n")
        for group, stats in permdisp_results['group_summary'].items():
            f.write(f"  {group}: {stats['mean_dist_to_centroid']:.4f} +/- {stats['sd']:.4f} (n={stats['n']})\n")
    print("Saved: results/beta_diversity_combined/permdisp_results.txt")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Gut vs Soil centroid separation (PERMANOVA): R²={permanova_results['r_squared']:.4f}, pseudo-F={permanova_results['pseudo_f']:.4f}, p={permanova_results['p_value']:.4f}")
    print(f"Gut vs Soil dispersion difference (PERMDISP): F={permdisp_results['f_statistic']:.4f}, p={permdisp_results['p_value']:.4f}")
    print("=" * 70)
