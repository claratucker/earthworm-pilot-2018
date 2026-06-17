#!/usr/bin/env python3
"""Beta diversity NMDS stratified by environment (Gut vs Soil separately).

PERMANOVA implemented per McArdle & Anderson (2001), the standard
distance-based formulation. Total and within-group sums of squares are
computed from squared pairwise distances directly, with a permutation
test (999 permutations by default) used to obtain a p-value. This
replaces an earlier implementation that approximated PERMANOVA using
raw group-mean distances rather than the correct sum-of-squares
decomposition, and that did not compute a p-value at all.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from sklearn.manifold import MDS
import warnings
warnings.filterwarnings('ignore')

from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms


def parse_sample_metadata(sample_names):
    """Parse metadata from sample names."""
    metadata_list = []
    for sample in sample_names:
        treatment = 'Control' if sample[0] == 'C' else 'Roundup'
        environment = 'Gut' if '.In' in sample else 'Soil'
        replicate_str = sample.split('.')[0][1:]
        replicate = int(replicate_str)

        metadata_list.append({
            'sample': sample,
            'treatment': treatment,
            'environment': environment,
            'replicate': replicate
        })

    return pd.DataFrame(metadata_list).set_index('sample')


def load_feature_table(filepath):
    """Load QIIME2 exported feature table TSV."""
    feature_table = pd.read_csv(filepath, sep='\t', index_col=0, skiprows=1)
    return feature_table


def bray_curtis_distance(sample1, sample2):
    """Bray-Curtis dissimilarity."""
    numerator = np.sum(np.abs(sample1 - sample2))
    denominator = np.sum(sample1 + sample2)
    return numerator / denominator if denominator > 0 else 0.0


def compute_distance_matrix(feature_table):
    """Compute pairwise Bray-Curtis distances."""
    samples = feature_table.columns.tolist()
    n_samples = len(samples)
    distance_matrix = np.zeros((n_samples, n_samples))

    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            sample_i = feature_table.iloc[:, i].values
            sample_j = feature_table.iloc[:, j].values
            dist = bray_curtis_distance(sample_i, sample_j)
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist

    return pd.DataFrame(distance_matrix, index=samples, columns=samples)


def nmds_ordination(distance_matrix, n_dimensions=2, random_state=42):
    """NMDS with stress."""
    mds = MDS(n_components=n_dimensions, dissimilarity='precomputed',
              random_state=random_state)
    scores = mds.fit_transform(distance_matrix)

    ordination_distances = pdist(scores)
    original_distances = pdist(distance_matrix, metric='euclidean')

    stress = np.sqrt(np.sum((original_distances - ordination_distances) ** 2) /
                    np.sum(original_distances ** 2))

    nmds_df = pd.DataFrame(scores, index=distance_matrix.index,
                          columns=[f'NMDS{i+1}' for i in range(n_dimensions)])

    return nmds_df, stress


def confidence_ellipse(x, y, n_std=2.0):
    """95% confidence ellipse."""
    cov = np.cov(x, y)
    if len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        return 0.1, 0.1, np.mean(x), np.mean(y), 0

    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_x = np.mean(x)
    mean_y = np.mean(y)

    return scale_x, scale_y, mean_x, mean_y, pearson


def permanova_sum_of_squares(dist_sq_matrix, group_labels):
    """Compute PERMANOVA SS components per McArdle & Anderson (2001).

    dist_sq_matrix: square matrix of squared distances (n x n)
    group_labels: array of group assignments, length n
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


def permanova_test(distance_matrix, metadata, treatment_column='treatment',
                   n_permutations=999, random_state=42):
    """PERMANOVA (McArdle & Anderson 2001) with permutation-based p-value."""
    groups = metadata[treatment_column].values
    dist_sq_matrix = distance_matrix.values ** 2
    n = len(groups)
    unique_groups = np.unique(groups)
    a = len(unique_groups)

    ss_within_obs, ss_total = permanova_sum_of_squares(dist_sq_matrix, groups)
    ss_between_obs = ss_total - ss_within_obs
    r_squared = ss_between_obs / ss_total if ss_total > 0 else 0.0

    df_between = a - 1
    df_within = n - a

    if df_within <= 0 or ss_within_obs == 0:
        return {'r_squared': r_squared, 'pseudo_f': np.nan, 'p_value': np.nan,
               'n_permutations': n_permutations}

    pseudo_f_obs = (ss_between_obs / df_between) / (ss_within_obs / df_within)

    rng = np.random.RandomState(random_state)
    permuted_f = np.zeros(n_permutations)
    for i in range(n_permutations):
        permuted_labels = rng.permutation(groups)
        ss_within_perm, _ = permanova_sum_of_squares(dist_sq_matrix, permuted_labels)
        ss_between_perm = ss_total - ss_within_perm
        permuted_f[i] = (ss_between_perm / df_between) / (ss_within_perm / df_within)

    p_value = (np.sum(permuted_f >= pseudo_f_obs) + 1) / (n_permutations + 1)

    return {'r_squared': r_squared, 'pseudo_f': pseudo_f_obs, 'p_value': p_value,
           'n_permutations': n_permutations}


def plot_nmds_environment(nmds_df, metadata, stress_value, environment, permanova_results, output_file=None):
    """Create NMDS plot for single environment with stats."""

    plot_data = nmds_df.reset_index()
    plot_data = plot_data.merge(metadata.reset_index(), left_on='index', right_on='sample')
    plot_data = plot_data.rename(columns={'index': 'sample'})

    fig, ax = plt.subplots(figsize=(12, 10))

    if environment == 'Gut':
        colors = {'Control': '#90EE90', 'Roundup': '#0B3D0B'}
    else:
        colors = {'Control': '#F5DEB3', 'Roundup': '#654321'}

    for treatment in ['Control', 'Roundup']:
        subset = plot_data[plot_data['treatment'] == treatment]
        ax.scatter(subset['NMDS1'], subset['NMDS2'],
                  c=colors[treatment], s=150, alpha=0.6,
                  edgecolors='black', linewidth=0.8, label=treatment)

    for treatment in ['Control', 'Roundup']:
        subset = plot_data[plot_data['treatment'] == treatment]

        if len(subset) < 2:
            continue

        centroid_x = subset['NMDS1'].mean()
        centroid_y = subset['NMDS2'].mean()

        ax.scatter(centroid_x, centroid_y, c=colors[treatment],
                  marker='+', s=400, linewidth=3, edgecolors='black', zorder=5)

        scale_x, scale_y, mean_x, mean_y, pearson = confidence_ellipse(
            subset['NMDS1'].values, subset['NMDS2'].values, n_std=2.0)

        if scale_x > 0 and scale_y > 0:
            ellipse = Ellipse((0, 0), width=scale_x*2, height=scale_y*2,
                            facecolor=colors[treatment], alpha=0.15,
                            edgecolor=colors[treatment], linewidth=2, linestyle='--')

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

    ax.set_title(f'NMDS: {environment} Environment - Bray-Curtis Distances\n(Stress = {stress_value:.3f}, {stress_interp})',
                fontsize=13, fontweight='bold')

    ax.legend(loc='best', fontsize=11, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle=':')

    n_samples = len(plot_data)
    f_str = f"{permanova_results['pseudo_f']:.4f}" if not np.isnan(permanova_results['pseudo_f']) else "NA"
    p_str = f"{permanova_results['p_value']:.4f}" if not np.isnan(permanova_results['p_value']) else "NA"
    stats_text = (f'{environment} (n={n_samples})\n'
                 f"R² (treatment) = {permanova_results['r_squared']:.4f}\n"
                 f"pseudo-F = {f_str}\n"
                 f"p = {p_str} ({permanova_results['n_permutations']} permutations)")
    ax.text(0.02, 0.98, stats_text,
           transform=ax.transAxes, verticalalignment='top', fontsize=10,
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved: {output_file}")

    return fig


def run_environment_analysis(feature_table_path, environment, output_dir='results/beta_diversity_by_environment'):
    """Run NMDS for single environment."""
    import os
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 70)
    print(f"BETA DIVERSITY - {environment.upper()} ENVIRONMENT")
    print("=" * 70)
    print()

    print("Loading feature table...")
    feature_table = load_feature_table(feature_table_path)

    print("Parsing metadata...")
    metadata = parse_sample_metadata(feature_table.columns)

    environment_samples = metadata[metadata['environment'] == environment].index.tolist()
    feature_table_subset = feature_table[environment_samples]
    metadata_subset = metadata.loc[environment_samples]

    print(f"Samples in {environment}: {len(environment_samples)}")
    print(f"Replicates per treatment:")
    print(metadata_subset.groupby('treatment').size())
    print()

    print("Computing Bray-Curtis distances...")
    distance_matrix = compute_distance_matrix(feature_table_subset)
    print()

    print("Running NMDS (2D)...")
    nmds_scores, stress = nmds_ordination(distance_matrix, n_dimensions=2)
    print(f"Stress: {stress:.4f}")
    print()

    print("Running PERMANOVA (999 permutations)...")
    perm_results = permanova_test(distance_matrix, metadata_subset, 'treatment', n_permutations=999)
    print(f"R² (treatment): {perm_results['r_squared']:.4f}")
    print(f"pseudo-F: {perm_results['pseudo_f']:.4f}")
    print(f"p-value: {perm_results['p_value']:.4f}")
    print()

    print("Saving results...")
    nmds_scores.to_csv(f'{output_dir}/{environment.lower()}_nmds_scores.csv')
    distance_matrix.to_csv(f'{output_dir}/{environment.lower()}_bray_curtis_distances.csv')
    metadata_subset.to_csv(f'{output_dir}/{environment.lower()}_metadata.csv')

    print("Creating plot...")
    fig = plot_nmds_environment(nmds_scores, metadata_subset, stress, environment,
                               perm_results,
                               output_file=f'{output_dir}/{environment.lower()}_nmds_ordination.pdf')

    caption = (f"NMDS ordination of {environment} microbiome showing Bray-Curtis distances "
              f"between control and Roundup-treated samples (stress = {stress:.3f}). "
              f"Points colored by treatment. Centroids marked with '+' symbols. "
              f"95% confidence ellipses (dashed lines) show treatment group spread. "
              f"R²={perm_results['r_squared']:.4f}, pseudo-F={perm_results['pseudo_f']:.4f}, "
              f"p={perm_results['p_value']:.4f} ({perm_results['n_permutations']} permutations). "
              f"n={len(environment_samples)}.")

    with open(f'{output_dir}/{environment.lower()}_figure_caption.txt', 'w') as f:
        f.write(caption)

    print()
    print("=" * 70)
    print(f"{environment.upper()} RESULTS")
    print(f"  Stress: {stress:.4f}")
    print(f"  R² (treatment): {perm_results['r_squared']:.4f}")
    print(f"  pseudo-F: {perm_results['pseudo_f']:.4f}")
    print(f"  p-value: {perm_results['p_value']:.4f}")
    print(f"  n={len(environment_samples)}")
    print("=" * 70)
    print()

    return nmds_scores, distance_matrix, metadata_subset, stress, perm_results


if __name__ == '__main__':

    gut_nmds, gut_dist, gut_meta, gut_stress, gut_perm = run_environment_analysis(
        'r/exported/feature-table.tsv', 'Gut')

    soil_nmds, soil_dist, soil_meta, soil_stress, soil_perm = run_environment_analysis(
        'r/exported/feature-table.tsv', 'Soil')

    print("\n" + "=" * 70)
    print("SUMMARY: TREATMENT EFFECT BY ENVIRONMENT (corrected PERMANOVA)")
    print("=" * 70)
    print(f"Gut   (n=16): Stress={gut_stress:.3f}, R²={gut_perm['r_squared']:.4f}, pseudo-F={gut_perm['pseudo_f']:.4f}, p={gut_perm['p_value']:.4f}")
    print(f"Soil  (n=6):  Stress={soil_stress:.3f}, R²={soil_perm['r_squared']:.4f}, pseudo-F={soil_perm['pseudo_f']:.4f}, p={soil_perm['p_value']:.4f}")
    print("=" * 70)
