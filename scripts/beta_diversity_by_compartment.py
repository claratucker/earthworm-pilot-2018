#!/usr/bin/env python3
"""Beta diversity NMDS stratified by compartment (Gut vs Soil separately)."""

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
        compartment = 'Gut' if '.In' in sample else 'Soil'
        replicate_str = sample.split('.')[0][1:]
        replicate = int(replicate_str)
        
        metadata_list.append({
            'sample': sample,
            'treatment': treatment,
            'compartment': compartment,
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


def plot_nmds_compartment(nmds_df, metadata, stress_value, compartment, output_file=None):
    """Create NMDS plot for single compartment."""
    
    plot_data = nmds_df.reset_index()
    plot_data = plot_data.merge(metadata.reset_index(), left_on='index', right_on='sample')
    plot_data = plot_data.rename(columns={'index': 'sample'})
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    colors = {'Control': '#DEEBF7', 'Roundup': '#08519C'}
    
    # Plot individual samples
    for treatment in ['Control', 'Roundup']:
        subset = plot_data[plot_data['treatment'] == treatment]
        ax.scatter(subset['NMDS1'], subset['NMDS2'],
                  c=colors[treatment], s=150, alpha=0.6,
                  edgecolors='black', linewidth=0.8, label=treatment)
    
    # Plot centroids and ellipses
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
    
    ax.set_title(f'NMDS: {compartment} Compartment - Bray-Curtis Distances\n(Stress = {stress_value:.3f}, {stress_interp})',
                fontsize=13, fontweight='bold')
    
    ax.legend(loc='best', fontsize=11, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle=':')
    
    n_samples = len(plot_data)
    ax.text(0.02, 0.98, f'{compartment} (n={n_samples})\nLight = Control, Dark = Roundup\nCentroids: + symbols\n95% ellipses: dashed lines',
           transform=ax.transAxes, verticalalignment='top', fontsize=9,
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved: {output_file}")
    
    return fig


def permanova_test(distance_matrix, metadata, treatment_column='treatment'):
    """Simplified PERMANOVA."""
    groups = metadata[treatment_column].values
    unique_groups = np.unique(groups)
    
    grand_mean = distance_matrix.values.mean()
    
    ss_between = 0
    for group in unique_groups:
        group_indices = np.where(groups == group)[0]
        group_distances = distance_matrix.values[np.ix_(group_indices, group_indices)]
        group_mean = group_distances.mean()
        ss_between += len(group_indices) * (group_mean - grand_mean) ** 2
    
    ss_total = np.sum((distance_matrix.values - grand_mean) ** 2)
    r_squared = ss_between / ss_total
    
    return {'r_squared': r_squared}


def run_compartment_analysis(feature_table_path, compartment, output_dir='results/beta_diversity_by_compartment'):
    """Run NMDS for single compartment."""
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 70)
    print(f"BETA DIVERSITY - {compartment.upper()} COMPARTMENT")
    print("=" * 70)
    print()
    
    print("Loading feature table...")
    feature_table = load_feature_table(feature_table_path)
    
    print("Parsing metadata...")
    metadata = parse_sample_metadata(feature_table.columns)
    
    # Filter to compartment
    compartment_samples = metadata[metadata['compartment'] == compartment].index.tolist()
    feature_table_subset = feature_table[compartment_samples]
    metadata_subset = metadata.loc[compartment_samples]
    
    print(f"Samples in {compartment}: {len(compartment_samples)}")
    print(f"Replicates per treatment:")
    print(metadata_subset.groupby('treatment').size())
    print()
    
    print("Computing Bray-Curtis distances...")
    distance_matrix = compute_distance_matrix(feature_table_subset)
    print(f"Distance matrix: {distance_matrix.shape}")
    print()
    
    print("Running NMDS (2D)...")
    nmds_scores, stress = nmds_ordination(distance_matrix, n_dimensions=2)
    print(f"Stress: {stress:.4f}")
    print()
    
    print("Running PERMANOVA...")
    perm_results = permanova_test(distance_matrix, metadata_subset, 'treatment')
    print(f"R² (treatment): {perm_results['r_squared']:.4f}")
    print()
    
    print("Saving results...")
    nmds_scores.to_csv(f'{output_dir}/{compartment.lower()}_nmds_scores.csv')
    distance_matrix.to_csv(f'{output_dir}/{compartment.lower()}_bray_curtis_distances.csv')
    metadata_subset.to_csv(f'{output_dir}/{compartment.lower()}_metadata.csv')
    
    print("Creating plot...")
    fig = plot_nmds_compartment(nmds_scores, metadata_subset, stress,
                               compartment, output_file=f'{output_dir}/{compartment.lower()}_nmds_ordination.pdf')
    
    caption = f"""NMDS ordination of {compartment} microbiome showing Bray-Curtis distances between control and Roundup-treated samples (stress = {stress:.3f}). Points colored by treatment (light blue = control, dark blue = Roundup). Centroids marked with '+' symbols. 95% confidence ellipses (dashed lines) show treatment group spread. n={len(compartment_samples)}."""
    
    with open(f'{output_dir}/{compartment.lower()}_figure_caption.txt', 'w') as f:
        f.write(caption)
    
    print()
    print("=" * 70)
    print(f"{compartment.upper()} RESULTS")
    print(f"  Stress: {stress:.4f}")
    print(f"  R² (treatment): {perm_results['r_squared']:.4f}")
    print(f"  n={len(compartment_samples)}")
    print("=" * 70)
    print()
    
    return nmds_scores, distance_matrix, metadata_subset, stress, perm_results


if __name__ == '__main__':
    
    gut_nmds, gut_dist, gut_meta, gut_stress, gut_perm = run_compartment_analysis(
        'r/exported/feature-table.tsv', 'Gut')
    
    soil_nmds, soil_dist, soil_meta, soil_stress, soil_perm = run_compartment_analysis(
        'r/exported/feature-table.tsv', 'Soil')
    
    print("\n" + "=" * 70)
    print("SUMMARY: TREATMENT EFFECT BY COMPARTMENT")
    print("=" * 70)
    print(f"Gut   (n=16): Stress={gut_stress:.3f}, R²={gut_perm['r_squared']:.4f}")
    print(f"Soil  (n=6):  Stress={soil_stress:.3f}, R²={soil_perm['r_squared']:.4f}")
    print("=" * 70)
