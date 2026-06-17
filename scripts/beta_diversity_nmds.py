#!/usr/bin/env python3
"""Beta diversity NMDS for earthworm pilot (2x2 design)."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist
from sklearn.manifold import MDS
import warnings
warnings.filterwarnings('ignore')

from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms


def parse_sample_metadata(sample_names):
    """Parse metadata from sample names: C1.In, R3.Soil, etc."""
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
    """Bray-Curtis dissimilarity. Literature: Sorensen (1948)."""
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
    """NMDS with stress. Stress < 0.1 good, 0.1-0.2 acceptable, > 0.3 bad."""
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
    """95% confidence ellipse. Literature: Sokal & Rohlf (1995)."""
    cov = np.cov(x, y)
    if len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        return 0.1, 0.1, np.mean(x), np.mean(y), 0
    
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    
    return scale_x, scale_y, mean_x, mean_y, pearson


def plot_nmds_with_features(nmds_df, metadata, distance_matrix, stress_value,
                            figsize=(14, 10), output_file=None, n_dims=2):
    """Create NMDS plot with Pochron requirements."""
    
    plot_data = nmds_df.reset_index()
    plot_data = plot_data.merge(metadata.reset_index(), left_on='index', right_on='sample')
    plot_data = plot_data.rename(columns={'index': 'sample'})
    
    fig, ax = plt.subplots(figsize=figsize)
    
    colors = {'Control': '#DEEBF7', 'Roundup': '#08519C'}
    markers = {'Gut': 'o', 'Soil': 's'}
    
    for (treatment, compartment), group in plot_data.groupby(['treatment', 'compartment']):
        ax.scatter(group['NMDS1'], group['NMDS2'],
                  c=colors[treatment], marker=markers[compartment],
                  s=150, alpha=0.6, edgecolors='black', linewidth=0.8,
                  label=f'{treatment} - {compartment}')
    
    for treatment in ['Control', 'Roundup']:
        for compartment in ['Gut', 'Soil']:
            subset = plot_data[(plot_data['treatment'] == treatment) &
                             (plot_data['compartment'] == compartment)]
            
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
        stress_interp = "Unsatisfactory"
    
    ax.set_title(f'NMDS Ordination ({n_dims}D) - Bray-Curtis Distances\n(Stress = {stress_value:.3f}, {stress_interp})',
                fontsize=13, fontweight='bold')
    
    ax.legend(loc='best', fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle=':')
    
    ax.text(0.02, 0.98, f'Stress: {stress_interp}\nLight = Control, Dark = Roundup\nCircle = Gut, Square = Soil',
           transform=ax.transAxes, verticalalignment='top', fontsize=9,
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved: {output_file}")
    
    return fig


def permanova_test(distance_matrix, metadata, treatment_column='treatment'):
    """Simplified PERMANOVA. Literature: Anderson (2001)."""
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


def run_analysis(feature_table_path, output_dir='results/beta_diversity', n_dims=2):
    """Run complete beta diversity analysis pipeline."""
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 70)
    print("EARTHWORM PILOT STUDY - BETA DIVERSITY ANALYSIS")
    print("=" * 70)
    print()
    
    print("Loading feature table...")
    feature_table = load_feature_table(feature_table_path)
    print(f"Feature table: {feature_table.shape[0]} ASVs × {feature_table.shape[1]} samples")
    print()
    
    print("Parsing sample metadata from names...")
    metadata = parse_sample_metadata(feature_table.columns)
    print(f"Treatment: {metadata['treatment'].unique()}")
    print(f"Compartments: {metadata['compartment'].unique()}")
    print(f"Replicates per group:")
    print(metadata.groupby(['treatment', 'compartment']).size())
    print()
    
    print("Computing Bray-Curtis distances...")
    distance_matrix = compute_distance_matrix(feature_table)
    print(f"Distance matrix: {distance_matrix.shape}")
    print()
    
    print(f"Running NMDS ordination ({n_dims}D)...")
    nmds_scores, stress = nmds_ordination(distance_matrix, n_dimensions=n_dims)
    print(f"NMDS stress value: {stress:.4f}")
    if stress < 0.1:
        print("Interpretation: Good fit")
    elif stress < 0.2:
        print("Interpretation: Acceptable")
    else:
        print("Interpretation: Unsatisfactory (recommend 3D ordination)")
    print()
    
    print("Running PERMANOVA (simplified)...")
    perm_results = permanova_test(distance_matrix, metadata, 'treatment')
    print(f"R-squared (treatment): {perm_results['r_squared']:.4f}")
    print()
    
    print("Saving results...")
    nmds_scores.to_csv(f'{output_dir}/nmds_scores_{n_dims}d.csv')
    distance_matrix.to_csv(f'{output_dir}/bray_curtis_distances.csv')
    metadata.to_csv(f'{output_dir}/metadata.csv')
    
    print("Creating NMDS plot...")
    fig = plot_nmds_with_features(nmds_scores, metadata, distance_matrix, stress,
                                 output_file=f'{output_dir}/nmds_ordination_{n_dims}d.pdf', n_dims=n_dims)
    
    caption = f"""NMDS ordination based on Bray-Curtis distances showing microbial community composition differences between control and Roundup-exposed earthworms in gut and soil compartments (stress = {stress:.3f}). Points represent individual samples, colored by treatment (light blue = control, dark blue = Roundup) and shaped by compartment (circles = gut, squares = soil). Treatment centroids are marked with '+' symbols. 95% confidence ellipses (dashed lines) enclose each treatment-compartment group. Stress value {stress:.3f} indicates {'good' if stress < 0.1 else 'acceptable' if stress < 0.2 else 'unsatisfactory'} fit of dissimilarities to {n_dims}-dimensional ordination space."""
    
    with open(f'{output_dir}/nmds_figure_caption_{n_dims}d.txt', 'w') as f:
        f.write(caption)
    
    print()
    print("Figure caption:")
    print(caption)
    print()
    print("=" * 70)
    print(f"Analysis complete. Output files in {output_dir}:")
    print(f"  - nmds_ordination_{n_dims}d.pdf")
    print(f"  - nmds_scores_{n_dims}d.csv")
    print(f"  - bray_curtis_distances.csv")
    print(f"  - metadata.csv")
    print(f"  - nmds_figure_caption_{n_dims}d.txt")
    print("=" * 70)
    
    return nmds_scores, distance_matrix, metadata, stress


if __name__ == '__main__':
    print("\n*** Running 2D NMDS (stress likely high) ***\n")
    nmds_2d, dist_matrix, meta, stress_2d = run_analysis('r/exported/feature-table.tsv', n_dims=2)
    
    if stress_2d > 0.3:
        print("\n\n*** Stress unsatisfactory. Running 3D NMDS for better fit ***\n")
        nmds_3d, _, _, stress_3d = run_analysis('r/exported/feature-table.tsv', n_dims=3)
