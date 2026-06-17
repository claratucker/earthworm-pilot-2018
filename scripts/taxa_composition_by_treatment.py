#!/usr/bin/env python3
"""Taxonomic composition (top genera) by treatment and compartment."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

print("Loading data...")
feature_table = pd.read_csv('r/exported/feature-table.tsv', sep='\t', index_col=0, skiprows=1)
taxonomy = pd.read_csv('r/exported/taxonomy.tsv', sep='\t', index_col=0)

# Parse metadata
def parse_metadata(sample_names):
    metadata_list = []
    for sample in sample_names:
        treatment = 'Control' if sample[0] == 'C' else 'Roundup'
        compartment = 'Gut' if '.In' in sample else 'Soil'
        metadata_list.append({'sample': sample, 'treatment': treatment, 'compartment': compartment})
    return pd.DataFrame(metadata_list).set_index('sample')

metadata = parse_metadata(feature_table.columns)

# Extract genus from taxonomy
def extract_genus(taxon_string):
    if pd.isna(taxon_string):
        return 'Unclassified'
    parts = taxon_string.split(';')
    for part in parts:
        if part.startswith('g__'):
            genus = part.replace('g__', '').strip()
            return genus if genus else 'Unclassified'
    return 'Unclassified'

taxonomy['genus'] = taxonomy['Taxon'].apply(extract_genus)

# Aggregate feature table by genus
feature_by_genus = pd.DataFrame(index=taxonomy['genus'].unique())

for asv_id in feature_table.index:
    genus = taxonomy.loc[asv_id, 'genus']
    feature_by_genus.loc[genus, feature_table.columns] = feature_table.loc[asv_id].values

feature_by_genus = feature_by_genus.fillna(0)

# Normalize to relative abundance
feature_by_genus_rel = feature_by_genus.div(feature_by_genus.sum(axis=0), axis=1)

# Get top genera
top_n = 15
top_genera = feature_by_genus_rel.sum(axis=1).nlargest(top_n).index.tolist()

print(f"Top {top_n} genera by total abundance:")
for i, genus in enumerate(top_genera, 1):
    print(f"  {i}. {genus}")
print()

# Create composition matrix
composition = feature_by_genus_rel.loc[top_genera].copy()
composition = composition.T
composition['Other'] = 1 - composition.sum(axis=1)

# Add metadata
composition = composition.join(metadata)

# Plot by compartment and treatment
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

color_schemes = {
    'Gut': sns.color_palette("Greens", len(top_genera) + 1),
    'Soil': sns.color_palette("YlOrBr", len(top_genera) + 1)
}

for ax_idx, compartment in enumerate(['Gut', 'Soil']):
    ax = axes[ax_idx]
    
    comp_data = composition[composition['compartment'] == compartment]
    comp_data = comp_data.sort_values('treatment')
    
    # Order: Control first, Roundup second
    x_positions = []
    x_labels = []
    group_colors = color_schemes[compartment]
    
    sample_order = []
    for treatment in ['Control', 'Roundup']:
        treat_samples = comp_data[comp_data['treatment'] == treatment].index.tolist()
        sample_order.extend(sorted(treat_samples))
    
    comp_data = comp_data.loc[sample_order]
    
    # Stack bar plot
    bottom = np.zeros(len(comp_data))
    
    for genus_idx, genus in enumerate(top_genera + ['Other']):
        values = comp_data[genus].values
        ax.bar(range(len(comp_data)), values, bottom=bottom,
              label=genus, color=group_colors[genus_idx], edgecolor='white', linewidth=1)
        bottom += values
    
    ax.set_ylabel('Relative Abundance', fontsize=12)
    ax.set_xlabel('Sample', fontsize=12)
    ax.set_title(f'Taxonomic Composition - {compartment} Compartment', fontsize=13, fontweight='bold')
    ax.set_xticks(range(len(comp_data)))
    
    # Label x-axis with treatment info
    x_labels_formatted = [f"{sample}\n({comp_data.loc[sample, 'treatment'][0]})" 
                         for sample in comp_data.index]
    ax.set_xticklabels(x_labels_formatted, rotation=45, ha='right', fontsize=9)
    
    # Add treatment separator
    n_control = len(comp_data[comp_data['treatment'] == 'Control'])
    if n_control > 0:
        ax.axvline(n_control - 0.5, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
    
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3, axis='y')

# Shared legend
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center', bbox_to_anchor=(0.5, -0.05), ncol=6, fontsize=9, framealpha=0.9)

plt.tight_layout(rect=[0, 0.08, 1, 1])
plt.savefig('results/figures/figure_5_taxa_composition.pdf', dpi=300, bbox_inches='tight')
print("Saved: results/figures/figure_5_taxa_composition.pdf")

# Also save composition data
composition.to_csv('results/taxa_composition/taxa_composition_by_sample.csv')
print("Saved: results/taxa_composition/taxa_composition_by_sample.csv")

# Summary by treatment and compartment
print("\n" + "=" * 70)
print("MEAN RELATIVE ABUNDANCE BY TREATMENT AND COMPARTMENT")
print("=" * 70)
print()

for compartment in ['Gut', 'Soil']:
    print(f"\n{compartment.upper()}")
    comp_data = composition[composition['compartment'] == compartment]
    
    for treatment in ['Control', 'Roundup']:
        treat_data = comp_data[comp_data['treatment'] == treatment]
        print(f"\n  {treatment} (n={len(treat_data)}):")
        
        means = treat_data[top_genera].mean()
        means = means.sort_values(ascending=False)
        
        for genus, mean_val in means.head(5).items():
            print(f"    {genus}: {mean_val:.4f}")

