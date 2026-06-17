#!/usr/bin/env python3
"""Taxonomic composition (top genera) by treatment and compartment.

Genus colors are consistent across Gut and Soil panels so compositional
differences between compartments are directly comparable. Treatment
(Control vs Roundup) is indicated by x-axis tick label color, not bar color,
since bar color encodes genus identity.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

print("Loading data...")
feature_table = pd.read_csv('r/exported/feature-table.tsv', sep='\t', index_col=0, skiprows=1)
taxonomy = pd.read_csv('r/exported/taxonomy.tsv', sep='\t', index_col=0)

def parse_metadata(sample_names):
    metadata_list = []
    for sample in sample_names:
        treatment = 'Control' if sample[0] == 'C' else 'Roundup'
        compartment = 'Gut' if '.In' in sample else 'Soil'
        metadata_list.append({'sample': sample, 'treatment': treatment, 'compartment': compartment})
    return pd.DataFrame(metadata_list).set_index('sample')

metadata = parse_metadata(feature_table.columns)

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

feature_by_genus = pd.DataFrame(index=taxonomy['genus'].unique())
for asv_id in feature_table.index:
    genus = taxonomy.loc[asv_id, 'genus']
    feature_by_genus.loc[genus, feature_table.columns] = feature_table.loc[asv_id].values
feature_by_genus = feature_by_genus.fillna(0)

feature_by_genus_rel = feature_by_genus.div(feature_by_genus.sum(axis=0), axis=1)

top_n = 15
top_genera = feature_by_genus_rel.sum(axis=1).nlargest(top_n).index.tolist()

print(f"Top {top_n} genera by total abundance:")
for i, genus in enumerate(top_genera, 1):
    print(f"  {i}. {genus}")
print()

composition = feature_by_genus_rel.loc[top_genera].copy().T
composition['Other'] = 1 - composition.sum(axis=1)
composition = composition.join(metadata)

# Single fixed color map for genera, shared across both panels.
# Uses a qualitative palette since genus identity, not magnitude, is what
# needs to be visually distinguished and matched between panels.
genus_order = top_genera + ['Other']
qualitative_colors = sns.color_palette("tab20", len(genus_order))
genus_color_map = dict(zip(genus_order, qualitative_colors))
genus_color_map['Other'] = '#BBBBBB'

# Treatment label colors (matches established scheme: green=Gut, brown=Soil,
# light=Control, dark=Roundup). Applied to tick labels only, not bars.
treatment_label_colors = {
    ('Gut', 'Control'): '#3CB043',
    ('Gut', 'Roundup'): '#0B3D0B',
    ('Soil', 'Control'): '#C9A876',
    ('Soil', 'Roundup'): '#654321'
}

fig, axes = plt.subplots(1, 2, figsize=(17, 7))

for ax_idx, compartment in enumerate(['Gut', 'Soil']):
    ax = axes[ax_idx]

    comp_data = composition[composition['compartment'] == compartment]

    sample_order = []
    for treatment in ['Control', 'Roundup']:
        treat_samples = comp_data[comp_data['treatment'] == treatment].index.tolist()
        sample_order.extend(sorted(treat_samples))
    comp_data = comp_data.loc[sample_order]

    bottom = np.zeros(len(comp_data))
    for genus in genus_order:
        values = comp_data[genus].values
        ax.bar(range(len(comp_data)), values, bottom=bottom,
              label=genus, color=genus_color_map[genus], edgecolor='white', linewidth=1)
        bottom += values

    ax.set_ylabel('Relative Abundance', fontsize=12)
    ax.set_xlabel('Sample', fontsize=12)
    ax.set_title(f'Taxonomic Composition, {compartment} Compartment', fontsize=13, fontweight='bold')
    ax.set_xticks(range(len(comp_data)))
    ax.set_xticklabels(comp_data.index, rotation=45, ha='right', fontsize=9)

    for tick_idx, sample in enumerate(comp_data.index):
        treatment = comp_data.loc[sample, 'treatment']
        ax.get_xticklabels()[tick_idx].set_color(treatment_label_colors[(compartment, treatment)])
        ax.get_xticklabels()[tick_idx].set_fontweight('bold')

    n_control = len(comp_data[comp_data['treatment'] == 'Control'])
    if n_control > 0:
        ax.axvline(n_control - 0.5, color='black', linestyle='--', linewidth=1.5, alpha=0.5)

    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3, axis='y')

# One shared legend for genus identity (same colors both panels)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1.0, 0.5),
          fontsize=9, framealpha=0.9, title='Genus')

# Caption-style note explaining label colors, placed below the plots
fig.text(0.5, -0.02,
        'X-axis sample labels colored by treatment: green = Gut Control/Roundup, brown = Soil Control/Roundup (light = Control, dark = Roundup).',
        ha='center', fontsize=9, style='italic')

plt.tight_layout(rect=[0, 0.02, 0.88, 1])
plt.savefig('results/taxa_composition/taxa_composition_by_treatment.pdf', dpi=300, bbox_inches='tight')
plt.savefig('results/figures/figure_5_taxa_composition.pdf', dpi=300, bbox_inches='tight')
print("Saved: results/taxa_composition/taxa_composition_by_treatment.pdf")
print("Saved: results/figures/figure_5_taxa_composition.pdf")

composition.to_csv('results/taxa_composition/taxa_composition_by_sample.csv')
print("Saved: results/taxa_composition/taxa_composition_by_sample.csv")

print("\n" + "=" * 70)
print("MEAN RELATIVE ABUNDANCE BY TREATMENT AND COMPARTMENT")
print("=" * 70)

for compartment in ['Gut', 'Soil']:
    print(f"\n{compartment.upper()}")
    comp_data = composition[composition['compartment'] == compartment]
    for treatment in ['Control', 'Roundup']:
        treat_data = comp_data[comp_data['treatment'] == treatment]
        print(f"\n  {treatment} (n={len(treat_data)}):")
        means = treat_data[top_genera].mean().sort_values(ascending=False)
        for genus, mean_val in means.head(5).items():
            print(f"    {genus}: {mean_val:.4f}")
