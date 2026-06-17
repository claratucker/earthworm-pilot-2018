#!/usr/bin/env python3
"""Taxonomic composition (top genera) by treatment, gut vs soil.

Genus colors are consistent across both panels so compositional
differences are directly comparable. Treatment (Control vs Roundup) is
shown as a highlighted background box behind each x-axis tick label rather
than via bar color, since bar color encodes genus identity. Genera below
an explicit relative abundance threshold are pooled into "Other" and
labeled with that threshold.
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
        environment = 'Gut' if '.In' in sample else 'Soil'
        metadata_list.append({'sample': sample, 'treatment': treatment, 'environment': environment})
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

OTHER_THRESHOLD = 0.01  # 1% mean relative abundance across all samples

mean_abundance = feature_by_genus_rel.mean(axis=1)
top_genera = mean_abundance[mean_abundance >= OTHER_THRESHOLD].sort_values(ascending=False).index.tolist()

print(f"Genera at or above {OTHER_THRESHOLD*100:.0f}% mean relative abundance: {len(top_genera)}")
for i, genus in enumerate(top_genera, 1):
    print(f"  {i}. {genus} ({mean_abundance[genus]*100:.2f}%)")
print()

other_label = f'Other (<{OTHER_THRESHOLD*100:.0f}% mean abundance)'

composition = feature_by_genus_rel.loc[top_genera].copy().T
composition[other_label] = 1 - composition.sum(axis=1)
composition = composition.join(metadata)

genus_order = top_genera + [other_label]
qualitative_colors = sns.color_palette("tab20", len(top_genera))
genus_color_map = dict(zip(top_genera, qualitative_colors))
genus_color_map[other_label] = '#BBBBBB'

treatment_box_colors = {
    ('Gut', 'Control'): '#90EE90',
    ('Gut', 'Roundup'): '#0B3D0B',
    ('Soil', 'Control'): '#F5DEB3',
    ('Soil', 'Roundup'): '#654321'
}
treatment_text_colors = {
    ('Gut', 'Control'): 'black',
    ('Gut', 'Roundup'): 'white',
    ('Soil', 'Control'): 'black',
    ('Soil', 'Roundup'): 'white'
}

fig, axes = plt.subplots(1, 2, figsize=(18, 7))

for ax_idx, environment in enumerate(['Gut', 'Soil']):
    ax = axes[ax_idx]

    comp_data = composition[composition['environment'] == environment]

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
    ax.set_title(environment, fontsize=14, fontweight='bold')
    ax.set_xticks(range(len(comp_data)))
    ax.set_xticklabels(comp_data.index, rotation=45, ha='right', fontsize=9)
    ax.set_xlabel('')

    for tick_idx, sample in enumerate(comp_data.index):
        treatment = comp_data.loc[sample, 'treatment']
        box_color = treatment_box_colors[(environment, treatment)]
        text_color = treatment_text_colors[(environment, treatment)]
        label = ax.get_xticklabels()[tick_idx]
        label.set_color(text_color)
        label.set_fontweight('bold')
        label.set_bbox(dict(facecolor=box_color, edgecolor='none', boxstyle='round,pad=0.3'))

    n_control = len(comp_data[comp_data['treatment'] == 'Control'])
    if n_control > 0:
        ax.axvline(n_control - 0.5, color='black', linestyle='--', linewidth=1.5, alpha=0.5)

    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3, axis='y')

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center left', bbox_to_anchor=(0.93, 0.5),
          fontsize=9, framealpha=0.9, title='Genus')

fig.suptitle('Taxonomic Composition by Treatment', fontsize=15, fontweight='bold', y=1.02)

plt.tight_layout(rect=[0, 0, 0.89, 1])
plt.savefig('results/taxa_composition/taxa_composition_by_treatment.pdf', dpi=300, bbox_inches='tight')
plt.savefig('results/figures/figure_5_taxa_composition.pdf', dpi=300, bbox_inches='tight')
print("Saved: results/taxa_composition/taxa_composition_by_treatment.pdf")
print("Saved: results/figures/figure_5_taxa_composition.pdf")

composition.to_csv('results/taxa_composition/taxa_composition_by_sample.csv')
print("Saved: results/taxa_composition/taxa_composition_by_sample.csv")

print("\n" + "=" * 70)
print("MEAN RELATIVE ABUNDANCE BY TREATMENT")
print("=" * 70)

for environment in ['Gut', 'Soil']:
    print(f"\n{environment.upper()}")
    comp_data = composition[composition['environment'] == environment]
    for treatment in ['Control', 'Roundup']:
        treat_data = comp_data[comp_data['treatment'] == treatment]
        print(f"\n  {treatment} (n={len(treat_data)}):")
        means = treat_data[top_genera].mean().sort_values(ascending=False)
        for genus, mean_val in means.head(5).items():
            print(f"    {genus}: {mean_val:.4f}")
