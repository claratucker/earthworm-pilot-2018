#!/usr/bin/env python3
"""Analyze EPSPS Class I/II balance by treatment and compartment."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

print("Loading data...")
feature_table = pd.read_csv('r/exported/feature-table.tsv', sep='\t', index_col=0, skiprows=1)
taxonomy = pd.read_csv('r/exported/taxonomy.tsv', sep='\t', index_col=0)
epsps = pd.read_csv('epsps/epsps_classified.tsv', sep='\t', index_col=0)

def extract_genus_from_taxon(taxon_string):
    if pd.isna(taxon_string):
        return None
    parts = taxon_string.split(';')
    for part in parts:
        if part.startswith('g__'):
            genus = part.replace('g__', '').strip()
            return genus if genus else None
    return None

taxonomy['genus'] = taxonomy['Taxon'].apply(extract_genus_from_taxon)
epsps['genus'] = epsps.index.str.split('|').str[0]

def parse_metadata(sample_names):
    metadata_list = []
    for sample in sample_names:
        treatment = 'Control' if sample[0] == 'C' else 'Roundup'
        compartment = 'Gut' if '.In' in sample else 'Soil'
        metadata_list.append({'sample': sample, 'treatment': treatment, 'compartment': compartment})
    return pd.DataFrame(metadata_list).set_index('sample')

metadata = parse_metadata(feature_table.columns)

# Build genus->EPSPS class map
genus_to_class = {}
for genus in epsps['genus'].dropna().unique():
    genus_epsps = epsps[epsps['genus'] == genus].iloc[0]
    genus_to_class[genus] = genus_epsps['primary_class']

results = []

for sample in feature_table.columns:
    sample_counts = feature_table[sample]
    
    class_i_abundance = 0
    class_ii_abundance = 0
    total_classified = 0
    
    for asv_id in sample_counts.index:
        if asv_id not in taxonomy.index:
            continue
        
        genus = taxonomy.loc[asv_id, 'genus']
        if pd.isna(genus) or genus not in genus_to_class:
            continue
        
        primary_class = genus_to_class[genus]
        count = sample_counts[asv_id]
        
        if primary_class == 'I':
            class_i_abundance += count
        elif primary_class == 'II':
            class_ii_abundance += count
        
        total_classified += count
    
    if total_classified > 0:
        rel_class_i = class_i_abundance / total_classified
        rel_class_ii = class_ii_abundance / total_classified
    else:
        rel_class_i = 0
        rel_class_ii = 0
    
    treatment = metadata.loc[sample, 'treatment']
    compartment = metadata.loc[sample, 'compartment']
    
    results.append({
        'sample': sample,
        'treatment': treatment,
        'compartment': compartment,
        'class_I_sensitive': rel_class_i,
        'class_II_resistant': rel_class_ii
    })

results_df = pd.DataFrame(results)

print("=" * 70)
print("EPSPS CLASS RELATIVE ABUNDANCE BY TREATMENT AND COMPARTMENT")
print("=" * 70)
print()
print(results_df.groupby(['compartment', 'treatment'])[['class_I_sensitive', 'class_II_resistant']].agg(['mean', 'std', 'count']))
print()

print("=" * 70)
print("STATISTICAL TESTS (Mann-Whitney U)")
print("=" * 70)
print()

p_values = {}

for compartment in ['Gut', 'Soil']:
    print(f"\n{compartment.upper()} COMPARTMENT")
    print()
    
    comp_data = results_df[results_df['compartment'] == compartment]
    
    control_class_i = comp_data[comp_data['treatment'] == 'Control']['class_I_sensitive'].values
    roundup_class_i = comp_data[comp_data['treatment'] == 'Roundup']['class_I_sensitive'].values
    
    u_stat, p_value_i = mannwhitneyu(control_class_i, roundup_class_i)
    p_values[f'I_{compartment}'] = p_value_i
    print(f"Class I (Sensitive) by Treatment:")
    print(f"  Control (n={len(control_class_i)}): mean={control_class_i.mean():.4f}, sd={control_class_i.std():.4f}")
    print(f"  Roundup (n={len(roundup_class_i)}): mean={roundup_class_i.mean():.4f}, sd={roundup_class_i.std():.4f}")
    print(f"  Mann-Whitney U: U={u_stat:.1f}, p={p_value_i:.4f}")
    print()
    
    control_class_ii = comp_data[comp_data['treatment'] == 'Control']['class_II_resistant'].values
    roundup_class_ii = comp_data[comp_data['treatment'] == 'Roundup']['class_II_resistant'].values
    
    u_stat, p_value_ii = mannwhitneyu(control_class_ii, roundup_class_ii)
    p_values[f'II_{compartment}'] = p_value_ii
    print(f"Class II (Resistant) by Treatment:")
    print(f"  Control (n={len(control_class_ii)}): mean={control_class_ii.mean():.4f}, sd={control_class_ii.std():.4f}")
    print(f"  Roundup (n={len(roundup_class_ii)}): mean={roundup_class_ii.mean():.4f}, sd={roundup_class_ii.std():.4f}")
    print(f"  Mann-Whitney U: U={u_stat:.1f}, p={p_value_ii:.4f}")
    print()

# Plot 2x2: Class (I, II) x Compartment (Gut, Soil)
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

classes = [('class_I_sensitive', 'I'), ('class_II_resistant', 'II')]
class_titles = ['Class I (Sensitive)', 'Class II (Resistant)']
compartments = ['Gut', 'Soil']
color_schemes = {
    'Gut': ['#90EE90', '#0B3D0B'],
    'Soil': ['#F5DEB3', '#654321']
}

for row, ((class_col, class_key), class_title) in enumerate(zip(classes, class_titles)):
    for col, compartment in enumerate(compartments):
        ax = axes[row, col]
        
        comp_data = results_df[results_df['compartment'] == compartment]
        
        p_val = p_values[f'{class_key}_{compartment}']
        
        sns.boxplot(data=comp_data, x='treatment', y=class_col, ax=ax,
                   palette=color_schemes[compartment])
        sns.stripplot(data=comp_data, x='treatment', y=class_col, ax=ax,
                     color='black', alpha=0.4, size=6)
        
        ax.set_ylabel('Relative Abundance', fontsize=11)
        ax.set_xlabel('Treatment', fontsize=11)
        ax.set_title(f'{class_title} by Treatment\n{compartment} Compartment', 
                    fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        ax.text(0.5, 0.95, f'p = {p_val:.4f}', transform=ax.transAxes,
               ha='center', va='top', fontsize=11, fontweight='bold',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig('results/epsps_by_treatment/epsps_by_treatment.pdf', dpi=300, bbox_inches='tight')
print("Saved: results/epsps_by_treatment/epsps_by_treatment.pdf")

results_df.to_csv('results/epsps_by_treatment/epsps_by_treatment.csv', index=False)
print("Saved: results/epsps_by_treatment/epsps_by_treatment.csv")

print("=" * 70)
