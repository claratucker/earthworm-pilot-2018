#!/usr/bin/env python3
"""Generate EPSPS summary table for publication."""

import pandas as pd
import numpy as np
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
        'Compartment': compartment,
        'Treatment': treatment,
        'Class I (Sensitive)': rel_class_i,
        'Class II (Resistant)': rel_class_ii
    })

results_df = pd.DataFrame(results)

# Create summary table
summary_rows = []

for compartment in ['Gut', 'Soil']:
    for epsps_class in ['Class I (Sensitive)', 'Class II (Resistant)']:
        control_vals = results_df[(results_df['Compartment'] == compartment) & 
                                 (results_df['Treatment'] == 'Control')][epsps_class].values
        roundup_vals = results_df[(results_df['Compartment'] == compartment) & 
                                 (results_df['Treatment'] == 'Roundup')][epsps_class].values
        
        u_stat, p_val = mannwhitneyu(control_vals, roundup_vals)
        
        summary_rows.append({
            'Compartment': compartment,
            'EPSPS Class': epsps_class,
            'Control Mean': f"{control_vals.mean():.4f}",
            'Control SD': f"{control_vals.std():.4f}",
            'Control n': len(control_vals),
            'Roundup Mean': f"{roundup_vals.mean():.4f}",
            'Roundup SD': f"{roundup_vals.std():.4f}",
            'Roundup n': len(roundup_vals),
            'U Statistic': f"{u_stat:.1f}",
            'p-value': f"{p_val:.4f}",
            'Significance': 'NS' if p_val >= 0.05 else f'*'
        })

summary_table = pd.DataFrame(summary_rows)

print("\n" + "=" * 120)
print("TABLE: EPSPS Class Distribution by Treatment and Compartment")
print("=" * 120)
print()
print(summary_table.to_string(index=False))
print()
print("NS = not significant (p >= 0.05)")
print("=" * 120)

# Save as both CSV and formatted text
summary_table.to_csv('results/epsps_by_treatment/epsps_summary_table.csv', index=False)
print("\nSaved: results/epsps_by_treatment/epsps_summary_table.csv")

with open('results/epsps_by_treatment/epsps_summary_table.txt', 'w') as f:
    f.write("TABLE: EPSPS Class Distribution by Treatment and Compartment\n")
    f.write("=" * 120 + "\n\n")
    f.write(summary_table.to_string(index=False))
    f.write("\n\nNS = not significant (p >= 0.05)\n")

print("Saved: results/epsps_by_treatment/epsps_summary_table.txt")
