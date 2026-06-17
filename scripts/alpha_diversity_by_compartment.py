#!/usr/bin/env python3
"""Alpha diversity (Shannon/Richness) stratified by compartment."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

print("Loading alpha diversity data...")
alpha = pd.read_csv('r/alpha_diversity.csv', index_col=0)

alpha['compartment'] = alpha['compartment'].str.capitalize()
alpha['treatment'] = alpha['treatment'].str.capitalize()

print("=" * 70)
print("ALPHA DIVERSITY BY COMPARTMENT")
print("=" * 70)
print()

p_values = {}

for metric in ['Shannon', 'Observed']:
    print(f"{metric.upper()}")
    print()
    
    gut = alpha[alpha['compartment'] == 'Gut'][metric].values
    soil = alpha[alpha['compartment'] == 'Soil'][metric].values
    
    u_stat, p_value_comp = mannwhitneyu(gut, soil)
    p_values[f'{metric}_compartment'] = p_value_comp
    print(f"  Gut vs Soil: U={u_stat:.1f}, p={p_value_comp:.4f}")
    print(f"    Gut (n={len(gut)}): {gut.mean():.4f} ± {gut.std():.4f}")
    print(f"    Soil (n={len(soil)}): {soil.mean():.4f} ± {soil.std():.4f}")
    print()
    
    for compartment in ['Gut', 'Soil']:
        control = alpha[(alpha['compartment'] == compartment) & (alpha['treatment'] == 'Control')][metric].values
        roundup = alpha[(alpha['compartment'] == compartment) & (alpha['treatment'] == 'Roundup')][metric].values
        
        u_stat, p_value = mannwhitneyu(control, roundup)
        p_values[f'{metric}_{compartment}'] = p_value
        print(f"  {compartment}: Control vs Roundup: U={u_stat:.1f}, p={p_value:.4f}")
        print(f"    Control (n={len(control)}): {control.mean():.4f} ± {control.std():.4f}")
        print(f"    Roundup (n={len(roundup)}): {roundup.mean():.4f} ± {roundup.std():.4f}")
        print()

# Plot
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

metrics = [('Shannon', 'Shannon Diversity'), ('Observed', 'Observed Richness')]

for row, (metric, title) in enumerate(metrics):
    # By compartment
    ax = axes[row, 0]
    sns.boxplot(data=alpha, x='compartment', y=metric, ax=ax, palette=['#C9A876', '#8B5A3C'])
    sns.stripplot(data=alpha, x='compartment', y=metric, ax=ax, color='black', alpha=0.4, size=6)
    ax.set_xlabel('Compartment', fontsize=11)
    ax.set_ylabel(title, fontsize=11)
    ax.set_title(f'{title} by Compartment', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    p_val = p_values[f'{metric}_compartment']
    ax.text(0.5, 0.95, f'p = {p_val:.4f}', transform=ax.transAxes,
           ha='center', va='top', fontsize=10, fontweight='bold',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # By treatment within compartment
    ax = axes[row, 1]
    
    # Create custom palette for each compartment with light/dark colors
    gut_data = alpha[alpha['compartment'] == 'Gut']
    soil_data = alpha[alpha['compartment'] == 'Soil']
    
    # Plot separately to use correct colors
    for compartment, color_palette in [('Gut', ['#90EE90', '#0B3D0B']), ('Soil', ['#F5DEB3', '#654321'])]:
        comp_data = alpha[alpha['compartment'] == compartment]
        for treat_idx, treatment in enumerate(['Control', 'Roundup']):
            subset = comp_data[comp_data['treatment'] == treatment]
            ax.scatter([compartment]*len(subset), subset[metric],
                      color=color_palette[treat_idx], s=100, alpha=0.6,
                      edgecolors='black', linewidth=0.8, label=f'{compartment} {treatment}')
            
            # Add box for each group
            pos = ['Gut', 'Soil'].index(compartment)
            bp = ax.boxplot([subset[metric]], positions=[pos + (treat_idx-0.2)], widths=0.35,
                           patch_artist=True, showfliers=False)
            for patch in bp['boxes']:
                patch.set_facecolor(color_palette[treat_idx])
                patch.set_alpha(0.3)
    
    ax.set_xlabel('Compartment', fontsize=11)
    ax.set_ylabel(title, fontsize=11)
    ax.set_title(f'{title} by Compartment and Treatment', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Gut', 'Soil'])
    
    # Add p-values for treatment comparisons
    p_gut = p_values[f'{metric}_Gut']
    p_soil = p_values[f'{metric}_Soil']
    
    stats_text = f"Gut: p={p_gut:.4f}\nSoil: p={p_soil:.4f}"
    ax.text(0.98, 0.05, stats_text, transform=ax.transAxes,
           ha='right', va='bottom', fontsize=9,
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig('results/alpha_diversity_by_compartment.pdf', dpi=300, bbox_inches='tight')
print("\nSaved: results/alpha_diversity_by_compartment.pdf")

alpha.to_csv('results/alpha_diversity_by_compartment.csv')
print("Saved: results/alpha_diversity_by_compartment.csv")

print("=" * 70)
