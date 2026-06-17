#!/usr/bin/env python3
"""Alpha diversity (Shannon/Richness) stratified by compartment."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

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

# Plot 2x2: Metric (Shannon, Observed) x Comparison (by compartment, by treatment)
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

metrics = [('Shannon', 'Shannon Diversity'), ('Observed', 'Observed Richness')]

for row, (metric, metric_title) in enumerate(metrics):
    
    # Left: By compartment (both treatments combined)
    ax = axes[row, 0]
    sns.boxplot(data=alpha, x='compartment', y=metric, ax=ax, palette=['#C9A876', '#8B5A3C'])
    sns.stripplot(data=alpha, x='compartment', y=metric, ax=ax, color='black', alpha=0.4, size=6)
    ax.set_xlabel('Compartment', fontsize=11)
    ax.set_ylabel(metric_title, fontsize=11)
    ax.set_title(f'{metric_title} by Compartment', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    p_val = p_values[f'{metric}_compartment']
    ax.text(0.5, 0.95, f'p = {p_val:.4f}', transform=ax.transAxes,
           ha='center', va='top', fontsize=10, fontweight='bold',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Right: By treatment within compartment
    ax = axes[row, 1]
    
    color_schemes = {
        'Gut': ['#90EE90', '#0B3D0B'],
        'Soil': ['#F5DEB3', '#654321']
    }
    
    x_pos = 0
    positions = []
    labels = []
    tick_positions = []
    tick_labels = []
    
    for comp_idx, compartment in enumerate(['Gut', 'Soil']):
        comp_data = alpha[alpha['compartment'] == compartment]
        
        for treat_idx, treatment in enumerate(['Control', 'Roundup']):
            subset = comp_data[comp_data['treatment'] == treatment]
            ax.scatter([x_pos]*len(subset), subset[metric],
                      color=color_schemes[compartment][treat_idx], s=100, alpha=0.6,
                      edgecolors='black', linewidth=0.8)
            
            bp = ax.boxplot([subset[metric]], positions=[x_pos], widths=0.4,
                           patch_artist=True, showfliers=False)
            for patch in bp['boxes']:
                patch.set_facecolor(color_schemes[compartment][treat_idx])
                patch.set_alpha(0.4)
                patch.set_linewidth(1.5)
            
            x_pos += 1
        
        # Add separator and label
        if comp_idx == 0:
            tick_positions.append(0.5)
            tick_labels.append('Gut')
        else:
            tick_positions.append(2.5)
            tick_labels.append('Soil')
        
        x_pos += 0.5
    
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=11)
    ax.set_ylabel(metric_title, fontsize=11)
    ax.set_title(f'{metric_title} by Treatment\nwithin Compartment', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_xlim(-0.5, x_pos - 0.5)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#90EE90', edgecolor='black', label='Gut Control'),
        Patch(facecolor='#0B3D0B', edgecolor='black', label='Gut Roundup'),
        Patch(facecolor='#F5DEB3', edgecolor='black', label='Soil Control'),
        Patch(facecolor='#654321', edgecolor='black', label='Soil Roundup')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
    
    # Add p-values for treatment comparisons
    p_gut = p_values[f'{metric}_Gut']
    p_soil = p_values[f'{metric}_Soil']
    
    stats_text = f"Gut: p={p_gut:.4f}\nSoil: p={p_soil:.4f}"
    ax.text(0.98, 0.05, stats_text, transform=ax.transAxes,
           ha='right', va='bottom', fontsize=9,
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig('results/alpha_diversity_by_compartment/alpha_diversity_by_compartment.pdf', dpi=300, bbox_inches='tight')
print("\nSaved: results/alpha_diversity_by_compartment/alpha_diversity_by_compartment.pdf")

alpha.to_csv('results/alpha_diversity_by_compartment/alpha_diversity_by_compartment.csv')
print("Saved: results/alpha_diversity_by_compartment/alpha_diversity_by_compartment.csv")

print("=" * 70)
