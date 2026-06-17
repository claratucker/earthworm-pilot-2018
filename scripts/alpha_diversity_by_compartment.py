#!/usr/bin/env python3
"""Alpha diversity (Shannon/Richness) stratified by compartment."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

print("Loading alpha diversity data...")
alpha = pd.read_csv('r/alpha_diversity.csv', index_col=0)

# Capitalize for consistency
alpha['compartment'] = alpha['compartment'].str.capitalize()
alpha['treatment'] = alpha['treatment'].str.capitalize()

print("=" * 70)
print("ALPHA DIVERSITY BY COMPARTMENT")
print("=" * 70)
print()

# Statistical tests
print("STATISTICAL TESTS")
print()

p_values = {}

for metric in ['Shannon', 'Observed']:
    print(f"{metric.upper()}")
    print()
    
    # By compartment
    gut = alpha[alpha['compartment'] == 'Gut'][metric].values
    soil = alpha[alpha['compartment'] == 'Soil'][metric].values
    
    u_stat, p_value_comp = mannwhitneyu(gut, soil)
    p_values[f'{metric}_compartment'] = p_value_comp
    print(f"  Gut vs Soil: U={u_stat:.1f}, p={p_value_comp:.4f}")
    print(f"    Gut (n={len(gut)}): {gut.mean():.4f} ± {gut.std():.4f}")
    print(f"    Soil (n={len(soil)}): {soil.mean():.4f} ± {soil.std():.4f}")
    print()
    
    # By treatment within each compartment
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
    
    # Add p-value for compartment comparison
    p_val = p_values[f'{metric}_compartment']
    sig_text = f"p={p_val:.3f}*" if p_val < 0.05 else f"p={p_val:.3f}"
    ax.text(0.5, 0.95, sig_text, transform=ax.transAxes,
           ha='center', va='top', fontsize=10, fontweight='bold',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # By treatment within compartment
    ax = axes[row, 1]
    sns.boxplot(data=alpha, x='compartment', y=metric, hue='treatment', ax=ax, palette=['#D2B48C', '#654321'])
    ax.set_xlabel('Compartment', fontsize=11)
    ax.set_ylabel(title, fontsize=11)
    ax.set_title(f'{title} by Compartment and Treatment', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    ax.legend(title='Treatment', loc='best')
    
    # Add p-values for treatment comparisons
    p_gut = p_values[f'{metric}_Gut']
    p_soil = p_values[f'{metric}_Soil']
    
    stats_text = f"Gut: p={p_gut:.3f}\nSoil: p={p_soil:.3f}"
    ax.text(0.98, 0.05, stats_text, transform=ax.transAxes,
           ha='right', va='bottom', fontsize=9,
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig('results/alpha_diversity_by_compartment.pdf', dpi=300, bbox_inches='tight')
print("\nSaved: results/alpha_diversity_by_compartment.pdf")

alpha.to_csv('results/alpha_diversity_by_compartment.csv')
print("Saved: results/alpha_diversity_by_compartment.csv")

print("=" * 70)
