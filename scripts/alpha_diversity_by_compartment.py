#!/usr/bin/env python3
"""Alpha diversity (Shannon/Simpson/Chao1) stratified by compartment."""

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

print("SHANNON DIVERSITY")
print()
print(alpha.groupby(['compartment', 'treatment'])['Shannon'].agg(['mean', 'std', 'count']))
print()

print("OBSERVED RICHNESS")
print()
print(alpha.groupby(['compartment', 'treatment'])['Observed'].agg(['mean', 'std', 'count']))
print()

print("=" * 70)
print("STATISTICAL TESTS")
print("=" * 70)
print()

for metric in ['Shannon', 'Observed']:
    print(f"{metric.upper()}")
    print()
    
    # By compartment
    gut = alpha[alpha['compartment'] == 'Gut'][metric].values
    soil = alpha[alpha['compartment'] == 'Soil'][metric].values
    
    u_stat, p_value = mannwhitneyu(gut, soil)
    print(f"  Gut vs Soil: U={u_stat:.1f}, p={p_value:.4f}")
    print(f"    Gut (n={len(gut)}): {gut.mean():.4f} ± {gut.std():.4f}")
    print(f"    Soil (n={len(soil)}): {soil.mean():.4f} ± {soil.std():.4f}")
    print()
    
    # By treatment within each compartment
    for compartment in ['Gut', 'Soil']:
        control = alpha[(alpha['compartment'] == compartment) & (alpha['treatment'] == 'Control')][metric].values
        roundup = alpha[(alpha['compartment'] == compartment) & (alpha['treatment'] == 'Roundup')][metric].values
        
        u_stat, p_value = mannwhitneyu(control, roundup)
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
    sns.boxplot(data=alpha, x='compartment', y=metric, ax=ax, palette=['#9ECAE1', '#08519C'])
    sns.stripplot(data=alpha, x='compartment', y=metric, ax=ax, color='black', alpha=0.4, size=6)
    ax.set_xlabel('Compartment', fontsize=11)
    ax.set_ylabel(title, fontsize=11)
    ax.set_title(f'{title} by Compartment', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    # By treatment within compartment
    ax = axes[row, 1]
    sns.boxplot(data=alpha, x='compartment', y=metric, hue='treatment', ax=ax, palette=['#DEEBF7', '#08519C'])
    ax.set_xlabel('Compartment', fontsize=11)
    ax.set_ylabel(title, fontsize=11)
    ax.set_title(f'{title} by Compartment and Treatment', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    ax.legend(title='Treatment', loc='best')

plt.tight_layout()
plt.savefig('results/alpha_diversity_by_compartment.pdf', dpi=300, bbox_inches='tight')
print("Saved: results/alpha_diversity_by_compartment.pdf")

alpha.to_csv('results/alpha_diversity_by_compartment.csv')
print("Saved: results/alpha_diversity_by_compartment.csv")

print("=" * 70)
