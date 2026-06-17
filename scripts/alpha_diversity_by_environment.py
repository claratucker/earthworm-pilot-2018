#!/usr/bin/env python3
"""Alpha diversity (Shannon/Richness) stratified by environment."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

print("Loading alpha diversity data...")
alpha = pd.read_csv('r/alpha_diversity.csv', index_col=0)
alpha = alpha.rename(columns={'compartment': 'environment'})

alpha['environment'] = alpha['environment'].str.capitalize()
alpha['treatment'] = alpha['treatment'].str.capitalize()

print("=" * 70)
print("ALPHA DIVERSITY BY ENVIRONMENT")
print("=" * 70)
print()

p_values = {}

for metric in ['Shannon', 'Observed']:
    print(f"{metric.upper()}")
    print()
    
    gut = alpha[alpha['environment'] == 'Gut'][metric].values
    soil = alpha[alpha['environment'] == 'Soil'][metric].values
    
    u_stat, p_value_comp = mannwhitneyu(gut, soil)
    p_values[f'{metric}_environment'] = p_value_comp
    print(f"  Gut vs Soil: U={u_stat:.1f}, p={p_value_comp:.4f}")
    print(f"    Gut (n={len(gut)}): {gut.mean():.4f} ± {gut.std():.4f}")
    print(f"    Soil (n={len(soil)}): {soil.mean():.4f} ± {soil.std():.4f}")
    print()
    
    for environment in ['Gut', 'Soil']:
        control = alpha[(alpha['environment'] == environment) & (alpha['treatment'] == 'Control')][metric].values
        roundup = alpha[(alpha['environment'] == environment) & (alpha['treatment'] == 'Roundup')][metric].values
        
        u_stat, p_value = mannwhitneyu(control, roundup)
        p_values[f'{metric}_{environment}'] = p_value
        print(f"  {environment}: Control vs Roundup: U={u_stat:.1f}, p={p_value:.4f}")
        print(f"    Control (n={len(control)}): {control.mean():.4f} ± {control.std():.4f}")
        print(f"    Roundup (n={len(roundup)}): {roundup.mean():.4f} ± {roundup.std():.4f}")
        print()

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

metrics = [('Shannon', 'Shannon Diversity'), ('Observed', 'Observed Richness')]

for row, (metric, metric_title) in enumerate(metrics):
    
    # Left: By environment
    ax = axes[row, 0]
    sns.boxplot(data=alpha, x='environment', y=metric, ax=ax, palette=['#C9A876', '#8B5A3C'])
    sns.stripplot(data=alpha, x='environment', y=metric, ax=ax, color='black', alpha=0.4, size=6)
    ax.set_xlabel('Environment', fontsize=11)
    ax.set_ylabel(metric_title, fontsize=11)
    ax.set_title(f'{metric_title} by Environment', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    p_val = p_values[f'{metric}_environment']
    ax.text(0.5, 0.95, f'p = {p_val:.4f}', transform=ax.transAxes,
           ha='center', va='top', fontsize=10, fontweight='bold',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Right: By treatment within each environment - separate panels
    ax = axes[row, 1]
    
    # Create palette mapping for all combinations
    palette = {
        ('Gut', 'Control'): '#90EE90',
        ('Gut', 'Roundup'): '#0B3D0B',
        ('Soil', 'Control'): '#F5DEB3',
        ('Soil', 'Roundup'): '#654321'
    }
    
    # Create a hue order for the plot
    alpha['environment_treatment'] = alpha['environment'] + '\n' + alpha['treatment']
    plot_order = ['Gut\nControl', 'Gut\nRoundup', 'Soil\nControl', 'Soil\nRoundup']
    
    # Map colors
    colors = [palette[(row.split('\n')[0], row.split('\n')[1])] for row in plot_order]
    
    sns.boxplot(data=alpha, x='environment_treatment', y=metric, ax=ax,
               palette=colors, order=plot_order)
    sns.stripplot(data=alpha, x='environment_treatment', y=metric, ax=ax,
                 color='black', alpha=0.4, size=6, order=plot_order)
    
    ax.set_xlabel('Environment and Treatment', fontsize=11)
    ax.set_ylabel(metric_title, fontsize=11)
    ax.set_title(f'{metric_title} by Treatment\nwithin Environment', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add p-values for treatment comparisons
    p_gut = p_values[f'{metric}_Gut']
    p_soil = p_values[f'{metric}_Soil']
    
    stats_text = f"Gut: p={p_gut:.4f}\nSoil: p={p_soil:.4f}"
    ax.text(0.98, 0.05, stats_text, transform=ax.transAxes,
           ha='right', va='bottom', fontsize=9,
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig('results/alpha_diversity_by_environment/alpha_diversity_by_environment.pdf', dpi=300, bbox_inches='tight')
print("\nSaved: results/alpha_diversity_by_environment/alpha_diversity_by_environment.pdf")

alpha.to_csv('results/alpha_diversity_by_environment/alpha_diversity_by_environment.csv')
print("Saved: results/alpha_diversity_by_environment/alpha_diversity_by_environment.csv")

print("=" * 70)
