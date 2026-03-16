#!/usr/bin/env python3
"""
Script to create overlay plot of RMSD data from two CSV files:
- tables/rmsd_all.csv (RMSD_all column)
- tables/rmsd_results.csv (RMSD column)
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def create_rmsd_overlay_plot():
    # Read the CSV files
    rmsd_all_df = pd.read_csv('tables/rmsd_all.csv')
    rmsd_results_df = pd.read_csv('tables/rmsd_results.csv')
    
    # Extract RMSD values
    rmsd_all_values = rmsd_all_df['RMSD_all']
    rmsd_results_values = rmsd_results_df['RMSD']
    
    # Create the overlay plot
    plt.figure(figsize=(12, 8))
    
    # Plot histograms with transparency
    plt.hist(rmsd_all_values, bins=30, alpha=0.6, label='RMSD_all', color='blue', density=True)
    plt.hist(rmsd_results_values, bins=30, alpha=0.6, label='RMSD_results', color='red', density=True)
    
    # Add statistical information
    plt.axvline(rmsd_all_values.mean(), color='blue', linestyle='--', alpha=0.8, 
                label=f'RMSD_all mean: {rmsd_all_values.mean():.3f}')
    plt.axvline(rmsd_results_values.mean(), color='red', linestyle='--', alpha=0.8,
                label=f'RMSD_results mean: {rmsd_results_values.mean():.3f}')
    
    # Formatting
    plt.xlabel('RMSD Value (Å)')
    plt.ylabel('Density')
    plt.title('Overlay of RMSD Distributions')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Save and show
    plt.tight_layout()
    plt.savefig('rmsd_overlay_plot.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print summary statistics
    print(f"RMSD_all dataset:")
    print(f"  Count: {len(rmsd_all_values)}")
    print(f"  Mean: {rmsd_all_values.mean():.3f}")
    print(f"  Std: {rmsd_all_values.std():.3f}")
    print(f"  Min: {rmsd_all_values.min():.3f}")
    print(f"  Max: {rmsd_all_values.max():.3f}")
    
    print(f"\nRMSD_results dataset:")
    print(f"  Count: {len(rmsd_results_values)}")
    print(f"  Mean: {rmsd_results_values.mean():.3f}")
    print(f"  Std: {rmsd_results_values.std():.3f}")
    print(f"  Min: {rmsd_results_values.min():.3f}")
    print(f"  Max: {rmsd_results_values.max():.3f}")

if __name__ == "__main__":
    create_rmsd_overlay_plot()