# Iron Axial Distance Analysis Code Snippet
# From calculate_rmsd.py, function: create_distance_difference_histogram()

"""
ANALYSIS: This code is measuring DISTANCE DIFFERENCES, not absolute distances.

Key evidence from the code:
1. Function name: create_distance_difference_histogram()
2. Column used: 'Distance_Difference' 
3. Comment: "Create separate histograms of Fe-axial ligand distance differences"
4. X-axis label: 'Fe-X Distance Difference (Å)'

This means the plot shows the difference between the Fe-axial ligand distance
in the optimized structure compared to some reference (likely the crystal structure).

The distance difference is calculated as:
Distance_Difference = Distance_Optimized - Distance_Reference

Positive values = axial ligand moved away from Fe during optimization
Negative values = axial ligand moved closer to Fe during optimization
Zero = no change in distance
"""

def create_distance_difference_histogram(self, df, output_file="plots/iron_axial_distance_differences_stacked.png", exclude_suspicious=True):
    """
    Create separate histograms of Fe-axial ligand distance differences for each axial ligand, 
    where each bar is colored by axial ligand combinations for the most common 5 combinations.
    Optionally excludes suspicious PDB IDs from plotting.
    """
    if df.empty or 'Distance_Difference' not in df.columns:
        if self.verbose:
            print("No distance difference data to plot!")
        return

    # Filter out suspicious structures
    if exclude_suspicious:
        suspicious_ids = self.get_suspicious_pdb_ids()
        if suspicious_ids:
            df = df[~df['PDB_ID'].isin(suspicious_ids)]
            if self.verbose:
                print(f"Excluded {len(suspicious_ids)} suspicious structures")

    if df.empty:
        print("No data remaining after filtering suspicious structures!")
        return

    # Get axial ligand data for coloring
    plot_df = df.copy()
    use_axial_colors = False
    axial_color_map = None
    axial_combinations = None
    
    if 'PDB_ID' in plot_df.columns:
        # Get axial ligand information
        valid_pdb_ids = plot_df['PDB_ID'].unique()
        axial_df = self.process_iron_axial_distances(pdb_ids_filter=valid_pdb_ids)
        
        if not axial_df.empty:
            # Create axial combinations by pivoting data
            axial_pivot = axial_df.pivot_table(
                index='PDB_ID', 
                columns='Axial_Ligand', 
                values='Axial_Resname', 
                aggfunc='first'
            ).reset_index()
            
            # Create proper axial combinations (axial1-axial2) with HSD normalization
            if 1 in axial_pivot.columns and 2 in axial_pivot.columns:
                # Normalize ligand names (HSD -> HIS) and create sorted combination
                def create_normalized_combo(row):
                    if pd.notna(row[1]) and pd.notna(row[2]):
                        norm_axial1 = self.normalize_axial_ligand_name(row[1])
                        norm_axial2 = self.normalize_axial_ligand_name(row[2])
                        return '-'.join(sorted([norm_axial1, norm_axial2]))
                    else:
                        return 'Unknown-Unknown'
                axial_pivot['Axial_Combo'] = axial_pivot.apply(create_normalized_combo, axis=1)
                
    # Create figure with subplots for each axial ligand
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Function to create histogram for given axial ligand
    def create_single_histogram(ax, axial_ligand_num, title):
        axial_data = plot_df[plot_df['Axial_Ligand'] == axial_ligand_num]
        
        if axial_data.empty:
            ax.set_xlabel('Fe-X Distance Difference (Å)', fontsize=14)  # NOTE: "Difference" in label
            ax.set_ylabel('Frequency', fontsize=14)
            return
        
        distances = axial_data['Distance_Difference'].dropna()  # NOTE: Using 'Distance_Difference' column
        
        # Rest of plotting code...
        
# CONCLUSION: This plot shows Fe-axial ligand DISTANCE DIFFERENCES (changes upon optimization),
# NOT absolute distances. The filename should reflect this - it's correctly named as 
# "iron_axial_distance_differences_stacked.png"