#!/usr/bin/env python3
"""
Test script to demonstrate the new axial ligand coloring functionality 
in the comprehensive iron-plane analysis function.
"""

import pandas as pd
import numpy as np
from calculate_rmsd import RMSDAnalyzer

def test_comprehensive_function():
    """Test the enhanced comprehensive iron-plane analysis function."""
    
    # Create an analyzer instance
    analyzer = RMSDAnalyzer(verbose=True)
    
    # Create mock iron-plane data
    pdb_ids = ['1abc', '2def', '3ghi', '4jkl', '5mno']
    iron_plane_data = pd.DataFrame({
        'PDB_ID': pdb_ids,
        'PDB_Iron_Plane_Distance': np.random.uniform(0.1, 0.8, 5),
        'XYZ_Iron_Plane_Distance': np.random.uniform(0.2, 0.9, 5),
        'Plane_Distance_Difference': np.random.uniform(-0.3, 0.3, 5)
    })
    
    # Create mock axial ligand data
    axial_data = pd.DataFrame({
        'PDB_ID': ['1abc', '1abc', '2def', '2def', '3ghi', '3ghi', '4jkl', '4jkl', '5mno', '5mno'],
        'Axial_Ligand': [1, 2, 1, 2, 1, 2, 1, 2, 1, 2],
        'Axial_Resname': ['CYS', 'HIS', 'HIS', 'MET', 'CYS', 'HOH', 'HIS', 'OXY', 'CYS', 'HIS']
    })
    
    print("Testing comprehensive iron-plane analysis function...")
    print("=" * 60)
    
    # Test without axial ligand coloring (default behavior)
    print("\n1. Testing without axial ligand coloring:")
    analyzer.create_comprehensive_iron_plane_analysis(
        iron_plane_data, 
        output_file="plots/test_comprehensive_default.png",
        layout='horizontal'
    )
    
    # Test with axial ligand coloring
    print("\n2. Testing with axial ligand coloring:")
    analyzer.create_comprehensive_iron_plane_analysis(
        iron_plane_data, 
        output_file="plots/test_comprehensive_axial_colors.png",
        layout='horizontal',
        color_by_axial_ligands=True,
        axial_data=axial_data
    )
    
    # Test vertical layout with axial coloring
    print("\n3. Testing vertical layout with axial ligand coloring:")
    analyzer.create_comprehensive_iron_plane_analysis(
        iron_plane_data, 
        output_file="plots/test_comprehensive_vertical_axial.png",
        layout='vertical',
        color_by_axial_ligands=True,
        axial_data=axial_data
    )
    
    print("\n" + "=" * 60)
    print("Test completed! Check the generated plots in the plots/ directory.")
    print("The function now supports:")
    print("- color_by_axial_ligands=True/False flag")
    print("- Automatic axial ligand combination detection")
    print("- Color coding by axial ligand pairs (e.g., CYS-HIS, HIS-MET)")
    print("- Compatible with all layout modes (horizontal, vertical, single_overlay)")

if __name__ == "__main__":
    test_comprehensive_function()