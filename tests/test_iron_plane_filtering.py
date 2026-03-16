#!/usr/bin/env python3
"""
Test script to verify that iron plane distance analysis properly filters
axial ligand combinations and that the plots only show allowed combinations.
"""

import sys
import os
import pandas as pd

# Add current directory to path to import calculate_rmsd
sys.path.insert(0, os.getcwd())

from calculate_rmsd import RMSDAnalyzer

def test_iron_plane_filtering():
    """Test the iron plane distance analysis filtering."""
    
    print("Testing iron plane distance analysis with axial ligand filtering...")
    print("=" * 70)
    
    # Create analyzer with filtering enabled
    analyzer = RMSDAnalyzer(verbose=True, filter_axial_ligands=True)
    
    print(f"\nStep 1: Processing iron plane distances with filtering enabled...")
    df_plane = analyzer.process_iron_plane_distances()
    
    if not df_plane.empty:
        print(f"\nResults after filtering:")
        print(f"  Total structures analyzed: {len(df_plane)}")
        
        if 'Axial_Combination' in df_plane.columns:
            print(f"  Unique axial combinations found:")
            axial_counts = df_plane['Axial_Combination'].value_counts()
            for combo, count in axial_counts.items():
                print(f"    {combo}: {count} structures")
            
            # Verify that all combinations are allowed
            print(f"\nStep 2: Verifying all combinations are allowed...")
            all_allowed = True
            for combo in axial_counts.index:
                if combo != "Unknown":
                    # Split the combination and check if it's allowed
                    parts = combo.split('-')
                    if len(parts) == 2:
                        is_allowed = analyzer.is_allowed_axial_combination(parts[0], parts[1])
                        status = "✓" if is_allowed else "✗"
                        print(f"    {status} {combo}: {'Allowed' if is_allowed else 'NOT ALLOWED'}")
                        if not is_allowed:
                            all_allowed = False
            
            if all_allowed:
                print(f"\n✓ SUCCESS: All axial combinations in iron plane analysis are from the allowed list!")
            else:
                print(f"\n✗ ERROR: Some non-allowed axial combinations found in iron plane analysis!")
        
        # Save the results for manual inspection
        output_file = "tables/test_iron_plane_distances.csv"
        df_plane.to_csv(output_file, index=False)
        print(f"\nIron plane distance results saved to: {output_file}")
        
        # Test the comprehensive analysis method as well
        print(f"\nStep 3: Testing comprehensive iron plane analysis method...")
        
        # Get some axial data to test coloring
        axial_df = analyzer.process_iron_axial_distances(pdb_ids_filter=df_plane['PDB_ID'].unique()[:10])  # Test with first 10
        
        if not axial_df.empty:
            print(f"  Axial data for testing: {len(axial_df)} measurements")
            print(f"  Testing comprehensive analysis with axial coloring...")
            
            # Create a small test plot
            test_output = "plots/test_comprehensive_iron_plane_analysis.png"
            analyzer.create_comprehensive_iron_plane_analysis(
                df_plane.head(20),  # Use only first 20 for testing
                output_file=test_output,
                color_by_axial_ligands=True,
                axial_data=axial_df
            )
            print(f"  Test plot saved to: {test_output}")
        
    else:
        print(f"\n⚠ WARNING: No iron plane distance data found after filtering!")
        print(f"This could mean:")
        print(f"  - No structures have the required files")
        print(f"  - All structures have non-allowed axial combinations")
        print(f"  - There's an issue with the data processing")
    
    print(f"\nTest completed!")

if __name__ == "__main__":
    test_iron_plane_filtering()