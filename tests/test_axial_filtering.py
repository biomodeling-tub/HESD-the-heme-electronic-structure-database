#!/usr/bin/env python3
"""
Test script to verify that axial ligand filtering is working correctly
and that HSD is properly normalized to HIS.
"""

import sys
import os

# Add current directory to path to import calculate_rmsd
sys.path.insert(0, os.getcwd())

from calculate_rmsd import RMSDAnalyzer

def test_axial_filtering():
    """Test the axial ligand filtering functionality."""
    
    print("Testing axial ligand filtering and HSD normalization...")
    print("=" * 60)
    
    # Create analyzer with filtering enabled
    analyzer = RMSDAnalyzer(verbose=True, filter_axial_ligands=True)
    
    print(f"\nAllowed axial combinations:")
    for combo in analyzer.allowed_axial_combinations:
        print(f"  {sorted(combo)}")
    
    # Test HSD normalization
    print(f"\nTesting HSD normalization:")
    test_cases = [
        ("HSD", "HIS"),
        ("HIS", "HIS"), 
        ("CYS", "CYS"),
        ("HOH", "HOH"),
        ("MET", "MET"),
        ("OXY", "OXY"),
        (None, None)
    ]
    
    for input_val, expected in test_cases:
        result = analyzer.normalize_axial_ligand_name(input_val)
        status = "✓" if result == expected else "✗"
        print(f"  {status} {input_val} -> {result} (expected: {expected})")
    
    # Test axial combination checking
    print(f"\nTesting axial combination validation:")
    test_combinations = [
        ("HIS", "HIS", True),     # Should be allowed
        ("HSD", "HIS", True),     # Should be allowed (HSD -> HIS)
        ("HSD", "HSD", True),     # Should be allowed (both HSD -> HIS)
        ("CYS", "HOH", True),     # Should be allowed
        ("HIS", "MET", True),     # Should be allowed
        ("HIS", "OXY", True),     # Should be allowed
        ("HIS", "HOH", True),     # Should be allowed
        ("CYS", "CYS", False),    # Should NOT be allowed
        ("MET", "MET", False),    # Should NOT be allowed
        ("HIS", "LYS", False),    # Should NOT be allowed
        ("ARG", "ASP", False),    # Should NOT be allowed
    ]
    
    for axial1, axial2, expected in test_combinations:
        result = analyzer.is_allowed_axial_combination(axial1, axial2)
        status = "✓" if result == expected else "✗"
        # Show normalized combination
        norm1 = analyzer.normalize_axial_ligand_name(axial1)
        norm2 = analyzer.normalize_axial_ligand_name(axial2)
        norm_combo = '-'.join(sorted([norm1, norm2]))
        print(f"  {status} {axial1}-{axial2} -> {norm_combo} (allowed: {result}, expected: {expected})")
    
    # Generate and save filtered combinations report
    print(f"\nGenerating filtered combinations report...")
    try:
        filtered_df = analyzer.save_filtered_axial_combinations_report()
        if not filtered_df.empty:
            print(f"Found {len(filtered_df)} structures with non-allowed axial combinations")
            print(f"Report saved to: tables/filtered_axial_combinations_report.csv")
            
            # Show some examples
            print(f"\nTop 5 most common filtered combinations:")
            if 'Original_Combination' in filtered_df.columns:
                top_combos = filtered_df['Original_Combination'].value_counts().head(5)
                for combo, count in top_combos.items():
                    print(f"  {combo}: {count} structures")
        else:
            print("No filtered combinations found (all structures have allowed axial ligands)")
    except Exception as e:
        print(f"Error generating report: {e}")
    
    print(f"\nTest completed!")

if __name__ == "__main__":
    test_axial_filtering()