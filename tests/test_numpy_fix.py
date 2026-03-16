#!/usr/bin/env python3
"""
Simple test to verify the NumPy array boolean evaluation fix.
"""

import numpy as np

def test_numpy_fix():
    """Test the fixed NumPy array evaluation logic."""
    
    print("Testing NumPy Array Boolean Evaluation Fix")
    print("=" * 50)
    
    # Test cases that previously caused the ValueError
    test_cases = [
        None,                          # None case
        np.array([]),                  # Empty array
        np.array(['HIS-HIS', 'CYS-HOH']),  # Non-empty array
        []                             # Empty list
    ]
    
    for i, axial_combinations in enumerate(test_cases):
        print(f"Test case {i+1}: {type(axial_combinations).__name__} = {axial_combinations}")
        
        try:
            # This is the fixed logic from the function
            result = axial_combinations if axial_combinations is not None and len(axial_combinations) > 0 else []
            
            # Test the iteration that was causing the error
            combo_list = []
            for combo in sorted(result):
                combo_list.append(combo)
            
            print(f"  ✓ Success: {combo_list}")
            
        except Exception as e:
            print(f"  ✗ Error: {e}")
        
        print()
    
    print("Key fix:")
    print("- Replaced 'axial_combinations or []' with explicit null and length checks")
    print("- This prevents NumPy arrays from being evaluated in boolean context")
    print("- The fix handles None, empty arrays, and populated arrays correctly")

if __name__ == "__main__":
    test_numpy_fix()