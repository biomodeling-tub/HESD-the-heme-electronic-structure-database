#!/usr/bin/env python3
"""
Analyze the axial ligand encoding to understand why no swaps occurred.
"""

import pandas as pd

# Decoding dictionaries from preprocessor.py
axial1_decode = {0: 'CYS', 1: 'HIS'}
axial2_decode = {0: 'HIS', 1: 'HOH', 2: 'MET', 3: 'OXY'}

print("="*80)
print("AXIAL LIGAND ENCODING ANALYSIS")
print("="*80)

print("\nAxial1 encoding:")
for code, ligand in axial1_decode.items():
    print(f"  {code} → {ligand}")

print("\nAxial2 encoding:")
for code, ligand in axial2_decode.items():
    print(f"  {code} → {ligand}")

print("\n" + "="*80)
print("POSSIBLE ENCODED COMBINATIONS")
print("="*80)

print("\nAll possible combinations based on encoding:")
print(f"{'axial1 (code)':<20} {'axial2 (code)':<20} {'Combination':<20} {'Alphabetically Sorted?'}")
print("-"*90)

combinations_possible = []
for code1, ligand1 in axial1_decode.items():
    for code2, ligand2 in axial2_decode.items():
        combo = f"{ligand1}-{ligand2}"
        sorted_combo = '-'.join(sorted([ligand1, ligand2]))
        is_sorted = "✓ YES" if combo == sorted_combo else "✗ NO (would be " + sorted_combo + ")"

        combinations_possible.append({
            'axial1_code': code1,
            'axial2_code': code2,
            'ligand1': ligand1,
            'ligand2': ligand2,
            'combination': combo,
            'sorted_combination': sorted_combo,
            'is_sorted': combo == sorted_combo
        })

        print(f"{ligand1} ({code1}){' '*15} {ligand2} ({code2}){' '*15} {combo:<20} {is_sorted}")

print("\n" + "="*80)
print("ANALYSIS SUMMARY")
print("="*80)

num_possible = len(combinations_possible)
num_already_sorted = sum(1 for c in combinations_possible if c['is_sorted'])
num_needs_swap = num_possible - num_already_sorted

print(f"\nTotal possible combinations from encoding: {num_possible}")
print(f"Already in alphabetical order: {num_already_sorted} ({num_already_sorted/num_possible*100:.1f}%)")
print(f"Would need swapping: {num_needs_swap} ({num_needs_swap/num_possible*100:.1f}%)")

if num_needs_swap == 0:
    print("\n" + "⚠"*40)
    print("IMPORTANT FINDING:")
    print("The encoding scheme ALREADY enforces alphabetical ordering!")
    print("This is because:")
    print("  - axial1 can only be: CYS or HIS")
    print("  - axial2 can only be: HIS, HOH, MET, or OXY")
    print("  - Alphabetically: CYS < HIS < HOH < MET < OXY")
    print("  - Therefore, axial1 will ALWAYS be ≤ axial2 alphabetically")
    print("\nThis explains why:")
    print("  ✓ No swaps were detected in the test")
    print("  ✓ The plots haven't changed")
    print("  ✓ All combinations are already sorted")
    print("⚠"*40)
else:
    print(f"\n{'='*80}")
    print("COMBINATIONS THAT WOULD NEED SWAPPING:")
    print(f"{'='*80}")
    for c in combinations_possible:
        if not c['is_sorted']:
            print(f"{c['combination']} → {c['sorted_combination']}")

# Load actual data to verify
try:
    df = pd.read_csv("tables/processed_output.csv")
    if 'axial1' in df.columns and 'axial2' in df.columns:
        print(f"\n{'='*80}")
        print("ACTUAL DATA VERIFICATION")
        print(f"{'='*80}")

        # Check all combinations in the data
        actual_combinations = set()
        for _, row in df.iterrows():
            if pd.notna(row['axial1']) and pd.notna(row['axial2']):
                code1 = int(row['axial1'])
                code2 = int(row['axial2'])
                ligand1 = axial1_decode.get(code1, f'Unknown({code1})')
                ligand2 = axial2_decode.get(code2, f'Unknown({code2})')
                actual_combinations.add((code1, code2, f"{ligand1}-{ligand2}"))

        print(f"\nActual combinations found in data:")
        print(f"{'axial1':<10} {'axial2':<10} {'Combination':<20} {'Sorted?'}")
        print("-"*60)

        for code1, code2, combo in sorted(actual_combinations):
            ligands = combo.split('-')
            is_sorted = "✓" if len(ligands) == 2 and ligands[0] <= ligands[1] else "✗"
            print(f"{code1:<10} {code2:<10} {combo:<20} {is_sorted}")

except Exception as e:
    print(f"\nCould not load data: {e}")

print("\n" + "="*80)
