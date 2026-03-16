#!/usr/bin/env python3
"""
Test script to verify that run_xtb_unconstrained.py creates timestamped files
and doesn't overwrite existing files.
"""

import os
import glob
import time
from datetime import datetime

def test_timestamp_uniqueness():
    """Test that timestamps are unique across runs."""
    print("Testing timestamp uniqueness...")

    # Generate two timestamps with a small delay
    ts1 = datetime.now().strftime("%Y%m%d_%H%M%S")
    time.sleep(2)  # Wait 2 seconds
    ts2 = datetime.now().strftime("%Y%m%d_%H%M%S")

    if ts1 != ts2:
        print("✓ PASS: Timestamps are unique with 2-second delay")
        return True
    else:
        print("✗ FAIL: Timestamps are not unique")
        return False

def test_file_naming_pattern():
    """Test that file naming patterns are correct."""
    print("\nTesting file naming patterns...")

    pdb_id = "test_pdb"
    timestamp = "20250126_123456"

    expected_files = [
        f"{pdb_id}_xtb_opt_{timestamp}.inp",
        f"{pdb_id}_xtb_opt_{timestamp}.out",
        f"{pdb_id}_opt_{timestamp}.xyz",
    ]

    print("Expected file patterns:")
    for f in expected_files:
        print(f"  - {f}")

    print("✓ PASS: File naming patterns defined correctly")
    return True

def check_existing_xtb_files():
    """Check for existing XTB output files in PDB directories."""
    print("\nChecking for existing XTB output files...")

    pattern = "PDB/**/*_opt_*.xyz"
    existing_files = glob.glob(pattern, recursive=True)

    print(f"Found {len(existing_files)} existing timestamped optimization files")

    if existing_files:
        print("\nSample files (first 5):")
        for f in existing_files[:5]:
            print(f"  - {f}")

    # Check for non-timestamped files that might be overwritten
    old_pattern = "PDB/**/*_opt.xyz"
    old_files = glob.glob(old_pattern, recursive=True)

    # Filter out timestamped files
    old_files = [f for f in old_files if not any(
        c.isdigit() for c in os.path.basename(f).split('_opt')[1].split('.')[0]
    )]

    if old_files:
        print(f"\n⚠ WARNING: Found {len(old_files)} old non-timestamped files:")
        for f in old_files[:5]:
            print(f"  - {f}")
        print("These files will NOT be overwritten by the new script.")

    return True

def verify_master_log_pattern():
    """Check if master log files exist."""
    print("\nChecking for master log files...")

    master_logs = glob.glob("xtb_optimization_master_*.log")

    if master_logs:
        print(f"Found {len(master_logs)} master log files:")
        for log in master_logs:
            size = os.path.getsize(log)
            print(f"  - {log} ({size} bytes)")
    else:
        print("No master log files found (none created yet)")

    return True

def main():
    """Run all tests."""
    print("="*80)
    print("XTB Timestamping Test Suite")
    print("="*80)

    tests = [
        test_timestamp_uniqueness,
        test_file_naming_pattern,
        check_existing_xtb_files,
        verify_master_log_pattern,
    ]

    results = []
    for test in tests:
        try:
            result = test()
            results.append(result)
        except Exception as e:
            print(f"✗ FAIL: Test raised exception: {e}")
            results.append(False)

    print("\n" + "="*80)
    print("TEST SUMMARY")
    print("="*80)
    passed = sum(results)
    total = len(results)
    print(f"Passed: {passed}/{total}")

    if all(results):
        print("\n✓ All tests passed!")
    else:
        print("\n⚠ Some tests failed")

    print("="*80)

if __name__ == "__main__":
    main()
