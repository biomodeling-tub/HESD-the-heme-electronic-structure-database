#!/usr/bin/env python3
"""
Run XTB optimization on all heme structures without constraints.

This script finds all {pdb-id}_system_protonated.pdb files in the PDB directory
and its subdirectories, converts them to XYZ format if needed, and runs XTB
with full geometry optimization at the highest reasonable level of theory,
without any constraints.

RESUMABLE: The script automatically skips structures that have already been
optimized, allowing you to stop and restart the process without losing progress.

Optimization Strategy (with automatic fallback):
------------------------------------------------
1. Attempt vtight optimization (highest quality)
2. If fails, attempt crude optimization (lower quality, more robust)
3. If crude succeeds, refine with vtight from crude structure
4. Track which optimization level succeeded for each structure

Resume Capability:
------------------
- By default, the script checks for existing {pdb-id}_opt_*.xyz files
- If found, that structure is skipped (shown as "SKIPPED" in logs)
- This allows stopping overnight and resuming in the morning
- To force reprocessing: modify skip_existing=False in process_all_structures()

Output files (all timestamped to prevent overwriting):
----------------------------------------------------
For each structure in PDB/{pdb-id}/:
  - {pdb-id}_xtb_opt_{timestamp}.inp       : XTB input file
  - {pdb-id}_xtb_opt_{timestamp}.out       : Complete XTB log (all attempts)
  - {pdb-id}_opt_{timestamp}.xyz           : Optimized structure (checked for resume)
  - {pdb-id}_xtb_errors_{timestamp}.log    : Error log (if failures occurred)
  - {pdb-id}_xtbopt.log_{timestamp}.log    : XTB optimization convergence log
  - {pdb-id}_charges_{timestamp}           : Partial charges
  - {pdb-id}_wbo_{timestamp}               : Wiberg bond orders

In working directory:
  - xtb_optimization_master_{timestamp}.log   : Master log tracking all structures
  - xtb_optimization_history_{timestamp}.txt  : Optimization level history
  - xtb_optimization_stats_{timestamp}.txt    : Summary statistics

All files include timestamps (YYYYMMDD_HHMMSS) to preserve previous runs.
"""

import os
import subprocess
import glob
import argparse
from pathlib import Path
from datetime import datetime

def parse_pdb_id_list(pdb_list_arg):
    """
    Parse PDB-ID list from command-line argument.

    Parameters:
    -----------
    pdb_list_arg : str
        Either comma-separated PDB-IDs or path to a file containing PDB-IDs

    Returns:
    --------
    list : List of PDB-IDs
    """
    if not pdb_list_arg:
        return []

    # Check if it's a file
    if os.path.isfile(pdb_list_arg):
        with open(pdb_list_arg, 'r') as f:
            pdb_ids = [line.strip() for line in f if line.strip()]
    else:
        # Assume comma-separated list
        pdb_ids = [pid.strip() for pid in pdb_list_arg.split(',') if pid.strip()]

    return pdb_ids

def check_already_processed(pdb_id, output_dir):
    """
    Check if a structure has already been successfully optimized.

    Parameters:
    -----------
    pdb_id : str
        PDB ID of the structure
    output_dir : str
        Directory containing potential output files

    Returns:
    --------
    tuple : (is_processed: bool, existing_files: list)
        Returns True if optimized structure exists, along with list of found files
    """
    # Look for any optimized structure file with pattern {pdb_id}_opt_*.xyz
    pattern = os.path.join(output_dir, f"{pdb_id}_opt_*.xyz")
    existing_opt_files = glob.glob(pattern)

    if existing_opt_files:
        return True, existing_opt_files

    return False, []

def find_protonated_structures(base_dir="PDB"):
    """
    Find all {pdb-id}_system_protonated.pdb files in the directory tree.

    Parameters:
    -----------
    base_dir : str
        Base directory to search (default: "PDB")

    Returns:
    --------
    list : List of paths to protonated structure files
    """
    pattern = os.path.join(base_dir, "**", "*_system_protonated.pdb")
    pdb_files = glob.glob(pattern, recursive=True)
    print(f"Found {len(pdb_files)} protonated structure files")
    return sorted(pdb_files)

def pdb_to_xyz(pdb_file, xyz_file):
    """
    Convert PDB file to XYZ format using Open Babel.

    Parameters:
    -----------
    pdb_file : str
        Path to input PDB file
    xyz_file : str
        Path to output XYZ file

    Returns:
    --------
    bool : True if conversion successful, False otherwise
    """
    try:
        cmd = ['obabel', pdb_file, '-O', xyz_file]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        return os.path.exists(xyz_file)
    except subprocess.CalledProcessError as e:
        print(f"Error converting {pdb_file} to XYZ: {e}")
        return False
    except FileNotFoundError:
        print("Open Babel not found. Attempting to use XTB directly with PDB file.")
        return False

def create_xtb_input_unconstrained(pdb_id, output_dir, timestamp):
    """
    Create XTB input file for unconstrained optimization.

    Parameters:
    -----------
    pdb_id : str
        PDB ID of the structure
    output_dir : str
        Directory to save the input file
    timestamp : str
        Timestamp string for unique file naming

    Returns:
    --------
    str : Path to the created input file
    """
    input_file = os.path.join(output_dir, f"{pdb_id}_xtb_opt_{timestamp}.inp")

    # Create input file with no constraints, just optimization settings
    with open(input_file, 'w') as f:
        f.write("$opt\n")
        f.write("   optlevel=2\n")  # 2 = vtight, highest reasonable level
        f.write("$end\n")

    return input_file

def create_xtb_input_with_charge_mult(pdb_id, output_dir, timestamp, charge, multiplicity):
    """
    Create XTB input file for unconstrained optimization with specific charge and multiplicity.

    Parameters:
    -----------
    pdb_id : str
        PDB ID of the structure
    output_dir : str
        Directory to save the input file
    timestamp : str
        Timestamp string for unique file naming
    charge : int
        Molecular charge
    multiplicity : int
        Spin multiplicity (2S+1)

    Returns:
    --------
    str : Path to the created input file
    """
    input_file = os.path.join(output_dir, f"{pdb_id}_xtb_opt_q{charge}_m{multiplicity}_{timestamp}.inp")

    # Calculate number of unpaired electrons (spin state)
    unpaired_electrons = multiplicity - 1

    # Create input file with charge, multiplicity, and optimization settings
    with open(input_file, 'w') as f:
        f.write(f"$chrg {charge}\n")
        f.write(f"$spin {unpaired_electrons}\n")
        f.write("$opt\n")
        f.write("   optlevel=2\n")  # 2 = vtight, highest reasonable level
        f.write("$end\n")

    return input_file

def check_already_processed_charge_mult(pdb_id, output_dir, charge, multiplicity):
    """
    Check if a structure has already been optimized for a specific charge-multiplicity state.

    Parameters:
    -----------
    pdb_id : str
        PDB ID of the structure
    output_dir : str
        Directory containing potential output files
    charge : int
        Molecular charge
    multiplicity : int
        Spin multiplicity

    Returns:
    --------
    tuple : (is_processed: bool, existing_files: list)
        Returns True if optimized structure exists, along with list of found files
    """
    # Look for any optimized structure file with pattern {pdb_id}_opt_q{charge}_m{multiplicity}_*.xyz
    pattern = os.path.join(output_dir, f"{pdb_id}_opt_q{charge}_m{multiplicity}_*.xyz")
    existing_opt_files = glob.glob(pattern)

    if existing_opt_files:
        return True, existing_opt_files

    return False, []

def run_xtb_optimization(structure_file, pdb_id, output_dir, timestamp, use_xyz=True):
    """
    Run XTB optimization without constraints with fallback strategy.

    Strategy:
    1. Try vtight optimization
    2. If fails, try crude optimization
    3. If crude succeeds, refine with vtight from crude structure

    Parameters:
    -----------
    structure_file : str
        Path to structure file (XYZ or PDB)
    pdb_id : str
        PDB ID of the structure
    output_dir : str
        Directory containing the structure and for output
    timestamp : str
        Timestamp string for unique file naming
    use_xyz : bool
        Whether the structure file is in XYZ format

    Returns:
    --------
    tuple : (success: bool, output_file: str, log_file: str, opt_level: str)
    """
    # Create input file
    input_file = create_xtb_input_unconstrained(pdb_id, output_dir, timestamp)

    # Define output files with timestamp to prevent overwriting
    log_file = os.path.join(output_dir, f"{pdb_id}_xtb_opt_{timestamp}.out")
    error_file = os.path.join(output_dir, f"{pdb_id}_xtb_errors_{timestamp}.log")
    opt_structure = os.path.join(output_dir, f"{pdb_id}_opt_{timestamp}.xyz")

    # Change to the structure directory
    original_dir = os.getcwd()
    os.chdir(output_dir)

    # After changing directory, use only basenames for file operations
    log_file_basename = os.path.basename(log_file)
    error_file_basename = os.path.basename(error_file)
    opt_structure_basename = os.path.basename(opt_structure)

    success = False
    opt_level = "failed"
    structure_basename = os.path.basename(structure_file)

    try:
        # ========== ATTEMPT 1: VTIGHT OPTIMIZATION ==========
        print(f"Running XTB vtight optimization for {pdb_id}...")
        cmd = [
            'xtb',
            structure_basename,
            '--opt', 'vtight',
            '--gfn', '2',
            '--input', os.path.basename(input_file)
        ]
        print(f"Command: {' '.join(cmd)}")

        # Run XTB with vtight
        process = subprocess.Popen(
            ' '.join(cmd),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate()

        # Save output to log file
        with open(log_file_basename, 'w') as log:
            log.write("=== VTIGHT OPTIMIZATION ATTEMPT ===\n")
            log.write(f"Command: {' '.join(cmd)}\n")
            log.write("="*80 + "\n\n")
            log.write(stdout.decode())
            if stderr:
                log.write("\n=== STDERR ===\n")
                log.write(stderr.decode())

        # Check if vtight optimization succeeded
        if process.returncode == 0 and os.path.exists('xtbopt.xyz'):
            os.rename('xtbopt.xyz', os.path.basename(opt_structure))
            print(f"✓ Successfully optimized {pdb_id} with vtight")
            success = True
            opt_level = "vtight"
        else:
            # ========== ATTEMPT 2: CRUDE OPTIMIZATION ==========
            print(f"⚠ vtight failed, trying crude optimization for {pdb_id}...")

            # Log the vtight failure
            with open(error_file_basename, 'w') as err:
                err.write(f"vtight optimization failed with return code {process.returncode}\n")
                err.write("="*80 + "\n")
                err.write("STDOUT:\n")
                err.write(stdout.decode())
                err.write("\n" + "="*80 + "\n")
                err.write("STDERR:\n")
                err.write(stderr.decode())
                err.write("\n" + "="*80 + "\n\n")

            # Try crude optimization
            cmd_crude = [
                'xtb',
                structure_basename,
                '--opt', 'crude',
                '--gfn', '2',
                '--input', os.path.basename(input_file)
            ]
            print(f"Command: {' '.join(cmd_crude)}")

            process_crude = subprocess.Popen(
                ' '.join(cmd_crude),
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            stdout_crude, stderr_crude = process_crude.communicate()

            # Append crude output to log file
            with open(log_file_basename, 'a') as log:
                log.write("\n\n" + "="*80 + "\n")
                log.write("=== CRUDE OPTIMIZATION ATTEMPT ===\n")
                log.write(f"Command: {' '.join(cmd_crude)}\n")
                log.write("="*80 + "\n\n")
                log.write(stdout_crude.decode())
                if stderr_crude:
                    log.write("\n=== STDERR ===\n")
                    log.write(stderr_crude.decode())

            if process_crude.returncode == 0 and os.path.exists('xtbopt.xyz'):
                print(f"✓ Crude optimization succeeded for {pdb_id}")

                # ========== ATTEMPT 3: REFINE WITH VTIGHT ==========
                print(f"Refining with vtight from crude structure for {pdb_id}...")

                cmd_refine = [
                    'xtb',
                    'xtbopt.xyz',
                    '--opt', 'vtight',
                    '--gfn', '2',
                    '--input', os.path.basename(input_file)
                ]
                print(f"Command: {' '.join(cmd_refine)}")

                process_refine = subprocess.Popen(
                    ' '.join(cmd_refine),
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )
                stdout_refine, stderr_refine = process_refine.communicate()

                # Append refinement output to log file
                with open(log_file_basename, 'a') as log:
                    log.write("\n\n" + "="*80 + "\n")
                    log.write("=== VTIGHT REFINEMENT FROM CRUDE ===\n")
                    log.write(f"Command: {' '.join(cmd_refine)}\n")
                    log.write("="*80 + "\n\n")
                    log.write(stdout_refine.decode())
                    if stderr_refine:
                        log.write("\n=== STDERR ===\n")
                        log.write(stderr_refine.decode())

                if process_refine.returncode == 0 and os.path.exists('xtbopt.xyz'):
                    os.rename('xtbopt.xyz', os.path.basename(opt_structure))
                    print(f"✓ Successfully refined {pdb_id} with vtight from crude")
                    success = True
                    opt_level = "crude+vtight"
                else:
                    # Refinement failed, keep crude result
                    if os.path.exists('xtbopt.xyz'):
                        os.rename('xtbopt.xyz', os.path.basename(opt_structure))
                    print(f"⚠ Refinement failed, keeping crude optimization for {pdb_id}")
                    success = True
                    opt_level = "crude"

                    # Log refinement failure
                    with open(error_file_basename, 'a') as err:
                        err.write(f"vtight refinement failed with return code {process_refine.returncode}\n")
                        err.write("="*80 + "\n")
                        err.write("STDOUT:\n")
                        err.write(stdout_refine.decode())
                        err.write("\n" + "="*80 + "\n")
                        err.write("STDERR:\n")
                        err.write(stderr_refine.decode())
            else:
                # Crude optimization also failed
                print(f"✗ All optimization attempts failed for {pdb_id}")
                success = False
                opt_level = "failed"

                # Log crude failure
                with open(error_file_basename, 'a') as err:
                    err.write(f"crude optimization failed with return code {process_crude.returncode}\n")
                    err.write("="*80 + "\n")
                    err.write("STDOUT:\n")
                    err.write(stdout_crude.decode())
                    err.write("\n" + "="*80 + "\n")
                    err.write("STDERR:\n")
                    err.write(stderr_crude.decode())

    except Exception as e:
        print(f"✗ Error running XTB for {pdb_id}: {e}")
        success = False
        opt_structure = None
        opt_level = "error"

        # Log the exception
        with open(error_file_basename, 'a') as err:
            err.write(f"\nException occurred: {e}\n")

    finally:
        # Clean up temporary XTB files and rename important ones with timestamp
        temp_files = ['xtbrestart', 'xtbtopo.mol', 'wbo', 'charges',
                      'gfnff_topo', 'xtbopt.log', '.xtboptok']

        for temp_file in temp_files:
            # Files are in current directory (already changed to output_dir)
            if os.path.exists(temp_file):
                # Archive important files with timestamp, delete others
                if temp_file in ['xtbopt.log', 'charges', 'wbo']:
                    archived_name = f"{pdb_id}_{temp_file}_{timestamp}"
                    if temp_file == 'xtbopt.log':
                        archived_name += '.log'
                    try:
                        os.rename(temp_file, archived_name)
                    except:
                        pass  # If rename fails, file will be left as-is
                else:
                    try:
                        os.remove(temp_file)
                    except:
                        pass  # If deletion fails, file will be left as-is

        # Return to original directory
        os.chdir(original_dir)

    # Print summary of what happened
    if success:
        print(f"  Output: {opt_structure}")
        print(f"  Log: {log_file}")
        print(f"  Optimization level: {opt_level}")

    return success, opt_structure, log_file, opt_level

def run_xtb_charge_mult_exploration(structure_file, pdb_id, output_dir, timestamp, charge, multiplicity, use_xyz=True):
    """
    Run XTB optimization for a specific charge-multiplicity state.

    Parameters:
    -----------
    structure_file : str
        Path to structure file (XYZ or PDB)
    pdb_id : str
        PDB ID of the structure
    output_dir : str
        Directory containing the structure and for output
    timestamp : str
        Timestamp string for unique file naming
    charge : int
        Molecular charge
    multiplicity : int
        Spin multiplicity
    use_xyz : bool
        Whether the structure file is in XYZ format

    Returns:
    --------
    tuple : (success: bool, output_file: str, log_file: str, opt_level: str)
    """
    # Create input file with charge and multiplicity
    input_file = create_xtb_input_with_charge_mult(pdb_id, output_dir, timestamp, charge, multiplicity)

    # Define output files with timestamp and charge-mult labels
    log_file = os.path.join(output_dir, f"{pdb_id}_xtb_opt_q{charge}_m{multiplicity}_{timestamp}.out")
    error_file = os.path.join(output_dir, f"{pdb_id}_xtb_errors_q{charge}_m{multiplicity}_{timestamp}.log")
    opt_structure = os.path.join(output_dir, f"{pdb_id}_opt_q{charge}_m{multiplicity}_{timestamp}.xyz")

    # Change to the structure directory
    original_dir = os.getcwd()
    os.chdir(output_dir)

    # After changing directory, use only basenames for file operations
    log_file_basename = os.path.basename(log_file)
    error_file_basename = os.path.basename(error_file)
    opt_structure_basename = os.path.basename(opt_structure)

    success = False
    opt_level = "failed"
    structure_basename = os.path.basename(structure_file)

    try:
        # Run XTB with vtight optimization for this charge-multiplicity state
        print(f"Running XTB vtight optimization for {pdb_id} (q={charge}, m={multiplicity})...")
        cmd = [
            'xtb',
            structure_basename,
            '--opt', 'vtight',
            '--gfn', '2',
            '--input', os.path.basename(input_file)
        ]
        print(f"Command: {' '.join(cmd)}")

        # Run XTB
        process = subprocess.Popen(
            ' '.join(cmd),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate()

        # Save output to log file
        with open(log_file_basename, 'w') as log:
            log.write(f"=== VTIGHT OPTIMIZATION (q={charge}, m={multiplicity}) ===\n")
            log.write(f"Command: {' '.join(cmd)}\n")
            log.write("="*80 + "\n\n")
            log.write(stdout.decode())
            if stderr:
                log.write("\n=== STDERR ===\n")
                log.write(stderr.decode())

        # Check if optimization succeeded
        if process.returncode == 0 and os.path.exists('xtbopt.xyz'):
            os.rename('xtbopt.xyz', opt_structure_basename)
            print(f"✓ Successfully optimized {pdb_id} (q={charge}, m={multiplicity})")
            success = True
            opt_level = "vtight"
        else:
            print(f"✗ Optimization failed for {pdb_id} (q={charge}, m={multiplicity})")
            success = False
            opt_level = "failed"

            # Log the failure
            with open(error_file_basename, 'w') as err:
                err.write(f"vtight optimization failed with return code {process.returncode}\n")
                err.write(f"Charge: {charge}, Multiplicity: {multiplicity}\n")
                err.write("="*80 + "\n")
                err.write("STDOUT:\n")
                err.write(stdout.decode())
                err.write("\n" + "="*80 + "\n")
                err.write("STDERR:\n")
                err.write(stderr.decode())

    except Exception as e:
        print(f"✗ Error running XTB for {pdb_id} (q={charge}, m={multiplicity}): {e}")
        success = False
        opt_structure = None
        opt_level = "error"

        # Log the exception
        with open(error_file_basename, 'w') as err:
            err.write(f"Exception occurred: {e}\n")
            err.write(f"Charge: {charge}, Multiplicity: {multiplicity}\n")

    finally:
        # Clean up temporary XTB files and rename important ones with timestamp
        temp_files = ['xtbrestart', 'xtbtopo.mol', 'wbo', 'charges',
                      'gfnff_topo', 'xtbopt.log', '.xtboptok']

        for temp_file in temp_files:
            # Files are in current directory (already changed to output_dir)
            if os.path.exists(temp_file):
                # Archive important files with timestamp, delete others
                if temp_file in ['xtbopt.log', 'charges', 'wbo']:
                    archived_name = f"{pdb_id}_{temp_file}_q{charge}_m{multiplicity}_{timestamp}"
                    if temp_file == 'xtbopt.log':
                        archived_name += '.log'
                    try:
                        os.rename(temp_file, archived_name)
                    except:
                        pass  # If rename fails, file will be left as-is
                else:
                    try:
                        os.remove(temp_file)
                    except:
                        pass  # If deletion fails, file will be left as-is

        # Return to original directory
        os.chdir(original_dir)

    # Print summary of what happened
    if success:
        print(f"  Output: {opt_structure}")
        print(f"  Log: {log_file}")
        print(f"  Optimization level: {opt_level}")

    return success, opt_structure, log_file, opt_level

def process_all_structures(base_dir="PDB", convert_to_xyz=True, skip_existing=True):
    """
    Process all protonated structures with XTB optimization.

    Parameters:
    -----------
    base_dir : str
        Base directory containing PDB subdirectories
    convert_to_xyz : bool
        Whether to convert PDB to XYZ before running XTB
    skip_existing : bool
        If True, skip structures that have already been optimized (default: True)
        Set to False to force reprocessing of all structures

    Returns:
    --------
    dict : Statistics about the processing
    """
    # Generate timestamp for this run (all files will use the same timestamp)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Create master log file
    master_log = f"xtb_optimization_master_{timestamp}.log"

    # Find all protonated structures
    pdb_files = find_protonated_structures(base_dir)

    if not pdb_files:
        print("No protonated structure files found!")
        with open(master_log, 'w') as log:
            log.write(f"XTB Optimization Run - {timestamp}\n")
            log.write("="*80 + "\n")
            log.write("No protonated structure files found!\n")
        return {'total': 0, 'success': 0, 'failed': 0}

    # Statistics
    stats = {
        'total': len(pdb_files),
        'success': 0,
        'failed': 0,
        'skipped': 0,
        'vtight': 0,
        'crude+vtight': 0,
        'crude': 0,
        'timestamp': timestamp
    }

    # Create history file to track optimization levels
    history_file = f"xtb_optimization_history_{timestamp}.txt"

    # Write master log header
    with open(master_log, 'w') as log:
        log.write(f"XTB Optimization Run - {timestamp}\n")
        log.write("="*80 + "\n")
        log.write(f"Found {len(pdb_files)} protonated structure files\n")
        log.write("="*80 + "\n\n")

    # Write history file header
    with open(history_file, 'w') as hist:
        hist.write(f"XTB Optimization History - {timestamp}\n")
        hist.write("="*80 + "\n")
        hist.write(f"{'PDB ID':<15} {'Optimization Level':<20} {'Status'}\n")
        hist.write("="*80 + "\n")

    # Process each structure
    for i, pdb_file in enumerate(pdb_files, 1):
        print(f"\n{'='*80}")
        print(f"Processing {i}/{len(pdb_files)}: {pdb_file}")
        print(f"{'='*80}")

        # Extract PDB ID and directory
        pdb_path = Path(pdb_file)
        pdb_id = pdb_path.stem.replace('_system_protonated', '')
        output_dir = pdb_path.parent

        # Log to master log
        with open(master_log, 'a') as log:
            log.write(f"\n{'='*80}\n")
            log.write(f"Processing {i}/{len(pdb_files)}: {pdb_id}\n")
            log.write(f"Source: {pdb_file}\n")
            log.write(f"{'='*80}\n")

        # Check if already processed
        if skip_existing:
            already_done, existing_files = check_already_processed(pdb_id, output_dir)
            if already_done:
                print(f"⊙ SKIPPING {pdb_id}: Already optimized")
                print(f"  Existing file(s): {', '.join([os.path.basename(f) for f in existing_files])}")
                stats['skipped'] += 1

                # Log to master log
                with open(master_log, 'a') as log:
                    log.write(f"⊙ SKIPPED: Already optimized\n")
                    log.write(f"  Existing file(s): {', '.join([os.path.basename(f) for f in existing_files])}\n")

                # Log to history file
                with open(history_file, 'a') as hist:
                    hist.write(f"{pdb_id:<15} {'existing':<20} SKIPPED\n")

                continue

        # Decide whether to convert to XYZ
        if convert_to_xyz:
            xyz_file = os.path.join(output_dir, f"{pdb_id}_system_protonated.xyz")

            # Check if XYZ already exists
            if not os.path.exists(xyz_file):
                print(f"Converting {pdb_id} to XYZ format...")
                with open(master_log, 'a') as log:
                    log.write(f"Converting to XYZ format...\n")
                if not pdb_to_xyz(pdb_file, xyz_file):
                    print(f"Using PDB file directly for {pdb_id}")
                    with open(master_log, 'a') as log:
                        log.write(f"XYZ conversion failed, using PDB file directly\n")
                    structure_file = pdb_file
                    use_xyz = False
                else:
                    structure_file = xyz_file
                    use_xyz = True
                    with open(master_log, 'a') as log:
                        log.write(f"Successfully converted to XYZ: {xyz_file}\n")
            else:
                print(f"Using existing XYZ file for {pdb_id}")
                with open(master_log, 'a') as log:
                    log.write(f"Using existing XYZ file: {xyz_file}\n")
                structure_file = xyz_file
                use_xyz = True
        else:
            structure_file = pdb_file
            use_xyz = False
            with open(master_log, 'a') as log:
                log.write(f"Using PDB file directly: {pdb_file}\n")

        # Run XTB optimization
        success, opt_structure, log_file, opt_level = run_xtb_optimization(
            structure_file, pdb_id, output_dir, timestamp, use_xyz
        )

        # Update statistics
        if success:
            stats['success'] += 1
            if opt_level in stats:
                stats[opt_level] += 1
        else:
            stats['failed'] += 1

        # Log results to master log
        with open(master_log, 'a') as log:
            if success:
                log.write(f"✓ SUCCESS: Optimized structure saved to {opt_structure}\n")
                log.write(f"  XTB log: {log_file}\n")
                log.write(f"  Optimization level: {opt_level}\n")
            else:
                log.write(f"✗ FAILED: Optimization failed for {pdb_id}\n")
                log.write(f"  Check log: {log_file}\n")
                log.write(f"  Final status: {opt_level}\n")

        # Log to history file
        with open(history_file, 'a') as hist:
            status = "SUCCESS" if success else "FAILED"
            hist.write(f"{pdb_id:<15} {opt_level:<20} {status}\n")

    # Print summary
    print(f"\n{'='*80}")
    print("PROCESSING COMPLETE")
    print(f"{'='*80}")
    print(f"Total structures: {stats['total']}")
    print(f"Skipped (already done): {stats['skipped']}")
    print(f"Successfully optimized: {stats['success']}")
    print(f"  - vtight: {stats['vtight']}")
    print(f"  - crude+vtight: {stats['crude+vtight']}")
    print(f"  - crude only: {stats['crude']}")
    print(f"Failed: {stats['failed']}")
    print(f"Run timestamp: {timestamp}")
    print(f"Master log: {master_log}")
    print(f"History file: {history_file}")
    print(f"{'='*80}\n")

    # Write summary to master log
    with open(master_log, 'a') as log:
        log.write(f"\n{'='*80}\n")
        log.write("PROCESSING COMPLETE\n")
        log.write(f"{'='*80}\n")
        log.write(f"Total structures: {stats['total']}\n")
        log.write(f"Skipped (already done): {stats['skipped']}\n")
        log.write(f"Successfully optimized: {stats['success']}\n")
        log.write(f"  - vtight: {stats['vtight']}\n")
        log.write(f"  - crude+vtight: {stats['crude+vtight']}\n")
        log.write(f"  - crude only: {stats['crude']}\n")
        log.write(f"Failed: {stats['failed']}\n")
        log.write(f"Timestamp: {timestamp}\n")
        log.write(f"Master log: {master_log}\n")
        log.write(f"History file: {history_file}\n")
        log.write(f"{'='*80}\n")

    return stats

def process_charge_mult_structures(pdb_id_list, base_dir="PDB", convert_to_xyz=True, skip_existing=True):
    """
    Process specific structures with charge-multiplicity exploration.

    For each structure, runs XTB optimization for multiple charge-multiplicity states:
    - q=0, m=1 (neutral, singlet)
    - q=0, m=5 (neutral, quintet)
    - q=1, m=2 (cation, doublet)
    - q=1, m=6 (cation, sextet)

    Parameters:
    -----------
    pdb_id_list : list
        List of PDB-IDs to process
    base_dir : str
        Base directory containing PDB subdirectories
    convert_to_xyz : bool
        Whether to convert PDB to XYZ before running XTB
    skip_existing : bool
        If True, skip charge-mult states that have already been optimized

    Returns:
    --------
    dict : Statistics about the processing
    """
    # Define charge-multiplicity combinations
    charge_mult_combinations = [
        (0, 1),  # neutral, singlet
        (0, 5),  # neutral, quintet
        (1, 2),  # cation, doublet
        (1, 6)   # cation, sextet
    ]

    # Generate timestamp for this run
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Create master log file
    master_log = f"xtb_chargemult_master_{timestamp}.log"

    # Statistics
    stats = {
        'total_structures': len(pdb_id_list),
        'total_calculations': len(pdb_id_list) * len(charge_mult_combinations),
        'q0_m1': {'success': 0, 'failed': 0, 'skipped': 0},
        'q0_m5': {'success': 0, 'failed': 0, 'skipped': 0},
        'q1_m2': {'success': 0, 'failed': 0, 'skipped': 0},
        'q1_m6': {'success': 0, 'failed': 0, 'skipped': 0},
        'timestamp': timestamp
    }

    # Create history file
    history_file = f"xtb_chargemult_history_{timestamp}.txt"

    # Write master log header
    with open(master_log, 'w') as log:
        log.write(f"XTB Charge-Multiplicity Exploration - {timestamp}\n")
        log.write("="*80 + "\n")
        log.write(f"Processing {len(pdb_id_list)} structures\n")
        log.write(f"Each structure will be optimized for {len(charge_mult_combinations)} charge-multiplicity states\n")
        log.write("Charge-Multiplicity combinations:\n")
        for q, m in charge_mult_combinations:
            log.write(f"  - q={q}, m={m}\n")
        log.write("="*80 + "\n\n")

    # Write history file header
    with open(history_file, 'w') as hist:
        hist.write(f"XTB Charge-Multiplicity History - {timestamp}\n")
        hist.write("="*80 + "\n")
        hist.write(f"{'PDB ID':<15} {'Charge':<10} {'Mult':<10} {'Status':<15} {'Opt Level'}\n")
        hist.write("="*80 + "\n")

    # Process each PDB-ID
    for i, pdb_id in enumerate(pdb_id_list, 1):
        print(f"\n{'='*80}")
        print(f"Processing {i}/{len(pdb_id_list)}: {pdb_id}")
        print(f"{'='*80}")

        # Find the protonated structure file for this PDB-ID
        pattern = os.path.join(base_dir, "**", f"{pdb_id}_system_protonated.pdb")
        pdb_files = glob.glob(pattern, recursive=True)

        if not pdb_files:
            print(f"⚠ WARNING: No protonated structure found for {pdb_id}")
            with open(master_log, 'a') as log:
                log.write(f"\n{'='*80}\n")
                log.write(f"WARNING: No protonated structure found for {pdb_id}\n")
                log.write(f"{'='*80}\n")
            continue

        if len(pdb_files) > 1:
            print(f"⚠ WARNING: Multiple protonated structures found for {pdb_id}, using first one")

        pdb_file = pdb_files[0]
        pdb_path = Path(pdb_file)
        output_dir = pdb_path.parent

        # Log to master log
        with open(master_log, 'a') as log:
            log.write(f"\n{'='*80}\n")
            log.write(f"Processing {i}/{len(pdb_id_list)}: {pdb_id}\n")
            log.write(f"Source: {pdb_file}\n")
            log.write(f"{'='*80}\n")

        # Prepare structure file (convert to XYZ if needed)
        if convert_to_xyz:
            xyz_file = os.path.join(output_dir, f"{pdb_id}_system_protonated.xyz")

            if not os.path.exists(xyz_file):
                print(f"Converting {pdb_id} to XYZ format...")
                with open(master_log, 'a') as log:
                    log.write(f"Converting to XYZ format...\n")
                if not pdb_to_xyz(pdb_file, xyz_file):
                    print(f"Using PDB file directly for {pdb_id}")
                    with open(master_log, 'a') as log:
                        log.write(f"XYZ conversion failed, using PDB file directly\n")
                    structure_file = pdb_file
                    use_xyz = False
                else:
                    structure_file = xyz_file
                    use_xyz = True
                    with open(master_log, 'a') as log:
                        log.write(f"Successfully converted to XYZ: {xyz_file}\n")
            else:
                print(f"Using existing XYZ file for {pdb_id}")
                with open(master_log, 'a') as log:
                    log.write(f"Using existing XYZ file: {xyz_file}\n")
                structure_file = xyz_file
                use_xyz = True
        else:
            structure_file = pdb_file
            use_xyz = False
            with open(master_log, 'a') as log:
                log.write(f"Using PDB file directly: {pdb_file}\n")

        # Process each charge-multiplicity combination
        for charge, multiplicity in charge_mult_combinations:
            print(f"\n--- Charge: {charge}, Multiplicity: {multiplicity} ---")

            # Check if already processed
            if skip_existing:
                already_done, existing_files = check_already_processed_charge_mult(
                    pdb_id, output_dir, charge, multiplicity
                )
                if already_done:
                    print(f"⊙ SKIPPING {pdb_id} (q={charge}, m={multiplicity}): Already optimized")
                    print(f"  Existing file(s): {', '.join([os.path.basename(f) for f in existing_files])}")

                    state_key = f"q{charge}_m{multiplicity}"
                    stats[state_key]['skipped'] += 1

                    # Log to master log
                    with open(master_log, 'a') as log:
                        log.write(f"⊙ SKIPPED (q={charge}, m={multiplicity}): Already optimized\n")
                        log.write(f"  Existing file(s): {', '.join([os.path.basename(f) for f in existing_files])}\n")

                    # Log to history file
                    with open(history_file, 'a') as hist:
                        hist.write(f"{pdb_id:<15} {charge:<10} {multiplicity:<10} {'SKIPPED':<15} existing\n")

                    continue

            # Run XTB optimization for this charge-multiplicity state
            success, opt_structure, log_file, opt_level = run_xtb_charge_mult_exploration(
                structure_file, pdb_id, output_dir, timestamp, charge, multiplicity, use_xyz
            )

            # Update statistics
            state_key = f"q{charge}_m{multiplicity}"
            if success:
                stats[state_key]['success'] += 1
            else:
                stats[state_key]['failed'] += 1

            # Log results to master log
            with open(master_log, 'a') as log:
                if success:
                    log.write(f"✓ SUCCESS (q={charge}, m={multiplicity}): {opt_structure}\n")
                    log.write(f"  XTB log: {log_file}\n")
                    log.write(f"  Optimization level: {opt_level}\n")
                else:
                    log.write(f"✗ FAILED (q={charge}, m={multiplicity}): {pdb_id}\n")
                    log.write(f"  Check log: {log_file}\n")
                    log.write(f"  Final status: {opt_level}\n")

            # Log to history file
            with open(history_file, 'a') as hist:
                status = "SUCCESS" if success else "FAILED"
                hist.write(f"{pdb_id:<15} {charge:<10} {multiplicity:<10} {status:<15} {opt_level}\n")

    # Print summary
    print(f"\n{'='*80}")
    print("CHARGE-MULTIPLICITY EXPLORATION COMPLETE")
    print(f"{'='*80}")
    print(f"Total structures: {stats['total_structures']}")
    print(f"Total calculations: {stats['total_calculations']}")
    print()
    for q, m in charge_mult_combinations:
        state_key = f"q{q}_m{m}"
        print(f"q={q}, m={m}:")
        print(f"  Success: {stats[state_key]['success']}")
        print(f"  Failed: {stats[state_key]['failed']}")
        print(f"  Skipped: {stats[state_key]['skipped']}")
    print(f"\nRun timestamp: {timestamp}")
    print(f"Master log: {master_log}")
    print(f"History file: {history_file}")
    print(f"{'='*80}\n")

    # Write summary to master log
    with open(master_log, 'a') as log:
        log.write(f"\n{'='*80}\n")
        log.write("CHARGE-MULTIPLICITY EXPLORATION COMPLETE\n")
        log.write(f"{'='*80}\n")
        log.write(f"Total structures: {stats['total_structures']}\n")
        log.write(f"Total calculations: {stats['total_calculations']}\n\n")
        for q, m in charge_mult_combinations:
            state_key = f"q{q}_m{m}"
            log.write(f"q={q}, m={m}:\n")
            log.write(f"  Success: {stats[state_key]['success']}\n")
            log.write(f"  Failed: {stats[state_key]['failed']}\n")
            log.write(f"  Skipped: {stats[state_key]['skipped']}\n")
        log.write(f"\nTimestamp: {timestamp}\n")
        log.write(f"Master log: {master_log}\n")
        log.write(f"History file: {history_file}\n")
        log.write(f"{'='*80}\n")

    return stats

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='XTB Optimization Pipeline with optional charge-multiplicity exploration',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Standard mode (all structures)
  python run_xtb_unconstrained.py

  # Charge-multiplicity mode with comma-separated PDB-IDs
  python run_xtb_unconstrained.py --mode charge-mult --pdb-list "1A6G,2HHB,3RGK"

  # Charge-multiplicity mode with file containing PDB-IDs
  python run_xtb_unconstrained.py --mode charge-mult --pdb-list pdb_ids.txt
        '''
    )
    parser.add_argument(
        '--mode',
        choices=['standard', 'charge-mult'],
        default='standard',
        help='Optimization mode: "standard" runs normal optimization, "charge-mult" explores multiple charge-multiplicity states'
    )
    parser.add_argument(
        '--pdb-list',
        type=str,
        help='For charge-mult mode: comma-separated PDB-IDs or path to file with PDB-IDs (one per line)'
    )

    args = parser.parse_args()

    # Validate arguments
    if args.mode == 'charge-mult' and not args.pdb_list:
        parser.error("--pdb-list is required when using --mode charge-mult")

    # Check if XTB is available
    try:
        subprocess.run(['xtb', '--version'], capture_output=True, check=True)
        print("✓ XTB found\n")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("✗ XTB not found in PATH!")
        print("Please install XTB before running this script.")
        exit(1)

    # Run appropriate mode
    if args.mode == 'standard':
        # Standard mode: process all structures
        print("XTB Unconstrained Optimization Pipeline")
        print("=" * 80)
        print("This script will run XTB with full geometry optimization")
        print("at the highest reasonable level of theory (vtight, GFN2-xTB)")
        print("without any constraints on all protonated heme structures.")
        print()
        print("RESUMABLE: Structures already optimized will be skipped automatically.")
        print("           Safe to stop and restart anytime!")
        print()
        print("Optimization Strategy:")
        print("  1. Try vtight optimization (highest quality)")
        print("  2. If fails, try crude optimization (more robust)")
        print("  3. If crude succeeds, refine with vtight")
        print()
        print("All output files are timestamped to prevent overwriting.")
        print("=" * 80 + "\n")

        # Process all structures
        stats = process_all_structures(base_dir="PDB", convert_to_xyz=True)

        # Save statistics to timestamped file
        stats_file = f"xtb_optimization_stats_{stats['timestamp']}.txt"
        with open(stats_file, "w") as f:
            f.write(f"XTB Unconstrained Optimization Statistics\n")
            f.write(f"{'='*50}\n")
            f.write(f"Run Timestamp: {stats['timestamp']}\n")
            f.write(f"Total structures: {stats['total']}\n")
            f.write(f"Skipped (already done): {stats['skipped']}\n")
            f.write(f"Successfully optimized: {stats['success']}\n")
            f.write(f"  - vtight: {stats['vtight']}\n")
            f.write(f"  - crude+vtight: {stats['crude+vtight']}\n")
            f.write(f"  - crude only: {stats['crude']}\n")
            f.write(f"Failed: {stats['failed']}\n")
            f.write(f"{'='*50}\n")

        print(f"\nStatistics saved to: {stats_file}")
        print(f"Master log: xtb_optimization_master_{stats['timestamp']}.log")
        print(f"History file: xtb_optimization_history_{stats['timestamp']}.txt")

    elif args.mode == 'charge-mult':
        # Charge-multiplicity exploration mode
        print("XTB Charge-Multiplicity Exploration Pipeline")
        print("=" * 80)
        print("This script will run XTB optimization for multiple charge-multiplicity states:")
        print("  - q=0, m=1 (neutral, singlet)")
        print("  - q=0, m=5 (neutral, quintet)")
        print("  - q=1, m=2 (cation, doublet)")
        print("  - q=1, m=6 (cation, sextet)")
        print()
        print("RESUMABLE: Charge-mult states already optimized will be skipped automatically.")
        print("           Safe to stop and restart anytime!")
        print()
        print("All output files are timestamped to prevent overwriting.")
        print("=" * 80 + "\n")

        # Parse PDB-ID list
        pdb_ids = parse_pdb_id_list(args.pdb_list)
        if not pdb_ids:
            print("✗ No PDB-IDs found in the provided list!")
            exit(1)

        print(f"Found {len(pdb_ids)} PDB-IDs to process:")
        for pdb_id in pdb_ids:
            print(f"  - {pdb_id}")
        print()

        # Process charge-mult structures
        stats = process_charge_mult_structures(pdb_ids, base_dir="PDB", convert_to_xyz=True)

        # Save statistics to timestamped file
        stats_file = f"xtb_chargemult_stats_{stats['timestamp']}.txt"
        with open(stats_file, "w") as f:
            f.write(f"XTB Charge-Multiplicity Exploration Statistics\n")
            f.write(f"{'='*50}\n")
            f.write(f"Run Timestamp: {stats['timestamp']}\n")
            f.write(f"Total structures: {stats['total_structures']}\n")
            f.write(f"Total calculations: {stats['total_calculations']}\n\n")
            f.write("Results by charge-multiplicity state:\n")
            f.write("-"*50 + "\n")
            for q, m in [(0, 1), (0, 5), (1, 2), (1, 6)]:
                state_key = f"q{q}_m{m}"
                f.write(f"\nq={q}, m={m}:\n")
                f.write(f"  Success: {stats[state_key]['success']}\n")
                f.write(f"  Failed: {stats[state_key]['failed']}\n")
                f.write(f"  Skipped: {stats[state_key]['skipped']}\n")
            f.write(f"\n{'='*50}\n")

        print(f"\nStatistics saved to: {stats_file}")
        print(f"Master log: xtb_chargemult_master_{stats['timestamp']}.log")
        print(f"History file: xtb_chargemult_history_{stats['timestamp']}.txt")
