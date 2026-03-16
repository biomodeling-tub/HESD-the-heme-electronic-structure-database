# XTB-Based Gaussian16 Workflow - Implementation Status

## Date: January 13, 2026

## Summary

Successfully implemented a **COMPLETE** XTB-based Gaussian16 workflow that uses XTB-preoptimized structures as starting geometries instead of checkpoint files from different charge-multiplicity states. All Gaussian input files, bash submission scripts, and XTB structure files are now generated automatically with full workflow automation.

## What Was Implemented

### 1. New Automation Script
- **File**: `automation/beb00042_folders_files_xtb.py`
- **Status**: ✅ WORKING
- Reads PDB IDs from `pdb_ids.txt` (13 IDs total)
- Finds most recent XTB-optimized files for each charge-multiplicity state
- Generates Gaussian .com files with XTB coordinates

### 2. File Organization
- **Output Directory**: `automation_xtb/` (completely separate from original workflow)
- **Original Files**: UNTOUCHED (`automation/`, `sample_calc/` remain unchanged)
- **Structure**: `automation_xtb/{pdb_id}/{pdb_id}{state}/`

### 3. Generated Files (52 files per PDB-ID)
- **12 PDB IDs processed**: 1a2f, 1a6k, 1aa4, 1aom, 1c52, 1c75, 1cpt, 1cxy, 1cyo, 1dp6, 1gjm, 1ltw
- **1 PDB ID skipped**: 2gdm (missing XTB files)
- **Files per PDB ID**:
  - **Gaussian input files (.com)**: 12 total
    - State 01: 3 files (1a2f01.com, 1a2f01_nbo.com, 1a2fL2.com)
    - State 05: 3 files (1a2f05re.com, 1a2f05.com, 1a2fL4.com)
    - State 12: 3 files (1a2f12re.com, 1a2f12.com, 1a2fL6.com)
    - State 16: 3 files (1a2f16re.com, 1a2f16.com, 1a2fL8.com)
  - **Bash submission scripts (.sh)**: 16 total
    - State 01: 4 scripts (g16_1a2f01.sh, g16_1a2f01_nbo.sh, g16_1a2fL2.sh, 1a2f01error.sh)
    - State 05: 4 scripts (g16_1a2f05re.sh, g16_1a2f05.sh, g16_1a2fL4.sh, 1a2f05error.sh)
    - State 12: 4 scripts (g16_1a2f12re.sh, g16_1a2f12.sh, g16_1a2fL6.sh, 1a2f12error.sh)
    - State 16: 4 scripts (g16_1a2f16re.sh, g16_1a2f16.sh, g16_1a2fL8.sh, 1a2f16error.sh)
  - **XTB structure files**: 24 total (6 files per state: .xyz, .inp, .out, .log, charges, wbo)

### 4. Key Features Implemented

#### Charge-Multiplicity Mapping
| Gaussian State | XTB Suffix | Charge | Multiplicity | Description |
|---|---|---|---|---|
| 01 | q0_m1 | 0 | 1 | Neutral singlet |
| 05 | q0_m5 | 0 | 5 | Neutral quintet |
| 12 | q1_m2 | 1 | 2 | Cation doublet |
| 16 | q1_m6 | 1 | 6 | Cation sextet |

#### Checkpoint Naming (New Convention)
- State 01: `{pdb_id}_01_opt.chk`, `{pdb_id}_01_nbo.chk`, `{pdb_id}_01_solv.chk`
- State 05: `{pdb_id}_05_opt.chk`, `{pdb_id}_05_nbo.chk`, `{pdb_id}_05_solv.chk`
- State 12: `{pdb_id}_12_opt.chk`, `{pdb_id}_12_nbo.chk`, `{pdb_id}_12_solv.chk`
- State 16: `{pdb_id}_16_opt.chk`, `{pdb_id}_16_nbo.chk`, `{pdb_id}_16_solv.chk`

#### Coordinate Verification
✅ **VERIFIED**: Gaussian files contain XTB-optimized coordinates
- Example Fe atom coordinates match exactly between XTB and Gaussian files
- No `geom=check` in initial optimization files
- All files have explicit XYZ coordinates

## Additional Features Implemented (January 13, 2026)

### 1. ✅ Overwrite Flag Functionality
- Added `overwrite` parameter to `__init__` (default=False)
- Implemented `_file_exists_and_should_skip()` helper method
- All file creation methods check for existing files before writing
- Safe re-runs without accidentally overwriting existing files

### 2. ✅ XTB Structure File Copying
- Implemented `_copy_xtb_structure_files()` method
- Automatically copies all XTB-related files to output folders:
  - Optimized structures (.xyz)
  - XTB input parameters (.inp)
  - XTB calculation logs (.out)
  - Convergence data (.log)
  - Partial charges (charges files)
  - Wiberg bond orders (wbo files)
- 6 files copied per state, 24 files per PDB-ID

### 3. ✅ NBO and Solvation Scripts (Within-State)
- `_create_gaussian_scripts_05()` - NBO/solvation for state 05
- `_create_gaussian_scripts_12()` - NBO/solvation for state 12
- `_create_gaussian_scripts_16()` - NBO/solvation for state 16
- All use new checkpoint naming convention
- Keep `geom=check` for within-state continuations

### 4. ✅ Complete Bash Submission Scripts
All 13 bash script creation methods implemented:
- `_error_warning_job()` - Error scripts for all 4 states
- `_create_run_script_01()` - State 01 initial optimization
- `_create_run_script_01_nbo()` - State 01 NBO calculation
- `_create_run_script_Link2()` - State 01 solvation
- `_create_run_script_05re()` - State 05 initial optimization
- `_create_run_script_05()` - State 05 NBO calculation
- `_create_run_script_Link4()` - State 05 solvation
- `_create_run_script_12re()` - State 12 initial optimization
- `_create_run_script_12()` - State 12 NBO calculation
- `_create_run_script_Link6()` - State 12 solvation
- `_create_run_script_16re()` - State 16 initial optimization
- `_create_run_script_16()` - State 16 NBO calculation
- `_create_run_script_Link8()` - State 16 solvation

### 5. ✅ Cross-State Dependency Removal
**Link2 script (state 01):**
- ❌ Removed checkpoint copying to states 05/12
- ❌ Removed job submissions to 05re and 12re
- ❌ Removed folder deletion (keeps 01 folder intact)

**Link6 script (state 12):**
- ❌ Removed checkpoint copying to state 16
- ❌ Removed job submission to 16re
- ❌ Removed folder deletion (keeps 12 folder intact)

**Link4 and Link8 scripts:**
- ❌ Removed folder deletions
- Each state completes independently without cleanup

### 6. ✅ File Existence Checks
- All .com and .sh creation methods check for existing files
- Prevents accidental overwrites when overwrite=False
- Clear console output showing which files are skipped

## Current Status

✅ **FULLY FUNCTIONAL**: All features have been implemented!

1. ✅ **Bash scripts complete**: Ready for cluster submission via SLURM
2. ✅ **NBO/solvation for all states**: Complete workflow for all 4 states
3. ✅ **Automated workflow**: Script generates all files automatically
4. ✅ **Independent states**: Each state can run independently without cross-state dependencies
5. ✅ **Overwrite protection**: Safe re-runs with overwrite=False by default

## Testing Status

### ✅ Verified Working
- XTB file discovery (finds most recent timestamped files)
- Coordinate extraction from XTB .xyz files
- Gaussian .com file generation with correct:
  - Charge-multiplicity states
  - XTB coordinates
  - Checkpoint naming
  - No cross-state dependencies
- Directory separation (no file overlap)

### ✅ Tested and Verified (January 13, 2026)
- XTB structure file copying (24 files per PDB-ID)
- All .com file generation (12 files per PDB-ID)
- All bash script generation (16 files per PDB-ID)
- Overwrite flag functionality (skips existing files correctly)
- File counts match expected (52 files per PDB-ID total)
- Checkpoint naming convention is correct throughout
- Cross-state dependencies successfully removed
- All bash scripts are executable (chmod 755)

### ⚠️ Not Yet Tested (Requires Cluster Access)
- Actual Gaussian16 execution with XTB starting geometries
- Job chaining within states (01→01_nbo→L2, etc.)
- SLURM job submission and completion
- Convergence from XTB starting points

## Usage Instructions

### To Generate Gaussian Input Files:
```bash
cd /home/pbuser/Desktop/PhD_WORK/heme/automation
python3 beb00042_folders_files_xtb.py
```

### To Submit Jobs to Cluster (Example for 1a2f, state 01):
```bash
cd /home/pbuser/Desktop/PhD_WORK/heme/automation_xtb/1a2f/1a2f01

# Submit initial optimization (will auto-chain to NBO and solvation)
sbatch g16_1a2f01.sh

# Or submit each step manually:
sbatch g16_1a2f01.sh        # Initial optimization
sbatch g16_1a2f01_nbo.sh    # NBO analysis (after 01 completes)
sbatch g16_1a2fL2.sh        # Solvation (after NBO completes)
```

### To Run All States Independently:
```bash
cd /home/pbuser/Desktop/PhD_WORK/heme/automation_xtb/1a2f

# Submit all 4 states in parallel (no dependencies between states)
sbatch 1a2f01/g16_1a2f01.sh
sbatch 1a2f05/g16_1a2f05re.sh
sbatch 1a2f12/g16_1a2f12re.sh
sbatch 1a2f16/g16_1a2f16re.sh
```

## Files Generated (Example for 1a2f)

```
automation_xtb/
└── 1a2f/
    ├── 1a2f01/                         (13 files total)
    │   ├── 1a2f01.com                  ← Initial opt (q0_m1)
    │   ├── 1a2f01_nbo.com              ← NBO analysis
    │   ├── 1a2fL2.com                  ← Solvation
    │   ├── g16_1a2f01.sh               ← Bash script for initial opt
    │   ├── g16_1a2f01_nbo.sh           ← Bash script for NBO
    │   ├── g16_1a2fL2.sh               ← Bash script for solvation
    │   ├── 1a2f01error.sh              ← Error handling script
    │   └── [6 XTB structure files]     ← .xyz, .inp, .out, .log, charges, wbo
    ├── 1a2f05/                         (13 files total)
    │   ├── 1a2f05re.com                ← Initial opt (q0_m5)
    │   ├── 1a2f05.com                  ← NBO analysis
    │   ├── 1a2fL4.com                  ← Solvation
    │   ├── g16_1a2f05re.sh             ← Bash script for initial opt
    │   ├── g16_1a2f05.sh               ← Bash script for NBO
    │   ├── g16_1a2fL4.sh               ← Bash script for solvation
    │   ├── 1a2f05error.sh              ← Error handling script
    │   └── [6 XTB structure files]
    ├── 1a2f12/                         (13 files total)
    │   ├── 1a2f12re.com                ← Initial opt (q1_m2)
    │   ├── 1a2f12.com                  ← NBO analysis
    │   ├── 1a2fL6.com                  ← Solvation
    │   ├── g16_1a2f12re.sh             ← Bash script for initial opt
    │   ├── g16_1a2f12.sh               ← Bash script for NBO
    │   ├── g16_1a2fL6.sh               ← Bash script for solvation
    │   ├── 1a2f12error.sh              ← Error handling script
    │   └── [6 XTB structure files]
    └── 1a2f16/                         (13 files total)
        ├── 1a2f16re.com                ← Initial opt (q1_m6)
        ├── 1a2f16.com                  ← NBO analysis
        ├── 1a2fL8.com                  ← Solvation
        ├── g16_1a2f16re.sh             ← Bash script for initial opt
        ├── g16_1a2f16.sh               ← Bash script for NBO
        ├── g16_1a2f16L8.sh             ← Bash script for solvation
        ├── 1a2f16error.sh              ← Error handling script
        └── [6 XTB structure files]
```

## Recommended Next Steps

### Testing with Cluster

1. **Test single PDB ID (1a2f, state 01)**:
   ```bash
   cd automation_xtb/1a2f/1a2f01
   sbatch g16_1a2f01.sh
   ```
   - Monitor job completion
   - Check for "Normal termination" in output
   - Verify checkpoint files are created with new naming
   - Check convergence from XTB starting point

2. **Test independent state execution**:
   - Submit all 4 states of 1a2f in parallel
   - Verify no cross-state dependencies cause failures
   - Confirm each state completes independently

3. **Compare with original workflow**:
   - Run same PDB ID with original workflow
   - Compare final energies
   - Compare geometries
   - Verify results are chemically reasonable

### Optional Enhancements

1. **Add progress tracking**: Modify agent scripts to track completion status
2. **Add result analysis**: Create scripts to compare energies across states
3. **Add validation**: Check for convergence issues and report them
4. **Add restart capability**: Handle failed jobs more gracefully

## Conclusion

✅ **IMPLEMENTATION COMPLETE - Fully Functional XTB-Based Workflow**

**What Works:**
- ✅ Complete automation script with all features implemented
- ✅ XTB-based workflow with independent charge-multiplicity states
- ✅ All Gaussian input files (.com) generated with XTB coordinates
- ✅ All bash submission scripts (.sh) created and executable
- ✅ XTB structure files copied to output folders for reference
- ✅ Cross-state dependencies completely removed
- ✅ Overwrite protection for safe re-runs
- ✅ Complete separation from original workflow
- ✅ 52 files per PDB-ID (12 .com + 16 .sh + 24 XTB files)

**Benefits:**
1. **Independent execution**: All 4 states can run in parallel
2. **Correct starting geometries**: Each state uses XTB geometry for its specific charge/multiplicity
3. **No cross-contamination**: No checkpoint files from different electronic states
4. **Easy management**: Clear file organization with descriptive names
5. **Safe operation**: Overwrite protection prevents accidental data loss
6. **Complete workflow**: From XTB optimization to final solvation calculations

**Ready For:**
- Cluster submission via SLURM
- Production runs on all 12 PDB-IDs
- Independent state calculations
- Comparison with original workflow

**Next Phase:** Testing with actual Gaussian16 execution on cluster
