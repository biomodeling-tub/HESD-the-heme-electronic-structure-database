#!/bin/bash

# Verification script to test pipeline manager logic

echo "=========================================="
echo "Pipeline Manager Logic Verification"
echo "=========================================="
echo ""

# Check if script exists
if [[ ! -f "./gaussian_pipeline_manager.sh" ]]; then
    echo "ERROR: gaussian_pipeline_manager.sh not found in current directory"
    exit 1
fi

echo "✓ Script found: gaussian_pipeline_manager.sh"
echo ""

# Check for SLURM integration
echo "=== Checking SLURM Integration ==="
if grep -q "is_job_running" gaussian_pipeline_manager.sh; then
    echo "✓ SLURM queue checking function found"
else
    echo "✗ SLURM queue checking function NOT found"
fi

if grep -q "squeue -u \$USER" gaussian_pipeline_manager.sh; then
    echo "✓ SLURM queue command (squeue) found"
else
    echo "✗ SLURM queue command NOT found"
fi
echo ""

# Check that old restrictive checks were removed
echo "=== Checking for Old Restrictive Logic ==="
restrictive_check_count=$(grep 'if \[\[ ! -f.*\.log.*\]\]; then' gaussian_pipeline_manager.sh | wc -l)
if [[ $restrictive_check_count -eq 0 ]]; then
    echo "✓ Old restrictive log file checks removed"
else
    echo "⚠ Found $restrictive_check_count restrictive log file checks (should be 0)"
fi
echo ""

# Check SLURM queue
echo "=== Checking Your SLURM Queue ==="
if command -v squeue &> /dev/null; then
    running_jobs=$(squeue -u $USER 2>/dev/null | wc -l)
    running_jobs=$((running_jobs - 1))  # Subtract header line
    echo "✓ SLURM command available"
    echo "  Currently running jobs: $running_jobs"

    if [[ $running_jobs -gt 0 ]]; then
        echo ""
        echo "  Your running jobs:"
        squeue -u $USER -o "  %.18i %.9P %.50j %.8T"
    fi
else
    echo "⚠ SLURM command (squeue) not available"
    echo "  Note: Script will still run but can't check queue"
fi
echo ""

# Check for PDB directories
echo "=== Checking for PDB Directories ==="
pdb_count=0
for dir in */; do
    pdb=$(basename "$dir")
    if [[ ${#pdb} -eq 4 ]]; then
        pdb_count=$((pdb_count + 1))
        echo "  Found: $pdb/"
    fi
done

if [[ $pdb_count -gt 0 ]]; then
    echo "✓ Found $pdb_count PDB director(ies)"
else
    echo "⚠ No 4-character PDB directories found"
fi
echo ""

# Test run the script
echo "=== Test Run ==="
echo "Running: ./gaussian_pipeline_manager.sh"
echo ""
./gaussian_pipeline_manager.sh
echo ""

echo "=========================================="
echo "Verification Complete!"
echo "=========================================="
echo ""
echo "Summary:"
echo "  - Script checks for successful termination in log files"
echo "  - Jobs are submitted only if not already in SLURM queue"
echo "  - Failed jobs will be automatically resubmitted"
echo "  - States 05 and 12 can run in parallel"
echo ""
echo "To monitor your pipeline:"
echo "  1. Check SLURM queue: squeue -u \$USER"
echo "  2. Watch log files: tail -f {pdb}/{pdb}{calc}/{pdb}{calc}.log"
echo "  3. Run script periodically or via cron"
